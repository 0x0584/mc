#include "reader.hpp"

namespace mc {
namespace order_vertices {
order_by_neighbours by_neighbours;
order_by_neighbourhood by_neighbourhood;
order_by_degeneracy by_degeneracy;

adjacency_list_sorted
order_by_neighbours::operator()(adjacency_list &adj_lst_orig) const {
  adjacency_list_sorted orig;
  orig.reserve(adj_lst_orig.size());
  for (auto &[v, neighbours] : adj_lst_orig) {
    orig.emplace_back(v, std::move(neighbours));
  }
  std::sort(std::execution::par_unseq, orig.begin(), orig.end(),
            [](const auto &v, const auto &u) {
              return v.second.size() < u.second.size();
            });
  return orig;
}

adjacency_list_sorted
order_by_neighbourhood::operator()(adjacency_list &adj_lst_orig) const {
  std::unordered_map<vertex, size_type> neighbours_degree;
  neighbours_degree.max_load_factor(0.5);
  neighbours_degree.reserve(adj_lst_orig.size());

  adjacency_list_sorted orig;
  orig.reserve(adj_lst_orig.size());
  for (auto &[v, neighbours] : adj_lst_orig) {
    for (const auto &u : neighbours) {
      neighbours_degree[v] += adj_lst_orig.at(u).size();
    }
    orig.emplace_back(v, std::move(neighbours));
  }

  std::sort(std::execution::par_unseq, orig.begin(), orig.end(),
            [&](const auto &v, const auto &u) {
              return (v.second.size() < u.second.size() ||
                      (v.second.size() == u.second.size() &&
                       neighbours_degree.at(v.first) <
                           neighbours_degree.at(u.first)));
            });

  return orig;
}

adjacency_list_sorted
order_by_degeneracy::operator()(adjacency_list &adj_lst_orig) const {
  std::unordered_map<vertex, size_type> degrees;
  degrees.max_load_factor(0.5);
  degrees.reserve(adj_lst_orig.size());

  for (const auto &[v, neighbours] : adj_lst_orig) {
    degrees.emplace(v, neighbours.size());
  }

  // MAKE_TIMER(degeneracy);
  for (size_type todo = degrees.size(), level = 1; todo; ++level) {
    std::vector<vertex> current;
    current.reserve(todo);
    for (const auto &[v, degree] : degrees) {
      if (degree == level) {
        current.emplace_back(v);
      }
    }
    for (std::vector<vertex> next; not current.empty();
         current = std::move(next)) {
      todo -= current.size();
      next.reserve(current.size());
      for (const auto &v : current) {
        for (const auto &u : adj_lst_orig.at(v)) {
          size_type &deg = degrees[u];
          if (deg > level && --deg == level) {
            next.emplace_back(v);
          }
        }
      }
    }
  }
  // degeneracy.stop();

  adjacency_list_sorted orig;
  orig.reserve(adj_lst_orig.size());
  for (auto &[v, neighbours] : adj_lst_orig) {
    orig.emplace_back(v, std::move(neighbours));
  }

  // MAKE_TIMER(sort);
  std::sort(std::execution::par_unseq, orig.begin(), orig.end(),
            [&degrees](const auto &u, const auto &v) {
              return degrees.at(u.first) < degrees.at(v.first);
            });

  return orig;
}
} // namespace order_vertices

template <typename Type>
static inline bool read_pair(std::istream &in, Type &lvalue, Type &rvalue) {
  static std::size_t num_line = 1;
  std::string line;
  if (not std::getline(in, line)) {
    num_line = 1;
    return false;
  }
  if (std::istringstream str(line); not(str >> lvalue && str >> rvalue)) {
    std::string str_err = std::to_string(num_line) + " is ill-formated";
    throw std::runtime_error(str_err.c_str());
  }
  num_line++;
  return true;
}

reader::reader(const vertex_order &order, std::istream &in) {
  MAKE_TIMER(reader);

  size_type num_v, num_e;
  if (not read_pair(in, num_v, num_e)) {
    return;
  }

  assert(num_v != 0);
  assert(num_e != 0);

  g = graph(num_v, num_e);

  adjacency_list adj_lst_orig;
  adj_lst_orig.max_load_factor(0.5);
  adj_lst_orig.reserve(num_v);

  size_type n_edges = num_e;
  for (vertex u, v; read_pair(in, u, v) && n_edges; --n_edges) {
    adj_lst_orig[u].emplace(v);
    adj_lst_orig[v].emplace(u);
  }

  assert(n_edges == 0);
  assert(adj_lst_orig.size() == num_v);

#ifndef NDEBUG
  if (adj_lst_orig.size() != num_v) {
    std::cerr << "-- Nmber of Vertices mismatched " << adj_lst_orig.size()
              << "/" << num_v << "\n";
  }
  if (n_edges != 0) {
    std::cerr << "-- Number of Edges mismatched " << (num_e - n_edges) << "/"
              << num_e << "\n";
  }
#endif

  std::cerr << "Read graph(" << adj_lst_orig.size() << " vertices and " << num_e
            << " edges)\n";

  adjacency_list_sorted adj_lst_sorted =
      std::visit(order_vertices::dispatcher(adj_lst_orig), order);

  MAKE_TIMER(prepare);

  v2k.max_load_factor(0.5);
  v2k.reserve(num_v);
  k2v.reserve(num_v);
  for (const auto &[v, neighbours] : adj_lst_sorted) {
    v2k.emplace(v, g.adj_lst.size());
    k2v.emplace_back(v);
    g.adj_mtx.emplace_back(std::vector<bool>(num_v));
    g.adj_lst.emplace_back(std::vector<key>{});
    g.adj_lst.back().reserve(neighbours.size());
  }

  STOP_TIMER(prepare);

  MAKE_TIMER(fill);
  for (const auto &[v, neighbours_sorted] : adj_lst_sorted) {
    const key &v_key = v2k.at(v);
    std::vector<bool> &v_adj_mtx = g.adj_mtx[v_key];
    std::vector<key> &neighbours = g.adj_lst[v_key];
    for (const auto &u : neighbours_sorted) {
      const key &u_key = v2k.at(u);
      v_adj_mtx[u_key] = true;
      neighbours.emplace_back(u_key);
    }
    std::sort(std::execution::par_unseq, neighbours.begin(), neighbours.end());
  }
}

std::ostream &operator<<(std::ostream &os, const reader &gr) {
  const auto vertices_per_line = 25;
  const auto tab = [tab_size = 2u](unsigned num_tabs) {
    return std::string(num_tabs * tab_size, ' ');
  };

  for (size_type v = 0; v < gr.g.num_vertices(); ++v) {
    os << gr.k2v[v] << ": {";
    const auto &v_adj = gr.g.neighbours(v);
    for (size_type u = 0; u < v_adj.size(); ++u) {
      if (u % vertices_per_line == 0) {
        os << "\n" << tab(1);
      }
      os << gr.k2v[v_adj[u]] << " ";
    }
    os << "\n}\n\n";
  }
  return os;
}
} // namespace mc
