// graph.cpp
//
// Copyright (C) 2024  0x0584 (Anas)
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
// USA.

#include "graph.hpp"
#include "core.hpp"
#include "thread.hpp"

namespace mc {
void graph::print() const {
  std::ostringstream oss;
  oss << " Graph (vertices=" << A.size() << ", edges=" << _edge_count << ")\n";
  for (auto u = 0ul; u < V.size(); ++u) {
    oss << V[u] << " { ";
    for (key v : A[u]) {
      oss << V[v] << " ";
    }
    oss << "}\n";
  }
  logger::print(oss.str());
}

bool graph_builder::read_single_vertex(vertex &w, std::string::iterator &it,
                                       std::string::iterator end) {
  auto seek_token = [end](auto &it) {
    while (it != end && (*it == feed::sep || *it == feed::deli)) {
      ++it;
    }
  };
  auto read_vertex = [end](auto &it) {
    int limit = 0;
    while (limit++ < 12 && it != end && *it != feed::sep && *it != feed::deli) {
      ++it;
    }
    if (limit > 12)
      throw std::runtime_error("vertex key is too large! abort parsing.");
  };
  auto parse_vertex = [](auto &left, auto &right) {
    vertex::value_type v = 0;
    while (left != right && *left >= '0' && *left <= '9') {
      v = v * 10 + unsigned(*left - '0');
      ++left;
    }
    return vertex(v);
  };

  if (seek_token(it); it != end) {
    auto start = it;
    read_vertex(it);
    w = parse_vertex(start, it);
    return true;
  } else {
    return false;
  }
}

graph graph_builder::build(bool undirected) {
  auto begin = std::chrono::high_resolution_clock::now();

  logger::debug("Number of chunks", feed.estimate_chunks());
  logger::debug("Reading", feed.edges_per_chunk(),
                "edges per chunk (estimate)");
  logger::debug("Number of edges per vertex set estimate",
                feed.estimate_num_edges());

  std::mutex graph_mtx;

  std::pmr::unordered_map<vertex, std::pmr::unordered_set<vertex>> adj_lst(
      memory::pool());
  adj_lst.reserve(feed.num_vertices());

  graph G;

  G.V = std::pmr::vector<vertex>(memory::pool());
  G.V.reserve(feed.num_vertices());

  G.vertex_to_key = std::pmr::unordered_map<vertex, graph::key>(memory::pool());
  G.vertex_to_key.reserve(feed.num_vertices());

  auto get_neighbours =
      [&G, &adj_lst, edge_size_estimate = feed.estimate_num_edges()](
          const vertex &u) -> std::pmr::unordered_set<vertex> & {
    std::pmr::unordered_set<vertex> &neighs = adj_lst[u];
    if (neighs.empty()) {
      G.vertex_to_key[u] = G.V.size();
      G.V.emplace_back(u);
      neighs = std::pmr::unordered_set<vertex>(memory::pool());
      neighs.reserve(edge_size_estimate);
    }
    return neighs;
  };

  thread::pool pool(args::num_threads); // XXX switch to a thread pool
  using edge_pair = std::pair<vertex, vertex>;
  std::pmr::vector<std::pmr::vector<edge_pair>> local_edges(
      args::num_threads, std::pmr::vector<edge_pair>(memory::pool()),
      memory::pool());
  for (auto &local : local_edges) {
    local.reserve(feed.edges_per_chunk());
  }

  do {
    pool.exec([&get_neighbours, buff = feed.read_chunk(), &local_edges,
               &graph_mtx, &G, &undirected,
               this](std::uint16_t task_id) mutable {
      auto &local = local_edges[task_id];
      local.clear();

      std::string::iterator it = buff.begin();
      int vertices_read = 2;
      for (vertex u, v; vertices_read == 2;) {
        vertices_read = read_single_vertex(u, it, buff.end()) +
                        read_single_vertex(v, it, buff.end());
        if (vertices_read != 2) {
          break;
        } else if (u == v) {
          logger::error("found cycle for vertex", u);
          continue;
        }
        local.emplace_back(u, v);
      }
      logger::debug("read edges:", local.size());
      std::scoped_lock lock(graph_mtx);
      for (const auto &[u, v] : local) {
        std::pmr::unordered_set<vertex> &u_neighs = get_neighbours(u);
        if (!u_neighs.emplace(v).second) {
          logger::error("redundant edge from", u, "to", v);
        } else {
          G._edge_count++;
          if (undirected) {
            std::pmr::unordered_set<vertex> &v_neighs = get_neighbours(v);
            if (!v_neighs.emplace(u).second) {
              logger::error("redundant edge from", v, "to", u);
            }
          }
        }
      }
    });
  } while (feed);
  pool.join();

  if (G.vertex_count() != feed.num_vertices()) {
    logger::error("Number of vertices read mismatched");
  }
  if (G.edge_count() != feed.num_edges()) {
    logger::error("Number of edges read mismatched");
  }

  G.A.resize(feed.num_vertices());
  for (auto it = adj_lst.cbegin(); it != adj_lst.cend(); ++it) {
    pool.exec([it, &G](std::uint16_t) {
      const auto &[v, neighs] = *it;
      auto &key_neighs = G.A[G.vertex_to_key.at(v)] =
          std::pmr::unordered_set<graph::key>(memory::pool());
      key_neighs.reserve(neighs.size());
      for (const auto &v : neighs) {
        key_neighs.emplace(G.vertex_to_key.at(v));
      }
    });
  }
  pool.join();

  auto end = std::chrono::high_resolution_clock::now();
  logger::info("Graph with", adj_lst.size(), "vertices and", G.edge_count(),
               "edges was read in",
               logger::time_diff(begin, end, logger::bold));

  return G;
}
} // namespace mc
