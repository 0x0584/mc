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
#include "thread.hpp"
#include <list>

namespace mc {
void graph::print() const {
  std::ostringstream oss;
  oss << (undirected ? "Undirected" : "Directed")
      << " Graph (vertices=" << A.size() << ", edges=" << _edge_count << ")\n";
  for (const auto &[v, neighs] : A) {
    oss << v << " { ";
    for (vertex u : neighs) {
      oss << u << " ";
    }
    oss << "}\n";
  }
  log::print(oss.str());
}

bool graph::add_edge_undirected(vertex u, vertex v, std::size_t edge_set_size) {
  neighbours_set &u_neighs = A[u];
  if (u_neighs.empty()) {
    u_neighs = neighbours_set{&edges_pool};
    u_neighs.reserve(edge_set_size);
  }
  if (u_neighs.emplace(v).second) {
    ++_edge_count;
    if (undirected) {
      neighbours_set &v_neighs = A[v];
      if (v_neighs.empty()) {
        v_neighs = neighbours_set{&edges_pool};
        v_neighs.reserve(edge_set_size);
      }
      if (not v_neighs.emplace(u).second) {
        log::info("redundant edge from", v, "to", u);
        return false;
      }
    }
  } else {
    log::info("redundant edge from", u, "to", v);
    return false;
  }
  return true;
}

bool graph_builder::read_single_vertex(graph::vertex &w,
                                       std::string::iterator &it,
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
    graph::vertex v = 0;
    while (left != right && *left >= '0' && *left <= '9') {
      v = v * 10 + unsigned(*left - '0');
      ++left;
    }
    return v;
  };

  if (seek_token(it); it != end) {
    auto start = it;
    read_vertex(it);
    w = parse_vertex(start, it);
    //    w = std::atol(std::string{i, it}.c_str());
    return true;
  } else {
    return false;
  }
}

graph graph_builder::build(std::pmr::monotonic_buffer_resource &vertices_pool,
                           std::pmr::monotonic_buffer_resource &edges_pool) {
  auto begin = std::chrono::high_resolution_clock::now();
  std::mutex graph_mtx;
  graph G(vertices_pool, edges_pool);
  G.A.reserve(feed.num_vertices());
  G.undirected = args::undirected;

  std::atomic_uint tasks_count{0};
  std::mutex tasks_mtx;
  std::condition_variable tasks_barrier;
  const std::size_t estimate_num_edges =
      2.15 * feed.num_edges() / feed.num_vertices();
  const std::size_t edges_per_chunk =
      1.50 * feed.num_edges() / feed.estimate_chunks();
  log::info("edges per task", edges_per_chunk, "/ num chunks",
            feed.estimate_chunks(), "edge set estimate", estimate_num_edges);

  do {
    {
      std::unique_lock<std::mutex> lock(tasks_mtx);
      tasks_barrier.wait(lock, [&tasks_count] {
        return tasks_count.load(std::memory_order_acquire) <
               thread::num_threads;
      });
      tasks_count.fetch_add(1, std::memory_order_release);
    }
    std::string buffer = feed.read_chunk();
    Q.emplace_back(std::async(
        [&, this](std::string buffer) {
          // std::pmr::list<std::pair<graph::vertex, graph::vertex>>
          // local_edges;
          std::pmr::vector<std::pair<graph::vertex, graph::vertex>> local_edges;
          local_edges.reserve(edges_per_chunk);
          std::string::iterator it = buffer.begin();
          int vertices_read = 2;
          for (graph::vertex u = graph::nil_vertex, v = graph::nil_vertex;
               vertices_read == 2;) {
            vertices_read = read_single_vertex(u, it, buffer.end()) +
                            read_single_vertex(v, it, buffer.end());
            if (vertices_read != 2) {
              break;
            } else if (u == v) {
              log::info("found cycle for vertex", u);
              continue;
            }

            local_edges.emplace_back(u, v);
          }
          // log::info("edges read:", local_edges.size());
          {
            std::unique_lock lock(graph_mtx);
            for (auto [u, v] : local_edges) {
              G.add_edge_undirected(u, v, estimate_num_edges);
            }
          }
          {
            std::unique_lock<std::mutex> lock(tasks_mtx);
            tasks_count.fetch_sub(1, std::memory_order_release);
          }
          tasks_barrier.notify_one();
        },
        std::move(buffer)));
  } while (feed);
  std::for_each(Q.begin(), Q.end(), [](std::future<void> &f) { f.get(); });

  auto end = std::chrono::high_resolution_clock::now();

  log::info("Graph with", feed.num_vertices(), "vertices and", G._edge_count,
            "edges was read in", log::time_diff(begin, end, log::bold));

  assert(G.A.size() == feed.num_vertices());
  assert(G._edge_count == feed.num_edges());

  return G;
}

} // namespace mc
