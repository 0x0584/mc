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

namespace mc {
void graph::print() const {
  std::ostringstream oss;
  oss << (undirected ? "Undirected" : "Directed")
      << " Graph (vertices=" << A.size() << ", edges=" << edge_count << ")\n";
  for (const auto &[v, neighs] : A) {
    oss << v << " { ";
    for (vertex u : neighs) {
      oss << u << " ";
    }
    oss << "}\n";
  }
  log::print(oss.str());
}

bool graph::add_edge_undirected(vertex u, vertex v) {
  neighbours_set &u_neighs = A[u];
  if (u_neighs.empty()) {
    u_neighs = neighbours_set{&edges_pool};
  }
  if (u_neighs.emplace(v).second) {
    ++edge_count;
    if (undirected) {
      neighbours_set &v_neighs = A[v];
      if (v_neighs.empty()) {
        v_neighs = neighbours_set{&edges_pool};
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
  auto skip = [end](std::string::iterator &it) {
    while (it != end && (*it == feed::sep || *it == feed::deli)) {
      ++it;
    }
  };
  auto read_vertex = [end](std::string::iterator &it) {
    int limit = 0;
    while (limit++ < 12 && it != end && *it != feed::sep && *it != feed::deli) {
      ++it;
    }
  };

  if (skip(it); it != end) {
    std::string::iterator i = it;
    read_vertex(it);
    w = std::atol(std::string{i, it}.c_str());
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
  G.undirected = args::undirected;
  do {
    std::string buffer = feed.read_chunk();
    Q.emplace_back(std::async(
        std::launch::deferred,
        [&, this](std::string buffer) {
          std::string::iterator it = buffer.begin();
          int vertices_read = 2;
          for (graph::vertex u, v; vertices_read == 2;) {
            vertices_read = read_single_vertex(u, it, buffer.end()) +
                            read_single_vertex(v, it, buffer.end());
            if (vertices_read != 2) {
              break;
            } else if (u == v) {
              log::info("found cycle for vertex", u);
              continue;
            }

            std::unique_lock lock(graph_mtx);
            G.add_edge_undirected(u, v);
          }
        },
        std::move(buffer)));
  } while (feed);
  std::for_each(Q.begin(), Q.end(), [](std::future<void> &f) { f.get(); });

  auto end = std::chrono::high_resolution_clock::now();

  log::info("Graph with", feed.num_vertices(), "vertices and", G.edge_count,
            "edges was read in", log::time_diff(begin, end, log::bold));

  assert(G.A.size() == feed.num_vertices());
  assert(G.edge_count == feed.num_edges());

  return G;
}

} // namespace mc
