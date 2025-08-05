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
  logger::print(oss.str());
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
        logger::info("redundant edge from", v, "to", u);
        return false;
      }
    }
  } else {
    logger::info("redundant edge from", u, "to", v);
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
  std::mutex graph_mtx;
  graph G(vertices_pool, edges_pool);
  G.A.reserve(feed.num_vertices());
  G.undirected = args::undirected;

  logger::debug("Num of chunks", feed.estimate_chunks());
  logger::debug("Reading", feed.edges_per_chunk(),
                "edges per chunk (estimate)");
  logger::debug("Number of edges per vertex set estimate",
                feed.estimate_num_edges());

  thread::pool pool(args::num_threads); // XXX switch to a thread pool
  std::vector<std::vector<std::pair<graph::vertex, graph::vertex>>> local_edges(
      args::num_threads);
  for (auto &local : local_edges) {
    local.reserve(feed.edges_per_chunk());
  }

  auto begin = std::chrono::high_resolution_clock::now();
  do {
    pool.exec([buff = feed.read_chunk(), &local_edges, &G, &graph_mtx,
               this](std::uint16_t task_id) mutable {
      auto &local = local_edges[task_id];
      local.clear();
      std::string::iterator it = buff.begin();
      int vertices_read = 2;
      for (graph::vertex u = graph::nil_vertex, v = graph::nil_vertex;
           vertices_read == 2;) {
        vertices_read = read_single_vertex(u, it, buff.end()) +
                        read_single_vertex(v, it, buff.end());
        if (vertices_read != 2) {
          break;
        } else if (u == v) {
          logger::warn("found cycle for vertex", u);
          continue;
        }
        local.emplace_back(u, v);
      }
      logger::debug("edges read:", local_edges.size());
      std::unique_lock lock(graph_mtx);
      for (auto [u, v] : local) {
        G.add_edge_undirected(u, v, feed.estimate_num_edges());
      }
    });
  } while (feed);
  pool.join();
  auto end = std::chrono::high_resolution_clock::now();

  logger::info("Graph with", feed.num_vertices(), "vertices and", G._edge_count,
               "edges was read in",
               logger::time_diff(begin, end, logger::bold));

  if (G.A.size() != feed.num_vertices())
    logger::warn("Number of vertices read mismatched");
  if (G._edge_count != feed.num_edges())
    logger::warn("Number of edges read mismatched");

  assert(G.A.size() == feed.num_vertices());
  assert(G._edge_count == feed.num_edges());

  return G;
}
} // namespace mc
