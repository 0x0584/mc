struct fahle_max_clique {
  using index_vector = std::vector<unsigned>;

  struct solve {
    const graph::adjacency_matrix_t &adj_mtx;
    index_vector &sol, clique;
    solve(const graph::adjacency_matrix_t &adj_mtx, index_vector &sol)
      : adj_mtx(adj_mtx), sol(sol) {
      clique.reserve(adj_mtx.size());
      index_vector verts;
      verts.reserve(adj_mtx.size());
      for (unsigned i = 0; i < adj_mtx.size(); ++i) {
        verts.emplace_back(i);
      }
      expand(std::move(verts));
    }

    void expand(index_vector verts) {
      while (not verts.empty() && clique.size() + verts.size() > sol.size()) {
        unsigned v = verts.back();

        index_vector new_verts;
        for (auto u : verts) {
          if (adj_mtx.at(v).at(u)) {
            new_verts.emplace_back(u);
          }
        }

        clique.emplace_back(v);
        if (new_verts.empty() && clique.size() > sol.size()) {
          sol = clique;
        }
        if (not new_verts.empty()) {
          expand(std::move(new_verts));
        }

        clique.pop_back();
        verts.pop_back();
      }
    }
  };

  fahle_max_clique(const graph *g) : g(g) {
    struct solve solver(g->adjacency_matrix(), sol);
  }

  std::vector<vertex_t> solution() const {
    std::vector<vertex_t> sol_verts, verts = g->vertices();
    std::transform(sol.begin(), sol.end(),
                   std::back_inserter(sol_verts),
                   [&verts](unsigned i) { return verts.at(i); });
    return sol_verts;
  }

private:
  const graph *g;
  index_vector sol;
};

std::vector<vertex_t> graph::fahle_max_clique() const {
  struct fahle_max_clique algo(this);
  return algo.solution();
}

using color_t = graph::color_t;
std::pair<std::vector<vertex_t>, std::vector<color_t>> graph::color_welsh_powell() const {
  std::vector<vertex_t> verts = vertices();
  std::vector<color_t> coloring(verts.size());
  for (unsigned u = 0; u < verts.size(); ++u) {
    const auto &adjacent = adjacency_list().at(verts.at(u));
    std::set<color_t> colors;
    for (unsigned v = 0; v < adjacent.size(); ++v) {
      colors.emplace(coloring[v]);
    }
    color_t color = 0;
    while (colors.count(color)) {
      color++;
    }
    coloring[u] = color;
  }
  return std::make_pair(std::move(verts), std::move(coloring));
}

struct tomita_max_clique {
  tomita_max_clique(const graph *g) : colors(g->color_welsh_powell()) {

  }

  std::vector<vertex_t> solution() {
    return {};
  }

private:
  const std::pair<std::vector<vertex_t>, std::vector<color_t>> colors;
  std::vector<std::vector<vertex_t>> color_class;
  std::vector<vertex_t> sol;
};

std::vector<vertex_t> graph::tomita_max_clique() const {
  struct tomita_max_clique algo(this);
  return algo.solution();
}
