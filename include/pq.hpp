#ifndef PQ_HPP
#define PQ_HPP

#include "core.hpp"
#include <array>

namespace utils {
template <typename T, typename Cmp = std::less<T>> class pq {
  struct pq_node {
    pq_node() = default;

    pq_node(const pq_node &) = delete;
    pq_node(pq_node &&other) = default;

    template <typename... Args>
    pq_node(Args &&...args) : value_(std::forward<Args>(args)...) {}

    ~pq_node() { destruct(child_); }

    pq_node &operator=(const pq_node &) = delete;
    pq_node &operator=(pq_node &&other) = default;

    template <typename... Args> static pq_node *construct(Args &&...args) {
      pq_node *node = new pq_node(std::forward<Args>(args)...);
      node->right_ = node;
      node->left_ = node;
      return node;
    }

    static void destruct(pq_node *node) {
      if (node) {
        pq_node *walk = node->right_;
        while (walk != node) {
          pq_node *tmp = walk;
          walk = walk->right_;
          delete tmp;
        }
        delete node;
      }
    }

    bool singleton() const { return this == right_; }

    /**
     * @brief detach the node making it a singlteon.
     *
     * @note does not affect the child.
     * @note resets the parent
     */
    void detach() {
      if (singleton()) {
        return;
      }
      right_->left_ = left_;
      left_->right_ = right_;
      left_ = right_ = this;
    }

    void attach(pq_node *node) {
      if (node) {
        node->left_->right_ = right_;
        right_->left_ = node->left_;
        right_ = node;
        node->left_ = this;
      }
    }

    void absorb_child() {
      if (child_) {
        pq_node *child = child_;
        child_ = nullptr;
        attach(child);
      }
    }

    void adopt_child(pq_node *child) {
      if (child) {
        degree_++;
        child->detach();
        if (child_) {
          child_->attach(child);
        } else {
          child_ = child;
        }
      }
    }

    void log(const char *str = "") const {
      logger::debug(log_node(str, this));
      logger::flush();
    }

    static std::string log_node(const char *str, const pq_node *ptr) {
      std::ostringstream oss;
      oss << str << ptr->degree_ << ":[";
      const pq_node *walk = ptr;
      do {
        oss << walk->value_;
        if (walk->child_) {
          oss << log_node(" (", walk->child_) << ")";
        }
        walk = walk->right_;
        if (walk == ptr) {
          break;
        }
        oss << " -> ";
      } while (true);
      oss << "]";
      return oss.str();
    }

    T value_;
    pq_node *left_ = nullptr;
    pq_node *child_ = nullptr;
    pq_node *right_ = nullptr;
    std::size_t degree_ = 0;
  };

public:
  pq(const pq &) = delete;
  pq(pq &&) = default;

  pq &operator=(const pq &) = delete;
  pq &operator=(pq &&) = default;

  pq() : pq(Cmp()) {}

  ~pq() {
    if (root) {
      pq_node::destruct(root);
    }
  }

  template <typename Comparator>
  explicit pq(const Comparator &cmp_in) : cmp(cmp_in) {}

private:
  inline bool compare(const pq_node &a, const pq_node &b) const {
    return cmp(a.value_, b.value_);
  }

  inline bool compare(const pq_node *a, const pq_node *b) const {
    assert(a != nullptr, "this should not trigger");
    if (!a) {
      throw std::logic_error("cannot compare a nullptr");
    }
    assert(b != nullptr, "this should not trigger");
    if (!b) {
      throw std::logic_error("cannot compare a nullptr");
    }
    return compare(*a, *b);
  }

  const Cmp &cmp;

  pq_node *root = nullptr;
  std::size_t count = 0;
  static inline constexpr std::uint16_t PQ_HEIGHT_UPPER_BOUND = 64;
  std::array<pq_node *, PQ_HEIGHT_UPPER_BOUND> subtrees;

public:
  template <typename... Args> void emplace(Args &&...args) {
    count++;
    pq_node *node = pq_node::construct(std::forward<Args>(args)...);
    if (root == nullptr) {
      root = node;
    } else {
      if (compare(node, root)) {
        std::swap(root, node);
      }
      root->attach(node);
    }
  }

  void push(T val) { emplace(std::move(val)); }

  bool empty() const { return root == nullptr; }

  const T &top() const {
    if (empty()) {
      throw std::runtime_error("cannot pop an empty priority queue.");
    }
    return root->value_;
  }

  void pop() {
    if (empty()) {
      throw std::runtime_error("cannot pop an empty priority queue.");
    }
    count--;
    root->absorb_child();
    if (root->singleton()) {
      delete root;
      root = nullptr;
      return;
    }
    pq_node *new_root = root->right_;
    root->detach();
    delete root;
    root = new_root;
    subtrees.fill(nullptr);
    consolidate();
  }

private:
  void consolidate() {
    pq_node *new_root = root;

    while (true) {
      std::size_t d = root->degree_;
      if (d >= PQ_HEIGHT_UPPER_BOUND) {
        throw std::runtime_error("pq grew too much!");
      }
      pq_node *u = subtrees[d];
      if (u == root) {
        break;
      } else if (u) {
        subtrees[d] = nullptr;
        u->detach();
        if (compare(root, u)) {
          root->adopt_child(u);
        } else {
          pq_node *left = root->left_;
          root->detach();
          left->attach(u);
          u->adopt_child(root);
          root = u;
        }
      } else {
        if (compare(root, new_root)) {
          new_root = root;
        }
        subtrees[root->degree_] = root;
        root = root->right_;
      }
    }
    root = new_root;
  }
};
} // namespace utils
#endif
