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

    static void *operator new(std::size_t sz) {
      return memory::pool()->allocate(sz, alignof(pq_node));
    }

    static void operator delete(void *ptr) noexcept {
      memory::pool()->deallocate(ptr, sizeof(pq_node), alignof(pq_node));
    }

    static void operator delete(void *ptr, std::size_t sz) noexcept {
      memory::pool()->deallocate(ptr, sz, alignof(pq_node));
    }

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
        pq_node *node_head = node;
        pq_node *node_tail = node->left_;
        left_->right_ = node_head;
        node_head->left_ = left_;
        node_tail->right_ = this;
        left_ = node_tail;
      }
    }

    void absorb_child() {
      if (child_) {
        attach(child_);
        child_ = nullptr;
        degree_--;
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
      oss << str << "[";
      const pq_node *walk = ptr;
      do {
        oss << walk->value_ << ":" << walk->degree_;
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
  explicit pq(const Comparator &cmp_in) : cmp(cmp_in) {
    subtrees.fill(nullptr);
  }

private:
  inline bool compare(const pq_node &a, const pq_node &b) const {
    return cmp(a.value_, b.value_);
  }

  inline bool compare(const pq_node *a, const pq_node *b) const {
    if (!a) {
      throw std::logic_error("cannot compare a nullptr");
    }
    if (!b) {
      throw std::logic_error("cannot compare a nullptr");
    }
    return compare(*a, *b);
  }

  const Cmp &cmp;

  pq_node *root = nullptr;
  std::size_t count = 0;
  static inline constexpr std::uint16_t PQ_HEIGHT_UPPER_BOUND = 128;
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

  std::size_t size() const { return count; }

  const T &top() const {
    if (empty()) {
      throw std::runtime_error("cannot get top of an empty priority queue.");
    }
    return root->value_;
  }

  // void log_root(const char *str) const {
  //   if (root)
  //     root->log(str);
  //   else
  //     logger::warn(str, "logging an empty root!");
  // }

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
    consolidate();
  }

private:
  void consolidate() {
    pq_node sentinel;
    sentinel.left_ = root->left_;
    sentinel.right_ = root;
    sentinel.left_->right_ = &sentinel;
    sentinel.right_->left_ = &sentinel;

    pq_node *walk = sentinel.right_;
    do {
      pq_node *current = walk;
      walk = walk->right_;

      current->detach();

      std::size_t d = current->degree_;
      while (subtrees[d]) {
        pq_node *other = subtrees[d];
        subtrees[d] = nullptr;
        if (compare(current, other)) {
          current->adopt_child(other);
        } else {
          other->adopt_child(current);
          current = other;
        }
        d++;
      }

      subtrees[d] = current;
    } while (walk != &sentinel);

    sentinel.detach();

    root = nullptr;
    for (pq_node *&node : subtrees) {
      if (!node)
        continue;
      if (!root) {
        root = node;
        root->left_ = root->right_ = node;
      } else {
        root->attach(node);
        if (compare(node, root))
          root = node;
      }
      node = nullptr;
    }
  }
};
} // namespace utils
#endif
