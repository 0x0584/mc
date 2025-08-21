#ifndef PQ_HPP
#define PQ_HPP

#include "core.hpp"
#include <array>

namespace utils {
template <typename T, typename Cmp = std::less<T>> class pq {
public:
  ~pq() = default;

  pq(const pq &) = delete;
  pq(pq &&) = default;

  pq &operator=(const pq &) = delete;
  pq &operator=(pq &&) = default;

  pq() : pq(Cmp()) {}

  template <typename Comparator>
  explicit pq(const Comparator &cmp_in) : cmp(cmp_in) {}

private:
  struct pq_node : std::enable_shared_from_this<pq_node> {
    pq_node() = default;

    pq_node(const pq_node &) = delete;
    pq_node(pq_node &&other) = default;

    template <typename... Args>
    pq_node(Args &&...args) : value_(std::forward<Args>(args)...), degree_(0) {}

    ~pq_node() = default;

    pq_node &operator=(const pq_node &) = delete;
    pq_node &operator=(pq_node &&other) = default;

    template <typename... Args>
    static std::shared_ptr<pq_node> construct(Args &&...args) {
      auto node = gc.make_shared<pq_node>(std::forward<Args>(args)...);
      node->right_ = node;
      node->left_ = node;
      return node;
    }

    bool singleton() const { return this == right_.get(); }

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
      right()->left(left_);
      left()->right(right_);
      right_ = this->shared_from_this();
      left_ = right_;
    }

    void attach(std::shared_ptr<pq_node> node) {
      if (!node) {
        return;
      }
      right()->left(node->left());
      node->left()->right(std::move(right()));
      node->left(this->shared_from_this());
      right(std::move(node));
    }

    void absorb_child() {
      if (!child_) {
        return;
      }
      attach(std::move(child_));
    }

    void adopt_child(std::shared_ptr<pq_node> child) {
      child->detach();
      if (!child_) {
        child_ = std::move(child);
      } else {
        child_->attach(std::move(child));
      }
      degree_++;
    }

    void log(const char *str = "") {
      logger::debug(log_node(str, this->weak_from_this()));
    }

  private:
    std::string log_node(const char *str, const std::shared_ptr<pq_node> &ptr) {
      std::ostringstream oss;
      oss << str << ptr->degree() << ":[";
      std::weak_ptr<pq_node> walk = ptr;
      do {
        auto walk_ptr = walk.lock();
        oss << walk_ptr->value().str();
        if (walk_ptr->child_) {
          oss << child_->log_node(" (", walk_ptr->child_) << ")";
        }
        walk = walk_ptr->right();
        if (walk.lock() != ptr)
          oss << " -> ";
      } while (walk.lock() != ptr);
      oss << "]";
      return oss.str();
    }

    std::weak_ptr<pq_node> left_;
    std::shared_ptr<pq_node> child_;
    std::shared_ptr<pq_node> right_;
    T value_;
    std::size_t degree_;

  public:
    const T &value() const { return value_; }

    const std::shared_ptr<pq_node> &right() const { return right_; }
    void right(std::shared_ptr<pq_node> ptr) { right_ = std::move(ptr); }

    std::shared_ptr<pq_node> left() const { return left_.lock(); }
    void left(std::weak_ptr<pq_node> ptr) { left_ = std::move(ptr); }

    const std::shared_ptr<pq_node> &child() const { return child_; }
    void child(std::shared_ptr<pq_node> ptr) { child_ = std::move(ptr); }

    std::size_t degree() const { return degree_; }
  };

  inline bool compare(const pq_node &a, const pq_node &b) const {
    return cmp(a.value(), b.value());
  }

  inline bool compare(const std::weak_ptr<pq_node> &a,
                      const std::weak_ptr<pq_node> &b) const {
    return compare(a.lock(), b.lock());
  }

  inline bool compare(const std::shared_ptr<pq_node> &a,
                      const std::shared_ptr<pq_node> &b) const {
    assert(a != nullptr, "this should not trigger");
    assert(b != nullptr, "this should not trigger");
    return compare(*a, *b);
  }

  const Cmp &cmp;

  std::shared_ptr<pq_node> root;
  std::size_t count = 0;
  static inline constexpr std::uint16_t PQ_HEIGHT_UPPER_BOUND = 64;
  std::array<std::shared_ptr<pq_node>, PQ_HEIGHT_UPPER_BOUND> subtrees;

  static inline gc gc;

public:
  template <typename... Args> void emplace(Args &&...args) {
    count++;
    auto node = pq_node::construct(std::forward<Args>(args)...);
    if (root == nullptr) {
      root = std::move(node);
    } else {
      if (compare(node, root)) {
        std::swap(root, node);
      }
      root->attach(std::move(node));
    }
  }

  void push(T val) { emplace(std::move(val)); }

  bool empty() const { return root == nullptr; }

  const T &top() const {
    if (empty()) {
      throw std::runtime_error("cannot pop an empty priority queue.");
    }
    return root->value();
  }

  void pop() {
    if (empty()) {
      throw std::runtime_error("cannot pop an empty priority queue.");
    }
    count--;
    root->absorb_child();
    if (root->singleton()) {
      root = nullptr;
      return;
    }
    auto new_root = root->right();
    root->detach();
    root = std::move(new_root);
    subtrees.fill(nullptr);
    consolidate();
  }

private:
  void consolidate() {
    std::shared_ptr<pq_node> new_root = root;
    auto start = root;
    while (true) {
      auto d = root->degree();
      if (d >= PQ_HEIGHT_UPPER_BOUND) {
        throw std::runtime_error("pq grew too much!");
      }
      auto u = subtrees[d];
      if (u == start) {
        break;
      } else if (u) {
        subtrees[d] = nullptr;
        u->detach();
        if (compare(root, u)) {
          root->adopt_child(u);
        } else {
          auto left = root->left();
          root->detach();
          left->attach(u);
          u->adopt_child(root);
          root = u;
        }
      } else {
        if (compare(root, new_root)) {
          new_root = root;
        }
        subtrees[root->degree()] = root;
        root = root->right();
      }
    }

    root = new_root;
  }
};
} // namespace utils
#endif
