#ifndef PQ_HPP
#define PQ_HPP

#include "core.hpp"

namespace utils {
template <typename T, typename Cmp = std::less<T>> class pq {
public:
  using value_type = T;
  using comparator_type = Cmp;

  pq() = default;
  ~pq() = default;

  pq(const pq &) = delete;
  pq(pq &&) = default;

  pq &operator=(const pq &) = delete;
  pq &operator=(pq &&) = default;

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
      auto node = std::make_shared<pq_node>(std::forward<Args>(args)...);
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
    value_type value_;
    std::size_t degree_;

  public:
    const value_type &value() const { return value_; }

    const std::shared_ptr<pq_node> &right() const { return right_; }
    void right(std::shared_ptr<pq_node> ptr) { right_ = std::move(ptr); }

    std::shared_ptr<pq_node> left() const { return left_.lock(); }
    void left(std::weak_ptr<pq_node> ptr) { left_ = std::move(ptr); }

    const std::shared_ptr<pq_node> &child() const { return child_; }
    void child(std::shared_ptr<pq_node> ptr) { child_ = std::move(ptr); }

    std::size_t degree() const { return degree_; }
  };

  struct comparator {
    inline bool operator()(const pq_node &a, const pq_node &b) const {
      return cmp(a.value(), b.value());
    }

    inline bool operator()(const std::weak_ptr<pq_node> &a,
                           const std::weak_ptr<pq_node> &b) const {
      return operator()(a.lock(), b.lock());
    }

    inline bool operator()(const std::shared_ptr<pq_node> &a,
                           const std::shared_ptr<pq_node> &b) const {
      assert(a != nullptr, "this should not trigger");
      assert(b != nullptr, "this should not trigger");
      return operator()(*a, *b);
    }

  private:
    comparator_type cmp;
  } cmp;

  std::shared_ptr<pq_node> root;
  std::size_t count = 0;
  std::unordered_map<std::size_t, std::shared_ptr<pq_node>> subtrees;

public:
  template <typename... Args> void emplace(Args &&...args) {
    count++;
    auto node = pq_node::construct(std::forward<Args>(args)...);
    if (root == nullptr) {
      root = std::move(node);
    } else {
      if (cmp(node, root)) {
        std::swap(root, node);
      }
      root->attach(std::move(node));
    }
  }

  void push(value_type val) { emplace(std::move(val)); }

  bool empty() const { return root == nullptr; }

  const value_type &top() const {
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
    subtrees.clear();
    consolidate();
  }

private:
  void consolidate() {
    std::shared_ptr<pq_node> new_root = root;
    while (true) {
      auto d = root->degree();
      auto u = subtrees[d];
      if (u == root) {
        break;
      } else if (u) {
        subtrees.erase(d);
        u->detach();
        if (cmp(root, u)) {
          root->adopt_child(u);
        } else {
          auto left = root->left();
          root->detach();
          left->attach(u);
          u->adopt_child(root);
          root = u;
        }
      } else {
        if (cmp(root, new_root)) {
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
