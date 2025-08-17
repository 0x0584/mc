#ifndef PQ_HPP
#define PQ_HPP

#include "core.hpp"
#include <memory_resource>

struct Foo {
  Foo() { /*logger::debug("Foo()");*/ };
  Foo(std::string s) : s(std::move(s)) { /*logger::debug("Foo(s)");*/ }
  Foo(Foo &&o) : s(std::move(o.s)) { /*logger::debug("Foo(&&)"); */ }
  Foo(const Foo &o) : s(o.s) { /*logger::debug("Foo(const Foo &)");*/ }

  ~Foo() { /*logger::debug("~Foo()"); */ };

  Foo &operator=(Foo &&o) {
    s = std::move(o.s);
    // logger::debug("Foo::operator=(&&)");
    return *this;
  }

  Foo &operator=(const Foo &o) {
    s = o.s;
    // logger::debug("Foo::operator=(const Foo &)");
    return *this;
  }

  const std::string &str() const { return s; }

  bool operator<(const Foo &rhs) const { return s < rhs.s; }

  friend std::ostream &operator<<(std::ostream &oss, const Foo &foo) {
    return oss << foo.str();
  }

private:
  std::string s;
};

namespace utils {
// template <typename T, typename Cmp = std::less<T>>

using T = Foo;
using Cmp = std::less<T>;

class pq {
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
    pq_node() { logger::debug("pq_node()"); }

    pq_node(const pq_node &) = delete;
    pq_node(pq_node &&other)
        : parent_(std::move(other.parent_)), left_(std::move(other.left_)),
          child_(std::move(other.child_)), right_(std::move(other.right_)),
          value_(std::move(other.value_)), marked_(std::move(other.marked_)),
          degree_(std::move(other.degree_)) {
      // logger::debug("pq_node(&&)");
    }

    template <typename... Args>
    pq_node(Args &&...args)
        : value_(std::forward<Args>(args)...), marked_(false), degree_(0) {
      // logger::debug("pq_node(value)");
    }

    ~pq_node() { logger::debug("~pq_node()"); }

    pq_node &operator=(const pq_node &) = delete;

    pq_node &operator=(pq_node &&other) {
      parent_ = std::move(other.parent_);
      left_ = std::move(other.left_);
      child_ = std::move(other.child_);
      right_ = std::move(other.right_);
      value_ = std::move(other.value_);
      marked_ = std::move(other.marked_);
      degree_ = std::move(other.degree_);
      // logger::debug("pq_node::operator=(&&)");
      return *this;
    }

    template <typename... Args>
    static std::shared_ptr<pq_node> construct(Args &&...args) {
      auto node = std::make_shared<pq_node>(std::forward<Args>(args)...);
      node->right_ = node;
      node->left_ = node;
      // logger::debug("node pointer created");
      return node;
    }

    void mark() { marked_ = true; }
    void unmark() { marked_ = false; }

    void increase_degree() { ++degree_; }
    void decrease_degree() { --degree_; }

    bool singleton() const { return this == right_.get(); }

    /**
     * @brief detach the node making it a singlteon.
     *
     * @note does not affect the child.
     * @note resets the parent
     */
    bool detach() {
      if (singleton()) {
        return false;
      }
      right()->left(left_);
      left()->right(right_);
      parent_.reset();
      right_ = shared_from_this();
      left_ = right_;
      return true;
    }

    void attach(std::shared_ptr<pq_node> node) {
      right()->left(node->left());
      node->left()->right(std::move(right()));
      node->left(shared_from_this());
      right(std::move(node));
    }

    void absorb_child() {
      if (!child_) {
        return;
      }
      std::weak_ptr<pq_node> walk_weak = child_;
      do {
        auto &tmp = *walk_weak.lock();
        tmp.unmark();
        tmp.parent().reset();
        walk_weak = tmp.right();
      } while (walk_weak.lock() != child_);
      attach(std::move(child_));
    }

    void adopt_child(std::shared_ptr<pq_node> child) {
      child->detach();
      child->parent(shared_from_this());
      if (!child_) {
        child_ = std::move(child);
      } else {
        child_->attach(std::move(child));
      }
      degree_++;
    }

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

    void log(const char *str = "") {
      if constexpr (logger::current_level == logger::log_level::debug) {
        logger::debug(log_node(str, shared_from_this()));
        logger::flush();
      }
    }

  private:
    std::weak_ptr<pq_node> parent_;
    std::weak_ptr<pq_node> left_;
    std::shared_ptr<pq_node> child_;
    std::shared_ptr<pq_node> right_;

    value_type value_;
    bool marked_;
    std::size_t degree_;

  public:
    // FIXME: update the interface to weak_ptr instead
    const value_type &value() const { return value_; }

    const std::shared_ptr<pq_node> &right() const { return right_; }
    void right(std::shared_ptr<pq_node> ptr) { right_ = std::move(ptr); }

    std::shared_ptr<pq_node> left() const { return left_.lock(); }
    void left(std::weak_ptr<pq_node> ptr) { left_ = std::move(ptr); }

    std::shared_ptr<pq_node> parent() const { return parent_.lock(); }
    void parent(std::weak_ptr<pq_node> ptr) { parent_ = std::move(ptr); }

    const std::shared_ptr<pq_node> &child() const { return child_; }
    void child(std::shared_ptr<pq_node> ptr) { child_ = std::move(ptr); }

    std::size_t degree() const { return degree_; }

    bool marked() const { return marked_; }
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
    root =
        merge(pq_node::construct(std::forward<Args>(args)...), std::move(root));
    root->log("after emplace ");
    count++;
  }

  bool empty() const { return root == nullptr; }

  const value_type &top() const {
    if (empty()) {
      throw std::runtime_error("cannot pop an empty priority queue.");
    }

    logger::warn("top: ", root->value());
    logger::flush();

    return root->value();
  }

  /*
    XXX: reset child marks and parents
    XXX: extract root from the heap
    XXX: merge the roots child with the top level
    XXX: loop to add more child

    FIXME: FINISH THE IMPLMENTATION
   */
  void pop() {
    if (empty()) {
      throw std::runtime_error("cannot pop an empty priority queue.");
    }

    count--;
    logger::info("\n");
    logger::info("poping:", root->value());
    logger::flush();

    root->log("root before absord ");
    root->absorb_child();
    root->log("root after absord ");
    if (auto new_root = root->right();
        (root = root->detach() ? std::move(new_root) : nullptr)) {
      root->log("after detach ");
      consolidate();
      root->log("final root ");
    }
    logger::info("pop done");
  }

private:
  void consolidate() {
    subtrees.clear();
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

  std::shared_ptr<pq_node> merge(std::shared_ptr<pq_node> root_a,
                                 std::shared_ptr<pq_node> root_b) {
    if (root_a == nullptr) {
      return root_b;
    } else if (root_b == nullptr) {
      return root_a;
    }

    if (cmp(root_b, root_a)) {
      std::swap(root_a, root_b);
    }

    root_a->attach(std::move(root_b));

    return root_a;
  }
};
} // namespace utils
#endif
