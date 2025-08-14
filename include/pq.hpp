#ifndef PQ_HPP
#define PQ_HPP

#include "core.hpp"

struct Foo {
  Foo() { logger::debug("Foo()"); };
  Foo(std::string s) : s(std::move(s)) { logger::debug("Foo(s)"); }
  Foo(Foo &&o) : s(std::move(o.s)) { logger::debug("Foo(&&)"); }
  Foo(const Foo &o) : s(o.s) { logger::debug("Foo(const Foo &)"); }

  ~Foo() { logger::debug("~Foo()"); };

  Foo &operator=(Foo &&o) {
    s = std::move(o.s);
    logger::debug("Foo::operator=(&&)");
    return *this;
  }

  Foo &operator=(const Foo &o) {
    s = o.s;
    logger::debug("Foo::operator=(const Foo &)");
    return *this;
  }

  const std::string &str() const { return s; }

  bool operator<(const Foo &rhs) const { return s < rhs.s; }

private:
  std::string s;
};

namespace mc {
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
        : parent_(std::move(other.parent_)), prev_(std::move(other.prev_)),
          child_(std::move(other.child_)), next_(std::move(other.next_)),
          value_(std::move(other.value_)), marked_(std::move(other.marked_)),
          degree_(std::move(other.degree_)) {
      logger::debug("pq_node(&&)");
    }

    template <typename... Args>
    pq_node(Args &&...args)
        : value_(std::forward<Args>(args)...), marked_(false), degree_(0) {
      logger::debug("pq_node(value)");
    }

    ~pq_node() { logger::debug("~pq_node()"); }

    pq_node &operator=(const pq_node &) = delete;

    pq_node &operator=(pq_node &&other) {
      parent_ = std::move(other.parent_);
      prev_ = std::move(other.prev_);
      child_ = std::move(other.child_);
      next_ = std::move(other.next_);
      value_ = std::move(other.value_);
      marked_ = std::move(other.marked_);
      degree_ = std::move(other.degree_);
      logger::debug("pq_node::operator=(&&)");
      return *this;
    }

    template <typename... Args>
    static std::shared_ptr<pq_node> self(Args &&...args) {
      auto node = std::make_shared<pq_node>(std::forward<Args>(args)...);
      node->next_ = node;
      node->prev_ = node;
      logger::debug("node pointer created");
      return node;
    }

    std::shared_ptr<pq_node> self() { return shared_from_this(); }

    const value_type &value() const { return value_; }

    // FIXME: update the interface to weak_ptr instead
    std::shared_ptr<pq_node> next() const { return next_; }
    void next(std::shared_ptr<pq_node> ptr) { next_ = std::move(ptr); }

    std::shared_ptr<pq_node> prev() const { return prev_.lock(); }
    void prev(std::weak_ptr<pq_node> ptr) { prev_ = std::move(ptr); }

    std::shared_ptr<pq_node> parent() const { return parent_.lock(); }
    void parent(std::weak_ptr<pq_node> ptr) { parent_ = std::move(ptr); }

    std::shared_ptr<pq_node> child() const { return child_; }

    void mark() { marked_ = true; }
    void unmark() { marked_ = false; }

  private:
    std::weak_ptr<pq_node> parent_;
    std::weak_ptr<pq_node> prev_;
    std::shared_ptr<pq_node> child_;
    std::shared_ptr<pq_node> next_;

    value_type value_;
    bool marked_;
    std::size_t degree_;
  };

  std::shared_ptr<pq_node> root;
  comparator_type cmp;
  std::size_t count = 0;

public:
  template <typename... Args> void emplace(Args &&...args) {
    root = merge(std::move(root), pq_node::self(std::forward<Args>(args)...));
    count++;

    logger::debug("node was inserted");

    std::weak_ptr<pq_node> walk = root;
    do {
      auto walk_ptr = walk.lock();
      logger::warn(walk_ptr->value().str());
      walk = walk_ptr->next();
    } while (walk.lock() != root);

    logger::debug("node logging done");
  }

  bool empty() const { return root != nullptr; }

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

    for (std::weak_ptr<pq_node> walk = root->child();
         walk.lock() != root->child(); walk = walk.lock()->next()) {
      auto walk_ptr = walk.lock();
      assert(walk_ptr != nullptr, "this should never trigger");
      walk_ptr->mark();
      walk_ptr->parent(std::weak_ptr<pq_node>());
    }
  }

private:
  std::shared_ptr<pq_node> merge(std::shared_ptr<pq_node> root_a,
                                 std::shared_ptr<pq_node> root_b) {
    if (root_a == nullptr) {
      return root_b;
    } else if (root_b == nullptr) {
      return root_a;
    }

    assert(root_a->next() != nullptr, "this should not trigger!");
    assert(root_b->prev() != nullptr, "this should not trigger!");

    if (cmp(root_b->value(), root_a->value())) {
      std::swap(root_a, root_b);
    }

    root_a->next()->prev(root_b->prev());
    root_b->prev()->next(root_a->next());

    root_a->next(root_b);
    root_b->prev(root_a);

    return root_a;
  }
};
} // namespace utils
} // namespace mc
#endif
