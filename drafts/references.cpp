#include <cstdio>
#include <functional>
#include <vector>

struct S {
  S() { std::puts("S()"); }
  ~S() { std::puts("~S()"); }

  S(std::vector<int> v) : v(std::move(v)) { std::puts("S(v)"); }
  S(const S &s) : v(s.v) { std::puts("S(const S &)"); }
  S(S &&s) : v(std::move(s.v)) { std::puts("S(S &&)"); }

  S &operator=(const S &s) {
    std::puts("S::operator=(const S &)");
    v = s.v;
    return *this;
  }

  S &operator=(S &&s) {
    v = std::move(s.v);
    std::puts("S::operator=(S &&)");
    return *this;
  }

  void print() { std::printf("%d\n", v.empty() ? -1 : v.front()); }
  
  std::vector<int> v;
};

void foo(S &s) {
  std::puts("----------- begin ---------------");
  S t = std::move(s);
  std::puts("--------------------------");
  s = std::move(t);
  std::puts("----------- end ---------------");
}

void bar(S &s, std::reference_wrapper<S> s_ref) {
  std::puts("----------- begin ---------------");
  s = s_ref.get();
  
  std::puts("----------- end ---------------");
}

int main() {
  std::puts("--------------------------");

  // std::puts("------------ p1 --------------");
  // std::pair<S, S> p1;
  std::puts("------------ p2 --------------");
  std::pair<S, S> p2{std::vector<int>{1}, std::vector<int>{2}};
  // std::puts("------------ p3 ---------------");
  // std::pair<S, S> p3{S(std::vector<int>{1}), S(std::vector<int>{2})};
  // std::puts("------------ p4 ---------------");
  // std::pair<S, S> p4(std::vector<int>{1}, std::vector<int>{2});
  // std::puts("------------ p5 ---------------");
  // std::pair<S, S> p5(S(std::vector<int>{1}), S(std::vector<int>{2}));

  // std::puts("--------------------------");
  // foo(p2.first);
  // std::puts("--------------------------");

  std::puts("------------ bar --------------");
  bar(p2.first, p2.second);
  std::puts("----------- bar ---------------");

  return 0;
}
