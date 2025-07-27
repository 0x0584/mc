#include <cstdio>
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

void bar(S &s) { std::puts("bar(S &s)"); }

void foo1(S s) {
  std::puts("foo1(S s)");
  bar(s);
}

void foo2(S &&s) {
  std::puts("foo2(S &&s)");
  bar(s);
}

int main() {

  std::puts("begin");
  std::puts("");
  
  std::puts("rvalue");
  std::puts("-================");
  foo1(S(std::vector<int>{1}));
  std::puts("-================");
  foo1(std::move(S(std::vector<int>{1})));
  std::puts("-================");
  foo2(S(std::vector<int>{1}));
  std::puts("-================");
  foo2(std::move(S(std::vector<int>{1})));
  std::puts("-================");
  std::puts("");
  
  std::puts("lvalue");
  S s(std::vector<int>{1});
  std::puts("-================");
  foo1(s);
  std::puts("-================");
  foo2(std::move(s));
  std::puts("-================");

  std::puts("");
  std::puts("end");
  return 0;
}
