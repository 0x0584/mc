#include <iostream>

struct S {
  static inline int x = 10;
};

struct T {
  int i;
  T(int i = 0) : i(i) {}
  void foo(int x = S::x) { std::cout << x << " " << i << "\n"; }
};

int main() {

  T t1(1);

  t1.foo();

  S::x = 11;

  T t2(2);

  t2.foo();
}
