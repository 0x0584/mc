#include <iostream>

struct A {
  A(int x) : x(x) { std::cout << "default(A=" << x << ")\n"; }
  A(const A &other) : x(other.x + 1) { std::cout << "copy(A=" << x << ")\n"; }
  int x;

  virtual void print() { std::cout << what() << x << "\n"; }

  virtual const char *what() { return ">>"; }
};

struct B : A {
  B(int x) : A(x) { std::cout << "default(B=" << x << ")\n"; }

  const char *what() final { return "<<"; }
};

int main() {
  A a1(1);
  a1.print();
  std::cout << "\n";

  A a2 = a1;
  a2.print();
  std::cout << "\n";

  B b1(1);
  b1.print();
  std::cout << "\n";

  B b2 = b1;
  b2.print();
  std::cout << "\n";
}
