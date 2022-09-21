#include <iostream>
#include <numeric>
#include <vector>

struct S {
  const std::vector<int> &v;

  S(const std::vector<int> &v) : v(v) {}
  //  S(std::vector<int> v) : v(std::move(v)) {}

  void print() {
    for (auto i : v) {
      std::cout << i << " ";
    }
    std::cout << "\n";
  }
};

static std::vector<int> values(int start) {
  std::vector<int> v(0xffff);
  std::iota(v.begin(), v.end(), start);
  return v;
}

int main() {
  S s = values(-0xff);
  { s.print(); }

  return 0;
}
