#include <iostream>
#include <vector>
#include <numeric>

int main() {
  std::vector<int> v(100);
  std::iota(v.begin(), v.end(), 0);

  std::vector<int> v2(v.begin() + 10, v.end() - 2);

  for (auto i : v2) {
    std::cout << i << " ";
  }
  std::cout << "\n";
}
