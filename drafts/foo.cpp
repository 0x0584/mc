#include <algorithm>
#include <cassert>
#include <chrono>
#include <execution>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std::chrono_literals;

#define SIZE 100000000

int main() {
  std::cout << std::fixed << std::setprecision(3) << std::left;

  std::vector<std::size_t> v(SIZE);

  std::iota(v.begin(), v.end(), 0);

  assert(v.size() == SIZE);

  for (std::size_t i = 0; i < v.size(); ++i) {
    assert(v[i] == i);
  }

  std::vector<std::size_t> v2;
  v2.reserve(SIZE);
  {
    auto begin = std::chrono::high_resolution_clock::now();
    std::copy_n(std::execution::par, v.cbegin(), v.size(),
                std::back_inserter(v2));
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << " par duration: "
              << std::chrono::duration<double>(end - begin).count() << "s\n";
  }

  assert(v2.size() == SIZE);

  for (std::size_t i = 0; i < v.size(); ++i) {
    assert(v[i] == v2[i]);
  }

  v2.clear();

  {
    auto begin = std::chrono::high_resolution_clock::now();
    std::copy_n(v.cbegin(), v.size(), std::back_inserter(v2));
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << " seq duration: "
              << std::chrono::duration<double>(end - begin).count() << "s\n";
  }

  assert(v2.size() == SIZE);

  for (std::size_t i = 0; i < v.size(); ++i) {
    assert(v[i] == v2[i]);
  }

  return 0;
}
