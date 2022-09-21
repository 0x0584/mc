#include <cassert>
#include <iostream>
#include <set>
#include <string>
#include <string_view>

#define SIZE 1000000

int main() {
  std::set<std::string> s;
  std::set<std::string_view> sv;

  for (auto i = 0u; i < SIZE; ++i) {
    std::string str = std::to_string(size_t(SIZE * i));
    const auto &it_str = s.emplace(std::move(str));
    sv.emplace(*it_str.first);
  }

  auto it_str = s.begin();
  auto it_str_v = sv.begin();
  while (it_str != s.end() && it_str_v != sv.end()) {
    assert(*it_str == *it_str_v);
    it_str++;
    it_str_v++;
  }

  return 0;
}
