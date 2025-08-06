#include <cstdio>
#include <vector>
struct vector_pair : std::pair<std::vector<int>, std::vector<int>> {
  vector_pair(const vector_pair &) = default;
  vector_pair(vector_pair &&) = default;
  vector_pair &operator=(const vector_pair &) = default;
  vector_pair &operator=(vector_pair &&) = default;
  int getFirstBack() const { return first.back(); }
};
int main() {
  vector_pair p = std::make_pair<std::vector<int>, std::vector<int>>({1}, {2});
}
