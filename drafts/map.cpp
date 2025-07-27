#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <vector>

void print_vector(const std::vector<int> &v) {
  for (int i : v) {
    std::cout << i << " ";
  }
  std::cout << "\n";
}

struct vector_hash {
  std::size_t
  operator()(const std::reference_wrapper<std::vector<int>> &v_ref) const {
    std::puts("hash");
    std::size_t seed = v_ref.get().size();
    std::vector<int> tmp = v_ref.get();
    print_vector(tmp);
    std::sort(tmp.begin(), tmp.end());
    print_vector(tmp);
    for (int e : tmp) {
      seed ^= hash_func(e) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    std::cout << "seed=" << seed << "\n";
    return seed;
  }

private:
  std::hash<int> hash_func;
};

struct vector_equal {
  bool operator()(const std::reference_wrapper<std::vector<int>> &a,
                  const std::reference_wrapper<std::vector<int>> &b) const {
    std::puts("equal");
    print_vector(a.get());
    print_vector(b.get());
    std::vector<int> atmp = a.get();
    std::sort(atmp.begin(), atmp.end());
    std::vector<int> btmp = b.get();
    std::sort(btmp.begin(), btmp.end());
    return atmp == btmp;
  }
};

int main() {
  std::vector<int> v{1, 2, 3};
  std::vector<int> w{2, 3, 1};
  std::unordered_map<std::reference_wrapper<std::vector<int>>, int, vector_hash,
                     vector_equal>
      m;

  m.emplace(v, 1);
  if (m.contains(w)) {
    std::puts("OK");
  } else {
    std::puts("KO");
  }
  m[w] = 2;
  if (m[v] == 2) {
    std::puts("OK");
  } else {
    std::puts("KO");
  }
}
