#include <cassert>

#include <functional>
#include <iostream>

#include <list>
#include <optional>
#include <unordered_map>
#include <vector>

#include <memory_resource>
#include <mutex>
#include <thread>

struct cache {
  inline cache(std::size_t capacity)
      : keys(&buff), values(&buff), _capacity(capacity) {
    assert(capacity > 0);
    keys.reserve(capacity);
  }

  std::optional<int> get(int k) {
    std::unique_lock<std::mutex> lock(mtx);
    std::optional<int> value;
    if (keys.contains(k)) {
      std::pmr::list<std::pair<int, int>>::iterator it = keys[k];
      cache_hit(it);
      value.emplace(it->second);
    }
    return value;
  }

  bool set(int k, int v) {
    std::unique_lock<std::mutex> lock(mtx);
    if (keys.contains(k)) {
      std::pmr::list<std::pair<int, int>>::iterator it = keys[k];
      it->second = v;
      cache_hit(it);
      return false;
    } else {
      std::pmr::list<std::pair<int, int>>::iterator it;
      if (values.size() == _capacity) {
        it = values.end();
        --it;
        keys.erase(it->first);
        it->first = k;
        it->second = v;
        cache_hit(it);
      } else {
        it = values.emplace(values.begin(), k, v);
      }
      keys.emplace(k, it);
      return true;
    }
  }

  std::size_t size() const {
    std::unique_lock<std::mutex> lock(mtx);
    return values.size();
  }

  std::size_t capacity() const { return _capacity; }

  bool full() const { return size() == _capacity; }

private:
  inline void cache_hit(std::pmr::list<std::pair<int, int>>::iterator it) {
    values.splice(values.begin(), values, it);
  }

  mutable std::mutex mtx;

  std::pmr::monotonic_buffer_resource buff;
  std::pmr::unordered_map<int, std::pmr::list<std::pair<int, int>>::iterator>
      keys;
  std::pmr::list<std::pair<int, int>> values;
  const std::size_t _capacity;
};

int main() {
  const auto SIZE = 10'000u;
  const auto POOL_SIZE = 16u;
  cache lru(SIZE);

  std::mutex print_mtx;

  std::vector<std::thread> pool(POOL_SIZE);

  for (auto i = 0u; i < POOL_SIZE; i++) {
    pool[i] = std::thread([&print_mtx, &lru]() {
      for (auto i = 0u; i < SIZE; ++i) {
        lru.set(i, i * SIZE + i);
      }
    });
  }

  for (std::thread &th : pool) {
    if (th.joinable()) {
      th.join();
    }
  }

  for (auto i = 0u; i < SIZE; ++i) {
    std::cout << "GET " << i << "=" << lru.get(i).value_or(-1) << "\n";
  }
}
