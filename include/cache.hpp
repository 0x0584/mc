#ifndef CACHE_HPP
#define CACHE_HPP

#include <algorithm>
#include <list>
#include <memory_resource>
#include <mutex>
#include <unordered_map>

#include "logger.hpp"

namespace std {
template <typename T> struct hash<pmr::vector<T>> {
  [[nodiscard]] inline std::size_t
  operator()(const pmr::vector<T> &v_ref) const {
    std::size_t seed = v_ref.size();
    pmr::vector<T> tmp = v_ref;
    std::sort(tmp.begin(), tmp.end());
    for (T e : tmp) {
      seed ^= hash_func(e) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }

private:
  std::hash<T> hash_func;
};

template <typename T> struct hash<vector<T>> {
  [[nodiscard]] inline std::size_t operator()(const vector<T> &v_ref) const {
    std::size_t seed = v_ref.size();
    vector<T> tmp = v_ref;
    std::sort(tmp.begin(), tmp.end());
    for (T e : tmp) {
      seed ^= hash_func(e) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }

private:
  std::hash<T> hash_func;
};

template <typename T> struct equal_to<vector<T>> {
  [[nodiscard]] inline bool operator()(const vector<T> &a,
                                       const vector<T> &b) const {
    vector<T> atmp = a;
    std::sort(atmp.begin(), atmp.end());
    vector<T> btmp = b;
    std::sort(btmp.begin(), btmp.end());
    return atmp == btmp;
  }
};

template <typename T> struct equal_to<pmr::vector<T>> {
  [[nodiscard]] inline bool operator()(const pmr::vector<T> &a,
                                       const pmr::vector<T> &b) const {
    pmr::vector<T> atmp = a;
    std::sort(atmp.begin(), atmp.end());
    pmr::vector<T> btmp = b;
    std::sort(btmp.begin(), btmp.end());
    return atmp == btmp;
  }
};
} // namespace std

namespace mc {

template <typename Key, typename Value, typename Hash = std::hash<Key>,
          typename Equal = std::equal_to<Key>>
struct lru_cache {
  using key = Key;
  using value = Value;
  using entry = std::pair<key, value>;
  using cache_store = std::pmr::list<entry>;
  using key_hash = Hash;
  using key_equal = Equal;
  using key_store = std::pmr::unordered_map<std::reference_wrapper<const key>,
                                            typename cache_store::iterator,
                                            key_hash, key_equal>;
  explicit lru_cache(std::size_t capacity)
      : _capacity(capacity), keys(&buff), store(&buff) {
    Assert(capacity > 0, "cache capacity cannot be 0");
    logger::debug("Cache capacity is", capacity);
    keys.reserve(capacity);
  }

  std::optional<entry> get(const key &k) {
    std::unique_lock<std::mutex> lock(mtx);
    std::optional<entry> value;
    if (auto it = keys.find(k); it != keys.end()) {
      typename cache_store::iterator e = it->second;
      hit(e);
      value.emplace(*e);
    }
    return value;
  }

  entry set(const key &k, value v) {
    std::unique_lock<std::mutex> lock(mtx);
    typename cache_store::iterator e;
    if (auto it = keys.find(k); it != keys.end()) {
      e = it->second;
      e->second = std::move(v);
      hit(e);
    } else {
      if (store.size() == _capacity) {
        e = store.end();
        --e;
        keys.erase(e->first);
        e->first = std::move(k);
        e->second = std::move(v);
        hit(e);
      } else {
        e = store.emplace(store.begin(), std::move(k), std::move(v));
        if (store.size() == _capacity) {
          logger::debug("Cache is full!");
        }
      }
      keys.emplace(e->first, e);
    }
    return *e;
  }

  inline std::size_t size() const {
    std::unique_lock<std::mutex> lock(mtx);
    return store.size();
  }

  inline std::size_t capacity() const { return _capacity; }

  inline bool full() const { return size() == _capacity; }

private:
  inline void hit(cache_store::iterator entry_it) {
    store.splice(store.begin(), store, entry_it);
  }

  mutable std::mutex mtx;
  const std::size_t _capacity;

  std::pmr::monotonic_buffer_resource buff;
  key_store keys;
  cache_store store;
};

} // namespace mc
#endif
