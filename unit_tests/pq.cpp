#include "pq.hpp"

int main() {

  utils::pq pq;

  pq.emplace(std::string("4"));
  assert(pq.top().str() == "4");

  pq.emplace(std::string("2"));
  assert(pq.top().str() == "2");

  pq.emplace(std::string("1"));
  assert(pq.top().str() == "1");

  pq.emplace(std::string("3"));
  assert(pq.top().str() == "1");

  pq.pop();
  assert(pq.top().str() == "2");

  pq.pop();
  assert(pq.top().str() == "3");

  pq.emplace(std::string("6"));
  assert(pq.top().str() == "3");

  pq.emplace(std::string("5"));
  assert(pq.top().str() == "3");

  pq.emplace(std::string("7"));
  assert(pq.top().str() == "3");

  pq.pop();
  assert(pq.top().str() == "4");

  pq.emplace(std::string("8"));
  assert(pq.top().str() == "4");

  pq.pop();
  assert(pq.top().str() == "5");

  pq.pop();
  assert(pq.top().str() == "6");

  pq.emplace(std::string("9"));
  assert(pq.top().str() == "6");

  pq.pop();
  assert(pq.top().str() == "7");

  pq.pop();
  assert(pq.top().str() == "8");

  pq.pop();
  assert(pq.top().str() == "9");

  pq.pop();
  assert(pq.empty());
}
