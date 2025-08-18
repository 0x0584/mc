#include "pq.hpp"

int main() {
  utils::pq<std::string> pq;

  pq.emplace("4");
  assert(pq.top() == "4");

  pq.emplace("2");
  assert(pq.top() == "2");

  pq.emplace("1");
  assert(pq.top() == "1");

  pq.emplace("3");
  assert(pq.top() == "1");

  pq.pop();
  assert(pq.top() == "2");

  pq.pop();
  assert(pq.top() == "3");

  pq.emplace("6");
  assert(pq.top() == "3");

  pq.emplace("5");
  assert(pq.top() == "3");

  pq.emplace("7");
  assert(pq.top() == "3");

  pq.pop();
  assert(pq.top() == "4");

  pq.emplace("8");
  assert(pq.top() == "4");

  pq.pop();
  assert(pq.top() == "5");

  pq.pop();
  assert(pq.top() == "6");

  pq.emplace("9");
  assert(pq.top() == "6");

  pq.pop();
  assert(pq.top() == "7");

  pq.pop();
  assert(pq.top() == "8");

  pq.pop();
  assert(pq.top() == "9");

  pq.pop();
  assert(pq.empty());
}
