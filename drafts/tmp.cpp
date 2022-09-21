#include <chrono>
#include <cstring>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

void print(std::string str) {
  static std::mutex print_mtx;

  std::unique_lock print_lock(print_mtx);
  std::cout << str << std::endl;
}

void foo(int limit, int prio) {

  sched_param sch_params;
  sch_params.sched_priority += prio;

  if (pthread_setschedparam(pthread_self(), SCHED_RR, &sch_params)) {
    print(std::string("Failed to set thread scheduling : ") +
          std::strerror(errno));
  }

  std::this_thread::sleep_for(std::chrono::seconds(1));

  for (int l = limit; l--;) {
    // print(std::to_string(prio) + " " + std::to_string(l));
  }

  print(std::to_string(prio) + " done");
}

int main() {
  std::vector<std::thread> ths(64);
  int p = sched_get_priority_min(SCHED_RR);
  for (auto &th : ths) {
    th = std::thread(foo, 20000, p);
    if (p < sched_get_priority_max(SCHED_RR)) {
      p++;
    }
  }

  for (auto &th : ths) {
    th.join();
  }
}
