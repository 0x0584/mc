#include <cassert>
#include <cmath>
#include <cstring>
#include <pthread.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

#define NUM_TURNS 5000
#define NUM_THREADS 16

std::mutex iomutex;

using clk = std::chrono::high_resolution_clock;
using time_point = std::chrono::time_point<clk>;
using dur_double = std::chrono::duration<double>;
using std::chrono::duration_cast;

void f_FIFO(int th_id, int num) {
#ifndef _WIN32
  sched_param sch;
  sch.sched_priority = th_id;
  if (pthread_setschedparam(pthread_self(), SCHED_FIFO, &sch)) {
    std::cout << "Failed to setschedparam: " << std::strerror(errno) << '\n';
    exit(-1);
  }
#endif

  while (num--) {
  }
}

void f_RR(int th_id, int num) {
#ifndef _WIN32
  sched_param sch;
  int policy;
  sch.sched_priority = th_id;
  if (pthread_setschedparam(pthread_self(), SCHED_FIFO, &sch)) {
    std::cout << "Failed to setschedparam: " << std::strerror(errno) << '\n';
    exit(-1);
  }
#endif

  while (num--) {
  }
}

class Timer {
public:
  Timer() : _start{clk::now()} {};

  double time_ns() {
    auto duration = clk::now() - _start;
    auto elapsed_s = duration_cast<dur_double>(duration).count();
    return elapsed_s * 1000 * 1000 * 1000;
  }

private:
  time_point _start;
};

double mean(const std::vector<double> &v) {
  if (v.size() == 0) {
    return 0; // right.
  }

  double sum = 0;
  for (auto &x : v) {
    sum += x;
  }

  return sum / v.size();
}

double std_error(const std::vector<double> &v, double meanval) {
  if (v.size() == 0) {
    return 0; // right
  }

  double s{0};
  for (auto &x : v) {
    s += (x - meanval) * (x - meanval);
  }

  return std::sqrt(s / v.size());
}

template <class F> void printtime(F f) {
  std::vector<double> timings;

  for (size_t i = 0; i <= 1000; i += 1) {
    timings.push_back(f());
  }

  double meanval = mean(timings);
  double stde = std_error(timings, meanval);

  std::cout.precision(6);
  std::cout << meanval << " ns +/- " << stde << std::endl;
  std::cout << "min: " << *min_element(timings.begin(), timings.end())
            << std::endl;
  std::cout << "max: " << *max_element(timings.begin(), timings.end())
            << std::endl;
}

double test_FIFO_threads() {
  Timer t;
  std::vector<std::thread> ths(NUM_THREADS);

  for (int i = 0; i < ths.size(); ++i) {
    ths[i] = std::thread(f_FIFO, i + 1, NUM_TURNS * (ths.size() - i));
  }

  for (auto &th : ths) {
    if (th.joinable()) {
      th.join();
    }
  }

  return t.time_ns();
}

double test_RR_threads() {
  Timer t;
  std::vector<std::thread> ths(NUM_THREADS);

  for (int i = 0; i < ths.size(); ++i) {
    ths[i] = std::thread(f_RR, i + 1, NUM_TURNS * (ths.size() - i));
  }

  for (auto &th : ths) {
    if (th.joinable()) {
      th.join();
    }
  }

  return t.time_ns();
}

double test_MIXED_threads() {
  std::vector<std::thread> ths(NUM_THREADS);
  Timer t;

  for (int i = 0; i < ths.size(); ++i) {
    ths[i] = std::thread(i % 2 == 0 ? f_FIFO : f_RR, i + 1,
                         NUM_TURNS * (ths.size() - i));
  }

  for (auto &th : ths) {
    if (th.joinable()) {
      th.join();
    }
  }

  return t.time_ns();
}

int main() {
  printtime(test_FIFO_threads);
  printtime(test_RR_threads);
  printtime(test_MIXED_threads);
}
