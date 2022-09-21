#include <assert.h>
#include <algorithm>
#include <atomic>
#include <deque>
#include <chrono>
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

// struct local_thread {
//   void create(std::thread &&th) {
//     main_th = th;
//     wait_th = std::thread([this] {
//       main_th.join();
//       cv.notify_one();
//     });
//   }

//   bool cond;
//   std::condition_variable cv;
//   std::mutex mtx;
//   std::thread main_th, wait_th;
// };

int main() {
  std::mutex mutex;
  std::condition_variable cv;
  bool go = false;

  std::thread thread1([&](){
    std::unique_lock<std::mutex> lock(mutex);
    while (!go)
      cv.wait(lock);
    cv.notify_one();
    std::cout << "GO1!\n";
  });

  std::thread thread2([&](){
    std::unique_lock<std::mutex> lock(mutex);
      while (!go)
        cv.wait(lock);
      std::cout << "GO2!\n";
    cv.notify_one();
  });

  {
    std::unique_lock<std::mutex> lock(mutex);
    go = true;
    cv.notify_one(); // Something happened - the threads can now process it
  }

  thread1.join();
  thread2.join();

  std::size_t default_time_to_wait = 10, value = 25;
  unsigned int num_threads = std::thread::hardware_concurrency();

  while (num_threads < value) {
    std::size_t time_to_wait = default_time_to_wait;
    //    default_time_to_wait *= 1.5;

    std::vector<std::thread> threads(num_threads);
    // std::condition_variable cv;
    std::vector<bool> slot_available(num_threads, true);
    std::mutex mtx, print_mtx, value_mtx, slot_mtx;

    bool terminated = false;

    auto print = [&](std::string s, std::size_t thread_id = -1) {
      std::unique_lock<std::mutex> lock(print_mtx);
      std::cout << (thread_id == size_t(-1) ? "" : std::string((thread_id + 1), ' ')) << s << std::endl;
    };

    std::cout << "> num_threads=" << num_threads << " time_to_wait=" << time_to_wait << "\n";

    std::atomic_size_t j{0};
    std::size_t i = value;
    while (true) {
      print("> waiting for available slots");
      {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&] {
          // std::unique_lock<std::mutex> lock(slot_mtx);
          for (std::size_t i = 0; i < slot_available.size(); ++i) {
            if (slot_available[i]) {
              return true;
            }
          }
          return false;
        });
      }
      //      cv.notify_one();

      print("> new available slots");
      if (terminated) {
        break;
      }

      for (std::size_t slot = 0; slot < slot_available.size(); ++slot) {
        if (slot_available[slot]) {
          slot_available[slot] = false;

          if (threads[slot].joinable()) {
            threads[slot].join();
          }

          threads[slot] = std::thread([&](std::size_t thread_id, std::size_t time_to) {
            std::size_t pause_time = time_to;
            print(std::to_string(thread_id) + " waiting for " + std::to_string(pause_time),
                  thread_id);

            std::this_thread::sleep_for(std::chrono::milliseconds(pause_time));
            print(std::to_string(thread_id) + " starting after " + std::to_string(2 * pause_time),
                  thread_id);

            std::this_thread::sleep_for(std::chrono::milliseconds(2 * pause_time));
            print(std::to_string(thread_id) + " ready", thread_id);

            {
              std::unique_lock<std::mutex> lock(value_mtx);
              if (not terminated) {
                print(std::to_string(thread_id) + " has value " + std::to_string(i), thread_id);
                if (i == 0) {
                  terminated = true;
                } else {
                  i--;
                  j++;
                }
              }
            }

            print(std::to_string(thread_id) + " available after " + std::to_string(pause_time),
                  thread_id);
            std::this_thread::sleep_for(std::chrono::milliseconds(pause_time));

            {
              std::unique_lock<std::mutex> lock(mtx);
              // std::unique_lock<std::mutex> lock_slot(slot_mtx);
              slot_available[thread_id] = true;
            }

            print(std::to_string(thread_id) + " available ", thread_id);

            cv.notify_one();
          }, slot, time_to_wait += 50 * (slot + 1));
          //std::this_thread::sleep_for(std::chrono::milliseconds(150));
          print(std::to_string(slot) + " created", slot);
        }
      }

    }

    print("> stopped making threads\n");
    for (auto &th : threads) {
      if (th.joinable()) {
        th.join();
      }
    }

    std::cout << "\n> all threads joined i=" << i << " j=" << j <<"\n\n";
    std::this_thread::sleep_for(std::chrono::milliseconds(3000));
    assert(j == value);
    num_threads++;
  }
}
