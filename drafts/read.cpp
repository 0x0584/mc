#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <execution>
#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <thread>
#include <vector>

// XXX: read until newline appears in the buff

static inline std::ostream &operator<<(std::ostream &oss,
                                       const std::vector<char> &buff) {
  oss << "[";
  for (char c : buff) {
    oss << c;
  }
  return oss << "]";
}

static inline std::ostream &operator<<(std::ostream &oss,
                                       const std::pair<int, int> &word) {
  return oss << "(" << word.first << " " << word.second << ")";
}

static inline std::ostream &
operator<<(std::ostream &oss, const std::vector<std::pair<int, int>> &words) {
  oss << "{ ";
  for (const auto &w : words) {
    oss << w << " ";
  }
  return oss << "}";
}

constexpr std::size_t PORTION_SIZE = 16; // 16B
static_assert(PORTION_SIZE != 0);

constexpr std::size_t WORD_COUNT_PER_TASK = PORTION_SIZE / 10; // just a default

int main(int argc, char *argv[]) {
  // opening the file stream with default file
  std::ifstream ifs(argc != 2 ? "foo" : argv[1]);
  if (not ifs) {
    std::cerr << "run as " << argv[0] << " FILE";
    return 1;
  }

  // queue of tasks awaiting to launch and their buffers
  std::vector< // holding task_id, and vector of words
      std::future<std::pair<std::size_t, std::vector<std::pair<int, int>>>>>
      Q;
  // buffers of parts from incomplete word parsing of value pairs in a portion
  std::vector< // holding task_id, buffers and their leftover parts
      std::tuple<std::size_t, std::vector<char>,
                 std::vector<char>::const_iterator>>
      P;

  // I) reading the data from the stream
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  // first, read PORTION_SIZE or more from the stream. resulting in a whole
  // portion which is newlines (portion) and a probable extra parts.
  //
  // also, preserve the last incomplete portion in a temporary buffer. then,
  // launch a task to process the portion. after that, clear the buffer and
  // append the extra part
  //
  // II) processing the data
  // ~~~~~~~~~~~~~~~~~~~~~~~
  //
  // from the task buffer, read a word (two integers) that are:
  //
  //   1. might have 1 or more forward spaces or tabs
  //   2. separated by 1 or more white spaces or tabs
  //   3. might have 1 or more backward spaces or tabs
  //   4. delimited by newlines
  //
  // if any of the above constraints failed, the buffer is considered
  // ill-formatted, and an error shall be raised
  //
  std::mutex task_mtx;
  // count of each task when preparing and after the task have launched
  std::size_t queued = 0, incoming = 0;
  // TODO: tasks fail if the buffer is ill-formatted
  std::atomic_bool failed = false;
  // when a get launched, it starts processing the portion parsing words
  const auto launch_task = [&P, &failed, &task_mtx,
                            &incoming](std::vector<char> portion) {
    // XXX: skip spaces and newlines forward
    const auto seek = [](auto &it, const auto &end) {
      while (it != end && (*it == ' ' || *it == '\n')) {
        ++it;
      }
      return it != end;
    };

    // XXX: read single integer number
    const auto number = [](auto &value_ref, auto &it, const auto &end) {
      if (it == end || not std::isdigit(*it)) {
        return false;
      }
      auto count = 0; // maximum integer digit count
      while (count++ < 12 && it != end && std::isdigit(*it)) {
        ++it;
      }
      long tmp = 1;
      for (auto begin = it - count; begin != it;) {
        tmp += (*begin - '0') * ((tmp << 2 << 1) + (tmp << 1));
        //tmp += (tmp * 10) * (*begin - '0');
      }
      return (value_ref = static_cast<int>(tmp)) == tmp;
    };

    return std::async( // queue the task to process the portion immediately
        std::launch::deferred, // such so that it runs while waiting for IO
        [&](std::vector<char> portion) {
          std::size_t task;
          {
            std::unique_lock _(task_mtx);
            task = ++incoming;
            std::cerr << " > #" << incoming << " " << portion << "\n";
          }

          std::vector<std::pair<int, int>>
              words; // read two integer numbers N M
          words.reserve(WORD_COUNT_PER_TASK);

          auto it = portion.cbegin();
          const auto &end = portion.cend();
          seek(it, end);
          while (not failed && it != end) {
            std::pair<int, int> word{0, 0};
            if (number(word.first, it, end) && seek(it, end) &&
                number(word.second, it, end)) {
              words.emplace_back(std::move(word));
              seek(it, end);
            } else if (it == end) {
              P.emplace_back(task, std::move(portion), it);
              break;
            } else {
              failed = true;
            }
          }

          if (failed) {
            words.clear();
            std::unique_lock _(task_mtx);
            std::cerr << " ! #" << task
                      << " failed!! portion is ill-formatted\n";
          }

          return std::make_pair(task, std::move(words));
        },
        std::move(portion));
  };

  // buffer allocated for the task
  std::vector<char> portion(PORTION_SIZE);
  std::size_t size = 0;

  /* auto begin = std::chrono::high_resolution_clock::now(); */

  std::cerr << "BEGINNING OF STREAM\n";

  // break when reading stream ends or when a any task fails
  while (not failed && not ifs.eof() && not ifs.fail()) {
    // extract the part from the portion if it is not fit a whole
    std::vector<char> part;
    // skip if the buffer is already a portion
    while (not failed && not ifs.eof() && not ifs.fail()) {
      // FIXME: stack portions in case the input was too large for a buffer
      assert(portion.size() % PORTION_SIZE + 1 <=
             // only read N-1 whole potions where is the max value
             std::numeric_limits<std::size_t>::max() % PORTION_SIZE);

      // append the newly read buffer to the current portion
      if (ifs.read(portion.data() + size, PORTION_SIZE);
          ifs.gcount() == PORTION_SIZE) { // if we have read a full buffer
        // look for the part, which is everything after the last newline
        for (auto back = portion.begin() + size + ifs.gcount(), head = back,
                  front = portion.begin() + size;
             head != front;) {
          std::cout << (head - front) << "\n";
          //  in case there was an extra part
          if (--head; *head == '\n') {
            // resize the part to fit the tail
            part.resize(back - ++head + PORTION_SIZE);
            // and only if the iterator have moved from its initial position
            if (head != back) {
              // extract the part starting from the newline (excluded) on
              std::move(head, back, part.begin());
              // then discard the part
              portion.resize(head - front);
            }
            goto sync; // we have read a portion
          }
        }
      } else if (ifs.gcount() == 0) { // if nothing has bee read
        goto done;
      } else if (ifs.gcount() < PORTION_SIZE) { // queue the last read
        goto sync;
      }

      // XXX: reorder portion resizing
      // store current part size (after the last fill) to skip search afterwards
      portion.resize((size = portion.size()) + PORTION_SIZE);
    }

  sync:
    std::cerr << "q #" << ++queued << " " << portion
              << " size=" << portion.size() << "\n";
    Q.emplace_back(launch_task(std::move(portion)));

    portion = std::move(part);
  }
  std::cerr << "END OF STREAM\n";

done:

  /* auto end = std::chrono::high_resolution_clock::now(); */

  if (not failed) {
    for (auto &task : Q) {
      if (auto result = task.get(); not result.second.empty()) {
        std::cerr << " < #" << result.first << " " << result.second << "\n";
      }
    }
    std::sort(std::execution::par_unseq, P.begin(), P.end(),
              [](const auto &l, const auto &r) {
                return std::get<0>(l) < std::get<0>(r);
              });
    for (auto &task : Q) {
      // TODO: sequentially handle the remaining portion(s)
    }
  }

  // std::cerr << "done took "
  //           << std::chrono::duration<double>(end - begin).count() << "\n";

  return 0;
}
