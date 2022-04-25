#include <fstream>
#include <iomanip>

#include "max_clique.hpp"

struct stream {
  std::unique_ptr<std::ifstream> file;
  std::istream *stream_ = nullptr;

  stream(int argc, const char *argv[]) {
    if (argc > 1) {
      if (file.reset(new std::ifstream(argv[1])); not file->is_open()) {
        std::string str_err = std::string(argv[1]) + " not found";
        throw std::runtime_error(str_err.c_str());
      }
      stream_ = file.get();
    } else {
      stream_ = &std::cin;
    }
  }

  std::istream &get_stream() {
    stream_->clear();
    stream_->seekg(0, std::ios_base::beg);
    return *stream_;
  }
};

using namespace mc;

int main(int argc, const char *argv[]) {
  std::cout << std::fixed << std::setprecision(3);
  std::cerr << std::fixed << std::setprecision(3);

  stream strm(argc, argv);
  {
    MAKE_TIMER(by_neighbours);
    max_clique maxclique(
        reader(order_vertices::by_neighbours, strm.get_stream()));
  }
  std::cout << std::endl;
  {
    MAKE_TIMER(by_neighbourhood);
    max_clique maxclique(
        reader(order_vertices::by_neighbourhood, strm.get_stream()));
  }
  std::cout << std::endl;
  {
    MAKE_TIMER(by_degeneracy);
    max_clique maxclique(
        reader(order_vertices::by_degeneracy, strm.get_stream()));
  }
  return 0;
}
