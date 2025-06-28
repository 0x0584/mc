LOG ?= 0
RELEASE ?= 1

THREADS_PER_CORE ?= 8

MACPORTS_PATH = /opt/local/bin/

CXX = /opt/local/bin/clang++-mp-19
LLD = ld64.lld-mp-19

LIBCPP_PATH = /opt/local/libexec/llvm-19/lib/libc++/
LIBLLVM_PATH = /opt/local/libexec/llvm-19/lib/
LIBUNWIND_PATH = /opt/local/libexec/llvm-19/lib/libunwind

LIBCPP_FLAGS = -L $(LIBCPP_PATH) -Wl,-rpath,$(LIBCPP_PATH) \
 -L $(LIBLLVM_PATH) -Wl,-rpath,$(LIBLLVM_PATH) \
 -L $(LIBUNWIND_PATH) -Wl,-rpath,$(LIBUNWIND_PATH) \
 -lc++experimental -lc++

TBB_FLAGS = -L /opt/local/libexec/tbb/lib -ltbb

LDFLAGS = $(LIBCPP_FLAGS) #$(TBB_FLAGS) 

CXXFLAGS = -pthread -fexperimental-library -std=c++17 -Iinclude -Wformat=2 -pedantic -Wundef -Wall -Wextra \
 -Wdisabled-optimization -Woverloaded-virtual -Wsign-conversion -Wpessimizing-move

ifeq ($(RELEASE),1)
 CXXFLAGS += -DNDEBUG -O3
else
 CXXFLAGS += -DDEBUG -g3
endif

ifeq ($(LOG),1)
 CXXFLAGS += -DLOG
endif

CXXFLAGS += -DTHREADS_PER_CORE=$(THREADS_PER_CORE)

#LDFLAGS ?= -ltdd -pthread
SOURCE = src/mc.cpp src/graph.cpp src/main.cpp src/mc.cpp \
 src/enumerator.cpp src/flavour.cpp
HEADER = include/mc.hpp include/enumerator.hpp include/flavour.hpp \
 include/graph.hpp include/input.hpp include/log.hpp       \
 include/thread.hpp

OBJECT = $(patsubst %.cpp,%.o,$(SOURCE))

PROGRAM = max-clique

$(PROGRAM): $(OBJECT)
	@echo CXX $@
	@$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp $(HEADER)
	@echo CXX $@
	@$(CXX) $(CXXFLAGS) -c $< -o $@

info:
	@echo "LOG=$(LOG)"
	@echo "DEBUG=$(DEBUG)"
	@echo "RELEASE=$(RELEASE)"
	@echo "THREADS_PER_CORE=$(THREADS_PER_CORE)"
	@echo "CXX=$(CXX)"
	@echo "CXXFLAGS=$(CXXFLAGS)"

clean:
	@rm -f $(PROGRAM) $(OBJECT)

re: clean $(PROGRAM)
