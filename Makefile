LOG					?= 0
RELEASE				?= 0

THREADS_PER_CORE 	?= 8

CXX					?= g++
CXXFLAGS			?= -std=c++17 -Iinclude -Wformat=2 -pedantic -Wundef -Wall -Wextra \
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

SOURCE = src/mc.cpp src/graph.cpp src/main.cpp src/mc.cpp \
         src/enumerator.cpp src/flavour.cpp
HEADER = include/mc.hpp include/enumerator.hpp include/flavour.hpp \
         include/graph.hpp include/input.hpp include/log.hpp       \
		 include/thread.hpp

OBJECT = $(patsubst %.cpp,%.o,$(SOURCE))

PROGRAM = mc

$(PROGRAM): $(OBJECT)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp $(HEADER)
	$(CXX) $(CXXFLAGS) -c $< -o $@

info:
	@echo "LOG=$(LOG)"
	@echo "DEBUG=$(DEBUG)"
	@echo "RELEASE=$(RELEASE)"
	@echo "THREADS_PER_CORE=$(THREADS_PER_CORE)"
	@echo "CXX=$(CXX)"
	@echo "CXXFLAGS=$(CXXFLAGS)"

clean:
	rm -f $(PROGRAM) $(OBJECT)
