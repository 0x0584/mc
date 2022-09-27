LOG					?= 0
RELEASE				?= 0

THREADS_PER_CORE 	?= 8

CXX					?= g++
CXXFLAGS			?= -std=c++17 -Wformat=2 -pedantic -Wundef -Wall -Wextra -Wdisabled-optimization -Woverloaded-virtual -Wsign-conversion -Wpassimizing-move

ifeq ($(RELEASE),1)
	CXXFLAGS += -DNDEBUG -O3
else
	CXXFLAGS += -DDEBUG -g3
endif

ifeq ($(LOG),1)
	CXXFLAGS += -DLOG
endif

CXXFLAGS += -DTHREADS_PER_CORE=$(THREADS_PER_CORE)

SOURCE = mc.cpp
HEADER = mc.hpp

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
	rm $(PROGRAM) $(OBJECT)
