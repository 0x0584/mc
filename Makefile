PROGRAM = max-clique

LOG ?= 0
RELEASE ?= 1
THREADS_PER_CORE ?= 8

SOURCES = src/graph.cpp src/enumerator.cpp src/flavour.cpp \
 src/mc.cpp src/main.cpp

HEADERS = include/enumerator.hpp include/flavour.hpp include/graph.hpp \
 include/mc.hpp include/input.hpp include/log.hpp include/thread.hpp

OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))

CXX ?= g++

CXXFLAGS = -std=c++23 -Wformat=2 -Wpedantic -Wundef -Wall -Wextra -Iinclude
CXXFLAGS += -DTHREADS_PER_CORE=$(THREADS_PER_CORE)

TDDFLAGS = -ltbb #-L/opt/local/libexec/tbb/lib
PROFFLAGS = -lprofiler #-L/opt/local/lib

LDFLAGS ?= $(TDDFLAGS)

ifeq ($(RELEASE),1)
 CXXFLAGS += -DNDEBUG -O3 -Wdisabled-optimization -Woverloaded-virtual -Wsign-conversion -Wpessimizing-move
else
 CXXFLAGS += -DDEBUG -g3 -O2 -fno-omit-frame-pointer
 LDFLAGS += $(PROFFLAGS)
endif


ifeq ($(USE_CACHE),1)
 CXXFLAGS += -DUSE_CACHE
endif

ifeq ($(LOG),1)
 CXXFLAGS += -DLOG
endif

$(PROGRAM): $(OBJECTS) $(HEADERS)
	@echo CXX $@
	@$(CXX) -o $@ $(CXXFLAGS) $(OBJECTS) $(LDFLAGS)

%.o: %.cpp $(HEADERS)
	@echo CXX $@
	@$(CXX) $(CXXFLAGS) -c $< -o $@

info:
	@echo "LOG=$(LOG)"
	@echo "RELEASE=$(RELEASE)"
	@echo "THREADS_PER_CORE=$(THREADS_PER_CORE)"
	@echo "CXX=$(CXX)"
	@echo "CXXFLAGS=$(CXXFLAGS)"
	@echo "LDFLAGS=$(LDFLAGS)"

clean:
	@rm -f $(PROGRAM) $(OBJECTS)

re: clean $(PROGRAM)
