LOG					?= 0
DEBUG				?= 0
RELEASE				?= 1

THREADS_PER_CORE 	?= 8

CXX					?= g++
CXXFLAGS			?= -std=c++17 -Wformat=2 -pedantic -Wundef -Wall -Wextra -Wdisabled-optimization -Woverloaded-virtual -Wsign-conversion -Wpassimizing-move
CXXFLAGS 			+= -DTHREADS_PER_CORE=$(THREADS_PER_CORE)

ifeq ($(LOG),1)
	CXXFLAGS += -DLOG
endif

ifeq ($(DEBUG),1)
	CXXFLAGS += -DDEBUG -g -Og
else
	CXXFLAGS += -DNDEBUG
endif

ifeq ($(RELEASE),1)
	CXXFLAGS += -O3
endif

info:
	@echo "LOG=$(LOG)"
	@echo "DEBUG=$(DEBUG)"
	@echo "RELEASE=$(RELEASE)"
	@echo "THREADS_PER_CORE=$(THREADS_PER_CORE)"
	@echo "CXX=$(CXX)"
	@echo "CXXFLAGS=$(CXXFLAGS)"

%: %.cpp %.hpp
	@$(CXX) $(CXXFLAGS) $< -o $@
