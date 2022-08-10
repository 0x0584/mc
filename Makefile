LOG				?= 1
DEBUG			?= 1
PROFILER		?= 0
RELEASE			?= 1

THREADS_PER_CORE ?= 8

CXX				?= g++
CXXFLAGS		?= -std=c++17 -pedantic -Wall -Wextra -Wcast-align -Wcast-qual \
					-Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations \
					-Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow \
					-Wsign-conversion -Wsign-promo -Wstrict-overflow=5 -Wswitch-default -Wundef

CXXFLAGS +=  -DTHREADS_PER_CORE=$(THREADS_PER_CORE)

ifeq ($(LOG),1)
	CXXFLAGS += -DLOG
endif

ifeq ($(PROFILER),1)
	CXXFLAGS += -DPROFILER
endif

ifeq ($(DEBUG),1)
	CXXFLAGS += -DDEBUG
else
	CXXFLAGS += -DNDEBUG
endif

ifeq ($(RELEASE),1)
	CXXFLAGS += -O3
else
	CXXFLAGS += -Og -g
endif

%: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) $< -o $@
