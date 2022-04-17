LOG				?= 1
DEBUG			?= 1
PROFILER		?= 0
RELEASE			?= 1

THREADS_PER_CORE ?= 8

CXX				?= g++
CPPFLAGS		+= -std=c++17 -pedantic -Wall -Wextra -Wcast-align -Wcast-qual \
					-Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations \
					-Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow \
					-Wsign-conversion -Wsign-promo -Wstrict-overflow=5 -Wswitch-default -Wundef -DTHREADS_PER_CORE=$(THREADS_PER_CORE)

ifeq ($(LOG),1)
	CPPFLAGS += -DLOG
endif

ifeq ($(PROFILER),1)
	CPPFLAGS += -DPROFILER
endif

ifeq ($(DEBUG),1)
	CPPFLAGS += -DDEBUG
else
	CPPFLAGS += -DNDEBUG
endif

ifeq ($(RELEASE),1)
	CPPFLAGS += -O3
else
	CPPFLAGS += -Og -g
endif
