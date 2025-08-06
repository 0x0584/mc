BUILD ?= Release
RELEASE ?= 1
THREADS_PER_CORE ?= 2
LOG_LEVEL ?= info
ENABLE_PROFILING ?=OFF

SANITIZE ?=
STATIC_ANALYZER ?= OFF

JOBS ?= 4
TESTS ?= heuristic
BUILD_DIR ?= build

CXX = clang++

CXXFLAGS = -std=c++23 -march=native -mtune=native -Iinclude
CXXFLAGS += -Wformat=2 -Wall -Wextra -Wpedantic -Wundef -Wdisabled-optimization -Woverloaded-virtual -Wsign-conversion -Wpessimizing-move
CXXFLAGS += -DTHREADS_PER_CORE=$(THREADS_PER_CORE) -DLOG_LEVEL=$(LOG_LEVEL)

LDFLAGS = -ltbb -lprofiler

ifeq ($(RELEASE),1)
 CXXFLAGS += -DNDEBUG -O3
else
 CXXFLAGS += -DDEBUG -O2 -g3
endif

build:
	@cmake -B $(BUILD_DIR) \
		-DCMAKE_BUILD_TYPE=$(BUILD) \
		-DCMAKE_CXX_COMPILER=$(CXX) \
		-DTHREADS_PER_CORE=$(THREADS_PER_CORE) \
		-DLOG_LEVEL=$(LOG_LEVEL) \
		-DSANITIZE=$(SANITIZE) \
		-DSTATIC_ANALYZER=$(STATIC_ANALYZER) \
		-DENABLE_PROFILING=$(ENABLE_PROFILING)

build-release: BUILD=Release
build-release: LOG_LEVEL=info
build-release: build

build-debug: BUILD=Debug
build-debug: LOG_LEVEL=debug
build-debug: build

compile:
	@make -C $(BUILD_DIR) -j$(JOBS)

recompile:
	@make -C $(BUILD_DIR) -j$(JOBS) clean all

test: compile
	@ctest --test-dir $(BUILD_DIR) -L $(TESTS)

test-again: compile
	@ctest --test-dir $(BUILD_DIR) -L $(TESTS) --rerun-failed --output-on-failure

test-release: BUILD=Release
test-release: LOG_LEVEL=warn
test-release: build test

test-debug: BUILD=Debug
test-debug: LOG_LEVEL=debug
test-debug: build test

clean:
	@rm -rf $(BUILD_DIR)

re: clean build test-again

.PHONY: build build-release build-debug compile test test-release test-debug clean re
