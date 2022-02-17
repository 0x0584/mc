CXX=g++-mp-10
CPPFLAGS=-std=c++17 -Wall -Wextra -Wpedantic -Wpessimizing-move # -stdlib=macports-libstdc++
#LDFLAGS=-L/opt/local/lib/libcxx
LOG = 1
DEBUG = 0
ifeq ($(shell test $(LOG) -gt 0; echo $$?),0)
	CPPFLAGS += -DLOG
endif

ifeq ($(shell test $(LOG) -gt 1; echo $$?),0)
	CPPFLAGS += -DLOG_ALGORITHM
endif

ifeq ($(DEBUG),1)
	CPPFLAGS += -g -DDEBUG
else
	CPPFLAGS += -O3
endif
