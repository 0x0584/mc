CXX=clang++-13
CPPFLAGS=-std=c++17 -stdlib=macports-libstdc++ -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations -Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-overflow=5 -Wswitch-default -Wundef

 # -stdlib=macports-libstdc++
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
	CPPFLAGS += -Og -g -DDEBUG
else
	CPPFLAGS += -O3 #-DNDEBUG
endif
