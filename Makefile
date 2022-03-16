LOG				= 1
DEBUG			= 0
RELEASE			= 1

CXX				= clang++-13
CPPFLAGS		= -std=c++17 -stdlib=macports-libstdc++ -pedantic -Wall -Wextra -Wcast-align -Wcast-qual \
					-Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-declarations \
					-Wmissing-include-dirs -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow \
					-Wsign-conversion -Wsign-promo -Wstrict-overflow=5 -Wswitch-default -Wundef

ifeq ($(LOG),1)
	CPPFLAGS += -DLOG
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
