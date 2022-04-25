cmake_minimum_required(VERSION 3.5)

project(maxclique LANGUAGES CXX)

include_directories(include)

file(GLOB SOURCES "src/*.cpp")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  message("-- CMAKE_BUILD_TYPE is Debug")
elseif(CMAKE_BUILD_TYPE MATCHES Release)
  message("-- CMAKE_BUILD_TYPE is Release")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX17_STANDARD_COMPILE_OPTION "-std=c++17")
set(CMAKE_CXX17_EXTENSION_COMPILE_OPTION "-std=gnu++17")

if(NOT CMAKE_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER g++)
endif(NOT CMAKE_CXX_COMPILER)

if(NOT CMAKE_C_COMPILER)
  set(CMAKE_C_COMPILER gcc)
endif(NOT CMAKE_C_COMPILER)

if(NOT CMAKE_CXX_FLAGS)
set(CMAKE_CXX_FLAGS  "-Wall -Wextra -Wpedantic")
endif(NOT CMAKE_CXX_FLAGS)

#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -stdlib=libc++")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ -lc++abi")

add_executable(${PROJECT_NAME} ${SOURCES})
