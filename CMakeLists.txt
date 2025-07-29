cmake_minimum_required(VERSION 3.30)
project(max-clique CXX)

option(LOG "Enable logging" OFF)
set(THREADS_PER_CORE "2" CACHE STRING "Number of threads per core")
set(LOG_LEVEL "Info" CACHE STRING "LogLevel [Off|Debug|Info|Warn|Error]")
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

message(STATUS "Threads per core are set to ${THREADS_PER_CORE}")
message(STATUS "Log Level is set to ${LOG_LEVEL}")
message(STATUS "Using C++${CMAKE_CXX_STANDARD}")

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose build type" FORCE)
endif()

file(GLOB SOURCES
  "src/main.cpp"
  "src/graph.cpp"
  "src/enumerator.cpp"
  "src/flavour.cpp"
  "src/mc.cpp"
)
add_executable(${PROJECT_NAME} ${SOURCES})

set(LOCAL_INCLUDE_DIR "/opt/local/include")
set(LOCAL_LIB_DIR "/opt/local/lib")
set(LOCAL_EXEC_DIR "/opt/local/libexec")

set(TBB_INCLUDE_DIR "/opt/local/libexec/tbb/include/tbb")
set(TBB_LIB_DIR "/opt/local/libexec/tbb/lib")

target_include_directories(${PROJECT_NAME} PRIVATE
  include
  ${LOCAL_INCLUDE_DIR}
  ${TBB_INCLUDE_DIR}
)

target_link_directories(${PROJECT_NAME} PRIVATE
  ${LOCAL_LIB_DIR}
  ${LOCAL_EXEC_DIR}
  ${TBB_LIB_DIR}
)

target_link_libraries(${PROJECT_NAME} PRIVATE
  tbb
  pthread
  debug profiler
)

target_compile_options(${PROJECT_NAME} PRIVATE
  -Wformat=2
  -Wall
  -Wextra
  -Wpedantic
  -Wundef
  -Wdisabled-optimization
  -Woverloaded-virtual
  -Wsign-conversion
  -Wpessimizing-move
  -march=native
  -mtune=native
)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3 -O2")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -g3 -O2")
set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -Os")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

target_compile_definitions(${PROJECT_NAME} PRIVATE
  THREADS_PER_CORE=${THREADS_PER_CORE}
  LOG_LEVEL=${LOG_LEVEL}
)

if (LOG)
  target_compile_definitions(${PROJECT_NAME} PRIVATE LOG)
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(GCC14_ROOT "/opt/local/lib/gcc14/gcc/x86_64-apple-darwin19/14.2.0")
  set(GCC14_INCLUDE_DIR "/opt/local/include/gcc14/c++")

  target_include_directories(${PROJECT_NAME} PRIVATE
    ${GCC14_ROOT}/include
    ${GCC14_ROOT}/include-fixed
    ${GCC14_INCLUDE_DIR}
    ${GCC14_INCLUDE_DIR}/x86_64-apple-darwin19
    ${GCC14_INCLUDE_DIR}/backward
  )

elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  target_compile_options(${PROJECT_NAME} PRIVATE
    -fexperimental-library
  )
  target_link_libraries(${PROJECT_NAME} PRIVATE
    c++
    c++experimental
  )

  set(LLVM19_ROOT "/opt/local/libexec/llvm-19")
  set(LLVM19_LIB_DIR "${LLVM19_ROOT}/lib")
  set(LLVM19_INCLUDE_DIR "${LLVM19_LIB_DIR}/clang/19/include")
  set(LIBCXX_LIB_DIR "${LLVM19_LIB_DIR}/libc++")
  set(LIBUNWIND_LIB_DIR "${LLVM19_LIB_DIR}/libunwind")

  target_include_directories(${PROJECT_NAME} PRIVATE
    ${LLVM19_INCLUDE_DIR}
    ${LIBCXX_LIB_DIR}
    ${LLVM19_LIB_DIR}
    ${LIBUNWIND_LIB_DIR}
  )

  target_link_directories(${PROJECT_NAME} PRIVATE
    ${LIBCXX_LIB_DIR}
    ${LLVM19_LIB_DIR}
    ${LIBUNWIND_LIB_DIR}
  )

  target_link_options(${PROJECT_NAME} PRIVATE
    -Wl,-rpath,${LIBCXX_LIB_DIR}
    -Wl,-rpath,${LLVM19_LIB_DIR}
    -Wl,-rpath,${LIBUNWIND_LIB_DIR}
  )
endif()

include(tests.cmake)
