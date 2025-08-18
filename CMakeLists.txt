cmake_minimum_required(VERSION 3.30)
project(max-clique CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose build type" FORCE)
endif()

message(STATUS "Build Type ${CMAKE_BUILD_TYPE}")
message(STATUS "Using C++ Standard ${CMAKE_CXX_STANDARD}")

add_executable(${PROJECT_NAME})
add_custom_command(
  TARGET ${PROJECT_NAME} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E create_symlink
  "$<TARGET_FILE:${PROJECT_NAME}>"
  "${CMAKE_SOURCE_DIR}/${PROJECT_NAME}"
  COMMENT "Creating symlink for ${PROJECT_NAME} in source directory"
)

add_subdirectory(src)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O2")
set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -Os")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

set(THREADS_PER_CORE "2" CACHE STRING "Number of threads per core")
set(LOG_LEVEL "Info" CACHE STRING "LogLevel [Off|Debug|Info|Warn|Error]")
if(NOT ENABLE_PROFILING)
  set(ENABLE_PROFILING "OFF" CACHE STRING "Profiling")
endif()
target_compile_definitions(${PROJECT_NAME} PRIVATE
  THREADS_PER_CORE=${THREADS_PER_CORE}
  LOG_LEVEL=${LOG_LEVEL}
)

if(ENABLE_PROFILING STREQUAL "ON")
  message(STATUS "Profiling is enabled")
  target_compile_definitions(${PROJECT_NAME} PRIVATE
	ENABLE_PROFILING
  )
  target_link_libraries(${PROJECT_NAME} PRIVATE
	profiler
  )
endif()

set(LOCAL_INCLUDE_DIR "/opt/local/include")
set(TBB_INCLUDE_DIR "/opt/local/libexec/tbb/include/tbb")
target_include_directories(${PROJECT_NAME} PRIVATE
  include
  ${LOCAL_INCLUDE_DIR}
  ${TBB_INCLUDE_DIR}
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
  -Wno-gnu-zero-variadic-macro-arguments
  -march=native
  -mtune=native
)

set(LOCAL_LIB_DIR "/opt/local/lib")
set(LOCAL_EXEC_DIR "/opt/local/libexec")
set(TBB_LIB_DIR "/opt/local/libexec/tbb/lib")
target_link_directories(${PROJECT_NAME} PRIVATE
  ${LOCAL_LIB_DIR}
  ${LOCAL_EXEC_DIR}
  ${TBB_LIB_DIR}
)

target_link_libraries(${PROJECT_NAME} PRIVATE
  tbb
  pthread
)

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

  if(SANITIZE STREQUAL "Address")
    set(SANITIZER_FLAGS "-fsanitize=address")
  elseif(SANITIZE STREQUAL "Thread")
    set(SANITIZER_FLAGS "-fsanitize=thread")
  elseif(SANITIZE STREQUAL "Undefined")
    set(SANITIZER_FLAGS "-fsanitize=undefined")
  elseif(SANITIZE STREQUAL "Memory")
    set(SANITIZER_FLAGS "-fsanitize=memory")
  elseif(SANITIZE STREQUAL "Leak")
    set(SANITIZER_FLAGS "-fsanitize=leak")
  elseif(SANITIZE STREQUAL "Address_Undefined")
    set(SANITIZER_FLAGS "-fsanitize=address,undefined")
  elseif(SANITIZE STREQUAL "Address_Leak")
    set(SANITIZER_FLAGS "-fsanitize=address,leak")
  else()
	set(SANITIZER_FLAGS "")
  endif()

  if(NOT SANITIZER_FLAGS STREQUAL "")
	message(STATUS "Sanitizer enabled ${SANITIZE}")
  endif()

  target_compile_options(${PROJECT_NAME} PRIVATE
	-fexperimental-library
    $<$<CONFIG:Debug>:${SANITIZER_FLAGS}>
  )

  target_link_libraries(${PROJECT_NAME} PRIVATE
    c++
    c++experimental
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
    $<$<CONFIG:Debug>:${SANITIZER_FLAGS}>
  )

  if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND STATIC_ANALYZER STREQUAL "ON")
	find_program(CLANG_TIDY_EXE NAMES "clang-tidy")
	if(CLANG_TIDY_EXE)
      message(STATUS "Clang Tidy found at ${CLANG_TIDY_EXE}")
      set(CLANG_TIDY_COMMAND "${CLANG_TIDY_EXE}" "-checks=-*,modernize-*,readability-*,bugprone-*,performance-*")
      set_target_properties(${PROJECT_NAME} PROPERTIES
		CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND}"
      )
	endif()
  endif()
endif()

message(STATUS "Log Level is set to ${LOG_LEVEL}")
message(STATUS "Threads per core are set to ${THREADS_PER_CORE}")

include(tests.cmake)
