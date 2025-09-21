cmake_minimum_required(VERSION 3.30)
project(max-clique CXX)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

message(STATUS "Build Type ${CMAKE_BUILD_TYPE}")
message(STATUS "Using C++ Standard ${CMAKE_CXX_STANDARD}")

add_executable(${PROJECT_NAME})

add_custom_command(
  TARGET ${PROJECT_NAME}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E create_symlink
  "$<TARGET_FILE:${PROJECT_NAME}>"
  "${CMAKE_SOURCE_DIR}/${PROJECT_NAME}-${CMAKE_BUILD_TYPE}"
  COMMENT "Creating symlink ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}-${CMAKE_BUILD_TYPE}"
)

if(CMAKE_BUILD_TYPE STREQUAL Release)
  add_custom_command(
	TARGET ${PROJECT_NAME}
	POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E create_symlink
	"$<TARGET_FILE:${PROJECT_NAME}>"
	"${CMAKE_SOURCE_DIR}/${PROJECT_NAME}"
	COMMENT "Creating symlink ${CMAKE_SOURCE_DIR}/${PROJECT_NAME}"
  )
  add_custom_command(
	TARGET ${PROJECT_NAME}
	POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E create_symlink
	"$<TARGET_FILE:${PROJECT_NAME}>"
	"${CMAKE_SOURCE_DIR}/graphs/${PROJECT_NAME}"
	COMMENT "Creating symlink ${CMAKE_SOURCE_DIR}/graphs/${PROJECT_NAME}"
  )
endif()

# add_custom_target(copy_compile_commands
#   COMMAND ${CMAKE_COMMAND} -E copy_if_different
#   "${CMAKE_BINARY_DIR}/compile_commands.json"
#   "${CMAKE_SOURCE_DIR}/compile_commands.json"
#   COMMENT "Copying compile_commands.json to the source directory"
# )
# add_dependencies(${PROJECT_NAME} copy_compile_commands)

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
  target_compile_options(${PROJECT_NAME} PRIVATE
	-fopt-info-loop
  )
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
	target_compile_options(${PROJECT_NAME} PRIVATE
      -fno-plt
	)
  endif()

  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	target_compile_options(${PROJECT_NAME} PRIVATE
      -fopt-info-vec
      -fopt-info-missed
      -fopt-info-optimized
	)
  endif()

elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  set(MP_CLANG_CONFIG "${CMAKE_SOURCE_DIR}/clang-mp.cfg")

  if(EXISTS "${MP_CLANG_CONFIG}")
    message(STATUS "Clang configuration found: ${MP_CLANG_CONFIG}")
    target_compile_options(${PROJECT_NAME} PRIVATE
      --config ${MP_CLANG_CONFIG}
    )
    target_link_options(${PROJECT_NAME} PRIVATE
      --config ${MP_CLANG_CONFIG}
    )
  endif()

  target_compile_options(${PROJECT_NAME} PRIVATE
    $<$<CONFIG:Debug>:${SANITIZER_FLAGS}>
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

  if(NOT SANITIZER_FLAGS STREQUAL "" AND NOT CMAKE_BUILD_TYPE STREQUAL Release)
    message(STATUS "Sanitizer enabled ${SANITIZE}")
  endif()

  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(${PROJECT_NAME} PRIVATE
      # -Rpass=loop-unroll
      # -Rpass-missed=loop-unroll
      # -Rpass-analysis=loop-unroll
      # -Rpass=loop-vectorize
      # -Rpass-missed=loop-vectorize
      # -Rpass-analysis=loop-vectorize
      -fslp-vectorize
      -ffp-contract=fast
      -fno-rtti
    )
  endif()

  target_link_options(${PROJECT_NAME} PRIVATE
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

# set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Generate a compilation database for IDEs" FORCE)

#set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15")


add_subdirectory(src)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O1")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O2")
set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -Os")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

set(GC_MAX_NUM_BLOCKS "32" CACHE STRING "PMR Tweaking max_num_blocks")
set(GC_MAX_POOL_BLOCK_SIZE "2048" CACHE STRING "PMR Tweaking max_pool_block_size")

if(NOT GC_ENABLE_MEMORY_TRACKING)
  set(GC_ENABLE_MEMORY_TRACKING "OFF" CACHE STRING "Memory allocation telemetry")
endif()

if(NOT CMAKE_BUILD_TYPE STREQUAL Release)
  if(GC_ENABLE_MEMORY_TRACKING STREQUAL "ON")
	target_compile_definitions(${PROJECT_NAME} PRIVATE
	  GC_ENABLE_MEMORY_TRACKING
	)
	message(STATUS "GC_ENABLE_MEMORY_TRACKING Enabled")
  endif()
endif()

if(NOT ENABLE_PROFILING)
  set(ENABLE_PROFILING "OFF" CACHE STRING "Profiling")
endif()

if(ENABLE_PROFILING STREQUAL "ON")
  message(STATUS "Profiling is enabled")
  target_compile_definitions(${PROJECT_NAME} PRIVATE
	ENABLE_PROFILING
  )
  target_link_libraries(${PROJECT_NAME} PRIVATE
	profiler
  )
endif()


set(THREADS_PER_CORE "2" CACHE STRING "Number of threads per core")
set(CHUNK_SIZE "4096" CACHE STRING "Chunk Unit Size for every buffer of the Feed")
set(FEED_SCALE_JOBS "3" CACHE STRING "Number of available buffers per job")
set(LOG_LEVEL "Info" CACHE STRING "LogLevel [Off|Debug|Info|Warn|Error]")

set(LOCAL_INCLUDE_DIR "/opt/local/include")
set(TBB_INCLUDE_DIR "/opt/local/libexec/tbb/include/tbb")

target_compile_definitions(${PROJECT_NAME} PRIVATE
  THREADS_PER_CORE=${THREADS_PER_CORE}
  LOG_LEVEL=${LOG_LEVEL}
)

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
  -Wuninitialized  
  -Wno-gnu-zero-variadic-macro-arguments
  -march=native
  -mtune=native
)

if (CMAKE_BUILD_TYPE STREQUAL "Release")
  target_compile_options(${PROJECT_NAME} PRIVATE
    # -funroll-loops
    -ffast-math
    -fstrict-aliasing
    -fvisibility=hidden
    -flto
  )
  target_link_options(${PROJECT_NAME} PRIVATE
    -fno-semantic-interposition
  )
endif()

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

message(STATUS "Log Level is set to ${LOG_LEVEL}")
message(STATUS "Threads per core are set to ${THREADS_PER_CORE}")

#include(tests.cmake)
