enable_testing()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(TIMEOUT_SMALL 5)
  set(TIMEOUT_MEDIUM 10)
  set(TIMEOUT_LARGE 20)
  set(TIMEOUT_HUGE 30)
else()
  set(TIMEOUT_SMALL 0.1)
  set(TIMEOUT_MEDIUM 0.3)
  set(TIMEOUT_LARGE 2)
  set(TIMEOUT_HUGE 5)
endif()

file(GLOB TEST_GRAPHS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "graphs/*.mtx")

foreach(graph_file IN LISTS TEST_GRAPHS)
  get_filename_component(graph_name "${graph_file}" NAME_WE)
  execute_process(
    COMMAND wc -l "${CMAKE_CURRENT_SOURCE_DIR}/${graph_file}"
    OUTPUT_VARIABLE wc_output
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  string(REGEX MATCH "[0-9]+" line_count_str "${wc_output}")
  math(EXPR line_count "${line_count_str}")

  set(category "unknown")
  set(timeout 0)
  if(line_count LESS 10000)
    set(category "small")
    set(timeout ${TIMEOUT_SMALL})
  elseif(line_count LESS 100000)
    set(category "medium")
    set(timeout ${TIMEOUT_MEDIUM})
  elseif(line_count LESS 500000)
    set(category "large")
    set(timeout ${TIMEOUT_LARGE})
  else()
    set(category "huge")
    set(timeout ${TIMEOUT_HUGE})
  endif()

  # message(STATUS "Found ${category} test graph: ${graph_name} has ${line_count} lines")

  add_test(
    NAME "${graph_name}-heuristic"
    COMMAND $<TARGET_FILE:max-clique>
	-i "${CMAKE_CURRENT_SOURCE_DIR}/${graph_file}"
	-r 1
  )

  add_test(
    NAME "${graph_name}-hybrid"
    COMMAND $<TARGET_FILE:max-clique>
	-i "${CMAKE_CURRENT_SOURCE_DIR}/${graph_file}"
	-r 1
	-y
  )

  add_test(
    NAME "${graph_name}-exact"
    COMMAND $<TARGET_FILE:max-clique>
	-i "${CMAKE_CURRENT_SOURCE_DIR}/${graph_file}"
	-r 1
	-e
  )

  set_tests_properties("${graph_name}-exact" PROPERTIES
	LABELS "${category},exact"
	TIMEOUT "${timeout}"
  )
  set_tests_properties("${graph_name}-heuristic" PROPERTIES
	LABELS "${category},heuristic"
	TIMEOUT "${timeout}"
  )

  set_tests_properties("${graph_name}-hybrid" PROPERTIES
	LABELS "${category},hybrid"
	TIMEOUT "${timeout}"
  )
endforeach()
