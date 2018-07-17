include_directories (${CMAKE_CURRENT_SOURCE_DIR})
link_directories (${CMAKE_CURRENT_BINARY_DIR})

include_directories (${GAEDA_INCLUDE} ${CMAKE_CURRENT_SOURCE_DIR})
link_directories (${CMAKE_CURRENT_BINARY_DIR})

foreach (PROBLEM RANGE 1 20)
  list(APPEND SRCS "F${PROBLEM}.cpp")
endforeach(PROBLEM)
list(APPEND SRCS "Benchmarks.cpp")

add_library (cec2012 ${SRCS})
set_target_properties (cec2012 PROPERTIES COMPILE_FLAGS "-fPIC")

foreach (PROBLEM RANGE 1 20)
  add_custom_command (OUTPUT f${PROBLEM}.cc COMMAND sed s/XXX/${PROBLEM}/g ${CMAKE_CURRENT_SOURCE_DIR}/template.cc > f${PROBLEM}.cc COMMENT "Creating f${PROBLEM}.cc file" DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/template.cc)
  add_custom_target (cec2012_f${PROBLEM}.o f${PROBLEM}.cc)
  add_library (cec2012_f${PROBLEM} SHARED f${PROBLEM}.cc)
  target_link_libraries (cec2012_f${PROBLEM} ${GAEDA_LIBRARY})
  target_link_libraries (cec2012_f${PROBLEM} 'cec2012')
  set_target_properties (cec2012_f${PROBLEM} PROPERTIES SUFFIX ".so")
endforeach (PROBLEM)
