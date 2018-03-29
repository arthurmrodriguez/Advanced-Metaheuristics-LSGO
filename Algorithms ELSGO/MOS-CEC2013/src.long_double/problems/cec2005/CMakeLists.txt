set (COMPILE_CEC2005 OFF CACHE BOOL "Compile CEC 2005 Benchmark")

if (COMPILE_CEC2005)
  include_directories(${GAEDA_INCLUDE} ${LIBCONFIG_INCLUDE})
  link_directories (${CMAKE_CURRENT_SOURCE_DIR}/sprng/lib)

  set (CEC2005_SRCS aux.c benchmarks.c cec2005.cc functions.c initialize.c rand.c utils.c)
  add_library (cec2005 SHARED ${CEC2005_SRCS})

  target_link_libraries (cec2005 ${GAEDA_LIBRARY} ${LIBCONFIG} 'lcg')

  set_target_properties (cec2005 PROPERTIES SUFFIX ".so")

  add_dependencies (cec2005 sprng)

  add_subdirectory (sprng)
endif (COMPILE_CEC2005)
