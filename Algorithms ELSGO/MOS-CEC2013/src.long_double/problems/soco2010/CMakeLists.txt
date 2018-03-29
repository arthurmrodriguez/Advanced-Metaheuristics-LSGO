set (COMPILE_SOCO2010 OFF CACHE BOOL "Compile SOCO 2010 Benchmark")

if (COMPILE_SOCO2010)
  set (PROBLEMS_SRCS f12 f13 f14 f15 f16 f17 f18 f19 f20 f21)

  include_directories (${GAEDA_INCLUDE})

  foreach (PROBLEM ${PROBLEMS_SRCS})
    add_library (${PROBLEM} SHARED ${PROBLEM}.cc aux.cc)
    target_link_libraries (${PROBLEM} ${GAEDA_LIBRARY})
    set_target_properties (${PROBLEM} PROPERTIES SUFFIX ".so")
  endforeach (PROBLEM)
endif (COMPILE_SOCO2010)
