set (COMPILE_ISDA2009 OFF CACHE BOOL "Compile ISDA 2009 Benchmark")

if (COMPILE_ISDA2009)
  set (PROBLEMS_SRCS f7 f8 f9 f10 f11)

  include_directories (${GAEDA_INCLUDE})

  foreach (PROBLEM ${PROBLEMS_SRCS})
    add_library (${PROBLEM} SHARED ${PROBLEM}.cc)
    target_link_libraries (${PROBLEM} ${GAEDA_LIBRARY})
    set_target_properties (${PROBLEM} PROPERTIES SUFFIX ".so")
  endforeach (PROBLEM)
endif (COMPILE_ISDA2009)
