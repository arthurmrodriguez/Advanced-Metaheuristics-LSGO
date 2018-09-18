#Added commands
cmake_minimum_required (VERSION 2.8)
cmake_policy(SET CMP0042 NEW)

set (COMPILE_SOCO2010 OFF CACHE BOOL "Compile SOCO 2010 Benchmark")

  #set (PROBLEMS_SRCS f12 f13 f14 f15 f16 f17 f18 f19 f20 f21)
  #Add eeg_problem only
  set (PROBLEMS_SRCS eeg_problem)		
  include_directories (${GAEDA_INCLUDE})

  foreach (PROBLEM ${PROBLEMS_SRCS})
    add_library (${PROBLEM} SHARED ${PROBLEM}.cc)
    target_link_libraries (${PROBLEM} ${GAEDA_LIBRARY})
    set_target_properties (${PROBLEM} PROPERTIES SUFFIX ".so")
  endforeach (PROBLEM)

if (COMPILE_SOCO2010)
endif (COMPILE_SOCO2010)
