set (COMPILE_BBOB2010 OFF CACHE BOOL "Compile BBOB 2010 Benchmark")

if (COMPILE_BBOB2010)
  set (BBOB2010_EXP_SRCS  mos.cc mos_algorithm.cc benchmarksdeclare.cc dirOK-linux.cc benchmarkshelper.cc benchmarks.cc benchmarksnoisy.cc fgeneric.cc experiment.cc)
  set (BBOB2010_TIME_SRCS mos.cc mos_algorithm.cc benchmarksdeclare.cc dirOK-linux.cc benchmarkshelper.cc benchmarks.cc benchmarksnoisy.cc fgeneric.cc timing.cc)

  include_directories (${GAEDA_INCLUDE} ${LIBCONFIG_INCLUDE})

  add_executable (bbob2010_exp  ${BBOB2010_EXP_SRCS} )
  add_executable (bbob2010_time ${BBOB2010_TIME_SRCS})

  target_link_libraries (bbob2010_exp  ${GAEDA_LIBRARY} ${LIBCONFIG} ${Boost_PROGRAM_OPTIONS_LIBRARY})
  target_link_libraries (bbob2010_time ${GAEDA_LIBRARY} ${LIBCONFIG} ${Boost_PROGRAM_OPTIONS_LIBRARY})
endif (COMPILE_BBOB2010)
