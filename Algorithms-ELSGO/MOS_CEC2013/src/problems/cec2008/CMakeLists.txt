set (COMPILE_CEC2008 OFF CACHE BOOL "Compile CEC 2008 Benchmark")

if (COMPILE_CEC2008)
  set (PROBLEMS_SRCS schwefel_problem shifted_ackley shifted_griewank shifted_rastrigin shifted_rosembrock shifted_sphere)

  include_directories (${GAEDA_INCLUDE})

  foreach (PROBLEM ${PROBLEMS_SRCS})
    add_library (${PROBLEM} SHARED ${PROBLEM}.cc )
    target_link_libraries (${PROBLEM} ${GAEDA_LIBRARY})
    set_target_properties (${PROBLEM} PROPERTIES SUFFIX ".so")
  endforeach (PROBLEM)

  add_library (fastfractal SHARED fastfractal/DoubleDip.cc fastfractal/FastFractal.cc fastfractal/FractalFunction1D.cc fastfractal/RanQD1.cc fastfractal/RanTable.cc fastfractal/UnitFunction1D.cc fastfractal/FastFractal_dd.cc)
  target_link_libraries (fastfractal ${GAEDA_LIBRARY})
  set_target_properties (fastfractal PROPERTIES SUFFIX ".so")
endif (COMPILE_CEC2008)
