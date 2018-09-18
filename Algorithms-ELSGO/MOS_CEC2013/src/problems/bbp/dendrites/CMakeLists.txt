set (COMPILE_BBP_DENDRITES OFF CACHE BOOL "Compile Cajal BBP Dendrites Grow Problem")

if (COMPILE_BBP_DENDRITES)
  find_package (Boost 1.40 COMPONENTS regex REQUIRED)

  find_library(ALGLIB libalglib.a)
  find_path   (ALGLIB_INCLUDE xblas.h)

  set (PROBLEMS_SRCS Dendrite.cc instance.cc kolmogorov.cc weka.cc)
  set (PROBLEM dendrites)

  include_directories(${GAEDA_INCLUDE} ${ALGLIB_INCLUDE} ${LIBCONFIG_INCLUDE} ${Boost_INCLUDE_DIRS})

  add_library (${PROBLEM} SHARED ${PROBLEMS_SRCS} ${PROBLEM}.cc)
  target_link_libraries (${PROBLEM} ${GAEDA_LIBRARY} ${ALGLIB} ${LIBCONFIG} ${Boost_REGEX_LIBRARY})
  set_target_properties (${PROBLEM} PROPERTIES SUFFIX ".so")

  add_executable (test test.cc ${PROBLEMS_SRCS})
  target_link_libraries (test ${GAEDA_LIBRARY} ${Boost_REGEX_LIBRARY})

  add_executable (externalVal externalVal.cc ${PROBLEMS_SRCS})
  target_link_libraries (externalVal ${GAEDA_LIBRARY} ${LIBCONFIG} ${Boost_REGEX_LIBRARY})
endif (COMPILE_BBP_DENDRITES)
