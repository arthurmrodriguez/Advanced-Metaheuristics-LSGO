set (COMPILE_NEURALSOURCES OFF CACHE BOOL "Compile Neural Sources Problem")

if (COMPILE_NEURALSOURCES)
  find_path(UBLAS_INCLUDE "boost/numeric/ublas/matrix.hpp")
  if ("${UBLAS_INCLUDE}" STREQUAL "UBLAS_INCLUDE-NOTFOUND")
    MESSAGE(FATAL_ERROR "Couldn't find Boost UBLAS header files.")
  endif()

  include_directories(${GAEDA_INCLUDE} ${LIBCONFIG_INCLUDE} ${UBLAS_INCLUDE})

  set (SRCS neuralsources neuralsourcesoriginalobjf)

  foreach(TEST ${SRCS})
    add_library (${TEST} SHARED ${TEST}.cc)
    target_link_libraries (${TEST} ${GAEDA_LIBRARY} ${LIBCONFIG})
    set_target_properties (${TEST} PROPERTIES SUFFIX ".so")
  endforeach(TEST)
endif (COMPILE_NEURALSOURCES)
