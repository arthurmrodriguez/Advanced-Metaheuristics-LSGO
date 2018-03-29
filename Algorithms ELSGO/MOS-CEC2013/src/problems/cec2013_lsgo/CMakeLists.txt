set (COMPILE_CEC2013_LSGO_CXX    OFF CACHE BOOL "Compile CEC2013 LSGO C++ Code"      )
set (COMPILE_CEC2013_LSGO_MATLAB OFF CACHE BOOL "Compile CEC2013 LSGO Matlab Wrapper")
set (COMPILE_CEC2013_LSGO_OCTAVE OFF CACHE BOOL "Compile CEC2013 LSGO Octave Wrapper")

# ############## #
# Common actions #
# ############## #

set  (     MATLAB_DATA_CEC2013_LSGO "${CMAKE_CURRENT_SOURCE_DIR}/datafiles/")
file (GLOB MATLAB_SRC_CEC2013_LSGO  "${CMAKE_CURRENT_SOURCE_DIR}/matlab/*.m")

file (MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/datafiles/")
file (COPY ${MATLAB_DATA_CEC2013_LSGO} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/datafiles/")
file (COPY ${MATLAB_SRC_CEC2013_LSGO}  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/datafiles/")

# ############## #
# C++ actions    #
# ############## #

if (COMPILE_CEC2013_LSGO_CXX)
  # We must link directly to BLAS instead of using Armadillo wrappers
  # to preven getting unresolved symbols.
  # Wrappers are deactivated in include/config.hpp file of Armadillo distribution.
  find_package(Armadillo)
  find_package(HDF5)
  find_package(BLAS)

  find_path(MATIO_INCLUDE "matio.h")
  if ("${MATIO_INCLUDE}" STREQUAL "MATIO_INCLUDE-NOTFOUND")
    MESSAGE(FATAL_ERROR "Couldn't find MATIO header files.")
  endif()

  find_library(MATIO_LIBRARY matio)
  if ("${MATIO_LIBRARY}" STREQUAL "MATIO_LIBRARY-NOTFOUND")
    MESSAGE(FATAL_ERROR "Couldn't find MATIO library.")
  endif()

  # Workaround, when FindBLAS does not work
  if ("${BLAS_LIBRARIES}" STREQUAL "FALSE")
    if (NOT "${BLAS_blas_LIBRARY}" STREQUAL "BLAS_blas_LIBRARY-NOTFOUND")
       set(BLAS_LIBRARIES ${BLAS_blas_LIBRARY})
    else (NOT "${BLAS_blas_LIBRARY}" STREQUAL "BLAS_blas_LIBRARY-NOTFOUND")
        MESSAGE(FATAL_ERROR "Couldn't find none of the BLAS libraries.")
    endif()
  endif()

  include_directories(${GAEDA_INCLUDE} ${ARMADILLO_INCLUDE_DIR} ${MATIO_INCLUDE} ${CMAKE_CURRENT_SOURCE_DIR})

  add_executable(testCEC2013LSGO testCEC2013LSGO.cc cec2013_funcs.cc load_data.cc)
  target_link_libraries(testCEC2013LSGO ${ARMADILLO_LIBRARIES} ${HDF5_LIBRARIES} ${BLAS_LIBRARIES} ${MATIO_LIBRARY})

  foreach (PROBLEM RANGE 1 15)
    add_custom_command (OUTPUT f${PROBLEM}.cc COMMAND sed s/XXX/${PROBLEM}/g ${CMAKE_CURRENT_SOURCE_DIR}/template.cc > f${PROBLEM}.cc COMMENT "Creating f${PROBLEM}.cc file" DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/template.cc)
    add_custom_target (cec2013_f${PROBLEM}.o f${PROBLEM}.cc)
    add_library (cec2013_f${PROBLEM} SHARED f${PROBLEM}.cc cec2013_funcs.cc load_data.cc)
    target_link_libraries (cec2013_f${PROBLEM} ${GAEDA_LIBRARY} ${ARMADILLO_LIBRARIES} ${HDF5_LIBRARIES} ${BLAS_LIBRARIES} ${MATIO_LIBRARY})
    set_target_properties (cec2013_f${PROBLEM} PROPERTIES SUFFIX ".so")
  endforeach (PROBLEM)
endif (COMPILE_CEC2013_LSGO_CXX)

# ############## #
# Matlab actions #
# ############## #

if (COMPILE_CEC2013_LSGO_MATLAB)
  set (MATLAB_PATH $ENV{MATLAB_PATH})

  if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set (MATLAB_LINK_PATH "${MATLAB_PATH}/runtime/maci64/")
    set (MATLAB_LINK_EXT  "dylib")
  elseif (CMAKE_SYSTEM_NAME MATCHES "Linux")
    set (MATLAB_LINK_PATH "${MATLAB_PATH}/runtime/glnxa64/")
    set (MATLAB_LINK_EXT  "so")
  else (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    message (FATAL_ERROR "-- Error: Unsupported OS: ${CMAKE_SYSTEM_NAME} or architecture: ${CMAKE_SYSTEM_PROCESSOR}")
  endif (CMAKE_SYSTEM_NAME MATCHES "Darwin")

  find_library(matlablib mwmclmcrrt HINTS ${MATLAB_LINK_PATH})

  include_directories (${GAEDA_INCLUDE} "${MATLAB_PATH}/extern/include/" ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

  # Build Objectives for Matlab
  add_custom_command (OUTPUT libCEC2013LSGO.${MATLAB_LINK_EXT} COMMAND mcc -W lib:CEC2013LSGO -T link:lib ${MATLAB_SRC_CEC2013_LSGO} -d . > CEC2013LSGO.log 2> CEC2013LSGO.log COMMAND mv CEC2013LSGO.${MATLAB_LINK_EXT} libCEC2013LSGO.${MATLAB_LINK_EXT} COMMENT "Creating libCEC2013LSGO.${MATLAB_LINK_EXT} from Matlab code" DEPENDS ${MATLAB_SRC_CEC2013_LSGO} VERBATIM)

  add_custom_target (cec2013Matlablsgowrapper DEPENDS libCEC2013LSGO.${MATLAB_LINK_EXT})

  set (TEST_MATLAB_SRCS CEC2013LSGOWrapper.cc testCEC2013LSGO_Matlab.cc)

  add_executable   (testCEC2013LSGO_Matlab ${TEST_MATLAB_SRCS})
  add_dependencies (testCEC2013LSGO_Matlab cec2013Matlablsgowrapper)

  set (WRAPPER_LIB ${CMAKE_CURRENT_BINARY_DIR}/libCEC2013LSGO.${MATLAB_LINK_EXT})

  target_link_libraries (testCEC2013LSGO_Matlab ${WRAPPER_LIB} ${matlablib})

  foreach (PROBLEM RANGE 1 15)
    add_custom_command (OUTPUT f${PROBLEM}_matlab.cc COMMAND sed s/XXX/${PROBLEM}/g ${CMAKE_CURRENT_SOURCE_DIR}/template_matlab.cc > f${PROBLEM}_matlab.cc COMMENT "Creating f${PROBLEM}_matlab.cc file" DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/template_matlab.cc libCEC2013LSGO.${MATLAB_LINK_EXT})
    add_custom_target (cec2013_matlab_f${PROBLEM}.o f${PROBLEM}_matlab.cc)
    add_library (cec2013_matlab_f${PROBLEM} SHARED f${PROBLEM}_matlab.cc CEC2013LSGOWrapper.cc)
    target_link_libraries (cec2013_matlab_f${PROBLEM} ${GAEDA_LIBRARY} ${WRAPPER_LIB} ${matlablib})
    set_target_properties (cec2013_matlab_f${PROBLEM} PROPERTIES SUFFIX ".so")
  endforeach (PROBLEM)

  # Custom clean objectives
  set (CLEAN_FILES libCEC2013LSGO.${MATLAB_LINK_EXT} CEC2013LSGO.c CEC2013LSGO.exports CEC2013LSGO.h mccExcludedFiles.log CEC2013LSGO.log readme.txt)
  set_directory_properties (PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")
endif (COMPILE_CEC2013_LSGO_MATLAB)

# ############## #
# Octave actions #
# ############## #

if (COMPILE_CEC2013_LSGO_OCTAVE)
  list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
  find_package (Octave REQUIRED)

  include_directories (${GAEDA_INCLUDE} ${OCTAVE_INCLUDE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  link_directories    (${CMAKE_CURRENT_BINARY_DIR})

  # Build Objectives for Octave
  set (TEST_OCTAVE_SRCS testCEC2013LSGO_Octave.cc)

  add_executable        (testCEC2013LSGO_Octave ${TEST_OCTAVE_SRCS})
  target_link_libraries (testCEC2013LSGO_Octave ${OCTAVE_OCTAVE_LIBRARY} ${OCTAVE_OCTINTERP_LIBRARY})

  foreach (PROBLEM RANGE 1 15)
    add_custom_command (OUTPUT f${PROBLEM}_octave.cc COMMAND sed s/XXX/${PROBLEM}/g ${CMAKE_CURRENT_SOURCE_DIR}/template_octave.cc > f${PROBLEM}_octave.cc COMMENT "Creating f${PROBLEM}_octave.cc file" DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/template_octave.cc)
    add_custom_target (cec2013_octave_f${PROBLEM}.o f${PROBLEM}_octave.cc)
    add_library (cec2013_octave_f${PROBLEM} SHARED f${PROBLEM}_octave.cc)
    target_link_libraries (cec2013_octave_f${PROBLEM} ${GAEDA_LIBRARY} ${OCTAVE_OCTAVE_LIBRARY} ${OCTAVE_CRUFT_LIBRARY} ${OCTAVE_OCTINTERP_LIBRARY})
    set_target_properties (cec2013_octave_f${PROBLEM} PROPERTIES SUFFIX ".so")
  endforeach (PROBLEM)
endif (COMPILE_CEC2013_LSGO_OCTAVE)
