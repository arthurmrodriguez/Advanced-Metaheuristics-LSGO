set (COMPILE_CEC2011 OFF CACHE BOOL "Compile CEC 2011 Benchmark")

if (COMPILE_CEC2011)
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

  file (GLOB MATLAB_SRC_CEC2011  "${CMAKE_CURRENT_SOURCE_DIR}/matlab/*.m")
  file (GLOB MATLAB_DATA_CEC2011 "${CMAKE_CURRENT_SOURCE_DIR}/*.mat")

  foreach (F ${MATLAB_DATA_CEC2011})
    string (REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/" "" FM ${F})
    list (APPEND MATLAB_DATA_CEC2011_INST ${FM})
  endforeach()

  include_directories (${MATLAB_PATH}/extern/include/ ${CMAKE_CURRENT_BINARY_DIR} ${GAEDA_INCLUDE})
  link_directories    (${MATLAB_LINK_PATH}            ${CMAKE_CURRENT_BINARY_DIR})

  add_custom_command (OUTPUT libCEC2011.${MATLAB_LINK_EXT} COMMAND mcc -W lib:CEC2011 -T link:lib ${MATLAB_SRC_CEC2011} -d . > CEC2011.log 2> CEC2011.log COMMAND mv CEC2011.${MATLAB_LINK_EXT} libCEC2011.${MATLAB_LINK_EXT} COMMENT "Creating libCEC2011.${MATLAB_LINK_EXT} from Matlab code" DEPENDS ${MATLAB_SRC_CEC2011} VERBATIM)
  add_custom_command (OUTPUT cec2011data COMMAND cp ${MATLAB_DATA_CEC2011} ${CMAKE_CURRENT_BINARY_DIR} COMMENT "Installing data files" DEPENDS ${MATLAB_DATA_CEC2011})
  add_custom_target (cec2011 DEPENDS libCEC2011.${MATLAB_LINK_EXT} cec2011data)

  set (TEST_SRCS mCEC_Function.cc test_function.cc)

  add_executable   (test_function ${TEST_SRCS})
  add_dependencies (test_function cec2011)

  target_link_libraries (test_function CEC2011 ${matlablib})

  set (PROBLEMS_SRCS 1 2 3 4 5 6 7 8 9 10 11_5 11_10 12_6 12_13 12_15 12_40 12_140 13_1 13_2 13_3 14 15)

  include_directories (${CMAKE_CURRENT_SOURCE_DIR})

  foreach (PROBLEM ${PROBLEMS_SRCS})
    string (REPLACE "_" "" YNEW ${PROBLEM})
    add_custom_command (OUTPUT f${PROBLEM}.cc COMMAND sed s/XXX/${PROBLEM}/g ${CMAKE_CURRENT_SOURCE_DIR}/template.cc > f${PROBLEM}.cc COMMAND sed -i \"\" s/YYY/${YNEW}/g f${PROBLEM}.cc COMMENT "Creating f${PROBLEM}.cc file" DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/template.cc libCEC2011.${MATLAB_LINK_EXT})
    add_custom_target (cec2011_f${PROBLEM}.o f${PROBLEM}.cc)
    add_library (cec2011_f${PROBLEM} SHARED f${PROBLEM}.cc mCEC_Function.cc)
    target_link_libraries (cec2011_f${PROBLEM} ${GAEDA_LIBRARY} CEC2011 ${matlablib})
    set_target_properties (cec2011_f${PROBLEM} PROPERTIES SUFFIX ".so")
  endforeach (PROBLEM)

  set (CLEAN_FILES libCEC2011.${MATLAB_LINK_EXT} CEC2011.c CEC2011.exports CEC2011.h CEC2011_mcc_component_data.c CEC2011.log CEC2011.prj mccExcludedFiles.log readme.txt ${MATLAB_DATA_CEC2011_INST})
  set_directory_properties (PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")
endif (COMPILE_CEC2011)
