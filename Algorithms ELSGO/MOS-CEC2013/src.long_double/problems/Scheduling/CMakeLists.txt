set (COMPILE_Scheduling OFF CACHE BOOL "Compile Supercomputer Scheduling Problem")

if (COMPILE_Scheduling)
  include_directories(${GAEDA_INCLUDE} ${LIBCONFIG_INCLUDE} ${CMAKE_CURRENT_SOURCE_DIR}/Simulator)
  link_directories (${CMAKE_CURRENT_BINARY_DIR}/Simulator)

  add_library (scheduling SHARED scheduling.cc)

  target_link_libraries (scheduling ${GAEDA_LIBRARY} ${LIBCONFIG} 'simulator')

  set_target_properties (scheduling PROPERTIES SUFFIX ".so")

  add_dependencies (scheduling simulator)

  add_subdirectory (Simulator)
  add_subdirectory (DataGenerator)
endif (COMPILE_Scheduling)
