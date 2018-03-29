include_directories(${GAEDA_INCLUDE} ${LIBCONFIG_INCLUDE})
link_directories (${CMAKE_CURRENT_BINARY_DIR})

set (SIMULATOR_SRCS Classes.cc ClusterSimulator.cc EventContainer.cc EventVector.cc ClusterMachine.cc Config.cc EventQueue.cc Job.cc)

add_library (simulator ${SIMULATOR_SRCS})
set_target_properties (simulator PROPERTIES COMPILE_FLAGS "-fPIC")

add_executable (simul simulator.cc)
target_link_libraries (simul ${GAEDA_LIBRARY} ${LIBCONFIG} 'simulator')

add_executable (schedulers schedulers.cc)
target_link_libraries (schedulers ${GAEDA_LIBRARY} ${LIBCONFIG} 'simulator')

add_dependencies (simul simulator)
add_dependencies (schedulers simulator)
