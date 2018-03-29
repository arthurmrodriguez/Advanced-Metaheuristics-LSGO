set (COMPILE_DARP OFF CACHE BOOL "Compile DARP Problem")

if (COMPILE_DARP)
  include_directories(${GAEDA_INCLUDE} ${LIBCONFIG_INCLUDE})

  set (SRCS darp.cc DARPGenome.cc VNSDARPGenome.cc TSDARPGenome.cc darpOps.cc aux.cc DistMatrix.cc VerticesList.cc DARPCostMatrix.cc DARPEvaluator.cc darpInit.cc skybusInit/heuristicoInsercion.cc skybusInit/ordenacionPeticiones.cc skybusInit/peticionInsercion.cc skybusInit/vehiculo.cc RandomReqDistributor.cc ReqDistMatrix.cc ParReqDistMatrix.cc KMedoidsReqDistributor.cc extras/cluster.cc RouteDistMatrix.cc RandomRouteDistributor.cc KMedoidsRouteDistributor.cc LazyReqDistMatrix.cc)

  add_library (darp SHARED ${SRCS})
  target_link_libraries (darp ${GAEDA_LIBRARY} ${LIBCONFIG})
  set_target_properties (darp PROPERTIES SUFFIX ".so")

  add_subdirectory(tests)
endif(COMPILE_DARP)
