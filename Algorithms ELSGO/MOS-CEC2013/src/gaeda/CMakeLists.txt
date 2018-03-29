find_package (Boost 1.40 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set (GAEDA_SRCS Algorithm.cc ClassicElitism.cc DE.cc DEElitism.cc EvolutionStrategy.cc GABayesianNetwork.cc GADCrowdingGA.cc GAEDAConfig.cc GAGaussianNetwork.cc GAGenealogy.cc GAGenealogyMemory.cc GAGenealogyTracer.cc GenerationalElitism.cc GAGeneticAlgorithm.cc GAGraphModel.cc GAIncGA.cc GAIntOps.cc GAParameter.cc GAPopulation.cc GARealOps.cc GAScaling.cc GASelector.cc GASimpleEDA.cc GASimpleGA.cc GAStatistics.cc GASteadyStateGA.cc ImprovPerDimManager.cc IncElitism.cc LSElitism.cc MOSConversion.cc MOSEA.cc MOSEA2.cc MOSGenomeFactory.cc MOSTechnique.cc MOSTechniqueCMAES.cc MOSTechniqueDE.cc MOSTechniqueDEPopSize.cc MOSTechniqueAdapDE.cc MOSTechniqueAdapSTSDE.cc MOSTechniqueEDA.cc MOSTechniqueES.cc MOSTechniqueGA.cc MOSTechniqueLS.cc MOSTechniqueMemeticGA.cc MOSTechniqueNMA.cc MOSTechniquePopLS.cc MOSTechniqueSet.cc MOSTechniqueSTSDE.cc Matrix.cc PercentElitism.cc PopElitism.cc chi.c gaerror.cc gaparameters.cc garandom.cc normal.c quicksort.c NSC.cc MOSParticipationFunc.cc MOSParticipationFunction.cc MOSQualityFunction.cc MOSEARL.cc MOSEAMultiDeme.cc MOSConfig.cc RLImprovPerDimManager.cc RoutingAlg.cc VNS.cc VNSOpAction.cc VNSOp.cc VNSShakerResults.cc TSCordeau.cc cmaes.c solis.cc solisrand.cc MOSTechniqueLSObj.cc TimeConvergence.cc)

set (EXTRAS_SRCS extras/centroids.cc extras/distaux.cc extras/problemloader.cc extras/combinations.cc)

set (GENOMES_SRCS genomes/GA1DArrayGenome.cc genomes/GA1DBinStrGenome.cc genomes/GA2DArrayGenome.cc genomes/GA2DBinStrGenome.cc genomes/GA3DArrayGenome.cc genomes/GA3DBinStrGenome.cc genomes/GAAllele.cc genomes/GABin2DecGenome.cc genomes/GABinStr.cc genomes/GAGenome.cc genomes/GAList.cc genomes/GAListBASE.cc genomes/GAListGenome.cc genomes/GARealGenome.cc genomes/GAStringGenome.cc genomes/GATree.cc genomes/GATreeBASE.cc genomes/GATreeGenome.cc genomes/MOSGenome.cc genomes/PSOGenome.cc genomes/gabincvt.cc)

set (ISLANDS_SRCS islands/CentroidsGenerator.cc islands/CommManager.cc islands/EABestEmigrantsSelector.cc islands/EABestImmigrantsSelector.cc islands/EAEmigrantsSelector.cc islands/EAImmigrantsSelector.cc islands/EAIslandsModel.cc islands/EAIslandsModelClusteredPops.cc islands/EAIslandsModelSync.cc islands/EAIslandsTopology.cc islands/EARandomEmigrantsSelector.cc islands/GAIslandsTopologyA2A.cc islands/GAIslandsTopologyDynMedoidsBased.cc islands/GAIslandsTopologyFurthestNeighbor.cc islands/GAIslandsTopologyHyperCube.cc islands/GAIslandsTopologyNearestNeighbor.cc islands/GAIslandsTopologyRandom.cc islands/GAIslandsTopologyRing.cc islands/GAIslandsTopologyRing2.cc islands/RandomCentGenerator.cc islands/RandomHighSeparatedCentGenerator.cc islands/VoronoiIndInit.cc islands/islandsutils.cc islands/neighborconds.cc)

set (ROUTEDIST_SRCS routingdist/DistRoutingAlg.cc routingdist/RouteDistributor.cc)

# TODO: Still not included in the library
# islands/GAIslandsTopologyDisc.cc islands/GAIslandsTopologyHybrid.cc islands/GAIslandsTopologyMesh.cc islands/GAIslandsTopologyNNWithFixedNeighbors.cc islands/GAIslandsTopologyPairs.cc

set (LOGGER_SRCS logger/CollectionLogStat.cc logger/FSSLogStat.cc logger/GAFileLogger.cc logger/GALogger.cc logger/GANullLogger.cc logger/GenealogyLogStat.cc logger/ImprovementsLogStat.cc logger/OptDistLogStat.cc logger/OptNCompsLogStat.cc logger/ParticipationLogStat.cc logger/QualityLogStat.cc logger/SetLogStat.cc logger/SingleLogStat.cc logger/RoutingAlgScoreLogStat.cc)

# TODO: Still not included in the library
# logger/EntropyLogStat.cc logger/GenDivLogStat.cc logger/GrefenstetteLogStat.cc logger/NativePrcntLogStat.cc

set (LIBGAEDA_SRCS ${GAEDA_SRCS} ${CONFIG_SRCS} ${EXTRAS_SRCS} ${GENOMES_SRCS} ${ISLANDS_SRCS} ${ROUTEDIST_SRCS} ${LOGGER_SRCS})

add_library (gaeda SHARED ${LIBGAEDA_SRCS})

set_target_properties (gaeda PROPERTIES SUFFIX ".so")

target_link_libraries (gaeda ${Boost_PROGRAM_OPTIONS_LIBRARY} 'dl')

set (GAEDA_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
set (GAEDA_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/libgaeda.so" PARENT_SCOPE)
