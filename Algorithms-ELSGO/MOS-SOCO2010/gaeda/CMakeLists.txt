project (gaeda)

find_package (MPI REQUIRED)
find_package (Boost REQUIRED COMPONENTS program_options)

include (CMakeForceCompiler)
cmake_force_cxx_compiler (mpicxx "MPI C++ Compiler")
cmake_force_c_compiler (mpicc "MPI C Compiler")

set (GAEDA_SRCS Algorithm.cc DEElitism.cc GAEDAConfig.cc GAGeneticAlgorithm.cc GAParameter.cc GAPopulation.cc GARealOps.cc GAScaling.cc GASelector.cc GAStatistics.cc LSElitism.cc MOSConversion.cc MOSEA.cc MOSGenomeFactory.cc MOSTechnique.cc MOSTechniqueDE.cc MOSTechniqueLS.cc MOSTechniqueSet.cc gaerror.cc gaparameters.cc garandom.cc quicksort.c MOSParticipationFunction.cc MOSQualityFunction.cc MOSConfig.cc)

set (CONFIG_SRCS config/grammar.c config/libconfig.c config/libconfigcpp.cc config/scanner.c)

set (EXTRAS_SRCS extras/distaux.cc extras/problemloader.cc)

set (GENOMES_SRCS genomes/GA1DArrayGenome.cc genomes/GAAllele.cc genomes/GAGenome.cc genomes/MOSGenome.cc)

set (ISLANDS_SRCS islands/CommManager.cc)

set (LOGGER_SRCS logger/CollectionLogStat.cc logger/GAFileLogger.cc logger/ImprovementsLogStat.cc logger/GALogger.cc logger/GANullLogger.cc logger/ParticipationLogStat.cc logger/QualityLogStat.cc logger/SetLogStat.cc logger/SingleLogStat.cc)

set (LIBGAEDA_SRCS ${GAEDA_SRCS} ${CONFIG_SRCS} ${EXTRAS_SRCS} ${GENOMES_SRCS} ${ISLANDS_SRCS} ${LOGGER_SRCS})

add_library (gaeda SHARED ${LIBGAEDA_SRCS})

set_target_properties (gaeda PROPERTIES SUFFIX ".so")

target_link_libraries (gaeda ${Boost_LIBRARIES})

set (GAEDA_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}" PARENT_SCOPE)
set (GAEDA_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/libgaeda.so" PARENT_SCOPE)
