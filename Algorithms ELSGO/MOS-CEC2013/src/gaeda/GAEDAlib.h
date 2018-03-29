// $Header: /home/cvs/galib/ga/ga.h,v 1.3 2004/12/29 16:25:22 mwall Exp $
/*****************************************************************************
 * GAEDAlib: A C++ GA library with EDA and multiprocessor (MPI) support      *
 *                                                                           *
 * (C) 2005 Pedro Diaz (pdiaz@laurel.datsi.fi.upm.es)                        *
 *                                                                           *
 * GAEDAlib is distributed under the terms of the BSD software license       *
 *                                                                           *
 * GAEDAlib is heavily based on GAlib, a C++ GA library by Mathew Wall:      *
 * Copyright (c) 1995-1996 Massachusetts Institute of Technology (MIT)       *
 * Copyright (c) 1996-2000 Matthew Wall (author of GAlib)                    *
 *                                                                           *
 * Some portions of GAEDAlib's source code come from the GNU C++ compiler    *
 * library and therefore are covered under the terms of a different license, *
 * the GNU Public License.                                                   *
 *                                                                           *
 * You should have received a file named LICENSE along with this software.   *
 * This file contains more information about the licensing conditions of     *
 * GAEDAlib as well as the full text of each license involved.               *
 *                                                                           *
 * The file AUTHORS lists the people who have contributed (directly or       *
 * indirectly) to GAEDAlib                                                   *
 *****************************************************************************/


/**
 * @file
 * @brief Main GAEDAlib include file
 *
 * This is the main GAEDAlib include file. Unless you are using
 * very specific parts  of GAEDAlib, just include this file on your
 * project and everything should be fine
 */


#ifndef GAEDALIB_H
#define GAEDALIB_H

// Make sure that we get the configuration into each of the galib components
// that will be used.
#include "gaconfig.h"

#include "garandom.h"

// Several includes
/* #include "gaversion.h" */
/* #include "gaid.h" */
/* #include "gaparameters.h" */
/* #include "GAParameter.h" */
/* #include "GAPopulation.h" */

// These are the headers for all of the genetic algorithm classes.
#include "GASimpleGA.h"
#include "GASteadyStateGA.h"
#include "GAIncGA.h"
#include "GASimpleEDA.h"
#include "DE.h"
#include "EvolutionStrategy.h"
#include "VNS.h"
#include "VNSOpAction.h"
#include "VNSOp.h"
#include "VNSShakerResults.h"
#include "TSCordeau.h"
#include "Algorithm.h"
#include "GAGeneticAlgorithm.h"
#include "TimeConvergence.h"
#include "Recombinator.h"
#include "PopElitism.h"
#include "DEElitism.h"
#include "ClassicElitism.h"
#include "PercentElitism.h"
#include "GAPopulation.h"
#include "MOSEA.h"

// Genome factory
#include "MOSGenomeFactory.h"

// Selection Schemes
#include "GASelector.h"

// MOS Techniques
#include "MOSTechnique.h"
#include "MOSTechniqueGA.h"
#include "MOSTechniqueEDA.h"
#include "MOSTechniqueDE.h"
#include "MOSTechniqueSet.h"

#ifndef GALIB_USE_NO_TEMPLATES
#include "GADCrowdingGA.h"
#endif

// Here we include the headers for all of the various genome types.  We do
// *not* include the headers for template specializations.  This prevents
// unnecessary instantiations of template objects.

#include "genomes/GAGenome.h"
#include "genomes/GA1DBinStrGenome.h"
#include "genomes/GA2DBinStrGenome.h"
#include "genomes/GA3DBinStrGenome.h"
#include "genomes/GABin2DecGenome.h"
#include "genomes/GARealGenome.h"

#ifndef GALIB_USE_NO_TEMPLATES
#include "genomes/GA1DArrayGenome.h"
#include "genomes/GA2DArrayGenome.h"
#include "genomes/GA3DArrayGenome.h"
#include "genomes/GAListGenome.h"
#include "genomes/GATreeGenome.h"

// We do *not* include the headers for template specializations.  This prevents
// unnecessary instantiations of template objects which causes grief to some
// compilers.
//#include "genomes/GAStringGenome.h"
//#include "genomes/GARealGenome.h"
#endif

// The statistics and logger includes
#include "GAStatistics.h"
#include "logger/LogStat.h"
#include "logger/GALogger.h"
#include "logger/GAFileLogger.h"
#include "logger/GANullLogger.h"
#include "logger/AgeLogStat.h"
#include "logger/ScoreLogStat.h"
#include "logger/FitnessLogStat.h"
#include "logger/ParticipationLogStat.h"
#include "logger/QualityLogStat.h"
#include "logger/GenealogyLogStat.h"
#include "logger/OptNCompsLogStat.h"
#include "logger/OptDistLogStat.h"
#include "logger/GrefenstetteLogStat.h"
#include "logger/FitIncLogStat.h"
#include "logger/NumNewChildrenLogStat.h"
#include "logger/NativePrcntLogStat.h"
#include "logger/GenDivLogStat.h"
#include "logger/EntropyLogStat.h"
#include "logger/DistAvgPoint.h"
#include "logger/EvalsLogStat.h"
#include "logger/FSSLogStat.h"
#include "logger/ImprovementsLogStat.h"
#include "logger/RoutingAlgScoreLogStat.h"
#include "logger/RoutingAlgBestSolNonPenScore.h"
#include "logger/RoutingAlgBestSolPickupTimeLogStat.h"
#include "logger/RoutingAlgBestSolDeliverySumTimeLogStat.h"
#include "logger/RoutingAlgBestSolLoadV.h"
#include "logger/RoutingAlgBestSolRideV.h"
#include "logger/RoutingAlgBestSolTWV.h"
#include "logger/VNSSuccessfulShakerLogStat.h"
#include "logger/VNSShakerPosLogStat.h"
#include "logger/EvalRouteCallsLogStat.h"
#include "logger/RuntimeLogStat.h"
#include "logger/FreememLogStat.h"

// The island model
#include "islands/EAIslandsModel.h"
#include "islands/EAIslandsModelSync.h"
#include "islands/EAIslandsTopology.h"
#include "islands/GAIslandsTopologyRing.h"
#include "islands/GAIslandsTopologyRing2.h"
#include "islands/GAIslandsTopologyA2A.h"
#include "islands/GAIslandsTopologyMesh.h"
#include "islands/GAIslandsTopologyDisc.h"
#include "islands/GAIslandsTopologyPairs.h"
#include "islands/GAIslandsTopologyDynMedoidsBased.h"
#include "islands/GAIslandsTopologyNearestNeighbor.h"
#include "islands/GAIslandsTopologyFurthestNeighbor.h"
#include "islands/GAIslandsTopologyRandom.h"
#include "islands/GAIslandsTopologyHyperCube.h"
#include "islands/GAIslandsTopologyHybrid.h"
#include "islands/VoronoiIndInit.h"
#include "islands/RandomCentGenerator.h"
#include "islands/RandomHighSeparatedCentGenerator.h"
#include "islands/islandsutils.h"
#include "islands/CommManager.h"
#include "islands/EABestEmigrantsSelector.h"
#include "islands/EARandomEmigrantsSelector.h"
#include "islands/EABestImmigrantsSelector.h"
#include "islands/EAEmigrantsSelector.h"
#include "islands/EAImmigrantsSelector.h"

// Routing distributed algorithm
#include "routingdist/RouteDistributor.h"
#include "routingdist/DistRoutingAlg.h"

// The genealgoy class
#include "GAGenealogy.h"
#include "GAGenealogyMemory.h"
#include "GAGenealogyTracer.h"

// Extra stuff
#include "extras/centroids.h"
#include "extras/problemloader.h"
#include "extras/combinations.h"


#endif
