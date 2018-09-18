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
 * @brief Genetic algorithm parameters 
 *
 * This file contains the name of all the parameters accepted by the genetic algorithms
 * implemented in GAEDAlib, as well as their default values
 */

#ifndef GAPARAMETERS_H
#define GAPARAMETERS_H

/* INCLUDES */
#include "GABayesianNetwork.h"
#include "GAGaussianNetwork.h"
#include "GAStatistics.h"

/*
 * When specifying parameters for a GAlib object, you can use the fullname (the
 * name used in parameters data files) or the short name (the name typically 
 * used on the command line).  When specifying parameters in your code you can 
 * use a string, or use the predefined macros below (kind of like using the
 * resource/class names in Motif for you Xt jocks).
 */

#define gaNnGenerations          "number_of_generations"
#define gaSNnGenerations         "ngen"
#define gaNpConvergence          "convergence_percentage"
#define gaSNpConvergence         "pconv"
#define gaNnConvergence          "generations_to_convergence"
#define gaSNnConvergence         "nconv"
#define gaNpCrossover            "crossover_probability"
#define gaSNpCrossover           "pcross"
#define gaNpMutation             "mutation_probability"
#define gaSNpMutation            "pmut"
#define gaNpopulationSize        "population_size"
#define gaSNpopulationSize       "popsize"
#define gaNnPopulations          "number_of_populations"
#define gaSNnPopulations         "npop"
#define gaNpReplacement          "replacement_percentage"
#define gaSNpReplacement         "prepl"
#define gaNnReplacement          "replacement_number"
#define gaSNnReplacement         "nrepl"
#define gaNnBestGenomes          "number_of_best"
#define gaSNnBestGenomes         "nbest"
#define gaNscoreFrequency        "score_frequency"
#define gaSNscoreFrequency       "sfreq"
#define gaNflushFrequency        "flush_frequency"
#define gaSNflushFrequency       "ffreq"
#define gaNscoreFilename         "score_filename"
#define gaSNscoreFilename        "sfile"
#define gaNselectScores          "select_scores"
#define gaSNselectScores         "sscores"
#define gaNelitism               "elitism"
#define gaSNelitism              "el"
#define gaNnOffspring            "number_of_offspring"
#define gaSNnOffspring           "noffspr"
#define gaNrecordDiversity       "record_diversity"
#define gaSNrecordDiversity      "recdiv"
#define gaNpMigration            "migration_percentage"
#define gaSNpMigration           "pmig"
#define gaNnMigration            "migration_number"
#define gaSNnMigration           "nmig"
#define gaNminimaxi              "minimaxi"
#define gaSNminimaxi             "mm"
#define gaNseed                  "seed"
#define gaSNseed                 "seed"


/* GASimpleEDA */
#define gaNselectionPercentage         "selection_percentage"
#define gaSNselectionPercentage        "sperc"
#define gaNdiscreteLearningMethod      "discrete_learning_method"
#define gaSNdiscreteLearningMethod     "dlearn"
#define gaNdiscreteEBNAScoring         "discrete_ebna_scoring"
#define gaSNdiscreteEBNAScoring        "ebnasc"
#define gaNdiscreteSimulationMethod    "discrete_simulation_method"
#define gaSNdiscreteSimulationMethod   "dsim"
#define gaNcontinuousLearningMethod    "continuous_learning_method"
#define gaSNcontinuousLearningMethod   "clearn"
#define gaNcontinuousScoreMethod       "continuous_score_method"
#define gaSNcontinuousScoreMethod      "cscor"


/**
 * @brief Namespace for default parameters
 *
 * This namespace contains the default values of the GA parameters
 */

namespace DefaultParameters {

   // GCC 4.1: explicit qualifier en 'DefaultParameters::nGenerations;
   extern int       nGenerations;
   extern double     convergencePerc;
   extern int       convergenceGens;
   extern double     mutationProb;
   extern double     crossoverProb;
   extern int       populationSize;
   extern int       nPopulations;
   extern double     replacementPerc;
   extern int       replacementNum;
   extern int       scoreSelection;
   extern int       minimaxi;
   extern GABoolean recordDiversity;
   extern int       defaultSeed;

   /* GASimpleGA */
   extern GABoolean elitism;

   /* GAIncrementalGA */
   extern int       nOffsprings;

   /* GAIslandsGA */
   extern double     migrationPerc;
   extern int       migrationNum;

   /* GASimpleEDA */
   extern double    EDASelectionPerc;

   extern GABayesianNetwork::LearningMethod   bayesianLearningMethod;
   extern GABayesianNetwork::EBNALocalScoring bayesianEBNAScoring;
   extern GABayesianNetwork::SimulationMethod bayesianSimulationMethod;
   extern GAGaussianNetwork::LearningMethod   gaussianLearningMethod;
   extern GAGaussianNetwork::ScoreMethod      gaussianScoringMethod;

}

#endif
