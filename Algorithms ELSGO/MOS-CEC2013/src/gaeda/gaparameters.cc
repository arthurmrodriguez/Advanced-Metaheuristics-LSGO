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


/* INCLUDES */
#include "gaparameters.h"

/* GLOBAL VARIABLES */

namespace DefaultParameters {
	// Here we assign the default values of the GAlib default parameters.
	int   nGenerations			= 250;
	double convergencePerc		= 0.99;
	int   convergenceGens		=  30;
	double mutationProb			= 0.01;
	double crossoverProb			= 0.90;
	int   populationSize		= 50;
	int   nPopulations			= 10;
	double replacementPerc		= 0.5;
	int   replacementNum		= 10;
	int   nOffsprings			= 2;
	double migrationPerc			= 0.1;
	int   migrationNum			= 5;
	int   scoreSelection		= GAStatistics::Maximum;
	int   minimaxi				= 1;
	GABoolean recordDiversity	= gaFalse;
	GABoolean elitism			= gaTrue;
	int   defaultSeed           = 0;

	/* GASimpleEDA */
	double EDASelectionPerc		= 1.0;
	// Related to Bayesian networks
	GABayesianNetwork::LearningMethod bayesianLearningMethod = GABayesianNetwork::UMDA;
	GABayesianNetwork::EBNALocalScoring bayesianEBNAScoring= GABayesianNetwork::BIC_SCORE;
	GABayesianNetwork::SimulationMethod bayesianSimulationMethod = GABayesianNetwork::PLS;
	// Related to Gaussian networks
	GAGaussianNetwork::LearningMethod gaussianLearningMethod= GAGaussianNetwork::UMDA;
	GAGaussianNetwork::ScoreMethod gaussianScoringMethod = GAGaussianNetwork::BGe_SCORE;
}

