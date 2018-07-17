#ifndef GAEDALIB_H
#define GAEDALIB_H

// Make sure that we get the configuration into each of the galib components
// that will be used.
#include "gaconfig.h"

// These are the headers for all of the genetic algorithm classes.
#include "Algorithm.h"
#include "GAGeneticAlgorithm.h"
#include "Recombinator.h"
#include "DEElitism.h"
#include "GAPopulation.h"

// Genome factory
#include "MOSGenomeFactory.h"

// Selection Schemes
#include "GASelector.h"

// MOS Techniques
#include "MOSTechnique.h"
#include "MOSTechniqueDE.h"
#include "MOSTechniqueSet.h"

#ifndef GALIB_USE_NO_TEMPLATES
#include "genomes/GA1DArrayGenome.h"
#endif

// The statistics and logger includes
#include "GAStatistics.h"
#include "logger/LogStat.h"
#include "logger/GALogger.h"
#include "logger/GAFileLogger.h"
#include "logger/GANullLogger.h"
#include "logger/ScoreLogStat.h"
#include "logger/FitnessLogStat.h"
#include "logger/ParticipationLogStat.h"
#include "logger/QualityLogStat.h"
#include "logger/QualityFuncActiveLogStat.h"
#include "logger/FitIncLogStat.h"
#include "logger/EvalsLogStat.h"
#include "logger/ImprovementsLogStat.h"

// The island model
#include "islands/CommManager.h"

// Extra stuff
#include "extras/problemloader.h"


#endif
