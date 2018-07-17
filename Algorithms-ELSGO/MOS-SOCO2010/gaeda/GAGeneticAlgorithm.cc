// $Header: /home/cvs/galib/ga/GABaseGA.C,v 1.4 2004/12/29 16:24:43 mwall Exp $
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
 * @brief GAGeneticAlgorithm class impl.
 *
 * Source file for the base genetic algorithm object.
 */

/* INCLUDES */
#include <stdio.h>
#include <string.h>
#include <sstream>

#include "GAGeneticAlgorithm.h"

#include "GAEDAConfig.h"
#include "garandom.h"
#include "Recombinator.h"
#include "logger/GALogger.h"

// return the configuration string that identifies this build of the library.
static const char* rcsid = GALIB_LIBRARY_IDENTIFIER;
const char* GAConfig() { return rcsid; }


// Here are a few termination functions that you can use.  Terminators return
// gaTrue if the algorithm should finish, gaFalse otherwise.

/**
 * Termination function: Terminate after a fixed number of iterations
 *
 * @param ga Genetic algorithm that this termination function will be applied to
 */

GABoolean GAGeneticAlgorithm::TerminateUponGeneration (GAGeneticAlgorithm& ga) {

   return (ga.generation () < ga.nGenerations () ? gaFalse : gaTrue);

}


/**
 * Termination function:  If we are maximizing, then terminate when the
 * convergence has exceeded the specified convergence.  If we are
 * minimizing, then terminate when the convergence has dropped below
 * the specified convergence.
 *
 * @param ga Genetic Algorithm that this termination function will be applied to
 */

GABoolean GAGeneticAlgorithm::TerminateUponConvergence (GAGeneticAlgorithm& ga){

   GABoolean val = gaFalse;

   if (GAGenome::optCriterion() == GAGenome::MINIMIZATION)
      if(ga.convergence () == 0 || ga.convergence () > ga.pConvergence ())
         val = gaFalse;
      else
         val = gaTrue;
   else
      if (ga.convergence () < ga.pConvergence ())
         val = gaFalse;
      else
         val = gaTrue;

   return val;

}


GABoolean GAGeneticAlgorithm::TerminateUponPopConvergenceOrGen (GAGeneticAlgorithm & ga) {
  GABoolean val;
  if (TerminateUponPopConvergence(ga) || TerminateUponGeneration(ga)) val = gaTrue;
  else                                                                val = gaFalse;
  return val;
}

/**
 * Termination function: Use the ratio between the minimum and maximum to
 * determine whether the population has converged.  This method will
 * not work if the values cross zero!
 *
 * Note that this is significantly different than the definition (and the
 * bug-laden implementation) that was in version of GAlib prior to 2.4.5.
 *
 * For historical purposes, here is the old definition of this method:
 *
 * Use the ratio of the population mean divided by the population maximum to
 * determine whether the population has converged.  If we are maximizing, then
 * check to see if the ratio exceeds the convergence.  If we are minimizing,
 * then check to see if the ratio has dropped below the convergence.
 *
 * @param ga Genetic Algorithm that this termination function will be applied to
 */

GABoolean GAGeneticAlgorithm::TerminateUponPopConvergence (GAGeneticAlgorithm & ga) {

   GABoolean val = gaFalse;

   if (ga.statistics ().current (GAStatistics::Maximum) == 0)
      return val;

   long double ratio = ga.statistics ().current (GAStatistics::Minimum) /
                 ga.statistics ().current (GAStatistics::Maximum);

   stringstream message;
   message << "y el minimo al comprobar la convergencia es " << ga.statistics ().current (GAStatistics::Minimum) << " y el maximo es " << ga.statistics ().current (GAStatistics::Maximum) << endl;
   GALogger::instance()->appendLogMessage("TerminateUpoPopConvergence",message.str());

   if (GAGenome::optCriterion() == GAGenome::MINIMIZATION) {
      if (ratio <= ga.pConvergence ())
         val = gaTrue;
      else
         val = gaFalse;
   }
   else {
      if (ratio >= ga.pConvergence ())
         val = gaTrue;
      else
         val = gaFalse;
   }

   return val;

}


/*
 * New Terminator. It compares the genotypes with the reference point passed with the pConvergence. Example if
 * a 0.5 is passed as a parameter, this terminator function compares the fitness between the best individual and
 * the medium from the population. If both GENOTYPES are the same, this convergence function returns true
 */
GABoolean GAGeneticAlgorithm::TerminateUponPartialPopConvergence (GAGeneticAlgorithm & ga) {

   GABoolean val = gaFalse;
   int       ref = (int) round( (ga.pop->size()-1) * ga.pConvergence() );

  // old version with fitness values
  //if ( ga.pop->individual(0).score() == ga.pop->individual(ref).score() ) val = gaTrue;
   if ( ga.pop->individual(0).compare( ga.pop->individual(ref) ) ) val = gaTrue;

   return val;
}


GABoolean GAGeneticAlgorithm::TerminateUponScoreOrGen(GAGeneticAlgorithm & ga) {
  GABoolean val = gaFalse;

  if ( ga.statistics ().current(GAStatistics::Maximum) >= ga.pConvergence() || TerminateUponGeneration(ga) ) val = gaTrue;

  return val;
}


GABoolean GAGeneticAlgorithm::TerminateUponFitnessEvaluations(GAGeneticAlgorithm & ga) {
  return (ga.statistics().indEvals() >= ga.nEvaluations()) ? gaTrue : gaFalse;
}


/**
 * Each genetic algorithm defines this member function to declare the
 * parameters that work with it. Pass a parameter list to this function
 * and this function will configure the list with the default parameter
 * list and values for the genetic algorithm class from which you called it.
 * This is a statically defined function, so invoke it using the class name
 * of the genetic algorithm whose parameters you want to use, for example,
 * GASimpleGA::registerDefaultParameters(list).
 *
 * @param p Parameter list that will filled with the default parameters
 */

GAParameterList& GAGeneticAlgorithm::registerDefaultParameters (GAParameterList& p) {

   p.add (gaNseed,         gaSNseed,         GAParameter::INT, &DefaultParameters::defaultSeed);
   p.add (gaNminimaxi,     gaSNminimaxi,     GAParameter::INT, &DefaultParameters::minimaxi);
   p.add (gaNnGenerations, gaSNnGenerations, GAParameter::INT, &DefaultParameters::nGenerations);
   p.add (gaNnConvergence, gaSNnConvergence, GAParameter::INT, &DefaultParameters::convergenceGens);

   p.add (gaNpConvergence, gaSNpConvergence, GAParameter::DOUBLE, &DefaultParameters::convergencePerc);
   p.add (gaNpCrossover,   gaSNpCrossover,   GAParameter::DOUBLE, &DefaultParameters::crossoverProb);
   p.add (gaNpMutation,    gaSNpMutation,    GAParameter::DOUBLE, &DefaultParameters::mutationProb);

   p.add (gaNpopulationSize,  gaSNpopulationSize,  GAParameter::INT, &DefaultParameters::populationSize);
   p.add (gaNnBestGenomes,    gaSNnBestGenomes,    GAParameter::INT, &gaDefNumBestGenomes);
   p.add (gaNscoreFrequency,  gaSNscoreFrequency,  GAParameter::INT, &gaDefScoreFrequency1);
   p.add (gaNflushFrequency,  gaSNflushFrequency,  GAParameter::INT, &gaDefFlushFrequency);
   p.add (gaNrecordDiversity, gaSNrecordDiversity, GAParameter::INT, &DefaultParameters::recordDiversity);

   p.add (gaNscoreFilename,   gaSNscoreFilename,   GAParameter::STRING, gaDefScoreFilename);
   p.add (gaNselectScores,    gaSNselectScores,    GAParameter::INT,    &DefaultParameters::scoreSelection);

   return p;

}


/**
 * Constructor:  When we create a GA, we stuff the parameters with the
 * basics that will be needed by most genetic algorithms - num generations,
 * p convergence, etc.
 *
 * @param g Genome used to create the initial population
 */

GAGeneticAlgorithm::GAGeneticAlgorithm (const GAGenome& g) : stats (), params () {

   pop = new GAPopulation (g, DefaultParameters::populationSize);
   pop->geneticAlgorithm (*this);

   ud = (void*) 0;
   cf = GAGeneticAlgorithm::DEFAULT_TERMINATOR;

   d_seed = DefaultParameters::defaultSeed;
   params.add (gaNseed, gaSNseed, GAParameter::INT, &d_seed);

   ngen = DefaultParameters::nGenerations;
   params.add (gaNnGenerations, gaSNnGenerations, GAParameter::INT, &ngen);

   nconv = DefaultParameters::convergenceGens;
   stats.nConvergence (nconv);
   params.add (gaNnConvergence, gaSNnConvergence, GAParameter::INT, &nconv);

   pconv = DefaultParameters::convergencePerc;
   params.add (gaNpConvergence, gaSNpConvergence, GAParameter::DOUBLE, &pconv);

   pcross = DefaultParameters::crossoverProb;
   params.add (gaNpCrossover, gaSNpCrossover, GAParameter::DOUBLE, &pcross);

   pmut = DefaultParameters::mutationProb;
   params.add (gaNpMutation, gaSNpMutation, GAParameter::DOUBLE, &pmut);

   int psize = pop->size ();
   params.add (gaNpopulationSize, gaSNpopulationSize, GAParameter::INT, &psize);

   stats.scoreFrequency (gaDefScoreFrequency1);
   params.add (gaNscoreFrequency, gaSNscoreFrequency, GAParameter::INT, &gaDefScoreFrequency1);

   stats.flushFrequency (gaDefFlushFrequency);
   params.add (gaNflushFrequency, gaSNflushFrequency, GAParameter::INT, &gaDefFlushFrequency);

   stats.recordDiversity (DefaultParameters::recordDiversity);
   params.add (gaNrecordDiversity, gaSNrecordDiversity, GAParameter::INT, &DefaultParameters::recordDiversity);

   stats.scoreFilename (gaDefScoreFilename);
   params.add (gaNscoreFilename, gaSNscoreFilename, GAParameter::STRING, gaDefScoreFilename);

   stats.selectScores (DefaultParameters::scoreSelection);
   params.add (gaNselectScores, gaSNselectScores, GAParameter::INT, &DefaultParameters::scoreSelection);

   stats.nBestGenomes (g, gaDefNumBestGenomes);
   params.add (gaNnBestGenomes, gaSNnBestGenomes, GAParameter::INT, &gaDefNumBestGenomes);

   scross = g.sexual  ();
   across = g.asexual ();

   nevals = 5000; // Safe default value

   recombinator_ = NULL;

}


/**
 * Constructor:  When we create a GA, we stuff the parameters with the
 * basics that will be needed by most genetic algorithms - num generations,
 * p convergence, etc.
 *
 * @param p Initial population of the algorithm
 */

GAGeneticAlgorithm::GAGeneticAlgorithm (const GAPopulation& p) : stats (), params () {

   pop = new GAPopulation (p);
   pop->geneticAlgorithm (*this);

   ud = (void*) 0;
   cf = GAGeneticAlgorithm::DEFAULT_TERMINATOR;

   d_seed = DefaultParameters::defaultSeed;
   params.add (gaNseed, gaSNseed, GAParameter::INT, &d_seed);

   ngen = DefaultParameters::nGenerations;
   params.add (gaNnGenerations, gaSNnGenerations, GAParameter::INT, &ngen);

   nconv = DefaultParameters::convergenceGens;
   stats.nConvergence (nconv);
   params.add (gaNnConvergence, gaSNnConvergence, GAParameter::INT, &nconv);

   pconv = DefaultParameters::convergencePerc;
   params.add (gaNpConvergence, gaSNpConvergence, GAParameter::DOUBLE, &pconv);

   pcross = DefaultParameters::crossoverProb;
   params.add (gaNpCrossover, gaSNpCrossover, GAParameter::DOUBLE, &pcross);

   pmut = DefaultParameters::mutationProb;
   params.add (gaNpMutation, gaSNpMutation, GAParameter::DOUBLE, &pmut);

   int psize = pop->size ();
   params.add (gaNpopulationSize, gaSNpopulationSize, GAParameter::INT, &psize);

   stats.scoreFrequency (gaDefScoreFrequency1);
   params.add (gaNscoreFrequency, gaSNscoreFrequency, GAParameter::INT, &gaDefScoreFrequency1);

   stats.flushFrequency (gaDefFlushFrequency);
   params.add (gaNflushFrequency, gaSNflushFrequency, GAParameter::INT, &gaDefFlushFrequency);

   stats.recordDiversity (DefaultParameters::recordDiversity);
   params.add (gaNrecordDiversity, gaSNrecordDiversity, GAParameter::INT, &DefaultParameters::recordDiversity);

   stats.scoreFilename (gaDefScoreFilename);
   params.add (gaNscoreFilename, gaSNscoreFilename, GAParameter::STRING, gaDefScoreFilename);

   stats.selectScores (DefaultParameters::scoreSelection);
   params.add (gaNselectScores, gaSNselectScores, GAParameter::INT, &DefaultParameters::scoreSelection);

   stats.nBestGenomes (p.individual (0), gaDefNumBestGenomes);
   params.add (gaNnBestGenomes, gaSNnBestGenomes, GAParameter::INT, &gaDefNumBestGenomes);

   scross = p.individual (0).sexual  ();
   across = p.individual (0).asexual ();

   nevals = 5000; // Safe default value

   recombinator_ = NULL;

}


/**
 * Constructor:  When we create a GA, we stuff the parameters with the
 * basics that will be needed by most genetic algorithms - num generations,
 * p convergence, etc.
 *
 * @param ga Genetic algorithm from which the new population will be cloned
 */

GAGeneticAlgorithm::GAGeneticAlgorithm (const GAGeneticAlgorithm& ga) : stats (ga.stats), params (ga.params) {

   pop = ga.pop->clone ();
   pop->geneticAlgorithm (*this);

   cf = ga.cf;
   ud = ga.ud;

   ngen  = ga.ngen;
   nconv = ga.nconv;
   pconv = ga.pconv;

   pcross = ga.pcross;
   pmut   = ga.pmut;

   scross = ga.scross;
   across = ga.across;

   d_seed = ga.d_seed;

   nevals = ga.nevals;

   recombinator_ = ga.recombinator_;

}


/**
 * Generic destructor: Destroys the algorithm population
 */

GAGeneticAlgorithm::~GAGeneticAlgorithm () {

   delete pop;
   if (recombinator_)
     delete recombinator_;

}


/**
 * Copy function: Copies the population,stats and parameters from a given
 * genetic algorithm to the current population
 *
 * @param ga Genetic algorithm from which the data will be copied
 */

void GAGeneticAlgorithm::copy (const GAGeneticAlgorithm& ga) {

   if (pop)
      pop->copy (*(ga.pop));
   else
      pop = ga.pop->clone ();

   pop->geneticAlgorithm (*this);

   stats  = ga.stats;
   params = ga.params;

   cf = ga.cf;
   ud = ga.ud;

   ngen  = ga.ngen;
   nconv = ga.nconv;
   pconv = ga.pconv;

   pcross = ga.pcross;
   pmut   = ga.pmut;

   scross = ga.scross;
   across = ga.across;

   d_seed = ga.d_seed;

   nevals = ga.nevals;

   recombinator_ = ga.recombinator_;

}


/**
 * Returns the algorithm's parameters
 *
 * @return Algorithm's parameters
 */

const GAParameterList& GAGeneticAlgorithm::getParameters () const {

   return params;

}


/**
 * Sets the parameters of the algorithm given a parameter list
 *
 * @param list Parameter list to set the parameters from
 * @return The algorithm's parameter list after if has been updated
 */

const GAParameterList& GAGeneticAlgorithm::setParameters (const GAParameterList& list) {

   for (int i = 0; i < list.size (); i++)
      setptr (list [i].fullname (), list [i].value ());

   return params;

}


/**
 * Sets the parameters of the algorithm, given C/C++ argv/argc
 *
 * @param argc Number of elements in argv
 * @param argv String with command line arguments
 * @param flag Used by GAParameter::parse()
 * @return The algorithm's parameter list after if has been updated
 */

const GAParameterList& GAGeneticAlgorithm::setParameters (int& argc, char** argv, GABoolean flag) {

   params.parse (argc, argv, flag); // get the args we understand

   for (int i = 0; i < params.size (); i++)
      setptr (params [i].fullname (), params [i].value ());

   return params;

}


#ifdef GALIB_USE_STREAMS

/**
 * Sets the parameters of the algorithm from a file
 *
 * @param filename Path to the configuration file
 * @param flag Used by GAParameter::read()
 * @return The algorithm's parameter list after if has been updated
 */

const GAParameterList& GAGeneticAlgorithm::setParameters (const char* filename, GABoolean flag) {

   params.read (filename, flag);

   for (int i = 0; i < params.size (); i++)
      setptr (params [i].fullname (), params [i].value ());

   return params;

}


/**
 * Sets the parameters of the algorithm from a istream
 *
 * @param is Stream used to read the parameters
 * @param flag Used by GAParameter::read()
 * @return The algorithm's parameter list after if has been updated
 */

const GAParameterList& GAGeneticAlgorithm::setParameters (STD_ISTREAM& is, GABoolean flag) {

   params.read (is, flag);

   for (int i = 0; i < params.size (); i++)
      setptr (params [i].fullname (), params [i].value ());

   return params;

}

#endif


/**
 * Set a parameter given a pointer to the value.
 * Return 0 if everything is OK, non-zero if error.  If we did not set
 * anything then we return non-zero (this is not an error, but we indicate
 * that we did not do anything).  The set method must set both the GA's
 * parameter and the value in the parameter list (kind of stupid to maintain two
 * copies of the same data, but oh well).  The call to set on params is redundant
 * for the times when this method is called *after* the parameter list has been
 * updated, but it is necessary when this method is called directly by the user.
 *
 * @param name Name of the parameter
 * @param value Pointer to the parameter's value
 * @return 0 if Everything is OK, non-zero if error or if nothing was set
 */

int GAGeneticAlgorithm::setptr (const char* name, const void* value) {

   int status = 0;

   params.set (name, value); // redundant for some cases, but not others

   if (strcmp (name, gaNnBestGenomes ) == 0 ||
       strcmp (name, gaSNnBestGenomes) == 0    ) {

      nBestGenomes (*((int*) value));

   }
   else if (strcmp (name, gaNpopulationSize ) == 0 ||
            strcmp (name, gaSNpopulationSize) == 0    ) {

      populationSize (*((int*) value));

   }
   else if (strcmp (name, gaNnGenerations ) == 0 ||
            strcmp (name, gaSNnGenerations) == 0    ) {

      nGenerations (*((int*) value));

   }
   else if (strcmp (name, gaNpConvergence ) == 0 ||
            strcmp (name, gaSNpConvergence) == 0    ) {

      pConvergence (*((long double*) value));

   }
   else if (strcmp (name, gaNnConvergence ) == 0 ||
            strcmp (name, gaSNnConvergence) == 0    ) {

      nConvergence (*((int*) value));

   }
   else if (strcmp (name, gaNpCrossover ) == 0 ||
            strcmp (name, gaSNpCrossover) == 0    ) {

      pCrossover (*((long double*) value));

   }
   else if (strcmp (name, gaNpMutation ) == 0 ||
            strcmp (name, gaSNpMutation) == 0    ) {

      pMutation (*((long double*) value));

   }
   else if (strcmp (name,gaNscoreFrequency ) == 0 ||
            strcmp (name,gaSNscoreFrequency) == 0    ) {

      stats.scoreFrequency (*((int*) value));

   }
   else if (strcmp (name,gaNflushFrequency ) == 0 ||
            strcmp (name,gaSNflushFrequency) == 0    ) {

      stats.flushFrequency (*((int*) value));

   }
   else if (strcmp (name,gaNrecordDiversity ) == 0 ||
            strcmp (name,gaSNrecordDiversity) == 0    ) {

      stats.recordDiversity (*((int*) value) ? gaTrue : gaFalse);

   }
   else if (strcmp (name,gaNselectScores ) == 0 ||
            strcmp (name,gaSNselectScores) == 0    ) {

      stats.selectScores(*((int*)value));

   }
   else if (strcmp (name,gaNscoreFilename ) == 0 ||
            strcmp (name,gaSNscoreFilename) == 0    ) {

      char tmpname [64];
      memcpy (tmpname, value, strlen ((char*) value) + 1);
      stats.scoreFilename (tmpname);

   }
   else
      status = 1;

   return status;

}


// This is a pretty ugly little hack to make long doubles/long doubles work transparently.
int GAGeneticAlgorithm::set (const char* name, long double v) {

   int status = 1;

   for (int i = 0; i < params.size (); i++) {
      if (strcmp (name, params [i].fullname ()) == 0 || strcmp (name, params [i].shrtname ()) == 0) {
         if (params [i].type () == GAParameter::FLOAT) {

            long double fval = (long double) v;
            status = setptr (name, (void*) &fval);

         }
         else
            status = setptr (name, (void*) &v);
      }
   }

   return status;

}


/**
 * Get a parameter from the algorithm
 *
 * @param name name of the parameter
 * @param value pointer to the buffer in which to store the parameter's value
 */

int GAGeneticAlgorithm::get (const char* name, void* value) const {

   int status = 0;

   if (strcmp (name, gaNseed)  == 0 ||
       strcmp (name, gaSNseed) == 0    ) {

      *((int*) value) = d_seed;

   }
   else if (strcmp (name, gaNnBestGenomes ) == 0 ||
            strcmp (name, gaSNnBestGenomes) == 0    ) {

      *((int*) value) = stats.nBestGenomes ();

   }
   else if (strcmp (name, gaNpopulationSize ) == 0 ||
            strcmp (name, gaSNpopulationSize) == 0    ) {

      *((int*) value) = pop->size ();

   }
   else if (strcmp (name, gaNnGenerations ) == 0 ||
            strcmp (name, gaSNnGenerations) == 0    ) {

      *((int*) value) = ngen;

   }
   else if (strcmp (name, gaNpConvergence ) == 0 ||
            strcmp (name, gaSNpConvergence) == 0    ) {

      *((long double*) value) = pconv;

   }
   else if (strcmp (name, gaNnConvergence ) == 0 ||
            strcmp (name, gaSNnConvergence) == 0    ) {

      *((int*) value) = nconv;

   }
   else if (strcmp (name, gaNpCrossover ) == 0 ||
            strcmp (name, gaSNpCrossover) == 0    ) {

      *((long double*) value) = pcross;

   }
   else if (strcmp (name, gaNpMutation ) == 0 ||
            strcmp (name, gaSNpMutation) == 0    ) {

      *((long double*) value) = pmut;

   }
   else if (strcmp (name, gaNscoreFrequency ) == 0 ||
            strcmp (name, gaSNscoreFrequency) == 0    ) {

      *((int*) value) = stats.scoreFrequency ();

   }
   else if (strcmp (name, gaNflushFrequency ) == 0 ||
            strcmp (name, gaSNflushFrequency) == 0    ) {

      *((int*) value) = stats.flushFrequency ();

   }
   else if (strcmp (name, gaNrecordDiversity ) == 0 ||
            strcmp (name, gaSNrecordDiversity) == 0    ) {

      *((int*) value) = stats.recordDiversity ();

   }
   else if (strcmp (name, gaNselectScores ) == 0 ||
            strcmp (name, gaSNselectScores) == 0    ) {

      *((int*) value) = stats.selectScores ();

   }
   else if (strcmp (name, gaNscoreFilename ) == 0 ||
            strcmp (name, gaSNscoreFilename) == 0    ) {

      *((const char**) value) = stats.scoreFilename ();

   }
   else
      status = 1;

   return status;

}


/**
 * Set the objective function on all individuals used by the genetic algorithm.
 * This can be changed during the course of an evolution.
 *
 * @param f Objective function
 */

void GAGeneticAlgorithm::objectiveFunction (GAGenome::Evaluator f) {

   for (unsigned int i = 0; i < pop->size (); i++)
      pop->individual (i).evaluator (f);

}


/**
 * Set the objective data member on all individuals used by the genetic algorithm.
 * This can be changed during the course of an evolution.
 *
 * @param v Objective function data
 */

void GAGeneticAlgorithm::objectiveData (const GAEvalData& v) {

   for (unsigned int i = 0; i < pop->size (); i++)
      pop->individual (i).evalData (v);

}


/**
 * Set the population. Returns a reference to the current population.
 *
 * @param p Source population
 * @return Reference to the current popultion
 */

const GAPopulation& GAGeneticAlgorithm::population (const GAPopulation& p) {

   if (p.size () < 1) {

      GAErr (GA_LOC, className (), "population", gaErrNoIndividuals);
      return *pop;

   }

   pop->copy (p);
   pop->geneticAlgorithm (*this);

   return *pop;

}


/**
 * Set the population size. This can be changed during the course of an evolution.
 *
 * @param value New population size
 * @return The new population size
 */

unsigned int GAGeneticAlgorithm::populationSize (unsigned int value) {

   unsigned int ps = value;

   params.set (gaNpopulationSize, value);

   return pop->size (ps);

}


/**
 * Returns a pointer to the user data member
 *
 * @return The algorithm's user data pointer
 */

void* GAGeneticAlgorithm::getUserData () const {

   return ud;

}


/**
 * Sets the user data member. This is a general purpose
 * void pointer, application-dependent
 *
 * @param d User data pointer
 */

void* GAGeneticAlgorithm::setUserData (void* d) {

   return ud = d;

}

void GAGeneticAlgorithm::updateNoStepsStats(){
  stats.updateNoStepsStats(*pop);
}

/**
 * Returns the terminator function
 *
 * @return The terminator function used in this algorithm
 */

GAGeneticAlgorithm::Terminator GAGeneticAlgorithm::getTerminator () const {

   return cf;

}


/**
 * Sets the terminator function
 *
 * @param f New terminator function
 */

GAGeneticAlgorithm::Terminator GAGeneticAlgorithm::setTerminator (Terminator f) {

   return cf = f;

}


void GAGeneticAlgorithm::printStats(string message){
  GALogger::instance()->appendStats(message, population());
}


/**
 * Initialize the genetic algorithm then evolve it until the termination
 * criteria have been satisfied. This function first calls initialize() then
 * calls the step member function until the done member function returns
 * gaTrue. It calls the flushScores() member as needed when the evolution
 * is complete. You may pass a seed to evolve if you want to specify
 * your own random seed.
 *
 * @param seed Random seed, as used in initialize()
 */

inline void GAGeneticAlgorithm::evolve (unsigned int seed) {
  // initialized call has been moved to Algorithm.cc
  //printStats("Initial Stats");
  while (!done ()){
    step ();
    //printStats("End of Step Stats");
  }

  if (stats.flushFrequency () > 0)
    stats.flushScores ();

  // This was added so that algorithms such as MTS that do not add nor remove
  // any individuals to the population display the correct best individual
  // (and not the first one, as it would do without this line)
  pop->touch();

}


/*
 * For storing the statistics that come from the application of variation operator (crossover, mutation in gas).
 * This cannot be done later with the population since some of the individuals generated are eliminated by the elitism
 * procedure. There is a relation 1:1 between parents and children, in the way that parent(i) and parent(i+1) generated
 * (by crossover or by copying) children i and i+1. With odd population sizes, children(n) and children(n-1) are the same
 * (since in this case only one child gets generated by the two parents). The reason of taking out this code from the
 * population is to avoid an extra dependency between population and the Communitcation Manager.
 *
 */
void GAGeneticAlgorithm::updateVarOperatorsStats(const GAPopulation& old_pop,
						 const vector<GAGenome*>& parents,
                                                 vector<GAGenome*>& children) {
  // Maybe this method should change its name or create a new method for the lines below
  for (vector<GAGenome*>::iterator it=children.begin(); it!=children.end(); it++){
    (*it)->genBirth(stats.generation());
    (*it)->setNewIndAge();
    (*it)->origin(rank);
  }
  stats.updateFitnessIncrement(old_pop,parents,children);
}



// Maybe a better idea with island_rank_ would be to move it to the GAGenome class and with a class method
// initialize this value within the main method (and have something similar to the setNewIndAge method)
void GAGeneticAlgorithm::setBirthData(GAPopulation& children){
  int curgen = stats.generation();

  for (unsigned int i = 0; i < children.size(); i++){
    children.individual(i).genBirth     ( curgen);
    children.individual(i).setNewIndAge ();
    children.individual(i).origin       ( rank );
  }
}

void GAGeneticAlgorithm::recombinator(Recombinator* recombinator){
  recombinator_ = recombinator;
}
