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
 * @brief GAGeneticAlgorithm class hdr.
 *
 * Header for the base genetic algorithm class
 */

#ifndef GABASEGA_H
#define GABASEGA_H

/* INCLUDES */
#include <stdio.h>
#include <string.h>

#include "Algorithm.h"
#include "gaversion.h"  // gets the RCS string in for ident purposes
#include "gaconfig.h"
#include "gaparameters.h"
#include "GAParameter.h"
#include "GAStatistics.h"
#include "GAPopulation.h"
#include "genomes/GAGenome.h"

class Recombinator;

/**
 * @brief Genetic Algorithms base class (abstract)
 *
 * This is an abstract class that cannot be instantiated. Each genetic
 * algorithm, when instantiated, will have default operators defined for
 * it. See the documentation for the specific genetic algorithm type for
 * details.
 *
 * The base genetic algorithm class keeps track of evolution statistics
 * such as number of mutations, number of crossovers, number of
 * evaluations, best/mean/worst in each generation, and initial/current
 * population statistics. It also defines the terminator, a member function
 * that specifies the stopping criterion for the algorithm.
 *
 * You can maximize or minimize by calling the appropriate member
 * function. If you derive your own genetic algorithm, remember that users
 * of your algorithm may need either type of optimization.
 *
 * Statistics can be written to file each generation or periodically by
 * specifying a flush frequency. Generational scores can be recorded each
 * generation or less frequently by specifying a score frequency.
 *
 * Parameters such as generations-to-completion, crossover probability and
 * mutation probability can be set by member functions, command-line, or
 * from file.
 *
 * The evolve member function first calls initialize then calls the step
 * member function until the done member function returns gaTrue. It calls
 * the flushScores member as needed when the evolution is complete. If you
 * evolve the genetic algorithm without using the evolve member function,
 * be sure to call initialize before stepping through the evolution. You
 * can use the step member function to evolve a single generation. You
 * should call flushScores when the evolution is finished so that any
 * buffered scores are flushed.
 *
 * The names of the individual parameter member functions correspond to
 * the #defined string names. You may set the parameters on a genetic
 * algorithm one at a time (for example, using the nGenerations member
 * function), using a parameter list (for example, using the parameters
 * member function with a GAParameterList), by parsing the command line
 * (for example, using the parameters member function with argc and argv),
 * by name-value pairs (for example, using the set member function with a
 * parameter nam
 *
 */

class GAGeneticAlgorithm : public Algorithm {

   public:

      GADefineIdentity ("GAIncrementalGA", GAID::BaseGA);

      typedef GABoolean (*Terminator) (GAGeneticAlgorithm& ga);

      enum {MINIMIZE = -1, MAXIMIZE = 1 };

      static GAParameterList& registerDefaultParameters (GAParameterList&);

      static GABoolean TerminateUponGeneration     (GAGeneticAlgorithm&);
      static GABoolean TerminateUponConvergence    (GAGeneticAlgorithm&);
      static GABoolean TerminateUponPopConvergence (GAGeneticAlgorithm&);

      static GABoolean TerminateUponPopConvergenceOrGen   (GAGeneticAlgorithm&);
      static GABoolean TerminateUponPartialPopConvergence (GAGeneticAlgorithm&);
      static GABoolean TerminateUponScoreOrGen            (GAGeneticAlgorithm&);
      static GABoolean TerminateUponFitnessEvaluations    (GAGeneticAlgorithm&);

   public:

      GAGeneticAlgorithm (const GAGenome&);
      GAGeneticAlgorithm (const GAPopulation&);
      GAGeneticAlgorithm (const GAGeneticAlgorithm&);

      GABoolean done ();

      // Virtual functions - documented here

      /**
       * Destructor
       */
      virtual ~GAGeneticAlgorithm ();

      virtual void copy (const GAGeneticAlgorithm&);


      /**
       * Initialize the genetic algorithm. If you specify a seed, this function
       * calls GARandomSeed with that value. If you do not specify a seed, GAEDAlib will
       * choose one for you. It then initializes the population and does the first
       * population evaluation.
       *
       * @param seed  Random seed
       */
      virtual void initialize (unsigned int seed = 0) = 0;

      virtual void offspring (GAPopulation* offpop) = 0;
      //virtual void offspring (GAPopulation* offpop) {};

      /**
       * Evolve the genetic algorithm for one generation.
       */
      virtual void step() = 0;

      virtual void evolve (unsigned int seed = 0);


#ifdef GALIB_USE_STREAMS
      virtual int write (const char*) const {return 0;}
      virtual int write (STD_OSTREAM&) const {return 0;}

      virtual int read (const char*) {return 0;}
      virtual int read (STD_ISTREAM&) {return 0;}
#endif

      void* getUserData () const;
      void* setUserData (void* d);

      Terminator getTerminator () const;
      Terminator setTerminator (Terminator f);

      const GAParameterList& getParameters () const;
      const GAParameterList& setParameters (const GAParameterList&);
      const GAParameterList& setParameters (int&, char**, GABoolean flag = gaFalse);


#ifdef GALIB_USE_STREAMS
      const GAParameterList& setParameters (const char* filename, GABoolean f = gaFalse);
      const GAParameterList& setParameters (STD_ISTREAM&, GABoolean flag = gaFalse);
#endif


      virtual int get (const char*, void*) const;
      virtual int setptr (const char*, const void*);

      int set (const char* s, int v) {return setptr (s, (void*) &v);}
      int set (const char* s, unsigned int v) {return setptr (s, (void*) &v);}
      int set (const char* s, char v) {return setptr (s, (void*) &v);}
      int set (const char* s, const char* v) {return setptr (s, (void*) v);}
      int set (const char* s, const void* v) {return setptr (s, v);}
      int set (const char* s, long double v);

      int nGenerations () const;
      int nGenerations (unsigned int n);
      int nConvergence () const;
      int nConvergence (unsigned int n);
      int nEvaluations () const;
      int nEvaluations (unsigned int n);

      long double pConvergence () const;
      long double pConvergence (long double p);
      long double pCrossover () const;
      long double pCrossover (long double p);
      long double pMutation () const;
      long double pMutation (long double p);

      GAGenome::SexualCrossover crossover (GAGenome::SexualCrossover f);
      GAGenome::SexualCrossover sexual () const {return scross;}

      GAGenome::AsexualCrossover crossover (GAGenome::AsexualCrossover f) {return across = f;}
      GAGenome::AsexualCrossover asexual () const {return across;}

      const GAStatistics& statistics () const {return stats;}
      long double convergence () const {return stats.convergence ();}
      int   generation  () const {return stats.generation  ();}

      void  flushScores () {

         if (stats.flushFrequency () > 0)
            stats.flushScores ();

      }

      int scoreFrequency () const {return stats.scoreFrequency ();}
      int scoreFrequency (unsigned int x) {

         params.set (gaNscoreFrequency, x);
         return stats.scoreFrequency (x);

      }

      int flushFrequency () const {return stats.flushFrequency ();}
      int flushFrequency (unsigned int x) {

         params.set (gaNflushFrequency, x);
         return stats.flushFrequency (x);

      }

      const char* scoreFilename () const {return stats.scoreFilename ();}
      const char* scoreFilename (const char* fn) {

         params.set (gaNscoreFilename, fn);
         return stats.scoreFilename (fn);

      }

      int selectScores () {return stats.selectScores ();}
      int selectScores (int w) {

         params.set (gaNselectScores, w);
         return stats.selectScores (w);

      }

      GABoolean recordDiversity () const {return stats.recordDiversity ();}
      GABoolean recordDiversity (GABoolean f) {

         params.set (gaNrecordDiversity, (int) f);
         return stats.recordDiversity (f);

      }

      virtual const GAPopulation& population () const {return *pop;}
      virtual const GAPopulation& population (const GAPopulation&);

      virtual unsigned populationSize () const {return pop->size ();}
      virtual unsigned populationSize (unsigned int n);

      virtual int nBestGenomes () const {return stats.nBestGenomes ();}
      virtual int nBestGenomes (unsigned int n) {

         params.set (gaNnBestGenomes, n);
         return stats.nBestGenomes (pop->individual (0), n);

      }

      virtual GAScalingScheme& scaling () const {return pop->scaling ();}
      virtual GAScalingScheme& scaling (const GAScalingScheme& s) {return pop->scaling (s);}

      virtual GASelectionScheme& selector () const {return pop->selector ();}
      virtual GASelectionScheme& selector (const GASelectionScheme& s) {return pop->selector (s);}

      virtual void objectiveFunction (GAGenome::Evaluator f);
      virtual void objectiveData     (const GAEvalData& v);

     void updateNoStepsStats();

      void updateVarOperatorsStats (const GAPopulation& oldpop,
				    const vector<GAGenome*>& parents,
				    vector<GAGenome*>& children );

      GAGenome* best() const { return pop->best().clone();}

      void printStats(string message);

      //void islandRank(unsigned island_rank);

      int getRank() {return rank;}
      void setRank(int newRank) {rank=newRank;} //   void islandRank(unsigned island_rank); did the same

      void setBirthData (GAPopulation& children);

      void recombinator(Recombinator* recombinator);

   protected:

      GAStatistics    stats;   ///< Runtime statistics
      GAParameterList params;  ///< Algorithm parameters
      GAPopulation*   pop;     ///< Genetic population
      Terminator      cf;      ///< Function for determining done-ness
      void*           ud;      ///< Pointer to user data structure

      int             d_seed;
      unsigned int    ngen;    ///< Number of generations
      unsigned int    nconv;   ///< Generations to convergence
      long double           pconv;   ///< Convergence percentage
      long double           pcross;  ///< Crossover probability
      long double           pmut;    ///< Mutation probability

      GAGenome::SexualCrossover  scross; ///< Sexual crossover to use
      GAGenome::AsexualCrossover across; ///< Asexual crossover to use

      int rank; //island_rank_

      unsigned int nevals;      ///< Maximum number of fitness evaluations

      Recombinator*     recombinator_;  // GAGeneticAlgorithm gets ownership, i.e it deletes in the destructor
};


/* INLINE FUNCTIONS */

// Inline functions should be defined in the header file. That's why
// we also document them here


/**
 * Checks the end of evolution
 */

inline GABoolean GAGeneticAlgorithm::done (){

   GABoolean result;

   result = (*cf) (*this);

   return result;

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

/* inline void GAGeneticAlgorithm::evolve (unsigned int seed) { */
/*   // initialized call has been moved to Algorithm.cc */
/*   printStats("Initial Stats"); */
/*   //   initialize (seed); */

/*   while (!done ()){ */
/*     step (); */
/*     printStats("End of Step Stats"); */
/*   } */

/*    if (stats.flushFrequency () > 0) */
/*       stats.flushScores (); */

/* } */


/**
 * Get the number of generations.
 */

inline int GAGeneticAlgorithm::nGenerations () const {

   return ngen;

}


/**
 * Set the number of generations.
 *
 * @param n number of iterations
 */

inline int GAGeneticAlgorithm::nGenerations (unsigned int n) {

   params.set (gaNnConvergence, n);

   return ngen = stats.nConvergence (n);

}


/**
 * Get the number of evaluations.
 */

inline int GAGeneticAlgorithm::nEvaluations () const {

   return nevals;

}


/**
 * Set the number of evaluations.
 *
 * @param n number of evaluations
 */

inline int GAGeneticAlgorithm::nEvaluations (unsigned int n) {

   return nevals = n;

}


/**
 * Get the number of generations used for the convergence test.
 */

inline int GAGeneticAlgorithm:: nConvergence () const {

   return nconv;

}


/**
 * Set the number of generations used for the convergence test.
 *
 * @param n Number of generations
 */

inline int GAGeneticAlgorithm::nConvergence (unsigned int n) {

   params.set (gaNnConvergence, n);

   return nconv = stats.nConvergence (n);

}


/**
 * Get the convergence percentage. The convergence is defined as the
 * ratio of the Nth previous best-of-generation score to the current best-of-generation score.
 * N is defined by the nConvergence() member function.
 */

inline long double GAGeneticAlgorithm::pConvergence () const {

   return pconv;

}


/**
 * Set the convergence percentage
 *
 * @param p Convergence percentage
 */

inline long double GAGeneticAlgorithm::pConvergence (long double p) {

   params.set (gaNpConvergence, p);

   return pconv = p;

}


/**
 * Get the crossover probability.
 */

inline long double GAGeneticAlgorithm::pCrossover () const {

   return pcross;

}


/**
 * Set the crossover probability.
 *
 * @param p Crossover probability
 */

inline long double GAGeneticAlgorithm::pCrossover (long double p) {

   params.set (gaNpCrossover, p);

   return pcross = p;

}


/**
 * Get the mutation probability.
 */

inline long double GAGeneticAlgorithm::pMutation () const {

   return pmut;

}


/**
 * Set the mutation probability.
 *
 * @param p Mutation probability
 */

inline long double GAGeneticAlgorithm::pMutation (long double p) {

   params.set (gaNpMutation, p);

   return pmut = p;

}


/**
 * Specify the mating method to use for evolution. This can be
 * changed during the course of an evolution. This genetic algorithm uses only
 * sexual crossover.
 *
 * @param f Mating method
 */

inline GAGenome::SexualCrossover GAGeneticAlgorithm::crossover (GAGenome::SexualCrossover f) {

   return scross = f;

}


#endif
