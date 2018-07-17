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

#ifndef CONFIG_H
#define CONFIG_H

#include "MOSConfig.h"
#include "logger/GALogger.h"
#include "extras/problemloader.h"
#include "genomes/GAGenome.h"

class CommManager;
class GAScalingScheme;
class GASelectionScheme;
class MOSQuality;

class GAEDAConfig {

   public:

      typedef enum {MOS} algorithmtype;
      typedef enum {ITERS, CONV, POPCONV,  PARTIALPOPCONV, POPCONVITERS, SCOREORGEN, EVALS} convergencecriterion;
      typedef enum {FullPrinting, DefaultPrinting} mos_genome_printing_t;

      enum logstat_t {Score, Fitness, FitInc, Participation, Quality, QualityFuncActive, Evals, Improve};

      enum selector_t {RandomSelec, Roulette, Tournament, Unif_Tournament, Rank, Stochastic, Deterministic};
      enum scaling_t  {None, Standard, Linear, Sigma_Truncation, Power_Law};

      static GAEDAConfig* handle(CommManager& comm_manager);
      static GAEDAConfig* handle();

      static void destroy();

      // Parse command line arguments
      bool parse (int argc, char **argv, bool embedded = false);

      // Show command line usage
      void showUsage (char *prgname);

      // General options
      bool                   quiet               (void) const {return mQuiet;        }
      char*                  getLogFile          (void) const {return mLogFile;      }
      char*                  setLogFile          (const char* logFile) {return mLogFile = strdup(logFile);}
      GALogger::LogLevel     getLogLevel         (void) const {return mLogLevel;     }
      vector<logstat_t>      getLogStatsNames    (void) const {return logstats_names;}
      unsigned int           getStatsDisplayFreq (void) const {return mStatsDispFreq;}

      algorithmtype          getAlgorithm        (void) const {return mAlg;}

      int                    getPopSize          (void) const {return mPopSize;}

      int                    getProblemSize      (void) const {return mProbSize;}
      int                    setProblemSize      (int sz) {return mProbSize = sz;}

      GAGenome::OptCriterion getOptCriterion     (void) const;
      char*                  getProblemData      (void) const {return mProblemData;}

      // Convergence options
      convergencecriterion getConvCrit (void) const {return mConvCrit;}

      int            getIters      (void) const {return mIters;     }
      long double    getConv       (void) const {return mConv;      }
      int            getConvIters  (void) const {return mConvIters; }
      unsigned long  getEvals      (void) const {return mEvals;     }
      unsigned long  setEvals      (unsigned long evals) {return mEvals = evals;}

      // Genetic algorithms options
      long double        getMutProb   (void) const {return mMutProb;}
      long double        getCrossProb (void) const {return mCrossProb;}
      GASelectionScheme* getSelector  (void) const;
      GAScalingScheme*   getScaling   (void) const;

      // MOS options
      const char*           getTechsRepository   (void) const {return techsRepository;}
      const vector<string>* getTechsToUse        (void) const {return &techsToUse;}
      long double           getMinPart           (void) const {return minPart;}
      mos_genome_printing_t getMOSGenomePrinting (void) const {return mosGenomePrinting;}
      std::string           getQualityMeasure    (void) const {return mQualityMeasure;}
      long double           getBonus             (void) const {return mBonus;}
      unsigned              getSharedEvals       (void) const {return mSharedEvals;}
      unsigned              setSharedEvals       (unsigned sharedEvals) {return mSharedEvals = sharedEvals;}
      long double           getSharedEvalsPercent(void) const {return mSharedEvalsPercent;}
      long double           setSharedEvalsPercent(long double sharedEvalsPercent) {return mSharedEvalsPercent = sharedEvalsPercent;}
      MOSQuality*           getQualityFunction   (void) const;

      // Debug options
      unsigned getSeed  (void) const {return seed;}
      int      getDebug (void) const {return loop;}

      // Problem file, struct and options
      char*             getProblemFile     (void) const {return mProblemFile;}
      GAProblemStruct*  getProblemStruct   (void) const {return mProblemStruct;}

      GAGenome*             getOptimum         (void) const;

      GAGenome::Initializer getIndInitFunction (void) const {return mProblemStruct->individualInit;}
      distToOptFuncPtr      getDistToOptFunc   (void) const {return mProblemStruct->distToOpt;}
      nComponentsOptFuncPtr getNComponentsOpt  (void) const {return mProblemStruct->nComponentsOpt;}

   protected:

      GAEDAConfig (CommManager& comm_manager);
      virtual ~GAEDAConfig();

      static GAEDAConfig* ptrSelf;

   private:

      // General options
      bool                  mQuiet;             // Useful for not printing information when embedding the algorithm
      char*                 mLogFile;           // Path to the log file
      GALogger::LogLevel    mLogLevel;          // Log level used bu the tracer
      vector<logstat_t>     logstats_names;
      unsigned              mStatsDispFreq;     // The frequency of generations for displaying the stats

      scaling_t             scaling_;

      algorithmtype         mAlg;               // Algorithm

      int                   mPopSize;           // Population size

      int                   mProbSize;          // Problem size

      char*                 mProblemData;       // Additional problem data

      // Convergence options
      convergencecriterion  mConvCrit;          // Convergence criterion

      int                   mIters;             // Iteration limit
      long double           mConv;              // Convergence limit
      int                   mConvIters;         // Iterations used to calculate the convergence
      unsigned long         mEvals;             // Evaluations limit

      // Genetic algorithms options
      long double           mMutProb;           // Mutation probability
      long double           mCrossProb;         // Crossover probability
      selector_t            selector_;          // The Selection Scheme
      int                   selector_degree_;

      // MOS options
      char*                 techsRepository;
      vector<string>        techsToUse;
      long double           minPart;            // MOS min participation
      mos_genome_printing_t mosGenomePrinting;  // Type of printing for MOS genome
      std::string           mQualityMeasure;    // Quality measure used for dynamic adjustment of participation
      long double           mBonus;
      unsigned              mSharedEvals;       // Number of shared evaluations in a MOS HRH algorithm
      long double           mSharedEvalsPercent;// Percentage of shared evaluations in a MOS HRH algorithm

      // Debug options
      unsigned              seed;               // the seed of the random method
      int                   loop;

      // Problem structure
      char*                 mProblemFile;       // Proble file
      GAProblemStruct*      mProblemStruct;     // Problem struct

      // Communication manager
      CommManager&           comm_manager_;

};

#endif
