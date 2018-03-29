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

#include "GAPopulation.h"
#include "MOSConfig.h"
#include "logger/GALogger.h"
#include "logger/GAFileLogger.h"
#include "extras/problemloader.h"
#include "islands/neighborconds.h"
#include "routingdist/DistRoutingAlg.h"

class CommManager;
class GAGenealogy;
class Recombinator;
class CentroidsGenerator;
class GAScalingScheme;
class GASelectionScheme;
class EAIslandsTopology;
class EAEmigrantsSelector;
class MOSParticipation;
class MOSQuality;
class VNSOp;
class VNSShaker;

class GAEDAConfig {

   public:

      typedef enum {SYNC, ASYNC, ASYNCFAST, CLUST, SEQ, ROUTINGDISTSYNC} modeltype;
      typedef enum {GA, SSGA, EDA, DE, ES, MOS, MOS2, MOSMultiDeme, MOSRL, VNS, TSCordeau} algorithmtype;
      typedef enum {ITERS, CONV, POPCONV,  PARTIALPOPCONV, POPCONVITERS, POPCONVEVALS, SCOREORGEN, EVALS, EVALROUTECALLS, TIMECONV} convergencecriterion;
      typedef enum {Mesh, a2a, Ring, Ring2, Pairs, NearestNeighbor, FurthestNeighbor,
                    RandomTopology, HyperCubeTopology, HybridTopology} topology;
      typedef enum {NoGenealogy, Tracer, Memory, All, Full} genealogy;
      typedef enum {PMAX, PHC, WOLF} mosrl_policy;
      typedef enum {FullPrinting, DefaultPrinting} mos_genome_printing_t;

      enum neighbor_conds    {Gabriel, Rel_neighbors, Square};
      enum centroids_gen_t   {HighSeparatedRandom, Random};
      enum islands_init_t    {VoronoiIslandInit, RandomIslandInit};
      enum emigrants_selec_t {BestEmigrantsSelector, RandomEmigrantsSelector};

      enum logstat_t {Score, Fitness, Age, OptNcomps, OptDist, FitInc, NumNewChild, NativePrcnt, GenDiv,
                      Entropy, GrefBias, DistAvgPoint, MaxDivDynCoevLogStat, SupDivDynCoevLogStat,
                      InfDivDynCoevLogStat, AvgDivIncDynCoevLogStat, CoevPhaseLogStat, Epistasis,
                      Participation, Quality, Genealogy, Evals, FSSAvg, FSSBest, FSSWorst, Improve, RoutingAlgScore,
                      VNSSuccessfulShaker, VNSShakerpos, EvalRouteCalls, RoutingAlgBestSolPickupTime,
                      RoutingAlgBestSolDeliveryTime, RoutingAlgBestSolNonPenScore, Runtime, Freemem, RoutingAlgBestSolTWV, RoutingAlgBestSolLoadV,
                      RoutingAlgBestSolRideV};

      enum recombinator_t {SingleElitism, FullElitism, DE_Elitism, Percent_Elitism};
      enum selector_t     {RandomSelec, Roulette, Tournament, Unif_Tournament, Rank, Stochastic, Deterministic};
      enum scaling_t      {None, Standard, Linear, Sigma_Truncation, Power_Law};

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

      long double                 getElisitsmPercent  (void) const {return elitism_percentage;}
      Recombinator*          getRecombinator     (void) const;

      algorithmtype          getAlgorithm        (void) const {return mAlg;}

      int                    getProblemSize      (void) const {return mProbSize;}
      int                    setProblemSize      (int sz) {return mProbSize = sz;}

      GAGenome::OptCriterion getOptCriterion     (void) const;

      char*                  getProblemData      (void) const {return mProblemData;}

      int                    getRepetitions      (void) const {return mRepetitions;}

      long double                 getOptimumFitness   (void) const {return mOptimumFitness;}
      long double                 setOptimumFitness   (long double fit) {return mOptimumFitness = fit;}

      int                    getFunction         (void) const {return mFunction;}

      int                    getCMAESMaxPopSize  (void) const {return mCMAESMaxPopSize;}
      std::string            getCMAESInitials    (void) const {return mCMAESInitials;}
      std::string            getBBOBOutputPath   (void) const {return mBBOBOutputPath;}
      std::string            getBBOBConfigName   (void) const {return mBBOBConfigName;}

      long double                 getResetPercentage  (void) const {return mResetPerc;}
      long double                 getRerunPercentage  (void) const {return mRerunPerc;}

      bool                   getPrintBestSol     (void) const {return mPrintBestSol;}

      // Convergence options
      convergencecriterion getConvCrit (void) const {return mConvCrit;}

      int            getIters             (void) const {return mIters;     }
      long double         getConv              (void) const {return mConv;      }
      int            getConvIters         (void) const {return mConvIters; }
      unsigned long  getEvals             (void) const {return mEvals;     }
      unsigned long  setEvals             (unsigned long evals) {return mEvals = evals;}
      int            getMaxEvalRouteCalls (void) const {return mEvalRouteCalls; }
      long double         getMaxExecSeconds    (void) const {return maxExecSeconds; }
      long double         getPrecission        (void) const {return mPrecission;}
      long double         setPrecission        (long double precission) {return mPrecission = precission;}
      bool           usePrecission        (void) const {return mUsePrecission;}

      // Parallel model options
      int                  getPopSize               (void) const {return mPopSize;}
      EAIslandsTopology*   getTopology              (void) const;
      modeltype            getModel                 (void) const {return mModel;}
      int                  getMigrationFrequency    (void) const {return mMF;}
      int                  getMigrationPopSize      (void) const {return mMigPopSize;}

      EAEmigrantsSelector* getEmigrantsSelector     (GAGeneticAlgorithm& ea) const;

      unsigned             getTopDegree             (void) const;

      neighbor_cond_func   getNeighborConditionFunc (void) const;
      unsigned             getMinNeighborsNum       (void) const;
      unsigned             getMaxNeighborsNum       (void) const;

      unsigned             getClusteringPeriod      (void) const {return clustering_period_;}

      CentroidsGenerator*  getCentroidsGenerator    (GAGenome& gen) const;

      int                  getIslandInit            (void) const;

      bool                 useCentroids             (void) const {return mUseCentroids;}
      long double               getCentroidProb          (void) const {return mCentProb;}
      long double               getCentroidMinDist       (void) const {return mMinCentDist;}
      int                  getCentroidMaxIters      (void) const {return mMaxCentIters;}

      // Distributed routing algorithm options
      DistRoutingAlg::emigrantRouteType getEmigrantRouteType()    const { return emigrant_route_t_; }


      // Genetic algorithms options
      long double             getMutProb   (void) const {return mMutProb;}
      long double             getCrossProb (void) const {return mCrossProb;}
      GASelectionScheme* getSelector  (void) const;
      GAScalingScheme*   getScaling   (void) const;

      // Diferential evolution options
      long double getFFactor  () const {return DE_F_factor;}
      long double getCRFactor () const {return DE_CR_factor;}
      GAGenome::DECrossover getDECrossover() const {return CrossOp;}

      // Evolution strategies options
      unsigned getMu     () const {return ES_mu;    }
      unsigned getRo     () const {return ES_ro;    }
      unsigned getLambda () const {return ES_lambda;}

      // VNS options
      VNSOp*                    getRoutingLocalSearch()         const;
      std::vector< VNSShaker* > getVNSShakersList()             const;
      long double                    getVNSTinitRatio()              const;
      long double                    getVNSTinitProb()               const;
      long double                    getVNSLSRatio1()                const;
      long double                    getVNSLSProb()                  const;
      long double                    getVNSLSRatio2()                const;
      int                       getVNSMaxEvalsShaker()          const;
      bool                      getVNSUseAllShakersInEachStep() const;
      int                       getVNSMaxShakersRes()           const;
      bool                      getVNSUpdatePenalizations()     const;

      // MTS options
      long double getMTSAdjustFailed() const {return mts_adjustFailed;}
      long double getMTSAdjustMin   () const {return mts_adjustMin;}
      long double getMTSMoveLeft    () const {return mts_moveLeft;}
      long double getMTSMoveRight   () const {return mts_moveRight;}

      // MTS Reduced options
      long double getMTSRedSearchProb() const {return mtsred_searchProb;}
      long double getMTSRedMinProb   () const {return mtsred_minProb; }

      // Solis-Wets options
      long double getSWMaxSuccess   () const {return sw_maxSuccess;}
      long double getSWMaxFailed    () const {return sw_maxFailed;}
      long double getSWAdjustSuccess() const {return sw_adjustSuccess;}
      long double getSWAdjustFailed () const {return sw_adjustFailed;}
      long double getSWDelta        () const {return sw_delta;}

      // Adaptive DE options
      long double getAdapDEFl   () const {return adapde_Fl;}
      long double getAdapDEFu   () const {return adapde_Fu;}
      long double getAdapDETauF () const {return adapde_tauF;}
      long double getAdapDETauCR() const {return adapde_tauCR;}

      // Adaptive STSDE options
      long double getSTSDEProb() const {return adapde_Fl;}

      // MOS options
      const char*           getTechsRepository   (void) const {return techsRepository;}
      const vector<string>* getTechsToUse        (void) const {return &techsToUse;}
      ParticipationFunction getPartFunction      (void) const {return partFunction;}
      long double                getMinPart           (void) const {return minPart;}
      evolutionType         getEvolutiveApproach (void) const {return evolutiveApproach;}
      mos_genome_printing_t getMOSGenomePrinting (void) const {return mosGenomePrinting;}
      long double                getMOSQLPMax         (void) const {return mMOSQLPmax;}
      long double                getMOSRLAlfa         (void) const {return mMOSRLAlfa;}
      long double                getMOSRLBeta         (void) const {return mMOSRLBeta;}
      long double                getMOSRLGamma        (void) const {return mMOSRLGamma;}
      long double                getMOSRLDelta        (void) const {return mMOSRLDelta;}
      long double                getMOSRLDeltaL       (void) const {return mMOSRLDeltaL;}
      long double                getMOSRLDeltaW       (void) const {return mMOSRLDeltaW;}
      mosrl_policy          getMOSRLPolicy       (void) const {return mMOSRLPolicy;}
      GAGenealogy*          getGenealogy         (void) const;
      bool                  printGenealogy       (void) const {return mPrintGenealogy;}
      std::string           getQualityMeasure    (void) const {return mQualityMeasure;}
      bool                  weightedAutonomic    (void) const {return mWeightedAutonomic;}
      long double                getBonus             (void) const {return mBonus;}
      unsigned              getSharedEvals       (void) const {return mSharedEvals;}
      unsigned              setSharedEvals       (unsigned sharedEvals) {return mSharedEvals = sharedEvals;}
      long double                getSharedEvalsPercent(void) const {return mSharedEvalsPercent;}
      long double                setSharedEvalsPercent(long double sharedEvalsPercent) {return mSharedEvalsPercent = sharedEvalsPercent;}
      long double                getSharedEvalsFactor (void) const {return mSharedEvalsFactor;}
      long double                setSharedEvalsFactor (unsigned sharedEvalsFactor) {return mSharedEvalsFactor = sharedEvalsFactor;}
      std::string           getHetMOSCfgFile     (void) const {return mHetMOSCfgFile;}
      MOSParticipation*     getParticipationFunction (void) const;
      MOSQuality*           getQualityFunction   (void) const;

      // Debug options
      unsigned getSeed  (void) const {return seed;}
      int      getDebug (void) const {return loop;}

      // Problem file, struct and options
      char*                 getProblemFile     (void) const {return mProblemFile;}
      GAProblemStruct*      getProblemStruct   (void) const {return mProblemStruct;}

      GAGenome*             getOptimum         (void) const;

      GAGenome::Initializer     getIndInitFunction (void) const {return mProblemStruct->individualInit;}
      GAPopulation::Initializer getPopInitFunction (void) const {return mProblemStruct->populationInit;}
      distToOptFuncPtr      getDistToOptFunc   (void) const {return mProblemStruct->distToOpt;}
      nComponentsOptFuncPtr getNComponentsOpt  (void) const {return mProblemStruct->nComponentsOpt;}

      // Utility functions
      bool isSequential (void) const {return (mModel == GAEDAConfig::SEQ);}
      unsigned getRestarts (void) const {return mRestarts;}
      unsigned increaseRestarts (void) {return mRestarts++;}
      unsigned setRestarts (unsigned restarts) {return mRestarts = restarts;}

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

      long double                elitism_percentage;
      recombinator_t        recombinator_;
      scaling_t             scaling_;

      algorithmtype         mAlg;               // Algorithm

      int                   mProbSize;          // Problem size

      char*                 mProblemData;       // Additional problem data

      int                   mRepetitions;       // Number of repetitions of the execution

      long double                mOptimumFitness;

      int                   mFunction;

      int                   mCMAESMaxPopSize;
      std::string           mCMAESInitials;
      std::string           mBBOBOutputPath;
      std::string           mBBOBConfigName;

      long double                mResetPerc;         // Percentage of the population that will be re-initialized after a restart
      long double                mRerunPerc;         // Percentage of the population that will be re-initialized after in multiple runs

      bool                  mPrintBestSol;      // Print or not the best solution at the end of the execution

      // Convergence options
      convergencecriterion  mConvCrit;          // Convergence criterion

      int                   mIters;             // Iteration limit
      long double                mConv;              // Convergence limit
      int                   mConvIters;         // Iterations used to calculate the convergence
      unsigned long         mEvals;             // Evaluations limit
      int                   mEvalRouteCalls;    // Number of calls to evalRoute (Routing Algs only)
      long double                maxExecSeconds;     // Maximum number of seconds of execution
      long double                mPrecission;        // Precission required to consider convergence
      bool                  mUsePrecission;     // Use the precission value to stop the algorithm?

      // Parallel model options
      int                   mPopSize;           // Population size
      topology              mTop;               // Topology
      modeltype             mModel;             // Model
      int                   mMF;                // Migration frequency
      int                   mMigPopSize;        // Emigrant population size

      emigrants_selec_t     mEmigrantsSelector; // The Emigrants Selector Type

      int                   mTopDegree;         // Topology Degree

      neighbor_conds        mNeighbor_cond;     // Neighbor condition
      int                   mMinNeighborsNum;   // Minimum number of neighbors for each medoid
      int                   mMaxNeighborsNum;   // Maximum number of neighbors for each medoid

      unsigned int          clustering_period_;

      centroids_gen_t       mCentGen;           // Centroids initializer

      islands_init_t        mIslandsInit;       // Islands Initializer

      bool                  mUseCentroids;      // Use centroids initialization?
      long double                mCentProb;          // Centroid mutation probability
      long double                mMinCentDist;       // Minimum distance in centroid generation
      int                   mMaxCentIters;      // Maximum iterations in centroid generation

      //Distributed routing algorithm options
      DistRoutingAlg::emigrantRouteType emigrant_route_t_;

      // Genetic algorithms options
      long double                mMutProb;           // Mutation probability
      long double                mCrossProb;         // Crossover probability
      selector_t            selector_;          // The Selection Scheme
      int                   selector_degree_;

      // Differential Evolution options
      long double                DE_F_factor;        // The mutator factor of the DE algorithm
      long double                DE_CR_factor;       // The recombinator probaility of the DE algorithm
      GAGenome::DECrossover CrossOp;            // crossover option

      // Evolution Strategies options
      unsigned              ES_mu;              // Mu: parents population size
      unsigned              ES_ro;              // Ro: mixing degree (number of parents for recombination)
      unsigned              ES_lambda;          // Lambda: offspring size

      // VNS Options
      long double vnstinitratio_;
      long double vnstinitprob_;
      long double vnslsratio1_;
      long double vnslsprob_;
      long double vnslsratio2_;
      int    vnsmaxevalsshaker_;
      bool   vnsuseallshakersineachstep_;
      int    vnsmaxshakersres_;
      bool   vnsupdatepenalizations_;

      // MTS options
      long double mts_adjustFailed;
      long double mts_adjustMin;
      long double mts_moveLeft;
      long double mts_moveRight;

      // MTS Reduced options
      long double mtsred_searchProb;
      long double mtsred_minProb;

      // Solis-Wets options
      int    sw_maxSuccess;
      int    sw_maxFailed;
      long double sw_adjustSuccess;
      long double sw_adjustFailed;
      long double sw_delta;

      // Adaptive DE options
      long double adapde_Fl;
      long double adapde_Fu;
      long double adapde_tauF;
      long double adapde_tauCR;

      // STSDE options
      long double stsde_prob;

      // MOS options
      char*                 techsRepository;
      vector<string>        techsToUse;
      long double                minPart;            // MOS min participation
      ParticipationFunction partFunction;       // MOS participation function
      evolutionType         evolutiveApproach;  // MOS evolutive approach
      mos_genome_printing_t mosGenomePrinting;  // Type of printing for MOS genome
      long double                mMOSQLPmax;
      long double                mMOSRLAlfa;
      long double                mMOSRLBeta;
      long double                mMOSRLGamma;
      long double                mMOSRLDelta;
      long double                mMOSRLDeltaW;
      long double                mMOSRLDeltaL;
      mosrl_policy          mMOSRLPolicy;
      genealogy             geneal;
      bool                  mPrintGenealogy;    // Should we print the genealogy at the end?
      std::string           mQualityMeasure;    // Quality measure used for dynamic adjustment of participation
      bool                  mWeightedAutonomic; // Weighted or arithmetic average in MOS autonomic
      long double                mBonus;
      std::string           mHetMOSCfgFile;     // Config file for the heterogeneous MOS
      unsigned              mSharedEvals;       // Number of shared evaluations in a MOS HRH algorithm
      long double                mSharedEvalsPercent;// Percentage of shared evaluations in a MOS HRH algorithm
      long double                mSharedEvalsFactor; // Factor to compute the number of shared evaluations in a MOS HRH algorithm (factor * popSize)

      // Debug options
      unsigned              seed;               // the seed of the random method
      int                   loop;

      // Problem structure
      char*                 mProblemFile;       // Proble file
      GAProblemStruct*      mProblemStruct;     // Problem struct

      // Communication manager
      CommManager&           comm_manager_;

      unsigned mRestarts;
};

#endif
