/**
  * @file
  * @brief Genetic operators for integer problems.
  *
  */

#ifndef GARealOps_H
#define GARealOps_H

#include <vector>

#include "genomes/GAGenome.h"

class GAPopulation;

// Initializers

extern "C" void RealUniformInitializer (GAGenome& g);
extern "C" void RealOrderedInitializer (GAGenome& g);

// Comparators

extern "C" double RealElementComparator             (const GAGenome& g1, const GAGenome& g2);
extern "C" double RealEuclideanComparator           (const GAGenome& g1, const GAGenome& g2);
extern "C" double RealEuclideanLinearizedComparator (const GAGenome& g1, const GAGenome& g2);
extern "C" double RealChebyshevComparator           (const GAGenome& g1, const GAGenome& g2);
extern "C" double RealManhattanComparator           (const GAGenome& g1, const GAGenome& g2);

// Crossovers
extern "C" int RealOnePointCrossover     (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2);
extern "C" int RealTwoPointsCrossover    (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2);
extern "C" int RealEvenOddCrossover      (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2);
extern "C" int RealPartialMatchCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2);
extern "C" int RealUniformCrossover      (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2);
extern "C" int RealOrderCrossover        (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2);
extern "C" int RealCycleCrossover        (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2);

extern "C" int RealAlternativeOrderCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2);
extern "C" int RealAlternativeCycleCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2);

extern "C" int RealBlendCrossover       (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2);
extern "C" int RealArithmethicCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2);

// DE Crossovers
extern "C" int RealExponentialCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, const double prob);
extern "C" int RealBinomialCrossover    (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, const double prob);

// Mutators
extern "C" int RealSwapMutator             (GAGenome& g, double pmut);
extern "C" int RealFlipMutator             (GAGenome& g, double pmut);
extern "C" int RealRepeatedExchangeMutator (GAGenome& g, double pmut);
extern "C" int RealSimpleInversionMutator  (GAGenome& g, double pmut);

extern "C" int RealUniformMutator  (GAGenome& g, double pmut);
extern "C" int RealGaussianMutator (GAGenome& g, double pmut);

// ES Crossovers
extern "C" int RealIntermediateCrossover           (const std::vector<GAGenome*> parents, GAGenome* child);
extern "C" int RealIntermediateCrossoverUpdateOnly (const std::vector<GAGenome*> parents, GAGenome* child);
extern "C" int RealDominantCrossover               (const std::vector<GAGenome*> parents, GAGenome* child);
extern "C" int RealDominantCrossoverUpdateOnly     (const std::vector<GAGenome*> parents, GAGenome* child);

// ES Mutators
extern "C" int RealIsotropicMutator              (GAGenome& g);
extern "C" int RealIsotropicMutatorUpdateOnly    (GAGenome& g);
extern "C" int RealNonIsotropicMutator           (GAGenome& g);
extern "C" int RealNonIsotropicMutatorUpdateOnly (GAGenome& g);

// ES Initializers
extern "C" void RealESUniformInitializer (GAGenome& g);

// MTS Local Searches
const unsigned BONUS1 = 10;
const unsigned BONUS2 = 1;

typedef double (*LocalSearch) (GAGenome& Xk, double& fitBest, unsigned& evals, double& fit_inc_acum, unsigned& improvements, double& SR, unsigned maxIters);

extern "C" double MTS_LS1 (GAGenome& Xk, double& fitBest, unsigned& evals, double& fit_inc_acum, unsigned& improvements, double& SR, unsigned maxIters);
extern "C" double MTS_LS2 (GAGenome& Xk, double& fitBest, unsigned& evals, double& fit_inc_acum, unsigned maxIters);
extern "C" double MTS_LS3 (GAGenome& Xk, double& fitBest, unsigned& evals, double& fit_inc_acum, unsigned maxIters);

extern "C" double MTS_LS1_Reduced_Dims (GAGenome& Xk, double& fitBest, unsigned& evals, double& fit_inc_acum, unsigned& improvements, unsigned maxIters);

extern "C" double RandomGreedyLS (GAGenome& Xk, double& fitBest, unsigned& evals, double& fit_inc_acum, unsigned& improvements, double& SR, unsigned maxIters);

typedef bool (*PopLocalSearch) (GAPopulation& pop, unsigned& used_evals, double& fit_inc_acum, unsigned& improvements, unsigned max_evals);

extern "C" bool DiffsLS (GAPopulation& pop, unsigned& used_evals, double& fit_inc_acum, unsigned& improvements, unsigned max_evals);

class LSType {
  public:
    GAGenome::Mutator ls;
    unsigned strategy;
    unsigned freq;

    typedef enum {BEST = 1, RANDOM = 2, ALL = 3} LSTypes;
};



#endif
