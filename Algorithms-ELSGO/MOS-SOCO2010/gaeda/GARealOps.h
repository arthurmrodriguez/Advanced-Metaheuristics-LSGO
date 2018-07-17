/**
  * @file
  * @brief Genetic operators for integer problems.
  *
  */

#ifndef GARealOps_H
#define GARealOps_H

class GAGenome;

// Initializers
extern "C" void RealUniformInitializer (GAGenome& g);

// Comparators
extern "C" long double RealEuclideanComparator (const GAGenome& g1, const GAGenome& g2);

// GA Crossovers
extern "C" int RealBlendCrossover       (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2);

// GA Mutators
extern "C" int RealGaussianMutator (GAGenome& g, long double pmut);

// DE Crossovers
extern "C" int RealExponentialCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, const long double prob);

// MTS Local Searches
const unsigned BONUS1 = 10;
const unsigned BONUS2 = 1;

typedef long double (*LocalSearch) (GAGenome& Xk, long double& fitBest, unsigned& evals, long double& fit_inc_acum, unsigned& improvements, long double& SR, unsigned maxIters);

extern "C" long double MTS_LS1 (GAGenome& Xk, long double& fitBest, unsigned& evals, long double& fit_inc_acum, unsigned& improvements, long double& SR, unsigned maxIters);

#endif
