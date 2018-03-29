/**
  * @file
  * @brief Genetic operators for integer problems.
  *
  */

#ifndef GAIntOps_H
#define GAIntOps_H

#include "genomes/GAGenome.h"
#include "GACommonOps.h"

// Initializers
extern "C" void IntegerUniformInitializer (GAGenome& g) {return UniformInitializer<int> (g);}
extern "C" void IntegerOrderedInitializer (GAGenome& g) {return OrderedInitializer<int> (g);}

// Comparators
extern "C" double IntegerElementComparator (const GAGenome& g1, const GAGenome& g2) {return ElementComparator<int> (g1, g2);}

// Crossovers
extern "C" int IntegerOnePointCrossover     (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return OnePointCrossover<int>     (p1, p2, c1, c2);}
extern "C" int IntegerTwoPointsCrossover    (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return TwoPointsCrossover<int>    (p1, p2, c1, c2);}
extern "C" int IntegerEvenOddCrossover      (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return EvenOddCrossover<int>      (p1, p2, c1, c2);}
extern "C" int IntegerPartialMatchCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return PartialMatchCrossover<int> (p1, p2, c1, c2);}
extern "C" int IntegerUniformCrossover      (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return UniformCrossover<int>      (p1, p2, c1, c2);}
extern "C" int IntegerOrderCrossover        (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return OrderCrossover<int>        (p1, p2, c1, c2);}
extern "C" int IntegerCycleCrossover        (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return CycleCrossover<int>        (p1, p2, c1, c2);}

extern "C" int IntegerAlternativeOrderCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return AlternativeOrderCrossover<int> (p1, p2, c1, c2);}
extern "C" int IntegerAlternativeCycleCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return AlternativeCycleCrossover<int> (p1, p2, c1, c2);}

// DE Crossovers
extern "C" int IntegerExponentialCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, const double prob) {return ExponentialCrossover<int> (p1, p2, c1, prob);}
extern "C" int IntegerBinomialCrossover    (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, const double prob) {return BinomialCrossover<int>    (p1, p2, c1, prob);}

// Mutators
extern "C" int IntegerSwapMutator             (GAGenome& g, double pmut) {return SwapMutator<int>             (g, pmut);}
extern "C" int IntegerFlipMutator             (GAGenome& g, double pmut) {return FlipMutator<int>             (g, pmut);}
extern "C" int IntegerRepeatedExchangeMutator (GAGenome& g, double pmut) {return RepeatedExchangeMutator<int> (g, pmut);}
extern "C" int IntegerSimpleInversionMutator  (GAGenome& g, double pmut) {return SimpleInversionMutator<int>  (g, pmut);}

// ES Crossovers
extern "C" int BinaryDominantCrossover (const std::vector<GAGenome*> parents, GAGenome* child);
extern "C" int BinaryDominantCrossoverUpdateOnly (const std::vector<GAGenome*> parents, GAGenome* child);

// ES Mutators
extern "C" int BinaryFlipMutator (GAGenome& g);
extern "C" int BinaryFlipMutatorUpdateOnly (GAGenome& g);
#endif
