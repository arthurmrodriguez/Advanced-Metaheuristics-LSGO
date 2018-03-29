/**
  * @file
  * @brief Genetic operators for string problems.
  *
  */

#ifndef GAStringOps_H
#define GAStringOps_H

#include "genomes/GAGenome.h"
#include "GACommonOps.h"

// Initializers
extern "C" void StringUniformInitializer (GAGenome& g) {return UniformInitializer<char> (g);}

// Comparators

extern "C" double StringElementComparator (const GAGenome& g1, const GAGenome& g2) {return ElementComparator<char> (g1, g2);};

// Crossovers
extern "C" int StringOnePointCrossover     (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return OnePointCrossover<char>     (p1, p2, c1, c2);}
extern "C" int StringTwoPointsCrossover    (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return TwoPointsCrossover<char>    (p1, p2, c1, c2);}
extern "C" int StringEvenOddCrossover      (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return EvenOddCrossover<char>      (p1, p2, c1, c2);}
extern "C" int StringPartialMatchCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return PartialMatchCrossover<char> (p1, p2, c1, c2);}
extern "C" int StringUniformCrossover      (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return UniformCrossover<char>      (p1, p2, c1, c2);}
extern "C" int StringOrderCrossover        (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return OrderCrossover<char>        (p1, p2, c1, c2);}
extern "C" int StringCycleCrossover        (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return CycleCrossover<char>        (p1, p2, c1, c2);}

extern "C" int StringAlternativeOrderCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return AlternativeOrderCrossover<char> (p1, p2, c1, c2);}
extern "C" int StringAlternativeCycleCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2) {return AlternativeCycleCrossover<char> (p1, p2, c1, c2);}

// Mutators
extern "C" int StringSwapMutator             (GAGenome& g, double pmut) {return SwapMutator<char>             (g, pmut);}
extern "C" int StringFlipMutator             (GAGenome& g, double pmut) {return FlipMutator<char>             (g, pmut);}
extern "C" int StringRepeatedExchangeMutator (GAGenome& g, double pmut) {return RepeatedExchangeMutator<char> (g, pmut);}
extern "C" int StringSimpleInversionMutator  (GAGenome& g, double pmut) {return SimpleInversionMutator<char>  (g, pmut);}

#endif