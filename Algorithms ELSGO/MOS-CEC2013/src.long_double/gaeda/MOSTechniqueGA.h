/**
 * @file
 * @brief MOSTechniqueGA class hdr.
 *
 */

#ifndef MOSTechniqueGA_H
#define MOSTechniqueGA_H

#include "MOSTechnique.h"

/**
 * @brief Clase que representa las técnicas de tipo GA
 */
class MOSTechniqueGA : public MOSTechnique {

 public:

  GADefineIdentity("MOSTechniqueGA", GAID::TechniqueGA);

  /* Constructor */
  MOSTechniqueGA(techIdType id, std::string description,
		   GAGenome::Mutator mutator, GAGenome::SexualCrossover crossover,
		   GAGenome::Comparator comparator, GAGenome::Initializer init,
		   GAGenome::Evaluator evaluator,
		   long double crossProb, long double mutProb,
		   encodingType encoding, GAGenome* genomeBase,
                   GASelectionScheme* selector);

  /* Produces a new offspring from a population */
  void offspring (GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset);

  /* Produces two children from two parents */
  void offspring(MOSGenome& dad, MOSGenome& mom, GAPopulation* destPop, unsigned size, unsigned offset);
  void offspring(MOSGenome& dad, MOSGenome& mom, GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset);

  unsigned offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged);

 protected:

  void offspring_internal(MOSGenome& dad, MOSGenome& mom, GAPopulation* destPop, unsigned offset, unsigned nChilds);

  /* Select a genoma from de population */
  MOSGenome* selectParent(GAPopulation* pop);

  // Atributos
  GAGenome::Mutator         _mutator;       // Operador de mutación
  GAGenome::SexualCrossover _crossover;     // Operador de cruce
  GAGenome::Comparator      _comparator;    // Operador de comparación
  long double                    _crossProb;     // Probabilidad de cruce
  long double                    _mutProb;       // Probabilidad de mutación
};

#endif
