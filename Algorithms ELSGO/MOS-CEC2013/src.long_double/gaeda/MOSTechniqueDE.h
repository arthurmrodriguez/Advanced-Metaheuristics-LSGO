#ifndef MOSTechniqueDE_H
#define MOSTechniqueDE_H

#include "MOSTechnique.h"

/**
 * @brief Clase que representa las técnicas de tipo DE
 */
class MOSTechniqueDE : public MOSTechnique {

 public:

  GADefineIdentity("MOSTechniqueDE", GAID::TechniqueDE);

  /* Constructor */
  MOSTechniqueDE(techIdType id, std::string description,
                 GAGenome::Initializer init, GAGenome::Evaluator evaluator,
                 encodingType encoding, GAGenome* genomeBase,
                 GAGenome::DECrossover crossover,
                 long double F, long double CR, GASelectionScheme* selector);

  /* Produces a new offspring from a population */
  void offspring (GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset);

  /* Produces two children from two parents */
  void offspring(MOSGenome& dad, MOSGenome& mom, GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset);

  unsigned offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged);

 protected:

  virtual void offspring_internal(MOSGenome& genome, MOSGenome& newGenome, MOSGenome& genome1, MOSGenome& genome2, MOSGenome& genome3);

  bool selectParents(GAPopulation* oldPop, MOSGenome* x_i, MOSGenome*& t1, MOSGenome*& t2, MOSGenome*& t3);

  // Atributos
  GAGenome::DECrossover _crossover;     // Operador de cruce
  long double _F;   // Factor FI de mutacion del DE
  long double _CR;  // Factor de la probabilidad de cruce del DE

  unsigned _curTarget;

  unsigned attempts;

};

#endif
