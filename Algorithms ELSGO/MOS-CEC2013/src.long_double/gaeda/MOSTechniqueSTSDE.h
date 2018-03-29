#ifndef MOSTechniqueSTSDE_H
#define MOSTechniqueSTSDE_H

#include "MOSTechnique.h"

/**
 * @brief Clase que representa las tecnicas de tipo STSDE
 */
class MOSTechniqueSTSDE : public MOSTechnique {

 public:

  GADefineIdentity("MOSTechniqueSTSDE", GAID::TechniqueSTSDE);

  /* Constructor */
  MOSTechniqueSTSDE(techIdType id, std::string description,
                    GAGenome::Initializer init, GAGenome::Evaluator evaluator,
                    encodingType encoding, GAGenome* genomeBase,
                    GAGenome::DECrossover crossover,
                    long double F, long double CR, GASelectionScheme* selector);

  /* Produces a new offspring from a population */
  void offspring (GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset);

  /* Produces two children from two parents */
  void offspring(MOSGenome& dad, MOSGenome& mom, GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset) {}

  unsigned offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged);

 protected:

  void updateRanges(GAPopulation* pop);
  void genSTSPop   (GAPopulation* pop, GAPopulation* stsPop);

  virtual void offspring_internal(MOSGenome& genomeM, MOSGenome& newGenomeM, MOSGenome& genome1M, MOSGenome& genome2M, MOSGenome& genome3M);
  bool selectParents(GAPopulation* oldPop, MOSGenome* x_i, MOSGenome*& t1, MOSGenome*& t2, MOSGenome*& t3);

  // Atributos
  GAGenome::DECrossover _crossover;     // Operador de cruce
  long double _F;   // Factor FI de mutacion del DE
  long double _CR;  // Factor de la probabilidad de cruce del DE

  long double*       min_range_;
  long double*       max_range_;

  unsigned _curTarget;

  unsigned attempts;

};

#endif
