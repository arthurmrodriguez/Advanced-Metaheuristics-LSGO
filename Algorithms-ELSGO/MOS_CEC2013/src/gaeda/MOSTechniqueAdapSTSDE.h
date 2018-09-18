#ifndef MOSTechniqueAdapSTSDE_H
#define MOSTechniqueAdapSTSDE_H

#include "MOSTechniqueSTSDE.h"

class MOSTechniqueAdapSTSDE : public MOSTechniqueSTSDE {
 public:
  GADefineIdentity("MOSTechniqueAdapSTSDE", GAID::TechniqueAdapSTSDE);

  MOSTechniqueAdapSTSDE(techIdType id, std::string description,
                     GAGenome::Initializer init, GAGenome::Evaluator evaluator,
                     encodingType encoding, GAGenome* genomeBase,
                     GAGenome::DECrossover crossover,
                     double F, double CR, GASelectionScheme* selector);

  void offspring (GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset);

  void offspring(MOSGenome& dad, MOSGenome& mom, GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset);

  unsigned offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged);

 protected:
  virtual void offspring_internal(MOSGenome& genomeM, MOSGenome& newGenomeM, MOSGenome& genome1M, MOSGenome& genome2M, MOSGenome& genome3M);
  void selfAdjustDEParams(const GAPopulation& pop);

  std::vector<double> _Fv;
  std::vector<double> _CRv;
};

#endif
