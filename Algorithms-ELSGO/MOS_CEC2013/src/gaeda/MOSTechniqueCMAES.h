#ifndef MOSTechniqueCMAES_H
#define MOSTechniqueCMAES_H

#include "MOSTechnique.h"

#include "cmaes.h"

/**
 * @brief Clase que representa las técnicas de tipo CMAES
 */
class MOSTechniqueCMAES : public MOSTechnique {

 public:

  GADefineIdentity("MOSTechniqueCMAES", GAID::TechniqueCMAES);

  /* Constructor */
  MOSTechniqueCMAES(techIdType id, std::string description,
                    GAGenome::Initializer init, GAGenome::Evaluator evaluator,
                    encodingType encoding, GAGenome* genomeBase,
                    GASelectionScheme* selector);

  void offspring (GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset) {}
  void offspring(MOSGenome& dad, MOSGenome& mom, GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset) {};

  unsigned offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged);

  bool restartRequired();
  bool restartInnerData(GAPopulation* pop);

 protected:

  double* _fitvals;
  cmaes_t _cmaes_evo;

  double _old_best;
  double _old_worst;

};

#endif
