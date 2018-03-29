/**
 * @file
 * @brief MOSTechniqueMemeticGA class hdr.
 *
 */

#ifndef MOSTechniqueMemeticGA_H
#define MOSTechniqueMemeticGA_H

#include "MOSTechnique.h"
#include "GARealOps.h"

/**
 * @brief Clase que representa las técnicas de tipo GA
 */
class MOSTechniqueMemeticGA : public MOSTechnique {

 public:

  GADefineIdentity("MOSTechniqueMemeticGA", GAID::TechniqueGA);

  /* Constructor */
  MOSTechniqueMemeticGA(techIdType id, std::string description,
		   GAGenome::Mutator mutator, GAGenome::SexualCrossover crossover,
		   GAGenome::Comparator comparator, GAGenome::Initializer init,
		   GAGenome::Evaluator evaluator,
		   double crossProb, double mutProb,
		   encodingType encoding, GAGenome* genomeBase,
       GASelectionScheme* selector, std::vector<LSType>& lss);

  unsigned offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged);

  bool restartRequired();
  bool restartInnerData(GAPopulation* pop);

 protected:

  void offspring_internal(MOSGenome& dad, MOSGenome& mom, GAPopulation* destPop, unsigned offset, unsigned nChilds);
  void applyLS(MOSGenome& genome, GAPopulation* pop, unsigned k);

  /* Select a genoma from de population */
  MOSGenome* selectParent(GAPopulation* pop);

  // Atributos
  GAGenome::Mutator         _mutator;       // Operador de mutación
  GAGenome::SexualCrossover _crossover;     // Operador de cruce
  GAGenome::Comparator      _comparator;    // Operador de comparación
  double                    _crossProb;     // Probabilidad de cruce
  double                    _mutProb;       // Probabilidad de mutación
  std::vector<LSType> _lss;
  unsigned _gensCount;
};

#endif
