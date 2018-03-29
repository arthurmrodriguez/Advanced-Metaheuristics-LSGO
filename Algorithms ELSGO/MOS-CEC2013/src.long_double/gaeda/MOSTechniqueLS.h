/**
 * @file
 * @brief MOSTechniqueLS class hdr.
 *
 */

#ifndef MOSTechniqueLS_H
#define MOSTechniqueLS_H

#include "MOSTechnique.h"
#include "GARealOps.h"

/**
 * @brief Class for MTS based LS techniques
 */
class MOSTechniqueLS : public MOSTechnique {

 public:

  GADefineIdentity("MOSTechniqueLS", GAID::TechniqueLS);

  // Constructor
  MOSTechniqueLS(techIdType id, std::string description, LocalSearch ls,
                 GAGenome::Comparator comparator, GAGenome::Initializer init,
                 GAGenome::Evaluator evaluator, encodingType encoding,
                 GAGenome* genomeBase, GASelectionScheme* selector);

  ~MOSTechniqueLS();

  // Produces a new offspring from a population
  void offspring (GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset);

  unsigned offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged);

 protected:

  // Atributos
  LocalSearch     _ls;              // LS method
  GAGenome::Comparator _comparator; // Comparison operator

  long double _SR;
  MOSGenome* _currentGenome;

};

#endif
