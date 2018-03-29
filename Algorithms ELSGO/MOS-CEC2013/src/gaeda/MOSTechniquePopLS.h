/**
 * @file
 * @brief MOSTechniquePopLS class hdr.
 *
 */

#ifndef MOSTechniquePopLS_H
#define MOSTechniquePopLS_H

#include "MOSTechnique.h"
#include "GARealOps.h"

/**
 * @brief Class for MTS based LS techniques
 */
class MOSTechniquePopLS : public MOSTechnique {

 public:

  GADefineIdentity("MOSTechniquePopLS", GAID::TechniquePopLS);

  // Constructor
  MOSTechniquePopLS(techIdType id, std::string description, PopLocalSearch ls,
                   GAGenome::Comparator comparator, GAGenome::Initializer init,
                   GAGenome::Evaluator evaluator, encodingType encoding,
                   GAGenome* genomeBase, GASelectionScheme* selector);

  void offspring (GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset) {std::cout << "Not implemented for this technique." << std::endl; exit(-1);}

  unsigned offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged);

 protected:

  // Atributos
  PopLocalSearch     _ls;           // LS method
  GAGenome::Comparator _comparator; // Comparison operator

  MOSGenome* _currentGenome;

};

#endif
