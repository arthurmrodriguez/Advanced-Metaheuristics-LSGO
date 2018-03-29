/**
 * @file
 * @brief MOSTechniqueLSObj class hdr.
 *
 */

#ifndef MOSTechniqueLSObj_H
#define MOSTechniqueLSObj_H

#include "MOSTechnique.h"
#include "GARealOps.h"
#include "solis.h"

class SolisParams;

/**
 * @brief Class for MTS based LS techniques
 */
class MOSTechniqueLSObj : public MOSTechnique {

 public:

  GADefineIdentity("MOSTechniqueLSObj", GAID::TechniqueLSObj);

  // Constructor
  MOSTechniqueLSObj(techIdType id, std::string description, std::string ls,
                 GAGenome::Comparator comparator, GAGenome::Initializer init,
                 GAGenome::Evaluator evaluator, encodingType encoding,
                 GAGenome* genomeBase, GASelectionScheme* selector);

  ~MOSTechniqueLSObj();

  // Produces a new offspring from a population
  void offspring (GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset);

  unsigned offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged);

 protected:

  // Atributos
  GAGenome::Comparator _comparator; // Comparison operator

  MOSGenome* _currentGenome;

  SolisWets _sw;
  SolisParams* _swparams;

};

#endif
