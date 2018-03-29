/**
 * @file
 * @brief MOSTechniqueES class hdr.
 *
 */

#ifndef MOSTechniqueES_H
#define MOSTechniqueES_H


#include "MOSTechnique.h"

/**
 * @brief Class for Evolution Strategies Techniques
 */

class MOSTechniqueES : public MOSTechnique {

 public:

  GADefineIdentity("MOSTechniqueES", GAID::TechniqueES);

  /* Constructor */
  MOSTechniqueES(techIdType id, std::string desc, GAGenome::ESMutator mut, GAGenome::ESMutator mut_upd,
                 GAGenome::ESCrossover cx, GAGenome::ESCrossover cx_upd,
                 GAGenome::Comparator cmp, GAGenome::Initializer init, GAGenome::Evaluator eval,
                 unsigned mu, unsigned ro, unsigned lambda, encodingType encoding, GAGenome* genomeBase,
                 GASelectionScheme* slct);

  /* Produces a new offspring from a population */
  void offspring (GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset);

  /* Produces two children from two parents */
  void offspring(MOSGenome& dad, MOSGenome& mom, GAPopulation* destPop, unsigned size, unsigned offset);
  void offspring(MOSGenome& dad, MOSGenome& mom, GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset);

  unsigned offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged) {return 0;}

  bool updateStrategies(const std::vector<GAGenome*>& parents, GAGenome* child) const;

 protected:

  void offspring_internal(const std::vector<MOSGenome*> parents, GAPopulation* destPop, unsigned offset);

  /* Select a genome from thee population */
  MOSGenome* selectParent(GAPopulation* pop);

  // Attributes
  GAGenome::ESMutator       _mutator;       // Mutation operator
  GAGenome::ESMutator       _mutator_upd;   // Mutation operator for update
  GAGenome::ESCrossover     _crossover;     // Crossover operator
  GAGenome::ESCrossover     _crossover_upd; // Crossover operator for update
  GAGenome::Comparator      _comparator;    // Compration operator
  unsigned                  _mu;            // Mu (parents population size)
  unsigned                  _ro;            // Ro (mixing ratio)
  unsigned                  _lambda;        // Lambda (children population size)
  long double                    _tau;           // Tau   (Factor for mutation)
  long double                    _tau0;          // Tau_0 (Factor for mutation)

  bool                      _message_printed;
};

#endif
