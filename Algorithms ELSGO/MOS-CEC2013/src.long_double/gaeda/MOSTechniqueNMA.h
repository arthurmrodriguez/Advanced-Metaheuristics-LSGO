/**
 * @file
 * @brief MOSTechniqueNMA class hdr.
 *
 */

#ifndef MOSTechniqueNMA_H
#define MOSTechniqueNMA_H

#include "MOSTechnique.h"
#include "GARealOps.h"
#include "genomes/GA1DArrayGenome.h"

/**
 * @brief Class for MTS based LS techniques
 */
class MOSTechniqueNMA : public MOSTechnique {

 public:

  typedef enum {Reflection, Expansion, ContractionInOK, ContractionOutOK, Shrink} Movements;

  GADefineIdentity("MOSTechniqueNMA", GAID::TechniqueNMA);

  // Constructor
  MOSTechniqueNMA(techIdType id, std::string description,
                  GAGenome::Comparator comparator, GAGenome::Initializer init,
                  GAGenome::Evaluator evaluator, encodingType encoding,
                  GAGenome* genomeBase, GASelectionScheme* selector);

  // Produces a new offspring from a population
  void offspring (GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset);

  unsigned offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged);

 protected:

  MOSGenome* centroidOf        (const GAPopulation& pop) const;
  MOSGenome* reflection        (const MOSGenome& centroid, const MOSGenome& worst_ind, long double rho) const;
  void        expansion        (const MOSGenome& centroid, const MOSGenome& reflx_ind, MOSGenome& worst_ind, long double chi, stringstream& msg) const;
  MOSGenome* outsideContraction(const MOSGenome& centroid, const MOSGenome& reflx_ind, long double gamma) const;
  MOSGenome* insideContraction (const MOSGenome& centroid, const MOSGenome& worst_ind, long double gamma) const;
  void       shrinkPop         (const GAPopulation& pop, unsigned& used_evals, unsigned max_evals, long double sigma) const;
  void       contraction       (const MOSGenome& centroid, const MOSGenome& reflx_ind, MOSGenome& worst_ind, const GAPopulation& pop, unsigned& used_evals, unsigned max_evals, unsigned& branch, long double gamma, long double sigma, stringstream& msg) const;

  bool       NMA               (GAPopulation& pop, /* out */ unsigned& used_evals, unsigned max_evals, unsigned& branch, long double rho=1.0, long double chi=2.0, long double gamma=0.5, long double sigma=0.5);

  // Atributos
  GAGenome::Comparator _comparator; // Comparison operator

};

#endif
