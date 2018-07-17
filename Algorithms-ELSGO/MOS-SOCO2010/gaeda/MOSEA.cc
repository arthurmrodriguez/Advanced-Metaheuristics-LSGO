/**
 * @file
 * @brief MOSEA class impl.
 *
 * Implementation of the MOSEA class
 */

#include "MOSEA.h"

#include "garandom.h"
#include "GAEDAConfig.h"
#include "GAPopulation.h"
#include "GARealOps.h"
#include "MOSGenomeFactory.h"
#include "MOSTechnique.h"
#include "MOSTechniqueSet.h"
#include "MOSParticipationFunction.h"
#include "MOSQualityFunction.h"

#include <dlfcn.h>
#include <algorithm>

/**
 * Constructors
 */
MOSEA::MOSEA (const GAGenome& genome, MOSParticipation* part, MOSQuality* qual)
  : GAGeneticAlgorithm(genome),
    _partFunction(part),
    _qualFunction(qual)
{}

MOSEA::MOSEA (const GAPopulation& population, MOSParticipation* part, MOSQuality* qual)
  : GAGeneticAlgorithm(population),
    _partFunction(part),
    _qualFunction(qual)
{}

MOSEA::MOSEA (const GAGeneticAlgorithm& alg, MOSParticipation* part, MOSQuality* qual)
  : GAGeneticAlgorithm(alg),
    _partFunction(part),
    _qualFunction(qual)
{}

/**
 * Destructor
 */
MOSEA::~MOSEA() {
  if (_partFunction)
    delete _partFunction;
  if(_qualFunction)
    delete _qualFunction;
  return;
}


/**
 * Initialization of the algorithm
 * @param seed Random seed
 */
void MOSEA::initialize (unsigned seed) {

  GARandomSeed(seed);

  _partFunction->setBaseQual(*this);

  pop->initialize();
  pop->evaluate(gaTrue);
  pop->scale();

  stats.scoreFrequency(gaDefScoreFrequency2);
  stats.reset(*pop);

  printStats("Initial Stats");

  return;
}


/**
 * Evolves a generation
 */
void MOSEA::step () {

  unsigned sharedEvals = GAEDAConfig::handle()->getSharedEvals();

  MOSTechniqueSet& techniqueSet = *MOSTechniqueSet::handle();

  // Compute maximum number of evaluations allowed for this offspring
  unsigned maxEvals     = nevals - stats.indEvals();
  unsigned evalsToShare = sharedEvals < maxEvals ? sharedEvals : maxEvals;

  // Update participation ratios
  _partFunction->update(*this);

  MOSTechniqueSet& techSet = *MOSTechniqueSet::handle();
  unsigned nTechs = techSet.nTechniques();

  long double   ratiosAcum = techSet.sumPartRatios();
  unsigned totalEvals = 0;
  unsigned evalsPerTech[nTechs];

  MOSTechniqueSet::MOSTechniqueSetIterator it;
  unsigned i;

  // Distribute the iterations according to their participation ratio
  for (it = techSet.begin(), i = 0; it != techSet.end(); it++, i++) {
    evalsPerTech[i] = (unsigned) (techSet.getPartRatio(it->first) * evalsToShare / ratiosAcum);
    totalEvals += evalsPerTech[i];
  }

  // Distribute the remaining iterations in the case that the division is not exact
  while (totalEvals < evalsToShare) {
    (*(std::min_element(evalsPerTech, evalsPerTech + nTechs)))++;
    totalEvals++;
  }

  if (stats.current(GAStatistics::Maximum) == stats.current(GAStatistics::Minimum)) {
    for (unsigned i = 1; i < pop->size(); i++) {
      pop->individual(i).initialize();
      pop->individual(i).evaluate();
    }
  }


  bool converged = false;
  for (it = techSet.begin(), i = 0; it != techSet.end() && !converged; it++, i++) {
    long double quality = it->second->evolve(pop, evalsPerTech[i], _qualFunction, converged);
    techSet.setTechQuality(it->first, quality);
  }

  // Add stats from each technique to the global stats object
  stats.numsel  = techniqueSet.getAllSelections();
  stats.numcro  = techniqueSet.getAllCrossovers();
  stats.nummut  = techniqueSet.getAllMutations();
  stats.numrep  = techniqueSet.getAllReplacements();
  stats.numeval = techniqueSet.getAllEvals();

  // Update the statistics one generation
  stats.update(*pop);

  // Print stats of this generation
  printStats("End of Step Stats");

  return;

}
