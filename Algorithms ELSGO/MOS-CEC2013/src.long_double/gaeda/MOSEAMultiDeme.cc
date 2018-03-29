/**
 * @file
 * @brief MOSEAMultiDeme class impl.
 *
 * Implementation of the MOSEAMultiDeme class
 */

#include "MOSEAMultiDeme.h"

#include "garandom.h"
#include "GAEDAConfig.h"
#include "GAPopulation.h"
#include "GAGenealogy.h"
#include "GAGenealogyMemory.h"
#include "Recombinator.h"
#include "MOSConversion.h"
#include "MOSTechnique.h"
#include "MOSTechniqueSet.h"
#include "NSC.h"
#include "genomes/MOSGenome.h"

/**
 * Constructors
 */
MOSEAMultiDeme::MOSEAMultiDeme (const std::vector<MOSGenome*>& genomes, const long double elitismPercent, const long double minPart, const evolutionType evolutiveApproach)
  : MOSEA (*genomes[0], elitismPercent, minPart, evolutiveApproach),
    _pops     (MOSTechniqueSet::handle()->nTechniques(), (GAPopulation*)0),
    _auxPops  (MOSTechniqueSet::handle()->nTechniques(), (GAPopulation*)0)
{

  unsigned nTechs  = _techniqueSet->nTechniques ();
  unsigned popSize = GAEDAConfig::handle()->getPopSize();

  std::vector<unsigned> sizes (nTechs, 0);

  unsigned sizePerTech = popSize / nTechs;
  int extraIndivs = popSize - (sizePerTech * nTechs);

  MOSTechniqueSet::MOSTechniqueSetIterator it;
  for (it=_techniqueSet->begin(); it!=_techniqueSet->end(); it++, extraIndivs--) {
    sizes[it->first] = sizePerTech;
    if (extraIndivs > 0)
      sizes[it->first]++;
  }

  for (unsigned i = 0; i < _pops.size(); i++) {
    _pops[i]    = new GAPopulation (*genomes[i], sizes[i]);
    _auxPops[i] = new GAPopulation (*genomes[i], sizes[i]);
  }

}


/**
 * Destructor
 */
MOSEAMultiDeme::~MOSEAMultiDeme() {

  return;

}


/**
 * Initialization of the algorithm
 * @param seed Random seed
 */
void MOSEAMultiDeme::initialize () {
  qualityMeasure = GAEDAConfig::handle()->getQualityMeasure();

  _techniqueSet->initPartRatios();

  unsigned offset = 0;
  for (unsigned i = 0; i < _pops.size(); i++) {
    _pops[i]->initialize();
    _pops[i]->evaluate(gaTrue);
    _pops[i]->scale();

    for (unsigned j = 0; j < _pops[i]->size(); j++) {
      delete pop->replace(_pops[i]->individual(j).clone(), offset);
      offset++;
    }
  }

  pop->scale();
  pop->sort();

  stats.scoreFrequency(gaDefScoreFrequency2);
  stats.reset(*pop);

  // BEGIN: Genealogy
  if (GAGenealogy::handle()) {

    // Add first individuals to the genealogy
    for (unsigned i=0; i < pop->size(); i++)
      GAGenealogy::handle()->addNode(pop->individual(i));

    if (GAGenealogy::handle()->isGenealogyMemory()) {

      // Put the ranking of each genome in the current generation in the genealogy
      GAGenealogyMemory* genealmem = dynamic_cast<GAGenealogyMemory*> (GAGenealogy::handle());
      genealmem->updateRankings(pop);

      if (qualityMeasure == "NSC") { // Prepare the NSC object
        _nsc = new NSC (pop->size());
        _nsc->addPopulation (pop);
      }

    }

  }
  // END: Genealogy

  printStats("Initial Stats");

  return;
}


/**
 * Evolves a generation
 */
void MOSEAMultiDeme::step () {

  // BEGIN: Genealogy
  GAGenealogyMemory *genealmem = NULL;

  // New step in the Genealogy
  if (GAGenealogy::handle()) {
    if (GAGenealogy::handle()->isGenealogyMemory())
      genealmem = dynamic_cast<GAGenealogyMemory*> (GAGenealogy::handle());
    GAGenealogy::handle()->newGeneration();
  }
  // END: Genealogy

  // Update participation ratios
  updatePartRatios();

  // Create new offspring
  _techniqueSet->offspringAll(_pops, _auxPops);

  // BEGIN: Genealogy
  // Take current population to prepare the NSC data
  if (GAGenealogy::handle())
    if (GAGenealogy::handle()->isGenealogyMemory() && qualityMeasure == "NSC")
      _nsc->prepareNSCData(pop);
  // END: Genealogy

  // Copy individuals from aux populations to _auxPop
  // This is necessary to be able to evaluate the quality of
  // the whole offspring
  unsigned offset = 0;

  for (unsigned i = 0; i < _auxPops.size(); i++) {
    _auxPops[i]->evaluate();
    _auxPops[i]->scale();

    for (unsigned j = 0; j < _auxPops[i]->size(); j++) {
      delete _auxPop->replace(_auxPops[i]->individual(j).clone(), offset);
      offset++;
    }
  }

  // Compute quality of offspring population
  evalTechQuality(_auxPop);

  _auxPop->evaluate();
  _auxPop->scale();

  // Apply elitism to old and new populations
  MOSTechniqueSet::MOSTechniqueSetIterator it;

  for (it=_techniqueSet->begin(); it!=_techniqueSet->end(); it++) {
    recombinator_->recombine (*_pops[it->first], *_auxPops[it->first]);
    std::swap (_pops[it->first], _auxPops[it->first]);
    _pops[it->first]->sort();
  }

  // Exchange individuals every N generations
  std::vector< std::vector<GAGenome*> > genomes (_pops.size());
  std::vector<unsigned> origPopSizes (_pops.size());

  // First, store individuals to migrate and old populations sizes
  if (this->generation() % 20 == 0) {

    for (unsigned i = 0; i < _pops.size(); i++) {

      origPopSizes [i] = _pops[i]->size();

      for (unsigned j = 0; j < 0.10 * _pops[i]->size(); j++)
        genomes[i].push_back(&_pops[i]->individual(j));

    }

    // Second, copy individuals from other populations and restore
    // old population size
    for (unsigned i = 0; i < _pops.size(); i++) {

      for (unsigned j = 0; j < _pops.size(); j++)
        if (i != j)
          for (unsigned k = 0; k < genomes[j].size(); k++)
            _pops[i]->add(*genomes[j][k]);

      _pops[i]->size(origPopSizes[i]);

    }

  }

  // Create pop population from multiple populations
  offset = 0;
  for (unsigned i = 0; i < _pops.size(); i++) {
    for (unsigned j = 0; j < _pops[i]->size(); j++) {
      delete pop->replace(_pops[i]->individual(j).clone(), offset);
      offset++;
    }
  }

  // FIXME: This should not be necessary
  // Evalua el fitness de los genomas de la poblacion
  pop->evaluate();

  // Escala el valor del fitness, para fijarlo entre 0 y 1
  pop->scale();

  // Add stats from each technique to the global stats object
  stats.numsel  = _techniqueSet->getAllSelections();
  stats.numcro  = _techniqueSet->getAllCrossovers();
  stats.nummut  = _techniqueSet->getAllMutations();
  stats.numrep  = _techniqueSet->getAllReplacements();
  stats.numeval = _techniqueSet->getAllEvals();

  // Update the statistics one generation
  stats.update(*pop);

  // BEGIN: Genealogy
  //Put ranking of each genoma for that generation in the genealogy
  if (GAGenealogy::handle())
    if (GAGenealogy::handle()->isGenealogyMemory())
      genealmem->updateRankings(pop);
  // END: Genealogy

  // Print stats of this generation
  printStats("End of Step Stats");

  return;

}
