/**
 * @file
 * @brief MOSEA class impl.
 *
 * Implementation of the MOSEA class
 */

#include "MOSEA.h"

#include "extras/distaux.h"
#include "garandom.h"
#include "NSC.h"
#include "GAEDAConfig.h"
#include "GAPopulation.h"
#include "GAGenealogy.h"
#include "GAGenealogyMemory.h"
#include "Recombinator.h"
#include "MOSConversion.h"
#include "MOSTechnique.h"
#include "MOSTechniqueSet.h"
#include "MOSParticipationFunc.h"
#include "genomes/MOSGenome.h"

/**
 * Constructors
 */
MOSEA::MOSEA (const GAGenome& genome, const double elitismPercent, const double minPart, const evolutionType evolutiveApproach)
  : GAGeneticAlgorithm(genome), _bonus(GAEDAConfig::handle()->getBonus()),
    _auxPop(NULL),
    _evolutiveApproach(evolutiveApproach),
    _partRestoreGen(0), _partFunction(constantPF),
    _PFData(NULL), _minPart(minPart),
    _elitismPercent(elitismPercent),
    _nsc (NULL),
    _actualPopSize(pop->size())
{
  _auxPop = new GAPopulation(*pop);
  _techniqueSet = MOSTechniqueSet::handle();
  _conversion = MOSConversion::handle();
}

MOSEA::MOSEA (const GAPopulation& population, const double elitismPercent, const double minPart, const evolutionType evolutiveApproach)
  : GAGeneticAlgorithm(population), _bonus(GAEDAConfig::handle()->getBonus()),
    _auxPop(NULL),
    _evolutiveApproach(evolutiveApproach),
    _partRestoreGen(0), _partFunction(constantPF),
    _PFData(NULL), _minPart(minPart),
    _elitismPercent(elitismPercent),
    _nsc (NULL),
    _actualPopSize(pop->size())
{
  _auxPop = new GAPopulation(*pop);
  _techniqueSet = MOSTechniqueSet::handle();
  _conversion = MOSConversion::handle();
}

MOSEA::MOSEA (const GAGeneticAlgorithm& alg, const double elitismPercent, const double minPart, const evolutionType evolutiveApproach)
  : GAGeneticAlgorithm(alg), _bonus(GAEDAConfig::handle()->getBonus()),
    _auxPop(NULL),
    _evolutiveApproach(evolutiveApproach),
    _partRestoreGen(0), _partFunction(constantPF),
    _PFData(NULL), _minPart(minPart),
    _elitismPercent(elitismPercent),
    _nsc (NULL),
    _actualPopSize(pop->size())
{
  _auxPop = new GAPopulation(*pop);
  _techniqueSet = MOSTechniqueSet::handle();
  _conversion = MOSConversion::handle();
}

/**
 * Destructor
 */
MOSEA::~MOSEA() {

  // BEGIN: Genealogy
  if (GAGenealogy::handle()) {

    // Decease the population of the last generation
    for (unsigned i=0; i<pop->n; i++)
      GAGenealogy::handle()->deceased(pop->individual(i));

    if (GAGenealogy::handle()->isGenealogyTracer())
      GAGenealogy::handle()->bestGenome(pop->best());

    // Print the hash table
    if (GAGenealogy::handle()->isGenealogyMemory() &&
        GAEDAConfig::handle()->printGenealogy   ()    ) {
      GAGenealogyMemory* genealmem = dynamic_cast<GAGenealogyMemory*> (GAGenealogy::handle());
      genealmem->print();
    }

  }

  if (_nsc)
    delete _nsc;
  // END: Genealogy

  if (_auxPop) {
    delete _auxPop;
  }

  return;

}


/**
 * Initialization of the algorithm
 * @param seed Random seed
 */
void MOSEA::initialize () {

  qualityMeasure = GAEDAConfig::handle()->getQualityMeasure();

  _techniqueSet->initPartRatios();

  pop->initialize();
  pop->evaluate(gaTrue);
  pop->scale();

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
 * Updates the participation ratio of every technique
 */
void MOSEA::updatePartRatios() {
  _partFunction(*this);
  return;
}


/**
 * Evolves a generation
 */
void MOSEA::step () {

  // BEGIN: Genealogy
  GAGenealogyMemory *genealmem = NULL;

  // New step in the Genealogy
  if (GAGenealogy::handle()) {

    if (GAGenealogy::handle()->isGenealogyMemory())
      genealmem = dynamic_cast<GAGenealogyMemory*> (GAGenealogy::handle());

    GAGenealogy::handle()->newGeneration();

  }
  // END: Genealogy

  // Compute maximum number of evaluations allowed for this offspring
  bool shouldMaxEvals = GAEDAConfig::handle()->getConvCrit() == GAEDAConfig::EVALS;
  int maxEvals = shouldMaxEvals ? nevals - stats.indEvals() : pop->size();

  switch (_evolutiveApproach) {
  case CentralEvolution:
    _auxPop = stepCentral(maxEvals);
    break;

  case AutonomicEvolution:
    _auxPop = stepAutonomic(maxEvals);
    break;

  default:
    return;
  }

  // Nullify the score of the genomes generated by the DE techniques
  if (_techniqueSet->hasDE()) {
    for (unsigned i=0; i<pop->size(); i++) {
      MOSGenome* mosgenome = dynamic_cast<MOSGenome*> (&pop->individual(i));
      if (mosgenome->mustBeNulled()) {
        mosgenome->setWorstScore();
        mosgenome->mustBeNulled(false);
      }
    }
    for (unsigned i=0; i<_auxPop->size(); i++) {
      MOSGenome* mosgenome = dynamic_cast<MOSGenome*> (&_auxPop->individual(i));
      if (mosgenome->mustBeNulled()) {
        mosgenome->setWorstScore();
        mosgenome->mustBeNulled(false);
      }
    }
  }

  // Evaluate fitness of new population
  _auxPop->evaluate();

  // Scale fitness of new individuals between 0 and 1
  _auxPop->scale();

  // Compute quality of offspring population
  if (_evolutiveApproach == CentralEvolution)
    evalTechQuality(_auxPop);

  // Apply elitism to old and new populations
  recombinator_->recombine (*pop, *_auxPop);
  std::swap (pop, _auxPop);

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


/**
 * Computes quality value by calling the actual quality measure used
 */
void MOSEA::evalTechQuality(GAPopulation* evalPop) {

  if      (qualityMeasure == "fAvg")      return fitnessAverageQuality(evalPop);
  else if (qualityMeasure == "fitIncAvg") return fitnessIncrementAverageQuality(evalPop);
  else if (qualityMeasure == "NSC")       return NSCQuality(evalPop);
  else if (qualityMeasure == "dAvg")      return diversityAvgQuality(evalPop);
  else if (qualityMeasure == "Compass")   return compassQuality(evalPop);

  return;

}


/**
 * Quality function: Fitness Average
 */
void MOSEA::fitnessAverageQuality(GAPopulation* evalPop) {

  std::map<techIdType,unsigned> nIndividualVector;
  std::map<techIdType,double>   acumQualityVector;

  // Initialization of both vectors
  MOSTechniqueSet::MOSTechniqueSetIterator it;
  for (it=_techniqueSet->begin(); it!=_techniqueSet->end(); it++) {
    nIndividualVector[it->first] = 0;
    acumQualityVector[it->first] = 0.0;
  }

  // Analize population to retrieve the total number of individuals generated by each technique
  // and their scores
  for (unsigned i = 0; i < evalPop->size (); i++) {
    MOSGenome& g = dynamic_cast<MOSGenome&> (evalPop->individual(i));
    techIdType techniqueId = g.getTechniqueId();
    nIndividualVector[techniqueId] = nIndividualVector[techniqueId] + 1;
    acumQualityVector[techniqueId] = acumQualityVector[techniqueId] + g.score();
  }

  // Compute the average fitness of each technique and store this value as
  // the quality value of the technique
  for (it=_techniqueSet->begin(); it!=_techniqueSet->end(); it++) {
    if (nIndividualVector[it->first] > 0)
      _techniqueSet->setTechQuality(it->first, (acumQualityVector[it->first] / nIndividualVector[it->first]));
    else
      _techniqueSet->setTechQuality(it->first, 0.0);
  }

  return;

}


/**
 * Quality function: Fitness Increment Average
 */
void MOSEA::fitnessIncrementAverageQuality(GAPopulation* evalPop) {

  std::map<techIdType,unsigned> nIndividualVector;
  std::map<techIdType,double>   acumQualityVector;

  // Initialization of both vectors
  MOSTechniqueSet::MOSTechniqueSetIterator it;
  for (it=_techniqueSet->begin(); it!=_techniqueSet->end(); it++) {
    nIndividualVector[it->first] = 0;
    acumQualityVector[it->first] = 0.0;
  }

  // Analize population to retrieve the total number of individuals generated by each technique
  // and their fitness increments
  for (unsigned i = 0; i < evalPop->size (); i++) {
    MOSGenome& g = dynamic_cast<MOSGenome&> (evalPop->individual(i));
    techIdType techniqueId = g.getTechniqueId();
    nIndividualVector[techniqueId] = nIndividualVector[techniqueId] + 1;
    acumQualityVector[techniqueId] = acumQualityVector[techniqueId] + g.getFitnessIncrement();
  }

  // Compute the average fitness increment of each technique and store this value as
  // the quality value of the technique
  for (it=_techniqueSet->begin(); it!=_techniqueSet->end(); it++) {
    if (nIndividualVector[it->first] > 0)
      _techniqueSet->setTechQuality(it->first, (acumQualityVector[it->first] / nIndividualVector[it->first]));
    else
      _techniqueSet->setTechQuality(it->first, 0.0);
  }

  return;

}


void MOSEA::diversityAvgQuality(GAPopulation* evalPop) {

   std::map<techIdType,double> diversities_avg = computeDiversitiesAvg(evalPop);
   MOSTechniqueSet::MOSTechniqueSetIterator it;

   for (it=_techniqueSet->begin(); it!=_techniqueSet->end(); it++)
      _techniqueSet->setTechQuality(it->first, diversities_avg[it->first]);

   return;

}


void MOSEA::compassQuality(GAPopulation* evalPop) {
   // The fitness average is first computed and stored in the technique set quality vector
   fitnessAverageQuality(evalPop);
   std::map<techIdType,double> diversities_avg = computeDiversitiesAvg(evalPop);

   double tmp_fitavg, tmp_divavg;
   MOSTechniqueSet::MOSTechniqueSetIterator it;

   for (it=_techniqueSet->begin(); it!=_techniqueSet->end(); it++) {
      tmp_fitavg = _techniqueSet->getTechQuality(it->first);
      tmp_divavg = diversities_avg[it->first];

      double compassvalue = (tmp_fitavg + tmp_divavg) / 2.0;
      _techniqueSet->setTechQuality(it->first, compassvalue);
   }

   return;

}


/**
 * Quality function: NSC
 */
void MOSEA::NSCQuality(GAPopulation* evalPop) {

  // Compute the NSC of each technique and store this value as
  // the quality value of the technique
  MOSTechniqueSet::MOSTechniqueSetIterator it;
  for (it=_techniqueSet->begin(); it!=_techniqueSet->end(); it++) {
    double quality = _nsc->computeNSC(evalPop, it->first);
    _techniqueSet->setTechQuality(it->first, quality);
  }

  return;

}


std::map<techIdType,double> MOSEA::computeDiversitiesAvg(GAPopulation* evalPop) {

   std::map<techIdType,double>              diversities_avg;
   std::map<techIdType, vector<GAGenome*> > subpops;

   for (unsigned i=0; i<evalPop->size(); i++) {
      MOSGenome& ind         = (MOSGenome&)evalPop->individual(i);
      techIdType techniqueId = ind.getTechniqueId();
      subpops[techniqueId].push_back(&ind);
   }

   std::map< techIdType,vector<GAGenome*> >::iterator it;

   for (it=subpops.begin(); it!=subpops.end(); it++) {
      techIdType         techid = it->first;
      vector<GAGenome*>& subpop = it->second;

      diversities_avg[techid] = computeDiversityAvg(subpop);
   }

   return diversities_avg;

}


/**
 * Evolve one generation using the central approach
 */
GAPopulation* MOSEA::stepCentral(unsigned maxEvals) {

  // Update participation ratios
  updatePartRatios();

  // Create new offspring
  _techniqueSet->offspringAll(pop, _auxPop, maxEvals);

  // BEGIN: Genealogy
  // Take current population to prepare the NSC data
  if (GAGenealogy::handle())
    if (GAGenealogy::handle()->isGenealogyMemory() && qualityMeasure == "NSC")
      _nsc->prepareNSCData(pop);
  // END: Genealogy

  return _auxPop;

}


/**
 * Evolve one generation using the autonomic approach
 */
GAPopulation* MOSEA::stepAutonomic(unsigned maxEvals) {

  MOSGenome *dad, *mom;
  GAGenome* newGenome;
  MOSProbVector probVectorChild;

  std::vector<double> partialInds (_techniqueSet->nTechniques(), 0);

  if (maxEvals < _auxPop->size())
    _auxPop->size(maxEvals);

  assert (pop->size() >= _auxPop->size());

  for (unsigned i=0; i<pop->size(); ) {

    mom = (MOSGenome*)&pop->select();
    dad = (MOSGenome*)&pop->select();

    // Mix probability vectors of both parents
    std::vector< const MOSProbVector* > probs (2, (MOSProbVector*)0);
    probs[0] = &dad->getProbVector();
    probs[1] = &mom->getProbVector();

    std::vector<double> scores (2, 0.0);
    scores[0] = dad->score();
    scores[1] = mom->score();

    mixProbVector(probs, scores, probVectorChild);

    // Select the technique to be used based on the probability vector
    MOSTechnique* tech = selectTechnique(probVectorChild);

    int encoding = tech->getEncoding();

    // Convert parents to the encoding of the selected technique (if necessary)
    if (!dad->existEncoding(encoding)) {
      newGenome = tech->getGenome();
      dad->addEncoding(encoding, newGenome);
      _conversion->convertGenome(dad->getDefaultEncoding(), encoding, dad->getDefaultGenome(), newGenome);
    }

    if (!mom->existEncoding(encoding)) {
      newGenome = tech->getGenome();
      mom->addEncoding(encoding, newGenome);
      _conversion->convertGenome(mom->getDefaultEncoding(), encoding, mom->getDefaultGenome(), newGenome);
    }

    // Create offspring
    tech->offspring (*dad, *mom, pop, _auxPop, _auxPop->size(), i);

    if (tech->classID() == GAID::TechniqueDE)
      tech->offspring (*dad, *mom, pop, _auxPop, _auxPop->size(), i+1);

    i+=2;
    partialInds[tech->getId()]+=2;

    // Clear vector for next iteration
    probVectorChild.clear();

  }

  // Update quality measure
  fitnessAverageQuality(_auxPop);

  // Update participation ratio
  MOSTechniqueSet::MOSTechniqueSetIterator it;

  for (it = _techniqueSet->begin(); it != _techniqueSet->end(); it++)
    it->second->setPartRatio((double) partialInds[it->first] / (double) pop->size());

  return _auxPop;

}


/**
 * selectTechnique
 * @param probVector Probability vector
 * Selects a technique based on the probability vector passed as parameter
 * PRE: Sum(probability of each technique) = 1
 */
MOSTechnique* MOSEA::selectTechnique(MOSProbVector& probVector) {

  double randomValue = GARandomDouble(0.0, 1.0);
  double acumProb = 0.0;
  MOSProbVector::iterator it;

  for(it = probVector.begin(); it!=probVector.end(); it++) {
    acumProb += it->second;
    if (randomValue <= acumProb)
      return _techniqueSet->getTechnique(it->first);
  }

  return _techniqueSet->getTechnique(it->first);

}


/**
 * Sets the population size
 */
unsigned MOSEA::populationSize(unsigned size) {
  GAGeneticAlgorithm::populationSize(size);
  _auxPop->size(size);
  return size;
}


/**
 * Gets the population size
 */
unsigned MOSEA::populationSize () const {

  if (_partFunction == dynQualityMaxPartPF)
    return _actualPopSize;
  else
    return pop->size();

}


/**
 * Sets the scaling scheme
 *
 * @param s New scaling scheme
 */
GAScalingScheme& MOSEA::scaling(const GAScalingScheme& s) {
	_auxPop->scaling(s);
	return GAGeneticAlgorithm::scaling(s);
}

/**
 * Sets the selection scheme
 *
 * @param s New selection scheme
 */
GASelectionScheme& MOSEA::selector (const GASelectionScheme& s) {
  _auxPop->selector(s);
  return GAGeneticAlgorithm::selector(s);
}