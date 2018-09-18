/**
 * @file
 * @brief MOSTechniqueEDA class impl.
 *
 */

#include "MOSTechniqueEDA.h"

#include "MOSConversion.h"
#include "genomes/MOSGenome.h"

MOSTechniqueEDA::MOSTechniqueEDA(techIdType id, std::string description, double selPerc, GAGraphModel::ModelType networkType,
				     GABayesianNetwork::LearningMethod bayesLearnMethod, GABayesianNetwork::EBNALocalScoring localScoring,
				     GABayesianNetwork::SimulationMethod simMethod, encodingType encoding, GAGenome* genomeBase,
				     GAGenome::Initializer init, GAGenome::Evaluator evaluator) {

  _id          = id;
  _description = description;

  _encoding    = encoding;
  _genomeBase  = genomeBase;

  _initializer = init;
  _evaluator   = evaluator;
  _selector    = NULL;

  _selPerc     = selPerc;

  if (networkType == GAGraphModel::BAYESIAN_NETWORK)
    _network = (GAGraphModel*) new GABayesianNetwork(*genomeBase, bayesLearnMethod, localScoring, simMethod);
  else
    _network = NULL;

  return;

}

MOSTechniqueEDA::MOSTechniqueEDA(techIdType id, std::string description, double selPerc, GAGraphModel::ModelType networkType,
				     GAGaussianNetwork::LearningMethod gaussLearnMethod, GAGaussianNetwork::ScoreMethod scoreMethod, 
				     encodingType encoding, GAGenome* genomeBase, GAGenome::Initializer init, GAGenome::Evaluator evaluator) {

  _id          = id;
  _description = description;

  _encoding    = encoding;
  _genomeBase  = genomeBase;

  _initializer = init;
  _evaluator   = evaluator;
  _selector    = NULL;

  _selPerc     = selPerc;

  if (networkType == GAGraphModel::GAUSSIAN_NETWORK)
    _network = (GAGraphModel*) new GAGaussianNetwork(*genomeBase, gaussLearnMethod, scoreMethod);
  else
    _network = NULL;

  return;

}

MOSTechniqueEDA::~MOSTechniqueEDA() {
  delete _network;
}

/**
 * Creates offspring for this technique from onput population (used in central approach)
 * @param origPop Input population
 * @param destPop Output population
 * @param size Number of individuals to create
 * @param offset Offset to store individuals in ouput population
 */
void MOSTechniqueEDA::offspring (GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset) {

  MOSConversion* conversion = MOSConversion::handle();

  MOSGenome* mosGenome;
  GAGenome* newGenome;

  for (unsigned i = 0; i < origPop->size(); i++) {

    mosGenome = (MOSGenome*) &(origPop->individual(i));

    // Tell MOSGenome which is the next EDA technique to make wrappers for nominal and continuous genes work properly
    mosGenome->setNextEdaTech(this);

    // Ensure that every individual in orig population has the encoding used by this EDA technique
    if (!mosGenome->existEncoding(_encoding)) {
      newGenome = this->getGenome();
      mosGenome->addEncoding(_encoding, newGenome);
      conversion->convertGenome(mosGenome->getDefaultEncoding(), _encoding, mosGenome->getDefaultGenome(), newGenome);
    }

  }

  // Learn probabilistic model
  _network->learn(*origPop, _selPerc);

  // Sample as many individuals as requested
  for (unsigned i = 0; i < size; i++) {
      mosGenome = (MOSGenome*) &(destPop->individual(i+offset));

      // Convert individuals in dest population to EDA's encoding and clean any other encoding
      mosGenome->purgeGenome(this);

      // Create a new individual
      _network->simulate(*mosGenome->getGenome(_encoding));
  }

  mosGenome->setFitnessIncrement(mosGenome->computeFitnessIncrement(origPop->ave()));
  mosGenome->mustComputeQuality(true);

  _stats.numeval += size;
  _stats.numrep  += size;
  _stats.numcro  += size;
  _stats.numsel  += size;

}

