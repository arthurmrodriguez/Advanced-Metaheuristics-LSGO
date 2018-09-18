/**
 * @file
 * @brief MOSTechniqueES class impl.
 *
 */

#include "MOSTechniqueES.h"

#include "garandom.h"
#include "GAGenealogy.h"
#include "GAGenealogyMemory.h"
#include "MOSConversion.h"
#include "MOSTechniqueSet.h"
#include "logger/GALogger.h"
#include "genomes/MOSGenome.h"

/**
 * MOSTechniqueES Constructor
 */

MOSTechniqueES::MOSTechniqueES(techIdType id, std::string desc, GAGenome::ESMutator mut, GAGenome::ESMutator mut_upd,
                               GAGenome::ESCrossover cx, GAGenome::ESCrossover cx_upd, GAGenome::Comparator cmp,
			       GAGenome::Initializer init, GAGenome::Evaluator eval, unsigned mu, unsigned ro,
			       unsigned lambda, encodingType encoding, GAGenome* genomeBase, GASelectionScheme* slct) {

  _id            = id;
  _description   = desc;

  _encoding      = encoding;
  _genomeBase    = genomeBase;

  _initializer   = init;
  _evaluator     = eval;
  _selector      = slct;

  _quality       = 0.0;
  _partRatio     = 0.0;

  _mutator       = mut;
  _mutator_upd   = mut_upd;
  _crossover     = cx;
  _crossover_upd = cx_upd;
  _comparator    = cmp;

  _mu            = mu;
  _ro            = ro;
  _lambda        = lambda;

  _message_printed = false;

  return;

}

/**
 * Auxiliar offspring (the one that actually performs recombination)
 */
void MOSTechniqueES::offspring_internal(const std::vector<MOSGenome*> parents, GAPopulation* destPop, unsigned offset) {

  MOSGenome& child = (MOSGenome&)destPop->individual(offset);

  std::vector<GAGenome*> pars (parents.size(), (GAGenome*)0);

  double bestParentFit = parents[0]->score();

  for (unsigned i = 0; i < parents.size(); i++) {
    pars[i] = parents[i]->getGenome(_encoding);
    if (GAGenome::compareScores(parents[i]->score(), bestParentFit) == GAGenome::BETTER)
      bestParentFit = parents[i]->score();
  }

  // Clean obsolete encoding and create current one if not present
  child.purgeGenome(this);

  // Perform crossover
  _stats.numcro += (*_crossover)(pars, child.getGenome(_encoding));

  // BEGIN: Genealogy
  // Put crossover children in the genealogy
  if(GAGenealogy::handle() != NULL) {
    if (parents.size() == 2)
      GAGenealogy::handle()->familyRelationship(*parents[0], *parents[1], child, child);
    else if (!_message_printed) {
      std::cerr << "[MOSTechniqueES] Warning: ES technique with more than two parents. Impossible to "
		<< "add a family relationship to the Genealogy. Further errors will be ignored." << std::endl;
      _message_printed = true;
    }
  }
  // END: Genealogy

  GALogger::instance()->appendOperatorResult("Crossover result: ", *parents[0], *parents[1], child, child);
  GALogger::instance()->appendInd("Individuals before mutation: ", child);

  // Mutate
  _stats.nummut += _mutator(*child.getGenome(_encoding));

  GALogger::instance()->appendInd("Individuals after mutation: ", child);

  // BEGIN: Genealogy
  // Put children in the genealogy
  if (GAGenealogy::handle() != NULL)
    if (parents.size() == 2)
      GAGenealogy::handle()->mutated(*parents[0], child, true);
  // END: Genealogy

  _stats.numeval++;
  _stats.numrep++;
  _stats.numsel += _ro;

  child.setFitnessIncrement(child.computeFitnessIncrement(bestParentFit));
  child.mustComputeQuality(true);

  // No need to mark children individuals to be evaluated. purgeGenome did that for us...
  return;
}

/**
 * Generate two children from two parents (autonomic approach)
 */
void MOSTechniqueES::offspring(MOSGenome& dad, MOSGenome& mom, GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset) {

  MOSProbVector probVectorChild;

  // Selection of parents for recombination
  static std::vector<MOSGenome*> parents (_ro, (MOSGenome*)0);

  unsigned i = 1;

  // Always add dad to the set of parents
  parents[0] = &dad;

  // If mom is different than dad, we add it as well. We also increment
  // i in one, to select one less parent
  if (&dad != &mom) {
    parents[1] = &mom;
    i++;
  }

  // Select as many extra parents as needed
  for (; i < _ro; i++)
    parents[i] = selectParent(origPop);

  // Perform offspring
  offspring_internal(parents, destPop, offset);

  std::vector< const MOSProbVector* > probs (_ro, (MOSProbVector*)0);
  for (i = 0; i < _ro; i++)
    probs[i] = &parents[i]->getProbVector();

  std::vector<double> scores (_ro, 0.0);
  for (i = 0; i < _ro; i++)
    scores[i] = parents[i]->score();

  mixProbVector(probs, scores, probVectorChild);

  MOSTechniqueSet::handle()->bonusTechnique(probVectorChild, this->_id);

  ((MOSGenome&)destPop->individual(offset)).updateProbVector(probVectorChild);

  return;

}

/**
 * Generate two children from two parents (RL approach) (Not working)
 */

void MOSTechniqueES::offspring(MOSGenome& dad, MOSGenome& mom, GAPopulation* destPop, unsigned size, unsigned offset) {

  // TODO: Arreglar esto para que sea generico

  // Selection of parents for recombination
  static std::vector<MOSGenome*> parents (2, (MOSGenome*)0);

/*
  for (unsigned i = 0; i < 2; i++)
    parents[i] = selectParent(origPop);
*/

  parents[0] = &dad;
  parents[1] = &mom;

  offspring_internal(parents, destPop, offset);

  return;

}

/**
 * Creates offspring for this technique from onput population (used in central approach)
 * @param origPop Input population
 * @param destPop Output population
 * @param size Number of individuals to create
 * @param offset Offset to store individuals in ouput population
 */
void MOSTechniqueES::offspring(GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset) {

  // Initialization of selector. This is needed as populations are continuously swapped
  // within the main algorithm and the selector object stores a reference to the population
  // it works on.
  _selector->assign(*origPop);
  _selector->update();

  // Selection of parents for recombination
  static std::vector<MOSGenome*> parents (_ro, (MOSGenome*)0);

  for (unsigned i = 0; i < size; i++) {

    for (unsigned j = 0; j < _ro; j++)
      parents[j] = selectParent(origPop);

    this->offspring_internal(parents, destPop, offset+i);

  }

  _stats.numrep += size;

  return;

}

/**
 * Method to update endogenous parameters of individuals not generated with an ES technique
 */
bool MOSTechniqueES::updateStrategies(const std::vector<GAGenome*>& parents, GAGenome* child) const {

  (*_crossover_upd)(parents, child);
  (*_mutator_upd)(*child);

  return true;

}

/**
 * Selects a parent individual from the passed population and checks if it has the appropriate encoding
 * @param pop Input population
 */
MOSGenome* MOSTechniqueES::selectParent(GAPopulation* pop) {

  MOSConversion* conversion = MOSConversion::handle();

  GAGenome* newGenome;
  MOSGenome* mosGenome;

  int indiv = GARandomInt(0, pop->size()-1);

  mosGenome = DYN_CAST(MOSGenome*, &(pop->individual(indiv)));

  // Check for technique's encoding
  if(!mosGenome->existEncoding(_encoding)) {
    newGenome = this->getGenome();
    mosGenome->addEncoding(_encoding, newGenome);
    conversion->convertGenome(mosGenome->getDefaultEncoding(), _encoding, mosGenome->getDefaultGenome(), newGenome);
  }

  return mosGenome;

}
