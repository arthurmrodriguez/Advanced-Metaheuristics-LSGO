/**
 * @file
 * @brief MOSTechniqueGA class impl.
 *
 */

#include "MOSTechniqueGA.h"

#include "garandom.h"
#include "GAGenealogy.h"
#include "GAGenealogyMemory.h"
#include "MOSConversion.h"
#include "MOSTechniqueSet.h"
#include "PopElitism.h"
#include "logger/GALogger.h"
#include "genomes/MOSGenome.h"

/**
 * MOSTechniqueGA Constructor
 */
MOSTechniqueGA::MOSTechniqueGA(techIdType id, std::string description, GAGenome::Mutator mutator, GAGenome::SexualCrossover crossover,
			       GAGenome::Comparator comparator, GAGenome::Initializer init, GAGenome::Evaluator evaluator,
			       double crossProb, double mutProb, int encoding, GAGenome* genomeBase, GASelectionScheme* selector) {

  _id          = id;
  _description = description;

  _encoding    = encoding;
  _genomeBase  = genomeBase;

  _initializer = init;
  _evaluator   = evaluator;
  _selector    = selector;

  _quality     = 0.0;
  _partRatio   = 0.0;

  _mutator     = mutator;
  _crossover   = crossover;
  _comparator  = comparator;
  _crossProb   = crossProb;
  _mutProb     = mutProb;

  _recombinator = new PopElitism();

  return;

}

/**
 * Aux offspring for internal use
 */
void MOSTechniqueGA::offspring_internal(MOSGenome& dad, MOSGenome& mom, GAPopulation* destPop, unsigned offset, unsigned nChilds) {

  int c1 = 0, c2 = 0, mut = 0;

  MOSGenome* child1 = (MOSGenome*)&destPop->individual(offset);
  MOSGenome* child2 = NULL;

  if (nChilds > 1)
    child2 = (MOSGenome*)&destPop->individual(offset+1);

  static std::vector<GAGenome*> parents (2, (GAGenome*)0);
  parents[0]=mom.getGenome(_encoding);
  parents[1]=dad.getGenome(_encoding);

  double bestParentFit = GAGenome::compareScores(dad.score(), mom.score()) == GAGenome::BETTER ? dad.score() : mom.score();

  // BEGIN: Update strategies
  MOSTechniqueSet::handle()->updateStrategies(parents, child1->getGenome(_encoding));
  if (nChilds>1)
    MOSTechniqueSet::handle()->updateStrategies(parents, child2->getGenome(_encoding));
  // END: Update strategies

  // BEGIN: Genealogy
  unsigned children = 0, mut1 = 0, mut2 = 0;
  // END: Genealogy

  _stats.numsel += 2;

  switch (nChilds) {
  case 2:

    // Clean old encodings (except that used by this technique)
    child1->purgeGenome(this);
    child2->purgeGenome(this);

    // Crossover
    if(GAFlipCoin(_crossProb)) {
      _stats.numcro += (*_crossover) (*mom.getGenome(_encoding), *dad.getGenome(_encoding), child1->getGenome(_encoding), child2->getGenome(_encoding));
      children = c1 = c2 = 1;
    }
    else {
      child1->inherit(mom);
      child2->inherit(dad);
    }

    // BEGIN: Genealogy
    // Put crossover children in the genealogy
    if(children)
      if(GAGenealogy::handle() != NULL)
        GAGenealogy::handle()->familyRelationship(dad, mom, *child1, *child2);
    // END: Genealogy

    GALogger::instance()->appendOperatorResult("Crossover result: ", dad, mom,  *child1, *child2);
    GALogger::instance()->appendOperatorResult("Individuals before mutation: ", *child1, *child2);

    // Mutations
    _stats.nummut += (mut = _mutator(*(child1->getGenome(_encoding)), _mutProb));
    if(mut > 0)
      mut1 = c1 = 1;

    _stats.nummut += (mut = _mutator(*(child2->getGenome(_encoding)), _mutProb));
    if(mut > 0)
      mut2 = c2 = 1;

    GALogger::instance()->appendOperatorResult("Individuals after mutation: ", *child1, *child2);

    // BEGIN: Genealogy
    // Put children in the genealogy
    if (GAGenealogy::handle() != NULL) {
      if (mut1) {
        if(children) // Individual was crossed and mutated
          GAGenealogy::handle()->mutated(mom, *child1, true);
        else // Individual was only mutated
          GAGenealogy::handle()->mutated(mom, *child1, false);
      }
      if (mut2) {
        if(children) // Individual was crossed and mutated
          GAGenealogy::handle()->mutated(dad, *child2, true);
        else // Individual was only mutated
          GAGenealogy::handle()->mutated(dad, *child2, false);
      }
    }
    // END: Genealogy

    _stats.numeval += c1 + c2;
    _stats.numrep += 2;
    break;

  case 1:
    // Clean old encodings (except that used by this technique)
    child1->purgeGenome(this);

    // Crossover
    if (GAFlipCoin(_crossProb)) {
      _stats.numcro += (*_crossover) (*mom.getGenome(_encoding), *dad.getGenome(_encoding), child1->getGenome(_encoding), (GAGenome*)0);
      children = c1 = 1;
    }
    else {
      if (GARandomBit()) {
        child1->inherit(mom);
        mut2 = 1;
      }
      else
        child1->inherit(dad);
    }

    // BEGIN: Genealogy
    // Put crossover children in the genealogy
    if(children)
      if(GAGenealogy::handle() != NULL)
        GAGenealogy::handle()->familyRelationship(dad, mom, *child1, *child1);
    // END: Genealogy

    GALogger::instance()->appendOperatorResult("Crossover result: ", dad, mom, *child1);
    GALogger::instance()->appendInd("Individuals before mutation: ", *child1);

    // Mutation
    _stats.nummut += (mut = _mutator(*(child1->getGenome(_encoding)), _mutProb));
    if(mut > 0)
      mut1 = c1 = 1;

    GALogger::instance()->appendInd("Individuals after mutation: ",  *child1);

    // BEGIN: Genealogy
    //Put mutated childrens in the genealogy
    if(GAGenealogy::handle() != NULL) {
      if(children && mut1 && mut2) // Individual was crossed and mutated
        GAGenealogy::handle()->mutated(mom, *child1, true);
      if(!children && mut1 && mut2) // Individual was only mutated
        GAGenealogy::handle()->mutated(mom, *child1, false);
      if(children && mut1 && !mut2) // Individual was crossed and mutated
        GAGenealogy::handle()->mutated(dad, *child1, true);
      if(!children && mut1 && !mut2) // Individual was only mutated
        GAGenealogy::handle()->mutated(dad, *child1, false);
    }
    // END: Genealogy

    _stats.numeval += c1;
    _stats.numrep += 1;
    break;

  default:
    break;

  }

  // Mark children individuals to be evaluated. This is necessary because, although purgeGenome
  // initially marked individuals for re-evaluation, a call to inherit when no crossover is performed
  // restores the _evaluated flag to true (as it should be). However, a mutation afterwards only
  // modifies the _evaluated flag of the internal genome, and not that of MOSGenome, the one that
  // will be checked when evaluating populations.
  if (c1)
      child1->evaluated(gaFalse);

  if (c2)
    child2->evaluated(gaFalse);

  if (c1) {
    child1->setFitnessIncrement(child1->computeFitnessIncrement(bestParentFit));
    child1->mustComputeQuality(true);
  }

  if (c2) {
    child2->setFitnessIncrement(child2->computeFitnessIncrement(bestParentFit));
    child2->mustComputeQuality(true);
  }

  return;

}

/**
 * Generate two children given two parents (used in autonomic approach)
 */
void MOSTechniqueGA::offspring(MOSGenome& dad, MOSGenome& mom, GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset) {

  MOSTechniqueSet* techSet = MOSTechniqueSet::handle();
  MOSProbVector probVectorChild;

  std::vector< const MOSProbVector* > probs (2, (MOSProbVector*)0);
  probs[0] = &dad.getProbVector();
  probs[1] = &mom.getProbVector();

  std::vector<double> scores (2, 0.0);
  scores[0] = dad.score();
  scores[1] = mom.score();

  mixProbVector(probs, scores, probVectorChild);

  techSet->bonusTechnique(probVectorChild, this->_id);

  if (size == offset + 1) {
    // If size is odd and we are generating the last individual
    offspring_internal(dad, mom, destPop, offset, 1);
    ((MOSGenome&)destPop->individual(offset)).updateProbVector(probVectorChild);
  } else {
    // Normal behavior: two children are crossed, mutated and updated
    offspring_internal(dad, mom, destPop, offset, 2);
    ((MOSGenome&)destPop->individual(offset)).updateProbVector(probVectorChild);
    ((MOSGenome&)destPop->individual(offset+1)).updateProbVector(probVectorChild);
  }

  return;

}


/**
 * Generate two children given two parents (used in RL approach)
 */
void MOSTechniqueGA::offspring(MOSGenome& dad, MOSGenome& mom, GAPopulation* destPop, unsigned size, unsigned offset) {

  if ((size % 2 != 0) && (size == offset+1))
    // If size is odd and we are generating the last individual
    offspring_internal(dad, mom, destPop, offset, 1);
  else
    // Normal behavior: two children are crossed, mutated and updated
    offspring_internal(dad, mom, destPop, offset, 2);

  return;

}

/**
 * Creates offspring for this technique from onput population (used in central approach)
 * @param origPop Input population
 * @param destPop Output population
 * @param size Number of individuals to create
 * @param offset Offset to store individuals in ouput population
 */
void MOSTechniqueGA::offspring(GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset) {

  // Initialization of selector. This is needed as populations are continuously swapped
  // within the main algorithm and the selector object stores a reference to the population
  // it works on.
  _selector->assign(*origPop);
  _selector->update();

  MOSGenome *mom, *dad;

  unsigned i;

  for (i = 0; i < size - 1; i += 2) {
    mom = selectParent(origPop);
    dad = selectParent(origPop);
    this->offspring_internal(*dad, *mom, destPop, offset + i, 2);
  }

  // If size is odd
  if (size % 2 != 0){
    dad = selectParent(origPop);
    mom = selectParent(origPop);
    this->offspring_internal(*dad, *mom, destPop, offset + i, 1);
  }

  return;

}

unsigned MOSTechniqueGA::offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged) {

  unsigned i;
  MOSGenome *mom, *dad;

  unsigned size = 0;

  if (maxEvals >= pop->size())
    if (pop->size() % 2 == 0)
      size = pop->size();
    else
      size = pop->size() - 1;
  else
    if (maxEvals % 2 == 0)
      size = maxEvals;
    else
      size = maxEvals - 1;

  for (i = 0; i < size && usedEvals < maxEvals && !converged; i+=2, usedEvals += 2) {
    mom = selectParent(pop);
    dad = selectParent(pop);
    this->offspring_internal(*dad, *mom, auxPop, i, 2);
    if (auxPop->individual(i).precissionReached() || auxPop->individual(i+1).precissionReached()) converged = true;
  }

  // If pop->size() is odd
  if (usedEvals < maxEvals && usedEvals < pop->size() && !converged) {
    dad = selectParent(pop);
    mom = selectParent(pop);
    this->offspring_internal(*dad, *mom, auxPop, i, 1);
    usedEvals++;
    if (auxPop->individual(i).precissionReached()) converged = true;
  }

  // We return usedEvals because, in this case, the number of evaluations
  // is the same that the number of new individuals
  return usedEvals;

}

/**
 * Selects a parent individual from the passed population and checks if it has the appropriate encoding
 * @param pop Input population
 */
MOSGenome* MOSTechniqueGA::selectParent(GAPopulation* pop) {

  MOSConversion* conversion = MOSConversion::handle();

  GAGenome*  newGenome;
  MOSGenome* mosGenome = (MOSGenome*)&(_selector->select());

  // Check for technique's encoding
  if (!mosGenome->existEncoding(_encoding)) {
    newGenome = this->getGenome();
    mosGenome->addEncoding(_encoding, newGenome);
    conversion->convertGenome(mosGenome->getDefaultEncoding(), _encoding, mosGenome->getDefaultGenome(), newGenome);
  }

  return mosGenome;

}
