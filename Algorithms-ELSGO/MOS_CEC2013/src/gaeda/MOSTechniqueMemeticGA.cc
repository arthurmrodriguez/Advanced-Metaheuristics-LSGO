/**
 * @file
 * @brief MOSTechniqueMemeticGA class impl.
 *
 */

#include "MOSTechniqueMemeticGA.h"

#include "garandom.h"
#include "GAGenealogy.h"
#include "GAGenealogyMemory.h"
#include "MOSConversion.h"
#include "MOSTechniqueSet.h"
#include "IncElitism.h"
#include "logger/GALogger.h"
#include "genomes/MOSGenome.h"

#include <stdexcept>

/**
 * MOSTechniqueMemeticGA Constructor
 */
MOSTechniqueMemeticGA::MOSTechniqueMemeticGA(techIdType id, std::string description, GAGenome::Mutator mutator, GAGenome::SexualCrossover crossover,
			       GAGenome::Comparator comparator, GAGenome::Initializer init, GAGenome::Evaluator evaluator,
			       double crossProb, double mutProb, int encoding, GAGenome* genomeBase, GASelectionScheme* selector, std::vector<LSType>& lss) {

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

  _recombinator = new IncElitism();

  _lss = lss;

  _gensCount = 0;

  return;

}


bool MOSTechniqueMemeticGA::restartRequired() {
  throw runtime_error("MOSTechniqueMemeticGA: restart method needs to be added.");
}


bool MOSTechniqueMemeticGA::restartInnerData(GAPopulation* pop) {
  return true;
}


/**
 * Aux offspring for internal use
 */
void MOSTechniqueMemeticGA::offspring_internal(MOSGenome& dad, MOSGenome& mom, GAPopulation* destPop, unsigned offset, unsigned nChilds) {

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

void MOSTechniqueMemeticGA::applyLS(MOSGenome& genome, GAPopulation* pop, unsigned k) {

  MOSGenome* gen = genome.clone();
  double   cost1 = gen->score();
  bool     mut1  = _lss[k].ls(*gen->getGenome(_encoding), _mutProb);

  if (mut1) {
    double cost2 = gen->evaluate(gaTrue);
    if (cost2 < cost1) {
      GAGenome* g = pop->replace(gen, GAPopulation::WORST);
      /* BEGIN: Genealogy */ if (GAGenealogy::handle()) GAGenealogy::handle()->deceased(*g); /* END: Genealogy */
      delete g;
    }
    else
      delete gen;
  }
  else
    delete gen;

}


unsigned MOSTechniqueMemeticGA::offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged) {

  while (usedEvals < maxEvals && !converged) {

    MOSGenome* mom = selectParent(pop);
    MOSGenome* dad = selectParent(pop);
    this->offspring_internal(*dad, *mom, auxPop, 0, 2);

    auxPop->individual(0).evaluate(gaTrue);
    auxPop->individual(1).evaluate(gaTrue);

    GAGenome* gen1 = pop->replace(auxPop->individual(0).clone(), GAPopulation::WORST);
    /* BEGIN: Genealogy */ if (GAGenealogy::handle()) GAGenealogy::handle()->deceased(*gen1); /* END: Genealogy */
    delete gen1;

    GAGenome* gen2 = pop->replace(auxPop->individual(1).clone(), GAPopulation::WORST);
    /* BEGIN: Genealogy */ if (GAGenealogy::handle()) GAGenealogy::handle()->deceased(*gen2); /* END: Genealogy */
    delete gen2;

    if (auxPop->individual(0).precissionReached() || auxPop->individual(1).precissionReached()) converged = true;
    usedEvals+=2;

    /******************* LOCAL SEARCHES *********************/

    for (unsigned k = 0; k < _lss.size(); k++) {
      // Apply the LS only each 'freq' iterations of the IncGA
      if (_gensCount % _lss[k].freq == 0) {
        // Mutate best individual in the population
        if (_lss[k].strategy == LSType::BEST || _lss[k].strategy == LSType::ALL) {
          MOSGenome& best = dynamic_cast<MOSGenome&> (pop->best());
          applyLS(best, pop, k);
        }

        // Mutate a random individual from the population
        if (_lss[k].strategy == LSType::RANDOM || _lss[k].strategy == LSType::ALL) {
          int r = GARandomInt(0, pop->size() - 1);
          MOSGenome& rnd = dynamic_cast<MOSGenome&> (pop->individual(r));
          applyLS(rnd, pop, k);
        }
      }
    }

    _gensCount++;

  }

  // We return usedEvals because, in this case, the number of evaluations
  // is the same that the number of new individuals
  return usedEvals;

}

/**
 * Selects a parent individual from the passed population and checks if it has the appropriate encoding
 * @param pop Input population
 */
MOSGenome* MOSTechniqueMemeticGA::selectParent(GAPopulation* pop) {

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
