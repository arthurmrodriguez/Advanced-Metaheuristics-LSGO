#include "MOSTechniqueDEPopSize.h"

#include "DEElitism.h"
#include "MOSConversion.h"
#include "MOSTechniqueSet.h"
#include "GAEDAConfig.h"
#include "GAGenealogy.h"
#include "GAGenealogyMemory.h"
#include "logger/GALogger.h"
#include "genomes/MOSGenome.h"
#include "genomes/GA1DArrayGenome.h"

#include <algorithm>
#include <vector>

/**
 * MOSTechniqueDEPopSize Constructor
 */
MOSTechniqueDEPopSize::MOSTechniqueDEPopSize(techIdType id, std::string description, GAGenome::Initializer init,
                               GAGenome::Evaluator evaluator, encodingType encoding, GAGenome* genomeBase,
                               GAGenome::DECrossover crossover, double F, double CR, GASelectionScheme* selector ) {

  _id          = id;
  _description = description;

  _encoding    = encoding;
  _genomeBase  = genomeBase;

  _quality     = 0.0;
  _partRatio   = 0.0;

  _initializer = init;
  _evaluator   = evaluator;
  _selector    = selector;

  _crossover   = crossover;
  _F           = F;
  _CR          = CR;

  _recombinator = new DEElitism();

  return;
}

/**
 * Aux offspring for internal use
 */
void MOSTechniqueDEPopSize::offspring_internal(MOSGenome& genomeM, MOSGenome& newGenomeM, MOSGenome& genome1M, MOSGenome& genome2M, MOSGenome& genome3M) {

  GA1DArrayAlleleGenome<double>* genome    = dynamic_cast <GA1DArrayAlleleGenome<double>*> (   genomeM.getGenome(_encoding));
  GA1DArrayAlleleGenome<double>* newGenome = dynamic_cast <GA1DArrayAlleleGenome<double>*> (newGenomeM.getGenome(_encoding));
  GA1DArrayAlleleGenome<double>* genome1   = dynamic_cast <GA1DArrayAlleleGenome<double>*> (  genome1M.getGenome(_encoding));
  GA1DArrayAlleleGenome<double>* genome2   = dynamic_cast <GA1DArrayAlleleGenome<double>*> (  genome2M.getGenome(_encoding));
  GA1DArrayAlleleGenome<double>* genome3   = dynamic_cast <GA1DArrayAlleleGenome<double>*> (  genome3M.getGenome(_encoding));

  static std::vector<GAGenome*> parents (4, (GAGenome*)0);
  parents[0]= genomeM.getGenome(_encoding);
  parents[1]=genome1M.getGenome(_encoding);
  parents[2]=genome2M.getGenome(_encoding);
  parents[3]=genome3M.getGenome(_encoding);

  double bestParentFit = genomeM.score();
  if (GAGenome::compareScores(genome1M.score(), bestParentFit) == GAGenome::BETTER) bestParentFit = genome1M.score();
  if (GAGenome::compareScores(genome2M.score(), bestParentFit) == GAGenome::BETTER) bestParentFit = genome2M.score();
  if (GAGenome::compareScores(genome3M.score(), bestParentFit) == GAGenome::BETTER) bestParentFit = genome3M.score();

  // BEGIN: Update strategies
  MOSTechniqueSet::handle()->updateStrategies(parents, newGenomeM.getGenome(_encoding));
  // END: Update strategies

  unsigned dim = genome->length();

  // Clean old encodings (except that used by this technique)
  newGenomeM.purgeGenome(this);

  // Mutation
  for (unsigned i = 0; i < dim; i++) {
    double new_val = genome1->gene(i) + _F * (genome2->gene(i) - genome3->gene(i));

    if (new_val > newGenome->alleleset(i).upper())
      new_val = newGenome->alleleset(i).upper();
    else if (new_val < newGenome->alleleset(i).lower())
      new_val = newGenome->alleleset(i).lower();

    newGenome->gene (i, new_val);

  }

  GALogger::instance()->appendOperatorResult("Mutation result: ", *genome, *newGenome);

  // Crossover
  _stats.numcro += (*_crossover) (*genome, *newGenome, newGenome, _CR);

  GALogger::instance()->appendOperatorResult("Crossover result: ", *genome, *newGenome);

  // BEGIN: Genealogy
  // Put crossover children in the genealogy
  if(GAGenealogy::handle() != NULL)
    if (GAGenealogy::handle()->isGenealogyMemory())
      GAGenealogy::handle()->familyRelationship(genomeM, genomeM, newGenomeM, newGenomeM);
    else {

      std::vector<GAGenome*> children, parents;

      children.push_back(&newGenomeM);
      parents.push_back(&genomeM);
      parents.push_back(&genome1M);
      parents.push_back(&genome2M);
      parents.push_back(&genome3M);

      GAGenealogy::handle()->familyRelationship(parents, children);

    }
  // END: Genealogy

  _stats.numsel  += 4;
  _stats.nummut  += 1;
  _stats.numrep  += 1;
  _stats.numeval += 1;

  double newScore = newGenomeM.score();

  newGenomeM.setFitnessIncrement(newGenomeM.computeFitnessIncrement(bestParentFit));
  newGenomeM.mustComputeQuality(true);
  if (GAGenome::compareScores(newScore, bestParentFit) == GAGenome::BETTER)
    _times_improved++;

  // No need to mark children individuals to be evaluated. purgeGenome did that for us...
  return;

}

/**
 * Generate two children given two parents (used in autonomic approach)
 */
void MOSTechniqueDEPopSize::offspring(MOSGenome& dad, MOSGenome& mom, GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset) {

  MOSTechniqueSet* techSet = MOSTechniqueSet::handle();
  MOSProbVector probVectorChild;

  MOSGenome *g1, *g2, *g3;

  // Check for duplicate initial parents
  g1 = (&dad == &mom) ? (MOSGenome*)&origPop->select() : &mom;

  // Select two more individuals
  g2 = (MOSGenome*)&origPop->select();
  g3 = (MOSGenome*)&origPop->select();

  std::vector< const MOSProbVector* > probs (4, (MOSProbVector*)0);
  probs[0] = &dad.getProbVector();
  probs[1] = &g1->getProbVector();
  probs[2] = &g2->getProbVector();
  probs[3] = &g3->getProbVector();

  std::vector<double> scores (4, 0.0);
  scores[0] = dad.score();
  scores[1] = g1->score();
  scores[2] = g2->score();
  scores[3] = g3->score();

  mixProbVector(probs, scores, probVectorChild);

  techSet->bonusTechnique(probVectorChild, this->_id);

  MOSGenome& newGenome = (MOSGenome&)destPop->individual(offset);

  // Generate offspring
  offspring_internal(dad, newGenome, *g1, *g2, *g3);
  ((MOSGenome&)destPop->individual(offset)).updateProbVector(probVectorChild);

  // Avoid selection of both parent and child individuals by the global elitism.
  // This is a tradeoff solution so that the DE functioning is not heavily penalized.
  if (GAGenome::compareScores(newGenome.score(), dad.score()) < GAGenome::BETTER)
    newGenome.mustBeNulled(true);
  else
    dad.mustBeNulled(true);

  return;

}

/**
 * Creates offspring for this technique from onput population (used in central approach)
 * @param origPop Input population
 * @param destPop Output population
 * @param size Number of individuals to create
 * @param offset Offset to store individuals in ouput population
 */
void MOSTechniqueDEPopSize::offspring(GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset) {

  // Initialization of selector. This is needed as populations are continuously swapped
  // within the main algorithm and the selector object stores a reference to the population
  // it works on.
  _selector->assign(*origPop);
  _selector->update();

  MOSGenome *genome, *newGenome, *g1, *g2, *g3;
  unsigned j = 0;

  for (unsigned i = offset; i < size + offset; i++) {

    genome = dynamic_cast <MOSGenome*> (&origPop->individual(j++));
    newGenome = (MOSGenome*) (& destPop->individual(i));

    // Select three individuals different among them
    selectParents (origPop, genome, g1, g2, g3);

    // Generate offspring
    this->offspring_internal(*genome, *newGenome, *g1, *g2, *g3);

    // Avoid selection of both parent and child individuals by the global elitism.
    // This is a tradeoff solution so that the DE functioning is not heavily penalized.
    if (GAGenome::compareScores(newGenome->score(), genome->score()) < GAGenome::BETTER)
      newGenome->mustBeNulled(true);
    else
      genome->mustBeNulled(true);

  }

  return;

}

unsigned MOSTechniqueDEPopSize::offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged) {

  MOSGenome *genome, *newGenome, *g1, *g2, *g3;

  assert(pop->size() == auxPop->size());

  if (pop->size() != _rndpos.size()) {
    _rndpos.resize(pop->size());

    for (unsigned i = 0; i < _rndpos.size(); i++)
      _rndpos[i] = i;
  }

  std::random_shuffle(_rndpos.begin(), _rndpos.end());

  _times_improved = 0;

  unsigned popsize = 50;

  assert(popsize <= pop->size());

  for (unsigned i = 0; i < popsize && usedEvals < maxEvals && !converged; i++, usedEvals++) {
    genome    = dynamic_cast<MOSGenome*> (&   pop->individual(_rndpos[i]));
    newGenome = dynamic_cast<MOSGenome*> (&auxPop->individual(_rndpos[i]));

    // Select three individuals different among them
    selectParents (pop, genome, g1, g2, g3);

    // Generate offspring
    this->offspring_internal(*genome, *newGenome, *g1, *g2, *g3);

    if (newGenome->precissionReached()) converged = true;
  }

  // We return usedEvals because, in this case, the number of evaluations
  // is the same that the number of new individuals
  return usedEvals;

}

bool MOSTechniqueDEPopSize::selectParents(GAPopulation* oldPop, MOSGenome* x_i, MOSGenome*& t1, MOSGenome*& t2, MOSGenome*& t3) {

  MOSConversion* conversion = MOSConversion::handle();

    unsigned attempts=0;
    do {
      t1 = dynamic_cast <MOSGenome*> (& oldPop->select() );
      attempts += 1;
      if (attempts >= oldPop->size() ){
        attempts=0;
        break;
      }
    }while (t1 == x_i);

    do {
       t1 = dynamic_cast <MOSGenome*> (& oldPop->individual(GARandomInt(0,oldPop->size()-1)) );
    }while (t1 == x_i);

    do{
      t2 = dynamic_cast <MOSGenome*> (& oldPop->select() );
      attempts+=1;
      if (attempts >= oldPop->size() ){
           attempts=0;
           break;
      }
    }while (t2 == x_i or t2 == t1);

    do{
       t2 = dynamic_cast <MOSGenome*> (& oldPop->individual(GARandomInt(0,oldPop->size()-1)) );
    }while (t2 == x_i or t2 == t1);

    do{
      t3 = dynamic_cast <MOSGenome*> (& oldPop->select() );
      attempts+=1;
      if (attempts >= oldPop->size() ){
           attempts=0;
           break;
         }
    }while (t3 == x_i or t3 == t2 or t3 == t1);

    do{
        t3 = dynamic_cast <MOSGenome*> (& oldPop->individual(GARandomInt(0,oldPop->size()-1)) );
    }while (t3 == x_i or t3 == t2 or t3 == t1);

  GAGenome* newGenome = NULL;

  // Check for technique's encoding
  if (!t1->existEncoding(_encoding)) {
    newGenome = this->getGenome();
    t1->addEncoding(_encoding, newGenome);
    conversion->convertGenome(t1->getDefaultEncoding(), _encoding, t1->getDefaultGenome(), newGenome);
  }

  if (!t2->existEncoding(_encoding)) {
    newGenome = this->getGenome();
    t2->addEncoding(_encoding, newGenome);
    conversion->convertGenome(t2->getDefaultEncoding(), _encoding, t2->getDefaultGenome(), newGenome);
  }

  if (!t3->existEncoding(_encoding)) {
    newGenome = this->getGenome();
    t3->addEncoding(_encoding, newGenome);
    conversion->convertGenome(t3->getDefaultEncoding(), _encoding, t3->getDefaultGenome(), newGenome);
  }

  return true;

}
