#include "MOSTechniqueDE.h"

#include "DEElitism.h"
#include "MOSConversion.h"
#include "MOSTechniqueSet.h"
#include "logger/GALogger.h"
#include "genomes/MOSGenome.h"
#include "genomes/GA1DArrayGenome.h"

#include <vector>

/**
 * MOSTechniqueDE Constructor
 */
MOSTechniqueDE::MOSTechniqueDE(techIdType id, std::string description, GAGenome::Initializer init,
                               GAGenome::Evaluator evaluator, encodingType encoding, GAGenome* genomeBase,
                               GAGenome::DECrossover crossover, long double F, long double CR, GASelectionScheme* selector ) {

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
void MOSTechniqueDE::offspring_internal(MOSGenome& genomeM, MOSGenome& newGenomeM, MOSGenome& genome1M, MOSGenome& genome2M, MOSGenome& genome3M) {

  GA1DArrayAlleleGenome<long double>* genome    = dynamic_cast <GA1DArrayAlleleGenome<long double>*> (   genomeM.getGenome(_encoding));
  GA1DArrayAlleleGenome<long double>* newGenome = dynamic_cast <GA1DArrayAlleleGenome<long double>*> (newGenomeM.getGenome(_encoding));
  GA1DArrayAlleleGenome<long double>* genome1   = dynamic_cast <GA1DArrayAlleleGenome<long double>*> (  genome1M.getGenome(_encoding));
  GA1DArrayAlleleGenome<long double>* genome2   = dynamic_cast <GA1DArrayAlleleGenome<long double>*> (  genome2M.getGenome(_encoding));
  GA1DArrayAlleleGenome<long double>* genome3   = dynamic_cast <GA1DArrayAlleleGenome<long double>*> (  genome3M.getGenome(_encoding));

  long double bestParentFit = genomeM.score();
  if (GAGenome::compareScores(genome1M.score(), bestParentFit) == GAGenome::BETTER) bestParentFit = genome1M.score();
  if (GAGenome::compareScores(genome2M.score(), bestParentFit) == GAGenome::BETTER) bestParentFit = genome2M.score();
  if (GAGenome::compareScores(genome3M.score(), bestParentFit) == GAGenome::BETTER) bestParentFit = genome3M.score();

  unsigned dim = genome->length();

  // Clean old encodings (except that used by this technique)
  newGenomeM.purgeGenome(this);

  // Mutation
  for (unsigned i = 0; i < dim; i++) {
    long double new_val = genome1->gene(i) + _F * (genome2->gene(i) - genome3->gene(i));

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

  _stats.numsel  += 4;
  _stats.nummut  += 1;
  _stats.numrep  += 1;
  _stats.numeval += 1;

  long double newScore = newGenomeM.score();

  newGenomeM.setFitnessIncrement(newGenomeM.computeFitnessIncrement(bestParentFit));
  newGenomeM.mustComputeQuality(true);
  if (GAGenome::compareScores(newScore, bestParentFit) == GAGenome::BETTER) {
    _times_improved++;
  }

  // No need to mark children individuals to be evaluated. purgeGenome did that for us...
  return;

}


unsigned MOSTechniqueDE::offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged) {

  MOSGenome *genome, *newGenome, *g1, *g2, *g3;

  assert(pop->size() == auxPop->size());

  _times_improved = 0;

  for (unsigned i = 0; i < pop->size() && usedEvals < maxEvals && !converged; i++, usedEvals++) {
    genome    = dynamic_cast<MOSGenome*> (&   pop->individual(i));
    newGenome = dynamic_cast<MOSGenome*> (&auxPop->individual(i));

    // Select three individuals different among them
    selectParents (pop, genome, g1, g2, g3);

    // Generate offspring
    this->offspring_internal(*genome, *newGenome, *g1, *g2, *g3);
  }

  // We return usedEvals because, in this case, the number of evaluations
  // is the same that the number of new individuals
  return usedEvals;

}

bool MOSTechniqueDE::selectParents(GAPopulation* oldPop, MOSGenome* x_i, MOSGenome*& t1, MOSGenome*& t2, MOSGenome*& t3) {

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
