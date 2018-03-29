#include "MOSTechniqueSTSDE.h"

#include "MOSConversion.h"
#include "MOSTechniqueSet.h"
#include "GAEDAConfig.h"
#include "GAGenealogy.h"
#include "GAGenealogyMemory.h"
#include "DEElitism.h"
#include "PopElitism.h"
#include "logger/GALogger.h"
#include "genomes/MOSGenome.h"
#include "genomes/GA1DArrayGenome.h"

#include <vector>

/**
 * MOSTechniqueSTSDE Constructor
 */
MOSTechniqueSTSDE::MOSTechniqueSTSDE(techIdType id, std::string description, GAGenome::Initializer init,
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

  const GA1DArrayAlleleGenome<double>* ind = dynamic_cast <const GA1DArrayAlleleGenome<double>* > (genomeBase); assert(ind != NULL);

  min_range_ = new double[ind->length()];
  max_range_ = new double[ind->length()];

  for (unsigned i=0; i<ind->length(); i++) {
    min_range_[i] = ind->alleleset(0).lower();
    max_range_[i] = ind->alleleset(0).upper();
  }

  return;
}


/**
 * Aux offspring for internal use
 */
void MOSTechniqueSTSDE::offspring_internal(MOSGenome& genomeM, MOSGenome& newGenomeM, MOSGenome& genome1M, MOSGenome& genome2M, MOSGenome& genome3M) {

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
 * Creates offspring for this technique from onput population (used in central approach)
 * @param origPop Input population
 * @param destPop Output population
 * @param size Number of individuals to create
 * @param offset Offset to store individuals in ouput population
 */
void MOSTechniqueSTSDE::offspring(GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset) {

  GAPopulation stsPop(*origPop);

  // Initialization of selector. This is needed as populations are continuously swapped
  // within the main algorithm and the selector object stores a reference to the population
  // it works on.
  _selector->assign(*origPop);
  _selector->update();

  if (GARandomDouble(0.0,1.0) <= GAEDAConfig::handle()->getSTSDEProb()) {
    updateRanges(origPop);
    genSTSPop   (origPop, &stsPop);

    PopElitism elitism;

    elitism.recombine(*origPop, stsPop);
  }
  else
    stsPop.copy(*origPop);

  MOSGenome *genome, *newGenome, *g1, *g2, *g3;
  unsigned j = 0;

  for (unsigned i = offset; i < size + offset; i++) {

    genome = dynamic_cast <MOSGenome*> (&stsPop.individual(j++));
    newGenome = (MOSGenome*) (& destPop->individual(i));

    // Select three individuals different among them
    selectParents (&stsPop, genome, g1, g2, g3);

    // Store current target vector to be used in the self-adaptive derived class
    _curTarget = i;

    // Generate offspring
    this->offspring_internal(*genome, *newGenome, *g1, *g2, *g3);

    // Avoid selection of both parent and child individuals by the global elitism.
    // This is a tradeoff solution so that the DE functioning is not heavily penalized.
    if (GAGenome::compareScores(newGenome->score(), genome->score()) < GAGenome::BETTER) {
      GAGenome* tmp;
      tmp = destPop->replace(genome, newGenome);
      assert (tmp == newGenome);
      tmp = stsPop.replace(newGenome, genome);
      assert (tmp == genome);
    }

  }

  return;

}


unsigned MOSTechniqueSTSDE::offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged) {

  GAPopulation stsPop(*pop);

  if (GARandomDouble(0.0, 1.0) <= GAEDAConfig::handle()->getSTSDEProb() && maxEvals > pop->size()) {
    updateRanges(pop);
    genSTSPop   (pop, &stsPop);
    usedEvals += pop->size();

    PopElitism elitism;
    elitism.recombine(*pop, stsPop);
  }
  else
    stsPop.copy(*pop);

  MOSGenome *genome, *newGenome, *g1, *g2, *g3;

  assert(pop->size() == auxPop->size());

  _times_improved = 0;

  for (unsigned i = 0; i < pop->size() && usedEvals < maxEvals && !converged; i++, usedEvals++) {

    genome = dynamic_cast <MOSGenome*> (&stsPop.individual(i));
    newGenome = (MOSGenome*) (& auxPop->individual(i));

    // Select three individuals different among them
    selectParents (&stsPop, genome, g1, g2, g3);

    // Store current target vector to be used in the self-adaptive derived class
    _curTarget = i;

    // Generate offspring
    this->offspring_internal(*genome, *newGenome, *g1, *g2, *g3);

    // Avoid selection of both parent and child individuals by the global elitism.
    // This is a tradeoff solution so that the DE functioning is not heavily penalized.
    if (GAGenome::compareScores(newGenome->score(), genome->score()) < GAGenome::BETTER) {
      GAGenome* tmp;
      tmp = auxPop->replace(genome, newGenome);
      assert (tmp == newGenome);
      tmp = stsPop.replace(newGenome, genome);
      assert (tmp == genome);
    }

    if (newGenome->precissionReached()) converged = true;

  }

  return usedEvals;

}

void MOSTechniqueSTSDE::updateRanges(GAPopulation* pop) {

  MOSGenome& tmp_indM = dynamic_cast<MOSGenome&>(pop->individual(0));
  GA1DArrayAlleleGenome<double> *tmp_ind = dynamic_cast<GA1DArrayAlleleGenome<double>*>(tmp_indM.getGenome(_encoding));
  assert(tmp_ind != NULL);
  int num_dims = tmp_ind->length();

  for (int dim_pos=0; dim_pos<num_dims; dim_pos++) {

    GA1DArrayAlleleGenome<double>& first_ind = dynamic_cast<GA1DArrayAlleleGenome<double>&>(*(dynamic_cast<MOSGenome&>(pop->individual(0)).getGenome(_encoding)));
    min_range_[dim_pos] = max_range_[dim_pos] = first_ind.gene(dim_pos);

    for (int ind_pos=1; ind_pos<pop->size(); ind_pos++) {
      GA1DArrayAlleleGenome<double>& ind_i = dynamic_cast<GA1DArrayAlleleGenome<double>&>(*(dynamic_cast<MOSGenome&>(pop->individual(ind_pos)).getGenome(_encoding)));

      if      (ind_i.gene(dim_pos) < min_range_[dim_pos]) min_range_[dim_pos] = ind_i.gene(dim_pos);
      else if (ind_i.gene(dim_pos) > max_range_[dim_pos]) max_range_[dim_pos] = ind_i.gene(dim_pos);
    }

  }

}


void MOSTechniqueSTSDE::genSTSPop(GAPopulation* pop, GAPopulation* stsPop) {

  double k = GARandomDouble(0.0,1.0);
  double new_value;

  for (int ind_pos=0; ind_pos<pop->size(); ind_pos++) {

    MOSGenome& indM     = dynamic_cast<MOSGenome&>(     pop->individual(ind_pos));
    MOSGenome& sts_indM = dynamic_cast<MOSGenome&>(stsPop->individual(ind_pos));

    // Clean old encodings (except that used by this technique)
    sts_indM.purgeGenome(this);

    GA1DArrayAlleleGenome<double>& ind     = dynamic_cast<GA1DArrayAlleleGenome<double>&>(    *indM.getGenome(_encoding));
    GA1DArrayAlleleGenome<double>& sts_ind = dynamic_cast<GA1DArrayAlleleGenome<double>&>(*sts_indM.getGenome(_encoding));

    for (int dim_pos=0; dim_pos<ind.length(); dim_pos++) {
      new_value = k * (min_range_[dim_pos] + max_range_[dim_pos]) - ind.gene(dim_pos);
      if (new_value < min_range_[dim_pos] ||
          new_value > max_range_[dim_pos]) {
        new_value = GARandomDouble(min_range_[dim_pos],max_range_[dim_pos]);
      }

      sts_ind.gene(dim_pos,new_value);
    }

  }

  stsPop->evaluate(gaTrue);
  _stats.numeval += stsPop->size();

}


bool MOSTechniqueSTSDE::selectParents(GAPopulation* oldPop, MOSGenome* x_i, MOSGenome*& t1, MOSGenome*& t2, MOSGenome*& t3) {

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
