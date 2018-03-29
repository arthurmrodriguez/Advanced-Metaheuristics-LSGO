/**
 * @file
 * @brief MOSTechniqueNMA class impl.
 *
 */

#include "MOSTechniqueNMA.h"

#include "garandom.h"
#include "GenerationalElitism.h"
#include "MOSConversion.h"
#include "genomes/MOSGenome.h"

#include <sstream>

/**
 * MOSTechniqueNMA Constructor
 */
MOSTechniqueNMA::MOSTechniqueNMA(techIdType id, std::string description, GAGenome::Comparator comparator,
                                 GAGenome::Initializer init, GAGenome::Evaluator evaluator, int encoding, GAGenome* genomeBase,
                                 GASelectionScheme* selector) {

  _id          = id;
  _description = description;

  _encoding    = encoding;
  _genomeBase  = genomeBase;

  _initializer = init;
  _evaluator   = evaluator;
  _selector    = selector;

  _quality     = 0.0;
  _partRatio   = 0.0;

  _comparator  = comparator;

  _recombinator = new GenerationalElitism();

  return;

}

/**
 * Creates offspring for this technique from input population (used in central approach)
 * @param origPop Input population
 * @param destPop Output population
 * @param size Number of individuals to create
 * @param offset Offset to store individuals in ouput population
 */
void MOSTechniqueNMA::offspring(GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset) {

  // Initialization of selector. This is needed as populations are continuously swapped
  // within the main algorithm and the selector object stores a reference to the population
  // it works on.
  _selector->assign(*origPop);
  _selector->update();

  GAPopulation* tmpPop = origPop->clone();
  unsigned totalEvals = 0;
  unsigned branch;

  for (unsigned i = 0; i < tmpPop->size(); i++)
    (dynamic_cast<MOSGenome&>(tmpPop->individual(i))).purgeGenome(this);

  do {
    NMA(*tmpPop, totalEvals, size, branch);
  } while (totalEvals < size);

  for (unsigned i = 0; i < size; i++)
    destPop->individual(offset+i).copy(tmpPop->individual(i));

  delete tmpPop;

  // Remind to update stats
  _stats.numsel+=1;
  _stats.numrep+=1;
  _stats.nummut+=0;
  _stats.numcro+=0;
  _stats.numeval+=totalEvals;

  return;

}


unsigned MOSTechniqueNMA::offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged) {

  unsigned totalEvals = 0;
  unsigned branch;

  for (unsigned i = 0; i < auxPop->size(); i++)
    (dynamic_cast<MOSGenome&>(auxPop->individual(i))).purgeGenome(this);

//  do {
    NMA(*auxPop, usedEvals, maxEvals, branch);
//  } while (totalEvals < maxEvals);

//  std::cout << "usedEvals: " << usedEvals << std::endl;

  // Remind to update stats
  _stats.numsel+=1;
  _stats.numrep+=1;
  _stats.nummut+=0;
  _stats.numcro+=0;
  _stats.numeval+=usedEvals;

  switch (branch) {
    case Reflection:
//      std::cout << "Reflection" << std::endl;
      return usedEvals;
      break;
    case Expansion:
//      std::cout << "Expansion" << std::endl;
      return usedEvals - 1;
      break;
    case ContractionInOK:
//      std::cout << "Contraction IN" << std::endl;
      return usedEvals - 1;
      break;
    case ContractionOutOK:
//      std::cout << "Contraction OUT" << std::endl;
      return usedEvals - 1;
      break;
    case Shrink:
      return usedEvals - 2;
      break;
    default:
//      std::cout << "Default" << std::endl;
      return usedEvals;
  }

}


// NMA helper functions

MOSGenome* MOSTechniqueNMA::centroidOf(const GAPopulation& pop) const {

  MOSGenome*                     centroid     = dynamic_cast<MOSGenome*>(pop.best().clone());
  GA1DArrayAlleleGenome<long double>& centroidReal = dynamic_cast<GA1DArrayAlleleGenome<long double>&>(*centroid->getGenome(_encoding));

  for (unsigned dim = 0; dim < centroidReal.length(); dim++) {

    long double new_val = 0.0;

    for (unsigned ind_pos = 0; ind_pos < pop.size(); ind_pos++) {
      MOSGenome& tmp_ind = dynamic_cast<MOSGenome&>(pop.individual(ind_pos));
      new_val           += dynamic_cast<GA1DArrayAlleleGenome<long double>&>(*tmp_ind.getGenome(_encoding)).gene(dim);
    }

    centroidReal.gene(dim, new_val / pop.size());

  }

  // Do not forget to evaluate the genome (as we work directly on the inner genome)
  centroid->evaluated(gaFalse);

  return centroid;

}


MOSGenome* MOSTechniqueNMA::reflection(const MOSGenome& centroid, const MOSGenome& worst_ind, long double rho) const {
  // Reflection: the NMA generates a trial solution xr by reflecting the worst individual through the opposite face of the polyhedron determined by all the points
  MOSGenome*                     reflx_ind      = dynamic_cast<MOSGenome*>(centroid.clone());
  GA1DArrayAlleleGenome<long double>& reflx_ind_real = dynamic_cast<GA1DArrayAlleleGenome<long double>&>(*reflx_ind->getGenome(_encoding));
  GA1DArrayAlleleGenome<long double>& centroid_real  = dynamic_cast<GA1DArrayAlleleGenome<long double>&>(  *centroid.getGenome(_encoding));
  GA1DArrayAlleleGenome<long double>& worst_ind_real = dynamic_cast<GA1DArrayAlleleGenome<long double>&>( *worst_ind.getGenome(_encoding));

  for (unsigned dim = 0; dim < reflx_ind_real.length(); dim++)
    reflx_ind_real.gene(dim, centroid_real.gene(dim) + rho * (centroid_real.gene(dim) - worst_ind_real.gene(dim)));

  // Do not forget to evaluate the genome (as we work directly on the inner genome)
  reflx_ind->evaluate(gaTrue);

  long double bestParentFit = worst_ind.score();
  reflx_ind->setFitnessIncrement(reflx_ind->computeFitnessIncrement(bestParentFit));
  reflx_ind->mustComputeQuality(true);

  return reflx_ind;
}


void MOSTechniqueNMA::expansion(const MOSGenome& centroid, const MOSGenome& reflx_ind, MOSGenome& worst_ind, long double chi, stringstream& msg) const {
  // Expansion: NMA further exploits a promising search direction by applying the expansion operation
  GA1DArrayAlleleGenome<long double>& reflx_ind_real = dynamic_cast<GA1DArrayAlleleGenome<long double>&>(*reflx_ind.getGenome(_encoding));
  GA1DArrayAlleleGenome<long double>& centroid_real  = dynamic_cast<GA1DArrayAlleleGenome<long double>&>( *centroid.getGenome(_encoding));

  MOSGenome*                     expan_ind      = dynamic_cast<MOSGenome*>(centroid.clone());
  GA1DArrayAlleleGenome<long double>& expan_ind_real = dynamic_cast<GA1DArrayAlleleGenome<long double>&>(*expan_ind->getGenome(_encoding));

  long double max_value = centroid_real.alleleset(0).upper();
  long double min_value = centroid_real.alleleset(0).lower();

  for (unsigned dim = 0; dim < expan_ind_real.length(); dim++) {
    long double new_gene = centroid_real.gene(dim) + chi * (reflx_ind_real.gene(dim) - centroid_real.gene(dim));
    if      (new_gene > max_value) new_gene = max_value;
    else if (new_gene < min_value) new_gene = min_value;
    expan_ind_real.gene(dim, new_gene);
  }

  // Do not forget to evaluate the genome (as we work directly on the inner genome)
  expan_ind->evaluate(gaTrue);

  long double bestParentFit = reflx_ind.score();
  expan_ind->setFitnessIncrement(expan_ind->computeFitnessIncrement(bestParentFit));
  expan_ind->mustComputeQuality(true);

  if (GAGenome::compareScores(expan_ind->score(), reflx_ind.score()) == GAGenome::BETTER) worst_ind.copy(*expan_ind);
  else                                                                                    worst_ind.copy( reflx_ind);

  delete expan_ind;
}


MOSGenome* MOSTechniqueNMA::outsideContraction(const MOSGenome& centroid, const MOSGenome& reflx_ind, long double gamma) const {
  MOSGenome*                     contr_ind      = dynamic_cast<MOSGenome*>(centroid.clone());
  GA1DArrayAlleleGenome<long double>& contr_ind_real = dynamic_cast<GA1DArrayAlleleGenome<long double>&>(*contr_ind->getGenome(_encoding));
  GA1DArrayAlleleGenome<long double>& centroid_real  = dynamic_cast<GA1DArrayAlleleGenome<long double>&>(  *centroid.getGenome(_encoding));
  GA1DArrayAlleleGenome<long double>& reflx_ind_real = dynamic_cast<GA1DArrayAlleleGenome<long double>&>( *reflx_ind.getGenome(_encoding));

  for (unsigned dim = 0; dim < contr_ind_real.length(); dim++)
    contr_ind_real.gene(dim, centroid_real.gene(dim) + gamma * (reflx_ind_real.gene(dim) - centroid_real.gene(dim)));

  // Do not forget to evaluate the genome (as we work directly on the inner genome)
  contr_ind->evaluate(gaTrue);

  long double bestParentFit = reflx_ind.score();
  contr_ind->setFitnessIncrement(contr_ind->computeFitnessIncrement(bestParentFit));
  contr_ind->mustComputeQuality(true);

  return contr_ind;
}


MOSGenome* MOSTechniqueNMA::insideContraction(const MOSGenome& centroid, const MOSGenome& worst_ind, long double gamma) const {
  MOSGenome*                     contr_ind      = dynamic_cast<MOSGenome*>(centroid.clone());
  GA1DArrayAlleleGenome<long double>& contr_ind_real = dynamic_cast<GA1DArrayAlleleGenome<long double>&>(*contr_ind->getGenome(_encoding));
  GA1DArrayAlleleGenome<long double>& centroid_real  = dynamic_cast<GA1DArrayAlleleGenome<long double>&>(  *centroid.getGenome(_encoding));
  GA1DArrayAlleleGenome<long double>& worst_ind_real = dynamic_cast<GA1DArrayAlleleGenome<long double>&>( *worst_ind.getGenome(_encoding));

  for (unsigned dim = 0; dim < contr_ind_real.length(); dim++)
    contr_ind_real.gene(dim, centroid_real.gene(dim) - gamma * (centroid_real.gene(dim) - worst_ind_real.gene(dim)));

  // Do not forget to evaluate the genome (as we work directly on the inner genome)
  contr_ind->evaluate(gaTrue);

  long double bestParentFit = worst_ind.score();
  contr_ind->setFitnessIncrement(contr_ind->computeFitnessIncrement(bestParentFit));
  contr_ind->mustComputeQuality(true);

  return contr_ind;
}


void MOSTechniqueNMA::shrinkPop(const GAPopulation& pop, unsigned& used_evals, unsigned max_evals, long double sigma) const {
  // /*LOG*/ GALogger::instance()->appendPopulation("NMA:shrink", "Before shrink pop is", pop);

  MOSGenome&                     best_ind      = dynamic_cast<MOSGenome&>(pop.individual(0));
  GA1DArrayAlleleGenome<long double>& best_ind_real = dynamic_cast<GA1DArrayAlleleGenome<long double>&>(*best_ind.getGenome(_encoding));

  for (int ind_pos=1; ind_pos<pop.size() && used_evals <= max_evals; ind_pos++, used_evals++) { // First pos = 1 not an error

    MOSGenome&                     tmp_ind      = dynamic_cast<MOSGenome&>(pop.individual(ind_pos));
    GA1DArrayAlleleGenome<long double>& tmp_ind_real = dynamic_cast<GA1DArrayAlleleGenome<long double>&>(*tmp_ind.getGenome(_encoding));

    long double bestParentFit = tmp_ind.score();

    for (unsigned dim = 0; dim < tmp_ind_real.length(); dim++)
      tmp_ind_real.gene(dim, best_ind_real.gene(dim) + sigma * (tmp_ind_real.gene(dim) - best_ind_real.gene(dim)));

    // Do not forget to evaluate the genome (as we work directly on the inner genome)
    tmp_ind.evaluate(gaTrue);

    tmp_ind.setFitnessIncrement(tmp_ind.computeFitnessIncrement(bestParentFit));
    tmp_ind.mustComputeQuality(true);

  }

}


void MOSTechniqueNMA::contraction(const MOSGenome& centroid, const MOSGenome& reflx_ind, MOSGenome& worst_ind, const GAPopulation& pop, unsigned& used_evals, unsigned max_evals, unsigned& branch, long double gamma, long double sigma, stringstream& msg) const {
  // Contraction between centroid and the best between worst and next_worst
  MOSGenome* contr_ind = NULL;

  // In/Out contraction
  used_evals++;

  if (GAGenome::compareScores(reflx_ind.score(), worst_ind.score()) == GAGenome::BETTER ) {

    // Outside Contraction
    contr_ind = outsideContraction(centroid, reflx_ind, gamma);

    // /*LOG*/ msg << "Outside contraction" << endl;

    if (GAGenome::compareScores(contr_ind->score(), reflx_ind.score()) == GAGenome::BETTER ||
        GAGenome::compareScores(contr_ind->score(), reflx_ind.score()) == GAGenome::EQUAL     ) {
      worst_ind.copy(*contr_ind);
      branch = MOSTechniqueNMA::ContractionOutOK;
    }
    else {
      shrinkPop(pop, used_evals, max_evals, sigma);
      branch = MOSTechniqueNMA::Shrink;
    }

  }
  else  {

    // Inside contraction
    contr_ind = insideContraction(centroid, worst_ind, gamma);

    // /*LOG*/ msg << "Inside contraction" << endl;

    if (GAGenome::compareScores(contr_ind->score(), worst_ind.score()) == GAGenome::BETTER) {
      worst_ind.copy(*contr_ind);
      branch = MOSTechniqueNMA::ContractionInOK;
    }
    else {
      shrinkPop(pop, used_evals, max_evals, sigma);
      branch = MOSTechniqueNMA::Shrink;
    }

  }

  // /*LOG*/ msg << " Contracted ind: " << GALogger::instance()->getIndLogStr(*contr_ind) << endl;
  // /*LOG*/ msg << " Worst new ind:"   << GALogger::instance()->getIndLogStr(worst_ind)  << endl;

  delete contr_ind;
}

bool MOSTechniqueNMA::NMA (GAPopulation& pop, /* out */ unsigned& used_evals, unsigned max_evals, unsigned& branch, long double rho, long double chi, long double gamma, long double sigma) {
  if (used_evals >= max_evals) return false;

  // Assuming a descending order, first pos is the best ind last pos is the worst ind
  assert(GAGenome::compareScores(pop.individual(0).score(), pop.individual(pop.size()-1).score()) == GAGenome::BETTER ||
         GAGenome::compareScores(pop.individual(0).score(), pop.individual(pop.size()-1).score()) == GAGenome::EQUAL);
  // /*LOG*/ GALogger::instance()->appendPopulation("NMA", "Before doing NMA pop is", pop);

  MOSGenome* centroid        = centroidOf(pop);
  MOSGenome& worst_ind       = dynamic_cast<MOSGenome&>(pop.individual(pop.size()-1)); // worst ind
  MOSGenome* reflx_ind       = reflection(*centroid, worst_ind, rho);
  MOSGenome& next_worst_ind  = dynamic_cast<MOSGenome&>(pop.individual(pop.size()-2)); // next worst ind
  long double     old_worst_score = worst_ind.score();

  used_evals++; // For reflection

     /*LOG*/ stringstream msg; msg << " Start NMA iteration used_evals = " << used_evals << " max evals=" << max_evals << endl;
  // /*LOG*/ msg << " rho= " << rho << " chi=" << chi << " gamma=" << gamma << " sigma=" << sigma << endl;
  // /*LOG*/ msg << " Centroid: "     << GALogger::instance()->getIndLogStr(*centroid)      << endl;
  // /*LOG*/ msg << " Worst: "        << GALogger::instance()->getIndLogStr(worst_ind )     << endl;
  // /*LOG*/ msg << " Reflx: "        << GALogger::instance()->getIndLogStr(*reflx_ind)     << endl;
  // /*LOG*/ msg << " Next Worst: "   << GALogger::instance()->getIndLogStr(next_worst_ind) << endl;

  if ((GAGenome::compareScores(reflx_ind->score(), pop.best().score()    ) == GAGenome::WORSE ||
       GAGenome::compareScores(reflx_ind->score(), pop.best().score()    ) == GAGenome::EQUAL    ) &&
       GAGenome::compareScores(reflx_ind->score(), next_worst_ind.score()) == GAGenome::BETTER        ) {
    // /*LOG*/ msg << " Accepting first reflected ind since f0 <= fr < fn" << endl;
    worst_ind.copy(*reflx_ind);
    branch = MOSTechniqueNMA::Reflection;
  }
  else if (GAGenome::compareScores(reflx_ind->score(), pop.individual(0).score()) == GAGenome::BETTER && used_evals <= max_evals) {
    // /*LOG*/ msg << " Checking Expansion due to fr < f0 " << endl;
    expansion(*centroid, *reflx_ind, worst_ind, chi, msg);
    branch = MOSTechniqueNMA::Expansion;
    used_evals++;
  }
  else if ((GAGenome::compareScores(reflx_ind->score(), next_worst_ind.score()) == GAGenome::WORSE ||
            GAGenome::compareScores(reflx_ind->score(), next_worst_ind.score()) == GAGenome::EQUAL    ) && used_evals <= max_evals ) {
    contraction(*centroid, *reflx_ind, worst_ind, pop, used_evals, max_evals, branch, gamma, sigma,msg);
  }

  // Final ordering, order of equivalents should be mantained. A Stable sort operation should be created in GAPopulation
  pop.sort(gaTrue);

  delete centroid;
  delete reflx_ind;

  // /*LOG*/ msg << "End NMA iteration used_evals=" << used_evals << " max evals=" << max_evals << endl;
  // /*LOG*/ GALogger::instance()->appendLogMessage("NMA", msg.str());
  // /*LOG*/ GALogger::instance()->appendPopulation("NMA", "after doing NMA pop is", pop);

  return ( GAGenome::compareScores(pop.individual(pop.size()-1).score(), old_worst_score) == GAGenome::BETTER );
}
