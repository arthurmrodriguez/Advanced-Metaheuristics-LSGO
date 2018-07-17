/**
 * @file
 * @brief MOSTechniqueLS class impl.
 *
 */

#include "MOSTechniqueLS.h"

#include "garandom.h"
#include "LSElitism.h"
#include "MOSConversion.h"
#include "genomes/GA1DArrayGenome.h"
#include "genomes/MOSGenome.h"

/**
 * MOSTechniqueLS Constructor
 */
MOSTechniqueLS::MOSTechniqueLS(techIdType id, std::string description, LocalSearch ls, GAGenome::Comparator comparator,
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

  _ls          = ls;
  _comparator  = comparator;

  _recombinator = new LSElitism();

  GA1DArrayAlleleGenome<long double>& g = dynamic_cast<GA1DArrayAlleleGenome<long double>&>(*genomeBase);
  _SR = (long double) (g.alleleset(0).upper() - g.alleleset(0).lower()) / 2.0;
  _currentGenome = NULL;

  return;

}


unsigned MOSTechniqueLS::offspring (unsigned maxEvals, unsigned& usedEvals, GAPopulation* pop, GAPopulation* auxPop, bool& converged) {

  if (!_currentGenome)
    _currentGenome = dynamic_cast<MOSGenome*>(pop->best().clone());

  // Do not forget to assign the individual to this technique
  _currentGenome->setTechnique(this);

  MOSGenome& g = dynamic_cast<MOSGenome&>(auxPop->individual(0));
  g.copy(pop->best());

  // Do not forget to assign the individual to this technique
  g.setTechnique(this);

  if (GAGenome::compareScores(g.score(), _currentGenome->score()) < GAGenome::BETTER)
    g.copy(*_currentGenome);

  // origFitBest is used to compute the fitness increment after the execution of the LS.
  // fitBest will be updated within the LS and will store the final fitness value of the
  // selected individual
  long double origFitBest = g.score();
  long double fitBest = origFitBest;

  long double fit_inc_acum_total = 0.0;

  // Iterate until the maximum number of FEs for this technique in this generation has
  // been exceeded
  do {
    long double fit_inc_acum = 0.0;
    unsigned evals = 0;
    long double grade = (*_ls)(*g.getGenome(_encoding), fitBest, evals, fit_inc_acum, _times_improved, _SR, maxEvals-usedEvals);
    g.score(g.getGenome(_encoding)->score());
    fit_inc_acum_total += fit_inc_acum;
    usedEvals+=evals;
  } while (usedEvals < maxEvals && !converged);

  // Compute fitness increment and store it in the modified individual.
  long double fit_inc = fit_inc_acum_total / (long double) usedEvals;
  g.setFitnessIncrement(fit_inc_acum_total);
  g.mustComputeQuality(true);

  _currentGenome->copy (g);

  // Remind to update stats
  _stats.numsel+=1;
  _stats.numrep+=1;
  _stats.nummut+=usedEvals;
  _stats.numcro+=0;
  _stats.numeval+=usedEvals;

  return 1;

}
