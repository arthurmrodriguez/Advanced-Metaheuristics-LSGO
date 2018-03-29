/**
 * @file
 * @brief MOSTechniqueLS class impl.
 *
 */

#include "MOSTechniqueLS.h"

#include "garandom.h"
#include "LSElitism.h"
// #include "GAGenealogy.h"
// #include "GAGenealogyMemory.h"
 #include "MOSConversion.h"
// #include "MOSTechniqueSet.h"
// #include "logger/GALogger.h"
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


MOSTechniqueLS::~MOSTechniqueLS() {
  if (_currentGenome)
    delete _currentGenome;
}


/**
 * Creates offspring for this technique from input population (used in central approach)
 * @param origPop Input population
 * @param destPop Output population
 * @param size Number of individuals to create
 * @param offset Offset to store individuals in ouput population
 */
void MOSTechniqueLS::offspring(GAPopulation* origPop, GAPopulation* destPop, unsigned size, unsigned offset) {

  if (!_currentGenome)
    _currentGenome = dynamic_cast<MOSGenome*>(origPop->best().clone());

  // Do not forget to assign the individual to this technique
  _currentGenome->setTechnique(this);

  // Initialization of selector. This is needed as populations are continuously swapped
  // within the main algorithm and the selector object stores a reference to the population
  // it works on.
  _selector->assign(*origPop);
  _selector->update();

  MOSGenome& g = dynamic_cast<MOSGenome&>(origPop->best());

  // Do not forget to assign the individual to this technique
  g.setTechnique(this);

  // Exchange the genome to me modified with the first assigned genome in the output
  // population. This way we avoid unnecesary copies and always apply the LS to the
  // same individual(s).
  GAGenome* tmp = &destPop->individual(offset);
  tmp->setWorstScore();
  tmp = destPop->replace(&g, tmp);
  tmp = origPop->replace(tmp, &g);

  if (GAGenome::compareScores(g.score(), _currentGenome->score()) < GAGenome::BETTER)
    g.copy(*_currentGenome);

  // origFitBest is used to compute the fitness increment after the execution of the LS.
  // fitBest will be updated within the LS and will store the final fitness value of the
  // selected individual
  long double origFitBest = g.score();
  long double fitBest = origFitBest;

  if (g.mustBeNulled())
    g.mustBeNulled(false);

  unsigned totalEvals = 0;
  long double fit_inc_acum_total = 0.0;

  _times_improved = 0;

  // Iterate until the maximum number of FEs for this technique in this generation has
  // been exceeded
  do {
    long double fit_inc_acum = 0.0;
    unsigned evals = 0;
    long double grade = (*_ls)(*g.getGenome(_encoding), fitBest, evals, fit_inc_acum, _times_improved, _SR, size-totalEvals);
    g.score(g.getGenome(_encoding)->score());
    fit_inc_acum_total += fit_inc_acum;
    totalEvals+=evals;
  } while (totalEvals < size);

  // Compute fitness increment and store it in the modified individual.
  // The fitness increment is scaled within [0, 1], as the first computation
  // returns a value within [-1, 1]
  long double fit_inc = fit_inc_acum_total / (long double) totalEvals;
  g.setFitnessIncrement(fit_inc);

  // Mark the rest of the offspring population with zero values, so that these individuals
  // are not selected by the elitism. Additionally, we assign to every individual the same
  // fitness increment value, so that the quality value of the whole population could be
  // computed properly
  for (unsigned i = 1; i < size; i++) {
    MOSGenome& g2 = dynamic_cast<MOSGenome&>(destPop->individual(offset+i));
    g2.purgeGenome(this);
    g2.setFitnessIncrement(fit_inc);
    destPop->individual(offset+i).setWorstScore();
  }

  _currentGenome->copy (g);

  // Remind to update stats
  _stats.numsel+=1;
  _stats.numrep+=1;
  _stats.nummut+=totalEvals;
  _stats.numcro+=0;
  _stats.numeval+=totalEvals;

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
    if (g.precissionReached()) converged = true;
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
