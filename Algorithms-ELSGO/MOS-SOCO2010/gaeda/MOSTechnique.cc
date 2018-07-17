/**
 * @file
 * @brief MOSTechnique class impl.
 *
 */
#include "MOSTechnique.h"

#include "MOSQualityFunction.h"

bool MOSTechnique::improvement_override = false;

MOSTechnique::MOSTechnique () : _times_improved (0), _evals_last_gen (0) {}

long double MOSTechnique::evolve(GAPopulation*& pop, unsigned maxEvals, MOSQuality* qualityFunction, bool& converged) {

  unsigned totalEvals = 0;
  unsigned iters      = 0;
  long double   qualAcum   = 0.0;
  GAPopulation* auxPop = pop->clone();

  // To allow the computation of the ratio of improvements
  _evals_last_gen = maxEvals;

  _times_improved = 0;

  long double total_improvements = 0;

  while (totalEvals < maxEvals && !converged) {
    // Initialization of selector. This is needed as populations are continuously swapped
    // within the main algorithm and the selector object stores a reference to the population
    // it works on.
    _selector->assign(*pop);
    _selector->update();

    unsigned usedEvals = 0;
    unsigned newInds = offspring(maxEvals - totalEvals, usedEvals, pop, auxPop, converged);

    totalEvals += usedEvals;
    long double qualAct = qualityFunction->update(*auxPop, _id, usedEvals);
    qualAcum += qualAct;

    iters++;
    total_improvements += _times_improved;

    _recombinator->recombine(*pop, *auxPop);

    auxPop->sort();
    auxPop->scale();
    std::swap(pop, auxPop);
  }

  delete auxPop;

  _times_improved = total_improvements;

  return (qualAcum / (long double) iters);

}

