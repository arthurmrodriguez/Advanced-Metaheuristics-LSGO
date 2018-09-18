#include "MOSConfig.h"

#include "MOSTechniqueSet.h"
#include "GAEDAConfig.h"

#include <assert.h>

/**
 * mixProbVector
 * Computes the average value of two probability vectors (from parents)
 */
void mixProbVector(const std::vector< const MOSProbVector* >& probs, const std::vector<double>& scores, MOSProbVector& pvChild) {

  bool weighted = GAEDAConfig::handle()->weightedAutonomic();

  MOSTechniqueSet* techniqueSet = MOSTechniqueSet::handle();
  MOSTechniqueSet::MOSTechniqueSetIterator it;
  MOSProbVector::const_iterator itp;

  double quotient = weighted ? 1.0 : probs.size();

  for (it=techniqueSet->begin(); it!=techniqueSet->end(); it++) {

    double sum_tech = 0.0;

    for (unsigned i = 0; i < scores.size(); i++)
      sum_tech += scores[i];

    pvChild[it->first] = 0.0;

    for (unsigned i = 0; i < probs.size(); i++) {
      double factor = weighted ? (scores[i] / sum_tech) : 1.0;

      itp = probs[i]->find(it->first);
      assert (itp != probs[i]->end());

      pvChild[it->first] = pvChild[it->first] + (factor * itp->second / quotient);
    }

  }

  return;

}
