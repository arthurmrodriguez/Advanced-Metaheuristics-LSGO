#include "MOSConfig.h"

#include "MOSTechniqueSet.h"
#include "GAEDAConfig.h"

#include <assert.h>

/**
 * mixProbVector
 * Computes the average value of two probability vectors (from parents)
 */
void mixProbVector(const std::vector< const MOSProbVector* >& probs, const std::vector<long double>& scores, MOSProbVector& pvChild) {

  MOSTechniqueSet* techniqueSet = MOSTechniqueSet::handle();
  MOSTechniqueSet::MOSTechniqueSetIterator it;
  MOSProbVector::const_iterator itp;

  for (it=techniqueSet->begin(); it!=techniqueSet->end(); it++) {

    long double sum_tech = 0.0;

    for (unsigned i = 0; i < scores.size(); i++)
      sum_tech += scores[i];

    pvChild[it->first] = 0.0;

    for (unsigned i = 0; i < probs.size(); i++) {
      itp = probs[i]->find(it->first);
      assert (itp != probs[i]->end());

      pvChild[it->first] = pvChild[it->first] + (itp->second / probs.size());
    }

  }

  return;

}
