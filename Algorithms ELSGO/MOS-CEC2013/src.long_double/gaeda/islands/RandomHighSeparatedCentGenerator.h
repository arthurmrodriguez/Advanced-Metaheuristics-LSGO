#ifndef RANDOM_HIGH_SEP_CENTROIDS_GENERATOR__H
#define RANDOM_HIGH_SEP_CENTROIDS_GENERATOR__H 

#include "RandomCentGenerator.h"

/**
 * This class implements a centroids generator that creates N random
 *                  centroids distributions and gets the one that maximizes the sum of
 *                  the distances between the centroids
 */
class RandomHighSeparatedCentGenerator : public RandomCentGenerator {
  unsigned int max_tries;

  long double sumMinDistancesBetweenCentroids(vector<GAGenome*>& centroids) const;
  long double findMinDistanceToOtherCentroids(vector<GAGenome*>& centroids, unsigned int pos) const;
public:
  RandomHighSeparatedCentGenerator(unsigned int cent_number, GAGenome& sample_gen, GAGenome::Initializer pinit_func, unsigned int max_tries=1000);
  ~RandomHighSeparatedCentGenerator();

  vector<GAGenome*>* generateCentroids() const;
};

#endif

