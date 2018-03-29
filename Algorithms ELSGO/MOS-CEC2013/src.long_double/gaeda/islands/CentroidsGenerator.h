#ifndef CENTROIDS_GENERATOR__H
#define CENTROIDS_GENERATOR__H 

#include "../genomes/GAGenome.h"
#include <vector>

using namespace std;

/**
 * This is an abstract class which will be extended to implement the 
 *                  
 */
class CentroidsGenerator{
protected:
  unsigned int          cent_number;
  GAGenome&             sample_gen;
  GAGenome::Initializer init_func;
public:
  CentroidsGenerator(unsigned int cent_number,GAGenome& sample_gen, GAGenome::Initializer init_func);
  virtual ~CentroidsGenerator();
  virtual vector<GAGenome*>* generateCentroids() const = 0; 
};

#endif
