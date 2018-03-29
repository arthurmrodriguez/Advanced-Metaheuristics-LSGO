/*
 * =====================================================================================
 * 
 *       Filename:  RandomCentGenerator.h
 * 
 *    Description:  This class implements a centroids generator that creates N random
 *                  centroids distributions and gets the one that maximizes the sum of
 *                  the distances between the centroids
 * 
 *        Version:  1.0
 *        Created:  02/28/07 11:26:21 CET
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:   Santiago Muelas, 
 *        Company:  
 * 
 * =====================================================================================
 */

#ifndef RANDOM_CENTROIDS_GENERATOR__H
#define RANDOM_CENTROIDS_GENERATOR__H 

#include "CentroidsGenerator.h"

class RandomCentGenerator : public CentroidsGenerator{

protected:
  vector<GAGenome*>* createARandomCentroidsDistribution() const;
public:
  
  RandomCentGenerator(unsigned int cent_number, GAGenome& sample_gen, GAGenome::Initializer pinit_func);
  virtual ~RandomCentGenerator();

  virtual vector<GAGenome*>* generateCentroids() const;
};

#endif

