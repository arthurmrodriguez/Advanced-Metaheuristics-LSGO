#include "RandomCentGenerator.h"

//Methods
RandomCentGenerator::RandomCentGenerator(unsigned int cent_number, GAGenome& sample_gen, GAGenome::Initializer pinit_func) 
                                                  : CentroidsGenerator(cent_number,sample_gen, pinit_func){}
vector<GAGenome*>* RandomCentGenerator::generateCentroids() const {
  vector<GAGenome*>* centroids = createARandomCentroidsDistribution();
 return centroids; 
}

RandomCentGenerator::~RandomCentGenerator(){}

vector<GAGenome *>* RandomCentGenerator::createARandomCentroidsDistribution() const{
  vector<GAGenome*>* centroids = new vector<GAGenome*>(cent_number);
  for (unsigned int i=0; i<cent_number; i++){
    GAGenome* gen_temp = sample_gen.clone();
    init_func(*gen_temp);
    (*centroids)[i] = gen_temp;
  }
  return centroids;
}

