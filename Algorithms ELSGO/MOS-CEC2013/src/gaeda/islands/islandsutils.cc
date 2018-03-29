#include "islandsutils.h"

vector<GAGenome*>* convertPopToVector(GAPopulation* pop){
  vector<GAGenome*>*  inds = new vector<GAGenome *>; 
  for (unsigned i=0; i<pop->size(); i++) inds->push_back( pop->individual(i).clone() ); 
  return inds;
}
