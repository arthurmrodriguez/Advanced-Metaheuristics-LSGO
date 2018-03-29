#include "EAEmigrantsSelector.h"

/**
 * For effiency the emigrants pop is created only once so the new operators that extend this class
 * need only to copy the genome values to each individual of this emmigrant pop
 */
EAEmigrantsSelector::EAEmigrantsSelector(GAGeneticAlgorithm& ea, unsigned int mig_popsize) :
                                               ea_(ea), emigrants_(ea_.population().individual(0),mig_popsize){
  assert(mig_popsize>=0);
  assert(mig_popsize <= ea_.population().size() );
}


EAEmigrantsSelector::~EAEmigrantsSelector(){}

GAPopulation& EAEmigrantsSelector::getEmigrantsPop(){
  setEmigrantsPop(ea_.population());
  return emigrants_;
}
