#include "EABestEmigrantsSelector.h"

EABestEmigrantsSelector::EABestEmigrantsSelector(GAGeneticAlgorithm& ea, unsigned int mig_popsize): EAEmigrantsSelector(ea,mig_popsize){}

EABestEmigrantsSelector::~EABestEmigrantsSelector(){}

void EABestEmigrantsSelector::setEmigrantsPop(const GAPopulation& pop){
  pop.sort();	
  for (unsigned i=0; i<emigrants_.size(); i++) {
    emigrants_.individual(i).copy( pop.individual(i) );
  }
}
