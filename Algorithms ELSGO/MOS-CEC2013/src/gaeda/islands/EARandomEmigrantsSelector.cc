#include "EARandomEmigrantsSelector.h"
#include <stdlib.h>

EARandomEmigrantsSelector::EARandomEmigrantsSelector(GAGeneticAlgorithm& ea, unsigned int mig_popsize): EAEmigrantsSelector(ea,mig_popsize){}

EARandomEmigrantsSelector::~EARandomEmigrantsSelector(){}

void EARandomEmigrantsSelector::setEmigrantsPop(const GAPopulation& pop){
  for (unsigned i=0; i<emigrants_.size(); i++) {
    int random_pos = rand() % pop.size();
    emigrants_.individual(i).copy( pop.individual(random_pos) );
  }
}
