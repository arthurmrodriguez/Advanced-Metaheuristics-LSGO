#include "EABestImmigrantsSelector.h"
#include "../logger/GALogger.h"
#include <assert.h>

EABestImmigrantsSelector::EABestImmigrantsSelector(GAGeneticAlgorithm& ea): EAImmigrantsSelector(ea) {}

EABestImmigrantsSelector::~EABestImmigrantsSelector(){}

void EABestImmigrantsSelector::admitImmigrants(GAPopulation& immigrants){
  if (immigrants.size() > 0) {
    GAPopulation& pop               = const_cast<GAPopulation &> (ea_.population());
    int           original_pop_size = pop.size();     
    
    /*LOG*/  GALogger::instance()->appendPopulation("EABestImmigrantsSelector::run. ", "Population before immigrants addition", ea_.population() );
    for (unsigned i=0; i<immigrants.size(); i++) pop.add( immigrants.individual(i) );  
  
    /*LOG*/  GALogger::instance()->appendPopulation("EABestImmigrantsSelector::run. ", "Population with immigrants", ea_.population() );
    pop.sort(gaTrue);
  
    ea_.populationSize(original_pop_size);
    /*LOG*/  GALogger::instance()->appendPopulation("EABestImmigrantsSelector::run. ", "Population resized", ea_.population() );
  
    // Update stats
    GAStatistics& stat = const_cast<GAStatistics &> (ea_.statistics());
    pop.statistics(gaTrue);
    stat.setScore(pop);	
  }
}
