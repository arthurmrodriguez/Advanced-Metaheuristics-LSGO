#ifndef EAIMMIGRANTSSELECTOR_H_
#define EAIMMIGRANTSSELECTOR_H_

#include "../GAGeneticAlgorithm.h"
#include "../GAPopulation.h"

class EAImmigrantsSelector
{
protected:
	GAGeneticAlgorithm& ea_;           // This is neccesary since the GA has two populations and changes the actual pop from
	                                   // this two during evolution
public:
	EAImmigrantsSelector(GAGeneticAlgorithm& ea);
	virtual ~EAImmigrantsSelector();
	
	virtual void admitImmigrants(GAPopulation& immigrants)=0;
};

#endif /*EAIMMIGRANTSSELECTOR_H_*/
