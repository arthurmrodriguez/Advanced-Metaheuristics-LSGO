#ifndef EABESTIMMIGRANTSSELECTOR_H_
#define EABESTIMMIGRANTSSELECTOR_H_

#include "EAImmigrantsSelector.h"

class EABestImmigrantsSelector : public EAImmigrantsSelector
{
public:
	EABestImmigrantsSelector(GAGeneticAlgorithm& ea);
	~EABestImmigrantsSelector();
	
	void admitImmigrants(GAPopulation& immigrants);
};

#endif /*EABESTIMMIGRANTSSELECTOR_H_*/
