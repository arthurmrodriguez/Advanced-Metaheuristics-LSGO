#ifndef EARANDOMEMIGRANTSSELECTOR_H_
#define EARANDOMEMIGRANTSSELECTOR_H_

#include "EAEmigrantsSelector.h"

class EARandomEmigrantsSelector : public EAEmigrantsSelector
{
public:
	EARandomEmigrantsSelector(GAGeneticAlgorithm& ea, unsigned int mig_popsize);
	~EARandomEmigrantsSelector();
	
	void setEmigrantsPop(const GAPopulation& pop);
};

#endif /*EARANDOMEMIGRANTSSELECTOR_H_*/
