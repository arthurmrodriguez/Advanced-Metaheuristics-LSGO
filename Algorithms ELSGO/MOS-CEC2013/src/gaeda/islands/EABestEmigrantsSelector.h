#ifndef EABESTEMIGRANTSSELECTOR_H_
#define EABESTEMIGRANTSSELECTOR_H_

#include "EAEmigrantsSelector.h"

class EABestEmigrantsSelector : public EAEmigrantsSelector {
public:
	EABestEmigrantsSelector(GAGeneticAlgorithm& ea, unsigned int mig_popsize);
	~EABestEmigrantsSelector();

	void setEmigrantsPop(const GAPopulation& pop);
};

#endif /*EABESTEMIGRANTSSELECTOR_H_*/
