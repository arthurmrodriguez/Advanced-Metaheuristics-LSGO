#ifndef EAINMIGRANTSSELECTOR_H_
#define EAINMIGRANTSSELECTOR_H_

#include "../GAGeneticAlgorithm.h"
#include "../GAPopulation.h"

class EAEmigrantsSelector {
	GAGeneticAlgorithm& ea_;           // This is neccesary since the GA has two populations and changes the actual pop from
	                                   // this two during evolution
protected:
	GAPopulation        emigrants_;
public:
	EAEmigrantsSelector(GAGeneticAlgorithm& ea, unsigned int mig_popsize);
	virtual ~EAEmigrantsSelector();
	
	GAPopulation& getEmigrantsPop();
	virtual void  setEmigrantsPop(const GAPopulation& pop )=0; // This method is the one that needs to be redeclared in order to define new emigrants selectors
};

#endif /*EAINMIGRANTSSELECTOR_H_*/
