#ifndef EAISLANDSMODEL_H
#define EAISLANDSMODEL_H

#include "EAIslandsTopology.h"
#include "EAEmigrantsSelector.h"
#include "EAImmigrantsSelector.h"
#include "CommManager.h"
#include "../genomes/GAGenome.h"
#include "../GAGeneticAlgorithm.h"
#include "../Algorithm.h"


/**
 * @brief Island GA model 
 *
 * This is an abstract base class that allows building GA island models
 */

class EAIslandsModel : public Algorithm {
protected:
  /*
   * All the pointer are been aquired by the objects from this class, i.e. their destructor
   * is called with the destructor of this class 
   */
  EAIslandsTopology*    topology_; 
  GAGeneticAlgorithm*   ea_;       
  unsigned int          mig_period_;
  EAEmigrantsSelector*  emigrants_selector_;     
  EAImmigrantsSelector* immigrants_selector_;
  CommManager&          comm_manager_;
  
  virtual void evolve     ()=0;
  virtual void step       ();
  virtual void initialize ();
  
public:
  GADefineIdentity ("EAIslandsModel", GAID::EAIslandsModel);

  EAIslandsModel(GAGeneticAlgorithm*   ea, 
                 EAIslandsTopology*    top, 
                 unsigned int          mig_period, 
                 EAEmigrantsSelector*  emigrants_selector,
                 EAImmigrantsSelector* immigrants_selector,
                 CommManager&          comm_manager);

  virtual ~EAIslandsModel(void);

  GAGenome*           best       (               ) const;
  const GAStatistics& statistics (               ) const;
  const GAPopulation& population (               ) const;
  const GAGeneticAlgorithm* getEA (    ) const{return ea_;}
};

#endif
