/* INCLUDES */
#include "EAIslandsModel.h"
#include "../GAParameter.h"
#include "../gaparameters.h"
#include "../logger/GALogger.h"

#include <assert.h>
#include <memory>

EAIslandsModel::EAIslandsModel(GAGeneticAlgorithm*   ea, 
                               EAIslandsTopology*    top, 
                               unsigned int          mig_period, 
                               EAEmigrantsSelector*  emigrants_selector,
                               EAImmigrantsSelector* immigrants_selector,
                               CommManager&          comm_manager): 
                                                                  topology_(top), 
                                                                  ea_(ea), 
                                                                  mig_period_(mig_period), 
                                                                  emigrants_selector_(emigrants_selector),
                                                                  immigrants_selector_(immigrants_selector),
                                                                  comm_manager_(comm_manager){}

EAIslandsModel::~EAIslandsModel() {
	delete topology_;
	delete ea_;
	delete emigrants_selector_;
	delete immigrants_selector_;
}

void EAIslandsModel::initialize() {
  ea_->initialize();
}

GAGenome* EAIslandsModel::best() const{
  GAGenome* best_ind = NULL;
  // At the end of the algorithm, each slave island sends its best result to the master island that recollects the individuals 
  if ( comm_manager_.isIslandSlave() ) {             
    comm_manager_.sendIndToMaster(ea_->population().best());
  } else {
    auto_ptr< vector<GAGenome*> > inds ( comm_manager_.receiveOneIndFromAllSlaveIslands() );
    GAPopulation pop(*inds);
    pop.add( ea_->population().best() );
    pop.sort();
    best_ind = pop.best().clone();
  }
  return best_ind;
}

void EAIslandsModel::step() {
  ea_->step();
}

const GAStatistics& EAIslandsModel::statistics () const {return ea_->statistics();}
const GAPopulation& EAIslandsModel::population() const{ return ea_->population();}
