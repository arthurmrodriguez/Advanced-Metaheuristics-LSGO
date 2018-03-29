/* INCLUDES */
#include "EAIslandsModelSync.h"
#include "../logger/GALogger.h"
//#include "VoronoiIndInit.h"
#include <list>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <sstream>
#include <memory>

using namespace std; // For the STL stuff

// Auxiliary functions declaration
void deleteMedoids(vector<GAGenome*>* medoids);


EAIslandsModelSync::EAIslandsModelSync(GAGeneticAlgorithm*   ea,        
                                       EAIslandsTopology*    top,   
                                       unsigned int          mig_period,
                                       EAEmigrantsSelector*  emigrants_selector,
                                       EAImmigrantsSelector* immigrants_selector,
                                       CommManager&          comm_manager) :  
	                                                                    EAIslandsModel(ea, 
	                                                                    		             top, 
	                                                                    		             mig_period, 
      	                                                                    		       emigrants_selector,
      	                                                                    		       immigrants_selector,
      	                                                                    		       comm_manager){}

EAIslandsModelSync::~EAIslandsModelSync (){}

/* 
 * This method implements the sync algorithm. In a future, it should be refactored in order to take all the
 * common parts with the rest of the algorithms to the base class. This way, the base class should define the
 * main template with hook methods like prestepActions() poststepActions() and postIslandConvergedActions()
 */
void EAIslandsModelSync::evolve() {
  //ea_->printStats("Initial Stats");
  
  bool all_islands_converged = false;

  //TODO: refactor as in the DistRoutingAlg
  while (!ea_->done() || !all_islands_converged){
    /*LOG*/ GALogger::instance()->appendPopulation( "GAIslandsModelSync::run", "Before entering inner step", ea_->population() );

    if ( !ea_->done() ) {
      do {
        ea_->step();                                                     

        postStep();                                                             // Necessary for some derived models (cluster pop model). It doesnt do anything in normal scenario
        ea_->population().sort();                                               // Only meaningful if postStep changes the population
        //ea_->printStats("End of step");
      } while ( ea_->statistics().generation() % mig_period_ != 0);        
    } 
    else {
     ea_->updateNoStepsStats();                                                 // To update the number of gens and other stats
    }
    topology_->generateNewNeighbors( ea_->population() );                       // Only meaningful with dynamic topologies                                 
    
    doMigration();
    all_islands_converged = haveAllIslandsConverged();
  }
}

void EAIslandsModelSync::doMigration(){
  GAPopulation& emigrant_pop = emigrants_selector_->getEmigrantsPop();              
  
  auto_ptr< list<int>    > receivers_pos ( topology_->getDestinations() );
  assert (topology_->getDestinations()->size() > 0);
  auto_ptr< GAPopulation > inmigrants ( comm_manager_.sendAndReceiveInds(emigrant_pop, *receivers_pos) ); // sync sending
  
  immigrants_selector_->admitImmigrants(*inmigrants);             

  /*LOG*/  GALogger::instance()->appendPopulation("GAIslandsModelSync::run. ", "Population after migration", ea_->population() );
}

bool EAIslandsModelSync::haveAllIslandsConverged(){
  auto_ptr< list<int> > converged_islands_pos ( comm_manager_.sendAndReceiveConvergence( ea_->done() ) );

  {/*LOG*/ stringstream message; message << "n# converged islands=" << converged_islands_pos->size();
    GALogger::instance()->appendLogMessage("Convergence", message.str());
  }

  return ( (int) converged_islands_pos->size() == comm_manager_.getNumIslands() - 1 );
}



void EAIslandsModelSync::postStep(){}


