#ifndef ISLANDSMODELSYNC_H
#define ISLANDSMODELSYNC_H

#include "EAIslandsModel.h"
#include "EAIslandsTopology.h"
#include "../gatypes.h"
#include "../GAGeneticAlgorithm.h"
#include "../GAPopulation.h"

/**
 * @brief Islands EA model (synchronous communications)
 *
 * Island GA model that uses synchronous communications. This model is implemented
 * using MPI and assumes that the MPI subsystem has already been initialized
 */

class EAIslandsModelSync : public EAIslandsModel  {

  void doMigration             ();
  bool haveAllIslandsConverged ();

protected:
  virtual void postStep();
  virtual void evolve  ();
  
public:
  EAIslandsModelSync( GAGeneticAlgorithm*   ea,        
                      EAIslandsTopology*    top,   
                      unsigned int          mig_period,
                      EAEmigrantsSelector*  emigrants_selector,
                      EAImmigrantsSelector* immigrants_selector,
                      CommManager&          comm_manager);
		    
  virtual ~EAIslandsModelSync () ;
};

#endif
