#ifndef EAISLANDSMODELCLUSTEREDPOPS_H
#define EAISLANDSMODELCLUSTEREDPOPS_H

#include "EAIslandsModelSync.h"
#include "CentroidsGenerator.h"

class EAIslandsModelClusteredPops : public EAIslandsModelSync {
  typedef vector<GAPopulation*>* (*cluster_func_t)          (GAPopulation&       pop, 
                                                             unsigned int        num_islands, 
                                                             unsigned int        island_pop_size,
                                                             CentroidsGenerator* centroids_gen);
  typedef void                   (*popssizes_resize_func_t) (vector<GAPopulation*>& pops, unsigned int island_pop_size);
  
  vector<GAPopulation*>* obtainNewIslandsPops(GAPopulation& pop);

protected:  
  unsigned int            clustering_period_;
  cluster_func_t          cluster_func_;
  popssizes_resize_func_t popsizes_resize_func_;
  CentroidsGenerator*     centroids_gen_;     
  
  virtual void postStep();  
public:
  EAIslandsModelClusteredPops( GAGeneticAlgorithm*     ga,        
                               EAIslandsTopology*      top,   
                               unsigned int            mig_period,
                               EAEmigrantsSelector*    emigrants_selector,
                               EAImmigrantsSelector*   immigrants_selector,
                               unsigned int            clustering_period,
                               cluster_func_t          cluster_func,
                               popssizes_resize_func_t popsizes_resize_func,
                               CentroidsGenerator*     centroids_gen,         // It is going to be deleted in the destructor
                               CommManager&            comm_manager);
        
  virtual ~EAIslandsModelClusteredPops () ;  

  static vector<GAPopulation*>* kmeansCluster          (GAPopulation&       pop, 
                                                        unsigned int        num_islands, 
                                                        unsigned int        island_pop_size, 
                                                        CentroidsGenerator* cent_gen);
  static void                   resizeWithBestIndMutation (vector<GAPopulation*>& pops, unsigned int island_pop_size);
};



#endif /*EAISLANDSMODELCLUSTEREDPOPS_H*/
