#include "EAIslandsModelClusteredPops.h"
#include <assert.h>
#include <stdlib.h>
#include <memory>

// Auxiliary functions declarations
static unsigned int    getNearestMedoidPos            (vector<GAGenome*>& centroids, GAGenome& ind);
static void            repopulateEmptyPops            (vector<GAPopulation*>& pops);
vector<GAPopulation*>* createInitialkmeansPopulations (GAPopulation& pop, unsigned int num_islands, CentroidsGenerator* cent_gen);

EAIslandsModelClusteredPops::EAIslandsModelClusteredPops( GAGeneticAlgorithm*     ga,        
                                                          EAIslandsTopology*      top,   
                                                          unsigned int            mig_period,
                                                          EAEmigrantsSelector*    emigrants_selector,
                                                          EAImmigrantsSelector*   immigrants_selector,
                                                          unsigned int            clustering_period,
                                                          cluster_func_t          cluster_func,
                                                          popssizes_resize_func_t popsizes_resize_func,
                                                          CentroidsGenerator*     centroids_gen,         // It is going to be deleted in the destructor                                                          
                                                          CommManager&            comm_manager):
                                                               EAIslandsModelSync    (ga,top,mig_period,emigrants_selector,immigrants_selector,comm_manager),
                                                               clustering_period_    (clustering_period),
                                                               cluster_func_         (cluster_func),
                                                               popsizes_resize_func_ (popsizes_resize_func),
                                                               centroids_gen_        (centroids_gen) {}


EAIslandsModelClusteredPops::~EAIslandsModelClusteredPops(){
  delete centroids_gen_;
}

void EAIslandsModelClusteredPops::postStep(){
  if ( ea_->statistics().generation() % clustering_period_ == 0){
    GAPopulation& our_pop = const_cast<GAPopulation &> ( ea_->population() );
    
    if ( comm_manager_.isIslandSlave() ) {
      comm_manager_.sendPopToMaster(our_pop);
      auto_ptr<GAPopulation> pop_received ( comm_manager_.receivePopFromMaster() );
      ea_->population( *pop_received );
    }
    else {    // MASTER
      auto_ptr<GAPopulation> pop_received ( comm_manager_.receivePopsFromAllSlaveIslands() );
      assert (pop_received->size() > 0);
      for (unsigned i=0; i<our_pop.size(); i++) pop_received->add( our_pop.individual(i) );                   // We add the master population
      assert ( pop_received->size() > our_pop.size()  );
      
      auto_ptr< vector<GAPopulation*> > new_island_pops ( obtainNewIslandsPops ( *pop_received ) );
      for (unsigned int i=0; i<new_island_pops->size(); i++) assert( (*new_island_pops)[i]->size() > 0 );
      
      comm_manager_.sendPopsToSlaveIslands( *new_island_pops );
      
      ea_->population( *( (*new_island_pops)[comm_manager_.MASTER_RANK] ) );                             // The master pop is set
      
      for (unsigned int i=0; i<new_island_pops->size(); i++) delete (*new_island_pops)[i];               // the new pops need to be deleted
    }
  }
}
  
vector<GAPopulation*>* EAIslandsModelClusteredPops::obtainNewIslandsPops( GAPopulation& pop ){
  unsigned int num_islands     = comm_manager_.getNumIslands(); 
  unsigned int island_pop_size = pop.size() / num_islands;
  
  assert (pop.size() % num_islands == 0); // the size must be a multiple of the number of islands
  
  vector<GAPopulation*>* new_pops = cluster_func_(pop, num_islands, island_pop_size, centroids_gen_);
  assert(new_pops != NULL && new_pops->size() == num_islands);
  for (unsigned int i=0; i<new_pops->size(); i++) assert( ((*new_pops)[i]) != NULL );
  
  popsizes_resize_func_(*new_pops, island_pop_size);
  for(unsigned int i=0; i<new_pops->size(); i++) (*new_pops)[i]->sort();
  return new_pops;
}

vector<GAPopulation*>* EAIslandsModelClusteredPops::kmeansCluster( GAPopulation&       pop, 
                                                                   unsigned int        num_islands, 
                                                                   unsigned int        island_pop_size, 
                                                                   CentroidsGenerator* cent_gen ){
  
  vector<GAPopulation*>* new_pops = createInitialkmeansPopulations(pop, num_islands, cent_gen);
  
  bool inds_pos_changed;
  do {
    for (unsigned int i=0; i<new_pops->size(); i++) {
      assert( (*new_pops)[i]->size() <= pop.size() ); assert( (*new_pops)[i]->size() > 0 );
      for ( unsigned j=0; j <  (*new_pops)[i]->size(); j++){
        GAGenome& ind = (*new_pops)[i]->individual(j); assert ( (*new_pops)[i]->individual(j) != NULL ); assert (&ind != NULL);
      }
    }
           
    vector<GAGenome*> medoids( new_pops->size() );
    for (unsigned int i=0; i<new_pops->size(); i++) medoids[i] = &((*new_pops)[i]->medoid());      // We compute the medoids from the
                                                                                                   // populations and we reassign the inds
    inds_pos_changed = false;                                                                      // to the nearest medoid
    
    for (unsigned int pop_pos=0; pop_pos<new_pops->size(); pop_pos++){                             
      GAPopulation* island_pop = (*new_pops)[pop_pos];
      for (unsigned ind_pop_pos=0; ind_pop_pos<island_pop->size(); ind_pop_pos++){
        GAGenome& ind = island_pop->individual(ind_pop_pos); assert (&ind != NULL);
        
        unsigned int nearest_cent_pos  = getNearestMedoidPos(medoids,ind);
        
        if (nearest_cent_pos != pop_pos){
          GAGenome* old_genome = island_pop->remove(ind_pop_pos);
          (*new_pops)[nearest_cent_pos]->add(old_genome);
          inds_pos_changed = true;
        }       
      }
    }
    repopulateEmptyPops(*new_pops);    
  } while (inds_pos_changed); 
  
  return new_pops;
}

vector<GAPopulation*>* createInitialkmeansPopulations(GAPopulation& pop, unsigned int num_islands, CentroidsGenerator* cent_gen){
  vector<GAPopulation*>* new_pops = new vector<GAPopulation*>(num_islands);                        // creation of the new populations
  for (unsigned int i=0; i<new_pops->size(); i++) (*new_pops)[i] = new GAPopulation();
  
  auto_ptr< vector<GAGenome*> > centroids ( cent_gen->generateCentroids() );  
  for (unsigned i=0; i<pop.size(); i++) {                                                               // Initially we assign all the population
    GAGenome&    ind               = pop.individual(i);                                            // to the new populations based on the  
    unsigned int nearest_cent_pos  = getNearestMedoidPos( *centroids,ind);                         // proximity to the generated centroids
    (*new_pops)[nearest_cent_pos]->add(ind);
  }
  repopulateEmptyPops(*new_pops);    
  
  for (unsigned int i=0; i<centroids->size(); i++) delete (*centroids)[i];
  
  return new_pops;  
}

void repopulateEmptyPops(vector<GAPopulation*>& pops){
  for(unsigned int i=0; i<pops.size(); i++){
    int pop_pos_source;
    if (pops[i]->size() == 0){                                                  // An empty population  
      do {
        pop_pos_source = rand() % pops.size();                                  // We select a random island that is not empty
      } while (pop_pos_source == (int) i || pops[pop_pos_source]->size() < 2);  // and we take out the worst individual
      GAGenome* old_genome = pops[pop_pos_source]->remove();
      pops[i]->add(old_genome);      
    }
  }
}

unsigned int getNearestMedoidPos(vector<GAGenome*>& centroids, GAGenome& ind){
  unsigned int selected_cent_pos = 0;    
  double       selected_min_dist = ind.compare( *(centroids[selected_cent_pos]) );
  double       tmp_min_dist      = 0.0;
  
  for (unsigned int cent_pos=1; cent_pos<centroids.size(); cent_pos++){             // We have previously compared with the first case
    tmp_min_dist = ind.compare( *(centroids[cent_pos]) );
    
    if (tmp_min_dist < selected_min_dist){
      selected_min_dist = tmp_min_dist;
      selected_cent_pos = cent_pos;
    }
  }
  return selected_cent_pos;
}


void  EAIslandsModelClusteredPops::resizeWithBestIndMutation(vector<GAPopulation*>& all_pops, unsigned int island_pop_size){
  const float MUTATION_PROB = 0.5;
  
  for(unsigned int pop_pos=0; pop_pos<all_pops.size(); pop_pos++){
    GAPopulation* pop = all_pops[pop_pos]; 
    assert(pop != NULL);
    if ( pop->size() > island_pop_size ) {
      pop->sort();                                    // need to be sorted in order to
      pop->size(island_pop_size);                     // delete the worst ones
    }
    else if ( pop->size() < island_pop_size ) {
      int ninds_to_create = abs( (int) ( pop->size() - island_pop_size ) );
      assert ( ninds_to_create <= (int) all_pops.size() * (int) island_pop_size ); assert( pop->size() > 0);
        
      for (int i=0; i<ninds_to_create; i++){
        pop->sort();  
        GAGenome* new_ind = pop->best().clone();                                    // we add a new individual based on the mutation of the best one from the pop
        new_ind->mutate(MUTATION_PROB);                                             // the best one from the pop
        pop->add( new_ind ); 
      }      
    }
  }  
}

