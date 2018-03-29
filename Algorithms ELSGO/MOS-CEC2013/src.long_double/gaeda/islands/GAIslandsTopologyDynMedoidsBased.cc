#include "GAIslandsTopologyDynMedoidsBased.h"
#include "neighborconds.h"
#include "../logger/GALogger.h"
#include "islandsutils.h"
#include <list>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <memory>

using namespace std;

GAIslandsTopologyDynMedoidsBased::GAIslandsTopologyDynMedoidsBased( int                count, 
                                                                    int                my_island_pos,
                                                                    CommManager&       comm_manager ) : 
                                                                                                  EAIslandsTopology (count, my_island_pos), 
                                                                                                  comm_manager_(comm_manager),
                                                                                                  num_islands_(comm_manager.getNumIslands()){
  dist_matrix_ = new long double[num_islands_*num_islands_];                                                                      
}

GAIslandsTopologyDynMedoidsBased::~GAIslandsTopologyDynMedoidsBased() { 
  delete[] dist_matrix_;
}

void GAIslandsTopologyDynMedoidsBased::regenerateNeihborhoods (vector<GAGenome*>& medoids ){}

void GAIslandsTopologyDynMedoidsBased::calculateDistanceMatrix(vector<GAGenome*>& medoids) {
  assert(medoids.size() == num_islands_); 
  
  // Compute the distance matrix 
  for (unsigned i=0; i<num_islands_; i++){
    dist_matrix_[i*num_islands_+i] = 0;
    for (unsigned j=i+1; j<num_islands_; j++){
      GAGenome& medoid_i = *(medoids[i]);
      GAGenome& medoid_j = *(medoids[j]);
      dist_matrix_[i*num_islands_+j] = dist_matrix_[j*num_islands_+i] = medoid_i.compare(medoid_j);  
    }
  }
}

void GAIslandsTopologyDynMedoidsBased::generateNewNeighbors (const GAPopulation& our_pop) {
  GAGenome&           my_medoid = our_pop.medoid(); 
  vector<GAGenome*>*  medoids;

  if ( comm_manager_.isIslandSlave() ) {                  
    comm_manager_.sendIndToMaster(my_medoid);                                        // First, it sends the medoid to the master
    medoids = convertPopToVector( comm_manager_.receivePopFromMaster() );
    assert( (int) medoids->size() == comm_manager_.getNumIslands() );
  }
  else {                                                                              // The master node
    medoids = comm_manager_.receiveOneIndFromAllSlaveIslands();                       // First, it collects all the medoids
    medoids->insert( medoids->begin(),my_medoid.clone() );                            // We clone the medoid since at the end, the list is being deleted                     
    comm_manager_.sendIndsToSlaveIslands(*medoids);                                   // Then it resends them to all the islands
  }

  GALogger::instance()->appendMedoidsStats("GAIslandsTopologyDynMedoidsBased::run, ", our_pop, *medoids, my_medoid);
  
  setMedoids(medoids);
  delete medoids;
}

void GAIslandsTopologyDynMedoidsBased::setMedoids (vector<GAGenome*>* medoids ){
  resetAdjMatrix();                                                                   // Create the adjacency matrix
  calculateDistanceMatrix(*medoids);
  regenerateNeihborhoods(*medoids);
  
  {// LOGGING INITIAL NEIGHBORS
    stringstream message; message << "and the neighbors of island " <<  my_island_pos_ << " are ";
    auto_ptr< list<int> > destinations( getDestinations() );
    for( list<int>::iterator iter=destinations->begin(); iter!= destinations->end(); iter++) {message << *iter << ", ";} message << endl;
    stringstream title; title << getTopologyName() << " Topology initialization: neighbors are:";
    GALogger::instance()->appendLogMessage(title.str() ,message.str(), GALogger::normal);
  }
}


