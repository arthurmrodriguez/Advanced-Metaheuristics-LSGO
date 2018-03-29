#include "GAIslandsTopologyFurthestNeighbor.h"
#include "../logger/GALogger.h"
#include <list>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <memory>

using namespace std;

GAIslandsTopologyFurthestNeighbor::GAIslandsTopologyFurthestNeighbor( int              count, 
                                                                    unsigned int       min_nneighbors, 
                                                                    unsigned int       max_nneighbors,
                                                                    int                my_island_pos,
                                                                    CommManager&       comm_manager) : 
                                                                      GAIslandsTopologyDynMedoidsBased (count, my_island_pos, comm_manager ),
                                                                      min_num_neighbors_               (min_nneighbors),
                                                                      max_num_neighbors_               (max_nneighbors) {
  assert((int) min_nneighbors <= count);
  assert((int) max_nneighbors > 0);
}

GAIslandsTopologyFurthestNeighbor::~GAIslandsTopologyFurthestNeighbor() { }

void GAIslandsTopologyFurthestNeighbor::regenerateNeihborhoods(vector<GAGenome*>& medoids){
  calculateDistanceMatrix(medoids);  //Calculate the distance matrix

  addNeighborsUntilLimit(medoids);

  {// LOGGING INITIAL NEIGHBORS 
    stringstream message; message << "and the neighbors of island " <<  my_island_pos_ << " are ";
    auto_ptr< list<int> > destinations( getDestinations() );
    for( list<int>::iterator iter=destinations->begin(); iter!= destinations->end(); iter++) {message << *iter << ", ";} message << endl;
    GALogger::instance()->appendLogMessage("Furthest Neighbor Topology initialization initial neighbors are:",message.str(), GALogger::normal);
  }
}

void GAIslandsTopologyFurthestNeighbor::addNeighborsUntilLimit(vector<GAGenome*>& medoids){
  int  neighbors_to_add = min_num_neighbors_ - getNumNeighbors();      

  while ( neighbors_to_add-- > 0 ){
   unsigned int new_neighbor_pos = findFurthestNotNeighborMedoidPos(medoids);
   addNeighbor(new_neighbor_pos);
  } 
}


string GAIslandsTopologyFurthestNeighbor::getTopologyName (){
  return "Furthest Neighbor";
}

// This should be refactored to a more efficient version where some sorting is involved so that this function return not a single individual
// but a group of the nearest not neighbors.
unsigned int GAIslandsTopologyFurthestNeighbor::findFurthestNotNeighborMedoidPos(vector<GAGenome*>& medoids  ){

  unsigned first_pos;
  for (first_pos=0; first_pos <num_islands_; first_pos++){        // We obtain the first position of a not neighbor medoid
    if ( first_pos != my_island_pos_ && !areNeighbors(my_island_pos_,first_pos) ) break;
  }

  unsigned selected_pos = first_pos;                    
   
  for (unsigned i=first_pos; i<num_islands_; i++){
    // If the selected position is the same as the medoid to which we are adding neighbors
    // or if the selected position is already a selected neighbor, we continue the loop 
    if ( i==my_island_pos_ || areNeighbors(my_island_pos_,i) ) continue;
    
    if ( dist_matrix_[my_island_pos_*num_islands_+i] != 0 && dist_matrix_[my_island_pos_*num_islands_+i] > 
         dist_matrix_[my_island_pos_*num_islands_+selected_pos] ){
      selected_pos = i;
    } 
  }

  return selected_pos;
} 
