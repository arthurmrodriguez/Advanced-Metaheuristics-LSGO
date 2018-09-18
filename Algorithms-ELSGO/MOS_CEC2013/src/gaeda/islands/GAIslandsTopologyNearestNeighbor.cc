#include "GAIslandsTopologyNearestNeighbor.h"
#include "neighborconds.h"
#include "../logger/GALogger.h"
#include <list>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <memory>

using namespace std;


GAIslandsTopologyNearestNeighbor::GAIslandsTopologyNearestNeighbor( int                count, 
                                                                    neighbor_cond_func f, 
                                                                    unsigned int       min_nneighbors, 
                                                                    unsigned int       max_nneighbors,
                                                                    int                my_island_pos,
                                                                    CommManager&       comm_manager) : 
                                                               GAIslandsTopologyDynMedoidsBased(count,my_island_pos, comm_manager ), 
                                                               neighbor_cond_                  (f              ),
                                                               min_num_neighbors_              (min_nneighbors ),
                                                               max_num_neighbors_              (max_nneighbors ) { 
  assert((int) min_nneighbors <= count);
  assert((int) max_nneighbors > 0);
}

GAIslandsTopologyNearestNeighbor::~GAIslandsTopologyNearestNeighbor() { }


string GAIslandsTopologyNearestNeighbor::getTopologyName (){
  return "Nearest Neighbor";
}


/* 
 * It returns the positions of the neighbors that have been to be added in order to pass the minimum neighbors limit
 */
void GAIslandsTopologyNearestNeighbor::regenerateNeihborhoods(vector<GAGenome*>& medoids){
  calculateInitialNeighbors(medoids);

  {// LOGGING INITIAL NEIGHBORS 
    stringstream message; message << "and the neighbors of island " <<  my_island_pos_ << " are ";
    auto_ptr< list<int> > destinations( getDestinations() );
    for( list<int>::iterator iter=destinations->begin(); iter!= destinations->end(); iter++) {message << *iter << ", ";} message << endl;
    GALogger::instance()->appendLogMessage("Nearest Neighbor Topology initialization initial neighbors are:",message.str(), GALogger::normal);
  }

  int num_neighbors = static_cast<int> ( getNumNeighbors() );
  if      ( num_neighbors < min_num_neighbors_) addNeighborsUntilLimit(medoids);
  else if ( num_neighbors > max_num_neighbors_) removeNeighborsUntilLimit(); 

}

void GAIslandsTopologyNearestNeighbor::calculateInitialNeighbors(vector<GAGenome*>& medoids){
  // We calculate the initial neighbors, we take the medoids and we see if the neighbor condition is satisfied for
  // the rest of the medoids
  for (unsigned i=0; i<num_islands_; i++){
    
    //cout << "La distancia de mi pos " << my_island_pos_ <<  " " <<  *(medoids[my_island_pos_]) << " a la isla " << i << " " << *(medoids[i]) << " es " << dist_matrix_[my_island_pos_*num_islands_+i] << endl;
    //printf("La distancia de mi pos(%d) con valor %s a la isla %d %s es %20.20lf\n",my_island_pos_,*(medoids[my_island_pos_]),i,*(medoids[i]),dist_matrix_[my_island_pos_*num_islands_+i]);
    if ( i==my_island_pos_ || dist_matrix_[my_island_pos_*num_islands_+i] == 0 ) continue;
    
    bool are_my_island_pos_i_neighbors = true;
    for (unsigned k=0; k<num_islands_ && are_my_island_pos_i_neighbors; k++){
      // If the individual is the same or is at the same position, we continue the loop since it does not count for the
      // verification of neighborhood
      if (k == my_island_pos_                            || k == i || 
          dist_matrix_[my_island_pos_*num_islands_+k] == 0 || dist_matrix_[i*num_islands_+k] == 0 ) continue;
      if ( !neighbor_cond_(my_island_pos_,i,k, dist_matrix_, num_islands_) ){
        {//LOGGING the island who breaks the condition
          stringstream message; message << "the island " << k << " breaks the condition of neighborhood between islands " << my_island_pos_ << " and " << i << endl;
          GALogger::instance()->appendLogMessage("Nearest Neighbor Topology breaking condition:",message.str(), GALogger::normal);
        }
          are_my_island_pos_i_neighbors = false;
      }      
    }

    if (are_my_island_pos_i_neighbors) addNeighbor(i);
  }
}
 
void GAIslandsTopologyNearestNeighbor::addNeighborsUntilLimit(vector<GAGenome*>& medoids){
  int  neighbors_to_add = min_num_neighbors_ - getNumNeighbors();      

  stringstream message; message << "We need to add the following neighbors:" << endl; // For logging

  while ( neighbors_to_add-- > 0 ){
   unsigned int new_neighbor_pos = findNearestNotNeighborMedoidPos(medoids);
   addNeighbor(new_neighbor_pos);
   message << new_neighbor_pos << ", ";
  } 
  message << endl; GALogger::instance()->appendLogMessage("Nearest Neighbor Topology, adding neighbors", message.str(),GALogger::normal);
}

void GAIslandsTopologyNearestNeighbor::removeNeighborsUntilLimit(){
  int  neighbors_to_remove = getNumNeighbors() - max_num_neighbors_;      

  stringstream message; message << "We need to remove the following neighbors:" << endl; // For logging

  while ( neighbors_to_remove-- > 0 ){
    unsigned int neighbor_pos = findFurthestNeighborMedoidPos();
    removeNeighbor(neighbor_pos);
    message << neighbor_pos << ", ";
  }
  message << endl; GALogger::instance()->appendLogMessage("Nearest Neighbor Topology, removing neighbors", message.str(),GALogger::normal);
}


// This should be refactored to a more efficient version where some sorting is involved so that this function return not a single individual
// but a group of the nearest not neighbors.
unsigned int GAIslandsTopologyNearestNeighbor::findNearestNotNeighborMedoidPos(vector<GAGenome*>& medoids  ){
  unsigned first_pos;
  for (first_pos=0; first_pos<num_islands_; first_pos++){        // We obtain the first position of a not neighbor medoid
    if ( first_pos != my_island_pos_ && !areNeighbors(my_island_pos_,first_pos) ) break;
  }

  unsigned int selected_pos = first_pos;                    
   
  for (unsigned i=first_pos; i<num_islands_; i++){
    // If the selected position is the same as the medoid to which we are adding neighbors
    // or if the selected position is already a selected neighbor, we continue the loop 
    if ( i==my_island_pos_ || areNeighbors(my_island_pos_,i) ) continue;
    
    if ( dist_matrix_[my_island_pos_*num_islands_+i] != 0 && dist_matrix_[my_island_pos_*num_islands_+i] < 
         dist_matrix_[my_island_pos_*num_islands_+selected_pos] ){
      selected_pos = i;
    } 
  }

  return selected_pos;
} 

unsigned int GAIslandsTopologyNearestNeighbor::findFurthestNeighborMedoidPos(){
  auto_ptr< list<int> > neighbors_pos( getDestinations() );
  unsigned              selected_pos = static_cast<unsigned int> ( *( neighbors_pos->begin() ) );       //initialized with the first pos

  for( list<int>::iterator iter=neighbors_pos->begin(); iter!= neighbors_pos->end(); iter++) {
    if ( dist_matrix_[my_island_pos_*num_islands_+*iter] > dist_matrix_[my_island_pos_*num_islands_+selected_pos] ){
      selected_pos = static_cast<unsigned int> (*iter);
    }
  }
  
  return selected_pos;
} 

