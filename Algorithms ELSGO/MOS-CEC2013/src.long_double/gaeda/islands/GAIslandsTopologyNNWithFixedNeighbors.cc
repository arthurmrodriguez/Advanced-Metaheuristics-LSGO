#include "GAIslandsTopologyNNWithFixedNeighbors.h"
#include "neighborconds.h"
#include <list>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;


GAIslandsTopologyNNWithFixedNeighbors::GAIslandsTopologyNNWithFixedNeighbors( int                count, 
                                                                              neighbor_cond_func f, 
                                                                              unsigned int       min_nneighbors, 
                                                                              unsigned int       max_nneighbors,
                                                                              int       my_island_pos, 
                                                                              list<int>&         pfixed_neighbors_pos,
                                                                              CommManager&       comm_manager) : 
                                                                                GAIslandsTopologyNearestNeighbor(count,f,min_nneighbors,max_nneighbors,my_island_pos,comm_manager), 
                                                                                fixed_neighbors_pos(pfixed_neighbors_pos) { 
  stringstream message; message << "y en la topologia nn los venicos fijos con ";
  for (list<int>::iterator it =fixed_neighbors_pos.begin(); it!=fixed_neighbors_pos.end();it++){
    message << *it << ", ";
  }
  message << endl;
  cout << message.str();
}

GAIslandsTopologyNNWithFixedNeighbors::~GAIslandsTopologyNNWithFixedNeighbors() { }

void GAIslandsTopologyNNWithFixedNeighbors::resetAdjMatrix(){
  GAIslandsTopologyDynMedoidsBased::resetAdjMatrix();
  for (list<int>::iterator it=fixed_neighbors_pos.begin(); it != fixed_neighbors_pos.end(); it++) addNeighbor(*it);
}

string GAIslandsTopologyNNWithFixedNeighbors::getTopologyName (){
  stringstream message; message << "Nearest Neighbor with the following positions fixed ";
  for (list<int>::iterator it=fixed_neighbors_pos.begin(); it != fixed_neighbors_pos.end(); it++) message << *it << ", ";
  return message.str();
}

 
void GAIslandsTopologyNNWithFixedNeighbors::removeNeighborsUntilLimit(vector< vector<long double> >& distance_matrix){
  int neighbors_to_remove = getNumNeighbors() - max_num_neighbors_;      
  int max_possible_neighbors_to_remove = getNumNeighbors() - fixed_neighbors_pos.size();
  assert (neighbors_to_remove <= max_possible_neighbors_to_remove);
  
  GAIslandsTopologyNearestNeighbor::removeNeighborsUntilLimit(distance_matrix);
}


unsigned int GAIslandsTopologyNNWithFixedNeighbors::findFurthestNeighborMedoidPos(vector< vector<long double> >& distance_matrix ){
  auto_ptr< list<int> > neighbors_pos( getDestinations() );
  for (list<int>::iterator it=fixed_neighbors_pos.begin(); it!=fixed_neighbors_pos.end(); it++){
    neighbors_pos->remove(*it); // We remove all the fixed positions neighbors
  }
  
  unsigned int          selected_pos = static_cast<unsigned int> ( *( neighbors_pos->begin() ) );       //initialized with the first pos

  for( list<int>::iterator iter=neighbors_pos->begin(); iter!= neighbors_pos->end(); iter++) {
    if ( distance_matrix[my_island_pos_][*iter] > distance_matrix[my_island_pos_][selected_pos] ){
      selected_pos = static_cast<unsigned int> (*iter);
    }
  }
  
  return selected_pos;
} 
