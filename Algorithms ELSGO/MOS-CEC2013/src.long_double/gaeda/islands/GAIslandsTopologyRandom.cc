#include "GAIslandsTopologyRandom.h"
#include "../logger/GALogger.h"
#include "../garandom.h"
#include <list>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <time.h>


using namespace std;


GAIslandsTopologyRandom::GAIslandsTopologyRandom( int           count, 
                                                  int           my_island_pos, 
                                                  unsigned int  conn_degree ) : 
                                                      EAIslandsTopology (count, my_island_pos ),
                                                      conn_degree_ (conn_degree) {}

GAIslandsTopologyRandom::~GAIslandsTopologyRandom() { }


string GAIslandsTopologyRandom::getTopologyName ( ){return "Random";}
	
void  GAIslandsTopologyRandom::generateNewNeighbors (const GAPopulation& pop){
  resetAdjMatrix();
  // Needs to be refactored a bit
  vector<int> available_neighbors(island_count_);
  for (int i=0;i<island_count_; i++) { 
    if ( i != (int) my_island_pos_ ) {  available_neighbors.push_back(i);}
  }
  
  for (unsigned int i=0; i<conn_degree_; i++){
    assert(available_neighbors.size() > 0);
    int neighbor_pos = GARandomInt(0,available_neighbors.size()-1);
    addNeighbor(available_neighbors[neighbor_pos]);
    available_neighbors.erase(available_neighbors.begin()+neighbor_pos);
  }
  
  // CHECK that everything is correct
  bool         are_no_neighbors_repeated = false;
  list<int>*   dests = getDestinations();
  vector<bool> neighbors_assigned( getIslandsCount() );
  
  for (unsigned int i=0; i<neighbors_assigned.size(); i++) neighbors_assigned[i] = false;
  
  for (list<int>::iterator it=dests->begin(); it!=dests->end(); it++){
    if (neighbors_assigned[*it]){
      are_no_neighbors_repeated = true;
      break;
    }
    else {
      neighbors_assigned[*it] = true;
    }
  }
  assert(!are_no_neighbors_repeated);
}
