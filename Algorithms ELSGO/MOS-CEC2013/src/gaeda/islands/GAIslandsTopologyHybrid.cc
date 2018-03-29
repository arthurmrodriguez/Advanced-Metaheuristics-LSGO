#include "GAIslandsTopologyHybrid.h"
#include "neighborconds.h"
#include <list>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;


GAIslandsTopologyHybrid::GAIslandsTopologyHybrid( int count, 
                                                      neighbor_cond_func f, 
                                                      int                my_island_pos, 
                                                      unsigned int       conn_degree,
                                                      float              static_particip,
                                                      CommManager&       comm_manager) : 
                                                         GAIslandsTopologyDynMedoidsBased( count,my_island_pos, comm_manager ) {
  hypercube_top = new GAIslandsTopologyHyperCube ( count, my_island_pos,conn_degree);
   
  // We get all the hypercube neighbors and we only get the first ones until we exceed the participation ratio
  auto_ptr< list<int> > hypercube_neighbors ( hypercube_top->getDestinations() );
  unsigned int num_hypercube_neighbors_used = (unsigned int) round ( conn_degree * static_particip );
  list<int> hypercube_neighbors_passed;

  for (list<int>::iterator it=hypercube_neighbors->begin(); 
       it!=hypercube_neighbors->end(), num_hypercube_neighbors_used > 0; 
       it++, num_hypercube_neighbors_used--){ 
    hypercube_neighbors_passed.push_back(*it);
  }
  
  nearestn_top  = new GAIslandsTopologyNNWithFixedNeighbors ( count, 
                                                          f, 
                                                          (unsigned int) conn_degree, 
                                                          (unsigned int) conn_degree, 
                                                          my_island_pos,  
                                                          hypercube_neighbors_passed,
                                                          comm_manager); 

}

GAIslandsTopologyHybrid::~GAIslandsTopologyHybrid() { 
  delete nearestn_top;
  delete hypercube_top;
}

void GAIslandsTopologyHybrid::calculateNeighbors(vector<GAGenome*>& medoids){
  nearestn_top->setMedoids( medoids );

  auto_ptr< list<int> > nearestn_neighbors ( nearestn_top->getDestinations() );

  for (list<int>::iterator it=nearestn_neighbors->begin(); it!=nearestn_neighbors->end(); it++) addNeighbor(*it);
}

string GAIslandsTopologyHybrid::getTopologyName (){
  return "Hybrid";
}

