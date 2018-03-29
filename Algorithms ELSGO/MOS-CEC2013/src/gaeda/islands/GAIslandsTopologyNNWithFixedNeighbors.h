#include "GAIslandsTopologyNearestNeighbor.h"
#include "neighborconds.h"
#include <list>
#include <math.h>
#include <vector>

using namespace std; 

#ifndef GAISLANDSTOPOLOGYNNWITHFIXEDNEIGHBORS__H_
#define GAISLANDSTOPOLOGYNNWITHFIXEDNEIGHBORS__H_

/**
 * @brief Disconnected islands topology
 *
 * No edges on this topology
 */
class GAIslandsTopologyNNWithFixedNeighbors : public GAIslandsTopologyNearestNeighbor {
  list<int> fixed_neighbors_pos;

  void         removeNeighborsUntilLimit     (vector< vector<double> >& distance_matrix);
  unsigned int findFurthestNeighborMedoidPos (vector< vector<double> >& distance_matrix );
  string       getTopologyName               (                           );
  void         resetAdjMatrix                (                          );
public:
   GAIslandsTopologyNNWithFixedNeighbors( int                count, 
                                          neighbor_cond_func f, 
                                          unsigned int       min_nneighbors, 
                                          unsigned int       max_nneighbors,
                                          int                my_island_pos, 
                                          list<int>&         fixed_neighbors_pos,
                                          CommManager&       comm_manager);

	~GAIslandsTopologyNNWithFixedNeighbors ();

};

#endif
