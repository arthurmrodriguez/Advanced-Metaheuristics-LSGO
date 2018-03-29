#include "GAIslandsTopologyDynMedoidsBased.h"
#include "neighborconds.h"
#include <list>
#include <math.h>
#include <vector>

using namespace std; 

#ifndef GAISLANDSTOPOLOGYNEARESTNEIGHBOR__H_
#define GAISLANDSTOPOLOGYNEARESTNEIGHBOR__H_

/**
 * @brief Disconnected islands topology
 *
 * No edges on this topology
 */
class GAIslandsTopologyNearestNeighbor : public GAIslandsTopologyDynMedoidsBased {
protected:
  neighbor_cond_func neighbor_cond_;
  int                min_num_neighbors_;
  int                max_num_neighbors_;
  
  void                 calculateInitialNeighbors       (vector<GAGenome*>& medoids);
  void                 addNeighborsUntilLimit          (vector<GAGenome*>& medoids);
  virtual void         removeNeighborsUntilLimit       (                          );
  unsigned int         findNearestNotNeighborMedoidPos (vector<GAGenome*>& medoids);
  virtual unsigned int findFurthestNeighborMedoidPos   (                          );
  virtual string       getTopologyName                 (                          );
  void                 regenerateNeihborhoods          (vector<GAGenome*>& medoids); 

public:
   GAIslandsTopologyNearestNeighbor( int                count, 
                                     neighbor_cond_func f, 
                                     unsigned int       min_nneighbors, 
                                     unsigned int       max_nneighbors,
                                     int                my_island_pos,
                                     CommManager&       comm_manager);

	~GAIslandsTopologyNearestNeighbor ();	

};

#endif
