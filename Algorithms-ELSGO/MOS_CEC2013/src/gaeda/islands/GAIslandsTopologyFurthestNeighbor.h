#include "GAIslandsTopologyDynMedoidsBased.h"
#include "neighborconds.h"
#include <list>
#include <math.h>
#include <vector>

using namespace std; 

#ifndef GAISLANDSTOPOLOGYFURTHESTNEIGHBOR__H_
#define GAISLANDSTOPOLOGYFURTHESTNEIGHBOR__H_

class GAIslandsTopologyFurthestNeighbor : public GAIslandsTopologyDynMedoidsBased {
  int    min_num_neighbors_;
  int    max_num_neighbors_;

  string getTopologyName();
public:
   GAIslandsTopologyFurthestNeighbor( int                count, 
                                      unsigned int       min_nneighbors, 
                                      unsigned int       max_nneighbors,
                                      int                my_island_pos,
                                      CommManager&       comm_manager);

	~GAIslandsTopologyFurthestNeighbor ();	

  void         regenerateNeihborhoods (vector<GAGenome*>& medoids );
  unsigned int findFurthestNotNeighborMedoidPos(vector<GAGenome*>& medoids  );
  void         addNeighborsUntilLimit(vector<GAGenome*>& medoids);
};

#endif
