#include "GAIslandsTopologyDynMedoidsBased.h"
#include "GAIslandsTopologyHyperCube.h"
#include "GAIslandsTopologyNNWithFixedNeighbors.h"
#include "neighborconds.h"
#include <list>
#include <math.h>
#include <vector>

using namespace std; 

#ifndef GAISLANDSTOPOLOGYHYBRID__H_
#define GAISLANDSTOPOLOGYHYBRID__H_

class GAIslandsTopologyHybrid : public GAIslandsTopologyDynMedoidsBased {
  GAIslandsTopologyNNWithFixedNeighbors* nearestn_top;
  GAIslandsTopologyHyperCube*            hypercube_top;
    
  void   calculateNeighbors (vector<GAGenome*>& medoids );
  string getTopologyName    (                           );

public:

   GAIslandsTopologyHybrid( int count, 
                            neighbor_cond_func f, 
                            int                my_island_pos, 
                            unsigned int       conn_degree,
                            float              dyn_participation,
                            CommManager&       comm_manager);
      
	~GAIslandsTopologyHybrid ();	

};

#endif
