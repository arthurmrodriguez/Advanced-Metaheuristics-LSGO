#include "EAIslandsTopology.h"
#include "neighborconds.h"
#include <list>
#include <math.h>
#include <vector>

using namespace std; 

#ifndef GAISLANDSTOPOLOGYHYPERCUBE__H_
#define GAISLANDSTOPOLOGYHYPERCUBE__H_

class GAIslandsTopologyHyperCube : public EAIslandsTopology {
  unsigned int conn_degree_;

  void   calculateNeighbors   ();

public:
   GAIslandsTopologyHyperCube( int          count, 
                               int          my_island_pos, 
                               unsigned int conn_degree);

	~GAIslandsTopologyHyperCube ();	

  string getTopologyName      ();
};

#endif
