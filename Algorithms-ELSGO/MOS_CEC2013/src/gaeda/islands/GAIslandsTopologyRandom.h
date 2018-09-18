#include "GAIslandsTopologyDynMedoidsBased.h"
#include "neighborconds.h"
#include <list>
#include <math.h>
#include <vector>

using namespace std; 

#ifndef GAISLANDSTOPOLOGYRANDOM__H_
#define GAISLANDSTOPOLOGYRANDOM__H_

/**
 * @brief Disconnected islands topology
 *
 * No edges on this topology
 */
class GAIslandsTopologyRandom : public EAIslandsTopology {
  unsigned int conn_degree_;
  
public:
   GAIslandsTopologyRandom( int          count, 
                            int          my_island_pos,
                            unsigned int conn_degree);

	~GAIslandsTopologyRandom ();
	
	void   generateNewNeighbors (const GAPopulation& pop); // Only to be implemented by dynamic topologies
  string getTopologyName      ();

};

#endif
