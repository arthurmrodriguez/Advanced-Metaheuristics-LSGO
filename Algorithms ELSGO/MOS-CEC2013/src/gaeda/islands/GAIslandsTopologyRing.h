#ifndef GAISLANDSTOPOLOGYRING_H
#define GAISLANDSTOPOLOGYRING_H

#include <list>

#include "EAIslandsTopology.h"


using namespace std; 


/**
 * @brief Circular islands topology
 *
 * This topology arranges the islands on a ring, where an island "n" is connected
 * to the island "n+1 (mod) count"
 */
class GAIslandsTopologyRing : public EAIslandsTopology {

	public:
		/* Constructor and destructor */
		GAIslandsTopologyRing( int count, int my_island_pos );
		~GAIslandsTopologyRing();

		list<int>* getDestinations( int islandID );
    string getTopologyName ();
};

#endif
