#ifndef GAISLANDSTOPOLOGYRING2_H
#define GAISLANDSTOPOLOGYRING2_H

#include <list>

#include "EAIslandsTopology.h"


using namespace std; 


/**
 * @brief Doubly circular islands topology
 *
 * This topology arranges the islands on a long double ring, where an island "n" is connected
 * to the islands "n+1 (mod) count" and "n-1 (mod) count"
 */
class GAIslandsTopologyRing2 : public EAIslandsTopology {

	public:
		/* Constructor and destructor */
		GAIslandsTopologyRing2( int count, int my_island_pos );
		~GAIslandsTopologyRing2();


		/* Members */
		list<int>* getDestinations( int islandID );
		list<int>* getOrigins( int islandID );
		string     getTopologyName ();
	
};

#endif
