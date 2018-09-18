#ifndef GAISLANDSTOPOLOGYA2A_H
#define GAISLANDSTOPOLOGYA2A_H

#include <list>

#include "EAIslandsTopology.h"


using namespace std; 


/**
 * @brief All-to-all islands topology
 *
 * This topology arranges the islands on a complete graph , where an island "n" is connected
 * to every each other island ( [1..n-1] U [n+1..Count] )
 */
class GAIslandsTopologyA2A : public EAIslandsTopology {

	public:
		/* Constructor and destructor */
		GAIslandsTopologyA2A( int count, int my_island_pos );
		~GAIslandsTopologyA2A();
		string getTopologyName ();
};

#endif
