/*****************************************************************************
 * GAEDAlib: A C++ GA library with EDA and multiprocessor (MPI) support      *
 *                                                                           *
 * (C) 2005 Pedro Diaz (pdiaz@laurel.datsi.fi.upm.es)                        *
 *                                                                           *
 * GAEDAlib is distributed under the terms of the BSD software license       *
 *                                                                           *
 * GAEDAlib is heavily based on GAlib, a C++ GA library by Mathew Wall:      *
 * Copyright (c) 1995-1996 Massachusetts Institute of Technology (MIT)       *
 * Copyright (c) 1996-2000 Matthew Wall (author of GAlib)                    *
 *                                                                           *
 * Some portions of GAEDAlib's source code come from the GNU C++ compiler    *
 * library and therefore are covered under the terms of a different license, *
 * the GNU Public License.                                                   *
 *                                                                           *
 * You should have received a file named LICENSE along with this software.   *
 * This file contains more information about the licensing conditions of     *
 * GAEDAlib as well as the full text of each license involved.               *
 *                                                                           *
 * The file AUTHORS lists the people who have contributed (directly or       *
 * indirectly) to GAEDAlib                                                   *
 *****************************************************************************/


/**
 * @file
 * @brief GAIslandsTopologyPairs class impl.
 *
 * GAIslandsTopologyPairs class implementation. 
 */


/* INCLUDES */
#include <list>
#include "GAIslandsTopologyPairs.h"

using namespace std;

/**
 * Class constructor
 *
 * @param count Number of islands in the archipelago
 */
GAIslandsTopologyPairs::GAIslandsTopologyPairs( int count ) : EAIslandsTopology(count){ 
	/* We are not going to use the adjacency matrix, so we don't
	 * do anything here
	 */
}


/**
 * Class destructor
 */
GAIslandsTopologyPairs::~GAIslandsTopologyPairs() {
	/* Nothing to do here! */
}



list<int>* GAIslandsTopologyPairs::getDestinations( int islandID ) {

	/* No need to look at the adjacency matrix */
	list<int> *res;
	res = new list<int>;
	bool isOdd = islandID %2;

	if (islandID == getIslandsCount() -1 && !isOdd) {
		// Last island on an archipielago with odd number
		// of islands. This one is disconnected
	} else {
		if (isOdd) {
			res->push_back( islandID-1 );
		} else {
			res->push_back( islandID+1 );
		}
	}
	return res;
}


list<int>* GAIslandsTopologyPairs::getOrigins( int islandID ) {
	return getDestinations( islandID );
}

