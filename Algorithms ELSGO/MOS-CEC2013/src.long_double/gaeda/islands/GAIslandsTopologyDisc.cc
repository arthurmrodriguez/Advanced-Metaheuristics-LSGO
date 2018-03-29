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
 * @brief GAIslandsTopologyDisc class impl.
 *
 * GAIslandsTopologyDisc implementation. This class represents
 * a disconnected topology
 */


/* INCLUDES */
#include "GAIslandsTopologyDisc.h"
#include <list>
#include <assert.h>

using namespace std;

/**
 * Class constructor
 *
 * @param count Number of islands in the archipelago
 */
GAIslandsTopologyDisc::GAIslandsTopologyDisc( int count ) : EAIslandsTopology(count){}


/**
 * Class destructor
 */
GAIslandsTopologyDisc::~GAIslandsTopologyDisc() {
	/* Nothing to do here! */
}



list<int>* GAIslandsTopologyDisc::getDestinations( int islandID ) {
  assert (islandID < island_count_ );
	/* Look at the adjacency matrix */
	list<int> *res;
	res = new list<int>;
	for (int i=0; i < island_count_; i++) {
		if (adjmatrix_[islandID][i])
			res->push_back( i );
	}
	return res;
}


list<int>* GAIslandsTopologyDisc::getOrigins( int islandID ) {
	/* Look at the adjacency matrix */
	list<int> *res;
	res = new list<int>;
	for (int i=0; i < island_count_; i++) {
		if (adjmatrix_[i][islandID])
			res->push_back( i );
	}
	return res;
}



/* Additional methods for modifying this topology */

/**
 * Additional method for adding edges to this topology
 *
 * @param islandSource Source island
 * @param islandDest Destination island
 */
void GAIslandsTopologyDisc::addEdge( int islandSource, int islandDest ) {
  assert(islandSource < island_count_ && islandDest < island_count_);
	// Set it on the adjacency matrix
	adjmatrix_[islandSource][islandDest] = 1;
}
	
/**
 * Additional method for deleting edges in this topology
 *
 * @param islandSource Source island
 * @param islandDest Destination island
 */
void GAIslandsTopologyDisc::delEdge( int islandSource, int islandDest ) {
  assert(islandSource < island_count_ && islandDest < island_count_);
	// Set it on the adjacency matrix
	adjmatrix_[islandSource][islandDest] = 0;
}
	
int GAIslandsTopologyDisc::getNDestinations( int islandID ) {
  assert (islandID < island_count_ );
	/* Look at the adjacency matrix */
  int ndest = 0;
	for (int i=0; i < island_count_; i++) {
		if (adjmatrix_[islandID][i])  ndest++;
	}
	return ndest;
}
