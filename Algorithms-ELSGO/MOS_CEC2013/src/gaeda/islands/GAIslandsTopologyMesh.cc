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
 * @brief GAIslandsTopologyMesh class impl.
 *
 * GAIslandsTopologyMesh implementation. This class represents
 * a mesh topology
 */


/* INCLUDES */
#include <list>
#include <assert.h>
#include "GAIslandsTopologyMesh.h"

using namespace std;

/**
 * Class constructor
 *
 * @param count Number of islands in the archipelago
 */
GAIslandsTopologyMesh::GAIslandsTopologyMesh( int count ) : EAIslandsTopology(count){

	/* Create the mesg topology */

	
	int i,j;
	
	// We start with a ring
	for (i=0; i < count; i++) {
		j = (i+1) % count;
		adjmatrix_[i][j] = 1;
	}

	
	// We create the edges that go across
	for (i=1; i < count/2; i++) {
		j = (count -i -1);
		
		if (i % 2 ==0) {
			adjmatrix_[i][j] = 1;
		} else {
			adjmatrix_[j][i] = 1;
		}
	}
	
	return;
}


/**
 * Class destructor
 */
GAIslandsTopologyMesh::~GAIslandsTopologyMesh() {
	/* Nothing to do here! */
}



list<int>* GAIslandsTopologyMesh::getDestinations( int islandID ) {

	/* Look at the adjacency matrix */
	list<int> *res;
	res = new list<int>;
	for (int i=0; i < island_count_; i++) {
		if (adjmatrix_[islandID][i])
			res->push_back( i );
	}
	return res;
}


list<int>* GAIslandsTopologyMesh::getOrigins( int islandID ) {

	/* Look at the adjacency matrix */
	list<int> *res;
	res = new list<int>;
	for (int i=0; i < island_count_; i++) {
		if (adjmatrix_[i][islandID])
			res->push_back( i );
	}
	return res;

}

