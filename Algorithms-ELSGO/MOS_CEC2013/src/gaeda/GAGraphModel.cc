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

/* INCLUDES */

#include "GAGraphModel.h"

#include "GAPopulation.h"
#include "genomes/GAGenome.h"

/* FUNCTIONS */

/**
 * Static function that guesses which kind of graphical model should be used
 * given an population individual
 *
 * @param g Genome used to guess the network type
 *
 */

GAGraphModel::ModelType GAGraphModel::GuessModelType (const GAGenome& g) {

   switch (g.classID ()) {

      // Binary strings, integer and such -> discrete EDA -> Bayesian network
      case GAID::BinaryStringGenome:
      case GAID::Bin2DecGenome:
      case GAID::BinaryStringGenome2D:
      case GAID::BinaryStringGenome3D:
      case GAID::StringGenome:
      case GAID::IntGenome:

         return BAYESIAN_NETWORK;
         break;

      // Float and double variables -> cont. EDA -> Gaussian network
      case GAID::FloatGenome:
      case GAID::DoubleGenome:

         return GAUSSIAN_NETWORK;

      // Arrays -> messy stuff.
      case GAID::ArrayAlleleGenome:
      case GAID::ArrayAlleleGenome2D:
      case GAID::ArrayAlleleGenome3D:

         ModelType t = NONE;

         for (int i = 0; i < g.fixedSize (); i++)
            if (g.domain (i) > 0)
               if (t == GAUSSIAN_NETWORK) //Hybrid
                  return NONE;
               else
                  t = BAYESIAN_NETWORK;
            else
               if (t == BAYESIAN_NETWORK) //Hybrid
                  return NONE;
               else
                  t = GAUSSIAN_NETWORK;

         return t;

   }

   return NONE;

}


/**
 * Constructor. Creates an arcless network with no paths among
 * nodes
 */

GAGraphModel::GAGraphModel (const GAGenome& g) : mGenome (g) {

   // Set the individual size, type and such

   mIndividualSize = g.fixedSize ();
   mModelType = GuessModelType (g);

   // Create the state array
   if (mIndividualSize > 0) {

      mStates = new int [mIndividualSize];

      for (int i = 0; i < mIndividualSize; i++)
         mStates [i] = g.domain (i);

   }

   // An arcless network is created. All nodes have
   // no parents.
   mParents = new bool* [mIndividualSize];
   int i, j;

   for (i = 0; i < mIndividualSize; i++) {

      mParents [i] = new bool [mIndividualSize];

      for (j = 0; j < mIndividualSize; j++)
         mParents [i][j] = false;

   }

   // Thus, there are no paths among nodes.
   mPaths = new int* [mIndividualSize];

   for (i = 0; i < mIndividualSize; i++) {

      mPaths [i] = new int [mIndividualSize];

      for (j = 0; j < mIndividualSize; j++)
         if (i == j)
            mPaths [i][j] = 1;
         else
            mPaths [i][j] = 0;

   }

}


/**
 * Class destructor
 */

GAGraphModel::~GAGraphModel () {

   // Deallocate the states array
   delete [] mStates;

}
