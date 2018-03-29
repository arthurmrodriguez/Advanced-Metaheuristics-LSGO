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
 * @brief GAGraphModel class hdr.
 *
 * Contains the definition of the GAGraphModel class
 */

#ifndef GAGRAPHMODEL_H
#define GAGRAPHMODEL_H

/* INCLUDES */

#include "gaid.h"

class GAGenome;
class GAPopulation;

/* CLASS DECLARATION */

/**
 * @brief Base class for graphical models
 *
 * Abstract base class for the different graphical models
 * used in EDA algorithms. Currently two subclasses (BayesianNetwork
 * and GaussianNetwork) implement actual network models.
 */

class GAGraphModel {

   public:

      /* DATA TYPES */

      enum ModelType {NONE = 0, BAYESIAN_NETWORK = 1, GAUSSIAN_NETWORK = 2};


      /* MEMBER FUNCTIONS */

      // Constructor
      GAGraphModel (const GAGenome& g);

      // Destructor
      virtual ~GAGraphModel ();

      // Learn from the model
      virtual void learn (const GAPopulation& p, double sel_pctg) = 0;

      // Simulate
      virtual void simulate (GAGenome& g) = 0;

      // Get/set network-specific parameters
      virtual int get    (const char*, void*) const = 0;
      virtual int setptr (const char*, const void*) = 0;

      // Guess the graphical model from the genome
      static ModelType GuessModelType (const GAGenome& g);

      // Get the individual's size
      int getIndividualSize (void) const {return mIndividualSize;}

      // Get the value of one of the states
      int getState (int pos) const {return (pos < mIndividualSize) ? mStates [pos] : -1 ;}

      // Get the genome
      const GAGenome* getGenome (void) const {return &mGenome;}

      // Get the kind of model
      ModelType getModelType (void) const {return mModelType;}


   protected:

      virtual void initialize () = 0;

      /// Adjacency matrix representing the structure of
      /// the network. If mParents[i][j] is true
      /// the there is an arc j->i in the network.
      bool** mParents;

      /// The matrix which stores the number of paths between
      /// two nodes. mPaths[i][j] represents the number of
      /// paths from j to i. mPaths[i][i] is always 1.
      int** mPaths;


   private:

      // Size of the individual
      int mIndividualSize;

      // States
      int* mStates;

      // Genome
      const GAGenome& mGenome;

      // Model type
      ModelType mModelType;

};

#endif
