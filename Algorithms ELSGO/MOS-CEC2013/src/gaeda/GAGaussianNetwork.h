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
 * @brief GAGaussianNetwork class hdr.
 *
 * Contains the definition of the GAGaussianNetwork class
 */

#ifndef GAUSSIANNETWORK_H
#define GAUSSIANNETWORK_H

/* INCLUDES */
#include <set>

#include "GAGraphModel.h"

#ifndef CSet
#define CSet std::set<int>
#endif


/**
 * @brief A graphical model based on gaussian networks
 *
 * This class implements a gaussian network. Several parameters (such as learning
 * method, simulation algorithm, etc..) can be retrieved and modified using
 * the get() and setptr() functions
 *
 */

class GAGaussianNetwork : public GAGraphModel {

   public:

      /* DATA STRUCTURES */
      typedef enum {BGe_SCORE = 0,BIC_SCORE = 1} ScoreMethod;
      typedef enum {UMDA, EGNA_B, EGNA_LOCAL, MIMIC, EE, EMNA} LearningMethod;


      /* FUNCTIONS */

      // The constructor. It creates an arcless Gaussian network
      // which represents a uniform probability distribution.
      GAGaussianNetwork (const GAGenome& g, LearningMethod lm, ScoreMethod sm);

      // The destructor.
      virtual ~GAGaussianNetwork ();

      // Get/set network configuration options
      int get (const char* var, void* value) const;
      int setptr (const char* var, const void* value);

      // Learns the Gaussian network from the given cases.
      void learn (const GAPopulation& p, double sel_pctg);

      // Creates a new individual simulating the Gaussian network.
      void simulate (GAGenome& g);
      int* simulate ();


      // Functions for EMNA algorithms

      // Initialises the Gaussian Network for EMNA algorithms
      void initEMNAGaussianNetwork (const GAPopulation& cases);

      // To improve a new cycle in incremental and adaptative EMNA algorithms
      void EMNAAdapt(const GAPopulation& cases, const GAGenome& genes_new,
                     const GAGenome& genes_worse);


   protected:

      void initialize ();


   private:

      /* VARIABLES */

      // Network parameters
      ScoreMethod mScoringMethod;
      LearningMethod mLearningMethod;

      double* mMin;
      double* mMax;
      int mSelectionSize;
      double mNrPrecision;
      int mNrMax;

      // Related to the BIC score
      long double** m_Sxixj;
      long double*  m_S2xi;

      // reference to the database of cases.
      const GAPopulation* m_cases;

      // Unconditional means of the variables.
      long double* mMeans;

      // v values of the gaussian network.
      long double* mV;

      // b values of the gaussian network. m_b[i][j] stores
      // the b value corresponding to the arc j->i.
      long double** mB;

      // The topological order of the nodes. m_order[i]
      // represents the topologycal order of i.
      int* mTopSortedMap;

      // The nodes of the Gaussian network sorted according
      // to their topologycal order.
      int* mTopSorted;

      // Matrix storing the improvement of arc modifications.
      // mA[i][j] represents the gain which could be obtained
      // when adding/removing the arc j->i when it is absent/present.
      double** mA;

      // Variance-covariance matrix.
      long double** mSigma;

      // mT0
      double mT0;


      /* FUNCTIONS */

      double f (double x);
      double F (double x);

      double rejectionRegionBoundary ();
      void learnEE ();

      double logGamma (double num);
      double logC (int l, int alpha);

      double BGe (CSet& nodes);
      double BGe (int node);


      // functions for the BIC score
      double BIC (int node, CSet& nodes);
      double BIC (int node);
      void calculateStructuresforBIC (const GAPopulation& cases);
      void distributionParametersBIC ();

      // The empirical entropy of node given parent.
      double entropy (int node, int parent);

      // The empirical entropy of node.
      double entropy (int node);

      // It calculates the mA matrix of the current Gaussian networks.
      void calculateA ();

      // It calculates the mA[node] values of the current Gaussian network.
      void calculateANode (int node);

      // It sets the structure of the Gaussian network to
      // an arcless or complete graph.
      void initGaussianNetwork (bool complete = false);

      // It learns the structure of the Gaussian network using
      // the local search algorithm.
      void learnEGNALocal ();

      // It learns the structure of the Gaussian network using the B algorithm.
      void learnEGNAB ();

      // It learns the structure of the Gaussian network as in the MIMIC algorithm.
      void learnMIMIC ();

      // It calculates the probabilities of the Gaussian network.
      void learnDistributions ();

      // It returns the mean and variance of the conditional probability
      // distribution of the ordered_node-th gene conditioned to the values
      // that its parent nodes have in instance.
      void distribution (int ordered_node, GAGenome& instance,
                         double& mean, double& variance);

      // It adds the arc parent->node to the Gaussian network
      // and updates all the internal data.
      void addArc (int node, int parent);

      // It removes the arc parent->node from the Gaussian
      // network and updates all the internal data.
      void removeArc (int node, int parent);

      // It calculates the means and the covariance matrix.
      void calculateMeansAndCovariances (const GAPopulation& cases);

};

#endif
