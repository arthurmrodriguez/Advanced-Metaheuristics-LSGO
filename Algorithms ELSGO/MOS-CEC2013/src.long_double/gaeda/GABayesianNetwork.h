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
 * @brief GABayesianNetwork class hdr.
 *
 * Contains the definition of the GABayesianNetwork class
 */

#ifndef GABAYESIANNETWORK_H
#define GABAYESIANNETWORK_H

/* INCLUDES */
#include "GAGraphModel.h"

class GAPopulation;
class GAGenome;

// Needed when PBIL is executed
#define ALPHA (long double) 0.5


/**
 * @brief A graphical model based on bayesian networks
 *
 * This class implements a bayesian network. Several parameters (such as learning
 * method, simulation algorithm, etc..) can be retrieved and modified using
 * the get() and setptr() functions
 *
 */

class GABayesianNetwork : public GAGraphModel {

   public:

      /* DATA TYPES */

      /// Different kinds of learning methods for the bayesian network
      typedef enum {
         UMDA, EBNA_B, EBNA_LOCAL, PBIL,
         TREE, MIMIC, EBNA_K2
      } LearningMethod;

      /// Score methods for the EBNA_LOCAL learning method
      typedef enum {BIC_SCORE, K2_SCORE} EBNALocalScoring;

      /// Simulation method
      typedef enum {
         PLS, PLS_ALL_VALUES_1, PLS_ALL_VALUES_2,
         PLS_CORRECT, PENALIZATION
      } SimulationMethod;


      /* MEMBER FUNCTIONS */

      // Constructor. It creates an arcless Bayesian
      // network which represents a uniform probability
      // distribution.
      GABayesianNetwork (const GAGenome &g, LearningMethod lm, EBNALocalScoring es, SimulationMethod sm);

      // Destructor.
      ~GABayesianNetwork ();


      // Get/set parameters of the bayesian network
      int get (const char*, void*) const;
      int setptr (const char*, const void*);


      // Trains the Bayesian network from the given cases.
      void learn (const GAPopulation& p, long double sel_pctg);


      // Creates a new individual simulating the Bayesian network.
      void simulate (GAGenome& g);
      int * simulate ();


      // The BIC metric for the local structure (EBNA BIC) represented by
      // mParents [node].
      long double BICMetric (int node, int**& cases);

      // Function to execute the K2 metric for the BN (EBNA K2)
      long double K2Metric  (int node, int**& cases);

      long double deltaBIC  (int node, int nparent, int**& cases);
      long double deltaK2   (int node, int nparent, int**& cases);


   protected:

      void initialize ();


   private:

      /* VARIABLES */

      // The number of selected individuals.
      int mSelectedSize;

      // The learning type (i.e. the EDA type).
      LearningMethod mLearningMethod;

      // Scores for the EBNA-local learning types
      EBNALocalScoring mScoringMethod;

      // The simulation method
      SimulationMethod mSimulationMethod;


      // Matrixes storing the improvement of arc modifications.
      // mA[i][j] represents the gain which could be obtained
      // when adding/removing the arc j->i when it is absent/present.
      // mActualMetric[node] stores the actual value of the score of
      // the Bayesian network

      long double** mA;
      long double* mActualMetric;


      // The conditional probability distributions of the
      // bayesian network. mProbs[i][j*(r_i-1)+k] represents the
      // probability of i being in its k-th state conditiones
      // to its parent nodes being in their j-th configuration
      // (r_i is the number of possible states that i can take).

      long double** mProbs;


      // The topologycal order of the nodes. mOrder[i]
      // represents the topologycal order of i.

      int* mOrder;


      // The nodes of the Bayesian network ordered according
      // to their topologycal order.

      int* mOrdered;


      /* FUNCTIONS */

      // The logarithm of n!
      long double logFact (int n);


      // Calculates the mA matrix of the current Bayesian
      // networks.

      void calculateA (int**& cases);


      // Calculates the mA[node] values of the current
      // Bayesian network.

      void calculateANode (int node, int**& cases);


      // Sets the structure of the Bayesian network to
      // an arcless graph.

      void initBayesianNetwork ();


      // Learns the structure of the Bayesian network using
      // the local search algorithm.

      void learnEBNALocal (int**& cases);


      // Learns the structure of the Bayesian network using
      // the B algorithm.

      void learnEBNAB (int**& cases);


      // Learns the structure of the Tree by the MWST method
      // of Chow and Liu.

      void learnTree (int**& cases);


      // Learns the structure of the Bayesian network using
      // the MIMIC algorithm.

      void learnMIMIC (int**& cases);


      // Empirical entropy of node given parent.

      long double entropy (int node, int parent, int**& cases);


      // Empirical entropy of a node.

      long double entropy (int node, int**& cases);


      // Adds the arc parent->node to the Bayesian network
      // and updates all the internal data.

      void addArc (int node, int parent);


      // Removes the arc parent->node from the Bayesian
      // network and updates all the internal data.

      void removeArc (int node, int parent);


      // Compute the MI value between two variables in selected individuals

      long double computeMI (int bat, int bi, int**& cases);


      // Calculates the probabilities of the Bayesian network.

      void learnProbabilities (int**& cases /* ,long double **&values, long double sel_total */);


      // Returns the conditional probability distribution
      // of the orderede_node-th gene conditioned to the
      // states that its parent nodes have in instance.
      //    long double * Probabilities(int ordered_node,int * instance);

      long double* probabilities (int ordered_node, GAGenome& g);

};

#endif
