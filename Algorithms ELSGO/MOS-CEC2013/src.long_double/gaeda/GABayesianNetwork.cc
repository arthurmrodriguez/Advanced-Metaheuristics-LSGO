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

// GABayesianNetwork.cc: implementation of the GABayesianNetwork class.

/* INCLUDES */
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <assert.h>

#include "GABayesianNetwork.h"

#include "chi.h"
#include "garandom.h"
#include "gaparameters.h"
#include "GAPopulation.h"
#include "genomes/GAGenome.h"

// The constructor. It creates an arcless Bayesian
// network which represents a uniform probability
// distribution.

GABayesianNetwork::GABayesianNetwork (const GAGenome& g, LearningMethod lm,
                                      EBNALocalScoring es, SimulationMethod sm) : GAGraphModel (g) {

   // The number of selected individuals.
   mSelectedSize = 0;

   // Learning method
   mLearningMethod = lm;

   // Scores for the EBNA-local learning types
   mScoringMethod = es;

   // Simulation method
   mSimulationMethod = sm;

   initialize ();

}


// Get parameters of the bayesian network

int GABayesianNetwork::get (const char* var, void* value) const {

   int status = 1;

   if (strcmp (var, gaNdiscreteLearningMethod ) == 0 ||
       strcmp (var, gaSNdiscreteLearningMethod) == 0)   {

      LearningMethod *ptr = (LearningMethod*) value;
      *ptr = mLearningMethod;
      status = 0;

   }
   else if (strcmp (var, gaNdiscreteEBNAScoring ) == 0 ||
            strcmp (var, gaSNdiscreteEBNAScoring) == 0)   {

      EBNALocalScoring *ptr = (EBNALocalScoring*) value;
      *ptr = mScoringMethod;
      status = 0;

   }
   else if (strcmp (var, gaNdiscreteSimulationMethod ) == 0 ||
            strcmp (var, gaSNdiscreteSimulationMethod) == 0)   {

      SimulationMethod *ptr = (SimulationMethod*) value;
      *ptr = mSimulationMethod;
      status = 0;

   }

   return status;

}


// Set parameters of the bayesian network

int GABayesianNetwork::setptr (const char* var, const void* value) {

   int status = 1;

   if (strcmp (var, gaNdiscreteLearningMethod ) == 0 ||
       strcmp (var, gaSNdiscreteLearningMethod) == 0)   {

      LearningMethod *ptr = (LearningMethod*) value;
      mLearningMethod = *ptr;
      status = 0;

   }
   else if (strcmp (var, gaNdiscreteEBNAScoring ) == 0 ||
            strcmp (var, gaSNdiscreteEBNAScoring) == 0)   {

      EBNALocalScoring *ptr = (EBNALocalScoring*) value;
      mScoringMethod = *ptr;
      status = 0;

   }
   else if (strcmp (var, gaNdiscreteSimulationMethod ) == 0 ||
            strcmp (var, gaSNdiscreteSimulationMethod) == 0)   {

      SimulationMethod *ptr = (SimulationMethod*) value;
      mSimulationMethod = *ptr;
      status = 0;

   }

   return status;

}


void GABayesianNetwork::initialize () {

   int i,j;

   mSelectedSize = 0;

   // The probability distribution will be uniform.
   mProbs = new long double * [getIndividualSize ()];

   for (i = 0; i < getIndividualSize (); i++) {

      mProbs [i] = new long double [getState (i) - 1];

      // The final probability is found subtracting
      // all the others to 1.

      for (j = 0; j < getState (i) - 1; j++)
         mProbs [i][j] = 1 / (long double) getState (i);

   }

   // The order is trivial.
   mOrder = new int [getIndividualSize ()];
   for (i = 0; i < getIndividualSize (); i++)
      mOrder [i] = i;

   // The ordered nodes can be found from mOrder.
   mOrdered = new int [getIndividualSize ()];
   for (i = 0; i < getIndividualSize (); i++)
      mOrdered [mOrder [i]] = i;

   // memory for mA is allocated.
   mA = new long double * [getIndividualSize ()];
   for (i = 0; i < getIndividualSize (); i++) {

      mA [i] = new long double [getIndividualSize ()];

      for (j = 0; j < getIndividualSize (); j++)
         mA [i][j] = INT_MIN;

   }

   // Matrix to record the actual state of the Bayesian network.
   mActualMetric = new long double [getIndividualSize ()];

}


GABayesianNetwork::~GABayesianNetwork () {

   int i;

   for (i = 0; i < getIndividualSize (); i++)
      delete [] mProbs [i];

   delete [] mProbs;

   delete [] mOrder;
   delete [] mOrdered;

   for (i = 0; i < getIndividualSize (); i++)
      delete [] mParents [i];

   delete [] mParents;

   for (i = 0; i < getIndividualSize (); i++)
      delete [] mA [i];

   delete [] mA;

   for (i = 0; i < getIndividualSize (); i++)
      delete [] mPaths [i];

   delete [] mPaths;

   delete [] mActualMetric;

}


// Creates a new individual simulating the Bayesian network.

void GABayesianNetwork::simulate (GAGenome& g) {

   int i, j, v, variable;

   // Variables to modify the simulation step
   // we assume that, when using PLS_ALL_VALUES_x,
   // the number of states of all variables is always the same
   int* Indiv = new int [getState (0)];

   for (i = 0; i < getState (0); i++)
      Indiv [i] = 0;

   int variables_left = getIndividualSize ();
   int values_left = getState (0);

   long double P_indiv, K;
   long double sum;

   // table to store all the probabilities and access them more easily
   int max_states = 0; // variable to compute which is the maximum number of states of any of the variables

   for (i = 0; i < getIndividualSize (); i++)
      if (max_states < getState (i))
         max_states = getState (i);

   long double* Prob_table = new long double [max_states];


   // The individual will be generated simulating its genes according
   // to the order they have in the Bayesian network.

   for (variable = 0; variable < getIndividualSize (); variable++) {

      long double which = GARandomDouble ();

      // Calculate the probabilities and store them in the probability table
      long double* probs = probabilities (variable, g);

      sum = 0.0;
      for (i = 0; i < getState (variable) - 1; i++) {

         Prob_table [i] = probs [i];
         sum += probs [i];

      }

      Prob_table [getState (variable) - 1] = 1 - sum; // prob. of last variable = 1 - sum of the rest

      // Modify the probabilities to get correct individuals
      switch (mSimulationMethod) {

         case PLS:

            break; // PLS simple, no modifications

         case PLS_CORRECT:

            break; // PLS simple, modifications later

         case PLS_ALL_VALUES_1:

            if (variables_left == values_left) {

               P_indiv = 0.0;

               for (i = 0; i < getState (variable); i++)
                  if (Indiv [i] > 0)
                     P_indiv += Prob_table [i];

               for (i = 0; i < getState (variable); i++) {

                  if (Indiv [i] > 0)
                     Prob_table [i] = 0.0;
                  else
                     Prob_table [i] *= 1 / (1 - P_indiv);

               }

            }

            break;

         case PLS_ALL_VALUES_2:

            if (values_left > 0) {

               // calculate P_indiv
               P_indiv = 0.0;

               for (i = 0; i < getState (variable); i++)
                  if (Indiv [i] > 0)
                     P_indiv += Prob_table [i];

               // Change the probabilities of all values
               if (variables_left == values_left)
                  for (i = 0; i < getState (variable); i++)
                     if (Indiv [i] > 0)
                        Prob_table [i] = 0.0;
                     else
                        Prob_table [i] *= 1 / (1 - P_indiv);
               else {
                  // calculate K, the constant to divide the probabilities
                  K = int ((getIndividualSize () - variables_left) / (variables_left - values_left)) + 1;
                  // K = round (K);

                  for (i = 0; i < getState (variable); i++)
                     if (Indiv [i] > 0)
                        Prob_table [i] /= K;
                     else
                        Prob_table [i] *= (K - P_indiv) / (K * (1 - P_indiv));

               }

            }

         default:

            break;
      }

      int j;
      for (j = 0; j < (getState (mOrdered [variable]) - 1) && which > Prob_table [j]; j++)
         which -= Prob_table [j];

      // Correct the problem of the rounding in the computer
      if (Prob_table [j] == 0.0) {
         if (j == 0) // we will select the next variable that has a prob!=0
            do {
               j++;
            } while (Prob_table [j] == 0.0);
         else // we will select the previous variable that has a prob!=0
            do {
               j--;
            } while (Prob_table [j] == 0.0);
      }

      g.setValueOfNominalGene (mOrdered [variable], j);
      variables_left--;

      if (Indiv [j] == 0)
         values_left--;

      Indiv [j]++;

   }

   // Correct the individual if necessary and requested
   if ((mSimulationMethod == PLS_CORRECT) && (values_left > 0)) {

      while (values_left > 0) {

         // select a random position that contains a variable already appeared at least twice
         do {
            j = GARandomInt (0, getIndividualSize () - 1);
         } while (Indiv [g.getValueOfNominalGene (j)] < 2);

         // substitute this variable by next the missing value.
         //TODO: revisar esto
         for (v = 0; v < getIndividualSize (); v++)
            if (Indiv [v] == 0)
               break;

         //for (v=0; v<getState(variable); v++) {if (Indiv[v]==0) break;}

         Indiv [g.getValueOfNominalGene (j)]--;
         g.setValueOfNominalGene (j, v);
         Indiv [v]++;
         values_left--;

      }

   }

   delete [] Indiv;
   delete [] Prob_table;

}



// It learns the Bayesian network from the given cases.

void GABayesianNetwork::learn (const GAPopulation& p, long double sel_pctg) {

   // BEGIN: SUPER APISONATOR!!
   int i, I, g, G;

   I = p.size ();
   G = getIndividualSize ();
   mSelectedSize = (int) (I * sel_pctg);

   // GCC 4.1 (fuera los paréntesis del tipo)
   int** cases = new int* [mSelectedSize];

   for (i = 0; i < mSelectedSize; i++) {

      // GCC 4.1 (fuera los paréntesis del tipo)
      cases [i] = new int [G];

      for (g = 0; g < G; g++)
         cases [i][g] = p.individual (i).getValueOfNominalGene (g);

   }
   // END: SUPER APISONATOR!!

   // Learn the structure
   switch(mLearningMethod) {

      case UMDA:

         // No learning required.
         break;

      case EBNA_B:

         mScoringMethod = BIC_SCORE;
         learnEBNAB (cases);

         break;

      case EBNA_LOCAL:

         mScoringMethod = BIC_SCORE;
         learnEBNALocal (cases);

         break;

      case EBNA_K2:

         mScoringMethod = K2_SCORE;
         learnEBNALocal (cases);

         break;

      case PBIL:

         // No learning required.
         break;

      case TREE:

         learnTree (cases);
         break;

      case MIMIC:

         learnMIMIC (cases);
         break;

   }

   // Learn the probabilities
   learnProbabilities (cases);

   for (i = 0; i < mSelectedSize; i++)
      delete cases [i];

   delete [] cases;

}


long double* GABayesianNetwork::probabilities (int ordered_node, GAGenome& g) {

   // First the configuration of the parent nodes is found.
   int parent_configuration = 0;
   int i;

   for (i = 0; i < ordered_node; i++)
      if (mParents [mOrdered [ordered_node]][mOrdered [i]]) {
         parent_configuration *= getState (mOrdered [i]);
         parent_configuration += g.getValueOfNominalGene (mOrdered [i]);
      }

   // Then the conditional probability corresponding to that
   // parent configuration is returned.
   return &(mProbs [mOrdered [ordered_node]][parent_configuration * (getState (mOrdered [ordered_node]) - 1)]);

}


void GABayesianNetwork::learnProbabilities (int**& cases) {

   long i, j, k;
   int** nijk = NULL;

   for (i = 0; i < getIndividualSize (); i++) {

      // Calculate the number of parent configurations of i.
      long no_j = 1;

      for (j = 0; j < getIndividualSize (); j++)
         if (mParents [i][j]) no_j *= getState (j);

      // Allocate memory for all nijk-s.
      nijk = NULL;
      nijk = new int*[no_j];

      for (j = 0; j < no_j; j++) {

         nijk [j] = new int [getState (i)];

         for (k = 0; k < getState (i); k++)
            nijk [j][k] = 0;

      }

      // All nijk-s are calculated.
      for (j = 0; j < mSelectedSize; j++) {

         int parent_configuration = 0;

         for (int parent = 0; parent < getIndividualSize (); parent++)
            if (mParents [i][parent]) {

               parent_configuration *= getState (parent);
               parent_configuration += cases [j][parent];

            }

         nijk [parent_configuration][cases [j][i]]++;

      }

      // All probabilities are recalculated (as Buntine).

      if (mLearningMethod != PBIL) {

         delete [] mProbs [i];
         mProbs [i] = new long double [no_j * (getState (i) - 1)];

      }

      // PBIL needs the probabilities of the prior population.
      // Other learning algorithms do not need them.

      for (j=0;j<no_j;j++) {

         int nij = 0;

         for (k = 0; k < getState (i); k++)
            nij += nijk [j][k];

         long double BDeu = (long double) 1 / getState (i) / no_j;

         // Also when the BN is induced by EDA, conditional
         // probabilities are calculated in univariated=UMDA form
         if ((mLearningMethod != PBIL) )
            for (k = 0; k < getState (i) - 1; k++)
               mProbs [i][j * (getState (i) - 1) + k] = (long double) (nijk [j][k] + BDeu) / (long double) (nij + getState (i) * BDeu);

         if (mLearningMethod == PBIL)
            for (k = 0; k < getState (i) - 1; k++)
               mProbs [i][j * (getState (i) - 1) + k] = (1.0 - ALPHA) * (mProbs [i][j * (getState (i) - 1) + k]) +
                                                        ALPHA * ((long double) (nijk [j][k]) / (long double) (nij));

      }

      // Memory used for nijk-s is released.
      for (j = 0; j < no_j; j++)
         delete [] nijk [j];

      delete [] nijk;

   }

}


/* ******************************* *
 * EBNA-learning related functions *
 * ******************************* */

void GABayesianNetwork::learnEBNAB (int**& cases) {

   initBayesianNetwork ();

   learnEBNALocal (cases);

}


void GABayesianNetwork::learnEBNALocal (int**& cases) {

   int TimesArcsChanged = 0;
   long double max;

   do {

      TimesArcsChanged++; // To count the number of arcs removed/added

      // Find the arc which modification most increases
      // the BIC metric of the Bayesian network.

      int max_i = 0;
      int max_j = 0;
      max = INT_MIN;

      for (int i = 0; i < getIndividualSize (); i++)
         for (int j = 0; j < getIndividualSize (); j++)
            if (mPaths [j][i] == 0 && mA [i][j] > max) {

               // In order to modify the arc j->i
               // there cannot be a path from i to j.
               max_i = i;
               max_j = j;
               max = mA [i][j];

            }

      // If the BIC metric of the Bayesian network can
      // be improved the arc is modified.

      if (max > 0) {

         if (!mParents [max_i][max_j])
            // If there is no arc j->i, it is added.
            addArc (max_i, max_j);
         else
            // If there is an arc j->i, it is removed.
            removeArc (max_i, max_j);

         // Update the actual BIC or K2 metric
         mActualMetric [max_i] += mA [max_i][max_j];

         calculateANode (max_i, cases);

      }

   } while (max > 0);

}


/* ************* *
 * Tree-learning *
 * ************* */

void GABayesianNetwork::learnTree (int**& cases) {

   int i, j, k, l;

   int withnocycles = 0;
   int max_i = 0;
   int max_j = 0;
   long double max_mi = 0.0;
   bool no_info = false; // to check if the arc with the best MI has
                         // MI = 0.0--> in that case, no continue
                         // learning the tree structure.

   // An up triangular matrix which contains the MI (X,Y) values
   // of the variables of the program

   long double** mi = new long double* [getIndividualSize ()];

   for (j = 0; j < getIndividualSize (); j++) {

      mi [j] = new long double [getIndividualSize ()];

      for (int k = 0; k < getIndividualSize (); k++)
         mi [j][k] = (long double) 0.0;

   }

   // A boolean long double matrix containing if there is an undirected
   // path between two nodes between i and i there exists an undirected path

   bool** m_upaths = new bool* [getIndividualSize ()];

   for (j = 0; j < getIndividualSize (); j++) {

      m_upaths [j] = new bool [getIndividualSize ()];

      for (l = 0; l < getIndividualSize (); l++)
         m_upaths [j][l] = false;

   }

   for (l = 0; l < getIndividualSize (); l++)
      m_upaths [l][l] = true;

   for (i = 0; i < getIndividualSize (); i++)
      for (j = 0; j < getIndividualSize (); j++)
         if (j > i)
            mi [i][j] = computeMI (i, j, cases);

   initBayesianNetwork ();

   // Search for the biggest MI and go inducing the BN=Tree of
   // getIndividualSize()-1 arcs

   while((withnocycles<(getIndividualSize()-1)) && (no_info == false)) {

      max_mi = (long double)0.0;

      for (i = 0; i < getIndividualSize (); i++)
         // no undirected path between i and j
         for (j = 0; j < getIndividualSize (); j++)
            if (j > i && m_upaths [j][i] == false && m_upaths [i][j] == false && mi [i][j] >= max_mi) {

               max_i  = i;
               max_j  = j;
               max_mi = mi [i][j];

            }


      // The best arc has MI = 0 and the tree learn process must be stopped
      if (max_mi == (long double) 0.0)
         no_info = true;

      mi [max_i][max_j] = (long double) (0.0);

      // check if the MI of the 'best' arc is not 0.0;-->in this case, no continue learning the tree.
      // check if the edge can be added not to create cycles:
      // no directed paths from i to j, j is no parent of i;
      // no already arc j-->i.

      if ((no_info == false) && (!mParents [max_i][max_j]) && (mPaths [max_j][max_i] == 0)) {

         withnocycles++;

         // From j to i: j-->i
         mParents [max_i][max_j] = true;

         // The undirected paths of the Bayesian
         // network are updated. The nodes that can be accesed
         // from i, they can be now be accesed from j, and viceversa.

         m_upaths [max_i][max_j] = true;
         m_upaths [max_j][max_i] = true;

         for (l = 0; l < getIndividualSize (); l++) {

            if (m_upaths [l][max_i] == true || m_upaths [max_i][l] == true) {

               m_upaths [l][max_j] = true;
               m_upaths [max_j][l] = true;

            }

            if (m_upaths [l][max_j] == true || m_upaths [max_j][l] == true) {

               m_upaths [l][max_i] = true;
               m_upaths [max_i][l] = true;

            }

         }

         // Update the undirected paths of inderectly implicated nodes.
         for (l = 0; l < getIndividualSize (); l++) {

            if (m_upaths [l][max_i] == true || m_upaths [max_i][l] == true)
               for (int r = 0; r < getIndividualSize (); r++)
                  if (m_upaths [r][max_i] == true || m_upaths [max_i][r] == true) {

                     m_upaths [l][r] = true;
                     m_upaths [r][l] = true;

                  }

            if (m_upaths [l][max_j] == true || m_upaths [max_j][l] == true)
               for (int r = 0; r < getIndividualSize (); r++)
                  if (m_upaths [r][max_j] == true || m_upaths [max_j][r] == true) {

                     m_upaths [l][r] = true;
                     m_upaths [r][l] = true;

                  }


         }

         // The paths of the Bayesian network are updated.

         for (i = 0; i < getIndividualSize (); i++)
            for (j = 0; j < getIndividualSize (); j++)
               if (mPaths [max_j][j] > 0 && mPaths [i][max_i] > 0)
                  mPaths [i][j] += mPaths [max_j][j] * mPaths [i][max_i];


         // Update the ordering of the nodes.

         if (mOrder [max_i] < mOrder [max_j]) {

            // The order of max_i, its descendants, max_j
            // and its ancestors must be updated.

            int jump_j = 0; // How many positions the ancestors of max_j are moved.

            // Calculate how many positions the descendants of max_i must be moved.
            int jump_i = 0;

            for (k = mOrder [max_i]; k <= mOrder [max_j]; k++)
               if (mPaths [max_j][mOrdered [k]] > 0)
                  jump_i++;

            // Update the order of the nodes between max_i and max_j (both included).

            for (k = mOrder [max_i]; k <= mOrder [max_j]; k++)
               if (mPaths [max_j][mOrdered [k]] > 0) {

                  mOrder [mOrdered [k]] += jump_j;
                  jump_i--;

               }
               else {

                  mOrder [mOrdered [k]] += jump_i;
                  jump_j--;

               }


            // Update the ordered nodes.

            for (k = 0; k < getIndividualSize (); k++)
               mOrdered [mOrder [k]] = k;

         }

      } // del if de mPaths, mParents

   }

   // free memory

   for (i = 0; i < getIndividualSize (); i++)
      delete [] mi [i];

   delete [] mi;

   for (i = 0; i < getIndividualSize (); i++)
      delete [] m_upaths [i];

   delete [] m_upaths;

}


long double GABayesianNetwork::computeMI (int bat, int bi, int**& cases) {

   long double pbat   = (long double) 0.0;
   long double pbi    = (long double) 0.0;
   long double pbatbi = (long double) 0.0;
   long double mi     = (long double) 0.0;

   for (int i = 0; i < getState (i); i++)      // for bat
      for (int j = 0; j < getState (j); j++) { // for bi

         pbat   = (long double) 0.0;
         pbi    = (long double) 0.0;
         pbatbi = (long double) 0.0;

         for (int r = 0; r < mSelectedSize; r++) {

            if (cases [r][bat] == i)
               pbat = pbat + (long double) 1.0;

            if (cases [r][bi] == j)
               pbi = pbi + (long double) 1.0;

            if ((cases [r][bat] == i) && (cases [r][bi] == j))
               pbatbi = pbatbi + (long double) 1.0;

         }

         pbat   = pbat   / (long double) (mSelectedSize);
         pbi    = pbi    / (long double) (mSelectedSize);
         pbatbi = pbatbi / (long double) (mSelectedSize);

         if (pbatbi != (long double) 0.0 && pbat != (long double) 0.0 && pbi != (long double) 0.0)
            mi += pbatbi * (log ((long double) ((pbatbi) / (pbat * pbi))));

      }


   return mi;

}

/* ************** *
 * MIMIC learning *
 * ************** */

void GABayesianNetwork::addArc (int node, int parent) {

   int i, j, k;

   // Add the arc.
   mParents [node][parent] = true;

   // The paths of the Bayesian network are updated.
   for (i = 0; i < getIndividualSize (); i++)
      for (j = 0; j < getIndividualSize (); j++)
         if (mPaths [parent][j] > 0 && mPaths [i][node] > 0)
            mPaths [i][j] += mPaths [parent][j] * mPaths [i][node];

   // Update the ordering of the nodes.
   if (mOrder [node] < mOrder [parent]) {

      // The order of node, its descendants, parent
      // and its ancestors must be updated.

      int jump_parent = 0; // How many positions the ancestors of parent are moved.

      // Calculate how many positions the descendants of node must be moved.
      int jump_node = 0;

      for (k = mOrder [node]; k <= mOrder [parent]; k++)
         if (mPaths [parent][mOrdered [k]] > 0)
            jump_node++;

      // Update the order of the nodes between
      // node and parent (both included).

      for (k = mOrder [node]; k <= mOrder [parent]; k++)
         if (mPaths [parent][mOrdered [k]] > 0) {
            mOrder [mOrdered [k]] += jump_parent;
            jump_node--;
         }
         else {
            mOrder [mOrdered [k]] += jump_node;
            jump_parent--;
         }

      // Update the ordered nodes.
      for (k = 0; k < getIndividualSize (); k++)
         mOrdered [mOrder [k]] = k;

   }

}


long double GABayesianNetwork::entropy (int node, int**& cases) {

   int x, i;
   int* n_x = new int [getState (node)];

   for (x = 0; x < getState (node); x++)
      n_x [x] = 0;

   for (i = 0; i < mSelectedSize; i++)
      n_x [cases [i][node]]++;

   long double entropy = 0;

   for (x = 0; x < getState (node); x++)
      if (n_x [x] > 0)
         entropy += n_x [x] * log ((long double) n_x [x]);

   entropy /= -mSelectedSize;
   entropy += log ((long double) mSelectedSize);

   delete [] n_x;

   return entropy;

}


long double GABayesianNetwork::entropy (int node, int parent, int**& cases) {

   int x, i, y;
   int** n_xy = new int* [getState (node)];

   for (x = 0; x < getState (node); x++)
      n_xy[x] = new int [getState (parent)];

   for (x = 0; x < getState (node); x++)
      for (int y = 0; y < getState (parent); y++)
         n_xy [x][y] = 0;

   for (i = 0; i < mSelectedSize; i++)
      n_xy [cases [i][node]][cases [i][parent]]++;

   long double entropy = 0;

   for (x = 0; x < getState (node); x++) {

      long double entropy_y = 0;
      int n_x = 0;

      for (y = 0; y < getState (parent); y++)
         if (n_xy [x][y] > 0) {

            entropy_y += n_xy [x][y] * log ((long double) n_xy [x][y]);
            n_x += n_xy [x][y];

         }

      entropy -= entropy_y;

      if (n_x > 0)
         entropy += n_x * log ((long double) n_x);

   }

   entropy /= -mSelectedSize;
   entropy -= log ((long double) (mSelectedSize));

   for (x = 0; x < getState (node); x++)
      delete [] n_xy [x];

   delete [] n_xy;

   return entropy;

}


void GABayesianNetwork::learnMIMIC (int**& cases) {

   int i;
   initBayesianNetwork ();

   bool* linked = new bool [getIndividualSize ()];

   for (i = 0; i < getIndividualSize (); i++)
      linked [i] = false;

   int min_i = -1;
   long double min_entropy_i = INT_MAX;

   for (i = 0; i < getIndividualSize (); i++) {

      long double entropy_i = entropy (i, cases);

         if (entropy_i < min_entropy_i) {

            min_i = i;
            min_entropy_i = entropy_i;

         }

   }

   linked [min_i] = true;
   int head = min_i;

   for (i = getIndividualSize () - 2; i >= 0; i--) {

      int min_j = -1;
      long double min_entropy_j = INT_MAX;

      for (int j = 0; j < getIndividualSize (); j++)
         if (!linked [j]) {

            long double entropy_j = entropy (j, head, cases);

            if (entropy_j < min_entropy_j) {

               min_j = j;
               min_entropy_j = entropy_j;

            }

         }

      addArc (min_j, head);
      linked [min_j] = true;
      head = min_j;

   }

   delete [] linked;

}


void GABayesianNetwork::initBayesianNetwork () {

   int i, j;

   for (i = 0; i < getIndividualSize (); i++) {

      mOrder   [i] = i;
      mOrdered [i] = i;

      for (j = 0; j < getIndividualSize (); j++) {

         mParents [i][j] = false;

         if (i == j)
            mPaths [i][j] = 1;
         else
            mPaths [i][j] = 0;

      }

   }

}


void GABayesianNetwork::calculateA (int**& cases) {

   for (int node = 0; node < getIndividualSize (); node++) {

      switch (mScoringMethod) {

         case BIC_SCORE:

            mActualMetric [node] = BICMetric (node, cases);
            break;

         case K2_SCORE:

            mActualMetric [node] = K2Metric (node, cases);
            break;

      }

      calculateANode (node, cases);

   }

}


void GABayesianNetwork::calculateANode (int node, int**& cases) {

   for (int i = 0; i < getIndividualSize (); i++)
      if ((i != node)) {

         mParents [node][i] = !mParents [node][i];

         long double new_metric = 0;

         switch (mScoringMethod) {

            case BIC_SCORE:

               new_metric = BICMetric (node, cases);
               break;

            case K2_SCORE:

               new_metric = K2Metric (node, cases);
               break;

         }

         mParents [node][i] = !mParents [node][i];
         mA [node][i] = new_metric - mActualMetric [node];

      }
      else
         mA [node][i] = INT_MIN;

}


long double GABayesianNetwork::BICMetric (int node, int**& cases) {

   int j, k;

   // Calculate the number of parent configurations.
   int no_j = 1;
   for (j = 0; j < getIndividualSize (); j++)
      if (mParents [node][j]) no_j *= getState (j);

   // Allocate memory for all nijk-s.
   int** nijk = new int* [no_j];

   for (j = 0; j < no_j; j++) {

      nijk [j] = new int [getState (node)];

      for (k = 0; k < getState (node); k++)
         nijk [j][k] = 0;

   }

   // Calculate all nijk-s.
   for (j = 0; j < mSelectedSize; j++) {

      // Find the parent configuration for the j-th case.
      int parent_configuration = 0;

      for (int parent = 0;parent < getIndividualSize (); parent++)
         if (mParents [node][parent]) {
            parent_configuration *= getState (parent);
            parent_configuration += cases [j][parent];
         }

      // Update the corresponding nijk.
      nijk [parent_configuration][cases [j][node]]++;

   }

   // Calculate the BIC value.

   long double bic = 0;

   for (j = 0; j < no_j; j++) {

      int nij = 0;

      for (k = 0; k < getState (node); k++) {

         nij += nijk [j][k];

         // For rounding problems...
         if (nijk [j][k] != 0)
            bic += nijk [j][k] * log ((long double) (nijk [j][k]));

      }

      // For rounding problems...
      if (nij != 0)
         bic -= nij * log ((long double) nij);

   }

   bic -= log ((long double) mSelectedSize) * no_j / 2;

   // Free the memory allocated for the nijk-s.

   for (j = 0; j < no_j; j++)
      delete [] nijk [j];

   delete [] nijk;

   return bic;

}


// Function only to use with the Parallel version

long double GABayesianNetwork::deltaBIC (int node, int nparent, int**& cases) {

   int j, k;

   // Calculate the number of parent configurations.
   int no_j = 1;

   for (j = 0; j < getIndividualSize (); j++)
      if (j != nparent) {
         if (mParents [node][j])
            no_j *= getState (j);
      }
      else {
         if (!(mParents [node][nparent]))
            no_j *= getState (nparent);
      }

   // Allocate memory for all nijk-s.
   int** nijk = new int* [no_j];

   for (j = 0; j < no_j; j++) {

      nijk [j] = new int [getState (node)];

      for (k = 0; k < getState (node); k++)
         nijk [j][k] = 0;

   }

   // Calculate all nijk-s.
   for (j = 0; j < mSelectedSize; j++) {

      // Find the parent configuration for the j-th case.
      int parent_configuration = 0;

      for (int parent=0;parent<getIndividualSize();parent++)
         if (parent != nparent) {
            if (mParents [node][parent]) {

               parent_configuration *= getState (parent);
               parent_configuration += cases [j][parent];

            }
         }
         else {
            if (!(mParents [node][nparent])) {

               parent_configuration *= getState (nparent);
               parent_configuration += cases [j][nparent];

            }
         }

      // Update the corresponding nijk.
      nijk [parent_configuration][cases [j][node]]++;

   }

   // Calculate the BIC value.
   long double deltabic = 0;

   for (j = 0; j < no_j; j++) {

      int nij = 0;

      for (k = 0; k < getState (node); k++) {

         nij += nijk [j][k];

         // For rounding problems...
         if (nijk [j][k] != 0)
            deltabic += nijk [j][k] * log ((long double) (nijk [j][k]));

      }

      // For rounding problems...
      if (nij != 0)
         deltabic -= nij * log ((long double) nij);

   }

   deltabic -= log ((long double) mSelectedSize) * no_j / 2;


   // Free the memory allocated for the nijk-s.

   for (j = 0; j < no_j; j++)
      delete [] nijk [j];

   delete [] nijk;

   return deltabic;

}


void GABayesianNetwork::removeArc (int node, int parent) {

   // Remove the arc.

   mParents [node][parent] = false;

   // The paths of the Bayesian network are updated.

   for (int i = 0; i < getIndividualSize (); i++)
      for (int j = 0; j < getIndividualSize (); j++)
         if (mPaths [parent][j] > 0 && mPaths [i][node] > 0)
            mPaths [i][j] -= mPaths [parent][j] * mPaths [i][node];

}


long double GABayesianNetwork::K2Metric (int node, int**& cases) {

   int j;

   // Calculate the number of parent configurations.
   int no_j = 1;

   for (j = 0; j < getIndividualSize (); j++)
      if (mParents [node][j])
         no_j *= getState (j);

   // Allocate memory for all nijk-s.
   int** nijk = new int*[ no_j];

   for (j = 0; j < no_j; j++) {

      nijk [j] = new int [getState (node)];

      for (int k = 0; k < getState (node); k++)
         nijk [j][k] = 0;

   }

   // Calculate all nijk-s.
   for (j = 0; j < mSelectedSize; j++) {

      // Find the parent configuration for the j-th case.
      int parent_configuration = 0;

      for (int parent = 0; parent < getIndividualSize (); parent++)
         if (mParents [node][parent]) {

            parent_configuration *= getState (parent);
            parent_configuration += cases [j][parent];

         }

      // Update the corresponding nijk.
      nijk [parent_configuration][cases [j][node]]++;

   }

   // Calculate the K2 value.
   long double k2 = 0;

   for (j = 0; j < no_j; j++) {

      int nij = 0;

      for (int k = 0; k < getState (node); k++) {

         nij += nijk [j][k];
         k2 += logFact (nijk [j][k]);

      }

      k2 += logFact (getState (node) - 1);
      k2 -= logFact (nij + getState (node) - 1);

   }

   k2 -= log ((long double) mSelectedSize) * no_j / 2;

   // Free the memory allocated for the nijk-s.

   for (j = 0; j < no_j; j++)
      delete [] nijk [j];

   delete [] nijk;

   return k2;

}


long double GABayesianNetwork::deltaK2 (int node, int nparent, int**& cases) {

   int j;

   // Calculate the number of parent configurations.
   int no_j = 1;

   for (j = 0; j < getIndividualSize (); j++)
      if (j != nparent) {
         if (mParents [node][j])
            no_j *= getState (j);
      }
      else {
         if (!(mParents [node][nparent]))
            no_j *= getState (nparent);
      }

   // Allocate memory for all nijk-s.
   int** nijk = new int* [no_j];

   for (j = 0; j < no_j; j++) {

      nijk [j] = new int [getState (node)];

      for (int k = 0; k < getState (node); k++)
         nijk [j][k] = 0;

   }

   // Calculate all nijk-s.
   for (j = 0; j < mSelectedSize; j++) {

      // Find the parent configuration for the j-th case.
      int parent_configuration = 0;

      for (int parent=0;parent<getIndividualSize();parent++)
         if (parent != nparent) {
            if (mParents [node][parent]) {

               parent_configuration *= getState (parent);
               parent_configuration += cases [j][parent];

            }
         }
         else {
            if (!(mParents [node][nparent])) {

               parent_configuration *= getState (nparent);
               parent_configuration += cases [j][nparent];

            }
         }

      // Update the corresponding nijk.
      nijk [parent_configuration][cases [j][node]]++;

   }

   // Calculate the K2 value.
   long double k2 = 0;

   for (j = 0; j < no_j; j++) {

      int nij = 0;

      for (int k = 0; k < getState (node); k++) {

         nij += nijk [j][k];
         k2 += logFact (nijk [j][k]);

      }

      k2 += logFact (getState (node) - 1);
      k2 -= logFact (nij + getState (node) - 1);

   }

   k2 -= log ((long double) mSelectedSize) * no_j / 2;

   // Free the memory allocated for the nijk-s.

   for (j = 0; j < no_j; j++)
      delete [] nijk [j];

   delete [] nijk;

   return k2;

}


long double GABayesianNetwork::logFact (int n) {

   long double value = 0;

   for (int i = 2; i < n; i++)
      value += log ((long double) i);

   return value;

}
