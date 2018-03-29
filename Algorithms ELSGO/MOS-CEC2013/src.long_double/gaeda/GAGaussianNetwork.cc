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

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <vector>

#include "GAGaussianNetwork.h"

#include "garandom.h"
#include "chi.h"
#include "normal.h"
#include "Matrix.h"
#include "gaparameters.h"
#include "std_stream.h"

GAGaussianNetwork::GAGaussianNetwork (const GAGenome &g, LearningMethod lm, ScoreMethod sm ) : GAGraphModel (g) {

   // Set the network parameters
   this->mScoringMethod  = sm;
   this->mLearningMethod = lm;

   // Init the network
   initialize();

}


// Get parameters of the gaussian network
int GAGaussianNetwork::get (const char* var, void* value) const {

   int status = 1;

   if (strcmp (var, gaNcontinuousLearningMethod ) == 0 ||
       strcmp (var, gaSNcontinuousLearningMethod) == 0    ) {

      LearningMethod* ptr = (LearningMethod*) value;
      *ptr = mLearningMethod;
      status = 0;

   }
   else if (strcmp (var, gaNcontinuousScoreMethod ) == 0 ||
            strcmp (var, gaSNcontinuousScoreMethod) == 0    ) {

      ScoreMethod* ptr = (ScoreMethod*) value;
      *ptr = mScoringMethod;
      status = 0;

   }

   return status;

}


// Set parameters of the gaussian network
int GAGaussianNetwork::setptr (const char* var, const void* value) {

   int status = 1;

   if (strcmp (var, gaNcontinuousLearningMethod ) ==0 ||
       strcmp (var, gaSNcontinuousLearningMethod) ==0    ) {

      LearningMethod* ptr = (LearningMethod*) value;
      mLearningMethod = *ptr;
      status = 0;

   }
   else if (strcmp (var, gaNcontinuousScoreMethod ) == 0 ||
            strcmp (var, gaSNcontinuousScoreMethod) == 0    ) {

      ScoreMethod* ptr = (ScoreMethod*) value;
      mScoringMethod = *ptr;
      status = 0;

   }

   return status;

}


void GAGaussianNetwork::initialize () {

   int i, j;

   mMin = new long double [getIndividualSize ()];
   mMax = new long double [getIndividualSize ()];

   for (int x = 0; x < getIndividualSize (); x++) {

      mMin [x] = getGenome ()->min (x);
      mMax [x] = getGenome ()->max (x);

   }

   mSelectionSize = 0;
   mNrPrecision   = 0.0;
   mNrMax         = 0;

   // The order is trivial.
   mTopSortedMap = new int [getIndividualSize ()];

   for (i = 0; i < getIndividualSize (); i++)
      mTopSortedMap [i] = i;

   // The ordered nodes can be found from mTopSortedMap.
   mTopSorted = new int [getIndividualSize ()];

   for (i = 0; i < getIndividualSize (); i++)
      mTopSorted [mTopSortedMap [i]] = i;

   // Memory for mA is allocated.
   mA = new long double* [getIndividualSize ()];

   for (i = 0; i < getIndividualSize (); i++) {

      mA [i] = new long double [getIndividualSize ()];

      for (j = 0; j < getIndividualSize (); j++)
         mA [i][j] = INT_MIN;

   }

   // The reference to the cases' database is initialised
   m_cases = NULL;

   // Finally, the numerical parameters are created. All variables will be N(0,1).
   mV = new long double [getIndividualSize ()];

   for (i = 0; i < getIndividualSize (); i++)
      mV [i] = (mMax [i] - mMin [i]) / 2;


   mMeans = new long double [getIndividualSize ()];

   for (i = 0; i < getIndividualSize (); i++)
      mMeans [i] = mMin [i] + (mMax [i] - mMin [i]) / 2;


   mB = new long double* [getIndividualSize ()];

   for (i = 0; i < getIndividualSize (); i++) {

      mB [i] = new long double [getIndividualSize ()];

      for (j = 0; j < getIndividualSize (); j++)
         mB [i][j] = 0.0;

   }


   mSigma = new long double* [getIndividualSize ()];

   for (i = 0; i < getIndividualSize (); i++) {

      mSigma [i] = new long double [getIndividualSize ()];

      for (j = 0; j < getIndividualSize (); j++)
         if (i == j)
            mSigma [i][j] = mV [i];
         else
            mSigma [i][j] = 0.0;

   }


   mT0 = rejectionRegionBoundary ();


   // If score is BIC in the EGNA and EE cases, the structures have to be defined

   if (mScoringMethod == BIC_SCORE) {

      m_Sxixj = new long double* [getIndividualSize ()];

      for (i = 0; i < getIndividualSize (); i++) {

         m_Sxixj [i] = new long double [getIndividualSize ()];

         for (j = 0; j < getIndividualSize (); j++)
            m_Sxixj [i][j] = 0.0;

      }

      m_S2xi = new long double [getIndividualSize ()];

      for (i = 0; i < getIndividualSize (); i++)
         m_S2xi [i] = 0.0;

   }

}


GAGaussianNetwork::~GAGaussianNetwork () {

   int i;

   delete [] mTopSortedMap;
   delete [] mTopSorted;

   for (i = 0; i < getIndividualSize (); i++)
      delete [] mParents [i];

   delete [] mParents;

   for (i = 0; i < getIndividualSize (); i++)
      delete [] mA [i];

   delete [] mA;

   for (i = 0; i < getIndividualSize (); i++)
      delete [] mPaths [i];

   delete [] mPaths;

   delete [] mMeans;
   delete [] mV;

   for (i = 0; i < getIndividualSize (); i++)
      delete [] mB [i];

   delete [] mB;

   if (mScoringMethod == BIC_SCORE) {

      for (i = 0; i < getIndividualSize (); i++)
         delete [] m_Sxixj [i];

      delete [] m_Sxixj;
      delete [] m_S2xi;

   }

}

void GAGaussianNetwork::simulate (GAGenome& genes) {

  long double normal, mean, variance;

  // The individual will be generated simulating its genes according
  // to the order they have in the gaussian network.

  for (int order = 0; order < getIndividualSize (); order++) {

    distribution(order,genes,mean,variance);

      do
        normal = normal_random (mean, variance);
      while((normal < mMin [mTopSorted [order]]) || (normal>mMax [mTopSorted [order]]));
      /*   This bug was detected while testing with a deceptive function. It does not simulate correctly a truncated normal but
       *   the results are considerably better than the ones obtained with the previous function
      if (normal < mMin [mTopSorted [order]])
         genes.setValueOfContinuousGene (mTopSorted [order], mMin [mTopSorted [order]]);
      else if (normal>mMax [mTopSorted [order]])
         genes.setValueOfContinuousGene (mTopSorted [order], mMax [mTopSorted [order]]);
         else*/
      genes.setValueOfContinuousGene (mTopSorted [order], normal);

  }

}

void GAGaussianNetwork::learn (const GAPopulation& cases, long double sel_pctg) {

   // Learn the structure
   int I;

   I = cases.size ();
   mSelectionSize = (int) (I * sel_pctg);

   // Copy the cases' database (required for the BIC score)
   m_cases = &cases;

   // The given cases are transformed to the means and covariances. All
   // calculations will be made using those values.
   calculateMeansAndCovariances (cases);

   // If the score is BIC, we also need to calculate the Sxixj and S2xi values
   if (mScoringMethod == BIC_SCORE)
      calculateStructuresforBIC (cases);

   switch (mLearningMethod) {

      case UMDA:

         // No learning required.
         break;

      case EGNA_B:

         learnEGNAB ();
         break;

      case EGNA_LOCAL:

         learnEGNALocal ();
         break;

      case MIMIC:

         learnMIMIC ();
         break;

      case EE:

         learnEE ();
         break;

      case EMNA:

         // No learning required, the structure is
         // complete as initialised at the beginning
         break;

   }

   // Learn the probabilities
   if (mScoringMethod != BIC_SCORE)
      learnDistributions ();
   else
      distributionParametersBIC ();

}


void GAGaussianNetwork::distribution (int ordered_node, GAGenome& instance,
                                      long double& mean, long double& variance) {

   int j;
   long double k1, k2;

   mean = mMeans [mTopSorted [ordered_node]];

   for (int order = 0; order < ordered_node; order++)
      if(mParents [mTopSorted [ordered_node]][mTopSorted [order]]) {

         k1 = mean;
         k2 = mB [mTopSorted [ordered_node]][mTopSorted [order]]      *
              (instance.getValueOfContinuousGene (mTopSorted [order]) -
              mMeans [mTopSorted [order]]);

         mean += mB [mTopSorted [ordered_node]][mTopSorted [order]]      *
                 (instance.getValueOfContinuousGene (mTopSorted [order]) -
                 mMeans [mTopSorted [order]]);

      }

   variance = mV [mTopSorted [ordered_node]];

   if (variance < 0) {

      STD_CERR << "negative variance " << variance << " !!!!!" << STD_ENDL;
      STD_CERR << "index: " << mTopSorted [ordered_node] << " " << ordered_node << STD_ENDL;
      STD_CERR << "mV: { ";

      for (j = 0; j < getIndividualSize () - 1; j++)
         STD_CERR << mV[j] << ", ";

      STD_CERR << mV [getIndividualSize () - 1] << " }" << STD_ENDL;
      variance = 0.0;

   }

}


void GAGaussianNetwork::learnDistributions () {

   int i,j;
   // m values have already been calculated by LearnMeansAndCovariances.
   // b and v values are calculated according to DeGroot.

   for (i = 0; i < getIndividualSize (); i++) {

      Matrix sigma (mSigma, getIndividualSize (), getIndividualSize ());

      CSet x_1;
      x_1.insert (i);

      CSet x_2;
      for (j = 0; j < getIndividualSize (); j++)
         if (mParents [i][j])
            x_2.insert (j);

      Matrix sigma_11 = sigma.Partition (x_1, x_1);
      Matrix sigma_12 = sigma.Partition (x_1, x_2);
      Matrix sigma_21 = sigma.Partition (x_2, x_1);
      Matrix sigma_22 = sigma.Partition (x_2, x_2);

      Matrix i_sigma_22 = sigma_22.Invert ();
      Matrix b = sigma_12 * i_sigma_22;

      int idx = 0;

      for (j = 0; j < getIndividualSize (); j++)
         if (mParents [i][j]) {
            mB [i][j] = b.Item (0, idx);
            idx++;
         }
         else
            mB [i][j] = 0;


      if (x_2.size () == 0)
         mV [i] = sigma_11.Item (0, 0);
      else {

         Matrix i_t_11 = sigma_11 - sigma_12 * i_sigma_22 * sigma_21;
         mV [i] = i_t_11.Item (0, 0);

         // Comprobando varianza!
         if (mV [i] < 0)
            STD_CERR << "Errorea mV!!! " << STD_ENDL;

      }

   }

}


void GAGaussianNetwork::learnEGNAB () {

   initGaussianNetwork ();
   learnEGNALocal ();

}

void GAGaussianNetwork::learnEGNALocal () {

   int i, j;
   long double max;

   calculateA ();

   do {

      // Find the arc which modification most increases
      // the metric of the Gaussian network.

      int max_i = 0;
      int max_j = 0;
      max = INT_MIN;

      for (i = 0; i < getIndividualSize (); i++)
         for (j = 0; j < getIndividualSize (); j++)
            if (mPaths [j][i] == 0 && mA [i][j] > max) {
               // In order to modify the arc j->i there
               // cannot be a path from i to j.

               max_i = i;
               max_j = j;
               max = mA [i][j];

            }

      // If the metric of the Gaussian network can
      // be improved the arc is modified.

      if (max > 0) {

         if (!mParents [max_i][max_j]) // If there is no arc j->i, it is added.
            addArc (max_i, max_j);
         else // If there is an arc j->i, it is removed.
            removeArc (max_i, max_j);

         calculateANode(max_i);
      }

   } while (max > 0);

}


void GAGaussianNetwork::initGaussianNetwork (bool complete) {

   int i, j;

   for (i = 0; i < getIndividualSize (); i++) {

      mTopSortedMap [i] = i;
      mTopSorted    [i] = i;

      for (j = 0; j < getIndividualSize (); j++) {

         mParents [i][j] = false;

         if (i == j)
            mPaths [i][j] = 1;
         else
            mPaths [i][j] = 0;

      }

   }

   if (complete)
      for (i = 0; i < getIndividualSize (); i++)
         for (j = i + 1; j < getIndividualSize (); j++)
            addArc (j, i);

}


void GAGaussianNetwork::calculateA () {

   int i;

   for (i = 0; i < getIndividualSize (); i++)
      calculateANode (i);

}


void GAGaussianNetwork::calculateANode (int node) {

   int i;
   long double old_gh, new_gh;

   if (mScoringMethod == BGe_SCORE)
      old_gh = BGe (node);
   else
      old_gh = BIC (node);

   for (i = 0; i < getIndividualSize (); i++)
      if ((i != node) && (mPaths [i][node] == 0)) {

         mParents [node][i] = !mParents [node][i];

         if (mScoringMethod == BGe_SCORE)
            new_gh = BGe (node);
         else
            new_gh = BIC (node);

         mParents [node][i] = !mParents [node][i];
         mA [node][i] = new_gh - old_gh;

      }
      else
         mA [node][i] = INT_MIN;

}


void GAGaussianNetwork::removeArc (int node, int parent) {

   int i,j;

   // removeArc the arc.
   mParents [node][parent] = false;

   // The paths of the Gaussian network are updated.
   for (i = 0; i < getIndividualSize (); i++)
      for (j = 0; j < getIndividualSize (); j++)
         if (mPaths [parent][j] > 0 && mPaths [i][node] > 0)
            mPaths [i][j]-= mPaths [parent][j] * mPaths [i][node];

}


void GAGaussianNetwork::addArc (int node, int parent) {

   int i, j, k;

   // addArc the arc.
   mParents [node][parent] = true;

   // The paths of the Gaussian network are updated.
   for (i = 0; i < getIndividualSize (); i++)
      for (j = 0; j < getIndividualSize (); j++)
         if (mPaths [parent][j] > 0 && mPaths [i][node] > 0)
            mPaths [i][j] += mPaths [parent][j] * mPaths [i][node];


   // Update the ordering of the nodes.
   if (mTopSortedMap [node] < mTopSortedMap [parent]) {

      // The order of node, its descendants, parent
      // and its ancestors must be updated.

      // How many positions the ancestors of parent are moved.
      int jump_parent = 0;

      // Calculate how many positions the descendants of node must be moved.
      int jump_node = 0;

      for (k = mTopSortedMap [node]; k <= mTopSortedMap [parent]; k++)
         if (mPaths [parent][mTopSorted [k]] > 0)
            jump_node++;

      // Update the order of the nodes between
      // node and parent (both included).

      for (k = mTopSortedMap [node]; k <= mTopSortedMap [parent]; k++)
         if (mPaths [parent][mTopSorted [k]] > 0) {

            mTopSortedMap [mTopSorted [k]] += jump_parent;
            jump_node--;

         }
         else {

            mTopSortedMap [mTopSorted [k]] += jump_node;
            jump_parent--;

         }

      // Update the ordered nodes.

      for (k = 0; k < getIndividualSize (); k++)
         mTopSorted [mTopSortedMap [k]] = k;

   }

}


void GAGaussianNetwork::learnMIMIC () {

   int i, j;

   initGaussianNetwork ();

   bool* linked = new bool [getIndividualSize ()];

   for (i = 0; i < getIndividualSize (); i++)
      linked [i] = false;

   int min_i = -1;
   long double min_entropy_i = INT_MAX;

   for (i = 0; i < getIndividualSize (); i++) {

      long double entropy_i = entropy (i);

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

      for (j = 0; j < getIndividualSize (); j++)

         if (!linked [j]) {

            long double entropy_j = entropy (j, head);

            if(entropy_j<min_entropy_j) {

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


long double GAGaussianNetwork::entropy (int node) {

   return 0.5 * (1 + log (2 * M_PI)) + log (sqrt (mSigma [node][node]));

}


long double GAGaussianNetwork::entropy (int node, int parent) {

   return 0.5 * (1 + log (2 * M_PI)) +
          0.5 * log (mSigma [node][node] - mSigma [node][parent] / mSigma [parent][parent]);

}


void GAGaussianNetwork::calculateMeansAndCovariances (const GAPopulation& cases) {

   int i, j;

   for (i = 0; i < getIndividualSize (); i++) {

      mMeans [i] = 0;

      for (j = 0; j < getIndividualSize (); j++)
         mSigma [i][j] = 0;

   }

   for (i = 0; i < mSelectionSize; i++)
      for (j = 0; j < getIndividualSize (); j++)
         mMeans [j] += cases.individual (i).getValueOfContinuousGene (j);

   for (i = 0; i < getIndividualSize (); i++)
      mMeans [i] /= mSelectionSize;

   for (i = 0; i < mSelectionSize; i++)
      for (j = 0; j < getIndividualSize(); j++)
         for (int k = 0; k < getIndividualSize (); k++)
            mSigma [j][k] += (cases.individual (i).getValueOfContinuousGene (j) - mMeans [j]) *
                             (cases.individual (i).getValueOfContinuousGene (k) - mMeans [k]);

   for (i = 0; i < getIndividualSize (); i++)
      for (j = 0; j < getIndividualSize (); j++)
         mSigma [i][j] /= mSelectionSize;

}


long double GAGaussianNetwork::BGe (int node) {

   int i;
   CSet Y;
   long double gh_parents;

   for (i = 0; i < getIndividualSize (); i++)
      if (mParents [node][i])
         Y.insert (i);

   if (Y.size () == 0)
      gh_parents = 0;
   else
      gh_parents = BGe (Y);

   Y.insert (node);

   long double gh_node_parents = BGe (Y);

   return gh_node_parents - gh_parents;

}


long double GAGaussianNetwork::BGe (CSet& nodes) {

   int i;
   int m = mSelectionSize;
   int n = getIndividualSize ();

   int alpha_W = n + 2;
   int alpha_m = alpha_W;

   Matrix m_0 (getIndividualSize (), 1);

   for (i = 0; i < getIndividualSize (); i++)
      // Change introduced the 23-10-2000 to improve initialization on the
      // learning on EGNA algorithms
      m_0.Item (i, 0) = mMeans [i];


   Matrix W_0 (getIndividualSize (),getIndividualSize ());

   for (i = 0; i < getIndividualSize (); i++)
      W_0.Item (i, i) = 1;


   Matrix x_m (getIndividualSize (), 1);

   for (i = 0; i < getIndividualSize (); i++)
      x_m.Item (i, 0) = mMeans [i];


   Matrix xY_m = x_m.RowPartition (nodes);

   Matrix sigma (mSigma, getIndividualSize (), getIndividualSize ());
   Matrix SY_m = sigma.Partition (nodes, nodes) * mSelectionSize;

   int l = nodes.size ();

   Matrix mY_0 = m_0.RowPartition (nodes);
   Matrix WY_0 = W_0.Partition (nodes, nodes);

   int alphaY_W = alpha_W - n + l;

   Matrix dif = mY_0 - xY_m;
   Matrix t_dif = dif.Transpose ();

   Matrix WY_m = WY_0 + SY_m + (dif * t_dif) * (alpha_m * m / (alpha_m + m));

   long double detWY_0 = WY_0.Determinant ();
   long double detWY_m = WY_m.Determinant ();

   long double log_gh = 0;

   log_gh += -(long double) l * m / 2 * log (M_PI);
   log_gh +=  (long double) l / 2 * log ((long double) alpha_m / (alpha_m + m));

   log_gh += logC (l, alphaY_W + m);
   log_gh -= logC (l, alphaY_W);

   log_gh +=  (long double) alphaY_W / 2 * log (detWY_0);
   log_gh += -(long double) (alphaY_W + m) / 2 * log (detWY_m);

   return log_gh;

}


long double GAGaussianNetwork::logC (int l, int alpha) {

   int i;
   long double value = 0;

   for (i = 1; i <= l; i++)
      value += logGamma ((long double) (alpha + 1 - i) / 2);

   return value;

}


long double GAGaussianNetwork::logGamma (long double num) {

   long double value = 0;
   long double p = num;

   while (p > 1) {

      value += log (p - 1);
      p--;

   }

   if (p == 0.5)
      value += log (sqrt (M_PI));

   return value;

}


void GAGaussianNetwork::learnEE () {

   int i, j;
   initGaussianNetwork (true);

   Matrix sigma (mSigma, getIndividualSize (), getIndividualSize ());
   Matrix i_sigma = sigma.Invert ();

   for (i = 0; i < getIndividualSize (); i++)
      for (j = 0; j < getIndividualSize (); j++)
         if(mParents [i][j]) {

            // Calculate the likelihood ratio test statistic

            long double r_ij_rest = -i_sigma.Item (i, j) /
                                    sqrt (i_sigma.Item (i, i) * i_sigma.Item (j, j));

            long double T_i = -mSelectionSize * log (1 - r_ij_rest * r_ij_rest);

            if (T_i <= mT0)
               removeArc (i, j);

         }

}


long double GAGaussianNetwork::rejectionRegionBoundary () {

   int i;
   long double x = rand () % 5 + 1;

   for (i = 0; i < mNrMax && fabs (F (x)) >= mNrPrecision; i++) {

      x -= F (x) / f (x);

      if (x < 0)
         x = 1; // The rejection region boundary must be positive.

   }

   return x;

}


long double GAGaussianNetwork::F (long double x) {

   long double b = -0.5 * (2 * getIndividualSize () + 1) * x;

   return F_chisq (x) + b * f_chisq (x) / mSelectionSize - (1 - ALPHA);

}


long double GAGaussianNetwork::f (long double x) {

   long double a = 0.25 * (x - 1) * (2 * getIndividualSize () + 1);

   return f_chisq (x) + a * f_chisq (x) / mSelectionSize;

}


/****************************************************************
 ****************************************************************
 ****************************************************************
 FUNCTIONS FOR BIC SCORE
 ****************************************************************
 ****************************************************************
 ****************************************************************/

long double GAGaussianNetwork::BIC (int node) {

   int i;
   CSet Y;

   for (i = 0; i < getIndividualSize (); i++)
      if (mParents [node][i])
         Y.insert (i);

   long double gh_parents;

   if (Y.size () == 0)
      gh_parents = 0;
   else
      gh_parents = BIC (node, Y);

   Y.insert (node);

   long double gh_node_parents = BIC (node, Y);

   return gh_node_parents - gh_parents;

}


// This function also requires the use of the database (copied in m_cases)
long double GAGaussianNetwork::BIC (int node, CSet& ParentNodes) {

   int r, j;
   long double berredura;

   // Compute v [node]
   mV [node] = m_S2xi [node];

   for (j = 0; j < getIndividualSize (); j++)
      if (mParents[node][j]) {

         mV [node] -= m_Sxixj [j][node] / m_S2xi [j];

         for (int k = j + 1; k < getIndividualSize (); k++)
            if (mParents [node][k])
               mV [node] += (m_Sxixj [j][k] * m_Sxixj [j][node] * m_Sxixj [k][node]) / (m_S2xi [j] * m_S2xi [k]);

      }

   // Compute initial log
   long double log2pi = log ((long double) sqrt (2 * M_PI) * mV [node]);

   // Compute penalization
   long double penal = 0.5 * log ((long double) mSelectionSize) * (/*2*getIndividualSize()+*/ ParentNodes.size ());

   // Compute big parenthesis
   long double body = 0;

   for (r = 0; r < mSelectionSize; r++) {

      berredura = 0;

      for (j = 0; j < getIndividualSize (); j++)
         if (ParentNodes.count (j))
            berredura +=  mB [node][j] * (m_cases->individual (r).getValueOfContinuousGene (j) - mMeans [j]);

      berredura += m_cases->individual (r).getValueOfContinuousGene (node) - mMeans [node];
      body -= berredura * berredura * 0.5;

   }

   body -= log2pi * mSelectionSize;

   return (body - penal);

}


void GAGaussianNetwork::calculateStructuresforBIC (const GAPopulation& cases) {

   // All the mMeans[i] have already been calculated

   for (int i = 0; i < getIndividualSize (); i++) {

      // Sxixj
      for (int j = 0; j < getIndividualSize (); j++) {

         m_Sxixj [i][j] = 0.0;

         for (int r = 0; r < mSelectionSize; r++)
            m_Sxixj [i][j] += (cases.individual (r).getValueOfContinuousGene (j) - mMeans [j]) *
                              (cases.individual (r).getValueOfContinuousGene (i) - mMeans [i]) / (mSelectionSize);

      }

      // S2xi
      m_S2xi[i] = 0.0;

      for (int r = 0; r < mSelectionSize; r++)
         m_S2xi [i] += (cases.individual (r).getValueOfContinuousGene (i) - mMeans [i]) *
                       (cases.individual (r).getValueOfContinuousGene (i) - mMeans [i]) / (mSelectionSize);

   }

}


void GAGaussianNetwork::distributionParametersBIC () {

   int i;

   for (i=0; i<getIndividualSize(); i++)
      for (int j = 0; j < getIndividualSize (); j++) // b_ji
         mB [j][i] = m_Sxixj [j][i] / m_S2xi [j];

   // The v_i are computed here
   for (i = 0; i < getIndividualSize (); i++) {

      mV [i] = m_S2xi [i];

      for (int j=0; j<getIndividualSize(); j++)
         if (mParents [i][j]) {

            mV [i] -= m_Sxixj [j][i] / m_S2xi [j];

            for (int k = j + 1; k < getIndividualSize (); k++)
               if (mParents [i][k])
                  mV [i] += (m_Sxixj [j][k] * m_Sxixj [j][i] * m_Sxixj [k][i]) / (m_S2xi [j] * m_S2xi [k]);

         }

   }

}


/****************************************************************
 ****************************************************************
 ****************************************************************
 FUNCTIONS FOR EMNA ALGORITHMS
 ****************************************************************
 ****************************************************************
 ****************************************************************/

void GAGaussianNetwork::initEMNAGaussianNetwork (const GAPopulation& cases) {

   // Initialise the Gaussian Network. This will always remain == for EMNAs
   initGaussianNetwork (true);

   // Calculate the Means vector and Variance-Covariance matrix
   calculateMeansAndCovariances (cases);

   // Calculate mV[] and mB[]
   learnDistributions ();

}


void GAGaussianNetwork::EMNAAdapt (const GAPopulation& cases, const GAGenome& genes_new,
                                   const GAGenome& genes_worse) {

   int i, j;
   std::vector <long double> Cov_tmp (getIndividualSize ());

   // Pre-initializations - before changing mMeans[]
   Cov_tmp.resize (getIndividualSize ());

   for (i = 0; i < getIndividualSize (); i++) {

      Cov_tmp [i] = 0.0;

      for (int r = 0; r < mSelectionSize; r++)
         Cov_tmp [i] += (cases.individual (r).getValueOfContinuousGene (i) - mMeans [i]);

   }

   // Means
   for (i = 0; i < getIndividualSize (); i++)
      mMeans [i] += (genes_new.getValueOfContinuousGene (i) -
                     genes_worse.getValueOfContinuousGene (i)) / mSelectionSize;

   // Variances & Covariances
   for (i = 0; i < getIndividualSize (); i++)
      for (j = 0; j < getIndividualSize (); j++)
         mSigma [i][j] -= (genes_new.getValueOfContinuousGene (i) - genes_worse.getValueOfContinuousGene (i)) *
                          Cov_tmp [j] / (mSelectionSize * mSelectionSize) -
                          (genes_new.getValueOfContinuousGene (j) - genes_worse.getValueOfContinuousGene (j)) *
                          Cov_tmp [i] / (mSelectionSize * mSelectionSize) +
                          (genes_new.getValueOfContinuousGene (i) - genes_worse.getValueOfContinuousGene (i)) *
                          (genes_new.getValueOfContinuousGene (j) - genes_worse.getValueOfContinuousGene (j)) /
                          (mSelectionSize * mSelectionSize) - (genes_worse.getValueOfContinuousGene (i) - mMeans [i]) *
                          (genes_worse.getValueOfContinuousGene (j) - mMeans [j]) /
                          mSelectionSize + (genes_new.getValueOfContinuousGene (i) - mMeans [i]) *
                          (genes_new.getValueOfContinuousGene (j) - mMeans [j]) / mSelectionSize;

   // Calculate mV[] and mB[]
   learnDistributions ();

}
