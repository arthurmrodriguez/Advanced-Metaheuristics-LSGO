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
#include "centroids.h"

#include "../genomes/GAGenome.h"

/* FUNCTIONS */

// This function tries to generate a given number of centroids. These centroids are
// just genomes that are different enough from each other. The difference between
// genomes is measured by a distance function. The minimum distance between
// genomes is provided
//
// What is this useful for?. For island models, each centroid can be used to
// generate inidivuals for an island model, thus assuring that not two islands
// are working with similar individuals
//
// Arguments:
//
// number: Number of centroids
// minDist: Minimum distance between any two centroids
// sample: A sample genome. The initial centroids will be all a clone of this, and then
//         they will be mutated
// maxiters: Maximum number of iterations
// *success: Set to true if a solution is found before maxiters iterations. 0 Otherwise
// *minfound: Minimum distance between any two centroid
//
//

GACentroidVector* GAgenerateCentroids (int number, double minDist, GAGenome& sample, int maxiters, bool* success, double* minfound) {

   GACentroidVector* centroids = new GACentroidVector;
   double dist;
   int it = 0;
   bool done = false;

   *success = true;

   // Create an initial genomes vector
   for (int i = 0; i < number; i++)
      centroids->insert (centroids->end (), sample.clone ());

   // We now cycle until each genome is distant appart from the
   // others. Since we don't have specific knowledge about the
   // genomes or the distance function, this is just a trial and
   // error process.
   while (!done) {

      done = true;

      for (int i = 0; i < number; i++) {

         for (int j = i + 1; j < number; j++) {

            dist = (*centroids) [i]->compare (*(*centroids)[j]);

            if (dist < minDist) {

               (*centroids) [i]->mutate ((minDist - dist) / minDist);
               (*centroids) [j]->mutate ((minDist - dist) / minDist);

               done = false;

            }

         }

      }

      it++;

      if (it > maxiters) {

         *success = false;
          done    = true;
      }

   }

   // Calculate minimum distance between any two centroids
   *minfound = (*centroids) [0]->compare (*(*centroids) [1]);

   for (int i = 0; i < number; i++) {

      for (int j = i + 1; j < number; j++) {

         dist = (*centroids) [i]->compare (*(*centroids) [j]);

         if (dist < *minfound)
            *minfound = dist;

      }

   }

   // Just in case the mutations on last iteration find a correct solution
   if (*minfound >= minDist)
      *success = true;

   return centroids;

}
