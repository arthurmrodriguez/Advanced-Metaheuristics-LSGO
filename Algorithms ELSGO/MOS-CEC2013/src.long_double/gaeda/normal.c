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
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>

#include "normal.h"

/* Simulates a N(mean,variance) random variable */
long double normal_random (long double mean, long double variance) {

   int i;
   long double x = 0;
   long double normal_01;
   long double normal;

   // Sanity check
   assert (variance >= 0);

   // First, we generate a N(0,1) random number based on the Central

   // Average 12 U(0,1) variables.
   for (i = 0; i < 12; i++)
      x += rand ();

   x /= RAND_MAX;

   // Normalize.
   normal_01 = x - 6;
   normal = mean + normal_01 * sqrt (variance);

   return normal;

}
