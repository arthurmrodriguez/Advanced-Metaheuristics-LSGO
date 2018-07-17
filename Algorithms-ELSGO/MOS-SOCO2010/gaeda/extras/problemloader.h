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

#ifndef PROBLEMLOADER_H
#define PROBLEMLOADER_H

class GAGenome;
class GAPopulation;

typedef long double (*distToOptFuncPtr)      (GAGenome&);
typedef long double (*nComponentsOptFuncPtr) (GAGenome&);

/* DATA TYPES */
typedef struct {

   GAGenome* (*defineProblem) (void); // Pointer to the defineProblem function
   const char* (*describeProblem) (void);

   // Optional
   bool      (*postprocess)    (GAPopulation*, int);
   void      (*individualInit) (GAGenome&         );
   GAGenome* (*optimum)        (                  );
   int       optCriterion;

   distToOptFuncPtr      distToOpt;
   nComponentsOptFuncPtr nComponentsOpt;

} GAProblemStruct;

/* FUNCTION PROTOTYPES */
// Returns a pointer to the defineProblem function of the problem module
// located at libpath
GAProblemStruct* GAloadProblem (char* libpath);

#endif
