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
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "problemloader.h"

#include "../std_stream.h"
#include "../genomes/GAGenome.h"

GAProblemStruct* GAloadProblem (char* libpath) {

   void* handle;
   handle = dlopen (libpath, RTLD_NOW | RTLD_GLOBAL);

   GAProblemStruct* res = (GAProblemStruct*) malloc ((sizeof (GAProblemStruct)));

   if (handle == NULL){
      free (res);
      STD_CERR << "Problem handling the file." << STD_ENDL;
      STD_CERR << dlerror () << STD_ENDL;
      return NULL;
   }

   res->defineProblem = (GAGenome* (*) (void)) dlsym (handle, "defineProblem");

   if (res->defineProblem == NULL)  {
      free (res);
      STD_CERR << "Problem with function: 'defineProblem' or 'defineMOSproblem' (none provided)." << STD_ENDL;
      return NULL;
   }

   res->describeProblem = (const char* (*) (void)) dlsym (handle, "describeProblem");

   if (res->describeProblem == NULL){
      free (res);
      STD_CERR << "Problem with function: 'describeProblem'." << STD_ENDL;
      return NULL;
   }

   // These are optional
   res->individualInit = (void (*) (GAGenome&))     dlsym (handle, "individualInit");

   res->postprocess = (bool (*) (GAPopulation*,int)) dlsym (handle, "postprocess");

   // Chapuzillas rapidas que habra que cambiar en un futuro y no repetir codigo, por ej
   // con una macro
   void* sym  =  dlsym (handle, "getOptimumGenome");
   if (!sym) res->optimum = NULL;
   else      res->optimum =  (GAGenome* (*) () )sym;

   sym  =  dlsym (handle, "disToOpt");
   if (!sym) res->distToOpt = NULL;
   else      res->distToOpt = (distToOptFuncPtr )sym;

   sym  =  dlsym (handle, "nComponentsOpt");
   if (!sym) res->nComponentsOpt = NULL;
   else      res->nComponentsOpt = (nComponentsOptFuncPtr)sym;

   sym  =  dlsym (handle, "optCriterion");
   if (!sym) res->optCriterion = -1;
   else      res->optCriterion = ((GAGenome::OptCriterion (*) ())sym) ();

   return res;

}
