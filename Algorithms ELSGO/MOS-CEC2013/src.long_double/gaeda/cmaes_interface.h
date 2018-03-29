/* --------------------------------------------------------- */
/* --- File: cmaes_interface.h - Author: Nikolaus Hansen --- */
/* ---------------------- last modified:  IV 2007        --- */
/* --------------------------------- by: Nikolaus Hansen --- */
/* --------------------------------------------------------- */
/*   
     CMA-ES for non-linear function minimization. 

     Copyright (C) 1996, 2003, 2007 Nikolaus Hansen. 
     e-mail: hansen AT bionik.tu-berlin.de
             hansen AT lri.fr

     License: see file cmaes.c
*/
#ifndef CMAES_INTER_
#define CMAES_INTER_

#ifdef __cplusplus
extern "C" {
#endif

#include "cmaes.h"

/* --------------------------------------------------------- */
/* ------------------ Interface ---------------------------- */
/* --------------------------------------------------------- */

/* --- initialization, constructors, destructors --- */
long double * cmaes_init(cmaes_t *, int dimension , long double *xstart, 
		long double *stddev, long seed, int lambda, 
		const char *input_parameter_filename);
void cmaes_resume_distribution(cmaes_t *evo_ptr, char *filename);
void cmaes_exit(cmaes_t *);

/* --- core functions --- */
long double * const * cmaes_SamplePopulation(cmaes_t *);
long double *         cmaes_UpdateDistribution(cmaes_t *, 
					  const long double *rgFitnessValues);
const char *     cmaes_TestForTermination(cmaes_t *);

/* --- additional functions --- */
long double * const * cmaes_ReSampleSingle( cmaes_t *t, int index);
long double const *   cmaes_ReSampleSingle_old(cmaes_t *, long double *rgx); 
long double *         cmaes_SampleSingleInto( cmaes_t *t, long double *rgx);
void             cmaes_UpdateEigensystem(cmaes_t *, int flgforce);

/* --- getter functions --- */
long double         cmaes_Get(cmaes_t *, char const *keyword);
const long double * cmaes_GetPtr(cmaes_t *, char const *keyword); /* e.g. "xbestever" */
long double *       cmaes_GetNew( cmaes_t *t, char const *keyword); 
long double *       cmaes_GetInto( cmaes_t *t, char const *keyword, long double *mem); 

/* --- online control and output --- */
void           cmaes_ReadSignals(cmaes_t *, char const *filename);
void           cmaes_WriteToFile(cmaes_t *, const char *szKeyWord,
                                 const char *output_filename); 
char *         cmaes_SayHello(cmaes_t *);
/* --- misc --- */
long double *       cmaes_NewDouble(int n);  
void           cmaes_FATAL(char const *s1, char const *s2, char const *s3, 
			   char const *s4);

void cmaes_output_dir(cmaes_t * t, const char * dir);
long double const * cmaes_SetMean(cmaes_t *, const long double *xmean);

#ifdef __cplusplus
}
#endif

#endif
