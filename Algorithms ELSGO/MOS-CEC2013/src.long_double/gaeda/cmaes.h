/* --------------------------------------------------------- */
/* --- File: cmaes.h ----------- Author: Nikolaus Hansen --- */
/* ---------------------- last modified: VIII 2007       --- */
/* --------------------------------- by: Nikolaus Hansen --- */
/* --------------------------------------------------------- */
/*   
     CMA-ES for non-linear function minimization. 

     Copyright (C) 1996, 2003-2007  Nikolaus Hansen. 
     e-mail: hansen@bionik.tu-berlin.de
      
     License: see file cmaes.c
   
*/
#ifndef NH_cmaes_h /* only include ones */ 
#define NH_cmaes_h 

#ifdef __cplusplus
extern "C" {
#endif

#include <time.h>

typedef struct 
/* random_t 
 * sets up a pseudo random number generator instance 
 */
{
  /* Variables for Uniform() */
  long int startseed;
  long int aktseed;
  long int aktrand;
  long int *rgrand;
  
  /* Variables for Gauss() */
  short flgstored;
  long double hold;
} random_t;

typedef struct 
/* timings_t 
 * time measurement, used to time eigendecomposition 
 */
{
  /* for outside use */
  long double totaltime; 
  long double tictoctime; 
  long double lasttictoctime;
  
  /* local fields */
  clock_t lastclock;
  time_t lasttime;
  clock_t ticclock;
  time_t tictime;
  short istic;
  short isstarted; 

  long double lastdiff;
  long double tictoczwischensumme;
} timings_t;

typedef struct 
/* readpara_t
 * collects all parameters, in particular those that are read from 
 * a file before to start. This should split in future? 
 */
{
  /* input parameter */
  int N; /* problem dimension, must stay constant */
  unsigned int seed; 
  long double * xstart; 
  long double * typicalX; 
  int typicalXcase;
  long double * rgInitialStds;
  long double * rgDiffMinChange; 

  /* termination parameters */
  long double stopMaxFunEvals; 
  long double facmaxeval;
  long double stopMaxIter; 
  struct { int flg; long double val; } stStopFitness; 
  long double stopTolFun;
  long double stopTolFunHist;
  long double stopTolX;
  long double stopTolUpXFactor;

  /* internal evolution strategy parameters */
  int lambda;          /* -> mu, <- N */
  int mu;              /* -> weights, (lambda) */
  long double mucov, mueff; /* <- weights */
  long double *weights;     /* <- mu, -> mueff, mucov, ccov */
  long double damps;        /* <- cs, maxeval, lambda */
  long double cs;           /* -> damps, <- N */
  long double ccumcov;      /* <- N */
  long double ccov;         /* <- mucov, <- N */
  struct { int flgalways; long double modulo; long double maxtime; } updateCmode;
  long double facupdateCmode;

  /* supplementary variables */

  char *weigkey; 
  char resumefile[99];
  char **rgsformat;
  void **rgpadr;
  char **rgskeyar;
  long double ***rgp2adr;
  int n1para, n1outpara;
  int n2para;
} readpara_t;

typedef struct 
/* cmaes_t 
 * CMA-ES "object" 
 */
{
  char *version;
  readpara_t sp;
  random_t rand; /* random number generator */

  long double sigma;  /* step size */

  long double *rgxmean;  /* mean x vector, "parent" */
  long double *rgxbestever; 
  long double **rgrgx;   /* range of x-vectors, lambda offspring */
  int *index;       /* sorting index of sample pop. */
  long double *arFuncValueHist;

  short flgIniphase; /* not really in use anymore */
  short flgStop; 

  long double chiN; 
  long double **C;  /* lower triangular matrix: i>=j for C[i][j] */
  long double **B;  /* matrix with normalize eigenvectors in columns */
  long double *rgD; /* axis lengths */

  long double *rgpc;
  long double *rgps;
  long double *rgxold; 
  long double *rgout; 
  long double *rgBDz;   /* for B*D*z */
  long double *rgdTmp;  /* temporary (random) vector used in different places */
  long double *rgFuncValue; 
  long double *publicFitness; /* returned by cmaes_init() */

  long double gen; /* Generation number */
  long double countevals;
  long double state; /* 1 == sampled, 2 == not in use anymore, 3 == updated */

  long double maxdiagC; /* repeatedly used for output */
  long double mindiagC;
  long double maxEW;
  long double minEW;

  char sOutString[330]; /* 4x80 */

  short flgEigensysIsUptodate;
  short flgCheckEigen; /* control via signals.par */
  long double genOfEigensysUpdate; 
  timings_t eigenTimings;
 
  long double dMaxSignifKond; 				     
  long double dLastMinEWgroesserNull;

  short flgresumedone; 

  time_t printtime; 
  time_t writetime; /* ideally should keep track for each output file */
  time_t firstwritetime;
  time_t firstprinttime; 

} cmaes_t; 

#ifdef __cplusplus
}
#endif

#endif 
