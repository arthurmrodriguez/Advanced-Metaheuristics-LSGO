#include <stdio.h>
#include "tests.h"
#include <math.h>

#ifndef ANSI_ARGS
#ifdef __STDC__
#define ANSI_ARGS(args) args
#else
#define ANSI_ARGS(args) ()
#endif
#endif

long double xt ANSI_ARGS((long double x));
long double maxt ANSI_ARGS((long n, int t));
void set_t ANSI_ARGS((int t));
long double *V2;

int xt_t = 1;


#ifdef __STDC__
main(int argc, char *argv[])
#else
main(argc, argv)
int argc;
char *argv[];
#endif
{
  long ntests, n, i;
  long double *V, result;
  int t;
  
  if(argc != N_STREAM_PARAM + 3)
  {
    fprintf(stderr,"USAGE: %s (... %d arguments)\n",argv[0], N_STREAM_PARAM+2);
    exit(1);
  }
  
  ntests = init_tests(argc,argv);
  
  n = atol(argv[N_STREAM_PARAM+1]);
  t = atoi(argv[N_STREAM_PARAM+2]);
  
  V = (long double *) mymalloc(NTESTS*sizeof(long double));
  V2 = (long double *) mymalloc(n*sizeof(long double));
  
  for(i=0; i<ntests; i++)
  {
    V[i] = maxt(n,t);
    
    next_stream();
  }
  
  set_KS_n(NTESTS);

#if defined(SPRNG_MPI)
  getKSdata(V,NTESTS);
#endif
  
  if(proc_rank == 0)
  {
    result = KS(V,NTESTS,KSF);
    printf("\nResult: KS value = %f", result);
    result = KSpercent(result,NTESTS);
    printf("\t %% = %.2f\n\n", result*100.0);
  }
  

#if defined(SPRNG_MPI)
     MPI_Finalize();
#endif

}


#ifdef __STDC__
long double xt(long double x)
#else
long double xt(x)
long double x;
#endif
{
  return pow(x,(long double) xt_t);
}

#ifdef __STDC__
void set_t(int t)
#else
void set_t(t)
int t;
#endif
{
  xt_t = t;
}

#ifdef __STDC__
long double maxt(long n, int t)
#else
long double maxt(n, t)
int t;
long n;
#endif
{
  long double *V=V2, temp;
  int j;
  long i;
  
  for(i=0; i<n; i++)
  {
    V[i] = 0.0;
    for(j=0; j<t; j++)
    {
      temp=get_rn();
      
      if( temp > V[i])
	V[i] = temp;
    }
  }
  
  set_t(t);
  temp = KS(V,n,xt);
  /*printf("\tKS for stream = %f, %% = %f\n", temp, KSpercent(temp,n));*/
  
  /*free(V);*/
  
  return temp;
}
