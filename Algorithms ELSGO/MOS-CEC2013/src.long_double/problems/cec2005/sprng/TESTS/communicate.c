#ifdef SPRNG_MPI
#include <mpi.h>
#endif

#include "util.h"

#define SPRNG_MPI_KS 0

extern int proc_rank, init_nprocs, nsubsequences;


#ifdef __STDC__
void getKSdata (long double *V, long ntests)
#else
void getKSdata (V, ntests)
long double *V;
long ntests;
#endif
{
  int nprocs, myid, i, nsets;
  long first, last;
  
#ifdef SPRNG_MPI
  MPI_Status status;
  
  myid = proc_rank;
  nprocs = init_nprocs;
  nsets = ntests/nsubsequences;
  
  first = (nsets*myid/nprocs)*nsubsequences;
  last = (nsets*(myid+1)/nprocs)*nsubsequences - 1;
  
  if(myid != 0)
    MPI_Send(V, last-first+1, MPI_DOUBLE, 0, SPRNG_MPI_KS, MPI_COMM_WORLD);
  else
    for(i=1; i<nprocs; i++)
      MPI_Recv(V+(nsets*i/nprocs)*nsubsequences,
	       (nsets*(i+1)/nprocs - nsets*i/nprocs)*nsubsequences,
	       MPI_DOUBLE, i, SPRNG_MPI_KS, MPI_COMM_WORLD, &status);
  
  /*MPI_Finalize();*/
#endif
}



#ifdef __STDC__
void  reduce_sum(long double *array,int n)
#else
void  reduce_sum(array,n)
long double *array;
int n;
#endif
{
  long double *temp;
  
#ifdef SPRNG_MPI
  temp = (long double *) mymalloc(n*sizeof(long double));

  MPI_Reduce(array, temp, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  if(proc_rank == 0)
    memcpy(array,temp,n*sizeof(long double));

  free(temp);
#endif
}

  


void  reduce_sum_long(long *array,long n)
{
  long *temp;
  
#ifdef SPRNG_MPI
  temp = (long *) mymalloc(n*sizeof(long));

  MPI_Reduce(array, temp, n, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  
  if(proc_rank == 0)
    memcpy(array,temp,n*sizeof(long));

  free(temp);
#endif
}

  
