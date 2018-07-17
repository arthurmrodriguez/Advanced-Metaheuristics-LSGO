/**
  * @file
  * @brief Common genetic operators.
  *
  */

#ifndef GACommonOps_H
#define GACommonOps_H

#include "genomes/GA1DArrayGenome.h"

// Initializers

// The random initializer sets the elements of the array based on the alleles
// set.  We choose randomly the allele for each element.
template <typename T> void UniformInitializer (GAGenome& c) {

   GA1DArrayAlleleGenome<T>& child = DYN_CAST (GA1DArrayAlleleGenome<T>&, c);

   for (int i = child.length () - 1; i >= 0; i--)
      child.gene (i, child.alleleset (i).allele ());

}

// DE Crossovers

template <typename T> int ExponentialCrossover (const GAGenome& p1, const GAGenome& p2, GAGenome* c1, const long double prob) {

   const GA1DArrayGenome<T>& x_i    = DYN_CAST (const GA1DArrayGenome<T>&, p1);
   GA1DArrayGenome<T>&       x_new  = DYN_CAST (GA1DArrayGenome<T>&, *c1);

   int dim   = x_i.length();
   int nrand = GARandomInt(0,dim-1);
   int lrand = 0;

   while (GARandomDouble(0,1) < prob and (lrand<dim) )
    lrand = lrand+1;
   if (lrand==0)
    lrand=1;

   int initial_pos = ( nrand + lrand ) % dim;
   while (initial_pos != nrand) {
      x_new.gene(initial_pos,x_i.gene(initial_pos) );
      initial_pos = (initial_pos + 1) % dim;
   }

   return 1;

}


#endif
