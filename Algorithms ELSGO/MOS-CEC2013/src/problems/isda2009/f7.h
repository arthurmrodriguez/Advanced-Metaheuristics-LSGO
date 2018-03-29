#ifndef F7_H_
#define F7_H_

#include <genomes/GA1DArrayGenome.h>

// F7 also known as Schwefel's Problem 2.22

const double f_Schwefel2_22_BIAS = 0.0;

double f_Schwefel2_22(GAGenome& g) {

   GA1DArrayAlleleGenome<double>& genome = DYN_CAST (GA1DArrayAlleleGenome<double>&, g);

   unsigned dim = genome.length ();

   double sum, currentGen, prod;

   sum  = 0.0;
   prod = 1.0;

   for (unsigned i = 0; i < dim; i++) {
      currentGen = fabs (genome.gene (i));
      sum  += currentGen;
      prod *= currentGen;
   }

   return (sum + prod);

}

#endif /* F7_H_ */
