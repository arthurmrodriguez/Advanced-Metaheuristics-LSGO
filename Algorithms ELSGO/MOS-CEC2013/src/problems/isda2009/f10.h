#ifndef F10_H_
#define F10_H_

#include <genomes/GA1DArrayGenome.h>

// F10 also known as Bohachevsky

const double f_Bohachevsky_BIAS = 0.0;

double f_Bohachevsky(GAGenome& g) {

   GA1DArrayAlleleGenome<double>& genome = DYN_CAST (GA1DArrayAlleleGenome<double>&, g);

   unsigned dim = genome.length ();

   const double PI = 3.141592653589793;

   double sum = 0.0;
   double currentGen = genome.gene (0);
   double nextGen;

   for (unsigned i = 1; i < dim; i++) {
      nextGen = genome.gene (i);
      sum += currentGen * currentGen + 2.0 * nextGen * nextGen;
      sum += -0.3 * cos (3.0 * PI * currentGen) - 0.4 * cos (4.0 * PI * nextGen) + 0.7;
      currentGen = nextGen;
   }

   return sum;

}

#endif /* F10_H_ */
