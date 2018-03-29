#ifndef F8_H_
#define F8_H_

#include <genomes/GA1DArrayGenome.h>

// F8 also known as Schwefel's Problem 1.2

const double f_Schwefel1_2_BIAS = 0.0;

double f_Schwefel1_2(GAGenome& g) {

   GA1DArrayAlleleGenome<double>& genome = DYN_CAST (GA1DArrayAlleleGenome<double>&, g);

   unsigned dim = genome.length ();

   double sum = 0.0, val = 0.0;

   for (unsigned i = 0; i < dim; i++) {
      val += genome.gene(i);
      sum += val * val;
   }

   return sum;

}

#endif  /* F8_H_ */
