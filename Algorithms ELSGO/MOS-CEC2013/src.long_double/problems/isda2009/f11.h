#ifndef F11_H_
#define F11_H_

#include <genomes/GA1DArrayGenome.h>

// F11 also known as Schaffer

const long double f_Schaffer_BIAS = 0.0;

long double f_Schaffer(GAGenome& g) {

   GA1DArrayAlleleGenome<long double>& genome = DYN_CAST (GA1DArrayAlleleGenome<long double>&, g);

   unsigned dim = genome.length ();

   long double sum, aux, aux2;
   long double currentGen, nextGen;

   sum = 0.0;
   currentGen = genome.gene (0);
   currentGen = currentGen * currentGen;

   for (unsigned i = 1; i < dim; i++) {
      nextGen = genome.gene (i);
      nextGen = nextGen * nextGen;
      aux = currentGen + nextGen;
      currentGen = nextGen;
      aux2 = sin (50.0 * pow (aux, 0.1));
      sum += pow (aux, 0.25) * (aux2 * aux2 + 1.0);
   }

   return sum;

}

#endif /* F11_H_ */
