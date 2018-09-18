#include "GAIntOps.h"


// ES Crossovers

int BinaryDominantCrossover (const std::vector<GAGenome*> parents, GAGenome* child) {

  GA1DArrayAlleleGenome<int>& g = DYN_CAST(GA1DArrayAlleleGenome<int>&, *child);

  unsigned ro = parents.size();

  std::vector<GA1DArrayAlleleGenome<int>*> pars (ro, (GA1DArrayAlleleGenome<int>*)0);

  for (unsigned i = 0; i < ro; i++)
    pars [i] = DYN_CAST(GA1DArrayAlleleGenome<int>*, parents[i]);

  // Recombine both endogenous parameters and variables
  for (int i = 0; i < g.size(); i++) {

    unsigned p = GARandomInt (0, ro-1);
    g.gene  (i, pars[p]->gene  (i));

  }

  return 1;

}


int BinaryDominantCrossoverUpdateOnly (const std::vector<GAGenome*> parents, GAGenome* child) {
   return 0;
}


// ES Mutators

int BinaryFlipMutator (GAGenome& g) {

   GA1DArrayAlleleGenome<int>& child = DYN_CAST (GA1DArrayAlleleGenome<int>&, g);

   unsigned size = child.length();
   double pMut = 1.0 / (double) size;
   int nMut = 0;

   if (pMut <= 0.0)
      return 0;

   for (unsigned i = 0; i < size; i++) {

      if (GAFlipCoin(pMut)) {

         int val = child.gene (i);
         child.gene (i, val == 1 ? 0 : 1);

         nMut++;

      }

   }

   return nMut;

}

int BinaryFlipMutatorUpdateOnly (GAGenome& g) {
   return 0;
}
