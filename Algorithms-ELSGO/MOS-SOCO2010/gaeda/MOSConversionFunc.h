#ifndef MOSCONVERSIONFUNC_H
#define MOSCONVERSIONFUNC_H

#include "quicksort.h"

/**
 * Definicion de las funciones de conversion basicas para tsp
 */
extern "C" bool convertIntToReal (GAGenome* intGen, GAGenome* realGen) {
  GA1DArrayAlleleGenome<int>* _intGen = DYN_CAST(GA1DArrayAlleleGenome<int>*,intGen);
  GA1DArrayAlleleGenome<long double>* _realGen = DYN_CAST(GA1DArrayAlleleGenome<long double>*,realGen);

  unsigned len = _realGen->length();
  long double   inc = (_realGen->alleleset().upper () - _realGen->alleleset().lower ()) / (long double) len;

  for (register unsigned i = 0; i < len; i++) {
    long double low, up;
    low = _realGen->alleleset().lower() + (_intGen->gene(i) * inc);
    up  = _realGen->alleleset().lower() + ((_intGen->gene(i)+1) * inc);
    _realGen->gene(i, GARandomDouble (low, up));
  }
  return true;
}

extern "C" bool convertRealToInt (GAGenome* realGen, GAGenome* intGen) {
  GA1DArrayAlleleGenome<int>* _intGen = DYN_CAST(GA1DArrayAlleleGenome<int>*,intGen);
  GA1DArrayAlleleGenome<long double>* _realGen = DYN_CAST(GA1DArrayAlleleGenome<long double>*,realGen);

  unsigned len = _intGen->length();
  int      keys   [len];
  long double   values [len];

  // Prepare data
  for (register unsigned i=0; i<len; i++) {
    keys   [i] = i;
    values [i] = _realGen->gene(i);
  }

  // Sort all the values
  quicksort (keys, values, 0, len - 1);

  for (register unsigned i=0; i<len; i++) {
    _intGen->gene(keys[i], _intGen->alleleset().allele(i));
  }
  return true;
}

#endif
