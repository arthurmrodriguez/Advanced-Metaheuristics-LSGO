#include "RanQD1.h"

#include <math.h>
#include <iostream>
RanQD1::RanQD1()
{
}

RanQD1::RanQD1 (long seed)
  :idum_(seed)
{
  nextLong(); // one multiple to scatter seeds
}

RanQD1::~RanQD1()
{
}

void RanQD1::setSeed(long seed) {
  idum_ = seed;
  nextLong(); // one multiple to scatter seeds
}

long long RanQD1::nextLong() {
  idum_ = (A_ * idum_ + C_) & MASK_;
  return idum_;
}
 
long double RanQD1::nextDouble() {
  return nextLong()/MAX_INT_;
}

int RanQD1::nextInt(int min, int max) {
  return min + (int) floor(nextDouble()*(max-min+1));
}

long long   RanQD1::MASK_    = 0xffffffffl;  // lower order 32 bits of long
long double RanQD1::MAX_INT_ = 4294967295.0; // 2^32-1 (MASK as a long double)
long long   RanQD1::A_       = 1664525l;     // suggested by Knuth
long long   RanQD1::C_       = 1013904223l;  // suggested by Lewis
