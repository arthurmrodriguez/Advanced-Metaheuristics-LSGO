#ifndef RANTABLE_H_
#define RANTABLE_H_

#include "RanQD1.h"

class RanTable {

 private:

  long double* doubleTable_;
  int doubleTableIndex_;
  int doubleTableSize_;

  int* intTable_;
  int intTableIndex_;
  int intTableSize_;

  RanQD1* ran;

 public:

  RanTable(int doubleTableSize, int intTableSize, int aveInt, long index);

  virtual ~RanTable();

  void setSeed(long seed);

  long double nextDouble ();
  int         nextInteger();

};

#endif /*RANTABLE_H_*/
