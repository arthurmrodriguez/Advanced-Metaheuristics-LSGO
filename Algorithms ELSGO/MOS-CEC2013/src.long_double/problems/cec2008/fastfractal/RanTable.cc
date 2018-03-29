#include "RanTable.h"
#include <iostream>
RanTable::RanTable(int long doubleTableSize, int intTableSize, int aveInt, long index) {

  ran = new RanQD1(index);

  doubleTableSize_  = doubleTableSize;
  doubleTableIndex_ = 0;
  doubleTable_      = new long double[doubleTableSize_];

  for (int i = 0; i < doubleTableSize_; i++)
    doubleTable_[i] = ran->nextDouble();

  ran->setSeed(index);

  intTableSize_  = intTableSize;
  intTableIndex_ = 0;
  intTable_      = new int[intTableSize_];

  for (int i = 0; i < intTableSize_; i++)
    intTable_[i] = ran->nextInt(0, 2*aveInt);

  ran->setSeed(index);

}

RanTable::~RanTable()
{
  delete ran;
  delete [] doubleTable_;
  delete [] intTable_;
}

void RanTable::setSeed(long seed) {
  doubleTableIndex_ = (int) (seed & (doubleTableSize_-1));
  intTableIndex_    = (int) (seed & (intTableSize_-1));
}

long double RanTable::nextDouble() {
  doubleTableIndex_ = (doubleTableIndex_+1) & (doubleTableSize_-1);
  return doubleTable_[doubleTableIndex_];
}
int  RanTable::nextInteger() {
  intTableIndex_ = (intTableIndex_+1) & (intTableSize_-1);
  return intTable_[intTableIndex_];
}
