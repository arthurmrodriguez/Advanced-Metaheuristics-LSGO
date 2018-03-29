#ifndef DOUBLEDIP_H_
#define DOUBLEDIP_H_

#include "UnitFunction1D.h"

class DoubleDip : public UnitFunction1D {

 public:

  DoubleDip();
  DoubleDip(long double centre, long double scale);

  ~DoubleDip();

  virtual long double getValue(long double point);
  virtual long double twist(long double x, long double y);

};

#endif /*DOUBLEDIP_H_*/
