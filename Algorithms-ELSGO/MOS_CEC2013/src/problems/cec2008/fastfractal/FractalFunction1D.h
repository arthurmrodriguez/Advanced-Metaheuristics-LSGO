#ifndef FRACTALFUNCTION_H_
#define FRACTALFUNCTION_H_

#include "RanTable.h"
#include "UnitFunction1D.h"

class FractalFunction1D {

 private:

  static int DOUBLE_TABLE_SIZE_;   // 16384
  static int INT_TABLE_SIZE_;      //   256

  int fractalDepth_;  // maximum recursive depth of 40 suggested for 64-bit architecture
  int density_;
  long index_;
  bool func_;
  RanTable* ranTable_;
  UnitFunction1D* unitFunction_;

 public:

  FractalFunction1D();
  FractalFunction1D(UnitFunction1D* unitFunction, int fractalDepth, int density, long index);
  FractalFunction1D(UnitFunction1D* unitFunction, int density, long index);
  FractalFunction1D(UnitFunction1D* unitFunction, long index);
  FractalFunction1D(UnitFunction1D* unitFunction);

  virtual ~FractalFunction1D();

  void setIndex(long index);

  long double evaluate(long double x);
  long double getDepthLocal(long double x, int recDepth, long seed, long span);
  long double getDepthWRTSquare(long double x, long square, int recDepth, long seed, long span, long double scale);

};

#endif /*FRACTALFUNCTION1D_H_*/
