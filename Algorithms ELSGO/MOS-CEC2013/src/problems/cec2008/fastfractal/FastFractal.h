#ifndef FASTFRACTAL_H_
#define FASTFRACTAL_H_

#include "UnitFunction1D.h"
#include "FractalFunction1D.h"

#include <vector>

class FastFractal {

 private:

  int dimensions_;
  UnitFunction1D* unitFunction_;
  FractalFunction1D* ff_;

 public:

  FastFractal();
  FastFractal(const std::string& unitFunctionName, int fractalDepth,
	      int density, long index, int dimensions);

  virtual ~FastFractal();

  long double              evaluate(std::vector<long double>                point );
  std::vector<long double> evaluate(std::vector< std::vector<long double> > points);

};

#endif /*FASTFRACTAL_H_*/
