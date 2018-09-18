#include "FastFractal.h"
#include "DoubleDip.h"

#include <stdlib.h>

#include <stdexcept>
#include <iostream>

FastFractal::FastFractal()
{
}

FastFractal::FastFractal (const std::string& unitFunctionName, int fractalDepth,
			  int density, long index, int dimensions)
  : dimensions_ (dimensions)
{

  if (unitFunctionName == "DoubleDip")
    unitFunction_ = new DoubleDip();
  else {
    std::cerr << "Error (FastFractal): Unknown Unit1D function name." << std::endl;
    exit(1);
  }

  ff_ = new FractalFunction1D(unitFunction_, fractalDepth, density, index);

}

FastFractal::~FastFractal()
{
  delete unitFunction_;
  delete ff_;
}

long double FastFractal::evaluate (std::vector<long double> point) {

  if ((int) point.size() != dimensions_)
    throw std::runtime_error("Point does not have the proper number of dimensions.");

  long double depth = 0;
  long double x, lastx, dx;

  ff_->setIndex((6*dimensions_-1)+1);
  lastx = point[dimensions_-1];

  for (int i = 0; i < dimensions_; i++) {
    ff_->setIndex(6*i + 1);                   // spread to small "prime-ish" seeds for diversity
    x = point[i];
    dx = unitFunction_->twist(x, lastx);
    depth = depth + ff_->evaluate(x + dx);    // "twist" and evaluate
    lastx = x;
  }

  return depth;

}

std::vector<long double> FastFractal::evaluate (std::vector< std::vector<long double> > points) {

  if ((int) points[0].size() != dimensions_)
    throw std::runtime_error("Point does not have the proper number of dimensions.");

  std::vector<long double> results (points.size(), 0.0);

  for (unsigned i = 0; i < points.size(); i++)
    results[i] = evaluate(points[i]);

  return results;

}
