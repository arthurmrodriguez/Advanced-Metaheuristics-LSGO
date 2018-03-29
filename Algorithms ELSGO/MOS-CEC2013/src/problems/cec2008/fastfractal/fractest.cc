#include "FastFractal.h"

#include <iostream>
#include <vector>

int main (int argc, char* argv[]) {

  int dim = 1000;
  std::vector<long double> x (dim, 0.0);

  for (int i = 0; i < dim; i++)
    x[i] = 0.12323;

  FastFractal* ff = new FastFractal("DoubleDip", 3, 1, 1, dim);

  long double f = 0;

  f = ff->evaluate (x);
  std::cout << "f@0=" << f << std::endl;

  f = ff->evaluate (x);
  x[1] = 0.3;
  std::cout << "f@0=" << f << std::endl;

  f = ff->evaluate (x);
  std::cout << "f@0=" << f << std::endl;

  f = ff->evaluate (x);
  std::cout << "f@0=" << f << std::endl;

  x[1] = 0.5;
  f = ff->evaluate (x);
  std::cout << "f@0=" << f << std::endl;

  f = ff->evaluate (x);
  std::cout << "f@0=" << f << std::endl;

  return 0;
}
