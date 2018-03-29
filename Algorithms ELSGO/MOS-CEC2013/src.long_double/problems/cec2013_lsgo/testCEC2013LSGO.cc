#include <iostream>
#include <limits>
#include <stdlib.h>

#include <armadillo>
#include <matio.h>

#include "cec2013_funcs.h"

unsigned maxIters=1;

bool initial_flag=false;
mat_t* matfp = NULL;
arma::mat data;

void   init    (unsigned func);
long double evaluate(unsigned func);

int main(int argc, char** argv) {
  matfp=Mat_Open(argv[2], MAT_ACC_RDONLY);

  unsigned func=strtol(argv[1], NULL, 10);

  init(func);

  long double fit=0.0;
  for (unsigned i=0; i<maxIters; i++)
    fit=evaluate(func);

  std::cout << "f" << func << "(data)=" << fit << std::endl;

  Mat_Close(matfp);

  return 0;
}

void init(unsigned func) {
  switch(func) {
    case 1:
      init_f1(matfp);
      data=arma::zeros(1000,1);
      break;
    case 2:
      init_f2(matfp);
      data=arma::zeros(1000,1);
      break;
    case 3:
      init_f3(matfp);
      data=arma::zeros(1000,1);
      break;
    case 4:
      init_f4(matfp);
      data=arma::zeros(1000,1);
      break;
    case 5:
      init_f5(matfp);
      data=arma::zeros(1000,1);
      break;
    case 6:
      init_f6(matfp);
      data=arma::zeros(1000,1);
      break;
    case 7:
      init_f7(matfp);
      data=arma::zeros(1000,1);
      break;
    case 8:
      init_f8(matfp);
      data=arma::zeros(1000,1);
      break;
    case 9:
      init_f9(matfp);
      data=arma::zeros(1000,1);
      break;
    case 10:
      init_f10(matfp);
      data=arma::zeros(1000,1);
      break;
    case 11:
      init_f11(matfp);
      data=arma::zeros(1000,1);
      break;
    case 12:
      init_f12(matfp);
      data=arma::zeros(1000,1);
      break;
    case 13:
      init_f13(matfp);
      data=arma::zeros(905,1);
      break;
    case 14:
      init_f14(matfp);
      data=arma::zeros(905,1);
      break;
    case 15:
      init_f15(matfp);
      data=arma::zeros(1000,1);
      break;
    default:
      break;
  }
}

long double evaluate(unsigned func) {
  switch(func) {
    case 1:
      return f1(data);
    case 2:
      return f2(data);
    case 3:
      return f3(data);
    case 4:
      return f4(data);
    case 5:
      return f5(data);
    case 6:
      return f6(data);
    case 7:
      return f7(data);
    case 8:
      return f8(data);
    case 9:
      return f9(data);
    case 10:
      return f10(data);
    case 11:
      return f11(data);
    case 12:
      return f12(data);
    case 13:
      return f13(data);
    case 14:
      return f14(data);
    case 15:
      return f15(data);
    default:
      return 0.0;
      break;
  }
}
