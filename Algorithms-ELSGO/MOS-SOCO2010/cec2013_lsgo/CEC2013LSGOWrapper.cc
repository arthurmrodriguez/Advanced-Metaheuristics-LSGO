#include "CEC2013LSGOWrapper.h"
#include "mclmcrrt.h"
#include <CEC2013LSGO.h>

#include <string.h>
#include <iostream>

int Initialize_CEC2013_LSGO_Wrapper(const char *path) {
  if( !mclInitializeApplication(NULL, 0)) {
    std::cout << "Error: Could not initialize the application!!!" << std::endl;
    return _HAS_ERROR;
  }

  // initialize lib
  if( !CEC2013LSGOInitialize()) {
    std::cout << "Error: Could not initialize Cost function!!!" << std::endl;
    return _HAS_ERROR;
  }

  static mxArray *input = mxCreateString(path);

  mlfInitializeBenchmark(input);

  return _NO_ERROR;
}


void Terminate_CEC2013_LSGO_Wrapper(void) {
  CEC2013LSGOTerminate();
}


void bench_func(long double *x, long double *f, unsigned fun_num, unsigned dims) {
  // Allocate pointers in Matlab to store input arguments (plhs)
  // and return values (prhs)
  static mxArray *plhs[2] = {0, 0};
  static mxArray *prhs[1] = {0};

  // Pointers to store references to the previously created
  // arguments and return values: solution (_x), function
  // number (_fun_num) and fitness value (_f)
  long double *_x, *_fun_num, *_f;

  // Allocate memory for each of the parameters/return values
  if(plhs[0] == 0)
    plhs[0] = mxCreateDoubleMatrix(1, dims, mxREAL);
  if(plhs[1] == 0)
    plhs[1] = mxCreateDoubleMatrix(1, 1,    mxREAL);
  if(prhs[0] == 0)
    prhs[0] = mxCreateDoubleMatrix(1, 1,    mxREAL);

  // Copy input arguments into Matlab variables
  _x       = mxGetPr(plhs[0]);
  _fun_num = mxGetPr(plhs[1]);

  *_fun_num = fun_num;

  for(int i = 0; i < dims; i++) {
    _x[i] = x[i];
  }

  // Call Matlab wrapper
  mlfBenchmark_func(1, prhs, plhs[0], plhs[1]);

  // Store return value
  _f = mxGetPr(prhs[0]);
  *f = *_f;
}

