#ifndef _CEC2013_LSGO_WRAPPER_H
#define _CEC2013_LSGO_WRAPPER_H

const int _HAS_ERROR = 1;
const int  _NO_ERROR = 0;

int  Initialize_CEC2013_LSGO_Wrapper(const char *path);
void  Terminate_CEC2013_LSGO_Wrapper(void);

void bench_func(double *x, double *f, unsigned fun_num, unsigned dims);

#endif
