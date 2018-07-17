#ifndef _F1_H
#define _F1_H

#include "Benchmarks.h"

class F1:public Benchmarks{
public:
	F1();
	long double compute(long double* x) ;
	long double compute(vector<long double> x) ;
	~F1();
};
#endif
