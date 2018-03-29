#ifndef _F5_H
#define _F5_H

#include "Benchmarks.h"

class F5:public Benchmarks{
protected:
public:
	F5();
	long double compute(long double* x) ;
	long double compute(vector<long double> x) ;
	~F5();
};

#endif
