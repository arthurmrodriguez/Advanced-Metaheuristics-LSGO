#ifndef _F20_H
#define _F20_H

#include "Benchmarks.h"

class F20:public Benchmarks{
protected:
public:
	F20();
	long double compute(long double* x) ;
	long double compute(vector<long double> x) ;
	~F20();
};

#endif

