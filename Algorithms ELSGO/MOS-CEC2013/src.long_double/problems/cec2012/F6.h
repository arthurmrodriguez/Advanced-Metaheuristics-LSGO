#ifndef _F6_H
#define _F6_H

#include "Benchmarks.h"

class F6:public Benchmarks{
protected:
public:
	F6();
	long double compute(long double* x) ;
	long double compute(vector<long double> x) ;
	~F6();
};

#endif
