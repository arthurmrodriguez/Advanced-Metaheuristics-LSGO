#ifndef _F13_H
#define _F13_H

#include "Benchmarks.h"

class F13:public Benchmarks{
protected:
public:
	F13();
	long double compute(long double* x) ;
	long double compute(vector<long double> x) ;
	~F13();
};

#endif

