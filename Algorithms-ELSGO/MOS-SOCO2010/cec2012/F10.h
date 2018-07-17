#ifndef _F10_H
#define _F10_H

#include "Benchmarks.h"

class F10:public Benchmarks{
protected:
public:
	F10();
	long double compute(long double* x) ;
	long double compute(vector<long double> x) ;
	~F10();
};

#endif
