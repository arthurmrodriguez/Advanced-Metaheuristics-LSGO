#ifndef _F11_H
#define _F11_H

#include "Benchmarks.h"

class F11:public Benchmarks{
protected:
public:
	F11();
	long double compute(long double* x) ;
	long double compute(vector<long double> x) ;
	~F11();
};

#endif
