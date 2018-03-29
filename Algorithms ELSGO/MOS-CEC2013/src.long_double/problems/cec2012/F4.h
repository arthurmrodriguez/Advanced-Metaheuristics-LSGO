#ifndef _F4_H
#define _F4_H

#include "Benchmarks.h"

class F4:public Benchmarks{
protected:
public:
	F4();
	long double compute(long double* x) ;
	long double compute(vector<long double> x) ;
	~F4();
};

#endif
