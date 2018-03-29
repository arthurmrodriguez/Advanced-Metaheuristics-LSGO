#ifndef _F19_H
#define _F19_H

#include "Benchmarks.h"

class F19:public Benchmarks{
protected:
public:
	F19();
	long double compute(long double* x) ;
	long double compute(vector<long double> x) ;
	~F19();
};

#endif


