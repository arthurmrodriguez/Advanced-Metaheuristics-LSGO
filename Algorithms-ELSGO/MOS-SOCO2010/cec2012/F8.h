#ifndef _F8_H
#define _F8_H

#include "Benchmarks.h"

class F8:public Benchmarks{
protected:
public:
	F8();
	long double compute(long double* x) ;
	long double compute(vector<long double> x) ;
	~F8();
};

#endif
