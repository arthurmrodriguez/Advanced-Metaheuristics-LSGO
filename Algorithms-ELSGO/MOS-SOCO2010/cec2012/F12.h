#ifndef _F12_H
#define _F12_H

#include "Benchmarks.h"

class F12:public Benchmarks{
protected:
public:
	F12();
	long double compute(long double* x) ;
	long double compute(vector<long double> x) ;
	~F12();
};

#endif
