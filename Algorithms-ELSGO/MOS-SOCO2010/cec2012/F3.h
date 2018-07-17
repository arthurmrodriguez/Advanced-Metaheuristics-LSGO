#ifndef _F3_H
#define _F3_H

#include "Benchmarks.h"


class F3:public Benchmarks{
protected:

public:
	F3();
	long double compute(long double* x) ;
	long double compute(vector<long double> x) ;
	~F3();
};
#endif
