#ifndef _F2_H
#define _F2_H

#include "Benchmarks.h"


class F2:public Benchmarks{
protected:

public:
	F2();
	long double compute(long double* x) ;
	long double compute(vector<long double> x) ;
	~F2();
};
#endif
