#ifndef _F14_H
#define _F14_H

#include "Benchmarks.h"

class F14:public Benchmarks{
protected:
	void generateInterArray ( );
public:
	F14();
	long double compute(long double* x) ;
	long double compute(vector<long double> x) ;
	~F14();
};

#endif

