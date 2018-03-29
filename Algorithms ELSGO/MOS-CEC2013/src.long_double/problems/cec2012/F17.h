#ifndef _F17_H
#define _F17_H

#include "Benchmarks.h"

class F17:public Benchmarks{
protected:
	void generateInterArray ( );
public:
	F17();
	long double compute(long double* x) ;
	long double compute(vector<long double> x) ;
	~F17();
};

#endif

