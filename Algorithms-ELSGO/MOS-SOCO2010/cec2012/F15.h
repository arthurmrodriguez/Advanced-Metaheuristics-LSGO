#ifndef _F15_H
#define _F15_H

#include "Benchmarks.h"

class F15:public Benchmarks{
protected:
	void generateInterArray ( );
public:
	F15();
	long double compute(long double* x) ;
	long double compute(vector<long double> x) ;
	~F15();
};

#endif
