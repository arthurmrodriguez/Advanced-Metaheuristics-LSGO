#ifndef _F18_H
#define _F18_H

#include "Benchmarks.h"

class F18:public Benchmarks{
protected:
	void generateInterArray ( );
public:
	F18();
	long double compute(long double* x) ;
	long double compute(vector<long double> x) ;
	~F18();
};

#endif

