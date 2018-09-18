#ifndef ENTROPYSTATVALUE_H_
#define ENTROPYSTATVALUE_H_

#include "SingleLogStat.h"

class EntropyLogStat : public SingleLogStat {
public:
	EntropyLogStat(const Algorithm& alg) : SingleLogStat("entropy",alg){}
	~EntropyLogStat() {}
	
	double computeValue() {return alg_.population().entropy();}
};

#endif /*ENTROPYSTATVALUE_H_*/
