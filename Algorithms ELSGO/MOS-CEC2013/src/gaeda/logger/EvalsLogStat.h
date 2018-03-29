#ifndef EVALSSTATVALUE_H_
#define EVALSSTATVALUE_H_

#include "SingleLogStat.h"

class EvalsLogStat : public SingleLogStat {
public:
   EvalsLogStat(const Algorithm& alg) : SingleLogStat("evaluations",alg){}
   ~EvalsLogStat() {}

   double computeValue() {return alg_.statistics().indEvals();}
};

#endif /*EVALSSTATVALUE_H_*/
