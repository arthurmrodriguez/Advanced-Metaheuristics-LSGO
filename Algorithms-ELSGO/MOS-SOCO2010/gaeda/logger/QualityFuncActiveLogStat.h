#ifndef QUALITYFUNCACTIVESTATVALUE_H_
#define QUALITYFUNCACTIVESTATVALUE_H_

#include "../MOSTechnique.h"
#include "SingleLogStat.h"

class QualityFuncActiveLogStat : public SingleLogStat {
public:
   QualityFuncActiveLogStat(const Algorithm& alg) : SingleLogStat("qf_active",alg){}
   ~QualityFuncActiveLogStat() {}

   long double computeValue() {return (MOSTechnique::improvement_override ? 2.0 : 1.0);}
};

#endif /*QUALITYFUNCACTIVESTATVALUE_H_*/
