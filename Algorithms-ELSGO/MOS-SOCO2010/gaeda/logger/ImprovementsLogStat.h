#ifndef IMPROVEMENTSSTATVALUE_H_
#define IMPROVEMENTSSTATVALUE_H_

#include "SetLogStat.h"

class ImprovementsLogStat : public SetLogStat {
 public:
  ImprovementsLogStat(const Algorithm& alg);
  ~ImprovementsLogStat() {}

  void computeValues(vector<long double>*);
};

#endif /*IMPROVEMENTSSTATVALUE_H_*/
