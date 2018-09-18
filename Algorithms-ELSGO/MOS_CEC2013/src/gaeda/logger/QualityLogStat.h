#ifndef QUALITYSTATVALUE_H_
#define QUALITYSTATVALUE_H_

#include "SetLogStat.h"

class QualityLogStat : public SetLogStat {
 public:
  QualityLogStat(const Algorithm& alg);
  ~QualityLogStat() {}

  void computeValues(vector<double>*);
};

#endif /*QUALITYSTATVALUE_H_*/
