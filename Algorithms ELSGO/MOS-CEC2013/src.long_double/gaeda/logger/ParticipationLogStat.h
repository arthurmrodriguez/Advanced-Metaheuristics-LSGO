#ifndef PARTICIPATIONSTATVALUE_H_
#define PARTICIPATIONSTATVALUE_H_

#include "SetLogStat.h"

class ParticipationLogStat : public SetLogStat {
 public:
  ParticipationLogStat(const Algorithm& alg);
  ~ParticipationLogStat() {}

  void computeValues(vector<long double>*);
};

#endif /*PARTICIPATIONSTATVALUE_H_*/
