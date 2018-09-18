#ifndef ROUTINGALGSCORESTATVALUE_H_
#define ROUTINGALGSCORESTATVALUE_H_

#include "SetLogStat.h"

class RoutingAlgScoreLogStat : public SetLogStat {
 public:
  RoutingAlgScoreLogStat(const Algorithm& alg);
  virtual ~RoutingAlgScoreLogStat() {}

  void computeValues(vector<double>*);
};

#endif /*ROUTINGALGSCORESTATVALUE_H_*/
