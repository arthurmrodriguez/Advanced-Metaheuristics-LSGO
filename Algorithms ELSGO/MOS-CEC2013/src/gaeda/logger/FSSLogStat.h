#ifndef FSSSTATVALUE_H_
#define FSSSTATVALUE_H_

#include "SetLogStat.h"

class FSSLogStat : public SetLogStat {
 public:

  enum types {AVG, BEST, WORST};

  FSSLogStat(const Algorithm& alg, int type);
  ~FSSLogStat() {}

  void computeValues(vector<double>*);

 protected:
  int _algId;
  int _type;
};

#endif /*FSSSTATVALUE_H_*/
