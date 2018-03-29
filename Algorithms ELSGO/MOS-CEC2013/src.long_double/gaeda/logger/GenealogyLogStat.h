#ifndef NSCSTATVALUE_H_
#define NSCSTATVALUE_H_

#include "LogStat.h"

class GenealogyLogStat : public LogStat {
public:
  GenealogyLogStat(const Algorithm& alg);
  ~GenealogyLogStat() {}
  
  void computeValues(long double &max, long double &min, long double &avg, long double &dev, long double &imperc, long double &immax, long double &immin, long double &imavg, long double &imdev, long double &fdac);
  void update();

};

#endif /*NSCSTATVALUE_H_*/
