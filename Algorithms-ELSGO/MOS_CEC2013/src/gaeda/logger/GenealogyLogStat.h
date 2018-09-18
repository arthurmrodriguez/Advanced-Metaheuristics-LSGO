#ifndef NSCSTATVALUE_H_
#define NSCSTATVALUE_H_

#include "LogStat.h"

class GenealogyLogStat : public LogStat {
public:
  GenealogyLogStat(const Algorithm& alg);
  ~GenealogyLogStat() {}
  
  void computeValues(double &max, double &min, double &avg, double &dev, double &imperc, double &immax, double &immin, double &imavg, double &imdev, double &fdac);
  void update();

};

#endif /*NSCSTATVALUE_H_*/
