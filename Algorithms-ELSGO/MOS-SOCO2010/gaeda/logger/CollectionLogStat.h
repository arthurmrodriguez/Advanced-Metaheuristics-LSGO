#ifndef POPSTATVALUE_H_
#define POPSTATVALUE_H_

#include <string>

#include "LogStat.h"

class GAGenome;
class GAPopulation;
class Algorithm;

class CollectionLogStat : public LogStat {

protected:

  virtual void computeValues(long double& max, long double& min, long double& avg, long double& dev);

public:

  CollectionLogStat(std::string name, const Algorithm& alg);

  ~CollectionLogStat();

  void update();
  virtual long double indValue(GAGenome& ind) = 0;

};

#endif /*POPSTATVALUE_H_*/
