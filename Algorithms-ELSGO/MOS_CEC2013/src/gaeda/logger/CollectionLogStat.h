#ifndef POPSTATVALUE_H_
#define POPSTATVALUE_H_

#include <string>

#include "LogStat.h"

class GAGenome;
class GAPopulation;
class Algorithm;

class CollectionLogStat : public LogStat {

protected:

  virtual void computeValues(double& max, double& min, double& avg, double& dev);

public:

  CollectionLogStat(std::string name, const Algorithm& alg);

  ~CollectionLogStat();

  void update();
  virtual double indValue(GAGenome& ind) = 0;

};

#endif /*POPSTATVALUE_H_*/
