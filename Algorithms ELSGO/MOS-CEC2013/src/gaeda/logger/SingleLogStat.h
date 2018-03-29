#ifndef SINGLESTATVALUE_H_
#define SINGLESTATVALUE_H_

#include <string>

#include "LogStat.h"

using namespace std;

class SingleLogStat : public LogStat {
protected:
  virtual double computeValue() = 0;
public:
  SingleLogStat (string name,const Algorithm& alg);

  ~SingleLogStat();

  void update();
};

#endif /*SINGLESTATVALUE_H_*/
