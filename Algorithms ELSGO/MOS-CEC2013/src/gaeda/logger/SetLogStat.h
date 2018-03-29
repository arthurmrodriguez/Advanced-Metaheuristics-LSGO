#ifndef SETLOGSTAT_H_
#define SETLOGSTAT_H_

#include <vector>

#include "LogStat.h"

class Algorithm;

using namespace std;

class SetLogStat : public LogStat {
protected:
  virtual void computeValues(vector<double>*) = 0;
  unsigned _setSize;

public:
  SetLogStat (string name,const Algorithm& alg);

  ~SetLogStat();

  void update();
};

#endif /*SETLOGSTAT_H_*/
