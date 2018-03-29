#ifndef OPTNCOMPSSTATVALUE_H_
#define OPTNCOMPSSTATVALUE_H_

#include "CollectionLogStat.h"

class GAGenome;
class Algorithm;

class OptNCompsLogStat : public CollectionLogStat
{
  static GAGenome* opt_gen__;

  long double (*ncomps_func_)(GAGenome&);

  static long double normalOptNcompsFunc (GAGenome& gen);

public:
	OptNCompsLogStat(const Algorithm& alg, GAGenome* opt_gen, long double (*ncomps_func)(GAGenome&));
	~OptNCompsLogStat();

  long double indValue(GAGenome& ind);
};

#endif /*OPTNCOMPSSTATVALUE_H_*/
