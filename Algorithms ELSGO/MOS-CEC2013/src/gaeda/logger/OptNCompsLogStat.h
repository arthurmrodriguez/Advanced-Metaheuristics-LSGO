#ifndef OPTNCOMPSSTATVALUE_H_
#define OPTNCOMPSSTATVALUE_H_

#include "CollectionLogStat.h"

class GAGenome;
class Algorithm;

class OptNCompsLogStat : public CollectionLogStat
{
  static GAGenome* opt_gen__;

  double (*ncomps_func_)(GAGenome&);

  static double normalOptNcompsFunc (GAGenome& gen);

public:
	OptNCompsLogStat(const Algorithm& alg, GAGenome* opt_gen, double (*ncomps_func)(GAGenome&));
	~OptNCompsLogStat();

  double indValue(GAGenome& ind);
};

#endif /*OPTNCOMPSSTATVALUE_H_*/
