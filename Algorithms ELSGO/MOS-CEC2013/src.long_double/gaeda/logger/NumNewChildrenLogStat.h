#ifndef NUMNEWCHILDRENLOGSTAT_
#define NUMNEWCHILDRENLOGSTAT_

#include "SingleLogStat.h"
#include "../GAStatistics.h"

class NumNewChildrenLogStat : public SingleLogStat {
public:
  NumNewChildrenLogStat (const Algorithm& alg) : SingleLogStat("num_new_child",alg) {}
  ~NumNewChildrenLogStat(                    ){}

  long double computeValue(){ return alg_.statistics().numNewChild(); }
};

#endif /*NUMNEWCHILDRENLOGSTAT_*/
