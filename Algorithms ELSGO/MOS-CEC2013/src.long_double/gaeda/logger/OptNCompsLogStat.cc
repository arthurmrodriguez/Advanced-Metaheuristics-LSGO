#include "OptNCompsLogStat.h"

#include "../genomes/GAGenome.h"

GAGenome* OptNCompsLogStat::opt_gen__ = NULL;

long double OptNCompsLogStat::normalOptNcompsFunc (GAGenome& gen){
  return gen.compsCompare(*opt_gen__);
}

OptNCompsLogStat::OptNCompsLogStat(const Algorithm& alg,
                                   GAGenome*        opt_gen,
                                   long double           (*ncomps_func)(GAGenome&)) :
                                                               CollectionLogStat("opt_ncomps",alg),
                                                               ncomps_func_(ncomps_func) {
  assert(ncomps_func_ != NULL || opt_gen != NULL);
  opt_gen__ = opt_gen;
  if (ncomps_func_ == NULL) ncomps_func_ = normalOptNcompsFunc;
}

OptNCompsLogStat::~OptNCompsLogStat(){}

long double OptNCompsLogStat::indValue(GAGenome& ind){
  return ncomps_func_(ind);
}
