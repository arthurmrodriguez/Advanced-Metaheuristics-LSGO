#include "OptDistLogStat.h"

#include "../genomes/GAGenome.h"

GAGenome* OptDistLogStat::opt_gen__ = NULL;

long double OptDistLogStat::normalOptDistFunc (GAGenome& gen){
  return gen.compare(*opt_gen__);
}

OptDistLogStat::OptDistLogStat(const Algorithm& alg,
                               GAGenome*        opt_gen,
                               long double           (*dist_func)(GAGenome&)) :
                                                                    CollectionLogStat("opt_dist",alg),
                                                                    dist_func_       (dist_func ) {
  assert(dist_func_ != NULL || opt_gen != NULL);
  opt_gen__ = opt_gen;
  if (dist_func == NULL) dist_func_ = normalOptDistFunc;
}

OptDistLogStat::~OptDistLogStat(){}

long double OptDistLogStat::indValue(GAGenome& ind){
  return dist_func_(ind);
}
