#include "RoutingAlgScoreLogStat.h"

#include "../RoutingAlg.h"
#include "../MOSTechniqueSet.h"

RoutingAlgScoreLogStat::RoutingAlgScoreLogStat(const Algorithm& alg) : SetLogStat("sbest,s",alg) {
   _setSize = 2;
}

void RoutingAlgScoreLogStat::computeValues(vector<double>* v) {
   const RoutingAlg& alg = DYN_CAST(const RoutingAlg&, alg_);
   assert(&alg);

   (*v)[0] = alg.getBestSolScore();
   (*v)[1] = alg.getActualSolScore();
}
