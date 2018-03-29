#ifndef ROUTINGALGBESTSOLLOADV
#define ROUTINGALGBESTSOLLOADV

#include "SingleLogStat.h"
#include "../RoutingAlg.h"

class RoutingAlgBestSolLoadV: public SingleLogStat {
public:
  RoutingAlgBestSolLoadV(const Algorithm& alg) : SingleLogStat("routingalgbestsolloadv",alg) {}
  virtual ~RoutingAlgBestSolLoadV() {}

  virtual long double computeValue() {
    const RoutingAlg& vns     = dynamic_cast<const RoutingAlg&>(alg_); assert(&vns);
    RoutingGenome&    bestsol = vns.bestSol();
    return bestsol.loadV();
  }
};

#endif
