#ifndef ROUTINGALGBESTSOLRIDEV
#define ROUTINGALGBESTSOLRIDEV

#include "SingleLogStat.h"
#include "../RoutingAlg.h"

class RoutingAlgBestSolRideV: public SingleLogStat {
public:
  RoutingAlgBestSolRideV(const Algorithm& alg) : SingleLogStat("routingalgbestsolridev",alg) {}
  virtual ~RoutingAlgBestSolRideV() {}

  virtual long double computeValue() {
    const RoutingAlg& vns     = dynamic_cast<const RoutingAlg&>(alg_); assert(&vns);
    RoutingGenome&    bestsol = vns.bestSol();
    return bestsol.rideV();
  }
};

#endif
