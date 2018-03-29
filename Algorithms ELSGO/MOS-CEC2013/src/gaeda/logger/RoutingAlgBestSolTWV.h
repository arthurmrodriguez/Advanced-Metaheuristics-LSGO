#ifndef ROUTINGALGBESTSOLTWV
#define ROUTINGALGBESTSOLTWV

#include "SingleLogStat.h"
#include "../RoutingAlg.h"

class RoutingAlgBestSolTWV: public SingleLogStat {
public:
  RoutingAlgBestSolTWV(const Algorithm& alg) : SingleLogStat("routingalgbestsoltwv",alg) {}
  virtual ~RoutingAlgBestSolTWV() {}

  virtual double computeValue() {
    const RoutingAlg& vns     = dynamic_cast<const RoutingAlg&>(alg_); assert(&vns);
    RoutingGenome&    bestsol = vns.bestSol();
    return bestsol.TWV();
  }
};

#endif
