#ifndef ROUTINGALGBESTSOLNONPENSCORE
#define ROUTINGALGBESTSOLNONPENSCORE

#include "SingleLogStat.h"
#include "../RoutingAlg.h"

class RoutingAlgBestSolNonPenScore: public SingleLogStat {
public:
  RoutingAlgBestSolNonPenScore(const Algorithm& alg) : SingleLogStat("routingalgbestsolnonpenscore",alg) {}
  virtual ~RoutingAlgBestSolNonPenScore() {}

  virtual double computeValue() {
    const RoutingAlg& vns     = dynamic_cast<const RoutingAlg&>(alg_); assert(&vns);
    RoutingGenome&    bestsol = vns.bestSol();
    return bestsol.nonPenalizedScore();
  }
};

#endif
