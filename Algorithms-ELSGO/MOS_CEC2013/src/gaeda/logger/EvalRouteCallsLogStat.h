#ifndef EVALROUTECALLSLOGSTAT_H_
#define EVALROUTECALLSLOGSTAT_H_

#include "LogStat.h"

class EvalRouteCallsLogStat: public LogStat {
public:
  EvalRouteCallsLogStat(const Algorithm& alg) : LogStat("EvalRouteCallsLogStat",alg) {}
  virtual ~EvalRouteCallsLogStat() {}

  virtual void update() {
    const RoutingAlg& ralg = dynamic_cast<const RoutingAlg&>(alg_); assert(&ralg);
    sprintf(message_,"evalroutecalls=%d ",ralg.evalRouteCalls());

  }
};

#endif /* VNSBESTSHAKER_H_ */
