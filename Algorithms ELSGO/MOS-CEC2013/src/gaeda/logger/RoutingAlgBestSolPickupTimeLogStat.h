#ifndef ROUTINGALGBESTSOLPICKUPTIMELOGSTAT_H_
#define ROUTINGALGBESTSOLPICKUPTIMELOGSTAT_H_

#include "LogStat.h"
#include "../RoutingAlg.h"

class RoutingAlgBestSolPickupTimeLogStat: public LogStat {
public:
  RoutingAlgBestSolPickupTimeLogStat(const Algorithm& alg) : LogStat("RoutingAlgBestSolPickupTimeLogStat",alg) {}
  virtual ~RoutingAlgBestSolPickupTimeLogStat() {}

  virtual void update() {
    const RoutingAlg& alg = dynamic_cast<const RoutingAlg&>(alg_); assert(&alg);
    RoutingGenome& bestsol = alg.bestSol();
    sprintf(message_,"bestsolpickuptime=%f ", bestsol.pickupDelay());
  }
};

#endif
