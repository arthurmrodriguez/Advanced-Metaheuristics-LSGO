#ifndef ROUTINGALGBESTSOLDELIVERYTIMELOGSTAT
#define ROUTINGALGBESTSOLDELIVERYTIMELOGSTAT

#include "LogStat.h"
#include "../RoutingAlg.h"

class RoutingAlgBestSolDeliveryTimeLogStat: public LogStat {
public:
  RoutingAlgBestSolDeliveryTimeLogStat(const Algorithm& alg) : LogStat("RoutingAlgBestSolDeliveryTimeLogStat",alg) {}
  virtual ~RoutingAlgBestSolDeliveryTimeLogStat() {}

  virtual void update() {
    const RoutingAlg& vns = dynamic_cast<const RoutingAlg&>(alg_); assert(&vns);
    RoutingGenome& bestsol = vns.bestSol();
    sprintf(message_,"bestsoldeliverytime=%f ", bestsol.deliveryDelay());
  }
};

#endif
