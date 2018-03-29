#ifndef VNSSHAKERPOS_H_
#define VNSSHAKERPOS_H_

#include "LogStat.h"

class VNSShakerPosLogStat: public LogStat {
public:
  VNSShakerPosLogStat(const Algorithm& alg) : LogStat("VNSShakerPosLogStat",alg) {}
  virtual ~VNSShakerPosLogStat() {}

  virtual void update() {
    const VNS& vns = dynamic_cast<const VNS&>(alg_); assert(&vns);
    sprintf(message_,"shakerpos=%d ",vns.shakerPos());

  }
};

#endif
