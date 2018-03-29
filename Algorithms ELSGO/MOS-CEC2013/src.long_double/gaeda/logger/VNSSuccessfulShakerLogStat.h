/*
 * VNSBestShaker.h
 *
 *  Created on: Aug 4, 2011
 *      Author: santi
 */

#ifndef VNSSUCCESSFULSHAKER_H_
#define VNSSUCCESSFULSHAKER_H_

#include "LogStat.h"

class VNSSuccessfulShakerLogStat: public LogStat {
public:
  VNSSuccessfulShakerLogStat(const Algorithm& alg) : LogStat("VNSSuccessfulShakerLogStat",alg) {}
  virtual ~VNSSuccessfulShakerLogStat() {}

  virtual void update() {
    const VNS& vns = dynamic_cast<const VNS&>(alg_); assert(&vns);
    sprintf(message_,"succesfulshaker=%s ",vns.successfulShaker().c_str());

  }
};

#endif /* VNSBESTSHAKER_H_ */
