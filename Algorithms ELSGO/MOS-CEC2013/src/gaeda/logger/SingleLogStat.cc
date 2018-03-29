#include <stdio.h>

#include "SingleLogStat.h"

SingleLogStat::SingleLogStat(string name,const Algorithm& alg) : LogStat(name,alg){}

SingleLogStat::~SingleLogStat(){}

void SingleLogStat::update(){
  sprintf(message_,"%s=%.6E ",name_.c_str(),computeValue());
}
