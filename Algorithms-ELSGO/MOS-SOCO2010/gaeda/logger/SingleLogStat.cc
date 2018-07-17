#include "SingleLogStat.h"

#include <stdio.h>

SingleLogStat::SingleLogStat(string name,const Algorithm& alg) : LogStat(name,alg){}

SingleLogStat::~SingleLogStat(){}

void SingleLogStat::update(){
  sprintf(message_,"%s=%.16LE ",name_.c_str(),computeValue());
}
