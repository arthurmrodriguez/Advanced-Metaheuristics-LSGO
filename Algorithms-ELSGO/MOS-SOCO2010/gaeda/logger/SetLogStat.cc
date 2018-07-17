#include "SetLogStat.h"

#include <iomanip>
#include <math.h>
#include <string.h>
#include <sstream>
#include <stdexcept>

SetLogStat::SetLogStat(string name, const Algorithm& alg): LogStat(name,alg) {}

SetLogStat::~SetLogStat(){}

void SetLogStat::update(){

  vector<long double> stat_vector(_setSize);
  stringstream str;

  computeValues(&stat_vector);

  str << showbase << showpoint << showpos << scientific << setw(4) << setprecision(4) << right ;
  str << name_.c_str() << "=[";

  vector<long double>::iterator iter = stat_vector.begin();
  while(iter != stat_vector.end()) {
    if (iter == stat_vector.end () - 1)
      str << *iter;
    else
      str << *iter << ",";
    iter++;
  }

  str << "] ";

  if ((int) str.str().size() >= MESSAGE_SIZE)
    throw runtime_error ("Error: MESSAGE_SIZE exceeded.");

  strncpy (message_, str.str().c_str(), MESSAGE_SIZE);

}
