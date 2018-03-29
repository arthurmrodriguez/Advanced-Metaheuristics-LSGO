#ifndef PERCENTAGEELITISM_H_
#define PERCENTAGEELITISM_H_

#include "Recombinator.h"

class PercentElitism: public Recombinator {

 private:
  double _percentage;

 public:
  PercentElitism(double percent);
  ~PercentElitism();

  void recombine (const GAPopulation& old_pop, /*inout*/ GAPopulation& new_pop);

};

#endif
