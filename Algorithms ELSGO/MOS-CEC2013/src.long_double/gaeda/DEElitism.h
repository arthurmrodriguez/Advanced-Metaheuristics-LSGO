#ifndef DEELITISTRECOMBINATOR_H_
#define DEELITISTRECOMBINATOR_H_

#include "Recombinator.h"

class DEElitism : public Recombinator {
public:
  DEElitism();
  ~DEElitism();

  void recombine (const GAPopulation& old_pop, /*inout*/ GAPopulation& new_pop);
};

#endif /*DEELITISTRECOMBINATOR_H_*/
