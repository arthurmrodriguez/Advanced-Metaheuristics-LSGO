#ifndef INCRECOMBINATOR_H_
#define INCRECOMBINATOR_H_

#include "Recombinator.h"

class IncElitism : public Recombinator {
public:
	IncElitism();
	~IncElitism();

	void recombine (const GAPopulation& old_pop, /*inout*/ GAPopulation& new_pop);
};

#endif /*INCRECOMBINATOR_H_*/
