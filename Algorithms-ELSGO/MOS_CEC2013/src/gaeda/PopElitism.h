#ifndef ELITISTRECOMBINATOR_H_
#define ELITISTRECOMBINATOR_H_

#include "Recombinator.h"

class PopElitism : public Recombinator {
public:
	PopElitism();
	~PopElitism();

	void recombine (const GAPopulation& old_pop, /*inout*/ GAPopulation& new_pop);
};

#endif /*ELITISTRECOMBINATOR_H_*/
