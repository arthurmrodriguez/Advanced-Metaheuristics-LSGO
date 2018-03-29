#ifndef GENERATIONALRECOMBINATOR_H_
#define GENERATIONALRECOMBINATOR_H_

#include "Recombinator.h"

class GenerationalElitism : public Recombinator {
public:
	 GenerationalElitism();
	~GenerationalElitism();

	void recombine (const GAPopulation& old_pop, GAPopulation& new_pop);
};

#endif /*GENERATIONALRECOMBINATOR_H_*/
