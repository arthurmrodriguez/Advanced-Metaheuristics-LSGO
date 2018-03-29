#ifndef CLASSICELITISMRECOMBINATOR_H_
#define CLASSICELITISMRECOMBINATOR_H_

#include "Recombinator.h"

class ClassicElitism : public Recombinator
{
public:
	ClassicElitism();
	~ClassicElitism();

	void recombine (const GAPopulation& old_pop, /*inout*/ GAPopulation& new_pop);
};

#endif /*CLASSICELITISMRECOMBINATOR_H_*/
