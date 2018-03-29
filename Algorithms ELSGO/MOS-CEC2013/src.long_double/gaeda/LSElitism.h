#ifndef LSELITISMRECOMBINATOR_H_
#define LSELITISMRECOMBINATOR_H_

#include "Recombinator.h"

class LSElitism : public Recombinator
{
public:
	LSElitism();
	~LSElitism();

	void recombine (const GAPopulation& old_pop, GAPopulation& new_pop);
};

#endif /*LSELITISMRECOMBINATOR_H_*/
