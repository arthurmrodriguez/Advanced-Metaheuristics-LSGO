#ifndef RECOMBINATOR_H_
#define RECOMBINATOR_H_

class GAPopulation;

class Recombinator {
public:
	Recombinator() {}
	virtual ~Recombinator() {}

	virtual void recombine (const GAPopulation& pop, /*inout*/ GAPopulation& children) = 0;
};

#endif /*RECOMBINATOR_H_*/
