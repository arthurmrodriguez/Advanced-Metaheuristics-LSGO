#ifndef FASTFRACTAL_DD_H_
#define FASTFRACTAL_DD_H_
#include "../../../gaeda/genomes/GAGenome.h"
class FastFractal_dd
{
private:
  int dim_;
  GAGenome g_;
public:
	FastFractal_dd(int dim , GAGenome g);
	virtual ~FastFractal_dd();
};

#endif /*FASTFRACTAL_DD_H_*/
