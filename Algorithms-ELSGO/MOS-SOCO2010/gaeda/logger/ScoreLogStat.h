#ifndef SCORESTATVALUE_H_
#define SCORESTATVALUE_H_

#include "CollectionLogStat.h"

class ScoreLogStat : public CollectionLogStat
{
public:
	ScoreLogStat(const Algorithm& alg) : CollectionLogStat("score",alg) {}
	~ScoreLogStat(){}

	long double indValue(GAGenome& ind){return ind.score();}
};

#endif /*SCORESTATVALUE_H_*/
