#ifndef OPTDISTSTATVALUE_H_
#define OPTDISTSTATVALUE_H_

#include "CollectionLogStat.h"

class GAGenome;
class Algorithm;

/**
 * The object from this class computes the distance to optimum log stat.
 * Initially it received only the pointer to the optimum and the function invoqued
 * was always normalOptDistFunc that appeared on indValue. Later on, when a problem
 * with several global optimum appeared, it was decided that the user needed it to define
 * this function so a pointer to a distance function was incorporated to the constructor.
 * By default if this distance function is not defined, the objects from this class
 * execute the normal distance function.
 */
class OptDistLogStat : public CollectionLogStat {
  static GAGenome* opt_gen__ ;

  double (*dist_func_)(GAGenome&);

  static double normalOptDistFunc (GAGenome& gen);

public:
	OptDistLogStat(const Algorithm& alg,GAGenome* opt_gen, double (*dist_func)(GAGenome&) );
  ~OptDistLogStat();

  double indValue(GAGenome& ind);
};

#endif /*OPTDISTSTATVALUE_H_*/
