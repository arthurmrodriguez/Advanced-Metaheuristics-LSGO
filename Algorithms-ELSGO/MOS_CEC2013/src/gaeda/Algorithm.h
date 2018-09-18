#ifndef ALGORITHM_H_
#define ALGORITHM_H_

#include "gaid.h"

class GAGenome;
class GAStatistics;
class GAPopulation;

class Algorithm : public GAID {

protected:
  virtual void evolve () = 0;

public:

  GADefineIdentity ("Algorithm", GAID::Algorithm);

  Algorithm();
  virtual ~Algorithm();

  virtual void                initialize  () ;
  virtual void                run         (); // Hack to make sure that initialize is always called
  virtual GAGenome*           best        () const = 0;
  virtual const GAStatistics& statistics  () const = 0;
  virtual const GAPopulation& population  () const = 0;
  virtual void                step        () = 0;

};

#endif /*ALGORITHM_H_*/
