#ifndef DE_H
#define DE_H

#include "GAGeneticAlgorithm.h"

class DE : public GAGeneticAlgorithm {
protected:
  GAPopulation* oldPop;
  long double        F_;
  long double        CR_;
  GAGenome::DECrossover  deCross_;

public:

  DE(const GAGenome&, long double F, long double CR);
  DE(const GAPopulation&, long double F, long double CR);
  DE(const DE&);

  virtual ~DE();

  DE&                         operator=(const DE&);

  virtual void                copy(const DE&);
  virtual void                initialize();

  virtual void                step();

  virtual const GAPopulation& population(const GAPopulation&);
  virtual unsigned int        populationSize(unsigned int value);

  virtual GAScalingScheme&    scaling(const GAScalingScheme & s);

  virtual GASelectionScheme&  selector(const GASelectionScheme& s);

  virtual void                objectiveFunction(GAGenome::Evaluator f);
  virtual void                objectiveData(const GAEvalData& v);

  virtual void                offspring (GAPopulation* offpop);
};



#ifdef GALIB_USE_STREAMS
inline STD_OSTREAM& operator<< (STD_OSTREAM& os, DE & arg)
{ arg.write(os); return(os); }
inline STD_ISTREAM& operator>> (STD_ISTREAM& is, DE & arg)
{ arg.read(is); return(is); }
#endif

#endif
