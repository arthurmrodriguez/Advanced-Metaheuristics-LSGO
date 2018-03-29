/**
 * @file
 * @brief Evolution Strategy class hdr.
 *
 * Header file for the simple evolution strategy class
 */

#ifndef EVOLUTION_STRATEGY_H
#define EVOLUTION_STRATEGY_H

/* INCLUDES */
#include "GAGeneticAlgorithm.h"

template <class T> class GA1DArrayAlleleGenome;

class EvolutionStrategy : public GAGeneticAlgorithm {

 public:

  static GAParameterList& registerDefaultParameters(GAParameterList&);

  GADefineIdentity("EvolutionStrategy", GAID::ES);

  EvolutionStrategy(const GAGenome&);
  EvolutionStrategy(const GAPopulation&);
  EvolutionStrategy(const EvolutionStrategy&);

  virtual ~EvolutionStrategy();

  EvolutionStrategy& operator=(const EvolutionStrategy&);
  EvolutionStrategy& operator++();

  virtual void copy(const GAGeneticAlgorithm&);

  virtual void initialize();
  virtual void offspring (GAPopulation* offpop);

  virtual void step();

  virtual int setptr(const char* name, const void* value);
  virtual int get   (const char* name, void* value) const;

  virtual const GAPopulation& population(const GAPopulation&);
  virtual unsigned populationSize(unsigned value);

  virtual GAScalingScheme&   scaling (const GAScalingScheme&   s);
  virtual GASelectionScheme& selector(const GASelectionScheme& s);

  virtual void objectiveFunction(GAGenome::Evaluator f);
  virtual void objectiveData    (const GAEvalData& v);


  int IntermediateRecombination(const std::vector<GA1DArrayAlleleGenome<long double>*>& parents,
				GA1DArrayAlleleGenome<long double>& child) const;

 protected:
  GAPopulation* oldPop;
  unsigned _mu;
  unsigned _ro;
  unsigned _lambda;
  long double   _tau;
  long double   _tau0;

};

#ifdef GALIB_USE_STREAMS
inline STD_OSTREAM& operator<< (STD_OSTREAM& os, EvolutionStrategy& arg) {arg.write(os); return(os);}
inline STD_ISTREAM& operator>> (STD_ISTREAM& is, EvolutionStrategy& arg) {arg.read(is);  return(is);}
#endif

#endif
