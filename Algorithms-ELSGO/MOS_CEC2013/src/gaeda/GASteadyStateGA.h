// $Header: /home/cvs/galib/ga/GASStateGA.h,v 1.2 2004/12/28 00:12:12 mwall Exp $
/* ----------------------------------------------------------------------------
  gasteadystate.h
  mbwall 28jul94
  Copyright (c) 1995 Massachusetts Institute of Technology
                     all rights reserved

  Header file for the steady-state genetic algorithm class.
---------------------------------------------------------------------------- */
#ifndef GASSTATEGA_H
#define GASSTATEGA_H

/* INCLUDES */
#include "GAGeneticAlgorithm.h"

class GASteadyStateGA : public GAGeneticAlgorithm {
public:
  GADefineIdentity("GASteadyStateGA", GAID::SteadyStateGA);

  static GAParameterList& registerDefaultParameters(GAParameterList&);

public:
  GASteadyStateGA(const GAGenome&);
  GASteadyStateGA(const GAPopulation&);
  GASteadyStateGA(const GASteadyStateGA&);
  GASteadyStateGA& operator=(const GASteadyStateGA&);
  virtual ~GASteadyStateGA();
  virtual void copy(const GAGeneticAlgorithm&);

  virtual void initialize();
  virtual void offspring(GAPopulation* offpop);
  virtual void step();
  GASteadyStateGA & operator++() { step(); return *this; }

  virtual int setptr(const char* name, const void* value);
  virtual int get(const char* name, void* value) const;

  virtual const GAPopulation& population() const {return *pop;}
  virtual const GAPopulation& population(const GAPopulation&);
  virtual unsigned int populationSize() const {return pop->size();}
  virtual unsigned int populationSize(unsigned int n);
  virtual GAScalingScheme& scaling() const {return pop->scaling();}
  virtual GAScalingScheme& scaling(const GAScalingScheme & s)
    { /* tmpPop->scaling(s); */ return GAGeneticAlgorithm::scaling(s); }
  virtual GASelectionScheme& selector() const {return pop->selector(); }
  virtual GASelectionScheme& selector(const GASelectionScheme& s)
    { /* tmpPop->selector(s); */ return GAGeneticAlgorithm::selector(s); }
  virtual void objectiveFunction(GAGenome::Evaluator f);
  virtual void objectiveData(const GAEvalData& v);

  double pReplacement() const { return pRepl; }
  double pReplacement(double p);
  int nReplacement() const { return nRepl; }
  int nReplacement(unsigned int n);

protected:
  GAPopulation *tmpPop;		// temporary population for replacements
  double pRepl;			// percentage of population to replace each gen
  unsigned int nRepl;		// how many of each population to replace
  short which;			// 0 if prepl, 1 if nrepl
};



#ifdef GALIB_USE_STREAMS
inline STD_OSTREAM& operator<< (STD_OSTREAM& os, GASteadyStateGA & arg)
{ arg.write(os); return(os); }
inline STD_ISTREAM& operator>> (STD_ISTREAM& is, GASteadyStateGA & arg)
{ arg.read(is); return(is); }
#endif

#endif
