/**
 * @file
 * @brief GASimpleGA class hdr.
 *
 * Header file for the simple genetic algorithm class
 */


#ifndef GASIMPLEGA_H
#define GASIMPLEGA_H

/* INCLUDES */
#include "GAGeneticAlgorithm.h"


/**
 * @brief Simple genetic algorithm with non-overlapping populations
 *
 * This genetic algorithm is the 'simple' genetic algorithm that Goldberg
 * describes in his book. It uses non-overlapping populations. When you
 * create a simple genetic algorithm, you must specify either an individual
 * or a population of individuals. The new genetic algorithm will clone
 * the individual(s) that you specify to make its own population. You can
 * change most of the genetic algorithm behaviors after creation and during
 * the course of the evolution.
 *
 * The simple genetic algorithm creates an initial population by cloning
 * the individual or population you pass when you create it. Each
 * generation the algorithm creates an entirely new population of
 * individuals by selecting from the previous population then mating to
 * produce the new offspring for the new population. This process continues
 * until the stopping criteria are met (determined by the terminator).
 *
 * Elitism is optional. By default, elitism is on, meaning that the best
 * individual from each generation is carried over to the next generation.
 * To turn off elitism, pass gaFalse to the elitist member function.
 *
 * The score frequency for this genetic algorithm defaults to 1 (it records
 * the best-of-generation every generation). The default scaling is Linear,
 * the default selection is RouletteWheel.
 */
class GASimpleGA : public GAGeneticAlgorithm {
protected:
  GAPopulation* oldPop;
public:
    static GAParameterList& registerDefaultParameters(GAParameterList&);

    GADefineIdentity("GASimpleGA", GAID::SimpleGA);

		GASimpleGA(const GAGenome&);
		GASimpleGA(const GAPopulation&);
		GASimpleGA(const GASimpleGA&);

		virtual ~GASimpleGA();

    GASimpleGA&                 operator=(const GASimpleGA&);
    GASimpleGA&                 operator++();

		virtual void                copy(const GAGeneticAlgorithm&);
		virtual void                initialize();
		virtual void offspring(GAPopulation* offpop);

		virtual void                step();

		virtual int                 setptr(const char* name, const void* value);
		virtual int                 get(const char* name, void* value) const;

		virtual const GAPopulation& population(const GAPopulation&);
		virtual unsigned int                 populationSize(unsigned int value);

		virtual GAScalingScheme&    scaling(const GAScalingScheme & s);

		virtual GASelectionScheme&  selector(const GASelectionScheme& s);

		virtual void                objectiveFunction(GAGenome::Evaluator f);
		virtual void                objectiveData(const GAEvalData& v);
};



#ifdef GALIB_USE_STREAMS
inline STD_OSTREAM& operator<< (STD_OSTREAM& os, GASimpleGA & arg)
{ arg.write(os); return(os); }
inline STD_ISTREAM& operator>> (STD_ISTREAM& is, GASimpleGA & arg)
{ arg.read(is); return(is); }
#endif

#endif
