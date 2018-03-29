/*****************************************************************************
 * GAEDAlib: A C++ GA library with EDA and multiprocessor (MPI) support      *
 *                                                                           *
 * (C) 2005 Pedro Diaz (pdiaz@laurel.datsi.fi.upm.es)                        *
 *                                                                           *
 * GAEDAlib is distributed under the terms of the BSD software license       *
 *                                                                           *
 * GAEDAlib is heavily based on GAlib, a C++ GA library by Mathew Wall:      *
 * Copyright (c) 1995-1996 Massachusetts Institute of Technology (MIT)       *
 * Copyright (c) 1996-2000 Matthew Wall (author of GAlib)                    *
 *                                                                           *
 * Some portions of GAEDAlib's source code come from the GNU C++ compiler    *
 * library and therefore are covered under the terms of a different license, *
 * the GNU Public License.                                                   *
 *                                                                           *
 * You should have received a file named LICENSE along with this software.   *
 * This file contains more information about the licensing conditions of     *
 * GAEDAlib as well as the full text of each license involved.               *
 *                                                                           *
 * The file AUTHORS lists the people who have contributed (directly or       *
 * indirectly) to GAEDAlib                                                   *
 *****************************************************************************/

/**
 * @file
 * @brief GASimpleEDA class hdr.
 *
 * Header for the Pure EDA algorithm
 */

#ifndef GASIMPLEEDA_H
#define GASIMPLEEDA_H

/* INCLUDES */
#include "GAGeneticAlgorithm.h"

class GraphModel;

/* CLASS DECLARATION */

/**
 * @brief EDA (Estimation of Distribution Algorithm) implementation
 *
 *  GASimpleEDA provides a GA class that implements a pure Estimation of Distribution Algorithm.
 *  Bayesian networks are used to model discrete genomes, and Gaussian networks are used
 *  for continuous data.
 */
class GASimpleEDA : public GAGeneticAlgorithm {
	public:

		GADefineIdentity("GASimpleEDA", GAID::SimpleEDA);
		static GAParameterList& registerDefaultParameters(GAParameterList&);
		GASimpleEDA(GAGenome&);
		GASimpleEDA(const GAPopulation&);
		GASimpleEDA(const GASimpleEDA&);
		GASimpleEDA& operator=(const GASimpleEDA&);
		virtual ~GASimpleEDA();
		virtual void copy(const GAGeneticAlgorithm &ga);
		virtual void initialize();
		virtual void offspring(GAPopulation* offpop);
		virtual void step();
		GASimpleEDA & operator++() { step(); return *this; }

		virtual int setptr(const char* name, const void* value);
		virtual int get(const char* name, void* value) const;

		GABoolean elitist() const {return el;}
		GABoolean elitist(GABoolean flag) {params.set(gaNelitism, (int)flag); return el=flag;}

		virtual const GAPopulation& population() const {return *pop;}
		virtual const GAPopulation& population(const GAPopulation&);
		virtual unsigned int populationSize() const {return pop->size();}
		virtual unsigned int populationSize(unsigned int value);
		virtual GAScalingScheme& scaling() const {return pop->scaling();}
		virtual GAScalingScheme& scaling(const GAScalingScheme & s)
			{mOldPop->scaling(s); return GAGeneticAlgorithm::scaling(s);}
		virtual GASelectionScheme& selector() const {return pop->selector(); }
		virtual GASelectionScheme& selector(const GASelectionScheme& s)
			{mOldPop->selector(s); return GAGeneticAlgorithm::selector(s);}
		virtual void objectiveFunction(GAGenome::Evaluator f);
		virtual void objectiveData(const GAEvalData& v);

	private:
		GAPopulation *mOldPop;		///< Current and old populations
		GABoolean el;				///< Are we elitist?
		long double mSelPerc;			///< Selection percentage

		/// The Bayesian network estimated from the selected individuals.
		GAGraphModel* network;
};



#ifdef GALIB_USE_STREAMS
inline STD_OSTREAM& operator<< (STD_OSTREAM& os, GASimpleEDA & arg)
{ arg.write(os); return(os); }
inline STD_ISTREAM& operator>> (STD_ISTREAM& is, GASimpleEDA & arg)
{ arg.read(is); return(is); }
#endif

#endif
