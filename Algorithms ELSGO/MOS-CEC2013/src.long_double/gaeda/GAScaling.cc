// $Header: /home/cvs/galib/ga/GAScaling.C,v 1.1.1.1 1999/11/11 18:56:03 mbwall Exp $
/* ----------------------------------------------------------------------------
  scaling.C
  mbwall 10aug94
  Copyright (c) 1995 Massachusetts Institute of Technology
                     all rights reserved

 DESCRIPTION:
  Class definitions for the scaling objects.  These objects translate (scale)
the raw objective scores of a population into 'fitness' scores that are used to
determine which genomes are fit for mating/selection.
  In all of these routines we do a check to be sure that we have enough space
to do our evaluations.  We don't need to do this so long as the population
objects that call us do the test themselves.  I'll leave the redundancy for
now.
---------------------------------------------------------------------------- */
#include <math.h>

#include "GAScaling.h"

#include "gaerror.h"
#include "GAPopulation.h"
#include "GAGeneticAlgorithm.h"
#include "genomes/GAGenome.h"

long double gaDefLinearScalingMultiplier   = 1.2;
long double gaDefSigmaTruncationMultiplier = 2.0;
long double gaDefPowerScalingFactor        = 1.0005;
long double gaDefSharingCutoff             = 1.0;




/* ----------------------------------------------------------------------------
NoScaling
---------------------------------------------------------------------------- */
// Assign the fitness scores to be the same as the objective scores for all of
// the individuals in the population.
void
GANoScaling::evaluate(const GAPopulation& p) {
  for(unsigned i=0; i<p.size(); i++)
    p.individual(i).fitness(p.individual(i).score());
}


void GAStandardScaling::evaluate(const GAPopulation& p) {

  p.statistics(gaTrue);

  long double best  = GAGenome::optCriterion() == GAGenome::MAXIMIZATION ? p.max() : p.min();
  long double worst = GAGenome::optCriterion() == GAGenome::MAXIMIZATION ? p.min() : p.max();

  for (unsigned i = 0; i < p.size(); i++)
    p.individual(i).fitness((best - p.individual(i).score()) / (best - worst));

}


/* ----------------------------------------------------------------------------
LinearScaling
---------------------------------------------------------------------------- */
// Scale the objective scores in the population.  We assume that the raw
// evaluation has already taken place (the calculation of the objective scores
// of the genomes in the population).  We must use a precision higher than that
// of the genome scores so that we don't lose any information.

void
GALinearScaling::evaluate(const GAPopulation & p) {
// Here we calculate the slope and intercept using the multiplier and objective
// score ranges...

  long double pmin = p.min();
  long double pmax = p.max();
  long double pave = p.ave();

  long double delta, a, b;
  if(pave == pmax){	// no scaling - population is all the same
    a = 1.0;
    b = 0.0;
  }
  else if(pmin > ((long double)c * pave - pmax)/((long double)c - 1.0)){
    delta = pmax - pave;
    a = ((long double)c - 1.0) * pave / delta;
    b = pave * (pmax - (long double)c * pave) / delta;
  }
  else{				// stretch to make min be 0
    delta = pave - pmin;
    a = pave / delta;
    b = -pmin * pave / delta;
  }

// and now we calculate the scaled scaled values.  Negative scores are not
// allowed with this kind of scaling, so check for negative values.  If we get
// a negative value, dump an error message then set all of the scores to 0.

  for(unsigned i=0; i<p.size(); i++){
    long double f = p.individual(i).score();
    if(f < 0.0){
      GAErr(GA_LOC, className(), "evaluate", gaErrNegFitness);
      for(unsigned ii=0; ii<p.size(); ii++)
	p.individual(ii).fitness(0.0);
      return;
    }
    f = f * a + b;
    if(f < 0) f = 0.0;	// truncate if necessary (only due to roundoff error)
    p.individual(i).fitness((long double)f);       // might lose information here!
  }
}


// Set the multiplier for this selection type.  The fmultiplier must be greater
// than 1.0 or else we'll get a divide by zero error in our scaling operations.
long double
GALinearScaling::multiplier(long double fm) {
  if(fm <= 1.0){
    GAErr(GA_LOC, className(), "multiplier", gaErrBadLinearScalingMult);
    return c;
  }
  return c = fm;
}




/* ----------------------------------------------------------------------------
SigmaTruncationScaling
---------------------------------------------------------------------------- */
// This is an implementation of the sigma truncation scaling method descibed in
// goldberg p 124.  If the scaled fitness is less than zero, we arbitrarily set
// it to zero (thus the truncation part of 'sigma truncation').
void
GASigmaTruncationScaling::evaluate(const GAPopulation & p) {
  for(unsigned i=0; i<p.size(); i++){
    long double f = (long double)(p.individual(i).score()) - (long double)(p.ave());
    f += (long double)c * (long double)(p.dev());
    if(f < 0) f = 0.0;
    p.individual(i).fitness((long double)f);       // might lose information here!
  }
}


// Set the multiplier for this selection type.  It should be greater than or
// equal to zero.
long double
GASigmaTruncationScaling::multiplier(long double fm) {
  if(fm < 0.0){
    GAErr(GA_LOC, className(), "multiplier", gaErrBadSigmaTruncationMult);
    return c;
  }
  return c = fm;
}




/* ----------------------------------------------------------------------------
PowerLawScaling
---------------------------------------------------------------------------- */
// This is an implementation of the most basic form of power scaling, where the
// fitness is a function of the objective score raised to some power.  Negative
// objective scores are not allowed.  If we get one, we post an error and set
// all of the fitness scores to zero.
void
GAPowerLawScaling::evaluate(const GAPopulation & p) {
  for(unsigned i=0; i<p.size(); i++){
    long double f = p.individual(i).score();
    if(f < 0.0){
      GAErr(GA_LOC, className(), "evaluate", gaErrPowerNegFitness);
      for(unsigned ii=0; ii<p.size(); ii++)
	p.individual(ii).fitness(0.0);
      return;
    }
    f = pow(f,(long double)k);
    p.individual(i).fitness((long double)f);       // might lose information here!
  }
}



/* ----------------------------------------------------------------------------
Sharing
---------------------------------------------------------------------------- */
// This is an implementation of speciation using the sharing method described
// by goldberg in his book.  This requires a user-defined distance function in
// order to work.  The distance function returns a value between
// 0 and 1 inclusive to tell us how similar two genomes are to each other.
// A value of 0 means that the two genomes are identical to each other, a
// value of 1 means they are completely different.
//   A single genome is identical to itself, so d(i,i) is 0.
//   If alpha is 1 then we don't use pow().
//   If we have a comparator to use, use it.  If not, use the comparator of
// each genome.
//   We can cut in half the number of calls to the sharing function by keeping
// one half of the ixj matrix.  This is because d(i,j) is the same as d(j,i).
// We cache the distances in an upper right triangular matrix stored as a
// series of long doubles.
//   If the population is maximizing then we derate by dividing.  If the
// population is minimizing then we derate by multiplying.  First we check to
// see if there is a GA using the population.  If there is, we use its min/max
// flag to determine whether or not we should be minimizing or maximizing.  If
// there is not GA with the population, then we use the population's sort order
// as the basis for whether to minimize or maximize.
// *** This could be done with n*n/2 instead of n*n, to reduce storage, but we
// can't reduce computation any more...
// *** probably should use the diversity built-in to the population...
void
GASharing::evaluate(const GAPopulation& p) {
  if(p.size() > N){
    delete [] d;
    N = p.size();
    d = new long double[N*N];
  }
  int n = p.size();

  int i, j;
  if(df) {
    for(i=0; i<n; i++){		// calculate and cache the distances
      d[i*n+i] = 0.0;		// each genome is same as itself
      for(j=i+1; j<n; j++)
	d[i*n+j] = d[j*n+i] = (*df)(p.individual(i), p.individual(j));
    }
  }
  else {
    for(i=0; i<n; i++){		// calculate and cache the distances
      d[i*n+i] = 0.0;		// each genome is same as itself
      for(j=i+1; j<n; j++)
	d[i*n+j] = d[j*n+i] = p.individual(i).compare(p.individual(j));
    }
  }

  for(i=0; i<n; i++){		// now derate the fitness of each genome
    long double sum = 0.0;
    for(j=0; j<n; j++) {
      if(d[i*n+j] < _sigma) {
	if(_alpha == 1)
	  sum += ((d[i*n+j] >= _sigma) ? 0.0 : 1.0 - d[i*n+j]/_sigma);
	else
	  sum += ((d[i*n+j]>=_sigma) ? 0.0 : 1.0-pow(d[i*n+j]/_sigma,_alpha));
      }
    }
    long double f;
    if(GAGenome::optCriterion() == GAGenome::MINIMIZATION)
      f = p.individual(i).score() * sum;
    else
      f = p.individual(i).score() / sum;
    p.individual(i).fitness((long double)f);       // might lose information here!
  }
}

void
GASharing::copy(const GAScalingScheme & arg){
  if(&arg == this) return;

  GAScalingScheme::copy(arg);
  const GASharing& s = DYN_CAST(const GASharing&, arg);
  delete [] d;
  _sigma = s._sigma;
  _alpha = s._alpha;
  df = s.df;
  N = s.N;
  d = new long double[N*N];
  memcpy(d, s.d, N*N*sizeof(long double));
}


// The cutoff for triangular sharing must always be greater than 0
long double
GASharing::sigma(long double c) {
  if(c <= 0.0){
    GAErr(GA_LOC, className(), "sigma", gaErrBadSharingCutoff);
    return _sigma;
  }
  return _sigma = c;
}
