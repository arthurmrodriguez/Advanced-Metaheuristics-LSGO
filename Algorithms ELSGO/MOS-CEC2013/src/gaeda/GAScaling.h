// $Header: /home/cvs/galib/ga/GAScaling.h,v 1.2 2004/12/28 14:38:44 mwall Exp $
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
 * @brief GAScaling class hdr.
 *
 * Header file for the class responsible for doing the
 * fitness scaling
 */


#ifndef GASCALING_H
#define GASCALING_H

/* INCLUDES */
#include <string.h>

#include "gaid.h"
#include "gaconfig.h"
#include "gatypes.h"
#include "genomes/GAGenome.h"

class GAPopulation;

extern double gaDefLinearScalingMultiplier;
extern double gaDefSigmaTruncationMultiplier;
extern double gaDefPowerScalingFactor;
extern double gaDefSharingCutoff;


/**
 * @brief Class used to do scaling on a population
 *
 * This class is used to do scaling on a population.  When
 * a genome is evaluated, the user's objective function provides an overall
 * rating for each genome.  The GA is concerned with fitness, not objective
 * score (unless you do no scaling).  So this object is the container for the
 * scaled values.
 *
 * Examples of scaling include linear scaling, sigma truncation, and power law.
 * Goldberg's sharing function is also a type of scaling, and it is implemented
 * here as a unique class.

 * The scaling class is designed to be used with a population.  It is stupid -
 * it does know how to update itself, but it must be told when.
 */



/* ----------------------------------------------------------------------------
Scaling

  The scaling object is used to scale the objective scores of a population to
avoid clustering and premature convergence (among other things).  See golberg
for more about the theory.  This is basically just a container for any data
that the scaling object might need to do its thing.  The simplest scalings
don't store any data.
---------------------------------------------------------------------------- */
class GAScalingScheme : public GAID {
public:
  GADefineIdentity("GAScalingScheme", GAID::Scaling);

  GAScalingScheme() {}
  GAScalingScheme(const GAScalingScheme& s) { copy(s); }
  GAScalingScheme& operator=(const GAScalingScheme& s) {copy(s); return *this;}
  virtual ~GAScalingScheme() { }
  virtual GAScalingScheme * clone() const=0;
  virtual void copy(const GAScalingScheme &) {}
  virtual void evaluate(const GAPopulation & p)=0;
};


/* ----------------------------------------------------------------------------
NoScaling
---------------------------------------------------------------------------- */
class GANoScaling : public GAScalingScheme {
public:
  GADefineIdentity("GANoScaling", GAID::NoScaling);

  GANoScaling() : GAScalingScheme() {}
  GANoScaling(const GANoScaling &) : GAScalingScheme() {}
  GANoScaling& operator=(const GAScalingScheme&){ return *this; }
  virtual ~GANoScaling(){}
  virtual GAScalingScheme * clone() const {return new GANoScaling(*this);}
  virtual void evaluate(const GAPopulation & p);
};


class GAStandardScaling : public GAScalingScheme {
public:
  GADefineIdentity("GAStandardScaling", GAID::StandardScaling);

  GAStandardScaling() : GAScalingScheme() {}
  GAStandardScaling(const GAStandardScaling &) : GAScalingScheme() {}
  GAStandardScaling& operator=(const GAScalingScheme&){ return *this; }
  virtual ~GAStandardScaling(){}
  virtual GAScalingScheme * clone() const {return new GAStandardScaling(*this);}
  virtual void evaluate(const GAPopulation & p);
};


/* ----------------------------------------------------------------------------
LinearScaling

  This scaling object does linear scaling as described in goldberg pp 122-124.
---------------------------------------------------------------------------- */
class GALinearScaling : public GAScalingScheme {
public:
  GADefineIdentity("GALinearScaling", GAID::LinearScaling);

  GALinearScaling(double fm=gaDefLinearScalingMultiplier) {multiplier(fm);}
  GALinearScaling(const GALinearScaling & arg) {copy(arg);}
  GALinearScaling & operator=(const GAScalingScheme & arg)
    {copy(arg); return(*this);}
  virtual ~GALinearScaling(){}
  virtual GAScalingScheme * clone() const {return new GALinearScaling(*this);}
  virtual void evaluate(const GAPopulation & p);
  virtual void copy(const GAScalingScheme & arg){
    if(&arg != this){
      GAScalingScheme::copy(arg);
      c=(DYN_CAST(const GALinearScaling&,arg)).c;
    }
  }

  double multiplier(double fm);
  double multiplier() const {return c;}

protected:
  double c;			// linear scaling multiplier
};



/* ----------------------------------------------------------------------------
SigmaTruncationScaling

  This scaling object does sigma truncation as defined in goldberg p124.
---------------------------------------------------------------------------- */
class GASigmaTruncationScaling : public GAScalingScheme {
public:
  GADefineIdentity("GASigmaTruncationScaling", GAID::SigmaTruncationScaling);

  GASigmaTruncationScaling(double m=gaDefSigmaTruncationMultiplier)
    {multiplier(m);}
  GASigmaTruncationScaling(const GASigmaTruncationScaling & arg){copy(arg);}
  GASigmaTruncationScaling & operator=(const GAScalingScheme & arg)
    {copy(arg); return(*this);}
  virtual ~GASigmaTruncationScaling(){}
  virtual GAScalingScheme * clone() const
    {return new GASigmaTruncationScaling(*this);}
  virtual void evaluate(const GAPopulation & p);
  virtual void copy(const GAScalingScheme & arg){
    if(&arg != this){
      GAScalingScheme::copy(arg);
      c=(DYN_CAST(const GASigmaTruncationScaling&,arg)).c;
    }
  }

  double multiplier(double fm);
  double multiplier() const {return c;}

protected:
  double c;			// std deviation multiplier
};



/* ----------------------------------------------------------------------------
PowerLawScaling

  This scaling object does power law scaling as defined in goldberg p124.
---------------------------------------------------------------------------- */
class GAPowerLawScaling : public GAScalingScheme {
public:
  GADefineIdentity("GAPowerLawScaling", GAID::PowerLawScaling);

  GAPowerLawScaling(double f=gaDefPowerScalingFactor) {k = f;}
  GAPowerLawScaling(const GAPowerLawScaling & arg) {copy(arg);}
  GAPowerLawScaling & operator=(const GAScalingScheme & arg)
    {copy(arg); return(*this);}
  virtual ~GAPowerLawScaling(){}
  virtual GAScalingScheme * clone() const
    {return new GAPowerLawScaling(*this);}
  virtual void evaluate(const GAPopulation & p);
  virtual void copy(const GAScalingScheme & arg){
    if(&arg != this){
      GAScalingScheme::copy(arg);
      k=(DYN_CAST(const GAPowerLawScaling&,arg)).k;
    }
  }

  double power(double p){return k=p;}
  double power() const {return k;}

protected:
  double k;			// power scaling factor
};



/* ----------------------------------------------------------------------------
Sharing

  This scaling object does sharing as described in goldberg p 192.  This
implementation does triangular sharing with the (optional) alpha parameter for
changing the curvature of the sharing range and the (required) sigma parameter
for controlling the range of influence.  If you want a different type of
sharing function, then derive from this class and define a new (virtual)
evaluate method and add whatever member data you need to specify the shape
of your sharing function.
  We use the distance function to scale the objective scores of the
genomes so that if there are many similar genomes in a population
their scores will be decreased, thus giving sub-species a greater chance to
reproduce.
  This sharing function is defined as follows:

                     /
                    |  1 - (d(i,j)/sigma) ^ alpha       d(i,j) < sigma
        s(d(i,j)) = |
                    |  0                                d(i,j) >= sigma
                     \

  where d is the distance between any two individuals and is defined such
that d(i,j) = 0 means that individuals i and j are identical.  d() has no
upper bound, but it is never negative.
  The alpha controls the shape of the sharing function.  When alpha=1 then
the 'curve' is a straight line.  If alpha is < 1 then the curves are concave,
if alpha is > 1 then the curves are convex.
  If you decide to use this type of sharing, be careful of the sigma value
that you use.  It can make a HUGE difference, depending upon the objective.
  Distance functions are independent of the sharing functions (as defined in
Goldberg, that is).
  A similarity (distance) function is used with the sharing object.  It is a
type of speciation (similar in functionality to DeJong crowding, but this uses
fitness scaling rather than replacement strategy to affect the speciation).
If the genomes are identical, the similarity function should return a value of
0.0, if completely different then return a value of 1.0.  0 means less
diversity means all genomes are the same.
  You can specify whether the scaling should be maximize or minimize based.
If you are maximizing, the scaling will divide the raw score by the scaling
factor.  If you are minimizing, the scaling will multiply the score by the
scaling factor.  (By definition, the scaling factor will always be >= 1.0)
  By default, the scaling object uses the max/min settings that it contains
(so you can set the scaling independently of the GA).  If the scaling's min/max
was not set, then it tries to use the min/max settings in the GA that owns
the population to which the scaling object is attached.  If there is no GA,
then it bases its min/max on the sort order of the population.
  You can set the minimaxi to:

    GA::MINIMIZE - scale by multiplying the raw scores
    GA::MAXIMIZE - scale by dividing the raw scores
    0            - minimize or maximize based upon the GA's settings

*** This should be called TriangularSharing rather than simply Sharing.
---------------------------------------------------------------------------- */
class GASharing : public GAScalingScheme {
public:
  GADefineIdentity("GASharing", GAID::Sharing);

  GASharing(GAGenome::Comparator func,
	    double cut=gaDefSharingCutoff, double a=1.0)
    { N=0; d=(double*)0; df=func; _sigma = cut; _alpha = a; }
  GASharing(double cut=gaDefSharingCutoff, double a=1.0)
    { N=0; d=(double*)0; df=0; _sigma = cut; _alpha = a; }
  GASharing(const GASharing & arg) { N=0; d=(double*)0; copy(arg); }
  GASharing & operator=(const GAScalingScheme & arg){copy(arg); return(*this);}
  virtual ~GASharing(){ delete [] d;}
  virtual GAScalingScheme * clone() const {return new GASharing(*this);}
  virtual void copy(const GAScalingScheme & arg);
  virtual void evaluate(const GAPopulation & p);

  GAGenome::Comparator distanceFunction(GAGenome::Comparator f){return df=f;}
  GAGenome::Comparator distanceFunction() const {return df;}

  double sigma(double);
  double sigma() const { return _sigma; }

  double alpha(double c) { return _alpha = c; }
  double alpha() const { return _alpha; }

protected:
  GAGenome::Comparator df;		// the user-defined distance function
  unsigned int N;			// how many do we have? (n of n-by-n)
  double *d;				// the distances for each genome pair
  double _sigma;				// absolute cutoff from central point
  double _alpha;				// controls the curvature of sharing f
};


#endif
