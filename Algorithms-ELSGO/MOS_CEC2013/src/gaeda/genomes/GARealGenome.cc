// $Header: /nfs/dsi/cvs/galib/ga/GARealGenome.C,v 1.1.1.1 1999/11/11 18:56:03 mbwall Exp $
/* ----------------------------------------------------------------------------
  real.C
  mbwall 11nov95
  Copyright (c) 1995-1996 Massachusetts Institute of Technology
                          all rights reserved

 DESCRIPTION:
   Source file for the real number specialization of the array genome.
---------------------------------------------------------------------------- */
#include "GARealGenome.h"
#include "GA1DArrayGenome.h"
#include <math.h>

template <> double GAAlleleSet<double>::allele () const {

   double value = 0.0;

   if (core->type == GAAllele::ENUMERATED) {
      value = core->a [GARandomInt (0, core->sz-1)];
   }
   else if (core->type == GAAllele::DISCRETIZED) {
      double n = (core->a [1] - core->a [0]) / core->a [2];
      double m = (double) n;

      if(core->lowerb == GAAllele::EXCLUSIVE) 
        m -= 1.0;
      if(core->upperb == GAAllele::EXCLUSIVE) 
        m -= 1.0;

      value = core->a [0] + GARandomDouble (0, (double) m) * core->a [2];

      if(core->lowerb == GAAllele::EXCLUSIVE) 
        value += core->a[2];
   }
   else {
      if (core->a [0] == core->a [1]           &&
          core->lowerb == GAAllele::EXCLUSIVE  &&
          core->upperb == GAAllele::EXCLUSIVE     ) {
         value = core->a [0];
      }
      else {
        do value = GARandomDouble (core->a [0], core->a [1]);
        while ((core->lowerb == GAAllele::EXCLUSIVE && value == core->a [0]) ||
               (core->upperb == GAAllele::EXCLUSIVE && value == core->a [1])    );
      }
   }
   return value;
}

// If someone asks for a discretized item that is beyond the bounds, give them
// one of the bounds.  If they ask for allele item when there is no
// discretization or enumeration, then error and return lower bound.

template <> double GAAlleleSet<double>::allele (unsigned int i) const {
   double value = 0.0;

   if (core->type == GAAllele::ENUMERATED) 
     value = core->a [i % core->sz];
   else if (core->type == GAAllele::DISCRETIZED) {
      double n = (core->a [1] - core->a [0]) / core->a [2];
      double m = (double) n; // What about bogus limits?

      if (core->lowerb == GAAllele::EXCLUSIVE) 
        m -= 1.0;
      if (core->upperb == GAAllele::EXCLUSIVE) 
        m -= 1.0;

      if (i > m) 
        i = (unsigned int) m;

      value = core->a [0] + i * core->a [2];

      if (core->lowerb == GAAllele::EXCLUSIVE) 
        value += core->a [2];
   }
   else {
      GAErr (GA_LOC, "GAAlleleSet", "allele", gaErrNoAlleleIndex);
      value = core->a [0];
   }
   return value;
}

template <>
int GA1DArrayAlleleGenome<double>::compsCompare (const GAGenome& g) const{
  const GA1DArrayAlleleGenome<double>& gen = dynamic_cast< const GA1DArrayAlleleGenome<double>& > (g);
  
  double allowed_range = fabs( gen.alleleset(0).upper() - gen.alleleset(0).lower() ) * 0.01;
    
  int ncomps=0;
  for (unsigned i=0; i<nx; i++){
    if ( fabs( gen.gene(i) - gene(i) ) <= allowed_range ) 
      ncomps++;
  }
  return ncomps;
}


// now the specialization of the genome itself.

template <> const char * 
GA1DArrayAlleleGenome<double>::className() const {return "GARealGenome";}
template <> int
GA1DArrayAlleleGenome<double>::classID() const {return GAID::DoubleGenome;}


#ifdef GALIB_USE_STREAMS
// The read specialization takes in each number and stuffs it into the array.
template <> int
GA1DArrayAlleleGenome<double>::read(STD_ISTREAM & is) {
  unsigned int i=0;
  double val;
  do{
    is >> val;
    if(!is.fail()) gene(i++, val);
  } while(!is.fail() && !is.eof() && i < nx);

  if(is.eof() && i < nx){
    GAErr(GA_LOC, className(), "read", gaErrUnexpectedEOF);
    is.clear(STD_IOS_BADBIT | is.rdstate());
    return 1;
  }
  return 0;
}

// No need to specialize the write method.
#endif

// force instantiations of this genome type.
//
// These must be included _after_ the specializations because some compilers
// get all wigged out about the declaration/specialization order.  Note that
// some compilers require a syntax different than others when forcing the 
// instantiation (i.e. GNU wants the 'template class', borland does not).
#ifndef GALIB_USE_AUTO_INST
#include "GAAllele.cc"
#include "GA1DArrayGenome.cc"

#if defined(__BORLANDC__)
#define GALIB_REALGENOME_TEMPLATE_PREFACE
#else
#define GALIB_REALGENOME_TEMPLATE_PREFACE template class
#endif

GALIB_REALGENOME_TEMPLATE_PREFACE GAAlleleSet<double>;
GALIB_REALGENOME_TEMPLATE_PREFACE GAAlleleSetCore<double>;
GALIB_REALGENOME_TEMPLATE_PREFACE GAAlleleSetArray<double>;

GALIB_REALGENOME_TEMPLATE_PREFACE GAArray<double>;
GALIB_REALGENOME_TEMPLATE_PREFACE GA1DArrayGenome<double>;
GALIB_REALGENOME_TEMPLATE_PREFACE GA1DArrayAlleleGenome<double>;

#endif
