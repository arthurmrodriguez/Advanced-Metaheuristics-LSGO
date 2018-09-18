// $Header: /nfs/dsi/cvs/galib/ga/GARealGenome.h,v 1.1.1.1 1999/11/11 18:56:03 mbwall Exp $
/* ----------------------------------------------------------------------------
  real.h
  mbwall 25feb95
  Copyright (c) 1995 Massachusetts Institute of Technology
                     all rights reserved

 DESCRIPTION:
  This header defines the specialization of the array genome of type double
for the real number genome.
---------------------------------------------------------------------------- */
#ifndef _ga_real_h_
#define _ga_real_h_

#include "GAAllele.h"
#include "GA1DArrayGenome.h"

typedef GAAlleleSet<double>      GARealAlleleSet;
typedef GAAlleleSetArray<double> GARealAlleleSetArray;

typedef GA1DArrayAlleleGenome<double> GARealGenome;

// in one (and only one) place in the code that uses the string genome, you 
// should define INSTANTIATE_STRING_GENOME in order to force the specialization
// for this genome.
#if defined(INSTANTIATE_REAL_GENOME)
#include "GARealGenome.cc"
#endif

#endif
