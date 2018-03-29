// $Header: /home/cvs/galib/ga/GAStringGenome.h,v 1.2 2004/12/28 18:18:27 mwall Exp $
/* ----------------------------------------------------------------------------
  string.h
  mbwall 25feb95
  Copyright (c) 1995 Massachusetts Institute of Technology

 DESCRIPTION:
  This header defines the interface for the string genome.
---------------------------------------------------------------------------- */
#ifndef GASTRINGGENOME_H
#define GASTRINGGENOME_H

/* INCLUDES */
#include "GAAllele.h"
#include "GA1DArrayGenome.h"

typedef GAAlleleSet<char> GAStringAlleleSet;
typedef GAAlleleSet<char> GACharacterAlleleSet;
typedef GAAlleleSetArray<char> GAStringAlleleSetArray;

typedef GA1DArrayAlleleGenome<char> GAStringGenome;

// in one (and only one) place in the code that uses the string genome, you 
// should define INSTANTIATE_STRING_GENOME in order to force the specialization
// for this genome.
#if defined(INSTANTIATE_STRING_GENOME)
#include "GAStringGenome.cc"
#endif

#endif
