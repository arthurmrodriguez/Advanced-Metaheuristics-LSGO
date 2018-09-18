// $Header: /home/cvs/galib/ga/GAStringGenome.C,v 1.5 2004/12/28 22:17:30 mwall Exp $
/* ----------------------------------------------------------------------------
  string.C
  mbwall 21mar95
  Copyright (c) 1995-1996 Massachusetts Institute of Technology
                          all rights reserved

 DESCRIPTION:
   Source file for the string specialization of the array genome.
---------------------------------------------------------------------------- */
#include "GAStringGenome.h"

template <> const char * 
GA1DArrayAlleleGenome<char>::className() const {return "GAStringGenome";}
template <> int
GA1DArrayAlleleGenome<char>::classID() const {return GAID::StringGenome;}

#ifdef GALIB_USE_STREAMS
// The read specialization takes in each character whether it is whitespace or
// not and stuffs it into the genome.  This is unlike the default array read.
template <> int
GA1DArrayAlleleGenome<char>::read(STD_ISTREAM & is)
{
  unsigned int i=0;
  char c;
  do{
    is.get(c);
    if(!is.fail()) gene(i++, c);
  } while(!is.fail() && !is.eof() && i < nx);

  if(is.eof() && i < nx){
    GAErr(GA_LOC, className(), "read", gaErrUnexpectedEOF);
    is.clear(STD_IOS_BADBIT | is.rdstate());
    return 1;
  }
  return 0;
}

// Unlike the base array genome, here when we write out we don't put any
// whitespace between genes.  No newline at end of it all.
template <> int
GA1DArrayAlleleGenome<char>::write(STD_OSTREAM & os) const
{
  for(unsigned int i=0; i<nx; i++)
    os << gene(i);
  return 0;
}
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
#define GALIB_STRINGGENOME_TEMPLATE_PREFACE
#else
#define GALIB_STRINGGENOME_TEMPLATE_PREFACE template class
#endif

GALIB_STRINGGENOME_TEMPLATE_PREFACE GAAlleleSet<char>;
GALIB_STRINGGENOME_TEMPLATE_PREFACE GAAlleleSetCore<char>;
GALIB_STRINGGENOME_TEMPLATE_PREFACE GAAlleleSetArray<char>;

GALIB_STRINGGENOME_TEMPLATE_PREFACE GAArray<char>;
GALIB_STRINGGENOME_TEMPLATE_PREFACE GA1DArrayGenome<char>;
GALIB_STRINGGENOME_TEMPLATE_PREFACE GA1DArrayAlleleGenome<char>;

#endif
