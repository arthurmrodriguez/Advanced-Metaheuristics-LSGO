// $Header: /home/cvs/galib/ga/GAListGenome.h,v 1.3 2004/12/28 18:18:27 mwall Exp $
/* ----------------------------------------------------------------------------
  list.h
  mbwall 25feb95
  Copyright (c) 1995 Massachusetts Institute of Technology
                     all rights reserved

 DESCRIPTION:
  This header defines the interface for the list genome.
---------------------------------------------------------------------------- */
#ifndef GALISTGENOME_H
#define GALISTGENOME_H


/* INCLUDES */
#include "GAList.h"
#include "GAGenome.h"


/* ----------------------------------------------------------------------------
ListGenome
-------------------------------------------------------------------------------
---------------------------------------------------------------------------- */
template <class T>
class GAListGenome : public GAList<T>, public GAGenome {
public:
  GADeclareIdentity();

  static int DestructiveMutator(GAGenome &, double);
  static int SwapMutator(GAGenome &, double);
  static double NodeComparator(const GAGenome&, const GAGenome&);
  static int OnePointCrossover(const GAGenome&, const GAGenome&, 
			      GAGenome*, GAGenome*);
  static int PartialMatchCrossover(const GAGenome&, const GAGenome&, 
			      GAGenome*, GAGenome*);
  static int OrderCrossover(const GAGenome&, const GAGenome&, 
			      GAGenome*, GAGenome*);
  static int CycleCrossover(const GAGenome&, const GAGenome&, 
			      GAGenome*, GAGenome*);

public:
  GAListGenome(GAGenome::Evaluator f=0, void * u=(void *)0);
  GAListGenome(const GAListGenome<T> &);
  GAListGenome<T> & operator=(const GAGenome & orig)
    {copy(orig); return *this;}
  virtual ~GAListGenome();
  virtual GAGenome *clone(GAGenome::CloneMethod flag=CONTENTS) const;
  virtual void copy(const GAGenome &);

#ifdef GALIB_USE_STREAMS
  virtual int write (STD_OSTREAM &) const;
#endif

  virtual int equal(const GAGenome & c) const;

// Here we do inlined versions of the access members of the super class.  We
// do our own here so that we can set/unset the _evaluated flag appropriately.

  int destroy() { _evaluated = gaFalse; return GAList<T>::destroy(); }
  int swap(unsigned int i, unsigned int j)
    { _evaluated = gaFalse; return GAList<T>::swap(i,j); }
  T * remove() { _evaluated = gaFalse; return GAList<T>::remove(); }
  int insert(GAList<T> * t, GAListBASE::Location where=GAListBASE::AFTER)
    { _evaluated = gaFalse; return GAList<T>::insert(t, where); }
  int insert(const T & t, GAListBASE::Location where=GAListBASE::AFTER)
    { _evaluated = gaFalse; return GAList<T>::insert(t, where); }
};

#ifdef GALIB_USE_BORLAND_INST
#include "GAListGenome.cc"
#endif

#endif
