// $Header: /home/cvs/galib/ga/GATreeGenome.h,v 1.3 2004/12/28 18:18:27 mwall Exp $
/* ----------------------------------------------------------------------------
  tree.h
  mbwall 25feb95
  Copyright (c) 1995 Massachusetts Institute of Technology
                     all rights reserved

 DESCRIPTION:
  This header defines the interface for the tree genome.
---------------------------------------------------------------------------- */
#ifndef GATREEGENOME_H
#define GATREEGENOME_H

/* INCLUDES */
#include "GATree.h"
#include "GAGenome.h"


/* ----------------------------------------------------------------------------
TreeGenome
-------------------------------------------------------------------------------
  Beware that the tree genome can grow unbounded - there is no size limit
on the tree, so if you have an objective function that encourages size the tree
will grow until you run out of memory.
---------------------------------------------------------------------------- */
template <class T>
class GATreeGenome : public GATree<T>, public GAGenome {
public:
  GADeclareIdentity();

  static int DestructiveMutator(GAGenome &, long double);
  static int SwapNodeMutator(GAGenome &, long double);
  static int SwapSubtreeMutator(GAGenome &, long double);
  static int OnePointCrossover(const GAGenome&, const GAGenome&, 
			      GAGenome*, GAGenome*);
  static long double TopologyComparator(const GAGenome&, const GAGenome&);
//  static long double NodeComparator(const GAGenome&, const GAGenome&);

public:
  GATreeGenome(GAGenome::Evaluator f=NULL, void * u=NULL);
  GATreeGenome(const GATreeGenome<T> &);
  GATreeGenome<T> & operator=(const GAGenome & orig)
    {copy(orig); return *this;}
  virtual ~GATreeGenome();
  virtual GAGenome *clone(GAGenome::CloneMethod flag=CONTENTS) const;
  virtual void copy(const GAGenome &);

#ifdef GALIB_USE_STREAMS
  virtual int write (STD_OSTREAM &) const;
#endif

  virtual int equal(const GAGenome & c) const;

// Here we do inlined versions of the access members of the super class.  We
// do our own here so that we can set/unset the _evaluated flag appropriately.

  int destroy() { _evaluated = gaFalse; return GATree<T>::destroy(); }
  int swaptree(GATree<T> * t)
    { _evaluated = gaFalse; return GATree<T>::swaptree(t); }
  int swaptree(unsigned int i, unsigned int j)
    { _evaluated = gaFalse; return GATree<T>::swaptree(i,j); }
  int swap(unsigned int i, unsigned int j)
    { _evaluated = gaFalse; return GATree<T>::swap(i,j); }
  GATree<T> * remove() { _evaluated = gaFalse; return GATree<T>::remove(); }
  int insert(GATree<T> * t, GATreeBASE::Location where=GATreeBASE::BELOW)
    { _evaluated = gaFalse; return GATree<T>::insert(t, where); }
  int insert(const T & t, GATreeBASE::Location where=GATreeBASE::BELOW)
    { _evaluated = gaFalse; return GATree<T>::insert(t, where); }
};

#ifdef GALIB_USE_BORLAND_INST
#include "GATreeGenome.cc"
#endif

#endif
