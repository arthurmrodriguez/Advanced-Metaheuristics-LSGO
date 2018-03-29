// $Header: /home/cvs/galib/ga/GA1DArrayGenome.C,v 1.3 2004/12/28 22:17:30 mwall Exp $
/* ----------------------------------------------------------------------------
  array1.C
  mbwall 25feb95
  Copyright (c) 1995 Massachusetts Institute of Technology
                     all rights reserved

 DESCRIPTION:
  Source file for the 1D array genome.
---------------------------------------------------------------------------- */

#ifndef GA1DARRAYGENOME_CC
#define GA1DARRAYGENOME_CC

/* INCLUDES */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../garandom.h"
#include "GA1DArrayGenome.h"

/* ----------------------------------------------------------------------------
1DArrayGenome
---------------------------------------------------------------------------- */

template <class T> const char* GA1DArrayGenome<T>::className () const {return "GA1DArrayGenome";}

template <class T> int GA1DArrayGenome<T>::classID () const {return GAID::ArrayGenome;}

// Set all the initial values to NULL or zero, then allocate the space we'll
// need (using the resize method).  We do NOT call the initialize method at
// this point - initialization must be done explicitly by the user of the
// genome (eg when the population is created or reset).  If we called the
// initializer routine here then we could end up with multiple initializations
// and/or calls to dummy initializers (for example when the genome is
// created with a dummy initializer and the initializer is assigned later on).
// Besides, we default to the no-initialization initializer by calling the
// default genome constructor.
template <class T> GA1DArrayGenome<T>::GA1DArrayGenome (unsigned int length, GAGenome::Evaluator f, void* u) : GAArray<T> (length),
                                                                                                               GAGenome () {
   evaluator (f);
   userData (u);
   nx = minX = maxX = length;

   _stdDevs.resize(length);

   for (unsigned i = 0; i < length; i++)
     _stdDevs[i] = GAUnitGaussian();

}


// This is the copy initializer.  We set everything to the default values, then
// copy the original.  The Array creator takes care of zeroing the data.
template <class T> GA1DArrayGenome<T>::GA1DArrayGenome (const GA1DArrayGenome<T>& orig) : GAArray<T> (orig.sz), GAGenome () {

   GA1DArrayGenome<T>::copy (orig);

}


// Delete whatever we own.
template <class T> GA1DArrayGenome<T>::~GA1DArrayGenome () {}


// This is the class-specific copy method.  It will get called by the super
// class since the superclass operator= is set up to call ccopy (and that is
// what we define here - a virtual function).  We should check to be sure that
// both genomes are the same class and same dimension.  This function tries
// to be smart about they way it copies.  If we already have data, then we do
// a memcpy of the one we're supposed to copy.  If we don't or we're not the
// same size as the one we're supposed to copy, then we adjust ourselves.
//   The Array takes care of the resize in its copy method.
template <class T> void GA1DArrayGenome<T>::copy (const GAGenome& orig) {

   if (&orig == this)
      return;

   const GA1DArrayGenome<T>* c = DYN_CAST (const GA1DArrayGenome<T>*, &orig);

   if (c) {

      GAGenome::copy (*c);
      GAArray<T>::copy (*c);

      nx   = c->nx;
      minX = c->minX;
      maxX = c->maxX;

      _stdDevs.resize(c->_stdDevs.size());
      assert(c->size() == c->_stdDevs.size());

      _stdDevs = c->_stdDevs;
   }

}


template <class T> GAGenome* GA1DArrayGenome<T>::clone (GAGenome::CloneMethod flag) const {

   GA1DArrayGenome<T>* cpy = new GA1DArrayGenome<T> (nx);

   if (flag == GAGenome::CONTENTS)
    cpy->copy (*this);
   else {

    cpy->GAGenome::copy (*this);

    cpy->maxX = maxX;
    cpy->minX = minX;

   }

   return cpy;

}


//   Resize the genome.
//   A negative value for the length means that we should randomly set the
// length of the genome (if the resize behaviour is resizeable).  If
// someone tries to randomly set the length and the resize behaviour is fixed
// length, then we don't do anything.
//   We pay attention to the values of minX and maxX - they determine what kind
// of resizing we are allowed to do.  If a resize is requested with a length
// less than the min length specified by the behaviour, we set the minimum
// to the length.  If the length is longer than the max length specified by
// the behaviour, we set the max value to the length.
//   We return the total size (in bits) of the genome after resize.
//   We don't do anything to the new contents!
template <class T> int GA1DArrayGenome<T>::resize (int len) {

   if (len == STA_CAST (int, nx))
      return nx;

   if (len == GAGenome::ANY_SIZE)
      len = GARandomInt (minX, maxX);
   else if (len < 0)
      return nx;  // do nothing
   else if (minX == maxX)
      minX = maxX = len;
   else {

      if (len < STA_CAST (int,minX))
         len=minX;

      if (len > STA_CAST(int,maxX))
         len=maxX;

  }

   nx = GAArray<T>::size (len);

   _stdDevs.resize(len);

   _evaluated = gaFalse;
   return this->sz;

}


#ifdef GALIB_USE_STREAMS
// We don't define this one apriori.  Do it in a specialization.
template <class T> int GA1DArrayGenome<T>::read (STD_ISTREAM&) {

   GAErr (GA_LOC, className (), "read", gaErrOpUndef);
   return 1;

}


// When we write the data to a stream we do it with spaces between elements.
// Also, there is no newline at the end of the stream of digits.
template <class T> int GA1DArrayGenome<T>::write (STD_OSTREAM& os) const {

   for (unsigned int i = 0; i < nx; i++)
      os << this->gene (i) << " ";

   return 0;

}
#endif


//   Set the resize behaviour of the genome.  A genome can be fixed
// length, resizeable with a max and min limit, or resizeable with no limits
// (other than an implicit one that we use internally).
//   A value of 0 means no resize, a value less than zero mean unlimited
// resize, and a positive value means resize with that value as the limit.
template <class T> int GA1DArrayGenome<T>::resizeBehaviour (unsigned int lower, unsigned int upper) {

   if (upper < lower) {

      GAErr (GA_LOC, className (), "resizeBehaviour", gaErrBadResizeBehaviour);
      return resizeBehaviour ();

   }

   minX = lower;
   maxX = upper;

   if (nx > upper)
      GA1DArrayGenome<T>::resize (upper);

   if (nx < lower)
      GA1DArrayGenome<T>::resize (lower);

   return resizeBehaviour ();

}


template <class T> int GA1DArrayGenome<T>::resizeBehaviour () const {

   int val = maxX;

   if (maxX == minX)
      val = GAGenome::FIXED_SIZE;

   return val;

}


template <class T> int GA1DArrayGenome<T>::equal (const GAGenome& c) const {

   const GA1DArrayGenome<T>& b = DYN_CAST (const GA1DArrayGenome<T>&, c);

   return ((this == &c) ? 1 : ((nx != b.nx) ? 0 : GAArray<T>::equal (b, 0, 0, nx)));

}


template <class T> int GA1DArrayGenome<T>::compsCompare (const GAGenome& g) const{
  const GA1DArrayGenome<T>& gen = dynamic_cast< const GA1DArrayGenome<T>& > (g);

  int ncomps=0;
  for (unsigned i=0; i<nx; i++) if (gen.gene(i) == gene(i) ) ncomps++;
  return ncomps;
}

/* ----------------------------------------------------------------------------
1DArrayAlleleGenome

  These genomes contain an allele set.  When we create a new genome, it owns
its own, independent allele set.  If we clone a new genome, the new one gets a
link to our allele set (so we don't end up with zillions of allele sets).  Same
is true for the copy constructor.
  The array may have a single allele set or an array of allele sets, depending
on which creator was called.  Either way, the allele set cannot be changed
once the array is created.
---------------------------------------------------------------------------- */

template <class T> const char* GA1DArrayAlleleGenome<T>::className () const {return "GA1DArrayAlleleGenome";}


template <class T> int GA1DArrayAlleleGenome<T>::classID () const {return GAID::ArrayAlleleGenome;}


template <class T> GA1DArrayAlleleGenome<T>::GA1DArrayAlleleGenome (unsigned int length, const GAAlleleSet<T>& s,
                                                                   GAGenome::Evaluator f, void* u) : GA1DArrayGenome<T> (length, f, u) {

   naset = 1;
   aset = new GAAlleleSet<T> [1];
   aset [0] = s;

   _SR = (alleleset(0).upper() - alleleset(0).lower()) / 2;
   _improve = 1;
   _lastIterPos = 0;
   _step = 0;

}


template <class T> GA1DArrayAlleleGenome<T>:: GA1DArrayAlleleGenome (const GAAlleleSetArray<T> & sa,
                                                                     GAGenome::Evaluator f, void * u) : GA1DArrayGenome<T> (sa.size (), f, u) {

   naset = sa.size ();
   aset = new GAAlleleSet<T> [naset];

   for(int i = 0; i < naset; i++)
      aset [i] = sa.set (i);

   _SR = (alleleset(0).upper() - alleleset(0).lower()) / 2;
   _improve = 1;
   _lastIterPos = 0;
   _step = 0;

}


// The copy constructor creates a new genome whose allele set refers to the
// original's allele set.
template <class T> GA1DArrayAlleleGenome<T>::GA1DArrayAlleleGenome (const GA1DArrayAlleleGenome<T>& orig) : GA1DArrayGenome<T> (orig.sz) {

   naset = 0;
   aset = (GAAlleleSet<T>*) 0;
   GA1DArrayAlleleGenome<T>::copy (orig);

}


// Delete the allele set
template <class T> GA1DArrayAlleleGenome<T>::~GA1DArrayAlleleGenome () {

   delete [] aset;

}


// This implementation of clone does not make use of the contents/attributes
// capability because this whole interface isn't quite right yet...  Just
// clone the entire thing, contents and all.
template <class T> GAGenome* GA1DArrayAlleleGenome<T>::clone (GAGenome::CloneMethod) const {

   return new GA1DArrayAlleleGenome<T> (*this);

}


template <class T> void GA1DArrayAlleleGenome<T>::copy (const GAGenome& orig) {

   if (&orig == this)
      return;

   const GA1DArrayAlleleGenome<T>* c = DYN_CAST (const GA1DArrayAlleleGenome<T>*, &orig);

   if(c) {

      GA1DArrayGenome<T>::copy (*c);

      if (naset != c->naset) {
         delete [] aset;
         naset = c->naset;
         aset = new GAAlleleSet<T> [naset];
      }

      for (int i = 0; i < naset; i++)
         aset [i].link (c->aset [i]);

      _SR = c->_SR;
      _improve = c->_improve;
      _lastIterPos = c->_lastIterPos;
      _step = c->_step;

   }

}


// If we resize to a larger length then we need to set the contents to a valid
// value (ie one of our alleles).
template <class T> int GA1DArrayAlleleGenome<T>::resize (int len) {

   unsigned int oldx = this->nx;
   GA1DArrayGenome<T>::resize (len);

   if (this->nx > oldx)
      for (unsigned int i = oldx; i < this->nx; i++)
         this->a [i] = aset [i % naset].allele ();

  return len;

}


// Define these so they can easily be specialized as needed.
#ifdef GALIB_USE_STREAMS
template <class T> int GA1DArrayAlleleGenome<T>::read (STD_ISTREAM& is) {

   return GA1DArrayGenome<T>::read (is);

}


template <class T> int GA1DArrayAlleleGenome<T>::write (STD_OSTREAM& os) const {

   return GA1DArrayGenome<T>::write (os);

}
#endif


template <class T> int GA1DArrayAlleleGenome<T>::equal (const GAGenome& c) const {

  return GA1DArrayGenome<T>::equal (c);

}


/**
 * Necessary in order to have a different implementation for the double case. See GARealGenome
 * for the implementation of this method
 */

template <class T> int GA1DArrayAlleleGenome<T>::compsCompare (const GAGenome& g) const {

   return GA1DArrayGenome<T>::compsCompare (g);

}


template <class T> void GA1DArrayGenome<T>::writeObject (ostream& os) const {
   GAGenome::writeObject(os);

   os.write ( (char*) (&nx),  sizeof (nx) );
   os.write ( (char*) (&minX),sizeof (minX) );
   os.write ( (char*) (&maxX),sizeof (maxX) );

   for (unsigned int i = 0; i<nx; i++) {
      os.write ((char*) (&(gene(i))), sizeof (gene(i)));
   }

}


template <class T> void GA1DArrayGenome<T>::readObject(istream& is) {
   GAGenome::readObject(is);

   int tmpnx;
   is.read ( (char*) (&tmpnx), sizeof (tmpnx) );
   is.read ( (char*) (&minX),  sizeof (minX) );
   is.read ( (char*) (&maxX),  sizeof (maxX) );

   resize(tmpnx);
   assert(nx==tmpnx);

   T value;
   for (unsigned i = 0; i <nx; i++) {
     is.read ((char*) (&value), sizeof (value));
     gene (i, value);
  }

}


template <class T> void GA1DArrayAlleleGenome<T>::writeObject (ostream& os) const {

   GA1DArrayGenome<T>::writeObject(os);

   os.write((char*)&_SR,sizeof(double));
   os.write((char*)&_improve,sizeof(char));
   os.write((char*)&_lastIterPos,sizeof(unsigned));
   os.write((char*)&_step,sizeof(unsigned));

}


template <class T> void GA1DArrayAlleleGenome<T>::readObject(istream& is) {

   GA1DArrayGenome<T>::readObject(is);

   is.read((char*)&_SR,sizeof(double));
   is.read((char*)&_improve,sizeof(char));
   is.read((char*)&_lastIterPos,sizeof(unsigned));
   is.read((char*)&_step,sizeof(unsigned));

}


template <class T> double GA1DArrayGenome<T>::de_cross_prob__ = 0.0;

#endif
