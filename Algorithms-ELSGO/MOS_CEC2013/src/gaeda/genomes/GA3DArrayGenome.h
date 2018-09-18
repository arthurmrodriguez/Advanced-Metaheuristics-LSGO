// $Header: /home/cvs/galib/ga/GA3DArrayGenome.h,v 1.4 2004/12/28 22:17:30 mwall Exp $
/* ----------------------------------------------------------------------------
  array3.h
  mbwall 25feb95
  Copyright (c) 1995 Massachusetts Institute of Technology
                     all rights reserved

 DESCRIPTION:
  This header defines the interface for the 3D array genome.  See comments in
1D array file.
---------------------------------------------------------------------------- */
#ifndef GA3DARRAYGENOME_H
#define GA3DARRAYGENOME_H

/* INCLUDES */
#include "GAArray.h"
#include "GAGenome.h"
#include "GAAllele.h"


/* ----------------------------------------------------------------------------
3DArrayGenome
---------------------------------------------------------------------------- */
template <class T>
class GA3DArrayGenome : public GAArray<T>, public GAGenome {
public:
  GADeclareIdentity();

  static int SwapMutator(GAGenome&, double);
  static double ElementComparator(const GAGenome&, const GAGenome&);
  static int UniformCrossover(const GAGenome&, const GAGenome&, 
			      GAGenome*, GAGenome*);
  static int EvenOddCrossover(const GAGenome&, const GAGenome&, 
			      GAGenome*, GAGenome*);
  static int OnePointCrossover(const GAGenome&, const GAGenome&, 
			      GAGenome*, GAGenome*);

public:
  GA3DArrayGenome(unsigned int x, unsigned int y, unsigned int z,
		  GAGenome::Evaluator f=(GAGenome::Evaluator)0,
		  void * u=(void *)0);
  GA3DArrayGenome(const GA3DArrayGenome & orig);
  GA3DArrayGenome& operator=(const GAGenome& orig)
    {copy(orig); return *this;}
  GA3DArrayGenome& operator=(const T array []){
    for(unsigned int i=0; i<nx; i++)
      for(unsigned int j=0; j<ny; j++)
	for(unsigned int k=0; k<nz; k++)
	  gene(i,j,k, *(array+k*ny*nx+j*nx+i));
    return *this;
  }
  virtual ~GA3DArrayGenome();
  virtual GAGenome * clone(GAGenome::CloneMethod flag=CONTENTS) const;
  virtual void copy(const GAGenome & chrom);

#ifdef GALIB_USE_STREAMS
  virtual int read(STD_ISTREAM &);
  virtual int write (STD_OSTREAM &) const;
#endif

  virtual int equal(const GAGenome & c) const;

  const T & gene(unsigned int x, unsigned int y, unsigned int z) const {
    return this->a[z*ny*nx + y*nx + x];
  }
  T & gene(unsigned int x, unsigned int y, unsigned int z, const T & value){
    if(this->a[z*ny*nx + y*nx + x] != value){ 
      _evaluated = gaFalse;
      this->a[z*ny*nx + y*nx + x] = value;
    }
    return this->a[z*ny*nx + y*nx + x];
  }
  int width() const {return nx;}
  int width(int w){resize(w, ny, nz); return nx;}
  int height() const {return ny;}
  int height(int h){resize(nx, h, nz); return ny;}
  int depth() const {return nz;}
  int depth(int d){resize(nx, ny, d); return nz;}
  virtual int resize(int x, int y, int z);
  int resizeBehaviour(Dimension which) const ;
  int resizeBehaviour(Dimension which,
		      unsigned int lower, unsigned int upper);
  int resizeBehaviour(unsigned int lowerX, unsigned int upperX, 
		      unsigned int lowerY, unsigned int upperY, 
		      unsigned int lowerZ, unsigned int upperZ){
    return(resizeBehaviour(WIDTH,  lowerX, upperX) * 
	   resizeBehaviour(HEIGHT, lowerY, upperY) * 
	   resizeBehaviour(DEPTH,  lowerZ, upperZ));
  }
  void copy(const GA3DArrayGenome &, 
	    unsigned int, unsigned int, unsigned int,
	    unsigned int, unsigned int, unsigned int,
	    unsigned int, unsigned int, unsigned int);
  void swap(unsigned int a, unsigned int b, unsigned int c,
	    unsigned int d, unsigned int e, unsigned int f){
    GAArray<T>::swap(c*ny*nx+b*nx+a, f*ny*nx+e*nx+d);
    _evaluated = gaFalse; 
  }


protected:
  unsigned int nx, ny, nz, minX, minY, minZ, maxX, maxY, maxZ;
};











/* ----------------------------------------------------------------------------
3DArrayAlleleGenome
---------------------------------------------------------------------------- */
template <class T>
class GA3DArrayAlleleGenome : public GA3DArrayGenome<T> {
public:
  GADeclareIdentity();

  static void UniformInitializer(GAGenome&);
  static int FlipMutator(GAGenome&, double);

public:
  GA3DArrayAlleleGenome(unsigned int x, unsigned int y, unsigned int z,
			const GAAlleleSet<T> & a,
			GAGenome::Evaluator f=(GAGenome::Evaluator)0,
			void * u=(void *)0);
  GA3DArrayAlleleGenome(unsigned int x, unsigned int y, unsigned int z,
			const GAAlleleSetArray<T> & a,
			GAGenome::Evaluator f=(GAGenome::Evaluator)0,
			void * u=(void *)0);
  GA3DArrayAlleleGenome(const GA3DArrayAlleleGenome<T> & orig);
  GA3DArrayAlleleGenome<T>& operator=(const GAGenome& orig)
    {copy(orig); return *this;}
  GA3DArrayAlleleGenome<T>& operator=(const T array [])
    {GA3DArrayGenome<T>::operator=(array); return *this;}
  virtual ~GA3DArrayAlleleGenome();
  virtual GAGenome * clone(GAGenome::CloneMethod flag=GAGenome::CONTENTS) const;
  virtual void copy(const GAGenome &);

#ifdef GALIB_USE_STREAMS
  virtual int read(STD_ISTREAM& is);
  virtual int write (STD_OSTREAM& os) const;
#endif

  virtual int equal(const GAGenome & c) const ;
  virtual int resize(int x, int y, int z);

  const GAAlleleSet<T>& alleleset(unsigned int i=0) const 
    {return aset[i%naset];}
public:

  virtual int fixedSize() const
    { 
      if (this->resizeBehaviour(GAGenome::DEPTH)      == GAGenome::FIXED_SIZE
	  &&  this->resizeBehaviour(GAGenome::WIDTH) == GAGenome::FIXED_SIZE
	  &&  this->resizeBehaviour(GAGenome::HEIGHT) == GAGenome::FIXED_SIZE) return this->nx * this->ny * this->nz;
      GAErr(GA_LOC, this->className(), "fixedSize", gaErrSameLengthReqd); 
      return -1;
    }

virtual int domain(int gno) const
{ 
  if(aset[gno%naset].type()==GAAllele::ENUMERATED ||
     aset[gno%naset].type()==GAAllele::DISCRETIZED)
    return aset[gno%naset].size();
  return -1;
}

  virtual long getValueOfNominalGene(int gno) const
    { 
      if(aset[gno%naset].type()==GAAllele::ENUMERATED ||
	 aset[gno%naset].type()==GAAllele::DISCRETIZED)
	{
	  for(int i=0;i<aset[gno%naset].size();i++)
	    if(aset[gno%naset].allele(i)==this->gene(gno%this->nx, (gno/this->nx)%this->ny, gno/(this->nx*this->ny)))
	      return (long)i;
	}
      GAErr(GA_LOC, this->className(), "getValueOfNominalGene", gaErrOpUndef);  
      return (long)-1L;
    }

  virtual long setValueOfNominalGene(int gno, long value)
    { 
      if(aset[gno%naset].type()==GAAllele::ENUMERATED ||
	 aset[gno%naset].type()==GAAllele::DISCRETIZED)
	return (long)this->gene(gno%this->nx, (gno/this->nx)%this->ny, gno/(this->nx*this->ny),aset[gno%naset].allele(value));
      GAErr(GA_LOC, this->className(), "setValueOfNominalGene", gaErrOpUndef);  
      return (long) -1L;
    }


  virtual double getValueofContinuousGene(int gno) const
    { 
      if(aset[gno%naset].type()==GAAllele::BOUNDED)
	return (double)this->gene(gno%this->nx, (gno/this->nx)%this->ny, gno/(this->nx*this->ny));
      GAErr(GA_LOC, this->className(), "getValueOfContinuousGene", gaErrOpUndef);  
      return (double)-1.0;
    }

  virtual double setValueOfContinuousGene(int gno, double value)
    { 
      if(aset[gno%naset].type()==GAAllele::BOUNDED)
	return (double)this->gene(gno%this->nx, (gno/this->nx)%this->ny, gno/(this->nx*this->ny),(T)value);
      GAErr(GA_LOC, this->className(), "setValueOfContinuousGene", gaErrOpUndef);  
      return (double) -1.0;
    }


  virtual double min(int gno) const
    { 
      if(aset[gno%naset].type()==GAAllele::BOUNDED)
	return (double)aset[gno%naset].lower();
      GAErr(GA_LOC, this->className(), "min", gaErrOpUndef);  
      return -1.0;
    }

  virtual double max(int gno) const
    { 
      if(aset[gno%naset].type()==GAAllele::BOUNDED)
	return (double)aset[gno%naset].upper();
      GAErr(GA_LOC, this->className(), "max", gaErrOpUndef);  
      return -1.0;
    }


protected:
  int naset;
  GAAlleleSet<T> * aset;
};


#ifdef GALIB_USE_BORLAND_INST
#include "GA3DArrayGenome.cc"
#endif

#endif
