// $Header: /home/cvs/galib/ga/GA2DBinStrGenome.h,v 1.2 2004/12/28 00:12:11 mwall Exp $
/* ----------------------------------------------------------------------------
  binstr2.h
  mbwall 19apr95
  Copyright (c) 1995 Massachusetts Institute of Technology
                     all rights reserved

 DESCRIPTION:
  This header defines the interface for the 2D binary string genome, including
crossover objects and all the default and built-in operators.
---------------------------------------------------------------------------- */
#ifndef GA2DBINSTRGENOME_H
#define GA2DBINSTRGENOME_H

/* INCLUDES */
#include "GABinStr.h"
#include "GAGenome.h"

/* ----------------------------------------------------------------------------
2DBinaryStringGenome
-------------------------------------------------------------------------------
---------------------------------------------------------------------------- */
class GA2DBinaryStringGenome : public GABinaryString, public GAGenome {
	public:
		GADefineIdentity("GA2DBinaryStringGenome", GAID::BinaryStringGenome2D);

		static void UniformInitializer(GAGenome &);
		static void UnsetInitializer(GAGenome &);
		static void SetInitializer(GAGenome &);
		static int FlipMutator(GAGenome &, long double);
		static long double BitComparator(const GAGenome&, const GAGenome&);
		static int UniformCrossover(const GAGenome&, const GAGenome&,
				GAGenome*, GAGenome*);
		static int EvenOddCrossover(const GAGenome&, const GAGenome&,
				GAGenome*, GAGenome*);
		static int OnePointCrossover(const GAGenome&, const GAGenome&,
				GAGenome*, GAGenome*);

	public:
		GA2DBinaryStringGenome(unsigned int x, unsigned int y,
				GAGenome::Evaluator f=(GAGenome::Evaluator)0,
				void * u=(void *)0);
		GA2DBinaryStringGenome(const GA2DBinaryStringGenome & orig);
		GA2DBinaryStringGenome& operator=(const GAGenome& arg)
		{copy(arg); return *this;}
		GA2DBinaryStringGenome& operator=(const short array []){
			for(unsigned int i=0; i<nx; i++)
				for(unsigned int j=0; j<ny; j++)
					gene(i,j,*(array+j*nx+i));
			return *this;
		}
		GA2DBinaryStringGenome& operator=(const int array []){
			for(unsigned int i=0; i<nx; i++)
				for(unsigned int j=0; j<ny; j++)
					gene(i,j,*(array+j*nx+i));
			return *this;
		}
		virtual ~GA2DBinaryStringGenome();
		virtual GAGenome *clone(GAGenome::CloneMethod flag=CONTENTS) const;
		virtual void copy(const GAGenome & chrom);

#ifdef GALIB_USE_STREAMS
		virtual int read(STD_ISTREAM &);
		virtual int write (STD_OSTREAM &) const;
#endif

		virtual int equal(const GAGenome & c) const;

		// specific to this class
		short gene(unsigned int x, unsigned int y) const {return bit(x+nx*y);}
		short gene(unsigned int x, unsigned int y, short value) {
			_evaluated = gaFalse;
			return((bit(x+nx*y) == value) ? value : bit(x+nx*y, value));
		}
		int width() const {return nx;}
		int width(int w){resize(w, ny); return nx;}
		int height() const {return ny;}
		int height(int h){resize(nx, h); return ny;}
		int resize(int x, int y);
		int resizeBehaviour(Dimension which) const ;
		int resizeBehaviour(Dimension which,
				unsigned int lowerX, unsigned int upperX);
		int resizeBehaviour(unsigned int lowerX, unsigned int upperX,
				unsigned int lowerY, unsigned int upperY){
			return(resizeBehaviour(WIDTH, lowerX, upperX) *
					resizeBehaviour(HEIGHT, lowerY, upperY));
		}
		void copy(const GA2DBinaryStringGenome &,
				unsigned int, unsigned int,
				unsigned int, unsigned int,
				unsigned int, unsigned int);
		int equal(const GA2DBinaryStringGenome&,
				unsigned int, unsigned int,
				unsigned int, unsigned int,
				unsigned int, unsigned int) const;
		void set(unsigned int, unsigned int, unsigned int, unsigned int);
		void unset(unsigned int, unsigned int, unsigned int, unsigned int);
		void randomize(unsigned int, unsigned int, unsigned int, unsigned int);
		void randomize() { GABinaryString::randomize(); }
		void move(unsigned int, unsigned int,
				unsigned int, unsigned int,
				unsigned int, unsigned int);

	public:

		virtual int fixedSize() const
		{
			if (resizeBehaviour(WIDTH) == GAGenome::FIXED_SIZE
					&&  resizeBehaviour(HEIGHT) == GAGenome::FIXED_SIZE) return nx * ny;
			GAErr(GA_LOC, this->className(), "fixedSize", gaErrSameLengthReqd);
			return -1;
		}

		virtual int domain(int gno) const
		{
			return 2;
		}

		virtual long getValueOfNominalGene(int gno) const
		{
			return (long)gene(gno%nx, gno/nx);
		}

		virtual long setValueOfNominalGene(int gno, long value)
		{
			return (long)gene(gno%nx, gno/nx, (short)value);
		}

	protected:
		unsigned int nx, ny, minX, minY, maxX, maxY;
};

#endif
