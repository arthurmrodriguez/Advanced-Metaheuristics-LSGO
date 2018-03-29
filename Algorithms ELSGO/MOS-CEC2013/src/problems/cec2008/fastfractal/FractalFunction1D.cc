#include "FractalFunction1D.h"
#include "DoubleDip.h"

#include <math.h>

FractalFunction1D::FractalFunction1D()
  : fractalDepth_(3),  // maximum recursive depth of 40 suggested for 64-bit architecture
    density_(1),
    index_(1),
    func_(true),
    ranTable_(NULL),
    unitFunction_(NULL)
{
  unitFunction_ = new DoubleDip();
}

FractalFunction1D::FractalFunction1D (UnitFunction1D* unitFunction, int fractalDepth,
				      int density, long index)
  : fractalDepth_(fractalDepth),  // maximum recursive depth of 40 suggested for 64-bit architecture
    density_(density),
    index_(index),
    func_(false),
    ranTable_(NULL),
    unitFunction_(unitFunction)
{
  ranTable_ = new RanTable(DOUBLE_TABLE_SIZE_, INT_TABLE_SIZE_, density_, index_);
}

FractalFunction1D::FractalFunction1D (UnitFunction1D* unitFunction, int density, long index)
  : fractalDepth_(3),  // maximum recursive depth of 40 suggested for 64-bit architecture
    density_(density),
    index_(index),
    func_(false),
    ranTable_(NULL),
    unitFunction_(unitFunction)
{
  ranTable_ = new RanTable(DOUBLE_TABLE_SIZE_, INT_TABLE_SIZE_, density_, index_);
}

FractalFunction1D::FractalFunction1D (UnitFunction1D* unitFunction, long index)
  : fractalDepth_(3),  // maximum recursive depth of 40 suggested for 64-bit architecture
    density_(1),
    index_(index),
    func_(false),
    ranTable_(NULL),
    unitFunction_(unitFunction)
{
  ranTable_ = new RanTable(DOUBLE_TABLE_SIZE_, INT_TABLE_SIZE_, density_, index_);
}

FractalFunction1D::FractalFunction1D (UnitFunction1D* unitFunction)
  : fractalDepth_(3),  // maximum recursive depth of 40 suggested for 64-bit architecture
    density_(1),
    index_(1),
    func_(false),
    ranTable_(NULL),
    unitFunction_(unitFunction)
{
  ranTable_ = new RanTable(DOUBLE_TABLE_SIZE_, INT_TABLE_SIZE_, density_, index_);
}

FractalFunction1D::~FractalFunction1D()
{
  if (ranTable_)
    delete ranTable_;
  if (func_)
    delete unitFunction_;
}

void FractalFunction1D::setIndex(long index) {
  index_ = index;
  ranTable_->setSeed(index);
}

#include <iostream>
long double FractalFunction1D::evaluate(long double x) {

  x = fmod (x, 1);  // first map pos into (0,1], (0,1]
  // note in Java -4.3%1 is -0.3 not 0.7 ie Matlab 'rem' function not 'mod'

  if (x <= 0)
    x = x+1;        // 0 must move to 1, or will be in wrong "square"

  if (fractalDepth_ < 1)
    return 0;       // check for valid depth argument, should never fire
  else
    return getDepthLocal(x, 1, index_, 1);  // start recursion

}

long double FractalFunction1D::getDepthLocal(long double x, int recDepth, long seed, long span) {

  long double depth = 0;
  long double scale = 1.0 / span;
  long square =  (long) ceil(x * span);
  long newSeed, square1;
  long double x1;

  // get contribution from each of the 3 relevant squares...
  for (int offset = -1; offset < 2; offset++) {

    x1 = x;
    square1 = square + offset;

    if (square1 == 0) {           // wrap to r.h.s.
      square1 = span;
      x1 = x1 + 1;
    }
    else if (square1 > span) {    // wrap to l.h.s.
      square1 = 1;
      x1 = x1 - 1;
    }

    // accumulate contributions
    depth = depth + getDepthWRTSquare(x1, square1, recDepth, seed, span, scale);


  }

  // now fire recursion to next level down...
  if (recDepth < fractalDepth_) {
    newSeed = (span + seed) & (DOUBLE_TABLE_SIZE_-1); // unique seeds up to random table size
    long newSpan = span << 1;                         // newSpan = 2^(recDepth-1), bit shift faster
    depth = depth + getDepthLocal(x,recDepth+1,newSeed,newSpan);  // recur to next level
  }

  return depth;

}

long double FractalFunction1D::getDepthWRTSquare (long double x, long square, int recDepth,
					     long seed, long span, long double scale) {

  long double depth = 0;
  // unique seed for this square
  long squareSeed = (square-1);
  // unique seed for square at depth up to table size
  long localSeed = (seed + squareSeed) & (DOUBLE_TABLE_SIZE_-1);

  // apply seed for this square
  ranTable_->setSeed(localSeed);
  // choose number of unit functions from uniform dist whose average is 'density'
  int numUnits = ranTable_->nextInteger();

  for (int i = 1; i <= numUnits; i++) { // get contribution from each

    // get diameter from quadratic distribution
    long double diameter = 1 / (2 - ranTable_->nextDouble()) * scale;
    // get centre from uniform distribution
    long double centre = (square - ranTable_->nextDouble()) * scale;

    // save making unnecessary call if unit function is too far to affect point
    if ((x-centre)*(x-centre) < diameter*diameter / 4) {
      unitFunction_->setCentre(centre);         // faster to set individually (believe it or not)
      unitFunction_->setScale(diameter);
      depth = depth + unitFunction_->getValue(x);  // add this unit function's contribution
    }

  }

  return depth;

}

int FractalFunction1D::DOUBLE_TABLE_SIZE_ = 0x3fff+1;   // 16384
int FractalFunction1D::INT_TABLE_SIZE_    = 0xff+1;     //   256
