#include "VNSDARPGenome.h"
#include "DARPEvaluator.h"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <math.h>


/*
 * Parragh initial values
 * evalWeightsUpdateValueMin_ = 0.05
 * evalWeightsUpdateValueMax_ = 0.1
 * evalWeightsUpdateValue_ (delta)  = random value between the previous values
 */

VNSDARPGenome::VNSDARPGenome(unsigned x, long double minevalweightsvalue, long double maxevalsweightsvalue, GAGenome::Evaluator f, void* u)
  : DARPGenome(x,0,f,u),
    minEvalWeightsUpdateValue_(minevalweightsvalue),
    maxEvalWeightsUpdateValue_(maxevalsweightsvalue) {
  evalWeightsUpdateValue_ = GARandomDouble(minevalweightsvalue,maxevalsweightsvalue);
}

void VNSDARPGenome::copy (const GAGenome& orig) {
  DARPGenome::copy(orig);

  const VNSDARPGenome* other = dynamic_cast<const VNSDARPGenome*>(&orig); assert(other);

  minEvalWeightsUpdateValue_ = other->minEvalWeightsUpdateValue_;
  maxEvalWeightsUpdateValue_ = other->maxEvalWeightsUpdateValue_;
}

int VNSDARPGenome::write (STD_OSTREAM& os) const {
  os << endl;
  os << "minEvalWeightsUpdateValue_: " << minEvalWeightsUpdateValue_ << endl;
  os << "maxEvalWeightsUpdateValue_: " << maxEvalWeightsUpdateValue_ << endl;

  DARPGenome::write(os);

  return 0;
}

void VNSDARPGenome::writeObject(ostream& os) const{
  DARPGenome::writeObject(os);

  os.write ( (char*) (&minEvalWeightsUpdateValue_), sizeof(minEvalWeightsUpdateValue_) );
  os.write ( (char*) (&maxEvalWeightsUpdateValue_), sizeof(maxEvalWeightsUpdateValue_) );
}

void VNSDARPGenome::readObject (istream& is) {
  DARPGenome::readObject(is);

  is.read ( (char*) (&minEvalWeightsUpdateValue_), sizeof (minEvalWeightsUpdateValue_) );
  is.read ( (char*) (&maxEvalWeightsUpdateValue_), sizeof (maxEvalWeightsUpdateValue_) );
}

void VNSDARPGenome::updatePenalizations() {
  evalWeightsUpdateValue_ = GARandomDouble(minEvalWeightsUpdateValue_,maxEvalWeightsUpdateValue_);
  DARPGenome::updatePenalizations();
}
