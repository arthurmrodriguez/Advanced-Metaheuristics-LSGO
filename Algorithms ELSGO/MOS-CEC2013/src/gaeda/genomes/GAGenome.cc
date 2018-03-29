// $Header: /home/cvs/galib/ga/GAGenome.C,v 1.1.1.1 1999/11/11 18:56:03 mbwall Exp $
/* ----------------------------------------------------------------------------
   genome.C
   mbwall 19apr95
   Copyright (c) 1995 Massachusetts Institute of Technology
   all rights reserved

DESCRIPTION:
Definitions for genome base class.  See the header file for complete
documentation for deriving new classes.  Comments here are implementation-
specific details about base class member functions.
---------------------------------------------------------------------------- */
#include <limits>
#include <math.h>
#include <vector>

#include "GAGenome.h"
#include "../GAEDAConfig.h"

GAGenome::OptCriterion GAGenome::opt_criterion__ = MAXIMIZATION;

// These are the default genome operators.
// None does anything - they just post an error message to let you know that no
// method has been defined.  These are for the base class (which has no
// function by itself).

void GAGenome::NoInitializer (GAGenome & c) {

   GAErr (GA_LOC, c.className (), "initializer", gaErrOpUndef);

}


int GAGenome::NoSexualCrossover (const GAGenome & p1, const GAGenome & p2, GAGenome * c1, GAGenome * c2) {

   GAErr (GA_LOC, p1.className (), "crossover", gaErrOpUndef);
   return 0;

}


int GAGenome::NoMutator (GAGenome & c, double) {

   GAErr (GA_LOC, c.className (), "mutator", gaErrOpUndef);
   return 0;

}


double GAGenome::NoComparator (const GAGenome& c, const GAGenome&) {

   GAErr (GA_LOC, c.className (), "comparator", gaErrOpUndef);
   return -1.0;

}


GAGenome::GAGenome (Initializer i, Mutator m, Comparator c) {

   if (i == 0)
      i = NoInitializer;

   if (m == 0)
      m = NoMutator;

   if (c == 0)
      c = NoComparator;

   _score     = _fitness = 0.0;
   _evaluated = gaFalse;
   _neval     = 0;

   ga   = 0;
   ud   = 0;
   eval = 0;
   evd  = 0;

   init = i;
   mutr = m;
   cmp  = c;

   sexcross  = 0;
   asexcross = 0;

   age_          = 0;
   gen_of_birth_ = 0;
   origin_       = 0;

   id = 0;
   island = -1;
   myage = 0;
}


GAGenome::GAGenome (const GAGenome & orig) {

   evd    = 0;
   _neval = 0;

   GAGenome::copy (orig);

}


GAGenome::~GAGenome () {

   delete evd;

}


GAGenome* GAGenome::clone (CloneMethod) const {

   GAErr (GA_LOC, className (), "clone", gaErrOpUndef);
   return new GAGenome (*this);

}


// The eval count is not copied from the other genome - that would 
// inflate the count.

void GAGenome::copy (const GAGenome & orig) {

   if (&orig == this)
      return;

   _score     = orig._score;
   _fitness   = orig._fitness;
   _evaluated = orig._evaluated;

   ga   = orig.ga;
   ud   = orig.ud;
   eval = orig.eval;

   init = orig.init;
   mutr = orig.mutr;
   cmp  = orig.cmp;

   sexcross  = orig.sexcross;
   asexcross = orig.asexcross;

   age_          = orig.age_;
   origin_       = orig.origin_;
   gen_of_birth_ = orig.gen_of_birth_;

   id = orig.id;
   island = orig.island;
   myage = orig.myage;

   _neval = 0;

   if (orig.evd) { // don't delete if c doesn't have one
      if (evd)
         evd->copy (*orig.evd);
      else
         evd = orig.evd->clone ();
   }

}


double GAGenome::evaluate (GABoolean flag) const {

   if (_evaluated == gaFalse || flag == gaTrue) {

      GAGenome* This = (GAGenome*) this;

      if (eval) {

         This->_neval++;
         This->_score = (*eval) (*This);

      }

      This->_evaluated = gaTrue;

   }

   return _score;

}


int GAGenome::age () const {

   return age_;

}


void GAGenome::incAge () {

   age_++;

}


void GAGenome::setNewIndAge () {

  age_ = -1; // -1 because the next thing that we do with the population is to increment each individual age. With -1 we
             // assure that the generated children have 0 as it age.

}


int GAGenome::genBirth () const {

  return gen_of_birth_;

}


void GAGenome::genBirth (int gen) {

   gen_of_birth_ = gen;

}


int GAGenome::origin () const {

   return origin_;

}


void GAGenome::origin (int origin) {

  origin_ = origin;

}


void GAGenome::readObject (STD_ISTREAM& is) {

   is.read ((char*) (&_score),        sizeof (_score)       );
   is.read ((char*) (&age_),          sizeof (age_)         );
   is.read ((char*) (&gen_of_birth_), sizeof (gen_of_birth_));
   is.read ((char*) (&origin_),       sizeof (origin_)      );

   is.read ((char*) (&id),            sizeof (id)           );
   is.read ((char*) (&island),        sizeof (island)       );
   is.read ((char*) (&myage),         sizeof (myage)        );

   _evaluated = gaTrue; // Since we set the score we dont need to reevaluate it.

}


/**
 * Deserialization of the binary value enclosed in the stream. Note that the deserialization must be done in 
 * inverse order than the readObject.
 */

void GAGenome::writeObject(ostream& os) const{

   os.write ((char*) (&_score),        sizeof (_score)       );
   os.write ((char*) (&age_),          sizeof (age_)         );
   os.write ((char*) (&gen_of_birth_), sizeof (gen_of_birth_));
   os.write ((char*) (&origin_),       sizeof (origin_)      );

   os.write ((char*) (&id),            sizeof (id)           );
   os.write ((char*) (&island),        sizeof (island)       );
   os.write ((char*) (&myage),         sizeof (myage)        );

}


void GAGenome::setWorstScore() {
  score(worstPossibleScore());
}


bool GAGenome::hasBetterScoreThan ( const GAGenome& g) const {
  return (compareScores(score(),g.score())) == BETTER;
}


double GAGenome::computeFitnessIncrement (double old_score) const {
  // To be sure that the score value has been updated
  evaluate();
//  return GAGenome::compareScores(_score, old_score) == GAGenome::BETTER ? (fabs((_score - old_score) / old_score) + 1.0) / 2.0 : 0.0;
  return GAGenome::compareScores(_score, old_score) == GAGenome::BETTER ? fabs((_score - old_score) / old_score) : 0.0;
//  return GAGenome::compareScores(_score, old_score) == GAGenome::BETTER ? 1.0 : 0.0;
}


double GAGenome::bestPossibleScore() {
  return (opt_criterion__ == MINIMIZATION) ? std::numeric_limits<double>::min() : std::numeric_limits<double>::max();
}

double GAGenome::worstPossibleScore() {
  return (opt_criterion__ == MINIMIZATION) ? std::numeric_limits<double>::max() : std::numeric_limits<double>::min();
}


GAGenome::ScoreComparison GAGenome::compareScores(double score1, double score2) {
  if (score1 == score2) return EQUAL;
  else {

    switch(opt_criterion__) {
    case MINIMIZATION:
      return (score1 < score2) ? BETTER : WORSE;
      break;
    case MAXIMIZATION:
      return (score1 > score2) ? BETTER : WORSE;
      break;
    }
  }
}


GAGenome::OptCriterion GAGenome::optCriterion () {
  return opt_criterion__;
}


void GAGenome::optCriterion (OptCriterion crit) {
  opt_criterion__ = crit;
}


bool GAGenome::precissionReached () {
  if (GAEDAConfig::handle()->usePrecission())
    return (_score < GAEDAConfig::handle()->getOptimumFitness() + GAEDAConfig::handle()->getPrecission()) ? true : false;
  else
    return false;
}
