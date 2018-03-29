// $Header: /home/cvs/galib/ga/GADCrowdingGA.h,v 1.1.1.1 1999/11/11 18:56:03 mbwall Exp $
/* ----------------------------------------------------------------------------
  dcrowdingga.h
  mbwall 29mar99
  Copyright (c) 1999 Matthew Wall, all rights reserved

  Header file for the steady-state genetic algorithm class.
---------------------------------------------------------------------------- */

#ifndef GADCROWDINGGA_H
#define GADCROWDINGGA_H

/* INCLUDES */
#include "GAGeneticAlgorithm.h"

class GADCrowdingGA : public GAGeneticAlgorithm {

   public:

      GADefineIdentity ("GADeterministicCrowdingGA", 241);

      GADCrowdingGA (const GAGenome& g) : GAGeneticAlgorithm (g) {}
      virtual ~GADCrowdingGA () {}

      virtual void initialize ();
      virtual void step ();
      GADCrowdingGA& operator++ () {step (); return *this;}

};

#endif
