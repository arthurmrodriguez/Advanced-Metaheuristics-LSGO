#ifndef CEC2012_H_
#define CEC2012_H_

#include <stdexcept>

#include <GARealOps.h>
#include <MOSGenomeFactory.h>
#include <GAEDAConfig.h>
#include <genomes/GA1DArrayGenome.h>

#include "Header.h"

// Global variables
GAAlleleSet<long double>* alleles = 0;
Benchmarks *fp;

// Aux functions
void initializeAlleles (unsigned func);


extern "C" long double objective (GAGenome& g) {
   GA1DArrayAlleleGenome<long double>& gen = dynamic_cast<GA1DArrayAlleleGenome<long double>&> (g);
   return fp->compute((long double*) gen);
}

extern "C" void individualInit (GAGenome& g) {
   return RealUniformInitializer (g);
}

extern "C" GAGenome* defineProblem (int size, char *data) {

   GAEDAConfig* cfg = GAEDAConfig::handle();

   fp = new F19();

   initializeAlleles(19);
   GA1DArrayAlleleGenome<long double>* genome = new GA1DArrayAlleleGenome<long double> (cfg->getProblemSize(), *alleles, objective);

   // Common operators
   genome->initializer (RealUniformInitializer);
   genome->comparator  (RealEuclideanComparator);

   // Specific stuff for GAs
   genome->crossover   (RealBlendCrossover);
   genome->mutator     (RealGaussianMutator);

   // Specific stuff for DE
   genome->crossover   (RealExponentialCrossover);

   // Specific stuff for MOS
   MOSGenomeFactory::handle()->registerGenome (GAID::RealEncoding, genome);

   return genome;

}

extern "C" bool postprocess (GAPopulation* pop, int rank, int repetition) {
   delete alleles;
   delete fp;
   return true;
}

extern "C" const char *describeProblem () {
   return "CEC2012: F19";
}

extern "C" GAGenome::OptCriterion optCriterion(){
  return GAGenome::MINIMIZATION;
}

void initializeAlleles (unsigned func) {
  switch (func) {
    case 1:
    case 4:
    case 7:
    case 8:
    case 9:
    case 12:
    case 13:
    case 14:
    case 17:
    case 18:
    case 19:
    case 20:
      alleles = new GAAlleleSet<long double> (-100, 100);
      break;
    case 2:
    case 5:
    case 10:
    case 15:
      alleles = new GAAlleleSet<long double> (-5, 5);
      break;
    case 3:
    case 6:
    case 11:
    case 16:
      alleles = new GAAlleleSet<long double> (-32, 32);
      break;
    default:
      throw runtime_error("Error: alleles not defined for this function");
      break;
  }

  return;
}

#endif /* CEC2012_H_ */
