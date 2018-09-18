#ifndef CEC2013_LSGO_H_
#define CEC2013_LSGO_H_

#include <stdexcept>
#include <limits>

#include <armadillo>
#include <matio.h>

#include <GARealOps.h>
#include <MOSGenomeFactory.h>
#include <GAEDAConfig.h>
#include <GAGeneticAlgorithm.h>
#include <genomes/GA1DArrayGenome.h>

#include "cec2013_funcs.h"

GAAlleleSet<double>* alleles = 0;
mat_t* matfp = NULL;
arma::mat* x = NULL;

std::string pname;

unsigned DIMS = 0;

void initializeAlleles (unsigned func);

extern "C" double objective (GAGenome& g) {
   GA1DArrayAlleleGenome<double>& gen = dynamic_cast<GA1DArrayAlleleGenome<double>&> (g);

   for (unsigned i=0; i < DIMS; i++)
      (*x)(i, 0)=gen.gene(i);

   double fit=fXXX(*x);

   if (isnan(fit))
     return numeric_limits<double>::max();
   else
     return fit;
}

extern "C" void individualInit (GAGenome& g) {
   return RealUniformInitializer (g);
}

extern "C" GAGenome* defineProblem (int size, char *data) {
   GAEDAConfig* cfg = GAEDAConfig::handle();

   // Open Matlab Data file and initialize data
   std::stringstream inputFile;
   std::string padding=(XXX < 10) ? "0" : "";
   inputFile << cfg->getProblemData() << "/f" << padding << XXX << ".mat";
   std::cout << inputFile.str() << std::endl;
   matfp=Mat_Open(inputFile.str().c_str(), MAT_ACC_RDONLY);
   init_fXXX(matfp);

   // The problem size is 1000 for every function except for f13 and f14
   DIMS = 1000;
   if (XXX == 13 || XXX == 14) DIMS = 905;

   x=new arma::mat(DIMS, 1);

   initializeAlleles(XXX);
   GA1DArrayAlleleGenome<double>* genome = new GA1DArrayAlleleGenome<double> (DIMS, *alleles, objective);

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
   delete x;
   Mat_Close(matfp);
   return true;
}

extern "C" const char *describeProblem () {
   return pname.c_str();
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
    case 11:
    case 12:
    case 13:
    case 14:
    case 15:
      alleles = new GAAlleleSet<double>(-100.0, 100.0);
      break;

    case 2:
    case 5:
    case 9:
      alleles = new GAAlleleSet<double>(-5.0, 5.0);
      break;

    case 3:
    case 6:
    case 10:
      alleles = new GAAlleleSet<double>(-32.0, 32.0);
      break;

    default:
      throw runtime_error("Error: alleles not defined for this function");
      break;
  }

  return;
}

#endif /* CEC2013_LSGO_H_ */
