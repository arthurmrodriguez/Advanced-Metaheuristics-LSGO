#ifndef CEC2013_LSGO_H_
#define CEC2013_LSGO_H_

#include <stdexcept>
#include <limits>

#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/toplev.h>

#include <GARealOps.h>
#include <MOSGenomeFactory.h>
#include <GAEDAConfig.h>
#include <GAGeneticAlgorithm.h>
#include <genomes/GA1DArrayGenome.h>

GAAlleleSet<double>* alleles = 0;

std::string pname;

unsigned DIMS = 0;
double *sol;
double obj_val[10];

void initializeAlleles (unsigned func);

extern "C" double objective (GAGenome& g) {
   GA1DArrayAlleleGenome<double>& gen = dynamic_cast<GA1DArrayAlleleGenome<double>&> (g);

   octave_value_list evalArguments;

   Matrix sol (1, gen.size());

   for (unsigned i = 0; i < gen.size(); i++)
      sol(0, i) = gen.gene(i);

   evalArguments(0) = sol;
   evalArguments(1) = XXX;

   const octave_value_list fitness = feval("benchmark_func", evalArguments, 1);

   if (isnan(fitness(0).scalar_value()))
     return numeric_limits<double>::max();
   else
     return fitness(0).scalar_value();
}

extern "C" void individualInit (GAGenome& g) {
   return RealUniformInitializer (g);
}

extern "C" GAGenome* defineProblem (int size, char *data) {
   GAEDAConfig* cfg = GAEDAConfig::handle();

   // Initialize Octave Wrapper
   const char* argvv[] = {"", "--silent"};
   octave_main(2, (char**) argvv, true);

   octave_value_list pathArguments;
   pathArguments(0) = cfg->getProblemData();
   const octave_value_list resultPath = feval("addpath", pathArguments, 1);

   octave_value_list initializeArguments;
   initializeArguments(0) = cfg->getProblemData();
   const octave_value_list result = feval("initializeBenchmark", initializeArguments, 1);

   // The problem size is 1000 for every function except for f13
   DIMS = 1000;
   if (XXX == 13) DIMS = 905;

   sol = new double[DIMS];

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
   delete sol;
   do_octave_atexit ();
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
