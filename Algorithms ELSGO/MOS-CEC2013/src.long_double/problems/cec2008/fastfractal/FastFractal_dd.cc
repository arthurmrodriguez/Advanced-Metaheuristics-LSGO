#include <GARealOps.h>
#include <MOSGenomeFactory.h>
#include <genomes/GA1DArrayGenome.h>
#include <GAEDAConfig.h>
#include <GAGeneticAlgorithm.h>

#include "FastFractal.h"

const long double MIN_ALLELE_VALUE = -100;
const long double MAX_ALLELE_VALUE =  100;

GA1DArrayAlleleGenome<long double>* optimum_genome = NULL;

extern "C" long double objective (GAGenome& g) {

  GA1DArrayAlleleGenome<long double>& genome = DYN_CAST (GA1DArrayAlleleGenome<long double>&, g);

  unsigned dim = genome.length ();
  FastFractal* ff = new FastFractal ("DoubleDip", 3, 1, 1, dim);

  std::vector<long double> x (dim, 0.0);

  for (unsigned i = 0; i < dim; i++)
    x [i] = genome.gene (i);

  long double score = 100000 - ff->evaluate (x);

  delete ff;

  return score;

}

extern "C" void individualInit(GAGenome& g) {
  return RealUniformInitializer (g);
}

extern "C" GAGenome* defineProblem (int size, char *data) {

   GAEDAConfig* cfg = GAEDAConfig::handle();

   GAAlleleSet<long double> alleles (MIN_ALLELE_VALUE, MAX_ALLELE_VALUE);
   GA1DArrayAlleleGenome<long double>* genome = new GA1DArrayAlleleGenome<long double> (cfg->getProblemSize (), alleles, objective);

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

   optimum_genome = new GA1DArrayAlleleGenome<long double> (cfg->getProblemSize (), alleles, objective);
   optimum_genome->comparator (RealEuclideanComparator);
   optimum_genome->crossover  (RealBlendCrossover);

   for (int i = 0; i < optimum_genome->length (); i++)
     optimum_genome->gene (i, 0.0);

   return genome;

}

extern "C" bool postprocess (GAPopulation* pop, int rank, int repetition) {
  delete optimum_genome;
  return true;
}

extern "C" const char* describeProblem () {
  return "FastFractal DoubleDip";
}

extern "C" GAGenome* getOptimumGenome () {
  return optimum_genome;
}

extern "C" int optCriterion () {
  return GAGeneticAlgorithm::MINIMIZE;
}
