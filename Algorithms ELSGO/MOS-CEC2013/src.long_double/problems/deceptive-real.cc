#include <GARealOps.h>
#include <MOSGenomeFactory.h>
#include <genomes/GA1DArrayGenome.h>
#include <GAEDAConfig.h>
#include <GAGeneticAlgorithm.h>

const long double DEC_MIN = 0.0;
const long double DEC_MAX = 1.0;

int global_size;
GA1DArrayAlleleGenome<long double>* optimum_genome = NULL;

long double fdeceptive (long double x, long double y) {

   long double result = 0.0;
   long double temp = sqrt ((x*x + y*y) / 2.0);

   if (temp <= 0.8 )
      result = 0.8 - temp;
   else
      result = -4.0 + (5.0 * temp);

   return result;

}

extern "C" long double objective (GAGenome& g) {

   GA1DArrayAlleleGenome<long double>& genome = DYN_CAST (GA1DArrayAlleleGenome<long double>&, g);

   long double result = 0.0;

   for (int i = 0; i < genome.length () / 2; i++)
      result += fdeceptive (genome.gene (2 * i), genome.gene (2 * i + 1));

   return result;

}

extern "C" void individualInit (GAGenome& g) {
   return RealUniformInitializer (g);
}

extern "C" GAGenome* defineProblem () {

   GAEDAConfig* cfg = GAEDAConfig::handle();

   global_size = cfg->getProblemSize();

   GAAlleleSet<long double> alleles (DEC_MIN, DEC_MAX);
   GA1DArrayAlleleGenome<long double>* genome = new GA1DArrayAlleleGenome<long double> (cfg->getProblemSize(), alleles, objective);

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

   optimum_genome = new GA1DArrayAlleleGenome<long double> (cfg->getProblemSize(), alleles, objective);
   genome->initializer (RealUniformInitializer);
   genome->comparator  (RealEuclideanComparator);
   genome->crossover   (RealBlendCrossover);
   genome->mutator     (RealGaussianMutator);

   for (int i = 0; i < optimum_genome->length (); i++)
      optimum_genome->gene (i, 1.0);

   return genome;

}

extern "C" bool postprocess (GAPopulation* pop, int rank, int repetition) {
  delete optimum_genome;
  return true;
}

extern "C" GAGenome* getOptimumGenome () {
   return optimum_genome;
}

extern "C" int optCriterion() {
   return GAGeneticAlgorithm::MAXIMIZE;
}

extern "C" const char *describeProblem (void) {
   return "Deceptive Real";
}
