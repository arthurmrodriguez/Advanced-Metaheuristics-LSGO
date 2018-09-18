#include <GARealOps.h>
#include <MOSGenomeFactory.h>
#include <genomes/GA1DArrayGenome.h>
#include <GAEDAConfig.h>

const float RAST_MIN = -5.12;
const float RAST_MAX =  5.12;

extern "C" double objective (GAGenome& g) {

   GA1DArrayAlleleGenome<double>& genome = DYN_CAST (GA1DArrayAlleleGenome<double>&, g);

   int size = genome.length ();
   double result = 10.0 * size;

   for (int i=0; i<size; i++)
      result += (pow (genome.gene (i), 2) - 10 * cos (2 * M_PI * genome.gene (i)));

   return 1 / (result + 1);

}

extern "C" void individualInit (GAGenome& g) {
   return RealUniformInitializer (g);
}

extern "C" GAGenome* defineProblem () {

   GAEDAConfig* cfg = GAEDAConfig::handle();

   // Definition of the genome and the allele set
   GAAlleleSet<double> alleles (RAST_MIN, RAST_MAX);
   GA1DArrayAlleleGenome<double>* genome = new GA1DArrayAlleleGenome<double> (cfg->getProblemSize(), alleles, objective);

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

extern "C" const char* describeProblem (void) {
   return "Rastrigin problem";
}
