#include <GARealOps.h>
#include <MOSGenomeFactory.h>
#include <genomes/GA1DArrayGenome.h>
#include <GAEDAConfig.h>

const double GRI_MIN = -600.0;
const double GRI_MAX =  600.0;

extern "C" double objective (GAGenome& g) {

   GA1DArrayAlleleGenome<double>& genome = DYN_CAST (GA1DArrayAlleleGenome<double>&, g);

   double score = 1;
   double prod  = 1.0;

   for (int i = 0; i < genome.length (); i++) {
      score += (genome.gene (i) * genome.gene (i) / 4000.0);
      prod  *= cos (genome.gene (i) / sqrt (i + 1.0));
   }

   score -= prod;

   return (1 / (score + 1));

}

extern "C" void individualInit (GAGenome& g) {
   return RealUniformInitializer (g);
}

extern "C" GAGenome* defineProblem () {

   GAEDAConfig* cfg = GAEDAConfig::handle();

   GAAlleleSet<double> alleles (GRI_MIN, GRI_MAX);
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
   return "Griewank problem";
}
