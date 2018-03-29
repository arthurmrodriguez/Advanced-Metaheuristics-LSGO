#include <gaid.h>
#include <genomes/GA1DArrayGenome.h>
#include <genomes/MOSGenome.h>
#include <islands/CommManager.h>
#include <GAEDAConfig.h>
#include <GARealOps.h>
#include <MOSTechniqueLS.h>
#include <MOSGenomeFactory.h>
#include <garandom.h>

#include "../../src/problems/cec2008/shifted_rastrigin.h"

extern "C" double objective (GAGenome& g) {

  GA1DArrayAlleleGenome<double>& genome = DYN_CAST (GA1DArrayAlleleGenome<double>&, g);
  unsigned dim = genome.length ();

  double z = 0.0;
  double F = 0.0;
  double pi= acos (-1.0);
  double result = 0.0;

  for (unsigned i = 0; i < dim; i++) {
    z = genome.gene (i) - coords_opt [i];
    F = F + (pow (z, 2) - 10 * cos (2 * pi * z) + 10);
  }

  result = F /*+ bias*/;

  return 1 / (result + 1);
}

extern "C" GAGenome* defineProblem (int size, char *data) {

   GAAlleleSet<double> alleles (-5, 5);
   GA1DArrayAlleleGenome<double>* genome = new GA1DArrayAlleleGenome<double> (100, alleles, objective);

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


int main (int argc, char** argv) {

   CommManager cm (argc, argv);
   GAEDAConfig* cfg = GAEDAConfig::handle (cm);

   GARandomSeed (1234);

   defineProblem (0, (char*) "");

   MOSGenomeFactory* genFact = MOSGenomeFactory::handle();

   MOSTechniqueLS tech (0, "", MTS_LS1, RealEuclideanComparator, RealUniformInitializer, objective, GAID::RealEncoding, genFact->getGenome (GAID::RealEncoding), new GATournamentSelector(2));

   MOSGenome gen (&tech);

   GAPopulation* origPop = new GAPopulation (gen, 30);
   GAPopulation* destPop = new GAPopulation (gen, 30);

   unsigned totalEvals = 0;

   do {
      tech.offspring (origPop, destPop, 1, 0);
      std::swap (origPop, destPop);
      totalEvals++;
   } while (totalEvals < 500000);

   return 0;

}
