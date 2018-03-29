/*
 * The following initial values give the settings described by
 * John Holland at ICGA93 and in his posting to gadistr
 * (GA-List v7n22). See the DESCRIPTION.ps file for more details.
 */

#include <math.h>

#include <GARealOps.h>
#include <MOSGenomeFactory.h>
#include <genomes/GA1DArrayGenome.h>

const double MIN1 = -5.0;
const double MAX1 = 10.0;
const double MIN2 =  0.0;
const double MAX2 = 15.0;

const double EXTRA = 1000.0;

extern "C" double objective (GAGenome& g) {

   GA1DArrayAlleleGenome<double>& genome = DYN_CAST (GA1DArrayAlleleGenome<double>&, g);

   double score = EXTRA;

   double x1 = genome.gene (0),
          x2 = genome.gene (1);

   score -= pow (x2 - (5 * x1 * x1 / (4 * M_PI * M_PI)) + (5 * x1 / M_PI) - 6, 2) +
            10 * (1 - (1 / (8 * M_PI))) * cos (x1) +
            10;

   return score;

}

extern "C" void individualInit (GAGenome& g) {
   return RealUniformInitializer (g);
}

extern "C" GAGenome* defineProblem () {

   GAAlleleSetArray<double> alleles;
   alleles.add (MIN1, MAX1, GAAllele::EXCLUSIVE, GAAllele::EXCLUSIVE);
   alleles.add (MIN2, MAX2, GAAllele::EXCLUSIVE, GAAllele::EXCLUSIVE);
   GA1DArrayAlleleGenome<double>* genome = new GA1DArrayAlleleGenome<double> (alleles, objective);

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
   return "Branin problem";
}
