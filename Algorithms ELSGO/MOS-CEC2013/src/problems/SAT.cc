/* ----------------------------------------------------------------------------
   Resolución del problema SAT 4blocks con 758 atributos
   ---------------------------------------------------------------------------- */

#include <GAIntOps.h>
#include <MOSGenomeFactory.h>
#include <genomes/GA1DArrayGenome.h>

int GENE [758];

extern "C" double calculate_sat (int* genes);

extern "C" double objective (GAGenome& g) {

   GA1DArrayAlleleGenome<int>* genome = dynamic_cast<GA1DArrayAlleleGenome<int>*>(&g);

   for (int i = 0; i < genome->length (); i++)
      GENE [i] = genome->gene (i);

   return calculate_sat (GENE);

}

extern "C" void individualInit (GAGenome& g) {
   return IntegerUniformInitializer (g);
}

extern "C" GAGenome* defineProblem () {

   // Create the genome and the allele set
   GAAlleleSet<int> alleles;
   alleles.add (0);
   alleles.add (1);
   GA1DArrayAlleleGenome<int>* genome = new GA1DArrayAlleleGenome<int> (758, alleles, objective);

   // Common operators
   genome->initializer (IntegerUniformInitializer);
   genome->comparator  (IntegerElementComparator);

   // Specific stuff for GAs
   genome->crossover   (IntegerTwoPointsCrossover);
   genome->mutator     (IntegerSwapMutator);

   // Specific stuff for DE
   genome->crossover   (IntegerExponentialCrossover);

   // Specific stuff for MOS
   MOSGenomeFactory::handle()->registerGenome (GAID::IntegerEncoding, genome);

   return genome;

}

extern "C" const char* describeProblem (void) {
   return "SAT Problem 4blocks";
}
