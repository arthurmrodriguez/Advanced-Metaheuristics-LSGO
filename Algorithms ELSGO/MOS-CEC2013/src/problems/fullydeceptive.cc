/* ----------------------------------------------------------------------------
Resolución del problema 4-bit fully deceptive k

Para K = 2,3,4
---------------------------------------------------------------------------- */

// Tipo: puede ser F2 (tambien lo he visto con el nombre de F4) o F3
#ifdef TYPE_F2
int valores [16] = {28, 26, 24, 18, 22, 16, 14, 0, 20, 12, 10, 2, 8, 4, 6, 30};
char NAME [] = "4-bit fully deceptive problem k/f2";
#else
int valores [16] = {10, 25, 26, 5, 27, 5, 5, 0, 28, 5, 5, 0, 5, 0, 0, 30};
char NAME [] = "4-bit fully deceptive problem k/f3";
#endif

#include <GAIntOps.h>
#include <MOSGenomeFactory.h>
#include <genomes/GA1DArrayGenome.h>

extern "C" double objective (GAGenome& g) {

   int bloque, pos;
   int fitness = 0;

   GA1DArrayAlleleGenome<int>* genome = dynamic_cast<GA1DArrayAlleleGenome<int>*>(&g);

   for (bloque = 0; bloque < 10; bloque++) {

      pos = 0;
      pos = genome->gene (4 * bloque    ) * 8 +
            genome->gene (4 * bloque + 1) * 4 +
            genome->gene (4 * bloque + 2) * 2 +
            genome->gene (4 * bloque + 3);

      fitness += valores [pos];

   }

   return fitness;

}

extern "C" void individualInit (GAGenome& g) {
   return IntegerUniformInitializer (g);
}

extern "C" GAGenome* defineProblem () {

   // Create the genome and the allele set
   GAAlleleSet<int> alleles;
   alleles.add (0);
   alleles.add (1);

   GA1DArrayAlleleGenome<int>* genome = new GA1DArrayAlleleGenome<int> (40, alleles, objective);

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
   return NAME;
}
