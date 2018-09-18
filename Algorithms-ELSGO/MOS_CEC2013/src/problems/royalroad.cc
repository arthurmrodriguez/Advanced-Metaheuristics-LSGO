#include <GAIntOps.h>
#include <MOSGenomeFactory.h>
#include <genomes/GA1DArrayGenome.h>
#include <GAEDAConfig.h>

double optimum = 1.0;

GA1DArrayAlleleGenome<int>* optimum_genome = NULL;

// Royal road problem
static int    b     = 8;
static int    g     = 7;
static int    k     = 4;
static int    mstar = 4;
static double u     = 0.3;
static double ustar = 1.0;
static double v     = 0.02;

/*
 * Other variables used to compute fitnesses and keep track of levels.
 */
static int*    block_flags;
static int     nblocks;
static double* part_scores;
static int     region_length;

extern "C" double objective (GAGenome& g) {

   GA1DArrayAlleleGenome<int>& genome = DYN_CAST (GA1DArrayAlleleGenome<int>&, g);

   // Del original
   double fitness = 0.0;
   register int block;
   register int level;
   register int i;

   // Clear the array that indicates which blocks are complete.
   for (i = 0; i < nblocks; i++)
      block_flags [i] = 0;

   // PART calculation - rewards (or otherwise) for partial blocks.

   for (block = 0; block < nblocks; block++) {

      // Count number of bits in the block "base".
      register int base = block * region_length;
      register int ones_count = 0;
      register int limit = base + b;

      for (i = base; i < limit; i++)
         if (genome.gene (i) == 1)
            ones_count++;

      if (ones_count < b)
         fitness += part_scores [ones_count];
      else {

         /*
          * There is no PART reward for a complete block. But we
          * record the block complete for the BONUS computation.
          */

         block_flags [block] = 1;

      }

   }

   // BONUS calculations - rewards for complete blocks at all levels.

   for (level = 0; level <= k; level++) {

      register int nblocks_this_level = 1 << (k - level);
      register int complete_block_sets = 0;

      // Count the number of completed blocks at this level.
      for (i = 0; i < nblocks_this_level; i++)
         complete_block_sets += block_flags [i];

      // Give u* for the first and u for additional ones.
      if (complete_block_sets > 0)
         fitness += ustar + u * (double) (complete_block_sets - 1);

#ifdef NO_INTERMEDIATE_LEVELS
      break;
#endif

      /*
       * We can break if we don't have at least 2 completed block sets,
       * since there is no way we could have a target at the next level.
       */

      if (complete_block_sets < 2)
         break;

      // Setup next level
      for (i = 0; i < nblocks_this_level; i += 2)
         block_flags [i >> 1] = block_flags [i] & block_flags [i + 1];

   }

   return (10+fitness)/optimum;

   // Del original

}

void initEvaluator (void) {

   // The number of level zero blocks.
   nblocks = 1 << k; // O SEA = 16

   // An array to keep track of completed blocks.
   block_flags = (int*) malloc (nblocks * sizeof (int));

   // An array to save re-computing PART scores.
   part_scores = (double*) malloc ((b + 1) * sizeof (double));

   if (!block_flags || !part_scores) {

      fprintf (stderr, "Could not malloc() in jhrr_initialize().\n");
      exit (1);

   }

   /*
    * Initialize an array to contain the PART scores based on the
    * number of ones in the block.
    *
    * We don't need to initialize the second for loop, but it's clearer.
    */

   for (int i = 0; i <= mstar; i++)
      part_scores [i] = (double) i * v;

   for (int i = mstar + 1; i < b; i++)
      part_scores [i] = (double) (i - mstar) * -v;

    part_scores [b] = 0.0;

    // This saves a tiny bit of computation when evaluating.
    region_length = b + g;

}

extern "C" void individualInit (GAGenome& g) {
   return IntegerUniformInitializer (g);
}

extern "C" GAGenome* defineProblem () {

   GAEDAConfig* cfg = GAEDAConfig::handle();

   // Create the genome and the allele set
   GAAlleleSet<int> alleles;
   alleles.add (0);
   alleles.add (1);

   GA1DArrayAlleleGenome<int>* genome = new GA1DArrayAlleleGenome<int> (cfg->getProblemSize(), alleles, objective);

   // Create and evaluate the optimum genome
   initEvaluator ();

   optimum_genome = new GA1DArrayAlleleGenome<int> (cfg->getProblemSize(), alleles, objective);
   optimum_genome->resize(GAGenome::ANY_SIZE);

   for (int i = 0; i < optimum_genome->length (); i++)
      optimum_genome->gene (i, 1);

   optimum = objective (*optimum_genome);

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

extern "C" bool postprocess (GAPopulation* pop, int rank, int repetition) {
  free (block_flags);
  free (part_scores);
  delete optimum_genome;
  return true;
}

extern "C" const char* describeProblem (void) {
   return "Holland's Royal Road problem";
}
