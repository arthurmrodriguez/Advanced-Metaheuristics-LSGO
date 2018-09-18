#include <GAIntOps.h>
#include <MOSGenomeFactory.h>
#include <genomes/GA1DArrayGenome.h>
#include <GAEDAConfig.h>

// Data storage definition
typedef struct {
   int* literals;
   int  nlits;
} clause_t;

typedef struct {
   clause_t* clauses;
   int nclaus;
   int total_lits;
} sat_problem_t;

sat_problem_t problem;


// Parses a cnf file
bool parse_problem (FILE* f) {

   problem.nclaus  = 0;
   problem.clauses = NULL;

   // Parse the preamble
   char buffer [2000];

   while (fgets (buffer, 1999, f)) {

      // Check if this is a comment or the problem line
      if (buffer [0] == 'c')
         continue;
      else if (buffer [0] == 'p')
         break;
      else
         return false;

   }

   char problem_type [30];
   int nlits, nclaus;

   if (sscanf (buffer, "p %s %d %d", problem_type, &nlits, &nclaus) != 3)
      return false;

   if (strcmp (problem_type, "cnf") != 0)
      return false;

   problem.total_lits = nlits;
   problem.nclaus     = nclaus;

   // Allocate space for clauses
   problem.clauses = (clause_t *) malloc( (size_t) sizeof(clause_t) * nclaus );

   for (int i = 0; i < nclaus; i++) {

      int tmp;
      int claus_lits = 0;
      problem.clauses [i].literals = NULL;

      while (fscanf (f, "%d", &tmp) && tmp != 0) {
         claus_lits++;
         problem.clauses [i].literals = (int*) realloc (problem.clauses [i].literals, sizeof (int) * claus_lits);
         problem.clauses [i].literals [claus_lits - 1] = tmp;
      }

      problem.clauses [i].nlits = claus_lits;

   }

   // Consecutive lines: clauses
   // Format: number_of_literals literal_number1 literal_number2 ....
   // -literal_number means negation

   // The problem formula is in CNF format
   /*
   int nlits;

   while (fscanf (f, "%d", &nlits) == 1) {

      // Alloc this clause on the problem struct
      problem.clauses = (clause_t*) realloc (problem.clauses, (problem.nclaus + 1) * sizeof (clause_t));
      problem.nclaus++;
      clause_t* new_claus = &problem.clauses [problem.nclaus - 1];

      // Read the new clause
      new_claus->nlits = nlits;
      new_claus->literals = (int*) malloc ((size_t) sizeof (int) * nlits );

      for (int i = 0; i < nlits; i++)
         if (fscanf (f, "%d", &new_claus->literals [i] ) != 1)
            return false;

   }
   */

   /*
   printf ("-> %d\n", problem.nclaus);

   for (int i = 0; i < problem.nclaus; i++) {

      printf ("Claus %d: %d literals (", i, problem.clauses [i].nlits);

      for (int j = 0; j < problem.clauses [i].nlits; j++)
         printf ("%d ", problem.clauses [i].literals [j]);

      printf ( ")\n" );

   }
   */

   return true;

}

// The evaluator only counts the number of clauses that are true
extern "C" double objective (GAGenome& g) {

   int count = 0;
   GA1DArrayAlleleGenome<int>* gen = dynamic_cast<GA1DArrayAlleleGenome<int>*>(&g);

   for (int i = 0; i < problem.nclaus; i++) {

      // Evaluate this clause
      int res = 0;

      for (int j = 0; (!res) && (j < problem.clauses [i].nlits); j++) {

         int lit_num = problem.clauses [i].literals [j];
         int negate = 0;

         if (lit_num < 0 )  {

            negate = 1;
            lit_num = -lit_num;

         }

         // In the problem literals start with 1 (so we can have -1 also, not possible
         // with 0). Substract 1 to index the array
         lit_num--;

         if (negate)
            res += 1 - gen->gene (lit_num);
         else 
            res += gen->gene (lit_num);

      }

      count += res;

   }

   return ((double) count) / problem.nclaus;

}

extern "C" void individualInit (GAGenome& g) {
   return IntegerUniformInitializer (g);
}

extern "C" GAGenome* defineProblem () {

   GAEDAConfig* cfg = GAEDAConfig::handle();
   FILE* fprob;

   if (cfg->getProblemData () == NULL)
      return NULL;

   fprob = fopen (cfg->getProblemData (), "r");

   if (fprob == NULL)
      return NULL;

   // Very basic parsing of the problem file
   parse_problem (fprob);

   fclose (fprob);

   // Create the genome and the allele set
   GAAlleleSet<int> alleles;
   alleles.add (0);
   alleles.add (1);
   GA1DArrayAlleleGenome<int>* genome = new GA1DArrayAlleleGenome<int> (problem.total_lits, alleles, objective);

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
   return "SAT problem";
}
