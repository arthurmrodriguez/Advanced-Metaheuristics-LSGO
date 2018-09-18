#include <stdexcept>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <sstream>
#include <vector>
#include <iostream>
#include <iomanip>

#include <libconfig.h++>

#include <GARealOps.h>
#include <MOSGenomeFactory.h>
#include <genomes/GA1DArrayGenome.h>
#include <GAEDAConfig.h>
#include <GAPopulation.h>

#include "global.h"
#include "rand.h"

double ALLELE_MIN = 0.0;
double ALLELE_MAX = 0.0;

const double MIN_F7     =   0.0;
const double MAX_F7     = 600.0;
const double MIN_F25    =   2.0;
const double MAX_F25    =   5.0;

unsigned function = 1;
std::vector<double> opt_genome;
double opt_genome_score = 0.0;

GA1DArrayAlleleGenome<double>* genome_opt = NULL;

void initializeProblem () {

   GAEDAConfig* cfg = GAEDAConfig::handle ();

   // Extern variables from CEC 2005 aux code
   nfunc = 10;
   nreal = cfg->getProblemSize ();

	/* Call these routines to initialize random number generator */
	/* require for computing noise in some test problems */
	randomize();
   initrandomnormaldeviate ();

	/* nreal and nfunc need to be initialized before calling these routines */
	/* Routine to allocate memory to global variables */
   allocate_memory ();

	/* Routine the initalize global variables */
   if (function > 0 && function <= 25)
      initialize (function);

   switch (function) {

      case 1:
      case 2:
      case 3:
      case 4:
      case 5:
      case 6:
      case 14:
         ALLELE_MIN = -100.0;
         ALLELE_MAX =  100.0;
         break;

      case 7:
         ALLELE_MIN = -1000.0;
         ALLELE_MAX =  1000.0;
         break;

      case 8:
         ALLELE_MIN = -32.0;
         ALLELE_MAX =  32.0;
         break;

      case 9:
      case 10:
      case 13:
      case 15:
      case 16:
      case 17:
      case 18:
      case 19:
      case 20:
      case 21:
      case 22:
      case 23:
      case 24:
         ALLELE_MIN = -5.0;
         ALLELE_MAX =  5.0;
         break;

      case 11:
         ALLELE_MIN = -0.5;
         ALLELE_MAX =  0.5;
         break;

      case 12:
         ALLELE_MIN = -PI;
         ALLELE_MAX =  PI;
         break;

      case 25:
         ALLELE_MIN = -1000.0;
         ALLELE_MAX =  1000.0;
         break;

      default:
         printf("Error: The function number must be within the interval [1:25] (meaning F1-F25).\n");
         exit (-1);
         break;
   }

}

bool parse_config (char* fich) {

   libconfig::Config cfg;
   cfg.setAutoConvert (true);

   try {
      cfg.readFile (fich);
   }
   catch (libconfig::ParseException e) {
      std::cerr << "Error: parsing the configuration file in line (" << e.getLine () << "): " << e.getError() << "." << std::endl;
      exit (-1);
   }
   catch (libconfig::FileIOException e) {
      std::cerr << "Error: file not found or could not be open." << std::endl;
      exit (-1);
   }

   // First, we read mandatory parameters
   try {
      function = cfg.lookup ("function");
      strcpy (basedir, cfg.lookup ("datadir"));
   }
   catch (libconfig::SettingNotFoundException e) {
      std::cerr << "Error: one or more mandatory parameter(s) were not found." << std::endl;
      exit (-1);
   }

   return true;

}

void readAndEvalOptimumGenome () {

   std::string fname = basedir;
   fname += "global_optima.txt";

   std::ifstream f (fname.c_str ());

   char buffer [4096];

   for (unsigned i = 0; i < function; i++)
      f.getline (buffer, 4096);

   std::string str = buffer;
   istringstream ostr1 (str);

   for (unsigned i = 0; i < 100; i++) {
      double x;
      ostr1 >> x;
      opt_genome.push_back (x);
   }

}

extern "C" double objective (GAGenome& g) {

   GA1DArrayAlleleGenome<double>& genome = DYN_CAST (GA1DArrayAlleleGenome<double>&, g);
   long double x [nreal];
   double score = 0.0;

   for (int i = 0; i < nreal; i++)
      x [i] = (long double) genome.gene (i);

           calc_benchmark_norm (function);
   score = calc_benchmark_func (x, function);

   assert (score != INF);
   assert (!isnan (score));
   assert (!isinf (score));

//   return score;
   return (1 / (score + 1));

}

extern "C" void RealUniformInitializerF7 (GAGenome& g) {

   GA1DArrayAlleleGenome<double>& gen = DYN_CAST (GA1DArrayAlleleGenome<double>&, g);

   gen.resize (GAGenome::ANY_SIZE);

   for (int i = 0; i < gen.length (); i++)
      gen.gene (i, GARandomDouble (MIN_F7, MAX_F7));

}

extern "C" void RealUniformInitializerF25 (GAGenome& g) {

   GA1DArrayAlleleGenome<double>& gen = DYN_CAST (GA1DArrayAlleleGenome<double>&, g);

   gen.resize (GAGenome::ANY_SIZE);

   for (int i = 0; i < gen.length (); i++)
      gen.gene (i, GARandomDouble (MIN_F25, MAX_F25));

}

extern "C" void individualInit (GAGenome& g) {
   if (function == 7)
      return RealUniformInitializerF7 (g);
   else if (function == 25)
      return RealUniformInitializerF25 (g);
   else
      return RealUniformInitializer (g);
}

extern "C" GAGenome* defineProblem () {

   GAEDAConfig* cfg = GAEDAConfig::handle ();

   // Parse the config file for the problem
   if (!parse_config (cfg->getProblemData ()))
      throw runtime_error ("Error: A valid configuration file must be provided for the CEC 2005 problem.");

   initializeProblem ();
   readAndEvalOptimumGenome ();

   GAAlleleSet<double> alleles (ALLELE_MIN, ALLELE_MAX);
   GA1DArrayAlleleGenome<double>* genome = new GA1DArrayAlleleGenome<double> (cfg->getProblemSize (), alleles, objective);

   // Common operators
   if (function == 7)
      genome->initializer (RealUniformInitializerF7);
   else if (function == 25)
      genome->initializer (RealUniformInitializerF25);
   else
      genome->initializer (RealUniformInitializer);

   genome->comparator  (RealEuclideanComparator);

   // Specific stuff for GAs
   genome->crossover   (RealBlendCrossover);
   genome->mutator     (RealGaussianMutator);

   // Specific stuff for DE
   genome->crossover   (RealExponentialCrossover);

   // Specific stuff for MOS
   MOSGenomeFactory::handle()->registerGenome (GAID::RealEncoding, genome);

   // Evaluate Optimum genome
   genome_opt = new GA1DArrayAlleleGenome<double> (cfg->getProblemSize (), alleles, objective);
   for (int i = 0; i < cfg->getProblemSize (); i++)
      genome_opt->gene (i, opt_genome [i]);
   opt_genome_score = (1 / genome_opt->score ()) - 1;

   return genome;

}

extern "C" bool postprocess (GAPopulation* pop, int rank, int repetition) {

   if (rank == 0) {
      double best_score = (1 / pop->best ().score ()) - 1;
      std::cout << setiosflags (ios::fixed | ios::showpoint) << setprecision (10);
      std::cout << "-> BEST:  " << best_score << std::endl;
      std::cout << "-> OPTIM: " << opt_genome_score << std::endl;
      std::cout << "-> ERROR: " << best_score - opt_genome_score << std::endl;
   }

   // Routine to free the memory allocated at run time
   free_memory();

   delete genome_opt;

   return true;

}

extern "C" const char* describeProblem (void) {
   return "CEC 2005 problem.";
}
