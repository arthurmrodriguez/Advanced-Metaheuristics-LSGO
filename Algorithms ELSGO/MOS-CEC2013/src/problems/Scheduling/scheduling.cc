#include <stdexcept>
#include <iostream>
#include <string>

#include <libconfig.h++>

#include <GAEDAConfig.h>
#include <GARealOps.h>
#include <GAIntOps.h>
#include <MOSGenomeFactory.h>
#include <MOSConversion.h>
#include <MOSConversionFunc.h>
#include <genomes/GA1DArrayGenome.h>
#include <quicksort.h>
#include <GAPopulation.h>

#include <ClusterSimulator.h>
#include <EventVector.h>

/*
/////////////////////////////////////////////////////////////
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
*/

// Configuration files for the simulator
std::string jobfile, machinefile, queuefile, simulatorfile;

// Simulation type
std::string simulType;

// Number of jobs to schedule
unsigned nJobs;

// Simulator object
ClusterSimulator* simulator;

// Number of local search iterations
unsigned lsiters = 0;

bool print_summary_info = false;

// Previous scheduling: used for PopulationInitializer
std::vector<unsigned> previousScheduling;

/*
/////////////////////////////////////////////////////////////
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
*/

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
      simulType   = (const char*) cfg.lookup ("simultype"  );
      jobfile     = (const char*) cfg.lookup ("jobfile"    );
      machinefile = (const char*) cfg.lookup ("machinefile");

      // If simulation type is partial, these files are mandatory
      if (simulType == "partial") {
         queuefile     = (const char*) cfg.lookup ("queuefile"    );
         simulatorfile = (const char*) cfg.lookup ("simulatorfile");
      }
      // If not, they are optional
      else {
         cfg.lookupValue ("queuefile"    , queuefile    );
         cfg.lookupValue ("simulatorfile", simulatorfile);
      }
   }
   catch (libconfig::SettingNotFoundException e) {

      std::cerr << "Error: one or more mandatory parameter(s) were not found." << std::endl;
      exit (-1);

   }

   // Then, we read optional vars (the default value is not modified if the
   // parameter is not provided in the config file)
   try {
      libconfig::Setting& Sched = cfg.lookup ("prevsched");

      for (unsigned i = 0; i < (unsigned) Sched.getLength (); i++)
         previousScheduling.push_back ((unsigned) Sched [i]);
   }
   catch (libconfig::SettingNotFoundException e) {
      // Ignore this exception in this case as these options are not mandatory
   }

   cfg.lookupValue ("localsearch", lsiters);

   return true;

}


double localSearch (GA1DArrayAlleleGenome<int>& g, double old_score) {

   double score = old_score;
   unsigned sz = g.size ();

   for (unsigned i = 0; i < lsiters; i++) {

      unsigned a, b, aux;

      // Randomly select two positions from the individual
      a = GARandomInt (0, sz - 2);
      b = GARandomInt (0, sz - 2);

      // Switch positions if a > b
      if (a > b) {

         aux = a;
         a   = b;
         b   = aux;

      }

      // Abort this iteration if not valid positions were selected
      if ((a == b) || (a == b + 1))
         continue;

      // New jobs order vector with exchanged positions
      std::vector<unsigned> job_order;

      for (register unsigned i = 0; i < sz; i++) {

         if (i == a)
            job_order.push_back (g.gene (b));
         else if (i == b)
            job_order.push_back (g.gene (a));
         else
            job_order.push_back (g.gene (i));

      }

      // Re-simulate and get new scheduling and its score
      simulator->resetSimulation ();

      if (simulType == "partial")
         simulator->read_from_files ((const char*) jobfile.c_str     (), (const char*) queuefile.c_str     (),
                                     (const char*) machinefile.c_str (), (const char*) simulatorfile.c_str ());
      else
         simulator->read_from_files ((const char*) jobfile.c_str     (), (const char*) machinefile.c_str   ());

      int    cpus      = simulator->cluster_machine ().total_cpus ();
      double proc_time = simulator->total_processor_time ();
      double  sche_time = simulator->simulate (job_order);
      double  penalty   = simulator->get_total_penalty ();

      double  new_score = (proc_time) / ((sche_time + penalty) * cpus);

      // If the new score beats previous one, we apply changes to the individual
      if (new_score > score) {

         unsigned c1;

         score = new_score;

         c1 = g.gene (a);
         g.gene (a, g.gene (b));
         g.gene (b, c1);


         for (a = a + 1; a < b; a++, b--) {

            c1 = g.gene (a);
            g.gene (a, g.gene (b));
            g.gene (b, c1);

         }

      }

   }

   return score;

}


extern "C" double objectiveInt (GAGenome& g) {

   GA1DArrayAlleleGenome<int>& gen = dynamic_cast<GA1DArrayAlleleGenome<int>&> (g);
   std::vector<unsigned> job_order;

   for (register int i = 0; i < gen.length (); i++)
      job_order.push_back (gen.gene (i));

   simulator->resetSimulation ();

   if (simulType == "partial")
      simulator->read_from_files ((char*) jobfile.c_str     (), (char*) queuefile.c_str     (),
                                  (char*) machinefile.c_str (), (char*) simulatorfile.c_str ());
   else
      simulator->read_from_files ((char*) jobfile.c_str     (), (char*) machinefile.c_str   ());

   int     cpus      = simulator->cluster_machine ().total_cpus ();
   double  proc_time = simulator->total_processor_time ();
   double  sche_time = simulator->simulate (job_order);
   double  penalty   = simulator->get_total_penalty ();

   double  score = (proc_time) / ((sche_time + penalty) * cpus);

   if (print_summary_info) {
      // We enter this if when called by postProcess function and thus
      // we must disable local search
      lsiters = 0;

      std::cout << "           * Processor time: " << (proc_time) / (sche_time * cpus) << std::endl;
      std::cout << "           * Scheduled time: " << sche_time << std::endl;
      std::cout << "           * Accumulated penalty: " << penalty << std::endl;
      std::cout << "           * Total: " << sche_time + penalty << std::endl;
   }

   // If we must do a local search, we do it and return the new score
   if (lsiters > 0)
      score = localSearch (gen, score);

   if (score > 1.0)
      score = 1.0;

   return score;

}


extern "C" double objectiveReal (GAGenome& g) {

   GA1DArrayAlleleGenome<double>& gen = dynamic_cast<GA1DArrayAlleleGenome<double>&> (g);
   std::vector<unsigned> job_order;

   unsigned len = gen.length ();
   int      keys   [len];
   double   values [len];

   // Prepare data
   for (register unsigned i = 0; i < len; i++) {
      keys   [i] = i;
      values [i] = gen.gene (i);
   }

   // Sort all the values
   quicksort (keys, values, 0, len - 1);

   for (register unsigned i = 0; i < len; i++)
      job_order.push_back (keys [i]);

   simulator->resetSimulation ();

   if (simulType == "partial")
      simulator->read_from_files ((char*) jobfile.c_str     (), (char*) queuefile.c_str     (),
                                  (char*) machinefile.c_str (), (char*) simulatorfile.c_str ());
   else
      simulator->read_from_files ((char*) jobfile.c_str     (), (char*) machinefile.c_str   ());

   int     cpus      = simulator->cluster_machine ().total_cpus ();
   double  proc_time = simulator->total_processor_time ();
   double  sche_time = simulator->simulate (job_order);
   double  penalty   = simulator->get_total_penalty ();

   double  score = (proc_time) / ((sche_time + penalty) * cpus);

   if (print_summary_info) {
      // We enter this if when called by postProcess function and thus
      // we must disable local search
      lsiters = 0;

      std::cout << "           * Processor time: " << (proc_time) / (sche_time * cpus) << std::endl;
      std::cout << "           * Scheduled time: " << sche_time << std::endl;
      std::cout << "           * Accumulated penalty: " << penalty << std::endl;
      std::cout << "           * Total: " << sche_time + penalty << std::endl;
   }

   if (score > 1.0)
      score = 1.0;

   return score;

}


// Comparison functions to sort jobs vectors
bool compare_sjfc (unsigned a, unsigned b) {return (simulator->get_job_cpus (a) < simulator->get_job_cpus (b));}
bool compare_ljfc (unsigned a, unsigned b) {return (simulator->get_job_cpus (a) > simulator->get_job_cpus (b));}
bool compare_sjft (unsigned a, unsigned b) {return (simulator->get_job_time (a) < simulator->get_job_time (b));}
bool compare_ljft (unsigned a, unsigned b) {return (simulator->get_job_time (a) > simulator->get_job_time (b));}
bool compare_edf  (unsigned a, unsigned b) {return (simulator->get_job_deadline (a) < simulator->get_job_deadline (b));}
bool compare_eedf (unsigned a, unsigned b) {return (simulator->get_job_effective_deadline (a) < simulator->get_job_effective_deadline (b));}


extern "C" void PopInitializer (GAPopulation& pop) {

   int popSize = pop.size ();
   int j = 0;

   // Initialize first individuals (if possible) with trivial configurations
   if (popSize >= 6) {

      GA1DArrayAlleleGenome<int>& g1 = dynamic_cast<GA1DArrayAlleleGenome<int>&> (pop.individual (0));
      GA1DArrayAlleleGenome<int>& g2 = dynamic_cast<GA1DArrayAlleleGenome<int>&> (pop.individual (1));
      GA1DArrayAlleleGenome<int>& g3 = dynamic_cast<GA1DArrayAlleleGenome<int>&> (pop.individual (2));
      GA1DArrayAlleleGenome<int>& g4 = dynamic_cast<GA1DArrayAlleleGenome<int>&> (pop.individual (3));
      GA1DArrayAlleleGenome<int>& g5 = dynamic_cast<GA1DArrayAlleleGenome<int>&> (pop.individual (4));
      GA1DArrayAlleleGenome<int>& g6 = dynamic_cast<GA1DArrayAlleleGenome<int>&> (pop.individual (5));

      unsigned genSize = g1.length ();

      // Initialize first individual with FIFO order
      std::vector<unsigned> jobs = simulator->get_job_ids ();

      for (unsigned i = 0; i < genSize; i++)
         g1.gene (i, jobs [i]);

      // Initialize second individual with SJF (cpus) order
      sort (jobs.begin (), jobs.end (), compare_sjfc);

      for (unsigned i = 0; i < genSize; i++)
         g2.gene (i, jobs [i]);

      // Initialize third individual with LJF (cpus) order
      sort (jobs.begin (), jobs.end (), compare_ljfc);

      for (unsigned i = 0; i < genSize; i++)
         g3.gene (i, jobs [i]);

      // Initialize fourth individual with SJF (time) order
      sort (jobs.begin (), jobs.end (), compare_sjft);

      for (unsigned i = 0; i < genSize; i++)
         g4.gene (i, jobs [i]);

      // Initialize fifth individual with LJF (time) order
      sort (jobs.begin (), jobs.end (), compare_ljft);

      for (unsigned i = 0; i < genSize; i++)
         g5.gene (i, jobs [i]);

      // Initialize sixth individual with previous scheduling
      if (previousScheduling.size () == genSize) {
         for (unsigned i = 0; i < genSize; i++)
            g6.gene (i, previousScheduling [i]);

         j = 6;
      }
      else {
         j = 5;
      }

   }

   // Initialize the rest of the population uniformly
   for (;j < popSize; j++)
      pop.individual (j).initialize ();

   return;

}


extern "C" GAGenome* defineProblem () {

   GAEDAConfig* cfg = GAEDAConfig::handle ();

   // Parse the config file for the problem
   if (!parse_config (cfg->getProblemData ()))
     throw runtime_error ("Error: A valid configuration file must be provided for the Scheduling problem.");

   // Create the object ClusterSimulator (we need to read files
   // to know the jobs we have to schedule and build allele set)
   simulator = new ClusterSimulator ();

   if (simulType == "partial")
      simulator->read_from_files ((char*) jobfile.c_str     (), (char*) queuefile.c_str     (),
                                  (char*) machinefile.c_str (), (char*) simulatorfile.c_str ());
   else
      simulator->read_from_files ((char*) jobfile.c_str (), (char*) machinefile.c_str ());

   // Definition of the genomes and the allele sets
   GAAlleleSet<int>    int_alleles;
   GAAlleleSet<double> real_alleles (0, 1);

   std::vector<unsigned> ids = simulator->get_job_ids ();
   nJobs = ids.size ();

   for (unsigned i = 0; i < nJobs; i++)
      int_alleles.add (ids [i]);

   GA1DArrayAlleleGenome<double>* genomeReal = new GA1DArrayAlleleGenome<double> (nJobs, real_alleles, objectiveReal);
   GA1DArrayAlleleGenome<int>*    genomeInt  = new GA1DArrayAlleleGenome<int>    (nJobs, int_alleles,  objectiveInt );

   // Common operators
   genomeInt->initializer (IntegerOrderedInitializer);
   genomeInt->comparator  (IntegerElementComparator);

   // Specific stuff for GAs
   genomeInt->crossover   (IntegerAlternativeOrderCrossover);
   genomeInt->mutator     (IntegerRepeatedExchangeMutator);

   // Specific stuff for DE
   genomeInt->crossover   (IntegerExponentialCrossover);

   // Specific stuff for MOS
   MOSGenomeFactory::handle ()->registerGenome (GAID::RealEncoding,    genomeReal);
   MOSGenomeFactory::handle ()->registerGenome (GAID::IntegerEncoding, genomeInt );

   MOSConversion::handle ()->registerConvFunction (GAID::IntegerEncoding, GAID::RealEncoding,    convertIntToReal);
   MOSConversion::handle ()->registerConvFunction (GAID::RealEncoding,    GAID::IntegerEncoding, convertRealToInt);

   return genomeInt;

}


extern "C" const char *describeProblem (void) {
   return "Scheduling.";
}


extern "C" bool postprocess (GAPopulation* pop, int rank, int repetition) {

  std::cout << pop->best ().className() << std::endl;

   if (rank == 0) {
      GA1DArrayAlleleGenome<int>& g = dynamic_cast<GA1DArrayAlleleGenome<int>&> (pop->best ());

      print_summary_info = true;

      objectiveInt (g);
   }

   return true;

}
