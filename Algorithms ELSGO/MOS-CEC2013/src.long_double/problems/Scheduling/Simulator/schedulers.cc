#include <iostream>
#include <vector>

// For strcmp
#include <string.h>

#include "ClusterSimulator.h"
#include "Config.h"

ClusterSimulator* simulator;

// Comparison functions to sort jobs vectors
bool compare_sjfc (unsigned a, unsigned b) {return (simulator->get_job_cpus (a) < simulator->get_job_cpus (b));}
bool compare_ljfc (unsigned a, unsigned b) {return (simulator->get_job_cpus (a) > simulator->get_job_cpus (b));}
bool compare_sjft (unsigned a, unsigned b) {return (simulator->get_job_time (a) < simulator->get_job_time (b));}
bool compare_ljft (unsigned a, unsigned b) {return (simulator->get_job_time (a) > simulator->get_job_time (b));}
bool compare_edf  (unsigned a, unsigned b) {return (simulator->get_job_deadline (a) < simulator->get_job_deadline (b));}
bool compare_eedf (unsigned a, unsigned b) {return (simulator->get_job_effective_deadline (a) < simulator->get_job_effective_deadline (b));}

// Main program
int main (int argc, char* argv []) {

   // Structure to store configuration read from config file
   Config cfg;

   // Creation of the simulator object
   simulator = new ClusterSimulator ();

   // Check for correct number of arguments
   if (argc != 3) {

      std::cerr << "Error: wrong number of arguments." << std::endl;
      std::cerr << "Usage: " << argv [0] << " config_file scheduler_type" << std::endl;
      return 1;

   }

   // Check for correct configuration file
   if (!cfg.parse (argv [1])) {

      std::cerr << "Error: wrong configuration file." << std::endl;
      return 1;

   }

   if (!simulator->parse_classes (cfg.getClassesFile())) {

      std::cerr << "Error: wrong classes file." << std::endl;
      return 1;

   }

   // Check which read function must be called
   if (!strcmp (cfg.getStateSimulatorPath(), "NULL"))
      simulator->read_from_files (cfg.getWaitingVectorPath (), cfg.getStateMachinePath   ());
   else
      simulator->read_from_files (cfg.getWaitingVectorPath (), cfg.getExecutingQueuePath (),
                                  cfg.getStateMachinePath  (), cfg.getStateSimulatorPath ());

   // Get list of jobs
   std::vector<unsigned> job_order = simulator->get_job_ids ();

   // Get values to calculate the fitness
   int    cpus      = simulator->cluster_machine ().total_cpus ();
   long double proc_time = simulator->total_processor_time ();
   float  sche_time = 0.0;

   // Output debug information and sort jobs vector (if needed)
   if (strcmp (argv [2], "FIFO") == 0) {

      std::cout << "+ FIFO scheduling: " << std::endl;

   }
   else if (strcmp (argv [2], "SJFcpus") == 0) {

      std::cout << "+ SJF (cpus) scheduling: " << std::endl;
      sort (job_order.begin (), job_order.end (), compare_sjfc);

   }
   else if (strcmp (argv [2], "LJFcpus") == 0) {

      std::cout << "+ LJF (cpus) scheduling: " << std::endl;
      sort (job_order.begin (), job_order.end (), compare_ljfc);

   }
   else if (strcmp (argv [2], "SJFtime") == 0) {

      std::cout << "+ SJF (time) scheduling: " << std::endl;
      sort (job_order.begin (), job_order.end (), compare_sjft);

   }
   else if (strcmp (argv [2], "LJFtime") == 0) {

      std::cout << "+ LJF (time) scheduling: " << std::endl;
      sort (job_order.begin (), job_order.end (), compare_ljft);

   }
   else if (strcmp (argv [2], "EDF") == 0) {

      std::cout << "+ EDF scheduling: " << std::endl;
      sort (job_order.begin (), job_order.end (), compare_edf);

   }
   else if (strcmp (argv [2], "EEDF") == 0) {

      std::cout << "+ EEDF scheduling: " << std::endl;
      sort (job_order.begin (), job_order.end (), compare_eedf);

   }


   // Simulate with jobs vector with the appropriate simulator
   if (strcmp (argv [2], "BACK") == 0) {

      std::cout << "+ Backfilling: " << std::endl;
      sche_time = simulator->simulate_backfilling (job_order);

   }
   else if (strcmp (argv [2], "BACK_RES") == 0) {

      std::cout << "+ Backfilling with reservation: " << std::endl;
      sche_time = simulator->simulate_backfilling_res (job_order);

   }
   else {

      sche_time = simulator->simulate (job_order);

   }

   float penalty = simulator->get_total_penalty ();

   // Output fitness and new order for jobs
   std::cout << "Solution: " << (proc_time) / ((sche_time + penalty) * cpus) << std::endl;
   std::cout << "* Integer coding: ";

   for (unsigned i = 0; i < job_order.size (); i++)
      std::cout << job_order [i] << " ";

   std::cout << std::endl;

   std::cout << "* Processor time: " << (proc_time) / (sche_time * cpus) << std::endl;
   std::cout << "* Scheduled time: " << sche_time << std::endl;
   std::cout << "* Accumulated penalty: " << penalty << std::endl;
   std::cout << "* Total: " << sche_time + penalty << std::endl;

   return 0;

}
