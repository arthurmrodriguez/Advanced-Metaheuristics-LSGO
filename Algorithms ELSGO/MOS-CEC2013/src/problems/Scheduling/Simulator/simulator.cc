// For 'fork', 'exec', 'open', etc.
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>

// For streams
#include <iostream>
#include <sstream>
#include <fstream>

// For storing purposes
#include <vector>

// Project files
#include "Config.h"
#include "EventQueue.h"
#include "EventVector.h"
#include "ClusterSimulator.h"

// To associate a stream with a file
#include "fdstream.h"

// Aux functions
std::vector<unsigned> schedule (char** params);
bool createMOSConfigFile (const Config& cfg, const ClusterSimulator& simul, unsigned nMOS);

//-----------------------------------------------------------------------------------------------------------------

int main (int argc, char* argv []) {

   // Structure to store configuration read from config file
   Config cfg;

   // Check for correct number of arguments
   if (argc < 2) {

      std::cerr << "Error: wrong number of arguments." << std::endl;
      std::cerr << "Usage: " << argv [0] << " config_file" << std::endl;
      return 1;

   }

   // Check for correct configuration file
   if (!cfg.parse (argv [1])) {

      std::cerr << "Error: wrong configuration file." << std::endl;
      return 1;

   }

   // Parameters for different schedulers
   char* params_ga   [] = {(char*) "mpirun", (char*) "-np", (char*) cfg.getNodes (), (char*) "/home/alum/atorre/bin/mosexec", (char*) "-U", (char*) cfg.getFit (), (char*) "-p", (char*) cfg.getIndvs (), (char*) "-t", (char*) "mesh", (char*) "-m", (char*) "async", (char*) "-a", (char*) "mos", (char*) "-f", (char*) "10", (char*) "-A", (char*) cfg.getSchedulingConfigPath (), (char*) "-T", (char*) cfg.getSchedulingLogPath (), (char*) cfg.getLibsSchedulerPath (), NULL};
   char* params_fifo [] = {(char*) "./schedulers", argv[1], (char*) "FIFO", NULL};
   char* params_sjfc [] = {(char*) "./schedulers", argv[1], (char*) "SJFcpus", NULL};
   char* params_ljfc [] = {(char*) "./schedulers", argv[1], (char*) "LJFcpus", NULL};
   char* params_sjft [] = {(char*) "./schedulers", argv[1], (char*) "SJFtime", NULL};
   char* params_ljft [] = {(char*) "./schedulers", argv[1], (char*) "LJFtime", NULL};
   char* params_edf  [] = {(char*) "./schedulers", argv[1], (char*) "EDF", NULL};
   char* params_eedf [] = {(char*) "./schedulers", argv[1], (char*) "EEDF", NULL};

   // Vars controlling if a replanification must be done
   // (replanning_todo) and when it must be done (step)
   bool replanning_todo, step;

   // Aux vars to store events
   SimulatorEvent* event;
   EventVector* new_events;

   // Event queue to store jobs
   EventQueue event_queue;
   event_queue.read_from_file (cfg.getJobsFile (), in_event);

   // Simulator object
   ClusterSimulator* simulator = new ClusterSimulator (cfg.getMachineFile (), cfg.getClassesFile());
   //simulator->set_classes_file(cfg.getClassesFile ());

   // Accounting of MOS re-schedulings
   unsigned counter = 0;

   // Main loop
   while (!event_queue.empty () || simulator->num_events_waiting () > 0) {

      // Initialization of variables
      replanning_todo = false;
      step            = false;

      // Iterate while events continue arriving and none job has still finished
      while (!event_queue.empty () && !step) {

         // Get job at the top of the queue (next one)
         event = event_queue.first ();

         if (event->type () == in_event) { // An incoming job

            // Debug output
            std::cerr << "+ Step: current_time = " << event->m_job->timestamp () << std::endl;
            std::cerr << "+ Number of jobs executing: " << simulator->num_events_running () << std::endl;
            std::cerr << "+ Number of jobs waiting: " << simulator->num_events_waiting () << std::endl;
            std::cerr << "   + Searching jobs:" << std::endl;
            std::cerr << "     => Pop event of type 'IN':" << std::endl << "     "; event->print (std::cerr); std::cerr << std::endl;

            // Pop incoming event
            event_queue.remove ();

            // Add it to the waiting queue of the simulator
            simulator->add_event (event);

            // As a new job arrived, we must mark a replanification to be done
            // when the next executing job finishes
            replanning_todo = true;

            // Try to execute some of the waiting jobs
            new_events = simulator->step ();

            // Add returned output events associated to jobs that have
            // begun their execution after the call to 'step'
            while (!new_events->empty ()) {

               SimulatorEvent* e = new_events->first ();
               event_queue.add (e);
               new_events->remove ();

            }

            // Avoid memory leaks
            delete new_events;
            new_events = NULL;

         }
         else if (event->type () == out_event) { // A finished job

            // Debug output
            std::cerr << "+ Step: current_time = " << event->end_time () << std::endl;
            std::cerr << "+ Start penalty = " << event->penalty () << std::endl;
            std::cerr << "+ Number of jobs executing: " << simulator->num_events_running () << std::endl;
            std::cerr << "   + Searching jobs:" << std::endl;
            std::cerr << "     => Pop event of type 'OUT':" << std::endl << "     "; event->print (std::cerr); std::cerr << std::endl;

            // Pop finished job event and free its memory
            event_queue.remove ();
            delete event;

            // Finish job within the simulator and free its resources
            simulator->run_to_next_event ();

            // Mark the need to perform a re-schedule of waiting jobs
            step = true;

         }

      }

      // Check if a re-schedule must be done and if there are more
      // than one job in the waiting queue
      if (replanning_todo && simulator->num_events_waiting () > 1) {

         // Debug output
         std::cerr << "     => Re-scheduling." << std::endl;

         // Print the state of the simulator
         simulator->print_to_files (cfg.getWaitingVectorPath (), cfg.getExecutingQueuePath (),
                                    cfg.getStateMachinePath  (), cfg.getStateSimulatorPath ());

         // Always compute trivial schedulers
         std::cerr << "        + FIFO scheduling: " << std::endl;
         schedule (params_fifo);

         std::cerr << "        + SJF (cpus) scheduling: " << std::endl;
         schedule (params_sjfc);

         std::cerr << "        + LJF (cpus) scheduling: " << std::endl;
         schedule (params_ljfc);

         std::cerr << "        + SJF (time) scheduling: " << std::endl;
         schedule (params_sjft);

         std::cerr << "        + LJF (time) scheduling: " << std::endl;
         schedule (params_ljft);

         std::cerr << "        + EDF scheduling: " << std::endl;
         schedule (params_edf);

         std::cerr << "        + EEDF scheduling: " << std::endl;
         schedule (params_eedf);

         // Check if we have to launch a MOS scheduler
         if (cfg.getScheduler () == Config::MOS) {

            for (unsigned i = 0; i < cfg.getMOSConfigs (); i++) {

               char buff [32], logNameStr [1024];

               if (simulator->num_events_waiting () > 75)
                  sprintf (buff, "75");
               else if (simulator->num_events_waiting () > 25)
                  sprintf (buff, "%d", simulator->num_events_waiting ());
               else
                  sprintf (buff, "25");

               params_ga [7] = buff;

               std::stringstream logName;
               logName << cfg.getSchedulingConfigPath () << "-run" << counter++ << ".cfg";
               sprintf (logNameStr, "%s", logName.str ().c_str ());

               params_ga [19] = logNameStr;

               createMOSConfigFile (cfg, *simulator, i);

               // Debug output
               std::cerr << "        + MOS scheduling [" << i << "]: " << std::endl;

               std::vector<unsigned> order = schedule (params_ga);

               // Store new schedule into the simulator
               if (i == 0)
                  simulator->set_order (order);

            }

         }

         // Debug output
         std::cerr << "     => Re-scheduling done." << std::endl << std::endl;

      }

      // Try to execute some of the waiting jobs
      new_events = simulator->step ();

      // Add returned output events associated to jobs that have
      // begun their execution after the call to 'step'
      while (!new_events->empty ()) {

         SimulatorEvent* e = new_events->first ();
         event_queue.add (e);
         new_events->remove ();

      }

   }

   // simulator->print ();
   std::cout << "Final execution time: " << simulator->current_time () << std::endl;
   std::cout << "New time: " << simulator->m_new_time << std::endl; // Debug
   std::cout << "Total penalty: " << simulator->get_total_penalty () << std::endl;

   return 0;

}


// Executes the scheduler whose parameters are provided
std::vector<unsigned> schedule (char** params) {

   // Order returned by the scheduler
   std::vector<unsigned> order;

   // Vars needed by fork
   int pid = 0, status = 0;
   int fds [2];

   pipe (fds);

   // Fork new process to execute the scheduler
   pid = fork ();

   // Error in fork
   if (pid == -1) {

      std::cerr << "Error creating child process. Aborting..." << std::endl;
      exit (1);

   }
   else if (pid == 0) { // Child process

      // Redirect standard output of child process to the pipe
      close (1);
      dup (fds [1]);

      // Discard standard error output
      //int fp = open ("/dev/null", O_WRONLY);
      close (2);
      //dup (fp);
      dup (fds [1]);

      // Close unused descriptors
      close (fds [0]);

      // Launch scheduler process
      execvp (params [0], params);

   }
   else { // Parent process

      close (fds [1]);

      // Wait for scheduler termination
      wait (&status);

      // Associate the pipe with a stream buffer
      boost::fdistream f (fds [0]);

      if (f.bad ())
         std::cerr << "Error reading from pipe." << std::endl;
      else {

         // Buffer to store read lines
         char line [8192];

         // Execution time (for MOS)
         double exec_time = 0.0;

         // Parse output of the child (scheduler) process
         while (!f.eof ()) {

            // Read line from the stream associated to the pipe
            f.getline (line, 8192);

            // Convert read line into a string stream
            std::stringstream buff (line);

            // Aux var to extract discarded tokens
            std::string discard;

            // Fitness of the returned schedule
            double fitness;

            std::cout << line << std::endl;

            // Parse output of the child process
            if (buff.str ().find ("RES:") != std::string::npos) {

               // Discard two tokens ('->' and 'RES:') and then get the fitness
               buff >> discard;
               buff >> discard;
               buff >> fitness;

               std::cerr << "          * Fitness: " << fitness << std::endl;

            }
            if (buff.str ().find ("TIME:") != std::string::npos) {

               // Discard two tokens ('->' and 'TIME:') and then get the fitness
               buff >> discard;
               buff >> discard;
               buff >> exec_time;

            }
            else if (buff.str ().find ("Solution:") != std::string::npos) {

               // Discard one token ('Solution:') and then get the fitness
               buff >> discard;
               buff >> fitness;

               std::cerr << "          * Fitness: " << fitness << std::endl;

            }
            else if (buff.str ().find ("Integer coding") != std::string::npos) {

               // Discard first 3 items: '*', 'Integer' and 'coding:'
               for (unsigned i = 0; i < 3; i++)
                  buff >> discard;

               // Read job ids and store them in the order vector
               int item = 0, old = -1;

               while (1) {

                  old = item;
                  buff >> item;

                  if (old != item)
                     order.push_back (item);
                  else
                     break;

               }

               // Output received schedule
               std::cerr << "          * Best schedule found: ";

               for (unsigned i = 0; i < order.size (); i++)
                  std::cerr << order [i] << " ";

               std::cerr << std::endl;

            }
            else if (buff.str ().find ("Processor") != std::string::npos) {

               double proc_time;

               // Discard three tokens ('*', 'Processor' and 'time:') and then get the processor time
               buff >> discard;
               buff >> discard;
               buff >> discard;
               buff >> proc_time;

               std::cerr << "          * % of Processor time: " << proc_time << std::endl;

            }
            else if (buff.str ().find ("Scheduled") != std::string::npos) {

               double sche_time;

               // Discard three tokens ('*', 'Scheduled' and 'time:') and then get the scheduled time
               buff >> discard;
               buff >> discard;
               buff >> discard;
               buff >> sche_time;

               std::cerr << "          * Scheduled time: " << sche_time << std::endl;

            }
            else if (buff.str ().find ("Accumulated") != std::string::npos) {

               double penalty;

               // Discard three tokens ('*', 'Accumulated' and 'penalty:') and then get the penalty
               buff >> discard;
               buff >> discard;
               buff >> discard;
               buff >> penalty;

               std::cerr << "          * Accumulated penalty: " << penalty << std::endl;

            }
            else if (buff.str ().find ("Total") != std::string::npos) {

               double total;

               // Discard two tokens ('*' and 'total:') and then get the total time
               buff >> discard;
               buff >> discard;
               buff >> total;

               std::cerr << "          * Total (Scheduled time + Penalty time): " << total << std::endl;

               break;

            }
            else if (buff.str ().find ("failed") != std::string::npos ||
                     buff.str ().find ("Error" ) != std::string::npos    ) {

               std::cerr << " Error: the scheduler failed to execute. Aborting..." << std::endl;
               exit (-1);
               break;

            }
            else {

               // Ignore them if not debugging
               //std::cerr << line << std::endl;

            }

         }

         if (exec_time != 0.0)
            std::cerr << "          * Execution time: " << exec_time << std::endl;

      }

   }

   return order;

}


bool createMOSConfigFile (const Config& cfg, const ClusterSimulator& simul, unsigned nMOS) {

   std::ofstream of (cfg.getSchedulingConfigPath ());

   of << "simultype = \"partial\";" << std::endl;

   of << "jobfile = \""       << cfg.getWaitingVectorPath  () << "\";" << std::endl;
   of << "machinefile = \""   << cfg.getStateMachinePath   () << "\";" << std::endl;
   of << "queuefile = \""     << cfg.getExecutingQueuePath () << "\";" << std::endl;
   of << "simulatorfile = \"" << cfg.getStateSimulatorPath () << "\";" << std::endl;

   of << "techs = [";

   const std::vector<unsigned>& techs = cfg.getMOSTechniques (nMOS);

   for (unsigned i = 0; i < techs.size (); i++)
      if (i != (techs.size () - 1))
         of << techs [i] << ", ";
      else
         of << techs [i] << "];" << std::endl;

   of << "alg = \"mos\";"   << std::endl;
   of << "localsearch = 0;" << std::endl;

   const std::vector<unsigned> waitingJobs = simul.get_job_ids ();

   of << "prevsched = [";

   for (unsigned i = 0; i < waitingJobs.size (); i++)
      if (i != (waitingJobs.size () - 1))
         of << waitingJobs [i] << ", ";
      else
         of << waitingJobs [i] << "];" << std::endl;

   of.close ();

   return true;

}
