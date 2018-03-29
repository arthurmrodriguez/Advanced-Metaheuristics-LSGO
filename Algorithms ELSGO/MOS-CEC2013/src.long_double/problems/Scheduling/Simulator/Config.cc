#include <iostream>

#include <libconfig.h++>

#include "Config.h"

//-----------------------------------------------------------------------------------------------------------------

bool Config::parse (char* configFileName) {

   libconfig::Config cfg;

   cfg.setAutoConvert (true);

   try {

      cfg.readFile (configFileName);

   }
   catch (libconfig::ParseException e) {

      std::cerr << "Error: parsing the configuration file in line (" << e.getLine () << "): " << e.getError() << "." << std::endl;
      exit (-1);

   }
   catch (libconfig::FileIOException e) {

      std::cerr << "Error: file not found or could not be open." << std::endl;
      exit (-1);

   }


   try {

      mJobsFile    = strdup ((const char*) cfg.lookup ("JobsFile")   );
      mMachineFile = strdup ((const char*) cfg.lookup ("MachineFile"));
      mClassesFile = strdup ((const char*) cfg.lookup ("ClassesFile"));

      mFit   = strdup ((const char*) cfg.lookup ("Fit"  ));
      mIndvs = strdup ((const char*) cfg.lookup ("Indvs"));
      mNodes = strdup ((const char*) cfg.lookup ("Nodes"));

      long double fit = strtod (mFit, NULL);

      if (fit < 0.0 || fit > 1.0)
         return false;

      int indvs = atoi (mIndvs);

      if (indvs <= 0)
         return false;

      int nodes = atoi (mNodes);

      if (nodes <= 2)
         return false;

      mSchedulingConfigPath = strdup ((const char*) cfg.lookup ("SchedulingConfigPath"));

      libconfig::Setting& MOSTechs = cfg.lookup ("MOSTechniques");

      unsigned nMOSTechs = MOSTechs.getLength ();

      for (unsigned i = 0; i < nMOSTechs; i++) {

         libconfig::Setting& t = MOSTechs [i];

         std::vector<unsigned> techsv;

         for (unsigned j = 0; j < (unsigned) t.getLength (); j++)
            techsv.push_back ((unsigned) t [j]);

         mMOSTechs.push_back (techsv);

      }

      mSchedulingLogPath  = strdup ((const char*) cfg.lookup ("SchedulingLogPath" ));
      mLibsSchedulerPath  = strdup ((const char*) cfg.lookup ("LibsSchedulerPath" ));

      mWaitingVectorPath  = strdup ((const char*) cfg.lookup ("WaitingVectorPath" ));
      mExecutingQueuePath = strdup ((const char*) cfg.lookup ("ExecutingQueuePath"));
      mStateMachinePath   = strdup ((const char*) cfg.lookup ("StateMachinePath"  ));
      mStateSimulatorPath = strdup ((const char*) cfg.lookup ("StateSimulatorPath"));

      const char* scheduler = (const char*) cfg.lookup ("Scheduler");

      if (strcmp (scheduler, "mos") == 0)
         mScheduler = MOS;
      else if (strcmp (scheduler, "fifo") == 0)
         mScheduler = FIFO;
      else if (strcmp (scheduler, "sjfc") == 0)
         mScheduler = SJFcpus;
      else if (strcmp (scheduler, "ljfc") == 0)
         mScheduler = LJFcpus;
      else if (strcmp (scheduler, "sjft") == 0)
         mScheduler = SJFtime;
      else if (strcmp (scheduler, "ljft") == 0)
         mScheduler = LJFtime;
      else
         return false;

   }
   catch (libconfig::SettingNotFoundException e) {

      std::cerr << "Error: one or more mandatory parameter(s) were not found." << std::endl;
      exit (-1);

   }

   return true;

}

//-----------------------------------------------------------------------------------------------------------------
