#ifndef _CONFIG_H
#define _CONFIG_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <vector>

class Config {

   public:

      typedef enum {FIFO, MOS, SJFcpus, LJFcpus, SJFtime, LJFtime} Scheduler;

       Config () {}
      ~Config () {}

      bool parse (char* configFileName);

      const char* getJobsFile             () const {return mJobsFile;};
      const char* getMachineFile          () const {return mMachineFile;}
      const char* getClassesFile          () const {return mClassesFile;}
      const char* getFit                  () const {return mFit;}
      const char* getIndvs                () const {return mIndvs;}
      const char* getNodes                () const {return mNodes;}
      const char* getSchedulingLogPath    () const {return mSchedulingLogPath;}
      const char* getLibsSchedulerPath    () const {return mLibsSchedulerPath;}
      Scheduler   getScheduler            () const {return mScheduler;}
      const char* getWaitingVectorPath    () const {return mWaitingVectorPath;}
      const char* getExecutingQueuePath   () const {return mExecutingQueuePath;}
      const char* getStateMachinePath     () const {return mStateMachinePath;}
      const char* getStateSimulatorPath   () const {return mStateSimulatorPath;}
      unsigned    getMOSConfigs           () const {return mMOSTechs.size ();}
      const char* getSchedulingConfigPath () const {return mSchedulingConfigPath;}

      const std::vector<unsigned>& getMOSTechniques (unsigned i) const {return mMOSTechs [i];}


   private:

      char* mJobsFile;
      char* mMachineFile;
      char* mClassesFile;
      char* mFit;
      char* mIndvs;
      char* mNodes;
      char* mSchedulingConfigPath;
      char* mSchedulingLogPath;
      char* mLibsSchedulerPath;
      Scheduler mScheduler;
      char* mWaitingVectorPath;
      char* mExecutingQueuePath;
      char* mStateMachinePath;
      char* mStateSimulatorPath;

      std::vector< std::vector<unsigned> > mMOSTechs;

};

#endif
