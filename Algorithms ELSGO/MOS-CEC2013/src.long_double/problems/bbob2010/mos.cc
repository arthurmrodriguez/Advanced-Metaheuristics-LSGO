#include "mos.h"
#include "mos_algorithm.h"

#include <GAEDAConfig.h>
#include <garandom.h>
#include <MOSTechnique.h>
#include <islands/CommManager.h>

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <sstream>

#include "bbobStructures.h"

void MY_OPTIMIZER(long double(*fitnessfunction)(long double*), unsigned int dim, long double ftarget, long double maxfunevals, CommManager& comm_manager, std::string logFile)
{

  GAEDAConfig::handle()->setProblemSize(dim);
//  GAEDAConfig::handle()->setSharedEvals(maxfunevals*GAEDAConfig::handle()->getSharedEvalsPercent());
  GAEDAConfig::handle()->setSharedEvals(GAEDAConfig::handle()->getSharedEvalsFactor()*GAEDAConfig::handle()->getPopSize());
  GAEDAConfig::handle()->setEvals(maxfunevals);
  GAEDAConfig::handle()->setOptimumFitness(fgeneric_ftarget() - fgeneric_getDefaultPARAMS().precision);
  GAEDAConfig::handle()->setPrecission(fgeneric_getDefaultPARAMS().precision);
  GAEDAConfig::handle()->setLogFile(logFile.c_str());
  GAEDAConfig::handle()->setRestarts(0);
  MOSTechnique::improvement_override = false;

  mos_optimizer(fitnessfunction, comm_manager);

}


CommManager* initializeMOSConfig (int argc, char**argv) {

    CommManager* comm_manager = new CommManager(argc, argv);

    GAEDAConfig* cfg = GAEDAConfig::handle (*comm_manager);

    if (!cfg->parse (argc, argv, true)) {
       if (comm_manager->isIslandMaster ()) {
          std::cerr << "Try: " << argv[0] << " --help for more information." << std::endl;
          std::cout << "Seed: " << GAGetRandomSeed() << std::endl;
       }
       GAEDAConfig::destroy();
       exit(-1);
     }

  return comm_manager;

}

bool cleanMOSConfig (CommManager* comm_manager) {
  GAEDAConfig::destroy();
  delete comm_manager;
  return true;
}
