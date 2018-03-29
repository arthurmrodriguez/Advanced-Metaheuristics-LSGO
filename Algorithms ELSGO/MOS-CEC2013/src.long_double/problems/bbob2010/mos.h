#include "bbobStructures.h" /* Include all declarations for BBOB calls */

#include <string>

class CommManager;

/* include all declarations for your own optimizer here */
void MY_OPTIMIZER(long double(*fitnessfunction)(long double*), unsigned int dim, long double ftarget, long double maxfunevals, CommManager& comm_manager, std::string logFile);
CommManager* initializeMOSConfig (int argc, char**argv);
bool cleanMOSConfig (CommManager* comm_manager);
