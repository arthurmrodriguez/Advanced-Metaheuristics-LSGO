#include "bbobStructures.h" /* Include all declarations for BBOB calls */

#include <string>

class CommManager;

/* include all declarations for your own optimizer here */
void MY_OPTIMIZER(double(*fitnessfunction)(double*), unsigned int dim, double ftarget, double maxfunevals, CommManager& comm_manager, std::string logFile);
CommManager* initializeMOSConfig (int argc, char**argv);
bool cleanMOSConfig (CommManager* comm_manager);
