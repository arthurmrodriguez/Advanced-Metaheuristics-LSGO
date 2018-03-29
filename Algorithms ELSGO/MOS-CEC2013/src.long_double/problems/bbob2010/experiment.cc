/* runs an entire experiment for benchmarking MY_OPTIMIZER
* on the noise-free testbed
* or the noisy testbed (change the ifun loop in this case as given below).
*/

#include <sstream>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

#include <GAEDAConfig.h>
#include <islands/CommManager.h>

#include "bbobStructures.h" /* Include all declarations for BBOB calls */
#include "mos.h"

int main(int argc, char** argv)
{

    CommManager* comm_manager = initializeMOSConfig(argc, argv);

    // Original maximum number of FEs (it is overwritten in the inner loop, so we need a copy)
    unsigned origMax = GAEDAConfig::handle()->getEvals();

    unsigned int dim[6] = {2, 3, 5, 10, 20, 40};
    unsigned int instances[15] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    unsigned int idx_dim, ifun, idx_instances, seed;
    int independent_restarts;
    long double maxfunevals, minfunevals;

    clock_t t0 = clock();
    time_t Tval;
    /**************************************************
     *          BBOB Mandatory initialization         *
     *************************************************/
    /* retrieve all default parameters of BBOB calls  */
    ParamStruct params = fgeneric_getDefaultPARAMS();

    /* modify the following parameters, choosing a different setting
     * for each new experiment */
    strcpy(params.dataPath, (GAEDAConfig::handle()->getBBOBOutputPath() + "/" + GAEDAConfig::handle()->getBBOBConfigName()).c_str());  /* different folder for each experiment! */
    /* please beforehand run from the command-line 'python createfolders.py PUT_MY_BBOB_DATA_PATH'
     * to create the necessary folder structure to run an experiment. */
    strcpy(params.algName, GAEDAConfig::handle()->getBBOBConfigName().c_str());
    strcpy(params.comments, "PUT MORE DETAILED INFORMATION, PARAMETER SETTINGS ETC");

    seed = time(NULL);
    srand(seed); /* used by MY_OPTIMIZER */
    printf("random seed set to %d\n", seed);

    /* To make the noise deterministic. */
    /* fgeneric_noiseseed(30); printf("seed for the noise set to: 30\n"); */

    /* now the main loop */
    for (idx_dim = 0; idx_dim < 6; idx_dim++)
    {
        /* Function indices are from 1 to 24 (noiseless) or from 101 to 130 (noisy) */
        /* for the noisy functions exchange the for loop with */
        /* for (ifun = 101; ifun <= 130; ifun++) */
        int auxFuncNumber = GAEDAConfig::handle()->getFunction();

        if ((auxFuncNumber < 1 || auxFuncNumber > 24) && (auxFuncNumber < 101 || auxFuncNumber > 130)) {
          std::cerr << "Error: Wrong function number: " << auxFuncNumber << ". Aborting..." << std::endl;
          exit(-1);
        }

        for (ifun = auxFuncNumber; ifun <= auxFuncNumber; ifun++)
        {
            for (idx_instances = 0; idx_instances < 15; idx_instances++)
            {
                /* set DIM, funcId, instanceId to initialize BBOB fgeneric */
                params.DIM = dim[idx_dim];
                params.funcId = ifun;
                params.instanceId = instances[idx_instances];
                /* call the BBOB initialization */
                fgeneric_initialize(params);

                /* now call your optimizer so that it optimizes the function
                 * fgeneric_evaluate or
                 * fgeneric_evaluate_vector(long double * XX, unsigned int howMany,
                 *                          long double * result)
                 */

                /* The fgeneric interface can give some information:
                 *    e.g. fgeneric_ftarget() the target value only for termination
                 *         fgeneric_evaluations() the number of calls to fgeneric_evaluate 
                 *                 after fgeneric_initialization
                 *         fgeneric_best() the best value reached  
                 */
                maxfunevals = origMax * dim[idx_dim]; /* PUT APPROPRIATE MAX. NUMBER OF FEVALS */
                                                      /* 5. * dim should be fine to just check everything */
                minfunevals = dim[idx_dim] + 2;       /* PUT MINIMAL USEFUL NUMBER OF FEVALS */
                independent_restarts = 0;
                while (fgeneric_evaluations() + minfunevals <= maxfunevals)
                {
                    std::stringstream strbuff;
                    strbuff << GAEDAConfig::handle()->getBBOBOutputPath() + "/" + GAEDAConfig::handle()->getBBOBConfigName() + "/mos/mos-f" << ifun << "-" << dim[idx_dim] << "dims-r" << idx_instances+1;
                    MY_OPTIMIZER(&fgeneric_evaluate, dim[idx_dim], fgeneric_ftarget(),
                                 maxfunevals - fgeneric_evaluations(), *comm_manager, strbuff.str());
                    independent_restarts += GAEDAConfig::handle()->getRestarts();
                    if (fgeneric_best() < fgeneric_ftarget())
                        break;
                }

                printf("  f%d in %d-D, instance %d: FEs=%.0f with %d restarts,", ifun, dim[idx_dim],
                       instances[idx_instances], fgeneric_evaluations(), independent_restarts);
                printf(" fbest-ftarget=%.4e, elapsed time [h]: %.2f, ftarget=%.4e\n", 
                       fgeneric_best() - fgeneric_ftarget(), (long double)(clock()-t0)/CLOCKS_PER_SEC/60./60.,fgeneric_ftarget());
                /* call the BBOB closing function to wrap things up neatly */
                fgeneric_finalize();
            }
            Tval = time(NULL);
            printf("    date and time: %s", ctime(&Tval));
        }
        printf("---- dimension %d-D done ----\n", dim[idx_dim]);
    }

    cleanMOSConfig(comm_manager);

    return 0;
}
