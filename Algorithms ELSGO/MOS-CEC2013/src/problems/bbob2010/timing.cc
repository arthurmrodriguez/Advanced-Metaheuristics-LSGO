
/* runs the timing experiment for MY_OPTIMIZER.
 */

#include <sstream>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

#include <islands/CommManager.h>

#include "bbobStructures.h"
#include "mos.h"

int main(int argc, char** argv)
{
/* your own internal stuff */

    CommManager* comm_manager = initializeMOSConfig(argc, argv);

    unsigned int dim[6] = {2, 3, 5, 10, 20, 40};
    double timings[6];
    unsigned int runs[6]; /*not really needed*/
    unsigned int dims[6]; /*not really needed*/
    unsigned int i, idx_dim, nbrun;
    double ftarget;
    clock_t t0;

    ParamStruct params = fgeneric_getDefaultPARAMS();
    strcpy(params.dataPath, "tmp");
    params.funcId = 8;
    params.instanceId = 1;

    for (idx_dim = 0; idx_dim < 6; idx_dim++)
    {
        nbrun = 0;
        params.DIM = dim[idx_dim];
        ftarget = fgeneric_initialize(params);
        t0 = clock();
        while ((double)(clock()-t0)/CLOCKS_PER_SEC < 30.)
        {
            std::stringstream strbuff;
            strbuff << "mos-f" << params.funcId << "-" << dim[idx_dim] << "dims-r" << nbrun+1;

            MY_OPTIMIZER(&fgeneric_evaluate, dim[idx_dim], ftarget, 1e5, *comm_manager, strbuff.str()); /*adjust maxfunevals*/
            nbrun ++;
        }
        timings[idx_dim] = (double)(clock()-t0)/CLOCKS_PER_SEC / (double)fgeneric_evaluations();
        dims[idx_dim] = dim[idx_dim];
        runs[idx_dim] = nbrun;
        fgeneric_finalize();
        printf("Dimensions:");
        for (i = 0; i <= idx_dim; i++)
            printf(" %11d ", dims[i]);
        printf("\n      runs:");
        for (i = 0; i <= idx_dim; i++)
            printf(" %11d ", runs[i]);
        printf("\n times [s]:");
        for (i = 0; i <= idx_dim; i++)
            printf(" %11.1e ", timings[i]);
        printf("\n");
    }

    cleanMOSConfig(comm_manager);

    return 0;
}
