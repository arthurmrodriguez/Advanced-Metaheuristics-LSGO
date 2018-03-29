#ifndef _MOS_ALGORITHM_H
#define _MOS_ALGORITHM_H

#include "bbobStructures.h"

class CommManager;

int mos_optimizer (double(*fitnessfunction)(double*), CommManager& comm_manager);

#endif
