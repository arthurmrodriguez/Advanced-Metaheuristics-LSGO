#ifndef _MOS_ALGORITHM_H
#define _MOS_ALGORITHM_H

#include "bbobStructures.h"

class CommManager;

int mos_optimizer (long double(*fitnessfunction)(long double*), CommManager& comm_manager);

#endif
