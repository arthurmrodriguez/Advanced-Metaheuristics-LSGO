#include "aux.h"

#include "noTrans.h"
#include "shifted_rastrigin.h"
#include "f9NoTrans.h"

CREATEHYBRIDFUNC(F20, Extended_f_10NoTrans, Shifted_Rastrigin, 0.75, -5.0, 5.0, noTrans)
