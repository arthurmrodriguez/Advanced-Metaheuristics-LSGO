#include "aux.h"

#include "noTrans.h"
#include "shifted_rosembrock.h"
#include "f9NoTrans.h"

CREATEHYBRIDFUNC(F19, Extended_f_10NoTrans, Shifted_Rosenbrock, 0.75, -100.0, 100.0, noTrans)
