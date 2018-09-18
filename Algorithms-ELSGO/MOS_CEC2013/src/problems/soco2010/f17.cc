#include "aux.h"

#include "../cec2008/shifted_rastrigin.h"
#include "../cec2008/shifted_rosembrock.h"

CREATEHYBRIDFUNC(F17, Shifted_Rosenbrock, Shifted_Rastrigin, 0.5, -10.0, 10.0)
