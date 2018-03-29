#include <stdio.h>
#include <stdlib.h>

long double myrandom_()		/* remove _ before C compilation */
{
  return (long double) rand()/RAND_MAX;
}
