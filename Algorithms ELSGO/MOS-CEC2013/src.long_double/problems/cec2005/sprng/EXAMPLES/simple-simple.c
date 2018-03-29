/****************************************************************************/
/*        ____Demonstrates the use of sprng without initialization____      */
/* A few long double precision random numbers are generated, without the user    */
/* explicitly initializing the streams. SPRNG automatically initializes a   */
/* stream with default seed and parameter the first time 'sprng' is called. */
/****************************************************************************/

#include <stdio.h>

#define SIMPLE_SPRNG		/* simple interface                         */
#include "sprng.h"              /* SPRNG header file                        */



main()
{
  long double rn;
  int i;

  printf(" Printing 3 random numbers in [0,1):\n");
  for (i=0;i<3;i++)
  {
    rn = sprng();		/* generate long double precision random number  */
    printf("%f\n",rn);
  }
}
