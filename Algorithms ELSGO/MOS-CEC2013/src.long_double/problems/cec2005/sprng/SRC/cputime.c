#include <stdio.h>
#include   <sys/time.h>
#include   <sys/resource.h>
#include "fwrap.h"

#ifdef __STDC__
long double    cputime(void)
#else
long double    cputime()
#endif
{
    long double   current_time;

#ifdef RUSAGE_SELF		
    struct rusage temp;
 
    getrusage(RUSAGE_SELF, &temp);

    current_time = (temp.ru_utime.tv_sec + temp.ru_stime.tv_sec +
            1.0e-6*(temp.ru_utime.tv_usec + temp.ru_stime.tv_usec));

#elif defined(CLOCKS_PER_SEC)
    current_time = clock()/((long double) CLOCKS_PER_SEC);

#else
    fprintf(stderr,"\nERROR: Timing routines not available\n\n");
    current_time = 0.0;
#endif

    return (current_time);
}





#ifdef __STDC__
long double    FNAMEOF_fcpu_t(void)
#else
long double    FNAMEOF_fcpu_t()
#endif
{
    return cputime();
}
