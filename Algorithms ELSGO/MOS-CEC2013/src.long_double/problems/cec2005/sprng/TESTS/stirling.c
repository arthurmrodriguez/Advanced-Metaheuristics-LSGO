#include <stdio.h>
#include "util.h"


#ifdef __STDC__
long double stirling(int n, int m)
#else
long double stirling(n, m)
int n, m;
#endif
{
  long double **s, temp;
  int i, j;
  
  if(n < m)
    return 0.0;
  else if(n == m)
    return 1.0;
  else if(m == 0 && n > 0)
    return 0.0;
  else if(m == 1 && n > 0)
    return 1.0;
  
  s = (long double **) mymalloc(m*sizeof(long double *));
  for(i=0; i<m; i++)
    s[i] = (long double *) mymalloc((n-m+1)*sizeof(long double));
  
  for(j=0; j<=n-m; j++)
    s[0][j] = 1.0;
  
  for(i=1; i<m; i++)
  {
    s[i][0] = 1.0;
    for(j=1; j<=n-m; j++)
      s[i][j] = (i+1)*s[i][j-1] + s[i-1][j];
  }

  temp = s[m-1][n-m];
  
  for(i=0; i<m; i++)
    free(s[i]);
  
  return temp;
}


#if 0
main(int argc, char *argv[])
{
  int n, m;
  
  n = atoi(argv[1]);
  m = atoi(argv[2]);
  
  printf("S(%d,%d) = %f\n", n, m, stirling(n,m));
}
#endif

