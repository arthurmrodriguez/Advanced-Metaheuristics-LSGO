#ifndef _benchmarkshelper_H
#define _benchmarkshelper_H

long double round(long double a);
long double fmin(long double a, long double b);
long double fmax(long double a, long double b);
void unif(long double* r, int N, int inseed);
void gauss(long double * g, int N, int seed);
void computeXopt(int seed, int _DIM);
void monotoneTFosc(long double* f);
void freeStarStar(long double** M, int m);
long double** reshape(long double** B, long double* vector, int m, int n);
long double** computeRotation(long double ** B, int seed, int DIM);
long double myrand();
long double randn();
long double FGauss(long double Ftrue, long double beta);
long double FUniform(long double Ftrue, long double alpha, long double beta);
long double FCauchy(long double Ftrue, long double alpha, long double p);
int compare_long doubles (const void *a, const void *b);
void initbenchmarkshelper();
void finibenchmarkshelper();
long double computeFopt(int _funcId, int _trialId);
void setNoiseSeed(unsigned int _seed, unsigned int _seedn);

/* error handling routines - same arguments as printf, i.e. format first, then list of things to print */
/* this one exits after printing - severe error, not recoverable */
void ERROR(char *fmt, ...);
/* same, but returns to the caller, mild error */
void WARNING(char *fmt, ...);

/* Checks if sDir exists, 
   creates it if not
   checks if is writable thereafter
   Fatal ERROR if anything fails
*/
void dirOK(char *sDir);

/* create complete pathName from filename and dirname 
   is SYSTEM dependent (should be some #ifdef WINDOWS etc ...)
   fullFileName should already be allocated, at least 1024 bytes long
*/
void createFullFileName(char *fullFileName, char *dirName, char *fileName);

/* checks the existence of a file */
int existFile(char * fileName);

/* opens a file after checking it is there */
FILE * bbobOpenFile(char * fileName);

#endif
