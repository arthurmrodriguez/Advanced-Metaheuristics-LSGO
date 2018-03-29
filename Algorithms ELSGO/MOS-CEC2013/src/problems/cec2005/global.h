/* Global variable and function definitions */

# ifndef _GLOBAL_H
# define _GLOBAL_H

/* Global Constants */
# define INF 1.79e308
# define EPS 1.0e-10
# define E  2.7182818284590452353602874713526625
# define PI 3.1415926535897932384626433832795029

/* Global variables that you are required to initialize */
extern int nreal;					/* number of real variables */
extern int nfunc;					/* number of basic functions */
extern long double bound;		/* required for plotting the function profiles for nreal=2 */
extern int density;				/* density of grid points for plotting for nreal=2 */

/* Global variables being used in evaluation of various functions */
/* These are initalized in file def2.c */
extern long double C;
extern long double global_bias;
extern long double *trans_x;
extern long double *basic_f;
extern long double *temp_x1;
extern long double *temp_x2;
extern long double *temp_x3;
extern long double *temp_x4;
extern long double *weight;
extern long double *sigma;
extern long double *lambda;
extern long double *bias;
extern long double *norm_x;
extern long double *norm_f;
extern long double **o;
extern long double **g;
extern long double ***l;

extern long double **A_f5;
extern long double *B_f5;

extern long double **A_f12;
extern long double **B_f12;
extern long double *alpha;

extern char basedir [256];

#ifdef __cplusplus
extern "C" {
#endif

/* Auxillary function declarations */
long double maximum (long double, long double);
long double minimum (long double, long double);
long double modulus (long double*, int);
long double dot (long double*, long double*, int);
long double mean (long double*, int);

/* Basic funcion declarations */
long double calc_ackley (long double*);
long double calc_rastrigin (long double*);
long double calc_weierstrass (long double*);
long double calc_griewank (long double*);
long double calc_sphere (long double*);
long double calc_schwefel (long double*);
long double calc_rosenbrock (long double *x);
long double nc_schaffer (long double, long double);
long double nc_rastrigin (long double*);

/* Utility function declarations */
void allocate_memory();
void initialize(unsigned int);
void transform (long double*, int);
void transform_norm (int);
void calc_weight (long double*);
void free_memory();

/* Benchmark function declaration */
long double calc_benchmark_func (long double*, unsigned);
void calc_benchmark_norm(unsigned);

/* Code to construct the actual file name to open */
char* getFileFullPath (const char* basedir, const char* fname, char* dest);

#ifdef __cplusplus
}
#endif

# endif
