/* Global variables that you are required to initialize */
int nreal;					/* number of real variables */
int nfunc;					/* number of basic functions */
long double bound;		/* required for plotting the function profiles for nreal=2 */
int density;				/* density of grid points for plotting for nreal=2 */

/* Global variables being used in evaluation of various functions */
/* These are initalized in file def2.c */
long double C;
long double global_bias;
long double *trans_x;
long double *basic_f;
long double *temp_x1;
long double *temp_x2;
long double *temp_x3;
long double *temp_x4;
long double *weight;
long double *sigma;
long double *lambda;
long double *bias;
long double *norm_x;
long double *norm_f;
long double **o;
long double **g;
long double ***l;

long double **A_f5;
long double *B_f5;

long double **A_f12;
long double **B_f12;
long double *alpha;

char basedir  [ 256];
char fullName [1000];

/* Source file for custom initialization */
/* Hard-coded for every function based on the type and nature of input files */
/* At present hard-coded for D=2, 10, 30 and 50 */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

void initialize_f1()
{
    int i, j;
    FILE *fpt;
    fpt = fopen(getFileFullPath (basedir, "sphere_func_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
    }
    fclose(fpt);
    /*bias[0] = -450.0;*/
    return;
}

void initialize_f2()
{
    int i, j;
    FILE *fpt;
    fpt = fopen(getFileFullPath (basedir, "schwefel_102_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
    }
    fclose(fpt);
    /*bias[0] = -450.0;*/
    return;
}

void initialize_f3()
{
    int i, j;
    FILE *fpt;
    if (nreal==2) fpt = fopen(getFileFullPath (basedir, "elliptic_M_D2.txt", fullName),"r");
    if (nreal==10) fpt = fopen(getFileFullPath (basedir, "elliptic_M_D10.txt", fullName),"r");
    if (nreal==30) fpt = fopen(getFileFullPath (basedir, "elliptic_M_D30.txt", fullName),"r");
    if (nreal==50) fpt = fopen(getFileFullPath (basedir, "elliptic_M_D50.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nreal; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&g[i][j]);
            /* printf("\n G[%d][%d] = %LE",i+1,j+1,g[i][j]); */
        }
    }
    fclose(fpt);
    fpt = fopen(getFileFullPath (basedir, "high_cond_elliptic_rot_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
    }
    fclose(fpt);
    /*bias[0] = -450.0;*/
    return;
}

void initialize_f4()
{
    int i, j;
    FILE *fpt;
    fpt = fopen(getFileFullPath (basedir, "schwefel_102_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
    }
    fclose(fpt);
    /*bias[0] = -450.0;*/
    return;
}

void initialize_f5()
{
    int i, j;
    int index;
    FILE *fpt;
    char c;
    A_f5 = (long double **)malloc(nreal*sizeof(long double));
    for (i=0; i<nreal; i++)
    {
        A_f5[i] = (long double *)malloc(nreal*sizeof(long double));
    }
    B_f5 = (long double *)malloc(nreal*sizeof(long double));
    fpt = fopen(getFileFullPath (basedir, "schwefel_206_data.txt", fullName),"r");
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
        do
        {
            fscanf(fpt,"%c",&c);
        }
        while (c!='\n');
    }
    for (i=0; i<nreal; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&A_f5[i][j]);
            /* printf("\n A[%d][%d] = %LE",i+1,j+1,A_f5[i][j]); */
        }
        do
        {
            fscanf(fpt,"%c",&c);
        }
        while (c!='\n');
    }
    fclose(fpt);
    if (nreal%4==0)
    {
        index = nreal/4;
    }
    else
    {
        index = nreal/4 + 1;
    }
    for (i=0; i<index; i++)
    {
        o[0][i] = -100.0;
    }
    index = (3*nreal)/4 - 1;
    for (i=index; i<nreal; i++)
    {
        o[0][i] = 100.0;
    }
    for (i=0; i<nreal; i++)
    {
        B_f5[i] = 0.0;
        for (j=0; j<nreal; j++)
        {
            B_f5[i] += A_f5[i][j]*o[0][j];
        }
    }
    /*bias[0] = -310.0;*/
    return;
}

void initialize_f6()
{
    int i, j;
    FILE *fpt;
    fpt = fopen(getFileFullPath (basedir, "rosenbrock_func_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            o[i][j] -= 1.0;
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
    }
    fclose(fpt);
    bias[0] = 390.0;
    return;
}

void initialize_f7()
{
    int i, j;
    FILE *fpt;
    if (nreal==2)    fpt = fopen(getFileFullPath (basedir, "griewank_M_D2.txt", fullName),"r");
    if (nreal==10)    fpt = fopen(getFileFullPath (basedir, "griewank_M_D10.txt", fullName),"r");
    if (nreal==30)    fpt = fopen(getFileFullPath (basedir, "griewank_M_D30.txt", fullName),"r");
    if (nreal==50)    fpt = fopen(getFileFullPath (basedir, "griewank_M_D50.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nreal; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&g[i][j]);
            /* printf("\n G[%d][%d] = %LE",i+1,j+1,g[i][j]); */
        }
    }
    fclose(fpt);
    fpt = fopen(getFileFullPath (basedir, "griewank_func_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
    }
    fclose(fpt);
    /* bias[0] = -180.0; */
    return;
}

void initialize_f8()
{
    int i, j;
    int index;
    FILE *fpt;
    if (nreal==2)    fpt = fopen(getFileFullPath (basedir, "ackley_M_D2.txt", fullName),"r");
    if (nreal==10)    fpt = fopen(getFileFullPath (basedir, "ackley_M_D10.txt", fullName),"r");
    if (nreal==30)    fpt = fopen(getFileFullPath (basedir, "ackley_M_D30.txt", fullName),"r");
    if (nreal==50)    fpt = fopen(getFileFullPath (basedir, "ackley_M_D50.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nreal; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&g[i][j]);
            /* printf("\n M[%d][%d] = %LE",i+1,j+1,g[i][j]); */
        }
    }
    fclose(fpt);
    fpt = fopen(getFileFullPath (basedir, "ackley_func_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
    }
    fclose(fpt);
    index = nreal/2;
    for (i=1; i<=index; i++)
    {
        o[0][2*i-2] = -32.0;
    }
    /* bias[0] = -140.0; */
    return;
}

void initialize_f9()
{
    int i, j;
    FILE *fpt;
    fpt = fopen(getFileFullPath (basedir, "rastrigin_func_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
    }
    fclose(fpt);
    /* bias[0] = -330.0; */
    return;
}

void initialize_f10()
{
    int i, j;
    FILE *fpt;
    if (nreal==2)    fpt = fopen(getFileFullPath (basedir, "rastrigin_M_D2.txt", fullName),"r");
    if (nreal==10)    fpt = fopen(getFileFullPath (basedir, "rastrigin_M_D10.txt", fullName),"r");
    if (nreal==30)    fpt = fopen(getFileFullPath (basedir, "rastrigin_M_D30.txt", fullName),"r");
    if (nreal==50)    fpt = fopen(getFileFullPath (basedir, "rastrigin_M_D50.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nreal; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&g[i][j]);
            /* printf("\n M[%d][%d] = %LE",i+1,j+1,g[i][j]); */
        }
    }
    fclose(fpt);
    fpt = fopen(getFileFullPath (basedir, "rastrigin_func_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
    }
    fclose(fpt);
    /* bias[0] = -330.0; */
    return;
}

void initialize_f11()
{
    int i, j;
    FILE *fpt;
    if (nreal==2)    fpt = fopen(getFileFullPath (basedir, "weierstrass_M_D2.txt", fullName),"r");
    if (nreal==10)    fpt = fopen(getFileFullPath (basedir, "weierstrass_M_D10.txt", fullName),"r");
    if (nreal==30)    fpt = fopen(getFileFullPath (basedir, "weierstrass_M_D30.txt", fullName),"r");
    if (nreal==50)    fpt = fopen(getFileFullPath (basedir, "weierstrass_M_D50.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nreal; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&g[i][j]);
            /* printf("\n M[%d][%d] = %LE",i+1,j+1,g[i][j]); */
        }
    }
    fclose(fpt);
    fpt = fopen(getFileFullPath (basedir, "weierstrass_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
    }
    fclose(fpt);
    bias[0] = 90.0;
    return;
}

void initialize_f12()
{
    int i, j;
    FILE *fpt;
    char c;
    A_f12 = (long double **)malloc(nreal*sizeof(long double));
    B_f12 = (long double **)malloc(nreal*sizeof(long double));
    alpha = (long double *)malloc(nreal*sizeof(long double));
    for (i=0; i<nreal; i++)
    {
        A_f12[i] = (long double *)malloc(nreal*sizeof(long double));
        B_f12[i] = (long double *)malloc(nreal*sizeof(long double));
    }
    fpt = fopen(getFileFullPath (basedir, "schwefel_213_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    /* Reading A */
    for (i=0; i<nreal; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&A_f12[i][j]);
            /* printf("\n A[%d][%d] = %LE",i+1,j+1,A_f12[i][j]); */
        }
        do
        {
            fscanf(fpt,"%c",&c);
        }
        while (c!='\n');
    }
    if (i!=100)
    {
        for (i=nreal; i<100; i++)
        {
            do
            {
                fscanf(fpt,"%c",&c);
            }
            while(c!='\n');
        }
    }
    /* Reading B */
    for (i=0; i<nreal; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&B_f12[i][j]);
            /* printf("\n B[%d][%d] = %LE",i+1,j+1,B_f12[i][j]); */
        }
        do
        {
            fscanf(fpt,"%c",&c);
        }
        while (c!='\n');
    }
    if (i!=100)
    {
        for (i=nreal; i<100; i++)
        {
            do
            {
                fscanf(fpt,"%c",&c);
            }
            while(c!='\n');
        }
    }
    /* Reading alpha */
    for (i=0; i<nreal; i++)
    {
        fscanf(fpt,"%Lf",&alpha[i]);
        /* printf("\n alpha[%d] = %LE",i+1,alpha[i]); */
    }
    fclose(fpt);
    /* bias[0] = -460.0 */;
    return;
}

void initialize_f13()
{
    int i, j;
    FILE *fpt;
    fpt = fopen(getFileFullPath (basedir, "EF8F2_func_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            o[i][j] -= 1.0;
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
    }
    fclose(fpt);
    /* bias[0] = -130.0; */
    return;
}

void initialize_f14()
{
    int i, j;
    FILE *fpt;
    if (nreal==2)    fpt = fopen(getFileFullPath (basedir, "E_ScafferF6_M_D2.txt", fullName),"r");
    if (nreal==10)    fpt = fopen(getFileFullPath (basedir, "E_ScafferF6_M_D10.txt", fullName),"r");
    if (nreal==30)    fpt = fopen(getFileFullPath (basedir, "E_ScafferF6_M_D30.txt", fullName),"r");
    if (nreal==50)    fpt = fopen(getFileFullPath (basedir, "E_ScafferF6_M_D50.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nreal; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&g[i][j]);
            /* printf("\n M[%d][%d] = %LE",i+1,j+1,g[i][j]); */
        }
    }
    fclose(fpt);
    fpt = fopen(getFileFullPath (basedir, "E_ScafferF6_func_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
    }
    fclose(fpt);
    /* bias[0] = -300.0; */
    return;
}

void initialize_f15()
{
    int i, j;
    FILE *fpt;
    char c;
    fpt = fopen(getFileFullPath (basedir, "hybrid_func1_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
        do
        {
            fscanf(fpt,"%c",&c);
        }
        while (c!='\n');
        /* printf("\n"); */
    }
    fclose(fpt);
    lambda[0] = 1.0;
    lambda[1] = 1.0;
    lambda[2] = 10.0;
    lambda[3] = 10.0;
    lambda[4] = 1.0/12.0;
    lambda[5] = 1.0/12.0;
    lambda[6] = 5.0/32.0;
    lambda[7] = 5.0/32.0;
    lambda[8] = 1.0/20.0;
    lambda[9] = 1.0/20.0;
    global_bias = 120.0;
    return;
}

void initialize_f16()
{
    int i, j, k;
    FILE *fpt;
    char c;
    fpt = fopen(getFileFullPath (basedir, "hybrid_func1_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
        do
        {
            fscanf(fpt,"%c",&c);
        }
        while (c!='\n');
        /* printf("\n"); */
    }
    fclose(fpt);
    if (nreal==2)    fpt = fopen(getFileFullPath (basedir, "hybrid_func1_M_D2.txt", fullName),"r");
    if (nreal==10)    fpt = fopen(getFileFullPath (basedir, "hybrid_func1_M_D10.txt", fullName),"r");
    if (nreal==30)    fpt = fopen(getFileFullPath (basedir, "hybrid_func1_M_D30.txt", fullName),"r");
    if (nreal==50)    fpt = fopen(getFileFullPath (basedir, "hybrid_func1_M_D50.txt", fullName),"r");
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                fscanf(fpt,"%Lf",&l[i][j][k]);
                /* printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]); */
            }
            do
            {
                fscanf(fpt,"%c",&c);
            }
            while (c!='\n');
        }
        /* printf("\n"); */
    }
    lambda[0] = 1.0;
    lambda[1] = 1.0;
    lambda[2] = 10.0;
    lambda[3] = 10.0;
    lambda[4] = 1.0/12.0;
    lambda[5] = 1.0/12.0;
    lambda[6] = 5.0/32.0;
    lambda[7] = 5.0/32.0;
    lambda[8] = 1.0/20.0;
    lambda[9] = 1.0/20.0;
    global_bias = 120.0;
    return;
}

void initialize_f17()
{
    int i, j, k;
    FILE *fpt;
    char c;
    fpt = fopen(getFileFullPath (basedir, "hybrid_func1_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
        do
        {
            fscanf(fpt,"%c",&c);
        }
        while (c!='\n');
        /* printf("\n"); */
    }
    fclose(fpt);
    if (nreal==2)    fpt = fopen(getFileFullPath (basedir, "hybrid_func1_M_D2.txt", fullName),"r");
    if (nreal==10)    fpt = fopen(getFileFullPath (basedir, "hybrid_func1_M_D10.txt", fullName),"r");
    if (nreal==30)    fpt = fopen(getFileFullPath (basedir, "hybrid_func1_M_D30.txt", fullName),"r");
    if (nreal==50)    fpt = fopen(getFileFullPath (basedir, "hybrid_func1_M_D50.txt", fullName),"r");
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                fscanf(fpt,"%Lf",&l[i][j][k]);
                /* printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]); */
            }
            do
            {
                fscanf(fpt,"%c",&c);
            }
            while (c!='\n');
        }
        /* printf("\n"); */
    }
    lambda[0] = 1.0;
    lambda[1] = 1.0;
    lambda[2] = 10.0;
    lambda[3] = 10.0;
    lambda[4] = 1.0/12.0;
    lambda[5] = 1.0/12.0;
    lambda[6] = 5.0/32.0;
    lambda[7] = 5.0/32.0;
    lambda[8] = 1.0/20.0;
    lambda[9] = 1.0/20.0;
    global_bias = 120.0;
    return;
}

void initialize_f18()
{
    int i, j, k;
    FILE *fpt;
    char c;
    fpt = fopen(getFileFullPath (basedir, "hybrid_func2_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
        do
        {
            fscanf(fpt,"%c",&c);
        }
        while (c!='\n');
        /* printf("\n"); */
    }
    fclose(fpt);
    if (nreal==2)    fpt = fopen(getFileFullPath (basedir, "hybrid_func2_M_D2.txt", fullName),"r");
    if (nreal==10)    fpt = fopen(getFileFullPath (basedir, "hybrid_func2_M_D10.txt", fullName),"r");
    if (nreal==30)    fpt = fopen(getFileFullPath (basedir, "hybrid_func2_M_D30.txt", fullName),"r");
    if (nreal==50)    fpt = fopen(getFileFullPath (basedir, "hybrid_func2_M_D50.txt", fullName),"r");
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                fscanf(fpt,"%Lf",&l[i][j][k]);
                /* printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]); */
            }
            do
            {
                fscanf(fpt,"%c",&c);
            }
            while (c!='\n');
        }
        /* printf("\n"); */
    }
    for (i=0; i<nreal; i++)
    {
        o[nfunc-1][i] = 0.0;
    }
    sigma[0] = 1.0;
    sigma[1] = 2.0;
    sigma[2] = 1.5;
    sigma[3] = 1.5;
    sigma[4] = 1.0;
    sigma[5] = 1.0;
    sigma[6] = 1.5;
    sigma[7] = 1.5;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 5.0/16.0;
    lambda[1] = 5.0/32.0;
    lambda[2] = 2.0;
    lambda[3] = 1.0;
    lambda[4] = 1.0/10.0;
    lambda[5] = 1.0/20.0;
    lambda[6] = 20.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0/6.0;
    lambda[9] = 1.0/12.0;
    global_bias = 10.0;
    return;
}

void initialize_f19()
{
    int i, j, k;
    FILE *fpt;
    char c;
    fpt = fopen(getFileFullPath (basedir, "hybrid_func2_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
        do
        {
            fscanf(fpt,"%c",&c);
        }
        while (c!='\n');
        /* printf("\n"); */
    }
    fclose(fpt);
    if (nreal==2)    fpt = fopen(getFileFullPath (basedir, "hybrid_func2_M_D2.txt", fullName),"r");
    if (nreal==10)    fpt = fopen(getFileFullPath (basedir, "hybrid_func2_M_D10.txt", fullName),"r");
    if (nreal==30)    fpt = fopen(getFileFullPath (basedir, "hybrid_func2_M_D30.txt", fullName),"r");
    if (nreal==50)    fpt = fopen(getFileFullPath (basedir, "hybrid_func2_M_D50.txt", fullName),"r");
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                fscanf(fpt,"%Lf",&l[i][j][k]);
                /* printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]); */
            }
            do
            {
                fscanf(fpt,"%c",&c);
            }
            while (c!='\n');
        }
        /* printf("\n"); */
    }
    for (i=0; i<nreal; i++)
    {
        o[nfunc-1][i] = 0.0;
    }
    sigma[0] = 0.1;
    sigma[1] = 2.0;
    sigma[2] = 1.5;
    sigma[3] = 1.5;
    sigma[4] = 1.0;
    sigma[5] = 1.0;
    sigma[6] = 1.5;
    sigma[7] = 1.5;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 0.5/32.0;
    lambda[1] = 5.0/32.0;
    lambda[2] = 2.0;
    lambda[3] = 1.0;
    lambda[4] = 1.0/10.0;
    lambda[5] = 1.0/20.0;
    lambda[6] = 20.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0/6.0;
    lambda[9] = 1.0/12.0;
    global_bias = 10.0;
    return;
}

void initialize_f20()
{
    int i, j, k;
    int index;
    FILE *fpt;
    char c;
    fpt = fopen(getFileFullPath (basedir, "hybrid_func2_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
        do
        {
            fscanf(fpt,"%c",&c);
        }
        while (c!='\n');
        /* printf("\n"); */
    }
    fclose(fpt);
    index = nreal/2;
    for (i=1; i<=index; i++)
    {
        o[0][2*i-1] = 5.0;
    }
    if (nreal==2)    fpt = fopen(getFileFullPath (basedir, "hybrid_func2_M_D2.txt", fullName),"r");
    if (nreal==10)    fpt = fopen(getFileFullPath (basedir, "hybrid_func2_M_D10.txt", fullName),"r");
    if (nreal==30)    fpt = fopen(getFileFullPath (basedir, "hybrid_func2_M_D30.txt", fullName),"r");
    if (nreal==50)    fpt = fopen(getFileFullPath (basedir, "hybrid_func2_M_D50.txt", fullName),"r");
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                fscanf(fpt,"%Lf",&l[i][j][k]);
                /* printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]); */
            }
            do
            {
                fscanf(fpt,"%c",&c);
            }
            while (c!='\n');
        }
        /* printf("\n"); */
    }
    for (i=0; i<nreal; i++)
    {
        o[nfunc-1][i] = 0.0;
    }
    sigma[0] = 1.0;
    sigma[1] = 2.0;
    sigma[2] = 1.5;
    sigma[3] = 1.5;
    sigma[4] = 1.0;
    sigma[5] = 1.0;
    sigma[6] = 1.5;
    sigma[7] = 1.5;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 5.0/16.0;
    lambda[1] = 5.0/32.0;
    lambda[2] = 2.0;
    lambda[3] = 1.0;
    lambda[4] = 1.0/10.0;
    lambda[5] = 1.0/20.0;
    lambda[6] = 20.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0/6.0;
    lambda[9] = 1.0/12.0;
    global_bias = 10.0;
    return;
}

void initialize_f21()
{
    int i, j, k;
    FILE *fpt;
    char c;
    fpt = fopen(getFileFullPath (basedir, "hybrid_func3_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
        do
        {
            fscanf(fpt,"%c",&c);
        }
        while (c!='\n');
        /* printf("\n"); */
    }
    fclose(fpt);
    if (nreal==2)    fpt = fopen(getFileFullPath (basedir, "hybrid_func3_M_D2.txt", fullName),"r");
    if (nreal==10)    fpt = fopen(getFileFullPath (basedir, "hybrid_func3_M_D10.txt", fullName),"r");
    if (nreal==30)    fpt = fopen(getFileFullPath (basedir, "hybrid_func3_M_D30.txt", fullName),"r");
    if (nreal==50)    fpt = fopen(getFileFullPath (basedir, "hybrid_func3_M_D50.txt", fullName),"r");
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                fscanf(fpt,"%Lf",&l[i][j][k]);
                /* printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]); */
            }
            do
            {
                fscanf(fpt,"%c",&c);
            }
            while (c!='\n');
        }
        /* printf("\n"); */
    }
    sigma[0] = 1.0;
    sigma[1] = 1.0;
    sigma[2] = 1.0;
    sigma[3] = 1.0;
    sigma[4] = 1.0;
    sigma[5] = 2.0;
    sigma[6] = 2.0;
    sigma[7] = 2.0;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 1.0/4.0;
    lambda[1] = 1.0/20.0;
    lambda[2] = 5.0;
    lambda[3] = 1.0;
    lambda[4] = 5.0;
    lambda[5] = 1.0;
    lambda[6] = 50.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0/8.0;
    lambda[9] = 1.0/40.0;
    global_bias = 360.0;
    return;
}

void initialize_f22()
{
    int i, j, k;
    FILE *fpt;
    char c;
    fpt = fopen(getFileFullPath (basedir, "hybrid_func3_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
        do
        {
            fscanf(fpt,"%c",&c);
        }
        while (c!='\n');
        /* printf("\n"); */
    }
    fclose(fpt);
    if (nreal==2)    fpt = fopen(getFileFullPath (basedir, "hybrid_func3_HM_D2.txt", fullName),"r");
    if (nreal==10)    fpt = fopen(getFileFullPath (basedir, "hybrid_func3_HM_D10.txt", fullName),"r");
    if (nreal==30)    fpt = fopen(getFileFullPath (basedir, "hybrid_func3_HM_D30.txt", fullName),"r");
    if (nreal==50)    fpt = fopen(getFileFullPath (basedir, "hybrid_func3_HM_D50.txt", fullName),"r");
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                fscanf(fpt,"%Lf",&l[i][j][k]);
                /* printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]); */
            }
            do
            {
                fscanf(fpt,"%c",&c);
            }
            while (c!='\n');
        }
        /* printf("\n"); */
    }
    sigma[0] = 1.0;
    sigma[1] = 1.0;
    sigma[2] = 1.0;
    sigma[3] = 1.0;
    sigma[4] = 1.0;
    sigma[5] = 2.0;
    sigma[6] = 2.0;
    sigma[7] = 2.0;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 1.0/4.0;
    lambda[1] = 1.0/20.0;
    lambda[2] = 5.0;
    lambda[3] = 1.0;
    lambda[4] = 5.0;
    lambda[5] = 1.0;
    lambda[6] = 50.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0/8.0;
    lambda[9] = 1.0/40.0;
    global_bias = 360.0;
    return;
}

void initialize_f23()
{
    int i, j, k;
    FILE *fpt;
    char c;
    fpt = fopen(getFileFullPath (basedir, "hybrid_func3_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
        do
        {
            fscanf(fpt,"%c",&c);
        }
        while (c!='\n');
        /* printf("\n"); */
    }
    fclose(fpt);
    if (nreal==2)    fpt = fopen(getFileFullPath (basedir, "hybrid_func3_M_D2.txt", fullName),"r");
    if (nreal==10)    fpt = fopen(getFileFullPath (basedir, "hybrid_func3_M_D10.txt", fullName),"r");
    if (nreal==30)    fpt = fopen(getFileFullPath (basedir, "hybrid_func3_M_D30.txt", fullName),"r");
    if (nreal==50)    fpt = fopen(getFileFullPath (basedir, "hybrid_func3_M_D50.txt", fullName),"r");
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                fscanf(fpt,"%Lf",&l[i][j][k]);
                /* printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]); */
            }
            do
            {
                fscanf(fpt,"%c",&c);
                /*printf("\n got here \n");*/
            }
            while (c!='\n');
        }
        /* printf("\n"); */
    }
    sigma[0] = 1.0;
    sigma[1] = 1.0;
    sigma[2] = 1.0;
    sigma[3] = 1.0;
    sigma[4] = 1.0;
    sigma[5] = 2.0;
    sigma[6] = 2.0;
    sigma[7] = 2.0;
    sigma[8] = 2.0;
    sigma[9] = 2.0;
    lambda[0] = 1.0/4.0;
    lambda[1] = 1.0/20.0;
    lambda[2] = 5.0;
    lambda[3] = 1.0;
    lambda[4] = 5.0;
    lambda[5] = 1.0;
    lambda[6] = 50.0;
    lambda[7] = 10.0;
    lambda[8] = 1.0/8.0;
    lambda[9] = 1.0/40.0;
    global_bias = 360.0;
    return;
}

void initialize_f24f25()
{
    int i, j, k;
    FILE *fpt;
    char c;
    fpt = fopen(getFileFullPath (basedir, "hybrid_func4_data.txt", fullName),"r");
    if (fpt==NULL)
    {
        fprintf(stderr,"\n Error: Cannot open input file for reading \n");
        exit(0);
    }
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            fscanf(fpt,"%Lf",&o[i][j]);
            /* printf("\n O[%d][%d] = %LE",i+1,j+1,o[i][j]); */
        }
        do
        {
            fscanf(fpt,"%c",&c);
        }
        while (c!='\n');
        /* printf("\n"); */
    }
    fclose(fpt);
    if (nreal==2)    fpt = fopen(getFileFullPath (basedir, "hybrid_func4_M_D2.txt", fullName),"r");
    if (nreal==10)    fpt = fopen(getFileFullPath (basedir, "hybrid_func4_M_D10.txt", fullName),"r");
    if (nreal==30)    fpt = fopen(getFileFullPath (basedir, "hybrid_func4_M_D30.txt", fullName),"r");
    if (nreal==50)    fpt = fopen(getFileFullPath (basedir, "hybrid_func4_M_D50.txt", fullName),"r");
    for (i=0; i<nfunc; i++)
    {
        for (j=0; j<nreal; j++)
        {
            for (k=0; k<nreal; k++)
            {
                fscanf(fpt,"%Lf",&l[i][j][k]);
                /* printf("\n M[%d][%d][%d] = %LE",i+1,j+1,k+1,l[i][j][k]); */
            }
            do
            {
                fscanf(fpt,"%c",&c);
            }
            while (c!='\n');
        }
        /* printf("\n"); */
    }
    for (i=0; i<nfunc; i++)
    {
        sigma[i] = 2.0;
    }
    lambda[0] = 10.0;
    lambda[1] = 1.0/4.0;
    lambda[2] = 1.0;
    lambda[3] = 5.0/32.0;
    lambda[4] = 1.0;
    lambda[5] = 1.0/20.0;
    lambda[6] = 1.0/10.0;
    lambda[7] = 1.0;
    lambda[8] = 1.0/20.0;
    lambda[9] = 1.0/20.0;
    global_bias = 260.0;
    return;
}

void initialize (unsigned which) {

   switch (which) {

      case 1:
         initialize_f1();
         break;

      case 2:
         initialize_f2();
         break;

      case 3:
         initialize_f3();
         break;

      case 4:
         initialize_f4();
         break;

      case 5:
         initialize_f5();
         break;

      case 6:
         initialize_f6();
         break;

      case 7:
         initialize_f7();
         break;

      case 8:
         initialize_f8();
         break;

      case 9:
         initialize_f9();
         break;

      case 10:
         initialize_f10();
         break;

      case 11:
         initialize_f11();
         break;

      case 12:
         initialize_f12();
         break;

      case 13:
         initialize_f13();
         break;

      case 14:
         initialize_f14();
         break;

      case 15:
         initialize_f15();
         break;

      case 16:
         initialize_f16();
         break;

      case 17:
         initialize_f17();
         break;

      case 18:
         initialize_f18();
         break;

      case 19:
         initialize_f19();
         break;

      case 20:
         initialize_f20();
         break;

      case 21:
         initialize_f21();
         break;

      case 22:
         initialize_f22();
         break;

      case 23:
         initialize_f23();
         break;

      case 24:
      case 25:
         initialize_f24f25();
         break;

      default:
         fprintf(stderr, "Error: The function number must be within the interval [1:25] (meaning F1-F25).\n");
         break;
   }

   return;

}
