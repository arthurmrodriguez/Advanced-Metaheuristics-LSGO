#ifndef _BENCHMARKS_H
#define _BENCHMARKS_H

#include <vector>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include	<cmath>
#include	<ctime>
using namespace std;

#define PI (3.141592653589793238462643383279)
#define E  (2.718281828459045235360287471352)
#define L(i) ((int64_t)i)
#define D(i) ((long double)i)

struct IndexMap{
	unsigned arrIndex1;
	unsigned arrIndex2;
};

class Benchmarks{
protected:
	int next(int bits);
	int nextInt(int n);
	long double nextDouble();
	long double nextGaussian();
	unsigned ID;
	long double* createShiftVector(int dim, long double min,long double max);
	int* createPermVector(int dim);
	long double** createRotMatrix(int dim);
	long double* createRotMatrix1D(int dim);
	long double** createMultiRotateMatrix1D(int dim, int num);

	long double* lookupprepare(int dim);

	// Basic mathematical functions' declaration
	long double* multiply(long double*vector, long double*matrix,int dim);
	long double elliptic(long double*x,int dim);
	long double elliptic(long double*x, int dim, int k);
	long double rastrigin(long double*x,int dim);
	long double rastrigin(long double *x, int dim, int k); 
	long double ackley(long double*x,int dim);
	long double ackley(long double*x,int dim, int k);
	long double rot_elliptic(long double*x,int dim);
	long double rot_elliptic(long double*x,int dim, int k);
	long double rot_rastrigin(long double*x,int dim);
	long double rot_rastrigin(long double *x,int dim,int k);
	long double rot_ackley(long double*x,int dim);
	long double rot_ackley(long double*x,int dim,int k);
	long double schwefel(long double*x,int dim);
	long double schwefel(long double*x,int dim, int k);
	long double sphere(long double*x,int dim);
	long double sphere(long double*x,int dim, int k);
	long double rosenbrock(long double*x,int dim);
	long double rosenbrock(long double*x,int dim, int k);
	unsigned convertMatrixToArrayIndex ( unsigned i, unsigned j );
	void createIndexMapping (  ); 

	int64_t M;
	int64_t A;
	int64_t m_seed;
	int64_t MASK;
	long double m_nextGaussian;
	bool  m_havenextGaussian;
	bool setOvectorToZero;

	long double *Ovector;
	int*    Pvector;
	long double* RotMatrix;
	long double** MultiRotMatrix1D;
	long double *lookup;
	long double *lookup2;

	long double* anotherz;
	long double* anotherz1;
	long double* anotherz2;

	vector<bool> interArray;

	// running time setting for benchmarks
	int minX;
	int maxX;
	int dimension;
	int nonSeparableGroupSize;
	int64_t functionInitRandomSeed;
	struct IndexMap *indexMap;
	unsigned arrSize;

public:
	Benchmarks();
	virtual ~Benchmarks();
	virtual long double compute(long double* x){return 0;};
	virtual long double compute(vector<long double> x){return 0;};
	
	int getMinX();
	int getMaxX();
	unsigned getID();

	void setMinX(int);
	void setMaxX(int);
	void setSeed(int64_t);
	void setDimension(int);
	void setNonSeparableGroupSize(int);
	vector<bool> getInterArray (  );
	void ArrToMat ( unsigned I1, unsigned I2, unsigned &matIndex );
	void MatToArr ( unsigned &I1, unsigned &I2, unsigned matIndex );
};

#endif
