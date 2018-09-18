#include "F16.h"
#include <stdio.h>

/**
 * D/m-group Shifted and m-rotated Ackley’s Function
 *
 * as defined in "Benchmark Functions for the CEC'2010 Special Session
 * and Competition on Large-Scale Global Optimization" by Ke Tang,
 * Xiaodong Li, P. N. Suganthan, Zhenyu Yang, and Thomas Weise
 * published as technical report on January 8, 2010 at Nature Inspired
 * Computation and Applications Laboratory (NICAL), School of Computer
 * Science and Technology, University of Science and Technology of China,
 * Hefei, Anhui, China.
 */


F16::F16():Benchmarks(){
	m_havenextGaussian=0;
	Ovector = NULL;
	minX = -32;
	maxX = 32;
	ID = 16;

	Ovector=createShiftVector(dimension,minX,maxX);
	Pvector=createPermVector(dimension);
	RotMatrix=createRotMatrix1D(nonSeparableGroupSize);
}

F16::~F16(){
	delete[] Ovector;
	delete[] Pvector;
	delete[] RotMatrix;
}

double F16::compute(double*x){
	int i,k;
	double result=0.0;

	for(i=0;i<dimension;i++)
	{
		anotherz[i]=x[i]-Ovector[i];
	}
	for(k=1;k<=dimension/(nonSeparableGroupSize);k++)
	{
		result+=rot_ackley(anotherz,nonSeparableGroupSize,k);
	}
	return(result);
}

double F16::compute(vector<double> x){
	int i,k;
	double result=0.0;

	for(i=0;i<dimension;i++)
	{
		anotherz[i]=x[i]-Ovector[i];
	}
	for(k=1;k<=dimension/(nonSeparableGroupSize);k++)
	{
		result+=rot_ackley(anotherz,nonSeparableGroupSize,k);
	}
	return(result);
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  F16::generateInterArray
 *  Description:  
 * =====================================================================================
 */
	void
F16::generateInterArray ( )
{
	// initialize the basic structure
	for (unsigned i=0; i<(unsigned)dimension*(dimension-1)/2; i++){
		interArray.push_back(false);
	}

//	printf ( "Print P vector\n" );
//	for (unsigned i=0; i<(unsigned)dimension; i++){
//		printf ( "%d\t", Pvector[i] );
//	}

	// assign values
	unsigned baseIndex=0, compIndex=0;
	for (unsigned i=0; i<(unsigned)dimension/nonSeparableGroupSize; i++){
		for (unsigned j=0; j<(unsigned)nonSeparableGroupSize; j++){
			baseIndex =	Pvector[i*nonSeparableGroupSize+j];
			for (unsigned k=j+1; k<(unsigned)nonSeparableGroupSize; k++){
				compIndex = Pvector[i*nonSeparableGroupSize+k];
				if (baseIndex < compIndex){
//					printf ( "Mat: smallIndex %d, bigIndex %d; Arr: %d\n", baseIndex, compIndex, convertMatrixToArrayIndex(baseIndex, compIndex));
					interArray[convertMatrixToArrayIndex(baseIndex, compIndex)] = true;
				}else{
//					printf ( "Mat: smallIndex %d, bigIndex %d; Arr: %d\n", compIndex, baseIndex, convertMatrixToArrayIndex(compIndex, baseIndex));
//					printf ( "%d\n", convertMatrixToArrayIndex(compIndex, baseIndex));
					interArray[convertMatrixToArrayIndex( compIndex, baseIndex)] = true;
				}
			}
		}
	}
}		/* -----  end of function F16::generateMat  ----- */
