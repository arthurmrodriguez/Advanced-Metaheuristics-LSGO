//
// Created by Arthur Rodriguez on 26/7/18.
// Adaptation for SHADEILS algorithm
//
#ifndef _EEGPROBLEM_
#define _EEGPROBLEM_
#include "Benchmarks.h"

using namespace std;

// Create class wrapper for EGG_Function problem
class EEGProblem : public Benchmarks{

  public:
    EEGProblem(int dim, int noise);
    ~EEGProblem();

    // Data Members
    vector < vector < double > > Mixed;
    vector < vector < double > > Amatrix;
    vector < vector < double > > ICAcomponent;
    vector < vector < double > > optimumS1;

    long double f1max,f2max, f1min, f2min;

    // Put your absolute path to DATAin here
    string PATH_TO_DATA;

    // optimimum = ICAComponent but in a single vector
    vector < double> optimum;

    // Name of the problem
    string problemName;

    // Problem Dimension
    int problemDim, Dtype;


    //_______________________________________________________
    // FUNCTIONS
    //_______________________________________________________

    // Objective function
    double compute(double* x);

    // Function to load A, S & X matrix, as well as S1 optimum
    void load_data(string problem_id,string fname,int Dtype,int Dlength, vector< vector<double> >& dataMatrix);

    // Function to obtain the mean and stdev of a
    // set of values
    vector < double> newmeanstd(vector <  double >& v);

    // Calculates the COVAR of elements
    double mycorrelation2(vector <  double >& A, vector <  double >& B);

    // Creates the effective Covariance matrix of X and AxS1
    vector < vector <  double > > mycorrelation(vector < vector <  double > >& A, vector < vector <  double > >& B);

    // Calculate diagonal and non diagonal elements of COVAR matrix
    double myDiagonal(vector < vector < double> >& M);

    // Matrix multiplication
    vector < vector < double > > MultiplyWithOutAMP(vector < vector < double > >& A, vector < vector < double > >& B);

};
#endif
