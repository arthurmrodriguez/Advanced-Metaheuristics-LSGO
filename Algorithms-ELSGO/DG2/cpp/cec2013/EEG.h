//
// Created by Arthur Rodriguez on 23/8/18.
// Adaptation for DG2: Differential Grouping Algorithm
//

#ifndef EEG_DG2
#define EEG_DG2

#include "Benchmarks.h"

class EEG:public Benchmarks{

public:

    EEG(string problem);
    ~EEG();

    // Function to load A, S, X and OptS1 matrices
    void load_data(string problem_id,string fname,
      int Dtype,int Dlength, vector< vector<double> >& dataMatrix);

    // Mean and std for each vector
    vector <double> newmeanstd(vector <double>& v);

    // Calculates the COVAR of elements
    double mycorrelation2(vector <double >& A, vector <double >& B);

    // Creates the effective Covariance matrix of X and AxS1
    vector <vector <double > > mycorrelation(vector <vector <double> >& A, vector <vector <double > >& B);

    // Calculates both diagonal and nonDiagonal elements sum
    double myDiagonal(vector < vector < double> >& M);

    // Calculates the effective C element of equation
    vector <vector <double> > MultiplyWithOutAMP(vector <vector <double> >& A, vector<vector <double> >& B);

    // Current objective function of EEG problem
    double compute(double * g);

    // DATA MEMBERS
    string problemName, PATH_TO_DATA;
    double f1max, f2max, f1min, f2min;
    int DtypeG;

    // A, S, X and OptS1 matrices
    vector < vector < double > > Mixed;
    vector < vector < double > > Amatrix;
    vector < vector < double > > ICAcomponent;
    vector < vector < double > > optimumS1;

};


#endif //EEG_DG2
