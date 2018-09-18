//
// Created by Arthur Rodriguez on 19/7/18.
// Adaptation for SHADEILS algorithm
//
#include "EEGProblem.h"

EEGProblem::EEGProblem(int dim, int noise){

    // Upper & lower bound, dimension
    // and f1-f2 min-max values
    minX = -8;
    maxX = 8;
    problemDim = dim;
    ID = 16;

    f1max = f2max = -1000000;
    f1min = f2min = 10000000;

    // Put your absolute path to DATAin
    PATH_TO_DATA = "/Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/DATAin/";

    // Check type of problem in order to load proper files
    switch (dim) {
        case 1024:
            problemName = "D4";
            Dtype = 4;
            break;
        case 3072:
            problemName = "D12";
            Dtype = 12;
            break;
        case 4864:
            problemName = "D19";
            Dtype = 19;
            break;
        default:
            problemName = "D4";
            Dtype = 4;
            break;
    }

    // Add "N" at the end in case problem allows noise
    if(noise)
        problemName += "N";
    // Load A, S, X & S1 matrices
    load_data(problemName,"A.txt", Dtype,Dtype,Amatrix);
    load_data(problemName,"S.txt", Dtype,256,ICAcomponent);
    load_data(problemName,"X.txt", Dtype,256,Mixed);
    load_data(problemName,"S1.txt",Dtype,256,optimumS1);


}

EEGProblem::~EEGProblem(){
    delete[] Ovector;
    delete[] anotherz;
    delete[] Pvector;
}

void EEGProblem::load_data(string problem_id,string fname,
  int Dtype,int Dlength, vector< vector<double> >& dataMatrix){

    // Tempstring for file name
    char tempstring[256];
    vector < vector <double > > vec;
    sprintf(tempstring, "%s%s%s",PATH_TO_DATA.c_str(),problem_id.c_str(), fname.c_str());

    double val;
    fstream myfile;
    myfile.open(tempstring,ifstream::in);

    //cout<< "FILENAME IS " <<tempstring <<endl;

    vector<double> temp;
    if (myfile.is_open())
    {
        for (int i = 0; i < Dtype; i++) {
            for (int j = 0; j < Dlength; j++) {
                myfile >> val;
                temp.push_back(val);

            }
            dataMatrix.push_back(temp);
            temp.clear();
        }
        myfile.close();
    }
    else
      cout << "Unable to open file -- hahaha" << problem_id+fname << ".." << Dlength;

}

// Function to obtain the mean and stdev of a
// set of values
vector <double>  EEGProblem::newmeanstd(vector <double >& v){
    // Accumulates the values stored in v
    vector <double> result;
    double sum = 0;//accumulate(begin(v), end(v), 0.0);
    for (unsigned long i = 0;i<v.size();i++){
        sum += v.at(i);
    }

    // Obtains the mean
    double m = sum / v.size();

    double accum = 0.0;

    // Calculates the SUM(xi - m)^2
    for (unsigned long i = 0;i<v.size();i++){
        accum += (v.at(i)-m)*(v.at(i)-m);
    }

    // returns mean and stdev
    double stdev = sqrt(accum / (v.size() - 1));
    result.push_back(m);
    result.push_back(stdev);
    return result;
}

// Calculates the COVAR of elements
double EEGProblem::mycorrelation2(vector <double >& A, vector <double >& B){

    vector <double > A1 = newmeanstd(A);
    vector <double > B1 = newmeanstd(B);

    double c1 = 0;
    double temp1, temp2;

    double a = A1.at(1)*B1.at(1);
    if (abs(a) > 0.00001){

        for (unsigned long i = 0; i < A.size(); i++){
            temp1 = ((A.at(i) - A1.at(0)) );
            temp2 = ((B.at(i) - B1.at(0)) );

            c1 += temp1*temp2;

        }
        c1 /= (A.size()*a);
        return c1;
    }
    else
        return 0;

}

// Creates the effective Covariance matrix of X and AxS1
vector < vector <double > > EEGProblem::mycorrelation(vector < vector <double > >& A, vector < vector <double > >& B){
    vector < vector <double > > M;
    vector <double > Mtemp;


    for (unsigned long i = 0; i < A.size(); i++){
        for (unsigned long j = 0; j < B.size(); j++){
            Mtemp.push_back( mycorrelation2(A.at(i), B.at(j)));
        }
        M.push_back(Mtemp);
        Mtemp.clear();
    }
    return M;
}

// Calculate diagonal and non diagonal elements of COVAR matrix
double EEGProblem::myDiagonal(vector < vector <double> >& M){

    double diagonal = 0;
    double nonDiagonal = 0;

    for (unsigned long i = 0; i < M.size(); i++){
        for (unsigned long j = 0; j < M.size(); j++){
            if (i == j){
                diagonal = diagonal + pow((1 - M.at(i).at(j)), 2);

            }

            else{
                nonDiagonal = nonDiagonal + pow(M.at(i).at(j), 2);
            }

        }
    }

    double partialDiagonal = diagonal/M.size();
    double partialNonDiagonal = nonDiagonal/M.size()/(M.size()-1);

    return partialDiagonal + partialNonDiagonal;

}

// Matrix multiplication
vector < vector <double > > EEGProblem::MultiplyWithOutAMP(vector < vector <double > >& A, vector < vector <double > >& B) {

    vector < vector <double > > C;
    vector <double > Ctemp;
    for (unsigned long row = 0; row < A.size(); row++) {
        for (unsigned long col = 0; col < B.at(0).size(); col++) {
            Ctemp.push_back(0);
        }
        C.push_back(Ctemp);
        Ctemp.clear();
    }

    for (unsigned long row = 0; row < A.size(); row++) {
        for (unsigned long col = 0; col < B.at(row).size(); col++) {

            for (unsigned long inner = 0; inner < A.at(0).size(); inner++) {
                C.at(row).at(col) += A.at(row).at(inner) * B.at(inner).at(col);
            }

        }

    }

    return C;
}

// Objective function
double EEGProblem::compute(double* x) {

    // Adapt the 1D vector to a 2D matrix
    vector < vector <double > > s1;
    vector <double > s1temp;

    // Group the genome genes in a 2D vector to represent matrix S1
    for (int i = 0; i < Dtype; i++){
        for (unsigned long j = 0; j < ICAcomponent.at(0).size(); j++){
            s1temp.push_back(x[ (i*(ICAcomponent.at(0).size())+j) ]);
        }
        s1.push_back(s1temp);
        s1temp.clear();
    }

    // Calculate X1 = AxS1
    // A is nxn and S1 is nxm so X1 is nxm too
    vector < std::vector <double > > X1 = MultiplyWithOutAMP(Amatrix, s1);

    // Calculate COR1 = covar(X, X1)
    std::vector < std::vector <double > >  COR1 = mycorrelation(X1, Mixed);
    double musum = 0;

    // Calculates f2
    for (unsigned long i = 0; i < ICAcomponent.size(); i++){

        for (unsigned long j = 0; j < ICAcomponent.at(i).size(); j++){
            musum = musum + pow((ICAcomponent.at(i).at(j)-s1.at(i).at(j)),2);
        }
    }

    // mydiagonal1 calculates the diagonal elements and mydiagonal2 calculates
    // the non-diagonal elements, so it groups in a single function called myDiagonal
    double result = myDiagonal(COR1) + musum /(ICAcomponent.size() * ICAcomponent.at(0).size());

    // Check boundaries
    if (result > f1max)
        f1max = result;
    if (result < f1min)
        f1min = result;

    return result;

}
