//
// Created by Arthur Rodriguez on 23/8/18.
// Adaptation for DG2: Differential Grouping Algorithm
//

#include "EEG.h"

// Class constructor
EEG::EEG(string problem):Benchmarks(){

    Ovector = NULL;
    // Upper and lower bound for our problem
    minX = -8;
    maxX = 8;
    problemName = problem;

    f1max = f2max = -1000000;
    f1min = f2min = 10000000;

    // Put your absolute path to DATAin
    PATH_TO_DATA = "/";

    DtypeG = 4;
    // Dimension must be set properly, according
    // to problemName
    if (problemName == "D4" || problemName == "D4N"){
        dimension = 1024;
        DtypeG = 4;
        ID = 16;
    }
    else if(problemName == "D12" || problemName == "D12N"){
        dimension = 3072;
        DtypeG = 12;
        ID = 18;
    }
    else if(problemName == "D19" || problemName == "D19N"){
        dimension = 4864;
        DtypeG = 19;
        ID = 20;
    }
    else
        dimension = -1000;

    if(problemName == "D4N" || problemName == "D12N" || problemName == "D19N")
        ID+=1;

    // Effectively load problemName DATA
    // Load A, S, X & S1 matrices
    load_data(problemName,"A.txt", DtypeG,DtypeG,Amatrix);
    load_data(problemName,"S.txt", DtypeG,256,ICAcomponent);
    load_data(problemName,"X.txt", DtypeG,256,Mixed);
    load_data(problemName,"S1.txt",DtypeG,256,optimumS1);
    anotherz = new double[dimension];

}

EEG::~EEG(){
    delete[] Ovector;
    delete[] anotherz;
}

void EEG::load_data(string problem_id,string fname,
  int Dtype,int Dlength, vector< vector<double> >& dataMatrix){

    // Tempstring for file name
    char tempstring[256];
    vector < vector <double > > vec;
    sprintf(tempstring, "%s%s%s",PATH_TO_DATA.c_str(),problem_id.c_str(), fname.c_str());

    double val;
    fstream myfile;
    myfile.open(tempstring,ifstream::in);

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
      cout << "Unable to open file" << problem_id+fname << ".." << Dlength;

}


vector <double> EEG::newmeanstd(vector <double>& v){
    // Accumulates the values stored in v
    vector <double> result;
     double sum = 0;
    for (int i = 0;i<v.size();i++){
        sum += v.at(i);
    }

    // Obtains the mean
    double m = sum / v.size();

    double accum = 0.0;

    // Calculates the SUM(xi - m)^2
    for (int i = 0;i<v.size();i++){
        accum += (v.at(i)-m)*(v.at(i)-m);
    }

    // returns mean and stdev
    double stdev = sqrt(accum / (v.size() - 1));
    result.push_back(m);
    result.push_back(stdev);
    return result;
}

// Calculates the COVAR of elements
double EEG::mycorrelation2(vector <double >& A, vector <double >& B){

    vector <double > A1 = newmeanstd(A);
    vector <double > B1 = newmeanstd(B);

    double c1 = 0;
    double temp1, temp2;

    double a = A1.at(1)*B1.at(1);
    if (abs(a) > 0.00001){

        for (int i = 0; i < A.size(); i++){
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
vector <vector <double> > EEG::mycorrelation(vector <vector <double> >& A, vector <vector <double > >& B){
    vector <vector <double > > M;
    vector <double > Mtemp;


    for (int i = 0; i < A.size(); i++){
        for (int j = 0; j < B.size(); j++){
            Mtemp.push_back( mycorrelation2(A.at(i), B.at(j)));
        }
        M.push_back(Mtemp);
        Mtemp.clear();
    }
    return M;
}

// Calculates both diagonal and nonDiagonal elements sum
double EEG::myDiagonal(vector < vector < double> >& M){

    double diagonal = 0;
    double nonDiagonal = 0;

    for (int i = 0; i < M.size(); i++){
        for (int j = 0; j < M.size(); j++){
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

// Calculates the effective C element of equation
vector <vector <double> > EEG::MultiplyWithOutAMP(vector <vector <double> >& A, vector <vector <double > >& B) {


    vector <vector <double > > C;
    vector <double > Ctemp;
    for (int row = 0; row < A.size(); row++) {
        for (int col = 0; col < B.at(0).size(); col++) {
            Ctemp.push_back(0);
        }
        C.push_back(Ctemp);
        Ctemp.clear();
    }

    for (int row = 0; row < A.size(); row++) {
        for (int col = 0; col < B.at(row).size(); col++) {

            for (int inner = 0; inner < A.at(0).size(); inner++) {
                C.at(row).at(col) += A.at(row).at(inner) * B.at(inner).at(col);
            }

        }

    }

    return C;
}



// Current objective function of EEG problem
double EEG::compute(double* g) {

    // Adapt the 1D vector to a 2D matrix
    vector < vector < double > > s1;
    vector < double > s1temp;

    // Group the genome genes in a 2D vector to represent matrix S1
    for (int i = 0; i < DtypeG; i++){
        for (int j = 0; j < ICAcomponent.at(0).size(); j++){
            s1temp.push_back(g[ (i*(ICAcomponent.at(0).size())+j) ]);
        }

        s1.push_back(s1temp);
        s1temp.clear();
    }

    // Calculate X1 = AxS1
    // A is nxn and S1 is nxm so X1 is nxm too
    vector <vector <double > > X1 = MultiplyWithOutAMP(Amatrix, s1);

    // Calculate COR1 = covar(X, X1)
    vector <vector <double > >  COR1 = mycorrelation(X1, Mixed);
    double musum = 0;

    // Calculates f2
    for (int i = 0; i < ICAcomponent.size(); i++){

        for (int j = 0; j < ICAcomponent.at(i).size(); j++){
            musum = musum + pow((ICAcomponent.at(i).at(j)-s1.at(i).at(j)),2);
        }
    }

    double result = myDiagonal(COR1) + musum /(ICAcomponent.size() * ICAcomponent.at(0).size());

    // Check boundaries
    if (result > f1max)
        f1max = result;
    if (result < f1min)
        f1min = result;

    return result;

}
