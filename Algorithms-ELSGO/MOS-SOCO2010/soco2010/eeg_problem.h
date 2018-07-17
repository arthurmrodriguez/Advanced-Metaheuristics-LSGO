//
// Created by Arthur Rodriguez on 5/7/18.
//

#ifndef GAEDALIB_EEG_PROBLEM_H
#define GAEDALIB_EEG_PROBLEM_H

#include <genomes/GA1DArrayGenome.h>
#include "loaddata_bigopt.h"

// optimimum = ICAComponent but in a single vector
vector <long double> optimum;

// Name of the problem
string problemName;

// Function to obtain the mean and stdev of a
// set of values
vector <long double> newmeanstd(vector < long double > v){
    // Accumulates the values stored in v
    vector <long double> result;
    long double sum = 0;//accumulate(begin(v), end(v), 0.0);
    for (int i = 0;i<v.size();i++){
        sum += v.at(i);
    }

    // Obtains the mean
    long double m = sum / v.size();

    long double accum = 0.0;

    // Calculates the SUM(xi - m)^2
    for (int i = 0;i<v.size();i++){
        accum += (v.at(i)-m)*(v.at(i)-m);
    }

    // returns mean and stdev
    long double stdev = sqrt(accum / (v.size() - 1));
    result.push_back(m);
    result.push_back(stdev);
    return result;
}

// Calculates the COVAR of elements
long double mycorrelation2(vector < long double > A, vector < long double > B){

    vector < long double > A1 = newmeanstd(A);
    vector < long double > B1 = newmeanstd(B);

    long double c1 = 0;
    long double temp1, temp2;

    long double a = A1.at(1)*B1.at(1);
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
vector < vector < long double > > mycorrelation(vector < vector < long double > > A, vector < vector < long double > > B){
    vector < vector < long double > > M;
    vector < long double > Mtemp;


    for (int i = 0; i < A.size(); i++){
        for (int j = 0; j < B.size(); j++){
            Mtemp.push_back( mycorrelation2(A.at(i), B.at(j)));
        }
        M.push_back(Mtemp);
        Mtemp.clear();
    }
    return M;
}

long double mydiagonal1(vector < vector < long double > > M){
    long double ssum = 0;

    for (int i = 0; i < M.size(); i++){
        for (int j = 0; j < M.size(); j++){
            if (i == j){
                ssum = ssum + pow((1 - M.at(i).at(j)), 2);

            }

        }
    }
    return ssum/M.size();
}

long double mydiagonal2(vector < vector < long double > > M){
    long double ssum = 0;

    for (int i = 0; i < M.size(); i++){
        for (int j = 0; j < M.size(); j++){
            if (i == j){

            }
            else{
                ssum = ssum + pow(M.at(i).at(j), 2);

            }
        }
    }
    return ssum/M.size()/(M.size()-1);

}

long double myDiagonal(vector < vector <long double> > M){

    long double diagonal = 0;
    long double nonDiagonal = 0;

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

    long double partialDiagonal = diagonal/M.size();
    long double partialNonDiagonal = nonDiagonal/M.size()/(M.size()-1);

    return partialDiagonal + partialNonDiagonal;

}


vector < vector <long double > > MultiplyWithOutAMP(vector < vector <long double > > A, vector < vector <long double > > B) {


    vector < vector <long double > > C;
    vector <long double > Ctemp;
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




long double eeg_Function(GAGenome& g) {

    GA1DArrayAlleleGenome<long double>& genome = DYN_CAST (GA1DArrayAlleleGenome<long double>&, g);

    // Adapt the 1D vector to a 2D matrix
    vector < vector <long double > > s1;
    vector <long double > s1temp;

    // Group the genome genes in a 2D vector to represent matrix S1
    for (int i = 0; i < DtypeG; i++){
        for (int j = 0; j < ICAcomponent.at(0).size(); j++){
            s1temp.push_back(genome.gene((unsigned int)(i*(ICAcomponent.at(0).size())+j)));
        }
        s1.push_back(s1temp);
        s1temp.clear();
    }

    // Calculate X1 = AxS1
    // A is nxn and S1 is nxm so X1 is nxm too
    vector < std::vector <long double > > X1 = MultiplyWithOutAMP(Amatrix, s1);

    // Calculate COR1 = covar(X, X1)
    std::vector < std::vector <long double > >  COR1 = mycorrelation(X1, Mixed);
    long double musum = 0;

    // Calculates f2
    for (int i = 0; i < ICAcomponent.size(); i++){

        for (int j = 0; j < ICAcomponent.at(i).size(); j++){
            musum = musum + pow((ICAcomponent.at(i).at(j)-s1.at(i).at(j)),2);
        }
    }

    // mydiagonal1 calculates the diagonal elements and mydiagonal2 calculates
    // the non-diagonal elements, so it groups in a single function called myDiagonal
    long double result = myDiagonal(COR1) + musum /(ICAcomponent.size() * ICAcomponent.at(0).size());

    // Check boundaries
    if (result > f1max)
        f1max = result;
    if (result < f1min)
        f1min = result;

    return result;

}


#endif //GAEDALIB_EEG_PROBLEM_H
