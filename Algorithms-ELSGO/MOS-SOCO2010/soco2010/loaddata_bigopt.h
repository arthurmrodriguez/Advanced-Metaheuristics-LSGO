//
// Created by Arthur Rodriguez on 5/7/18.
//

#ifndef GAEDALIB_LOADDATA_BIGOPT_H
#define GAEDALIB_LOADDATA_BIGOPT_H

#include <vector>
#include <iostream>     // cout
#include <sstream>      // istringstream
#include <string>       // string

using namespace std;

char    tempstring[256];
int SingMul;
vector < vector <long double > > Mixed;
vector < vector <long double > > Amatrix;
vector < vector <long double > > ICAcomponent;
vector < vector <long double > > optimumS1;

long double f1max = -1000000;
long double f2max = -1000000;
long double f1min = 10000000;
long double f2min = 10000000;
int DtypeG;
int scaling;

// Put your absolute DAtain folder path in here
string PATH_TO_DATA = "/Users/Arthur18/Documents/UGR-GII/4toCurso/TFG-AdvancedMetaheuristicsForBigOpt/UGR-GII-TFG/Algorithms-ELSGO/MOS-SOCO2010/soco2010/DATAin/";

void printvec(vector < vector <long double > > ICAfreq){
//cout << ICAfreq.at(1).size()<<"size" << ICAfreq.size()<<"-l";
    for(int i=0;i<ICAfreq.size();i++){
        for(int j=0;j<ICAfreq.at(i).size();j++){
            cout << ICAfreq.at(i).at(j) << " ";
        }
        cout << endl;
    }
}



void load_data(string problem_id,string fname,int Dtype,int Dlength){
    vector < vector <long double > > vec;
    sprintf(tempstring, "%s%s%s",PATH_TO_DATA.c_str(),problem_id.c_str(), fname.c_str());

    long double val;
    fstream myfile;
    myfile.open(tempstring,ifstream::in);

    //cout<< "FILENAME IS " <<tempstring <<endl;

    vector<long double> temp;
    if (myfile.is_open())
    {


        for (int i = 0; i < Dtype; i++) {
            for (int j = 0; j < Dlength; j++) {
                myfile >> val;
                temp.push_back(val);

            }
            vec.push_back(temp);
            temp.clear();
        }
        myfile.close();
    }
    else cout << "Unable to open file -- hahaha" << fname << ".." << Dlength;

    if(fname == "X.txt"){
        Mixed = vec;
    }
    else if(fname == "S.txt"){
        ICAcomponent = vec;
    }
    else if(fname == "A.txt"){
        Amatrix = vec;
    }

    else if(fname == "S1.txt"){
        optimumS1 = vec;
    }
    else{
        cout << "wrong" << endl;
        cout << fname << endl;
    }

}


void loaddata(string problem_id){

    int Dtype;


    if(problem_id == "D4"){
        Dtype=4;
    }
    else if(problem_id == "D4N"){
        Dtype=4;
    }
    else if(problem_id == "D12"){
        Dtype=12;
    }
    else if(problem_id == "D12N"){
        Dtype=12;
    }
    else if(problem_id == "D19"){
        Dtype=19;
    }
    else if(problem_id == "D19N"){
        Dtype=19;
    }
    DtypeG=Dtype;

    // Load every file related to the problem
    load_data(problem_id,"X.txt", Dtype,256);
    load_data(problem_id,"S.txt", Dtype,256);
    load_data(problem_id,"A.txt", Dtype,Dtype);
    load_data(problem_id,"S1.txt",Dtype,256);


}


// Function to group the optimum matrix in a 1D vector
vector<long double> getOptimum(){

    vector<long double> optimum;
    for(int i = 0; i<optimumS1.size(); i++)
        for(int j = 0; j < optimumS1.at(i).size(); j++)
            optimum.push_back(optimumS1[i][j]);

    return optimum;

}


#endif //GAEDALIB_LOADDATA_BIGOPT_H
