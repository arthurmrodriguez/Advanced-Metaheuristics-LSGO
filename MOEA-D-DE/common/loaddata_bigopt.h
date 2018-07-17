#ifndef LOADDARA_H
#define LOADDARA_H


#include "../common/global.h"

#include <vector>
#include <iostream>     // std::cout
#include <sstream>      // std::istringstream
#include <string>       // std::string

char    tempstring[256];
int SingMul;
std::vector < std::vector <double > > Mixed;
std::vector < std::vector <double > > Amatrix;
std::vector < std::vector <double > > ICAcomponent;

double f1max = -1000000;
double f2max = -1000000;
double f1min = 10000000;
double f2min = 10000000;
int DtypeG;
int scaling;

void printvec(std::vector < std::vector <double > > ICAfreq){
//cout << ICAfreq.at(1).size()<<"size" << ICAfreq.size()<<"-l";
	for(int i=0;i<ICAfreq.size();i++){
		for(int j=0;j<ICAfreq.at(i).size();j++){
			cout << ICAfreq.at(i).at(j) << " ";
		}
		cout << endl;
	}
}



void load_data(char * problem_id,char * fname,int Dtype,int Dlength){
	std::vector < std::vector <double > > vec;
	sprintf(tempstring, "DATAin/%s%s",problem_id, fname);
	//cout << tempstring;

	double val;
	ifstream myfile(tempstring);
	std::vector<double> temp;
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

	if(!strcmp(fname,"X.txt")){
		Mixed = vec;
	}
	else if(!strcmp(fname,"S.txt")){
		ICAcomponent = vec;
	}
	else if(!strcmp(fname,"A.txt")){
		Amatrix = vec;
	}
	else{
		cout << "wrong" << endl;
		cout << fname << endl;
	}

}


void loaddata(char *problem_id){

	int Dtype;


	if(!strcmp(problem_id,"D4")){
		Dtype=4;
	}
	else if(!strcmp(problem_id,"D4N")){
		Dtype=4;
	}
	else if(!strcmp(problem_id,"D12")){
		Dtype=12;
	}
	else if(!strcmp(problem_id,"D12N")){
		Dtype=12;
	}
	else if(!strcmp(problem_id,"D19")){
		Dtype=19;
	}
	else if(!strcmp(problem_id,"D19N")){
		Dtype=19;
	}
	DtypeG=Dtype;

	load_data(problem_id,"X.txt", Dtype,256);
	load_data(problem_id,"S.txt", Dtype,256);
	load_data(problem_id,"A.txt", Dtype,Dtype);


}

#endif
