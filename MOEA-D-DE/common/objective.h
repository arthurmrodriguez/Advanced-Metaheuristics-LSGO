#ifndef __OBJECTIVE_bigopt_H_
#define __OBJECTIVE_bigopt_H_

#include "../common/global.h"
#include "../common/loaddata_bigopt.h"

vector <double> newmeanstd(vector < double > v){
	// Accumulates the values stored in v
	vector <double> result;
	double sum = 0;//std::accumulate(std::begin(v), std::end(v), 0.0);
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


double mycorrelation2(vector < double > A, vector < double > B){

	vector < double > A1 = newmeanstd(A);
	vector < double > B1 = newmeanstd(B);

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

vector < vector < double > > mycorrelation(vector < vector < double > > A, vector < vector < double > > B){
	vector < vector < double > > M;
	vector < double > Mtemp;


	for (int i = 0; i < A.size(); i++){
		for (int j = 0; j < B.size(); j++){
			Mtemp.push_back( mycorrelation2(A.at(i), B.at(j)));
		}
		M.push_back(Mtemp);
		Mtemp.clear();
	}
	return M;
}


double mydiagonal1(vector < vector < double > > M){
	double ssum = 0;

	for (int i = 0; i < M.size(); i++){
		for (int j = 0; j < M.size(); j++){
			if (i == j){
				ssum = ssum + pow((1 - M.at(i).at(j)), 2);

			}

		}
	}
	return ssum/M.size();
}

double mydiagonal2(vector < vector < double > > M){
	double ssum = 0;

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


std::vector < std::vector <double > > MultiplyWithOutAMP(std::vector < std::vector <double > > A, std::vector < std::vector <double > > B) {


	std::vector < std::vector <double > > C;
	std::vector <double > Ctemp;
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



void objective(vector<double> &x_var, vector <double> &y_obj){

	double* a = &x_var[0];

	std::vector < std::vector <double > > s1;
	std::vector <double > s1temp;

	for (int i = 0; i < DtypeG; i++){
		for (int j = 0; j < ICAcomponent.at(0).size(); j++){
			s1temp.push_back(a[i*(ICAcomponent.at(0).size())+j]);
		}
		s1.push_back(s1temp);
		s1temp.clear();
	}

	std::vector < std::vector <double > > X1;

	// Calculate X1 = AxS1
	// A is nxn and S1 is nxm so X1 is nxm too
	X1 = MultiplyWithOutAMP(Amatrix, s1);

	// Calculate COR1 = covar(X, X1)
	std::vector < std::vector <double > >  COR1 = mycorrelation(X1, Mixed);
	double musum = 0;


	for (int i = 0; i < ICAcomponent.size(); i++){

		for (int j = 0; j < ICAcomponent.at(i).size(); j++){
			musum = musum + pow((ICAcomponent.at(i).at(j)-s1.at(i).at(j)),2);
		}
	}



	switch (SingMul) {
	case 2:

		y_obj[0] = mydiagonal1(COR1) + mydiagonal2(COR1);
		y_obj[1] = musum /(ICAcomponent.size() * ICAcomponent.at(0).size());
		break;
	case 1:
		y_obj[0] = mydiagonal1(COR1) + mydiagonal2(COR1) + musum /(ICAcomponent.size() * ICAcomponent.at(0).size());
		y_obj[1] = 0;
		break;

	}


	if (y_obj[0] > f1max)
		f1max = y_obj[0];
	if (y_obj[0] < f1min)
		f1min = y_obj[0];
	if (y_obj[1] > f2max)
		f2max = y_obj[1];
	if (y_obj[1] < f2min)
		f2min = y_obj[1];
	if(scaling && SingMul==2){
		y_obj[1] = (y_obj[1]-f2min)*(f1max-f1min)/(f2max-f2min)+ f1min;
	}

}


#endif
