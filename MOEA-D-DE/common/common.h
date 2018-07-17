#ifndef __COMMON_H_
#define __COMMON_H_

#include "../common/global.h"
#include "../common/mylib.h"

template <class T>
void loadpfront(char *filename, vector<T> &ps)
{
	std::fstream fin;
	int line=0;
	char str[100]=" ";
	fin.open(filename,std::ios::in);
	if(fin.is_open())
	{
		const char* pch2;
		char  a[20],b[20],c[20],d[20];
		std::string str;
		while(!fin.eof())
		{
			T  data;
			std::getline(fin,str,'\n');
			pch2 = str.c_str();
			sscanf(pch2,"%s %s %s %s",a,b,c,d);
			data.y_obj[0] = atof(a);
			data.y_obj[1] = atof(b);
			if(nobj==3)
			{
			    data.y_obj[2] = atof(c);
			}
            if(nobj==4)
			{
			    data.y_obj[3] = atof(d);
			}
			line++;
			ps.push_back(data);
		}
	} //end if
	else
		std::cout<<"failed to open "<<filename<<endl;
    fin.close();
}


template <class T>
void loadpfrontB(char *filename, vector<T> &ps)
{
	std::fstream fin;
	int line=0;
	char str[100]=" ";
	fin.open(filename,std::ios::in);
	if(fin.is_open())
	{
		const char* pch2;
		char  a[20],b[20],c[20],d[20];
		std::string str;
		while(!fin.eof())
		{
			T  data;
			std::getline(fin,str,'\n');
			pch2 = str.c_str();
			sscanf(pch2,"%s %s %s %s",a,b,c,d);
			data.indiv.y_obj[0] = atof(a);
			data.indiv.y_obj[1] = atof(b);
			if(nobj==3)
			{
			    data.indiv.y_obj[2] = atof(c);
			}
            if(nobj==4)
			{
			    data.indiv.y_obj[3] = atof(d);
			}
			line++;
			ps.push_back(data);
		}
	} //end if
	else
		std::cout<<"failed to open "<<filename<<endl;
    fin.close();
}


template <class T>
double fitnessfunction(vector <double> &y_obj, vector <double> &namda, T* ind_arr)
{
    // Chebycheff Scalarizing Function
	double fitness = 0;

	if(!strcmp(strFunctionType,"_TCHE1"))
	{
		double max_fun = -1.0e+30;
		for(int n=0; n<nobj; n++)
		{
			//double diff = fabs(y_obj[n] - idealpoint[n] + scale[n]);
			//double diff = fabs(y_obj[n] - idealpoint[n] + 0.05);
		    double diff = fabs(y_obj[n] - idealpoint[n]);
			//double diff = fabs(y_obj[n] - 0);
			double feval;
			if(namda[n]==0)
				feval = 0.0001*diff;
			else
			    feval = diff*namda[n];
			if(feval>max_fun) max_fun = feval;

		}
		fitness = max_fun;
	}

	if(!strcmp(strFunctionType,"_TCHE2"))
	{
		// reference point in the CHIM
		double max_fun = -1.0e+30;
		for(int n=0; n<nobj; n++)
		{
			double diff = (y_obj[n] - idealpoint[n])/scale[n];
			double feval;
			if(namda[n]==0)
				feval = 0.0001*diff;
			else
			    feval = diff*namda[n];
			if(feval>max_fun) max_fun = feval;

		}
		fitness = max_fun;
	}



	// CHIM + Tchebycheff
	// CHIM is not available in 3 objectives
	if(!strcmp(strFunctionType,"_NBI1")){

		// quasi normal direction
		vector <double> norm;
		for(int i=0; i<nobj; i++)
		{
		   norm.push_back(0.0);
		   for(int j=0; j<nobj; j++){
			   norm[i]+= -ind_arr[j].y_obj[i];
		   }
		}

		// normalization
		double nd = norm_vector(norm);
		for(int i=0; i<nobj; i++)
			norm[i] = norm[i]/nd;


		// reference point in the CHIM
		vector <double> base;
		for(int i=0; i<nobj; i++)
		{
			double tp2 = 0;
			for(int j=0; j<nobj; j++)
				tp2+= ind_arr[j].y_obj[i]*namda[j];
			base.push_back(tp2);
		}

		// Tchebycheff function
		double max_fun = -1.0e+30;
		for(int n=0; n<nobj; n++)
		{
			double diff  = y_obj[n] - base[n];
			double feval = -diff*norm[n];
			if(feval>max_fun) max_fun = feval;

		}
		fitness = max_fun;
	}

	//* Boundary intersection approach
	//* reference point is chosen as the ideal point
	//* the direction is independent of CHIM
	if(!strcmp(strFunctionType,"_NBI2"))
	{

		double nd = norm_vector(namda);
		for(int i=0; i<nobj; i++)
			namda[i] = namda[i]/nd;

	    // penalty method
	    // temporary vectors NBI method
		vector <double> realA(nobj);
		vector <double> realB(nobj);

		// difference beween current point and reference point
		for(int n=0; n<nobj; n++)
			realA[n] = (y_obj[n] - idealpoint[n]);

		// distance along the search direction norm
		double d1 = fabs(prod_vector(realA,namda));
		//double d1 = prod_vector(realA,norm);

		// distance to the search direction norm
		for(int n=0; n<nobj; n++)
			realB[n] = (y_obj[n] - (idealpoint[n] + d1*namda[n]));
		double d2 = norm_vector(realB);

		fitness =  (d1 + 5*d2);

		//t2 = clock();
	    //total_sec+=(t2 - t1);
	}

	// NBI method
	if(!strcmp(strFunctionType,"_NBI3")){

		// quasi normal direction
		vector <double> norm;
		for(int i=0; i<nobj; i++)
		{
		   norm.push_back(0.0);
		   for(int j=0; j<nobj; j++){
			   norm[i]+= -ind_arr[j].y_obj[i];
		   }
		}

		// normalization
		double nd = norm_vector(norm);
		for(int i=0; i<nobj; i++){
			norm[i] = norm[i]/nd;
		}


		// reference point in the CHIM
		vector <double> base;
		for(int i=0; i<nobj; i++)
		{
			double tp2 = 0;
			for(int j=0; j<nobj; j++)
				tp2+= ind_arr[j].y_obj[i]*namda[j];
			base.push_back(tp2);
		}

	    // penalty method
	    // temporary vectors NBI method
		vector <double> realA;
		vector <double> realB;

		// difference beween current point and reference point
		for(int n=0; n<nobj; n++)
			realA.push_back(y_obj[n] - base[n]);

		// distance along the search direction norm
		double d1 = prod_vector(realA,norm);

		// distance to the search direction norm
		for(int n=0; n<nobj; n++)
			realB.push_back(y_obj[n] - (base[n] + d1*norm[n]));
		double d2 = norm_vector(realB);

		fitness =  -d1 + 2*d2;
	}


	return fitness;
}

void random_permutation(int *perm, int size)
{
	int *index = new int[size];
	bool *flag = new bool[size];
    for(int n=0; n<size; n++)  {
		index[n] = n;
		flag[n]  = true;
	}

	int num = 0;
	while(num<size){
	    int start = int(size*rnd_uni(&rnd_uni_init));
		while(1){
			if(flag[start]){
				perm[num] = index[start];
				flag[start] = false;
				num++;
				break;
			}
			if(start==(size-1))
				start = 0;
			else
			    start++;
		}
	}

	delete [] index;
	delete [] flag;

}

#endif
