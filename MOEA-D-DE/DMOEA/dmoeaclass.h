#ifndef __INDIVIDUAL_H_
#define __INDIVIDUAL_H_

#include "../common/global.h"
#include "../common/objective.h"

class CMOEADInd{
public:
	CMOEADInd();
	virtual ~CMOEADInd();

	vector <double> x_var;
	vector <double> y_obj;

	void   rnd_init();
	void   obj_eval();

    bool   operator<(const CMOEADInd &ind2);
	bool   operator<<(const CMOEADInd &ind2);
    bool   operator==(const CMOEADInd &ind2);
    void   operator=(const CMOEADInd &ind2);

	void show_objective();
	void show_variable();

	int    rank;

};

CMOEADInd::CMOEADInd()
{
	for(int i=0; i<nvar; i++)
		x_var.push_back(0.0);
	for(int n=0; n<nobj; n++)
        y_obj.push_back(0.0);
	rank = 0;
}

CMOEADInd::~CMOEADInd()
{

}

void CMOEADInd::rnd_init()
{
    for(int n=0;n<nvar;n++)
        x_var[n] = lowBound + rnd_uni(&rnd_uni_init)*(uppBound - lowBound);

}

void CMOEADInd::obj_eval()
{
	objective(x_var,y_obj);
}


void CMOEADInd::show_objective()
{
    for(int n=0; n<nobj; n++)
		printf("%f ",y_obj[n]);
	printf("\n");
}

void CMOEADInd::show_variable()
{
    for(int n=0; n<nvar; n++)
		printf("%f ",x_var[n]);
	printf("\n");
}

void CMOEADInd::operator=(const CMOEADInd &ind2)
{
    x_var = ind2.x_var;
	y_obj = ind2.y_obj;
	rank  = ind2.rank;
}

bool CMOEADInd::operator<(const CMOEADInd &ind2)
{
	bool dominated = true;
    for(int n=0; n<nobj; n++)
	{
		if(ind2.y_obj[n]<y_obj[n]) return false;
	}
	if(ind2.y_obj==y_obj) return false;
	return dominated;
}


bool CMOEADInd::operator<<(const CMOEADInd &ind2)
{
	bool dominated = true;
    for(int n=0; n<nobj; n++)
	{
		if(ind2.y_obj[n]<y_obj[n]  - 0.0001) return false;
	}
	if(ind2.y_obj==y_obj) return false;
	return dominated;
}

bool CMOEADInd::operator==(const CMOEADInd &ind2)
{
	if(ind2.y_obj==y_obj) return true;
	else return false;
}

class CSUB
{
public:
	CSUB();
	virtual ~CSUB();

	void show();

	CMOEADInd       indiv;     // best solution
	vector <int>    array;     // lattice point in a simplex
	vector <double> namda;     // weight vector
	vector <int>    table;     // neighbourhood table

	double          density, fitness;

    void  operator=(const CSUB &sub2);
};

CSUB::CSUB()
{
}

CSUB::~CSUB(){
}

void CSUB::show()
{
   for(int n=0; n<namda.size(); n++){
       printf("%f ",namda[n]);
   }
   printf("\n");
}

void CSUB::operator=(const CSUB &sub2){
    indiv  = sub2.indiv;
    array  = sub2.array;
	table  = sub2.table;
	namda  = sub2.namda;
}


#endif
