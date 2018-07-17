#ifndef __INDIVIDUAL__H_
#define __INDIVIDUAL__H_

#include "../common/global.h"
#include "../common/objective.h"

class CNSGA2Ind{
public:
    vector <double> x_var;
	vector <double> y_obj;

    int    rank, count;

    void   rnd_init();
    void   obj_eval();
    void   show_objective();

    CNSGA2Ind();
	~CNSGA2Ind();

    bool   operator<(const CNSGA2Ind& ind2);
    bool   operator==(const CNSGA2Ind& ind2);
    void   operator=(const CNSGA2Ind& ind2);
};

CNSGA2Ind::CNSGA2Ind()
{

	for(int i=0; i<nvar; i++)
		x_var.push_back(0.0);
	for(int j=0; j<nobj; j++)
	    y_obj.push_back(0.0);
	rank    = 0;
}

CNSGA2Ind::~CNSGA2Ind()
{
    x_var.clear();
	y_obj.clear();
}

void CNSGA2Ind::rnd_init()
{
    for(int n=0;n<nvar;n++)
        x_var[n] = lowBound + rnd_uni(&rnd_uni_init)*(uppBound - lowBound);
    obj_eval();
}

void CNSGA2Ind::obj_eval()
{
	objective(x_var,y_obj);
}

void CNSGA2Ind::show_objective()
{

	for(int n=0;n<nobj;n++)
		std::cout<<y_obj[n]<<" ";
	std::cout<<rank<<" ";
	std::cout<<"\n";

}


bool CNSGA2Ind::operator<(const CNSGA2Ind& ind2)
{
	int flag2 = 0;
	for(int n=0;n<nobj;n++)
	{
	    if(ind2.y_obj[n] < y_obj[n])
	        return false;
		if(ind2.y_obj[n] == y_obj[n])
			flag2++;
    }

    if(flag2==nobj) return false;

	return true;
}

bool CNSGA2Ind::operator==(const CNSGA2Ind& ind2)
{
	int flag = 0;
	for(int n=0;n<nobj;n++)
	{
	    if(ind2.y_obj[n] !=y_obj[n])
	        return false;
    }
    return true;
}

void CNSGA2Ind::operator=(const CNSGA2Ind& ind2)
{
	for(int n=0;n<nobj;n++)
	    y_obj[n] = ind2.y_obj[n];

	for(int n=0;n<nvar;n++)
	    x_var[n] = ind2.x_var[n];
    rank  = ind2.rank;
}

#endif
