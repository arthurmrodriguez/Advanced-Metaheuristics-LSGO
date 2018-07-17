/*==========================================================================
//  Implementation of MOEA/D Based on Differential Evolution (DE) for Continuous Multiobjective
//  Optimization Problems with Complicate Pareto Sets (2007)
//
//  See the details of MOEA/D-DE and test problems in the following paper
//  H. Li and Q. Zhang, Comparison Between NSGA-II and MOEA/D on a Set of Multiobjective Optimization
//  Problems with Complicated Pareto Sets, Technical Report CES-476, Department of Computer Science,
//  University of Essex, 2007
//
//  The component functions of each test instance can be found in "objective.h".
//
//  The source code of MOEA/D-DE and NSGA-II-DE were implemented by Hui Li and Qingfu Zhang
//
//  If you have any questions about the codes, please contact
//  Qingfu Zhang at qzhang@essex.ac.uk  or Hui Li at hzl@cs.nott.ac.uk
===========================================================================*/


#include "common/global.h"
#include "DMOEA/dmoeafunc.h"
#include "NSGA2/nsga2func.h"
#include "common/loaddata_bigopt.h"


void execute(char *alg);

void savescale(){
	std::fstream fout;
	fout.open("Scale.txt",std::ios::out);

	fout<<f1min<<endl;
	fout<<f1max<<endl;
	fout<<f2min<<endl;
	fout<<f2max<<endl;

	fout.close();
}

int main()
{

	//single objective =1, multiobjective =2
	SingMul = 2;

	//do you use scaling 0 don't use, 1 use
	scaling =1;

	//which data instance are you working on
	loaddata("D4");

	//cout << DtypeG;

	// The settings of test instances F1-F9
	// char *ins[] = {"F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9"}
	int pof[]      = {  21,   21,   21,   21,   21,   31,   21,   21,  22};       // type of pareto front
	int pos[]      = {  21,   22,   23,   24,   26,   32,   21,   21,  22};       // type of pareto set
	int dis[]      = {   1,    1,    1,    1,    1,    1,    3,    4,   1};       // type of non-negative function
	//int var[]      = {  30,   30,   30,   30,   30,   10,   10,   10,  30};       // dimensionality of search space
	//int obj[]      = {   2,    2,    2,    2,    2,    3,    2,    2,   2};       // number of objective functions

	int var[]      = {DtypeG*256};
	int obj[]      = {2};

	// The settings of algorithms
	int pop[] = { 20 };     // population size
	int gen[] = { 50 };     // number of generations

	for(int i=0; i<1; i++)
	{
		// the parameter setting of test instance
		dtype = dis[i];
		ptype = pof[i];
		ltype = pos[i];
		nvar  = var[i];
		nobj  = obj[i];


		// the parameter setting of MOEA/D-DE and
		pops    = pop[i];
		max_gen = gen[i];

		sprintf(strTestInstance,"P%dD%dL%d",ptype, dtype, ltype);
		//printf("Instances: pf shape %d  - distance %d, - ps shape %d \n ", ptype, dtype, ltype);
		//execute("DMOEA");

		execute("DMOEA");
	}
	savescale();

	int kin;
	cin>> kin;
}

void execute(char *alg)
{
	std::fstream fout;
	char filename[1024];
	// compute IGD-values
	if(strcmp(alg,"DMOEA")==0)
	    sprintf(filename,"GD/GD_DMOEA_%s.dat",strTestInstance);
	else
	    sprintf(filename,"GD/GD_NSGA2_%s.dat",strTestInstance);

	fout.open(filename,std::ios::out);

	for(int run=1; run<=max_run; run++)
	{
		vector<double> gd;
		if(strcmp(alg,"DMOEA")==0)
		{
			CMOEAD  MOEAD;
			gd = MOEAD.execute(run, "_TCHE1", "_DE");

		}
		else
		{
		    CNSGA2  NSGA2;
		    gd = NSGA2.execute(run);

		}
		for(int k=0; k<gd.size(); k++)
			fout<<gd[k]<<" ";
		fout<<"\n";
		gd.clear();
	}
	fout.close();
}
