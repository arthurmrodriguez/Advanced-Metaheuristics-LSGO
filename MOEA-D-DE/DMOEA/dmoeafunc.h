#ifndef __EVOLUTION_H_
#define __EVOLUTION_H_

#include "../common/global.h"
#include "../common/recombination.h"
#include "../common/common.h"
#include "dmoeaclass.h"


class CMOEAD
{

public:
	CMOEAD();
	virtual ~CMOEAD();


	void init_uniformweight();               // initialize the weights for subproblems
	void init_neighbourhood();               // calculate the neighbourhood of each subproblem
	void init_population();                  // initialize the population


	void update_reference(CMOEADInd &ind);                 // update ideal point which is used in Tchebycheff or NBI method
	void update_problem(CMOEADInd &ind, int id, int type); // update current solutions in the neighbourhood

	void diffevolution();                                      // DE-based recombination
	void matingselection(vector<int> &list, int cid, int size, int type);  // select mating parents

	// execute MOEAD
	vector<double> execute(int run, char *strfunc, char *stralg);

    void read_front(char filename[1024]);
    void calc_distance();

	void save_front(char savefilename[1024]);       // save the pareto front into files
	void save_ps(char savefilename[1024]);

    vector <CSUB>       population;
	vector <CMOEADInd>  ps;
	vector <CSUB>       offspring;
	vector <int>        array;
	CMOEADInd           *ind_arr;

	double              distance;                   // generational distance from PF to solutions found
	int                 popsize;

	void operator=(const CMOEAD &moea);
};

CMOEAD::CMOEAD()
{

	ind_arr = new CMOEADInd[nobj];

	// initialize ideal point
    for(int n=0; n<nobj; n++)
	{
		idealpoint.push_back(1.0e+30);
		ind_arr[n].rnd_init();
		ind_arr[n].obj_eval();
	}
}

CMOEAD::~CMOEAD()
{
	idealpoint.clear();
	delete [] ind_arr;
}


void CMOEAD::init_population()
{

    for(int i=0; i<population.size(); i++)
	{
		population[i].indiv.rnd_init();
		population[i].indiv.obj_eval();
		update_reference(population[i].indiv);
		nfes++;
	}
}

void CMOEAD::operator=(const CMOEAD &alg)
{

	population = alg.population;
	ps         = alg.ps;
	ind_arr    = alg.ind_arr;
	offspring  = alg.offspring;
	distance   = alg.distance;
	popsize    = alg.popsize;
}


// createt the weight vectors with uniform distribution
void CMOEAD::init_uniformweight()
{
	if(nobj==2)
	{
		//vector<CMOEADInd> ws;
		//loadpfront("F6Weight500.dat",ws);
		//pops = 500;

		for(int n=0; n<pops; n++)
		{
			CSUB sub;
			double a = 1.0*n/(pops - 1);
			sub.namda.push_back(a);
			sub.namda.push_back(1-a);

			//load weight vectors from file
			//sub.namda.push_back(ws[n].y_obj[0]);
			//sub.namda.push_back(ws[n].y_obj[1]);

			population.push_back(sub);
		}
		popsize = pops;
	}
	else
	{
		for(int i=0; i<=unit; i++)
		{
			for(int j=0; j<=unit; j++)
			{
				if(i+j<=unit)
				{
					CSUB sub;
					sub.array.push_back(i);
					sub.array.push_back(j);
					sub.array.push_back(unit-i-j);
					for(int k=0; k<sub.array.size(); k++)
						sub.namda.push_back(1.0*sub.array[k]/unit);
					population.push_back(sub);
				}
			}
		}

		popsize = population.size();
		pops    = popsize;
	}
}

void CMOEAD::init_neighbourhood()
{
    double *x   = new double[population.size()];
	int    *idx = new int[population.size()];
	for(int i=0; i<population.size(); i++)
	{
		// calculate the distances based on weight vectors
		for(int j=0; j<population.size(); j++)
		{
		    x[j]    = dist_vector(population[i].namda,population[j].namda);
			idx[j]  = j;
		}

		// find 'niche' nearest neighboring subproblems
		minfastsort(x,idx,population.size(),niche);
		for(int k=0; k<niche; k++)
		{
			population[i].table.push_back(idx[k]);
		}

	}
    delete [] x;
	delete [] idx;
}

void CMOEAD::update_problem(CMOEADInd &indiv, int id, int type)
{
	// indiv: child solution
	// id:   the id of current subproblem
	// type: update solutions in - neighborhood (1) or whole population (otherwise)
	int size, time = 0;
	if(type==1)	size = population[id].table.size();
	else        size = population.size();
	int *perm = new int[size];
	random_permutation(perm, size);
    for(int i=0; i<size; i++)
	{
		int k;
		if(type==1) k = population[id].table[perm[i]];
		else        k = perm[i];

		// calculate the values of objective function regarding the current subproblem
		double f1, f2;
		f1 = fitnessfunction(population[k].indiv.y_obj, population[k].namda, ind_arr);
		f2 = fitnessfunction(indiv.y_obj, population[k].namda, ind_arr);
		if(f2<f1)
		{
			population[k].indiv = indiv;
			time++;
		}
		// the maximal number of solutions updated is not allowed to exceed 'limit'
		if(time>=limit)
		{
			return;
		}
	}
	delete [] perm;
}

void CMOEAD::update_reference(CMOEADInd &ind)
{
	//ind: child solution
	for(int n=0; n<nobj; n++)
	{
		if(ind.y_obj[n]<idealpoint[n])
		{
			idealpoint[n] = ind.y_obj[n];
			ind_arr[n]    = ind;
		}
	}
}

void CMOEAD::matingselection(vector<int> &list, int cid, int size, int type){
	// list : the set of the indexes of selected mating parents
	// cid  : the id of current subproblem
	// size : the number of selected mating parents
	// type : 1 - neighborhood; otherwise - whole population
	int ss   = population[cid].table.size(), r, p;
    while(list.size()<size)
	{
		if(type==1){
		    r = int(ss*rnd_uni(&rnd_uni_init));
			p = population[cid].table[r];
		}
		else
			p = int(population.size()*rnd_uni(&rnd_uni_init));

		bool flag = true;
		for(int i=0; i<list.size(); i++)
		{
			if(list[i]==p) // p is in the list
			{
				flag = false;
				break;
			}
		}

		if(flag) list.push_back(p);
	}
}

void CMOEAD::diffevolution()
{
	pops = population.size();
	int *perm = new int[pops];
	random_permutation(perm, pops);

    for(int i=0; i<pops; i++)
	{
		int n = perm[i];
		// or int n = i;
		int type;
        double rnd = rnd_uni(&rnd_uni_init);

		// mating selection based on probability
		if(rnd<realb)    type = 1;   // neighborhood
		else             type = 2;   // whole population

		// select the indexes of mating parents
		vector<int> p;
		matingselection(p,n,2,type);  // neighborhood selection

		// produce a child solution
		CMOEADInd child;
		diff_evo_xover2(population[n].indiv,population[p[0]].indiv,population[p[1]].indiv,child);

		// apply polynomial mutation
		realmutation(child, 1.0/nvar);

		// evaluate the child solution
		child.obj_eval();

		// update the reference points and other solutions in the neighborhood or the whole population
		update_reference(child);
		update_problem(child, n, type);

		p.clear(); 	nfes++;
	}

    delete [] perm;
}


vector<double> CMOEAD::execute(int run, char *strfunc, char *stralg)
{
	seed = (seed + 23)%1377;
	rnd_uni_init = -(long)seed;

    strcpy(strFunctionType,strfunc);
	strcpy(strAlgorithmType,stralg);

	vector<double> gd;
    char filename[1024];

	// load the representative Pareto-optimal solutions
    sprintf(filename,"PF/pf_%s.dat",strTestInstance);
	loadpfront(filename,ps);

	// initialization
	int gen   = 1;
	nfes      = 0;
	init_uniformweight();
    init_neighbourhood();
	init_population();
    calc_distance();

	gd.push_back(0);  gd.push_back(distance);  // id and igd value

	// evolution
	for(gen=2; gen<=max_gen; gen++)
	{

		diffevolution();

		int dd = int(max_gen/25.0);

		// calculate igd-values
        if(gen%dd==0)
		{
	   	   calc_distance();
	   	   cout<<"gen = "<<gen<<"  gd = "<<distance<<"  "<<endl;
		   gd.push_back(int(1.0*gen/dd)); gd.push_back(distance);
		}

		// save the final population - F space
        if(gen%max_gen==0)
		{
	       sprintf(filename,"POF/PF_DMOEA_%s_R%d_G%d.dat",strTestInstance,run,gen);
	       save_front(filename);
		}

		// save the final population - X space
		if(gen%max_gen==0)
		{
		   sprintf(filename,"POS/PS_DMOEA_%s_R%d_G%d.dat",strTestInstance,run,gen);
		   save_ps(filename);
		}
	}

	population.clear();
	ps.clear();
	std::cout<<" Outcome of the "<<run<<"th run:  distance= "<<distance<<" nfes = "<<nfes<<endl;
	return gd;
}


void CMOEAD::save_front(char saveFilename[1024])
{
    std::fstream fout;
	fout.open(saveFilename,std::ios::out);
	for(int n=0; n<popsize; n++)
	{
		for(int k=0;k<nobj;k++)
			fout<<population[n].indiv.y_obj[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}

void CMOEAD::save_ps(char saveFilename[1024])
{
    std::fstream fout;
	fout.open(saveFilename,std::ios::out);
	for(int n=0; n<popsize; n++)
	{
		for(int k=0;k<nvar;k++)
			fout<<population[n].indiv.x_var[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}


void CMOEAD::calc_distance()
{
    distance = 0;
	for(int i=0; i<ps.size(); i++)
	{
	    double min_d = 1.0e+10;
		for(int j=0; j<population.size(); j++)
		{
            double d = dist_vector(ps[i].y_obj, population[j].indiv.y_obj);
			if(d<min_d)  min_d = d;
		}
		distance+= min_d;
	}
	distance/=ps.size();
}

#endif
