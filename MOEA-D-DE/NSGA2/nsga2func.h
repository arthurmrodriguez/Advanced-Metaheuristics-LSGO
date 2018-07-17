#ifndef __CNSGA2CLASS_H_
#define __CNSGA2CLASS_H_

#include "../common/global.h"
#include "../common/recombination.h"
#include "../common/common.h"
#include "nsga2class.h"

class CNSGA2
{
public:
	CNSGA2();
	virtual ~CNSGA2();

	vector<double> execute(int run);
	void init_population();
	void calc_distance();
	int  tour_selection();
	void fill_union(CNSGA2Ind &ind);
	void rank_popu();
	void eval_dens();
	void diffevolution();
	void save_front(char savefilename[1024]);                              // save the pareto front into files
	void save_ps(char savefilename[1024]);

    vector <CNSGA2Ind>  population;
	vector <CNSGA2Ind>  ps;
	vector <CNSGA2Ind>  offspring;
	double distance;

	int    nfes;

    void operator=(const CNSGA2 &alg);
};

CNSGA2::CNSGA2(){

}

CNSGA2::~CNSGA2(){

}

void CNSGA2::operator=(const CNSGA2 &alg)
{
	population = alg.population;
	offspring  = alg.offspring;
	ps         = alg.ps;
	nfes       = alg.nfes;
	distance   = alg.distance;
}

void CNSGA2::init_population()
{

	for(int n=0;n<pops;n++)
	{
		CNSGA2Ind ind;
		ind.rnd_init();
		population.push_back(ind);
		nfes++;
	}
}

int  CNSGA2::tour_selection()
{
	int p1 = int(rnd_uni(&rnd_uni_init)*pops);
	int p2 = int(rnd_uni(&rnd_uni_init)*pops);

	if(population[p1].rank<population[p2].rank)
		return p1;
	else
		return p2;
}


void CNSGA2::fill_union(CNSGA2Ind &ind)
{
    bool flag = true;
	int  size = offspring.size();
	for(int i=0; i<size; i++){
		if(ind==offspring[i]){
            flag = false;
			break;
		}
	}
	if(flag) offspring.push_back(ind);
}

void CNSGA2::diffevolution()
{
    for(int n=0; n<pops; n++)
	{
		int p1, p2, p3;
		p1 = tour_selection();
		while(1){  p2 = tour_selection();  	if(p2!=p1) break; }
		while(1){  p3 = tour_selection();  	if(p3!=p1&&p3!=p2) break; }
		CNSGA2Ind child;
		diff_evo_xover2(population[p1],population[p2],population[p3],child);
		realmutation(child, 1.0/nvar);
		child.obj_eval();

		// save both parents and offspring into the combined population
        fill_union(child);
		fill_union(population[n]);

		nfes++;
	}
}

void CNSGA2::rank_popu()
{
	int size = offspring.size();
	int** cset;

	int* rank = new int[size];

	cset = new int*[size];

    for(int i=0; i<size; i++)
        cset[i] = new int[size];

    for(int i=0; i<size; i++)
	{
		rank[i]  = 0;
		offspring[i].rank  = -1;
		offspring[i].count = 0;

	}

	for(int k=0; k<size; k++)
	    for(int j=0; j<size; j++)
		{
			if(k!=j)
			{
				if(offspring[j]<offspring[k]&&!(offspring[j]==offspring[k])) rank[k]++;

				if(offspring[k]<offspring[j]&&!(offspring[j]==offspring[k]))
				{
					offspring[k].count++;
					int m = offspring[k].count - 1;
					cset[k][m] = j;
				}
			}
		}

	int curr_rank = 0;
	while(1)
	{
		int stop_count = 0;
		int* rank2 = new int[size];
	    for(int k=0; k<size; k++)
			rank2[k] = rank[k];

        for(int k=0; k<size; k++)
		{
		    if(offspring[k].rank==-1&&rank[k]==0)
			{
			    offspring[k].rank = curr_rank;
				for(int j=0; j<offspring[k].count; j++)
				{
				   int id =	cset[k][j];
				   rank2[id]--;
				   stop_count++;
				}
			}
		}

	    for(int k=0; k<size; k++)
			rank[k] = rank2[k];

        delete [] rank2;
		curr_rank++;
		if(stop_count==0) break;
	}

    delete [] rank;

    for(int i=0; i<size; i++)
        delete cset[i];
	delete[] cset;
}

void CNSGA2::eval_dens()
{
	population.clear();
	int size = offspring.size();
	int rank = 0;
	while(1){
		int count = 0;
        for(int i=0; i<size; i++)
			if(offspring[i].rank==rank)
				count++;

		int size2 = population.size() + count;
		if(size2>pops) {
			break;
		}

        for(int i=0; i<size; i++)
  	        if(offspring[i].rank==rank)
			    population.push_back(offspring[i]);
		rank++;
		if(population.size()>=pops) break;
	}

	if(population.size()<pops){
	    vector<CNSGA2Ind> list;
		// save the individuals in the overflowed front
        for(int i=0; i<size; i++)
  	        if(offspring[i].rank==rank)
		        list.push_back(offspring[i]);
		int s2 = list.size();
		double *density = new double[s2];
		int    *idx     = new int[s2];

		for(int i=0; i<s2; i++){
			idx[i]     = i;
			density[i] = 0;
		}

		int    *idd     = new int[s2];
		double *obj     = new double[s2];
		for(int j=0; j<nobj; j++){
			for(int i=0; i<s2; i++){
			    idd[i] = i;
				obj[i] = list[i].y_obj[j];
			}
			minfastsort(obj,idd,s2,s2);
            density[idd[0]]    += -1.0e+30;
            density[idd[s2-1]] += -1.0e+30;
			for(int k=1; k<s2-1; k++)
				density[idd[k]]+= -(obj[k] - obj[k-1] + obj[k+1] - obj[k]);
		}
		delete [] idd;
		delete [] obj;

		int s3 = pops - population.size();

		minfastsort(density,idx,s2,s3);

		for(int i=0; i<s3; i++)
			population.push_back(list[idx[i]]);

		delete [] density;
		delete [] idx;
	}
	offspring.clear();
}

vector<double> CNSGA2::execute(int run)
{

	seed = (seed + 23)%1377;
	rnd_uni_init = -(long)seed;
	vector<double> gd;

	char filename[1024];

    sprintf(filename,"PF/pf_%s.dat",strTestInstance);
	loadpfront(filename,ps);
	nfes    = 0;


	// first generation
	int gen   = 1;
	init_population();
    calc_distance();
	gd.push_back(0);  gd.push_back(distance);


	for(int gen=2; gen<=max_gen; gen++)
	{
	    diffevolution();
		rank_popu();
		eval_dens();
		int dd = int(max_gen/25.0);
        if(gen%dd==0)
		{
	   	   calc_distance();
	   	   cout<<"gen = "<<gen<<"  gd = "<<distance<<"  "<<endl;
		   gd.push_back(int(1.0*gen/dd)); gd.push_back(distance);
		}

        // save the final population - F space
        if(gen%max_gen==0)
		{
	       sprintf(filename,"POF/PF_NSGA2_%s_R%d_G%d.dat",strTestInstance,run,gen);
	       save_front(filename);
		}

        // save the final population - X space
        if(gen%max_gen==0)
		{
		   sprintf(filename,"POS/PS_NSGA2_%s_R%d_G%d.dat",strTestInstance,run,gen);
		   save_ps(filename);
		}
	}

	population.clear();
	offspring.clear();
	ps.clear();
	std::cout<<" Outcome of the "<<run<<"th run:  distance= "<<distance<<" nfes = "<<nfes<<endl;
	return gd;

}

void CNSGA2::save_front(char saveFilename[1024])
{
    std::fstream fout;
	fout.open(saveFilename,std::ios::out);
	for(int n=0; n<pops; n++)
	{
		for(int k=0;k<nobj;k++)
			fout<<population[n].y_obj[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}

void CNSGA2::save_ps(char saveFilename[1024])
{
    std::fstream fout;
	fout.open(saveFilename,std::ios::out);
	for(int n=0; n<pops; n++)
	{
		for(int k=0;k<nvar;k++)
			fout<<population[n].x_var[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}


void CNSGA2::calc_distance()
{
    distance = 0;
	for(int i=0; i<ps.size(); i++)
	{
	    double min_d = 1.0e+10;
		for(int j=0; j<population.size(); j++)
		{
            double d = dist_vector(ps[i].y_obj, population[j].y_obj);
			if(d<min_d)  min_d = d;
		}
		distance+= min_d;
	}
	distance/=ps.size();
}

#endif
