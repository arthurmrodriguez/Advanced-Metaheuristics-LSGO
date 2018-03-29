#include <math.h>

#include "GenealogyLogStat.h"

#include "../gaconfig.h"
#include "../GAPopulation.h"
#include "../GAGenealogy.h"
#include "../GAGenealogyMemory.h"
#include "../genomes/GAGenome.h"

GenealogyLogStat::GenealogyLogStat(const Algorithm& alg) : LogStat("genealogy",alg) {
  return;
}

void GenealogyLogStat::computeValues(double &max, double &min, double &avg, double &dev, double &imperc, double &immax, double &immin, double &imavg, double &imdev, double &fac){
  const GAPopulation& pop = alg_.population();

  max = min = avg = dev = imperc = immax = immin = imavg = imdev = fac = 0.0;

  //if(geneal)
  if(GAGenealogy::handle()->isGenealogyMemory())
      {
	GAGenealogyMemory *genealmem = DYN_CAST(GAGenealogyMemory*, GAGenealogy::handle());

	GAGenome gen;
	int i, nimm, firstg, actgen = genealmem->getGeneration();
	double tmpfitage, tmpage, tmpdev, tmpimm = 0.0, tmpmax = 0.0, tmpmin = 30000;
	double tmpvar, tmpimax = 0.0, tmpimin = 30000;
	int island = genealmem->getIsland();

	//Get the immigrant percentage in the actual population [0-1]
	//Age Average of the population and the immigrant population
	for(tmpvar=tmpimm=0.0, nimm=i=0; i<(int)pop.size(); i++)
	  {
	    gen = pop.individual(i);
	    GAGenomeNode *genNode = genealmem->getGeneNode(gen);
	    firstg = genNode->getFirstG();
	    if(genNode->getFirstG() < 0)
	      tmpage = actgen + firstg;
	    else
	      tmpage = actgen - firstg;
	    tmpvar += tmpage;

	    if(tmpage > tmpmax)
	      tmpmax = tmpage;
	    if(tmpage < tmpmin)
	      tmpmin = tmpage;

	    if(gen.getIsland() != island)
	      {
		tmpimm += tmpage;
		nimm++;
		if(tmpage > tmpimax)
		  tmpimax = tmpage;
		if(tmpage < tmpimin)
		  tmpimin = tmpage;
	      }
	  }

	avg = tmpvar/(double) pop.size();
	max = tmpmax;
	min = tmpmin;
	if(nimm)
	  {
	    imperc = (double) nimm/(double) pop.size();
	    imavg = tmpimm/(double) nimm;
	    immax = tmpimax;
	    immin = tmpimin;
	  }
	else
	  imperc = imavg = immax = immin = 0.0;

	//Get the deviation age of the population and the immigrant population.
	//Also get the fitness-age correlation
	for(tmpfitage=tmpvar=tmpdev=0.0, i=0; i<(int)pop.size(); i++)
	  {
	    GAGenomeNode *genNode = genealmem->getGeneNode(pop.individual(i));
	    firstg = genNode->getFirstG();
	    if(firstg < 0)
	      tmpage = actgen + firstg;
	    else
	      tmpage = actgen - firstg;

	    tmpvar += pow((double) (tmpage) - avg, 2);
	    if(nimm)
	      if(pop.individual(i).getIsland() != island)
		tmpdev += pow((double) (tmpage) - imavg, 2);

	    tmpfitage += (pop.individual(i).fitness() - pop.fitave())*(tmpage - avg);
	  }

	dev = sqrt(tmpvar/((double) pop.size()-1));
	if(nimm && tmpdev != 0.0)
	  imdev = sqrt(tmpdev/((double) nimm-1));
	else
	  imdev = 0.0;

	if(pop.fitdev() != 0.0 && dev != 0.0)
	  fac = (tmpfitage/(double) pop.size())/(pop.fitdev()*dev);
	else
	  fac = 0.0;
      }
}


void GenealogyLogStat::update(){
  double max,min,avg,dev,ipr,ima,imi,iav,ide,fac;

  computeValues(max,min,avg,dev,ipr,ima,imi,iav,ide,fac);

  sprintf(message_,"%s=[max=%10.10lf,min=%10.10lf,avg=%10.10lf,dev=%10.10lf,ipr=%10.10lf,ima=%10.10lf,imi=%10.10lf,iav=%10.10lf,ide=%10.10lf,fac=%10.10lf] ",name_.c_str(),max,min,avg,dev,ipr,ima,imi,iav,ide,fac);
}
