#include "GenDivLogStat.h"
#include <math.h>
#include <assert.h>

GenDivLogStat::GenDivLogStat(const Algorithm& alg) : CollectionLogStat("gendiv",alg){}
GenDivLogStat::~GenDivLogStat(){}
  
void GenDivLogStat::computeValues(double& max, double& min, double& avg, double& dev){
  const GAPopulation& pop = alg_.population();
  pop.diversity();                            // Computes (if not already) the distance between the individuals
  
  max = alg_.statistics().gendistCurrent(GAStatistics::Maximum);
  min = alg_.statistics().gendistCurrent(GAStatistics::Minimum);
  avg = alg_.statistics().gendistCurrent(GAStatistics::Mean);
  dev = alg_.statistics().gendistCurrent(GAStatistics::Deviation);
}
   
double GenDivLogStat::indValue(GAGenome& ind){return 0.0;}
