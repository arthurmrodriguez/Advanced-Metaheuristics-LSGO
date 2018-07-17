#include "CollectionLogStat.h"

#include <math.h>
#include <assert.h>

#include "../Algorithm.h"
#include "../GAPopulation.h"
#include "../genomes/GAGenome.h"

CollectionLogStat::CollectionLogStat(string name, const Algorithm& alg): LogStat(name,alg) {}

CollectionLogStat::~CollectionLogStat(){}

void CollectionLogStat::update(){
  long double max,min,avg,dev;

  computeValues(max,min,avg,dev);

  sprintf(message_,"%s=[min=%.16LE,max=%.16LE,avg=%.16LE,dev=%.16LE] ",name_.c_str(),min,max,avg,dev );
}

void CollectionLogStat::computeValues(long double& max, long double& min, long double& avg, long double& dev){
  const GAPopulation& pop = alg_.population();
  long double              var, ind_value;

  long double sum = 0.0, sqr_sum = 0.0;

  unsigned pop_size = pop.size();

  if (pop_size == 1) {
    max = min = avg = pop.individual(0).score();
    dev = var = 0.0;
    return;
  }

  ind_value = indValue(pop.individual(0));

  max = min = sum = ind_value;
  sqr_sum = pow(ind_value,2.0);

  for (unsigned i=1; i<pop_size; i++) {
    ind_value = indValue(pop.individual(i));
    sum     += ind_value;
    sqr_sum += pow(ind_value ,2.0);
    max     =  GAMax( max,ind_value );
    min     =  GAMin( min,ind_value );
  }

  avg = sum / pop_size;
  var = (sqr_sum /( pop_size - 1.0) ) - ( ( (long double)pop_size/(pop_size - 1.0) ) * pow(avg,2.0) );

  if (var < 0.0 && var > -0.0000001 ) var = 0.0;

//  assert(!isnan(var)); assert(var >= 0 );

  dev = sqrt(var);
}
