#include "NativePrcntLogStat.h"
#include <math.h>
#include <assert.h>

NativePrcntLogStat::NativePrcntLogStat(int rank) : SingleLogStat("native_prcnt"), rank_(rank)  {}

NativePrcntLogStat::~NativePrcntLogStat(){}

void NativePrcntLogStat::update(GAPopulation& pop){
  long double value = 0.0;
  
  unsigned pop_size = pop.size();
        
  for (unsigned i=0; i<pop_size; i++) value  += (pop.individual(0).origin() == rank_ ) ? 1 : 0;
    
  value  /= pop_size;
     
  assert(!isnan(value)); 
  
  setMessage(value);
}