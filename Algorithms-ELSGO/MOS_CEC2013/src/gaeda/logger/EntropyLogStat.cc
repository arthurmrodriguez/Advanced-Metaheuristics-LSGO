#include "EntropyLogStat.h"
#include "../genomes/GA1DBinStrGenome.h"
#include <math.h>
#include <assert.h>
#include <stdexcept>

EntropyLogStat::EntropyLogStat() : SingleLogStat("entropy"){}

EntropyLogStat::~EntropyLogStat(){}

void EntropyLogStat::update(GAPopulation& pop){
  
  GA1DBinaryStringGenome* gen_tmp = dynamic_cast<GA1DBinaryStringGenome*>( &( pop.best() ) );
   
  if (!gen_tmp) throw runtime_error("In EntropyLogStat, genomes are not binary");
   
  unsigned pop_size = pop.size(); 
   
  double tmp_num_ones, tmp_num_zeros, entropy = 0.0;

  int num_genes = gen_tmp->length();
  for (int gen_pos=0; gen_pos<num_genes; gen_pos++){
    tmp_num_ones = 0;
    for (unsigned ind_pos=0; ind_pos<pop_size; ind_pos++){
      gen_tmp = dynamic_cast<GA1DBinaryStringGenome*>( &( pop.individual(ind_pos) ) );
      if (gen_tmp->gene(gen_pos) == 1 ) tmp_num_ones++;
    }
    tmp_num_zeros = pop_size - tmp_num_ones;
    
    // with no ones or zeros the logarithm returns nan so we dont consider this case
    if (tmp_num_zeros != 0 && tmp_num_ones != 0) {        
      entropy += -tmp_num_ones/pop_size * log2(tmp_num_ones/pop_size) - tmp_num_zeros/pop_size * log2(tmp_num_zeros/pop_size);
    }
  }
  entropy /= num_genes;
    
  assert(!isnan(entropy));
       
  setMessage(entropy);
}

