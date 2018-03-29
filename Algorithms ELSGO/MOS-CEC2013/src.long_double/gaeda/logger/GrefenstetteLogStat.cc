#include "GrefenstetteLogStat.h"
#include "genomes/GA1DBinStrGenome.h"
#include <math.h>
#include <assert.h>
#include <stdexcept>

GrefenstetteLogStat::GrefenstetteLogStat() : SingleLogStat("gref_bias"){}

GrefenstetteLogStat::~GrefenstetteLogStat(){}

void GrefenstetteLogStat::update(GAPopulation& pop){
  
  GA1DBinaryStringGenome* gen_tmp = dynamic_cast<GA1DBinaryStringGenome*>( &( pop.best() ) );
   
  if (!gen_tmp) throw runtime_error("In GrefenstetteLogStat, genomes are not binary");
   
  unsigned pop_size = pop.size(); 
   
  long double sum_genes, sum_compl_genes, gref_bias = 0.0;
  
  int num_genes = gen_tmp->length();
  for (int gen_pos=0; gen_pos<num_genes; gen_pos++){
    sum_genes = sum_compl_genes = 0;
    for (unsigned ind_pos=0; ind_pos<pop_size; ind_pos++){
      gen_tmp = dynamic_cast<GA1DBinaryStringGenome*>( &( pop.individual(ind_pos) ) );
      
      sum_genes       += gen_tmp->gene(gen_pos);
      sum_compl_genes += 1-gen_tmp->gene(gen_pos);
    }    
    gref_bias += (sum_genes > sum_compl_genes) ? sum_genes : sum_compl_genes;
  }
  gref_bias /= (num_genes * pop_size);
  
  gref_bias = 2*(1-gref_bias); // Corrected version so its range is in (0,1)
  
  assert(!isnan(gref_bias));
  assert(gref_bias <= 1.0 && gref_bias >= 0.0);
       
  setMessage(gref_bias);
}

