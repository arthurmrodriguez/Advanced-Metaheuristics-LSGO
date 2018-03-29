#include "NSC.h"

#include "GAPopulation.h"
#include "GAGenealogy.h"
#include "GAGenealogyMemory.h"
#include "genomes/MOSGenome.h"

NSC::NSC (unsigned size) : _sz (size), _nscn (size) {

  _nscd = (nscdata*) malloc (sizeof (nscdata) * _sz);

}

NSC::~NSC () {

  if (_nscd)
    free (_nscd);

}

void NSC::addPopulation (const GAPopulation* pop) {

  for (unsigned i=0; i < _sz; i++) {
    _nscd[i].id = pop->individual(i).getId();
    _nscd[i].island = pop->individual(i).getIsland();
    _nscd[i].fit = pop->individual(i).score();
  }

}

void NSC::prepareNSCData (const GAPopulation* pop) {

  GAGenealogyMemory* genealmem = dynamic_cast<GAGenealogyMemory*> (GAGenealogy::handle());

  unsigned long int idgen, lastid;
  int idis, lastis;

  GAGenome& g = pop->individual(0);
  lastid = _nscd[0].id = g.getId();
  lastis = _nscd[0].island = g.getIsland();
  _nscd[0].fit = g.score();
  _nscn = 1;

  // Copy individual and island ids (only from genomes that doesn't exist yet)
  for (unsigned i=1; i<pop->size(); i++) {
    GAGenome& gen = pop->individual(i);
    idgen = gen.getId();
    idis = gen.getIsland();

    if(lastid != idgen  || (lastid == idgen && lastis != idis)) {
       GAGenomeNode *genNode = genealmem->getGeneNode(gen);
       lastid = _nscd[_nscn].id = idgen;
       lastis = _nscd[_nscn].island = idis;
       _nscd[_nscn].fit = genNode->getScore();
       _nscn++;
    }
  }

}

long double NSC::computeNSC (const GAPopulation* pop, int tech) {

  GAGenome gen;
  unsigned i, j, k;
  unsigned part, npart, rest, ntotal;
  long double fitaux = 0.0, *dadfitlist = NULL, *childfitlist = NULL, *fitChild = NULL;
  long double negres = 0.0;   // Original NSC is with negative cummulated values
  long double totalres = 0.0; // Result with all cummulated values

  const GAGenealogyMemory *genealmem = dynamic_cast<GAGenealogyMemory*> (GAGenealogy::handle());

  // Reset children values
  fitChild = (long double*) malloc(_nscn*sizeof(long double));
  for (i=0; i<_nscn; i++)
    fitChild[i] = 0.0;

  // Put children values in the list
  for (i=0; i<pop->size(); i++) {
    MOSGenome& gen = dynamic_cast<MOSGenome&> (pop->individual(i));

    if (gen.getTechniqueId() != tech)
       continue;

    GAGenomeNode *genNode = genealmem->getGeneNode(gen);

    // If it was created in current generation and it's from this island
    if (genNode->getFirstG() == genealmem->getGeneration() &&
        genNode->getIsland() == genealmem->getIsland()          ) {

       // Update dad's value
       if (genNode->getIdDad())
         for (j=0; j<_nscn; j++)
           if (genNode->getIdDad() == _nscd[j].id && genNode->getIslandDad() == _nscd[j].island)
             if (fitChild[j] < genNode->getScore()) {
                fitChild[j] = genNode->getScore();
                break;
             }

       // Update mom's value
       if (genNode->getIdMom())
         for (j=0; j<_nscn; j++)
           if (genNode->getIdMom() == _nscd[j].id && genNode->getIslandMom() == _nscd[j].island)
             if (fitChild[j] < genNode->getScore()) {
                fitChild[j] = genNode->getScore();
                break;
             }

     }

  }

  // Number of gens in a group
  if(_nscn < 50)
    part = 2;
  else
    part = 10;

  npart = _nscn/part;
  rest  = _nscn - npart*part;

  // Number of partitions
  if(rest)
    ntotal=1+npart;
  else
    ntotal=npart;

  // Create the auxiliar arrays
  dadfitlist   = (long double*) malloc(ntotal*sizeof(long double));
  childfitlist = (long double*) malloc(ntotal*sizeof(long double));

  // Compute the mean in the auxiliar lists
  dadfitlist  [0] = 0.0;
  childfitlist[0] = 0.0;

  for (i=0, j=0, k=0; i<_nscn; i++, j++) {

    if (j==part) {

       dadfitlist  [k] = dadfitlist  [k]/(long double)part;
       childfitlist[k] = childfitlist[k]/(long double)part;

       k++;
       j=0;

       dadfitlist  [k] = 0.0;
       childfitlist[k] = 0.0;

    }

    // If a genome didn't have children, then put his fitness
    if (fitChild[i] == 0.0)
       fitChild[i] = _nscd[i].fit;

    dadfitlist  [k] += _nscd[i].fit;
    childfitlist[k] += fitChild[i];

  }

  dadfitlist  [k] = dadfitlist  [k]/(long double)j;
  childfitlist[k] = childfitlist[k]/(long double)j;

  // Compute the results
  for (i=0; i<ntotal-1; i++) {

    if (dadfitlist[i] - dadfitlist[i+1] == 0.0)
       fitaux = 0.0;
    else if (GAGenome::optCriterion() == GAGenome::MINIMIZATION)
       fitaux = (childfitlist[i] - childfitlist[i+1])/(dadfitlist[i] - dadfitlist[i+1]);
    else
       fitaux = (childfitlist[i+1] - childfitlist[i])/(dadfitlist[i+1] - dadfitlist[i]);

    if(fitaux < 0.0)
       negres += fitaux;

    totalres += fitaux;

  }

  free(fitChild);
  free(dadfitlist);
  free(childfitlist);

  negres = -1/(negres-1);

  // Original NSC show only negative difference
  return negres;

}
