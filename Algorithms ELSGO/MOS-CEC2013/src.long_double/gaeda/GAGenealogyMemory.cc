/*****************************************************************************
* GAEDAlib: A C++ GA library with EDA and multiprocessor (MPI) support      *
 *                                                                           *
* (C) 2007 Carlos Pascual (carlospmdn@gmail.com)                            *
 *                                                                           *
 * GAEDAlib is distributed under the terms of the BSD software license       *
 *                                                                           *
* GAEDAlib is heavily based on GAlib, a C++ GA library by Mathew Wall:      *
  * Copyright (c) 1995-1996 Massachusetts Institute of Technology (MIT)       *
  * Copyright (c) 1996-2000 Matthew Wall (author of GAlib)                    *
 *                                                                           *
 * Some portions of GAEDAlib's source code come from the GNU C++ compiler    *
 * library and therefore are covered under the terms of a different license, *
 * the GNU Public License.                                                   *
 *                                                                           *
 * You should have received a file named LICENSE along with this software.   *
 * This file contains more information about the licensing conditions of     *
 * GAEDAlib as well as the full text of each license involved.               *
 *                                                                           *
 * The file AUTHORS lists the people who have contributed (directly or       *
 * indirectly) to GAEDAlib                                                   *
 *****************************************************************************/

#include "GAGenealogyMemory.h"

#include "GAPopulation.h"
#include "logger/GAFileLogger.h"
#include "genomes/MOSGenome.h"

GAGenealogy* GAGenealogyMemory::create (int idis, executionType extype) {

  assert (instance_ptr == NULL);

  return instance_ptr = new GAGenealogyMemory (idis, extype);

}


GAGenealogyMemory::GAGenealogyMemory(int idis, executionType extype) : GAGenealogy(idis) {
  exType = extype;
  funcHeritage = NULL;
}


GAGenealogyMemory::~GAGenealogyMemory() {
  multimap<unsigned long int, GAGenomeNode*>::iterator it;

  //Delete all nodes in the genealogy
  for(it=genealogy.begin(); it != genealogy.end(); it++) {
    it->second->deleteGenomes();
    delete it->second;
  }

  genealogy.clear();
}


void GAGenealogyMemory::addNode(GAGenome &gen) {

  multimap<unsigned long int, GAGenomeNode*>::iterator it, itnext;
  GAGenomeNode *node;
  unsigned long int idgen = gen.getId();
  int idis = gen.getIsland();

  if (idgen == 0) {
    gen.setId(newId());
    gen.setIsland(island);
    node = new GAGenomeNode(gen, NULL, NULL, actGen, actGen, exType, false);
    genealogy.insert(pair<unsigned long int, GAGenomeNode*> (gen.getId(), node));
    return;
  }

  it = genealogy.find(idgen);

  if (it != genealogy.end()) {
    itnext = genealogy.upper_bound(idgen);

    while(it != itnext) {
      if(it->second->getIsland() == idis)
	return;
      it++;
    }
  }

  node = new GAGenomeNode(gen, NULL, NULL, actGen, actGen, exType, false);
  genealogy.insert(pair<unsigned long int, GAGenomeNode*> (idgen, node));
  return;

}


void GAGenealogyMemory::familyRelationship(GAGenome &dad, GAGenome &mom, GAGenome &child1, GAGenome &child2) {

  GAGenomeNode *node;
  unsigned long int idchild1 = 0, idchild2 = 0, idmom = 0, iddad = 0;

  iddad = dad.getId();
  idmom = mom.getId();

  child1.setId(idchild1 = newId());
  child1.setIsland(island);

  node = new GAGenomeNode(child1, &dad, &mom, actGen, actGen, exType, false);
  genealogy.insert(pair<unsigned long int, GAGenomeNode*> (idchild1, node));

  // In the case there isn't a selfcross
  if(&child1 != &child2) {
    child2.setId(idchild2 = newId());
    child2.setIsland(island);
    node = new GAGenomeNode(child2, &dad, &mom, actGen, actGen, exType, false);
    genealogy.insert(pair<unsigned long int, GAGenomeNode*> (idchild2, node));
  }

}


void GAGenealogyMemory::mutated(GAGenome &dad, GAGenome &child, bool crossAndMut) {

  GAGenomeNode *node;
  unsigned long int idchild;

  // If it was only mutated set ids and add it to the genealogy
  if(!crossAndMut) {
    child.setId(idchild = newId());
    child.setIsland(island);

    node = new GAGenomeNode(child, &dad, NULL, actGen, actGen, exType, true);
    genealogy.insert(pair<unsigned long int, GAGenomeNode*> (idchild, node));
  }
  else { // If genealogy is of type ALL or FULL set the mutated genome
    if(exType == GAGenealogyMemory::ALL || exType == GAGenealogyMemory::FULL) {
      node = getGeneNode(child); // Previous genome in genealogy has the same ids than child
      node->setMutatedGenome(child, exType);
    }
  }

}


void GAGenealogyMemory::deceased(GAGenome &gen) {

  multimap<unsigned long int, GAGenomeNode*>::iterator it, itnext;
  unsigned long int idgen = 0;
  int idis = gen.getIsland();
  stringstream mes;

  idgen = gen.getId();
  it = genealogy.find(idgen);

  if(it != genealogy.end()) {
    itnext = genealogy.upper_bound(idgen);
    while(it != itnext) {
      if(it->second->getIsland() == idis) {
	// Don't put the last generation because it's the last step of the algorithm, when merging all islands.
	if(it->second->getFirstG() > actGen)
	  return;

	it->second->setLastG(actGen);
	return;
      }
      it++;
    }
  }

  if(actGen != 0) {
    mes << "Decease error: couldn't find the genome with ids: " << idgen << " " << idis;
    GALogger::instance()->appendLogMessage("", mes.str(), GALogger::only_stats);
  }

}


//Add only the immigrants of this population to the genealogy
void GAGenealogyMemory::addImmigrant(GAGenome &gen) {

  int diffAge = 0;

  if(gen.getIsland() != island) {

    addNode(gen);
    GAGenomeNode *geneNode = getGeneNode(gen);

    // If this value is negative it means that the genome was previously in the
    // genealogy and that the algorithm is in the last step, when all populations
    // are merged. For this reason, we must not put a first generation, because it
    // will be greater than the last generation.
    if((actGen - gen.getAge()) > geneNode->getLastG())
      return;

    geneNode->setFirstG(actGen - gen.getAge());

    // Before immigration the rankings in that island don't exists
    diffAge = gen.getAge() - geneNode->getRankingList().size();

    if(diffAge > 0)
      for(int i = 0; i < diffAge; i++)
	geneNode->addRanking(-1);

  }

}


bool GAGenealogyMemory::exist(GAGenome &gen) {

  multimap<unsigned long int, GAGenomeNode*>::iterator it, itnext;
  unsigned long int idgen = gen.getId();
  int idis = gen.getIsland();

  // The genome was created without setting the generation
  if(!idgen || idis == -1)
    return false;

  it = genealogy.find(idgen);

  if(it != genealogy.end()) {

    itnext = genealogy.upper_bound(idgen);

    while(it != itnext) {
      if(it->second->getIsland() == idis)
	return true;
      it++;
    }

  }

  return false;

}


// Put the ranking of each genome  in the current population in the genealogy
void GAGenealogyMemory::updateRankings(const GAPopulation *pop) {

  if (GAGenome::optCriterion() == GAGenome::MINIMIZATION) {
    for (int i = 0; i < (int)pop->size(); i++) {
      int j = i-1;
      unsigned long int id = pop->individual(i).getId();
      int is = pop->individual(i).getIsland();

      //If it is a clone of a previous genome, don't insert the ranking
      for (; j>=0; j--)
	if (pop->individual(j).getId() == id && pop->individual(j).getIsland() == is)
	  break;

      if (j == -1)
	addGenomeRanking(id, is, i);
    }
  }
  else {
    for (int i = (int)pop->size()-1; i>=0; i--) {
      int j = i+1;
      unsigned long int id = pop->individual(i).getId();
      int is = pop->individual(i).getIsland();

      //If it is a clone of a previous genome don't insert the ranking
      for (; j<pop->size(); j++)
	if (pop->individual(j).getId() == id && pop->individual(j).getIsland() == is)
	  break;

      if (j == pop->size())
	addGenomeRanking(id, is, i);
    }
  }

}


void GAGenealogyMemory::addGenomeRanking(GAGenome &gen, unsigned int ranking) const {
  return addGenomeRanking(gen.getId(), gen.getIsland(), ranking);
}


void GAGenealogyMemory::addGenomeRanking(unsigned long int idGen, int idIs, unsigned int ranking) const {
  GAGenomeNode *node;
  node = getGeneNode(idGen, idIs);
  node->addRanking((signed short) ranking);
}


void GAGenealogyMemory::heritage(GAGenome &gen1, GAGenome &gen2) {
  return (*funcHeritage) (gen1, gen2);
}


GAGenomeNode *GAGenealogyMemory::getGeneNode(GAGenome &gen) const {
  return getGeneNode(gen.getId(), gen.getIsland());
}


GAGenomeNode *GAGenealogyMemory::getGeneNode(unsigned long int idgen, int idis) const {
  multimap<unsigned long int, GAGenomeNode*>::const_iterator it, itnext;

  stringstream mes;

  // The genome was created without setting the generation or the island.
  if(!idgen || idis == -1) {
    mes << "Get Node error: couldn't find the genome with ids: " << idgen << " " << idis;
    GALogger::instance()->appendLogMessage("", mes.str(), GALogger::only_stats);
    abort();
  }

  it = genealogy.find(idgen);

  if(it != genealogy.end()) {
    itnext = genealogy.upper_bound(idgen);

    while(it != itnext) {
      if(it->second->getIsland() == idis)
	return it->second;
      it++;
    }

  }

  mes << "Get Node error: couldn't find the genome with ids: " << idgen << " " << idis;
  GALogger::instance()->appendLogMessage("", mes.str(), GALogger::only_stats);
  abort();

}


// Structure to store technique and score for best individuals
typedef struct {
  int tech;
  long double fit;
} bestGen;


void GAGenealogyMemory::getLastBests(int *lastBests, int *techList, int nAlgs, int lastGen, int bestN) const {

  multimap<unsigned long int, GAGenomeNode*>::const_iterator it, itnext;
  GAGenomeNode *node;
  unsigned firstGen = actGen - lastGen, genAux = actGen-1;
  bestGen *bestList = (bestGen*) malloc(bestN*sizeof(bestGen));
  int i, j, bestListN = 0;
  long double fitAux;

  // Initialize to 0
  for (i=0; i<nAlgs; i++)
    lastBests[i] = 0;

  // If we have one generation to search the best genomes
  if (lastGen) {

    it = genealogy.end();
    it--;

    // Search within our generational period
    while (it->second->getFirstG() >= firstGen) {

      node = it->second;

      // Only with our genomes
      if (it->second->getIsland() == island) {

	// In a new generationwe  must get the number of best genomes per technique and reset the auxiliar vector
	if (genAux != node->getFirstG()) {

	  for (i=0; i<bestListN; i++)
	    for (j=0; j<nAlgs; j++)
	      if(techList[j] == bestList[i].tech) {
		lastBests[j]++;
		break;
	      }

	  bestListN = 1;
	  bestList[0].fit = node->getScore();
	  bestList[0].tech = node->getTechnique();
	  genAux = node->getFirstG();

	}
	else {

	  // We only have to add the technique and the score
	  if(bestListN < bestN) {
	    bestList[bestListN].fit = node->getScore();
	    bestList[bestListN].tech = node->getTechnique();
	    bestListN++;
	  }
	  else { // We must search the worst and exchange ir for the other one

	    fitAux = bestList[0].fit;
	    j = 0;

	    for(i=1; i<bestListN; i++)
	      if(GAGenome::compareScores(bestList[i].fit, fitAux) == GAGenome::WORSE) {
		fitAux = bestList[i].fit;
		j = i;
	      }

	    if(GAGenome::compareScores(node->getScore(), fitAux) == GAGenome::BETTER) {
	      bestList[j].fit = node->getScore();
	      bestList[j].tech = node->getTechnique();
	    }

	  }

	}

      }

      // If it is the first node of the multihash table
      if (it == genealogy.begin())
	break;
      else
	it--;

    }

  }

  // Add the number of genomes to the result
  for (i=0; i<bestListN; i++)
    for (j=0; j<nAlgs; j++)
      if (techList[j] == bestList[i].tech) {
	lastBests[j]++;
	break;
      }

  delete bestList;

  return;

}


void GAGenealogyMemory::print() {

  multimap<unsigned long int, GAGenomeNode*>::iterator it;
  list<signed short> rankList;
  list<signed short>::iterator rankIt;
  GAGenome *mutGen, *crossGen;

  stringstream mes;

  mes << "\n#GenealogyStats with " << actGen << " generations\n";
  mes << "#Type of execution: ";

  if(exType == GAGenealogyMemory::ALL)
    mes << "All\n";
  else if(exType == GAGenealogyMemory::FULL)
    mes << "Full\n";
  else
    mes << "Normal\n";

  mes << "# Id Island FirstGeneration LastGeneration  ";

  if(exType == GAGenealogyMemory::FULL)
    mes << "ScoreCross  ScoreMutation";
  else
    mes << "Score";

  mes << "  Technique   DadId DadIsland   MomId MomIsland";
  GALogger::instance()->appendLogMessage("", mes.str(), GALogger::only_stats);

  for (it=genealogy.begin(); it!=genealogy.end(); it++) {

    stringstream mesaux;

    crossGen = it->second->getCrossoverGenome();
    mutGen   = it->second->getMutationGenome ();

    mesaux << it->first << " " << it->second->getIsland() << " " << it->second->getFirstG() << " " << it->second->getLastG() << " ";

    if(exType == GAGenealogyMemory::FULL) {

      if(crossGen)
	mesaux << crossGen->score() << " ";
      else
	mesaux << "0.0 ";

      if(mutGen)
	mesaux << mutGen->score();
      else
	mesaux << "0.0";

    }
    else
      mesaux << it->second->getScore();

    mesaux << "  " << it->second->getTechnique() << "  " << it->second->getIdDad() << " ";
    mesaux << it->second->getIslandDad() << "  " << it->second->getIdMom() << " " << it->second->getIslandMom() << "\n";

    rankList = it->second->getRankingList();

    for(rankIt = rankList.begin(); rankIt != rankList.end(); rankIt++)
      mesaux << *rankIt << " ";

    mesaux << "\n";

    if(exType == GAGenealogyMemory::ALL) {
      if(crossGen)
	mesaux << *crossGen;
      else
	mesaux << *mutGen;
    }
    else if(exType == GAGenealogyMemory::FULL) {
      if(crossGen)
	mesaux << *crossGen;
      else
	mesaux << "#No genome crossover";

      mesaux << "\n";

      if(mutGen)
	mesaux << *mutGen;
      else
	mesaux << "#No genome mutation";
    }

    GALogger::instance()->appendLogMessage("", mesaux.str(), GALogger::only_stats);

  }

}


GAGenomeNode::GAGenomeNode(const GAGenome &g, GAGenome *d, GAGenome *m, int fg, int lg, GAGenealogyMemory::executionType extype, bool mutated) {

  // If the type of genealogy is ALL or FULL, then we must store a copy of
  // the genome (in the crossGen var if this is the result of a crossover or
  // in the mutGen if this is the result of a mutation)
  if(extype == GAGenealogyMemory::ALL || extype == GAGenealogyMemory::FULL) {
    if(mutated) {
      mutGen = g.clone();
      mutGen->evaluate();
      crossGen = NULL;
    }
    else {
      crossGen = g.clone();
      crossGen->evaluate();
      mutGen = NULL;
    }
  }
  else {
    crossGen = NULL;
    mutGen = NULL;
  }

  // Copy information about island of origin, score and
  // generations the individual is alive
  island = g.getIsland();
  score  = g.score();
  firstG = fg;
  lastG  = lg;

  // If we have a dad, copy its information
  if(d) {
    idDad = d->getId();
    isDad = d->getIsland();
  }
  else {
    idDad = 0;
    isDad = -1;
  }

  // If we have a mom, copy its information
  if(m && (m != d)) {
    idMom = m->getId();
    isMom = m->getIsland();
  }
  else {
    idMom = 0 ;
    isMom = -1;
  }

  // Retrieve the technique the individual has been created with
  const MOSGenome& mosg = dynamic_cast<const MOSGenome&> (g);
  tech = mosg.getTechniqueId();

}


// Delete the genomes created. It is here and not in the destructor because it brings problems
void GAGenomeNode::deleteGenomes() {

  if(crossGen)
    delete crossGen;

  if(mutGen)
    delete mutGen;

  crossGen = NULL;
  mutGen   = NULL;

}


void GAGenomeNode::setMutatedGenome(GAGenome &g, GAGenealogyMemory::executionType extype) {

  // In case genealogy is of type ALL, we must erase the intermmediate genome
  if(extype == GAGenealogyMemory::ALL)
    if(crossGen) {
      delete crossGen;
      crossGen = NULL;
    }

  mutGen = g.clone();
  mutGen->evaluate();

  score = mutGen->score();

}
