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

#ifndef GAGENEALOGYMEM_H
#define GAGENEALOGYMEM_H

#include <map>
#include <list>

#include "GAGenealogy.h"

class GAGenomeNode;
class GAGenome;
class GAPopulation;

class GAGenealogyMemory : public GAGenealogy {

 public:

  GADefineIdentity ("GenealogyMemory", GAID::GenealogyMemory);

  typedef enum {NORMAL, ALL, FULL} executionType;
  typedef void (*heritFunc) (GAGenome&, GAGenome&);

  static GAGenealogy* create (int idis, executionType extype = GAGenealogyMemory::NORMAL);

  ~GAGenealogyMemory();

  void addNode(GAGenome &gen);
  void familyRelationship(GAGenome &dad, GAGenome &mom, GAGenome &child1, GAGenome &child2);
  void familyRelationship(std::vector<GAGenome*>& parents, std::vector<GAGenome*>& children) {}
  void mutated(GAGenome &dad, GAGenome &child, bool crossAndMut);
  void deceased(GAGenome &gen);
  void addImmigrant(GAGenome &gen);
  void bestGenome(GAGenome &gen) {}

  void updateRankings(const GAPopulation *pop);
  void addGenomeRanking(GAGenome &gen, unsigned int ranking) const;
  void addGenomeRanking(unsigned long int idGen, int idIs, unsigned int ranking) const;

  bool exist(GAGenome &gen);

  void setExecutionType(executionType ext) {exType = ext;}

  heritFunc setHeritage(heritFunc hf){return funcHeritage=hf;}
  heritFunc getHeritage() {return funcHeritage;}
  void heritage(GAGenome &gen1, GAGenome &gen2);

  GAGenomeNode *getGeneNode(GAGenome &gen) const;
  GAGenomeNode *getGeneNode(unsigned long int idGen, int idIs) const;

  //Get the bestN number of genomes of the lastGen generations per technique. Return value is lastBests
  void getLastBests(int *lastBests, int *techList, int nAlgs, int lastGen, int bestN) const;

  //Print the whole genealogy
  void print();

 protected:

  GAGenealogyMemory(int idis, executionType extype = GAGenealogyMemory::NORMAL);

  executionType exType;                                 // Type of execution. Store a copy of the genome or not
  std::multimap<unsigned long int, GAGenomeNode*> genealogy; // Genealogy in a hash table
  heritFunc funcHeritage;                               // Heritage function

};


// Class: Node for the hash table
class GAGenomeNode {

 public:

  GAGenomeNode(const GAGenome &g, GAGenome *d, GAGenome *m, int fg, int lg, GAGenealogyMemory::executionType extype, bool mutated);
  ~GAGenomeNode() {rankList.clear();}

  // Delete the genomes created. It is here and not in the destructor because it brings problems
  void deleteGenomes();

  const GAGenome *getCrossoverGenome() const {return crossGen;}
  const GAGenome *getMutationGenome () const {return mutGen;}

  GAGenome *getCrossoverGenome() {return crossGen;}
  GAGenome *getMutationGenome () {return mutGen;}

  void setMutatedGenome(GAGenome &g, GAGenealogyMemory::executionType extype);

  const unsigned long int getIdMom() const {return idMom;}
  const unsigned long int getIdDad() const {return idDad;}

  const int getIslandDad() const {return isDad;}
  const int getIslandMom() const {return isMom;}

  const int getFirstG() const {return firstG;}
  const int getLastG() const {return lastG;}

  const long double getScore() const {return score;}
  const int getIsland() const {return island;}
  const int getTechnique() const {return tech;}

  const std::list<signed short> getRankingList() const {return rankList;}

  void setScore(long double sc) {score = sc;}
  void setFirstG(unsigned int first) {firstG = first;}
  void setLastG(unsigned int last) {lastG = last;}
  void addRanking(signed short ranking) {rankList.push_back(ranking);}

 private:

  GAGenome *crossGen;               // Genome after crossover
  GAGenome *mutGen;                 // Genome after mutation
  int island;                       // Island of creation of the genome
  unsigned long int idDad;          // Dad's ID
  unsigned long int idMom;          // Mom's ID
  int isDad;                        // Dad's island ID
  int isMom;                        // Mom's island ID
  long double score;                     // Score for the genome
  unsigned firstG;                  // First generation (when it was created)
  unsigned lastG;                   // Last generation (when it was deceased)
  int tech;                         // Technique that created the individual
  std::list<signed short> rankList; // List of rankings, one per generation

};

#endif
