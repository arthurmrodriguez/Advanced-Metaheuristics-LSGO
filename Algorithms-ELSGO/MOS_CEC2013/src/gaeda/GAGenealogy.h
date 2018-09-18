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

#ifndef GAGENEALOGY_H
#define GAGENEALOGY_H

#include <sstream>
#include <vector>

#include "gaid.h"

class GAGenome;

class GAGenealogy : public GAID {

 public:

  GADefineIdentity ("Genealogy", GAID::Genealogy);

  // Retrieves the pointer to the singleton object
  static GAGenealogy* handle();

  // Destructor
  virtual ~GAGenealogy();

  // Check for type of genealogy being used
  bool isGenealogyMemory() {return classID() == GAID::GenealogyMemory;}
  bool isGenealogyTracer() {return classID() == GAID::GenealogyTracer;}

  // Set and get current generation
  void  setGeneration(unsigned int gen) {actGen = gen;}
  const unsigned getGeneration() const {return actGen;}

  // Updates the generation counter
  void newGeneration() {actGen++;}

  // Returns current ID for a genome
  const unsigned long int getId() const {return id;}

  // Returns island ID
  const int getIsland() const {return island;}

  // Virtual methods to add information to the genealogy
  virtual void addNode(GAGenome &gen) = 0;
  virtual void familyRelationship(GAGenome &dad, GAGenome &mom, GAGenome &child1, GAGenome &child2) = 0;
  virtual void familyRelationship(std::vector<GAGenome*>& parents, std::vector<GAGenome*>& children) = 0;
  virtual void mutated(GAGenome &dad, GAGenome &child, bool crossAndMut) = 0;
  virtual void deceased(GAGenome &gen) = 0;
  virtual void addImmigrant(GAGenome &gen) = 0;
  virtual void bestGenome(GAGenome &gen) = 0;

 protected:

  // Protected constructor
  GAGenealogy(int idis);

  // Returns a new ID for a genome
  unsigned long int newId() {return id++;}

  static GAGenealogy* instance_ptr; // Protected pointer to the singleton instance
  unsigned actGen;                    // Current generation
  int island;                         // Island ID
  unsigned long int id;               // ID for genemes in the genealogy

};

#endif
