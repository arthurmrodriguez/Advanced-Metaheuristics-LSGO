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

#ifndef GAGENEALOGYTR_H
#define GAGENEALOGYTR_H

#include "GAGenealogy.h"

class GAGenealogyTracer : public GAGenealogy {

 public:

  GADefineIdentity ("GenealogyTracer", GAID::GenealogyTracer);

  static GAGenealogy* create (int idis);

  ~GAGenealogyTracer() {;}

  void addNode(GAGenome &gen);
  void familyRelationship(GAGenome &dad, GAGenome &mom, GAGenome &child1, GAGenome &child2);
  void familyRelationship(std::vector<GAGenome*>& parents, std::vector<GAGenome*>& children);
  void mutated(GAGenome &dad, GAGenome &child, bool crossAndMut);
  void deceased(GAGenome &gen);
  void addImmigrant(GAGenome &gen);
  void bestGenome(GAGenome &gen);

 protected:

  GAGenealogyTracer(int idis);

};

#endif
