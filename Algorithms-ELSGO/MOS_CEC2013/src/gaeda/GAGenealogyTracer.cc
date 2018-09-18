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

#include "GAGenealogyTracer.h"

#include "logger/GALogger.h"
#include "genomes/GAGenome.h"
#include "genomes/MOSGenome.h"

GAGenealogy* GAGenealogyTracer::create (int idis) {

  assert (instance_ptr == NULL);

  return instance_ptr = new GAGenealogyTracer (idis);

}

GAGenealogyTracer::GAGenealogyTracer(int idis) : GAGenealogy(idis) {;}

void GAGenealogyTracer::addNode(GAGenome &gen) {

  if(gen.getId() == 0) {
    gen.setId(newId());
    gen.setIsland(island);
  }

  // Retrieve the technique the individual has been created with
  unsigned tech = (dynamic_cast<const MOSGenome&> (gen)).getTechniqueId();

  stringstream mes;
  mes << "New Genome(" << actGen << "): " << gen.getId() << " " << gen.getIsland() << ", " << gen.score() << " " << tech;
  GALogger::instance()->appendLogMessage("", mes.str(), GALogger::only_stats);

}


void GAGenealogyTracer::familyRelationship(GAGenome &dad, GAGenome &mom, GAGenome &child1, GAGenome &child2) {

  child1.setId(newId());
  child1.setIsland(island);

  if(&child1 != &child2) {
    child2.setId(newId());
    child2.setIsland(island);
  }

  // Retrieve the technique the individuals have been created with
  unsigned tech_child1 = (dynamic_cast<const MOSGenome&> (child1)).getTechniqueId();
  unsigned tech_child2 = (dynamic_cast<const MOSGenome&> (child2)).getTechniqueId();
  unsigned tech_dad = (dynamic_cast<const MOSGenome&> (dad)).getTechniqueId();
  unsigned tech_mom = (dynamic_cast<const MOSGenome&> (mom)).getTechniqueId();

  stringstream mes;
  mes << "Children(" << actGen << "): " << child1.getId() << " " << child1.getIsland() << ", " << child1.score() << " " << tech_child1 << "; ";
  mes << child2.getId() << " " << child2.getIsland() << ", " << child2.score()  << " " << tech_child2;
  mes << " | dad: " << dad.getId() << " " << dad.getIsland() << ", " << dad.score() << " " << tech_dad;
  mes << "; mom: " << mom.getId() << " " << mom.getIsland() << ", " << mom.score()  << " " << tech_mom;
  GALogger::instance()->appendLogMessage("", mes.str(), GALogger::only_stats);

}


void GAGenealogyTracer::familyRelationship(std::vector<GAGenome*>& parents, std::vector<GAGenome*>& children) {

  stringstream mes;
  mes << "Children(" << actGen << "): ";

  for (unsigned i = 0; i < children.size(); i++) {

    children[i]->setId(newId());
    children[i]->setIsland(island);

    // Retrieve the technique the individual has been created with
    unsigned tech = (dynamic_cast<const MOSGenome&> (*(children[i]))).getTechniqueId();

    mes << children[i]->getId() << " " << children[i]->getIsland() << ", " << children[i]->score() << " " << tech;

    if (i != children.size() - 1)
      mes << "; ";

  }

  mes << " | ";

  for (unsigned i = 0; i < parents.size(); i++) {

    // Retrieve the technique the individual has been created with
    unsigned tech = (dynamic_cast<const MOSGenome&> (*(parents[i]))).getTechniqueId();

    mes << "parent: " << parents[i]->getId() << " " << parents[i]->getIsland() << ", " << parents[i]->score() << " " << tech;

    if (i != parents.size() - 1)
      mes << "; ";

  }

  GALogger::instance()->appendLogMessage("", mes.str(), GALogger::only_stats);

}


void GAGenealogyTracer::mutated(GAGenome &dad, GAGenome &child, bool crossAndMut) {

  if(!crossAndMut) {
    child.setId(newId());
    child.setIsland(island);
  }

  // Retrieve the technique the individuals have been created with
  unsigned tech_child = (dynamic_cast<const MOSGenome&> (child)).getTechniqueId();
  unsigned tech_dad = (dynamic_cast<const MOSGenome&> (dad)).getTechniqueId();

  stringstream mes;
  mes << "Mutated(" << actGen << "): " << child.getId() << " " << child.getIsland() << ", " << child.score() << " " << tech_child;
  mes << " | dad: " << dad.getId() << " " << dad.getIsland() << ", " << dad.score() << " " << tech_dad;
  GALogger::instance()->appendLogMessage("", mes.str(), GALogger::only_stats);

}


void GAGenealogyTracer::deceased(GAGenome &gen) {

  // Retrieve the technique the individual has been created with
  unsigned tech = (dynamic_cast<const MOSGenome&> (gen)).getTechniqueId();

  if(gen.getId() != 0) {
    stringstream mes;
    mes << "Deceased(" << actGen << "): " << gen.getId() << " " << gen.getIsland() << ", " << gen.score() << " " << tech;
    GALogger::instance()->appendLogMessage("", mes.str(), GALogger::only_stats);
  }

}


void GAGenealogyTracer::addImmigrant(GAGenome &gen) {

  // Retrieve the technique the individual has been created with
  unsigned tech = (dynamic_cast<const MOSGenome&> (gen)).getTechniqueId();

  stringstream mes;
  mes << "New Immigrant(" << actGen << "): " << gen.getId() << " " << gen.getIsland() << ", " << gen.score() << " " << tech;
  GALogger::instance()->appendLogMessage("", mes.str(), GALogger::only_stats);

}


void GAGenealogyTracer::bestGenome(GAGenome &gen) {

  // Retrieve the technique the individual has been created with
  unsigned tech = (dynamic_cast<const MOSGenome&> (gen)).getTechniqueId();

  if(gen.getId() != 0) {
    stringstream mes;
    mes << "Best(" << actGen << "): " << gen.getId() << " " << gen.getIsland() << ", " << gen.score() << " " << tech;
    GALogger::instance()->appendLogMessage("", mes.str(), GALogger::only_stats);
  }

}
