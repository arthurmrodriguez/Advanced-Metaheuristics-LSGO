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

#include "GAGenealogy.h"

GAGenealogy* GAGenealogy::instance_ptr = NULL;

GAGenealogy::GAGenealogy(int idis) : actGen(0), island(idis), id(1) {
}

GAGenealogy::~GAGenealogy() {

  if (instance_ptr)
    instance_ptr = NULL;

}

GAGenealogy* GAGenealogy::handle() {
  return instance_ptr;
}
