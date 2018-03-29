/*****************************************************************************
 * GAEDAlib: A C++ GA library with EDA and multiprocessor (MPI) support      *
 *                                                                           *
 * (C) 2005 Pedro Diaz (pdiaz@laurel.datsi.fi.upm.es)                        *
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

// Author: Ramon Etxeberria
// Date:   1999-12-09

#ifndef MATRIX_H
#define MATRIX_H

/* INCLUDES */

#include <iostream>
#include <set>

#include "std_stream.h"

#ifndef CSet
 #define CSet std::set<int>
#endif


/* CLASS DECLARATION */

class Matrix  {

   public:

      Matrix ();
      Matrix (int rows, int columns);
      Matrix (long double** matrix, int rows, int columns);
      Matrix (const Matrix& matrix);

      virtual ~Matrix ();

      Matrix operator+ (Matrix matrix);
      Matrix operator- (Matrix matrix);
      Matrix operator* (Matrix matrix);

      Matrix operator* (long double num);
      Matrix operator/ (long double num);

      long double& Item (int row, int column);
      long double* operator[] (int row);

      int rows () {return m_rows;}
      int cols () {return m_columns;}

      Matrix RowPartition (CSet& nodes);
      Matrix Partition    (CSet& x_set, CSet& y_set);

      long double Determinant ();
      long double Adjunt (int row,int column);
      Matrix      Adjunt ();

      Matrix Transpose ();
      Matrix Invert    ();

      friend STD_OSTREAM& operator<< (STD_OSTREAM& os, Matrix& matrix);
      friend STD_ISTREAM& operator>> (STD_ISTREAM& os, Matrix& matrix);


   protected:

      int m_rows;
      int m_columns;
      long double** m_data;

};

#endif
