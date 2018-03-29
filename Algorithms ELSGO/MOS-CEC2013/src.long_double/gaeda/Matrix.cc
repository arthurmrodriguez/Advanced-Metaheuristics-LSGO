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


/* Matrix.cc: implementation of the Matrix class */

/* INCLUDES */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "Matrix.h"


/* DEFINE */

// 0 = LU decomposition
// 1 = Adjunts.

#define DETERMINANT_CALCU 0


/* AUX FUNCTIONS */
void ludcmp (int nn, long double* matrix, long double* d);


Matrix::Matrix () : m_rows (0), m_columns (0), m_data (NULL) {

}


Matrix::Matrix (int rows, int columns): m_rows (rows), m_columns (columns) {

   if (m_rows == 0 || m_columns == 0) {

      m_data = NULL;
      return;

   }

   m_data = new long double* [m_rows];

   for (int i = 0; i < m_rows; i++) {

      m_data [i] = new long double [m_columns];

      for (int j = 0; j < m_columns; j++)
         m_data [i][j] = 0;

   }

}


Matrix::Matrix (long double** matrix, int rows, int columns): m_rows (rows), m_columns (columns) { // TODO

   // With this constructor a IND_SIZE x IND_SIZE matrix
   // is assumed.

   m_data = new long double*[ m_rows];

   for (int i = 0; i < m_rows; i++) {

      m_data [i] = new long double [m_columns];

      for (int j = 0; j < m_columns; j++)
         m_data [i][j] = matrix [i][j];

   }

}


Matrix::Matrix (const Matrix& matrix) : m_rows (matrix.m_rows), m_columns (matrix.m_columns) {

   if (m_rows == 0 || m_columns == 0) {

      m_data = NULL;
      return;

   }

   m_data = new long double* [m_rows];

   for (int i = 0; i < m_rows; i++) {

      m_data [i] = new long double [m_columns];

      for (int j = 0; j < m_columns; j++)
         m_data [i][j] = matrix.m_data [i][j];

   }

}


Matrix::~Matrix () {

   if (m_rows == 0 || m_columns == 0)
      return;

   for (int i = 0; i < m_rows; i++)
      delete [] m_data [i];

   delete [] m_data;

}


Matrix Matrix::operator+ (Matrix matrix) {

   Matrix result (m_rows, m_columns);

   for (int i = 0; i < m_rows; i++)
      for (int j = 0; j < m_columns; j++)
         result.m_data [i][j] = m_data [i][j] + matrix.m_data [i][j];

   return result;

}


Matrix Matrix::operator- (Matrix matrix) {

   Matrix result (m_rows, m_columns);

   for (int i = 0; i < m_rows; i++)
      for (int j = 0; j < m_columns; j++)
         result.m_data [i][j] = m_data [i][j] - matrix.m_data [i][j];

   return result;

}


Matrix Matrix::operator* (Matrix matrix) {

   Matrix result (m_rows, matrix.m_columns);

   for (int i = 0; i < m_rows; i++)
      for (int j = 0; j < matrix.m_columns; j++)
         for (int k = 0; k < m_columns; k++)
            result.m_data [i][j] += m_data [i][k] * matrix.m_data [k][j];

   return result;

}


Matrix Matrix::operator* (long double num) {

   Matrix multiplied (m_rows, m_columns);

   for (int i = 0; i < m_rows; i++)
      for (int j = 0; j < m_columns; j++)
         multiplied.m_data [i][j] = m_data [i][j] * num;

   return multiplied;

}


Matrix Matrix::operator/ (long double num) {

   // TODO: Comprobar divisiones por cero
   Matrix divided (m_rows, m_columns);

   for (int i = 0; i < m_rows; i++)
      for (int j = 0; j < m_columns; j++)
         divided.m_data [i][j] = m_data [i][j] / num;

   return divided;

}


Matrix Matrix::RowPartition (CSet& _set) {

   Matrix submatrix (_set.size (), m_columns);

   int row = 0;

   for (int i = 0; i < m_rows; i++)
      if(_set.count (i)) {

         for (int j = 0; j < m_columns; j++)
            submatrix.m_data [row][j] = m_data [i][j];

         row++;

      }

   return submatrix;

}


Matrix Matrix::Partition (CSet& x_set, CSet& y_set) {

   Matrix submatrix (x_set.size (), y_set.size ());

   int row = 0;

   for (int i = 0; i < m_rows; i++)
      if (x_set.count (i)) {

         int column = 0;

         for (int j = 0; j < m_columns; j++)
            if (y_set.count (j)) {

               submatrix.m_data [row][column] = m_data [i][j];
               column++;

            }

         row++;

      }

   return submatrix;

}


long double& Matrix::Item (int row, int column) {

   return m_data [row][column];

}


long double* Matrix::operator[] (int row) {

   return m_data [row];

}


long double Matrix::Determinant () {

   int i;

   if (m_columns == 0)
      return 0;

   if (m_columns == 1)
      return m_data [0][0];

   if (m_columns == 2)
      return m_data [0][0] * m_data [1][1] -
             m_data [0][1] * m_data [1][0];


   //LU decomposition.

   if (DETERMINANT_CALCU == 0) {

      long double* m_data2;
      long double aux;
      long double d;

      m_data2 = (long double*) calloc (m_rows * m_columns, sizeof (long double));

      if (m_data2 == NULL) {

         STD_CERR << STD_ENDL << "Error: No dynamic memory allocation done";
         exit (0);

      }

      for (i = 0; i < m_rows; i++)
         for (int j = 0; j < m_columns; j++)
            m_data2 [i * m_rows + j] = (long double) m_data [i][j];

      ludcmp (m_columns, m_data2, &d); // PLEASE note that matrix is modified by ludcmp() since
                                       // the decomposition of the matrix is returned in matrix itself.
      aux = (long double) d;

      for (i = 0; i < m_columns; i++) {

         if (m_data2 [i * m_rows + i] == (long double) 0.0) {

            STD_CERR << STD_ENDL << "Error: Matrix determinant is zero";
            exit (0);

         }

         aux *= m_data2 [i * m_rows + i];

      }

      free (m_data2);

      return (aux);

   }
   else {

      int column = 0;

      for (i = 0; i < m_columns && m_data [0][i] == 0; i++)
         column++;

      if (column == m_columns)
         return 0;

      for (i = column + 1; i < m_columns; i++) {

         long double k = m_data [0][i] / m_data [0][column];

         m_data [0][i] = 0;

         for (int j = 1; j < m_rows; j++) 
            m_data [j][i] -= m_data [j][column] * k;

      }

      return m_data [0][column] * Adjunt (0, column);

   }

}


long double Matrix::Adjunt (int row, int column) {

   int i;
   CSet rows;
   CSet columns;

   for (i = 0; i < m_rows; i++)
      if (i != row)
         rows.insert (i);

   for (i = 0; i < m_columns; i++)
      if (i != column)
         columns.insert (i);

   Matrix minor_matrix = Partition (rows, columns);
   long double minor = minor_matrix.Determinant ();

   if ((row + column) % 2 == 0)
      return minor;

   else return -minor;

}


Matrix Matrix::Adjunt () {

   Matrix adjunt (m_rows, m_columns);

   for (int i = 0; i < m_rows; i++)
      for (int j = 0; j < m_columns; j++)
         adjunt.m_data [i][j] = Adjunt (i, j);

   return adjunt;
}


Matrix Matrix::Transpose () {

   Matrix transposed (m_columns, m_rows);

   for (int i = 0; i < m_rows; i++)
      for (int j = 0; j < m_columns; j++)
         transposed.m_data [j][i] = m_data [i][j];

   return transposed;

}


Matrix Matrix::Invert () {

   Matrix adjunt = Adjunt ();
   Matrix t_adjunt = adjunt.Transpose ();
   long double determinant = Determinant ();

   //jmp if (determinant ==0) {
   //jmp           cout << "determinante a 0!!! sustituyendolo por 1.0E-30" << STD_ENDL;
   //jmp           determinant=1.0E-30;
   //jmp }

   Matrix inverted = t_adjunt / determinant;

   return inverted;

}


STD_OSTREAM& operator<< (STD_OSTREAM& os, Matrix& matrix) {

   os << " rows: " << matrix.m_rows << " columns: " << matrix.m_columns << STD_ENDL;

   for (int i = 0; i < matrix.m_rows; i++) {

      os << "{ " ;

      for (int j = 0; j < matrix.m_columns - 1; j++)
         os << matrix.m_data [i][j] << ", ";

      os << matrix.m_data [i][matrix.m_columns - 1] << " }" << STD_ENDL;

   }

   return os;

}


STD_ISTREAM& operator>> (STD_ISTREAM& os, Matrix& matrix) {

   int xdim, ydim;
   char TEXT [50];
   char KAR;

   os >> TEXT >> xdim >> TEXT >> ydim;

   //Create the matrix
   Matrix* tempmatrix = new Matrix (xdim, ydim);

   //matrix = tempmatrix;
   //matrix->m_rows = xdim;
   //matrix->m_columns = ydim;

   for (int i = 0; i < xdim; i++) {

      os >> KAR;

      for (int j = 0; j < ydim - 1; j++)
         os >> tempmatrix->m_data [i][j] >> KAR;

      os >> tempmatrix->m_data [i][ydim - 1] >> KAR;

   }

   matrix = (*tempmatrix);

   return os;

}


/* AUX FUNCTIONS */

void ludcmp (int nn, long double* matrix, long double* d) {

   int i, j, k, imax = 0;
   long double big, dum, sum, temp;
   long double* vv;
   long double TINY = (long double) 0.000001;

   vv = (long double*) calloc (nn, sizeof (long double));

   if (vv == NULL) {

      STD_CERR << STD_ENDL << "Error: No dynamic memory allocation done";
      exit (0);

   }

   *d = (long double) 1.0;

   for (i = 0; i < nn; i++) {

      big = (long double) 0.0;

      //It starts in j=i because the matrix is simmetric.
      for (j = i; j < nn; j++) {

         temp = (long double) fabs (matrix [i * nn + j]);

         if (temp > big)
            big = temp;

      }

      if (big == (long double) 0.0) {

         STD_CERR << STD_ENDL << "Error: singular matrix in routine ludcmp";
         exit (0);

      }

      vv [i] = (long double) ((long double) 1.0 / big);

   }

   for (j = 0; j < nn; j++) {

      for (i = 0; i < j; i++) {

         sum = matrix [i * nn + j];

         for (k = 0; k < i; k++)
            sum -= matrix [i * nn + k] * matrix [k * nn + j];

         matrix [i * nn + j] = sum;

      }

      big = (long double) 0.0;

      for (i = j; i < nn; i++) {

         sum = matrix [i * nn + j];

         for (k = 0; k < j; k++)
            sum -= matrix [i * nn + k] * matrix [k * nn + j];

         matrix [i * nn + j] = sum;

         dum = vv [i] * (long double) fabs (sum);

         if (dum >= big) {

            big = dum;
            imax = i;

         }

      }

      if (j != imax) {

         for (k = 0; k < nn; k++) {

            dum = matrix [imax * nn + k];
            matrix [imax * nn + k] = matrix [j * nn + k];
            matrix [j * nn + k] = dum;

         }

         *d = -(*d);
         vv [imax] = vv [j];

      }

      if (matrix [j * nn + j] == (long double) 0.0)
         matrix [j * nn + j] = TINY;

      if (j != nn) {

         dum = (long double) ((long double) 1.0 / matrix [j * nn + j]);

         for (i = j + 1; i < nn; i++)
            matrix [i * nn + j] *= dum;

      }

   }

   free (vv);

}
