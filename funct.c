/*****************************************************************************
 *  This file is part of the eccd_code program.                              *  
 *  Copyright (C) 2015 Xin Li <lixin.reco@gmail.com>                         *  
 *                                                                           *
 *  Filename:  funct.c                                                       *
 *  Function:  functions for memory handeling and matrix solving             *
 *  Updated:   2015-Mar-24                                                   *
 *  License:   GNU Public License, version 2                                 *
 *                                                                           *  
 *  This program is free software; you can redistribute it and/or modify     *  
 *  it under the terms of the GNU General Public License as published by     *  
 *  the Free Software Foundation; either version 2 of the License, or        *  
 *  (at your option) any later version.                                      *  
 *                                                                           *
 *  This program is distributed in the hope that it will be useful,          *  
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *  
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *  
 *  GNU General Public License for more details.                             *  
 *                                                                           *
 *  You should have received a copy of the GNU General Public License along  *  
 *  with this program; if not, write to the Free Software Foundation, Inc.,  *  
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.              *  
 *****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>

void* my_malloc(size_t bytes) 
{
   void* ptr = malloc(bytes);
   if(ptr == NULL) 
   {
      printf ("Error in allocating array!\n");
      exit(1);
   } 
   else 
   {
     return ptr;
   }
}

void my_eigen_symmv(gsl_matrix* data, int DIM,
                    gsl_vector* eval, gsl_matrix* evec)
{
   if ( DIM <= 0 )
   {
      printf ("Error: DIM is not larger than 0 in my_eigen_symmv!\n");
      exit(1);
   }
   else if ( data[0].size1 != DIM || data[0].size2 != DIM )
   {
      printf ("Error: dimension of matrix is not equal to DIM in my_eigen_symmv!\n");
      exit(1);
   }

   /* make a copy of 'data': 'matrix' */
   /* NOTE: 'data' will be destroyed after gsl_eigen_symmv */

   gsl_matrix *matrix = gsl_matrix_alloc (DIM, DIM);
   gsl_matrix_memcpy(matrix, data);

   /* calculate eigenvectors for real symmetric matrix */

   gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (DIM);
   gsl_eigen_symmv (data, eval, evec, w);
   gsl_eigen_symmv_free (w);

   /* sort eigenvalues and eigenvectors */

   gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
   
   /* NOTE: 'data' was destroyed after gsl_eigen_symmv */
   /* copy original elements back to 'data' */

   gsl_matrix_memcpy(data, matrix);

   /* check that 'data' can be diagonalized by 'evec' */

   {
      gsl_matrix *pd_1 = gsl_matrix_alloc (DIM, DIM);
      gsl_matrix_set_zero (pd_1);

      gsl_matrix *pd_2 = gsl_matrix_alloc (DIM, DIM);
      gsl_matrix_set_zero (pd_2);

      gsl_blas_dgemm (CblasTrans, CblasNoTrans,
                      1.0, evec, data,
                      0.0, pd_1);

      gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                      1.0, pd_1, evec,
                      0.0, pd_2);

      /* compare 'evec'(T) * 'data' * 'evec' with diagonal matrix */

      int i,j;
      gsl_matrix *diag = gsl_matrix_alloc (DIM, DIM);
      gsl_matrix_set_zero (diag);
      for (i=0; i<DIM; i++)
      {
         gsl_matrix_set(diag, i, i, gsl_vector_get(eval, i));
      }

      for (i=0; i<DIM; i++)
      for (j=0; j<DIM; j++)
      {
         if (fabs(gsl_matrix_get(pd_2, i, j) - 
                  gsl_matrix_get(diag, i, j)) > 1.0e-06)
         {
            printf ("Error: 'evec'(T) * 'data' * 'evec' is not diagonal!\n");
            exit(1);
         }
      }
   }
}

void my_prog_die(char* txt)
 {
   printf("%s\n", txt);
   exit(1);
 }

void my_norm_vec(gsl_vector* x, int DIM)
 {
   if ( x[0].size != DIM )
   {
      my_prog_die("Error: dimension of x is not equal to DIM in my_norm_vec!\n");
   }

 	int i;
 	double norm = gsl_blas_dnrm2(x);
   for (i=0; i<DIM; i++)
   {
      gsl_vector_set(x, i, gsl_vector_get(x, i)/norm);
   }
 }

void my_outer_prod(gsl_vector* x, gsl_vector* y, int DIM,
                   gsl_vector* z)
{
   if ( x[0].size != DIM || y[0].size != DIM )
   {
      my_prog_die("Error: dimension of x or y is not equal to DIM in my_outer_prod!\n");
   }

   int i;
   for (i=0; i<3; i++)
   {
      double x1 = gsl_vector_get(x, 0);
      double x2 = gsl_vector_get(x, 1);
      double x3 = gsl_vector_get(x, 2);

      double y1 = gsl_vector_get(y, 0);
      double y2 = gsl_vector_get(y, 1);
      double y3 = gsl_vector_get(y, 2);

      gsl_vector_set(z, 0, x2*y3-x3*y2);
      gsl_vector_set(z, 1, x3*y1-x1*y3);
      gsl_vector_set(z, 2, x1*y2-x2*y1);
   }
}
