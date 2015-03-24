/*****************************************************************************
 *  This file is part of the eccd_code program.                              *  
 *  Copyright (C) 2015 Xin Li <lixin.reco@gmail.com>                         *  
 *                                                                           *
 *  Filename:  tresp.h                                                       *
 *  Function:  fit transition electrostatic potential                        *
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
#include <gsl/gsl_linalg.h>

#include "funct.h"

int main( int argc, char *argv[] )
{
   /* define variables */

   FILE *file_inp, *file_cstr;
   char line[256],inpESP[128],inpCstr[128];
   double **A, *B;
   double **x, *q;
   int nAtoms, nPoints, nCstr, i, j, k;

   /* read esp file: xxx.esp */

   sprintf(inpESP, "%s", argv[1]);
   file_inp = fopen(inpESP, "r") ;
   if (file_inp == NULL)
   {
      my_prog_die("Error: Cannot open ESP input file !\n") ;
   }

   /* read first line: N_atoms */

   fgets( line, sizeof( line ), file_inp );
   sscanf( line, "%d", &nAtoms );
   printf ("nAtoms = %d\n",nAtoms);

   /* allocate arrays */

   x = my_malloc( sizeof(double *) * nAtoms );
   A = my_malloc( sizeof(double *) * nAtoms );
   B = my_malloc( sizeof(double) * nAtoms );
   q = my_malloc( sizeof(double) * nAtoms );

   for ( i=0; i<nAtoms; i++ ) 
   {
      x[i] = my_malloc( sizeof(double) * 3 );
      A[i] = my_malloc( sizeof(double) * nAtoms );

      /* initialize A and B */

      B[i] = 0.0;
      for (j=0; j<nAtoms; j++)
      {
         A[i][j] = 0.0;
      }
   }

   /* read molecule coordinates */

   for (i=0; i<nAtoms; i++)
   {
      fgets( line, sizeof( line ), file_inp );
      sscanf( line, "%lf%lf%lf", &x[i][0], &x[i][1], &x[i][2] );
   }

   /* read N_points */

   fgets( line, sizeof( line ), file_inp );
   sscanf( line, "%d", &nPoints );
   printf ("nPoints = %d\n",nPoints);

   /* read ESP coordinate and potential */
   /* j&k for atoms, i for points */

   for (i=0; i<nPoints; i++)
   {
      double *xi, esp;
      xi = my_malloc(sizeof(double) * 3);

      fgets( line, sizeof( line ), file_inp );
      sscanf( line, "%lf%lf%lf%le", &xi[0], &xi[1], &xi[2], &esp );

      /* contribution to B */

      for (k=0; k<nAtoms; k++)
      {
         double rik;
         rik = sqrt((xi[0] - x[k][0]) * (xi[0] - x[k][0]) + 
                    (xi[1] - x[k][1]) * (xi[1] - x[k][1]) + 
                    (xi[2] - x[k][2]) * (xi[2] - x[k][2]));
         B[k] += esp / rik;

         /* contribution to A */

         for (j=k; j<nAtoms; j++)
         {
            double rij;
            rij = sqrt((xi[0] - x[j][0]) * (xi[0] - x[j][0]) + 
                       (xi[1] - x[j][1]) * (xi[1] - x[j][1]) + 
                       (xi[2] - x[j][2]) * (xi[2] - x[j][2]));
            A[j][k] += 1.0 / (rij * rik);
         }
      }

      free(xi);
   }

   for (k=0; k<nAtoms; k++)
   {
      for (j=k; j<nAtoms; j++)
      {
         A[k][j] = A[j][k];
      }
   }

   fclose(file_inp);

   /* use GSL to solve linear equations */

   {
      /* read constraint file: xxx.cstr */
      
	   sprintf(inpCstr, "%s", argv[2]);
      file_cstr = fopen(inpCstr, "r") ;
      if (file_cstr == NULL)
      {
         my_prog_die("Error: Cannot open Cstr input file!\n") ;
      }

      /* read first line: N_atoms */

      fgets( line, sizeof( line ), file_cstr );
      sscanf( line, "%d", &nCstr );
      printf ("nConstraints = %d\n",nCstr );

      int DIM = nAtoms + nCstr + 1;

      gsl_matrix *matA = gsl_matrix_alloc(DIM, DIM);
      gsl_vector *vecB = gsl_vector_alloc(DIM);

      gsl_matrix_set_zero(matA);
      gsl_vector_set_zero(vecB);

      /* assign values from A and B */

      for (k=0; k<nAtoms; k++)
      {
         gsl_vector_set(vecB, k, B[k]);
         for (j=0; j<nAtoms; j++)
         {
            gsl_matrix_set(matA, k, j, A[k][j]);
         }
      }

      /* add constraints for qtot = 0 */

      for (k=0; k<nAtoms; k++)
      {
         gsl_matrix_set(matA, nAtoms, k, 1.0);
         gsl_matrix_set(matA, k, nAtoms, 1.0);
      }

      /* add further constraints */

      int iCstr;
      for (iCstr=0; iCstr<nCstr; iCstr++)
      {
         int ai, aj;
         fgets(line, sizeof(line), file_cstr);
         sscanf(line, "%d%d", &ai, &aj);

         gsl_matrix_set(matA, nAtoms+iCstr+1, ai-1,  1.0);
         gsl_matrix_set(matA, nAtoms+iCstr+1, aj-1, -1.0);
         gsl_matrix_set(matA, ai-1, nAtoms+iCstr+1,  1.0);
         gsl_matrix_set(matA, aj-1, nAtoms+iCstr+1, -1.0);
      }

      fclose(file_cstr);

      /* solve the linear equation system */

      int ss;
      gsl_vector *qfit = gsl_vector_alloc(DIM);
      gsl_permutation *pp = gsl_permutation_alloc(DIM);
  
      gsl_linalg_LU_decomp (matA, pp, &ss);
      gsl_linalg_LU_solve (matA, pp, vecB, qfit);
  
      for (i=0; i<nAtoms; i++)
      {
         q[i] = gsl_vector_get(qfit, i);
      }
  
      gsl_permutation_free(pp);
      gsl_vector_free(qfit);

      gsl_vector_free(vecB);
      gsl_matrix_free(matA);
   }

   /* check the quality of fitting */

   {
      double *dip, qsum;
      dip = my_malloc(sizeof(double) * 3);

      dip[0] = 0.0;
      dip[1] = 0.0;
      dip[2] = 0.0;
      qsum = 0.0;

      for (i=0; i<nAtoms; i++)
      {
         printf ("q[%2d] = %20.10f\n", i+1, q[i]);
         dip[0] += x[i][0] * q[i];
         dip[1] += x[i][1] * q[i];
         dip[2] += x[i][2] * q[i];
         qsum += q[i];
      }

      printf ("dip = %10.4f%10.4f%10.4f |%15.6e\n", 
               dip[0], dip[1], dip[2], qsum);

      free(dip);
   }
  
   /* free arrays */

   for ( i=0; i<nAtoms; i++ ) 
   {
      free(x[i]);
      free(A[i]);
   }

   free(x);
   free(A);
   free(B);
   free(q);

   /* end of main program */

   return 0;
}
