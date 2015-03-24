/*****************************************************************************
 *  This file is part of the eccd_code program.                              *  
 *  Copyright (C) 2015 Xin Li <lixin.reco@gmail.com>                         *  
 *                                                                           *
 *  Filename:  eccd.c                                                        *
 *  Function:  calculate exciton-coupled circular dichroism spectra          *
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

#include "funct.h"

int main( int argc, char *argv[] )
{
   /* user-defined parameter: FWHM in eV */
   double fwhm = 0.1;

   /* define variables */

   FILE *file_inp;
   char line[256];
   int i,j,mol,nAtoms,nMols,atID,iMol,jMol;
   int *molStart, *molEnd;
   double *eneExc, **excCoup;
   double **elecDip, **totEDip, **veloDip, **totVDip, **magnDip, **totMDip;
   double **x,*q;

   double box;

   /* Constants */
   /* pi = 3.141592653589793238462643 */
   /* epsilon_0 = 8.854187817e-12 F m^-1 */
   /* ea0 = 8.47835326e-30 C m */

   /* ea0 = 2.54174636369676 Debye */
   /*     = 2.54174636369676e-18 esu cm */

   /* Bohr magneton = 9.27400968e-24 J T^-1 */
   /*               = 9.27400968e-21 erg G^-1 */

   /* ea0 x Bohr magneton = 235.721803810286e-40 erg esu cm G^-1 */

   const double PI = 3.141592653589793238462643;
   const double EaBm2cgs = 235.721803810286;

   const double ea2CGS = 2.54174636369676e-18;
   const double BM2CGS = 9.27400968e-21;

   const double DDFactor  = 9.185929163775e-39; /* cgs unit */
   const double DBMFactor = 2.296482290944e-39; /* cgs unit */

   /* DBM2theta = 24.0 PI NA / h_bar c (in cgs unit) */
   /*           = 48.0*PI*PI * 6.02214129e+23 / (6.62606957e-27 * 2.99792458e+10) */

   const double DBM2theta = 1.43620101403094e+42;

   const double au2eps   = ea2CGS * ea2CGS / DDFactor; /* a.u. to M^-1 cm^-1 (abs) */
   const double au2theta = ea2CGS * BM2CGS * DBM2theta; /* a.u. to deg cm^2 dmol^-1 (ecd) */

   /* J. A. Schellman, Chem Rev 1975, 75, 323-331. */
   const double theta2eps = PI / (4.5 * log(10.0) * 1000.0);

   /* Physical constants */
   const double Eh2eV = 27.21138505;
   const double Ehxnm = 1.0e+02/2.194746313708;
   const double eVxnm = 1.0e+04/8.06554429;

   /* read input file: snap.inp */

   file_inp = fopen( "snap.inp", "r" ) ;
   if ( NULL==file_inp )
   {
      my_prog_die( "Error: Cannot open file snap.inp !\n") ;
   }

   /* read first line: N_molecule, N_atoms, box_size*/

   fgets( line, sizeof( line ), file_inp );
   sscanf( line, "%d%d%lf", &nMols, &nAtoms, &box );

   /* allocate arrays */

   eneExc = my_malloc( sizeof(double) * nMols );
   excCoup = my_malloc( sizeof(double *) * nMols );
   for ( mol=0; mol<nMols; mol++ ) 
   {
      excCoup[mol] = my_malloc( sizeof(double) * nMols );
   }

   molStart = my_malloc( sizeof(int) * nMols );
   molEnd   = my_malloc( sizeof(int) * nMols );

   elecDip = my_malloc( sizeof(double *) * nMols );
   totEDip = my_malloc( sizeof(double *) * nMols );
   veloDip = my_malloc( sizeof(double *) * nMols );
   totVDip = my_malloc( sizeof(double *) * nMols );
   magnDip = my_malloc( sizeof(double *) * nMols );
   totMDip = my_malloc( sizeof(double *) * nMols );
   for ( mol=0; mol<nMols; mol++ ) 
   {
      elecDip[mol] = my_malloc( sizeof(double) * 3 );
      totEDip[mol] = my_malloc( sizeof(double) * 3 );
      veloDip[mol] = my_malloc( sizeof(double) * 3 );
      totVDip[mol] = my_malloc( sizeof(double) * 3 );
      magnDip[mol] = my_malloc( sizeof(double) * 3 );
      totMDip[mol] = my_malloc( sizeof(double) * 3 );
   }

   x = my_malloc( sizeof(double *) * nAtoms );
   for ( i=0; i<nAtoms; i++ ) 
   {
      x[i] = my_malloc( sizeof(double) * 3 );
   }

   q = my_malloc( sizeof(double) * nAtoms );

   /* read molecule information */
   /* use atID to count the index of each atom */

   atID = 0;
   for (mol=0; mol<nMols; mol++)
   {
      /* read first line for each molecule block: */
      /* number of atoms per molecule and excitation energy (a.u.) */

      int apm;
      fgets( line, sizeof( line ), file_inp );
      sscanf( line, "%d%lf", &apm, &eneExc[mol] );

      fgets( line, sizeof( line ), file_inp );
      sscanf( line, "%lf%lf%lf", &elecDip[mol][0], 
                                 &elecDip[mol][1], 
                                 &elecDip[mol][2] );

      fgets( line, sizeof( line ), file_inp );
      sscanf( line, "%lf%lf%lf", &veloDip[mol][0], 
                                 &veloDip[mol][1], 
                                 &veloDip[mol][2] );

      fgets( line, sizeof( line ), file_inp );
      sscanf( line, "%lf%lf%lf", &magnDip[mol][0], 
                                 &magnDip[mol][1], 
                                 &magnDip[mol][2] );

      /* get start and end index of each molecule */

      molStart[mol] = atID;
      molEnd[mol] = atID + apm;

      /* coordinates, mass and charge of atoms in each molecule */

      for (i=0; i<apm; i++)
      {
         fgets( line, sizeof( line ), file_inp );
         sscanf( line, "%lf%lf%lf%lf", &x[atID][0],&x[atID][1],&x[atID][2],
                                       &q[atID]);
         atID++;
      }

   }

   /* close input file */

   fclose( file_inp );

   /* check consistency: number of atoms */

   if ( atID != nAtoms )
   {
      my_prog_die("Error: incorrect number of atoms!\n");
   }

   /* calculate Hamiltonian */

   /* diagonal elements: excitation energy of monomers */

   for (iMol=0; iMol<nMols; iMol++)
   {
      excCoup[iMol][iMol] = eneExc[iMol];
   }

   /* off-diagonal elements: Coulomb coupling between TrESP charges */

   for (iMol=0; iMol<nMols; iMol++)
   for (jMol=iMol+1; jMol<nMols; jMol++)
   {
      /* check the first atom to see if iMol and jMol are the same molecule */

      if ( fabs(x[molStart[iMol]][0] - x[molStart[jMol]][0]) +
           fabs(x[molStart[iMol]][1] - x[molStart[jMol]][1]) +
           fabs(x[molStart[iMol]][2] - x[molStart[jMol]][2]) < 1.0e-10 )
      {
         excCoup[iMol][jMol] = 0.0;
         excCoup[jMol][iMol] = 0.0;
      }

      else
      {
         double coulomb = 0.0;
         for (i=molStart[iMol]; i<molEnd[iMol]; i++)
         for (j=molStart[jMol]; j<molEnd[jMol]; j++)
         {
             double xij, yij, zij, rij;

             xij = x[i][0]-x[j][0];
             yij = x[i][1]-x[j][1];
             zij = x[i][2]-x[j][2];

             if      (xij < -0.5 * box) { xij += box; }
             else if (xij >  0.5 * box) { xij -= box; }
             if      (yij < -0.5 * box) { yij += box; }
             else if (yij >  0.5 * box) { yij -= box; }
             if      (zij < -0.5 * box) { zij += box; }
             else if (zij >  0.5 * box) { zij -= box; }

             rij = sqrt(xij*xij + yij*yij + zij*zij);

            coulomb += q[i]*q[j]/rij;
         }
         excCoup[iMol][jMol] = coulomb;
         excCoup[jMol][iMol] = coulomb;
      }
   }

   /* use GSL to diagnolize the Hamiltonian */

   {
      int DIM = nMols;

      /* save 2D array 'excCoup' as GSL matrix: 'data' */

      gsl_matrix *data = gsl_matrix_alloc (DIM, DIM);
      for (iMol=0; iMol<DIM; iMol++)
      for (jMol=0; jMol<DIM; jMol++)
      {
         gsl_matrix_set (data, iMol, jMol, excCoup[iMol][jMol]);
      }

      /* allocate eval and evec: eigenvalues and eigenvectors */

      gsl_vector *eval = gsl_vector_alloc (DIM);
      gsl_matrix *evec = gsl_matrix_alloc (DIM, DIM);

      /* calculate eigenvalues and eigenvectors */
 
      my_eigen_symmv(data, DIM, eval, evec);

      /* compute transition electric and magnetic dipole moment */

      for (i=0; i<DIM; i++)
      {
         for (j=0; j<3; j++)
         {
            totEDip[i][j] = 0.0;
            totVDip[i][j] = 0.0;
            totMDip[i][j] = 0.0;

            int k;
            for (k=0; k<DIM; k++)
            {
               totEDip[i][j] += gsl_matrix_get(evec, k, i) * elecDip[k][j];
               totVDip[i][j] += gsl_matrix_get(evec, k, i) * veloDip[k][j];
               totMDip[i][j] += gsl_matrix_get(evec, k, i) * magnDip[k][j];
            }
         }
      }

      /* get total excitation energy, oscillator strength and 
         rotational length */

      double *totEneExc, *totDipStr, *totRotLen, *totRotVel;
      totEneExc = my_malloc(sizeof(double) * DIM);
      totDipStr = my_malloc(sizeof(double) * DIM);
      totRotLen = my_malloc(sizeof(double) * DIM);
      totRotVel = my_malloc(sizeof(double) * DIM);

      for (i=0; i<DIM; i++)
      {
         totEneExc[i] = gsl_vector_get(eval, i);

         totDipStr[i] = totEDip[i][0]*totEDip[i][0] + 
                        totEDip[i][1]*totEDip[i][1] + 
                        totEDip[i][2]*totEDip[i][2];

         totRotLen[i] = ( totEDip[i][0] * totMDip[i][0] + 
                          totEDip[i][1] * totMDip[i][1] + 
                          totEDip[i][2] * totMDip[i][2] ) * (-1.0);

         totRotVel[i] = ( totVDip[i][0] * totMDip[i][0] + 
                          totVDip[i][1] * totMDip[i][1] + 
                          totVDip[i][2] * totMDip[i][2] ) / totEneExc[i];
      }

      FILE *file_txt, *file_dat;
      file_txt = fopen( "snap.txt", "w" ) ;
      file_dat = fopen( "snap.dat", "w" ) ;

      if (file_txt == NULL || file_dat == NULL )
      {
         my_prog_die( "Error: Cannot write to file snap.txt snap.dat snap.mat !\n") ;
      }

      int wavLen5;
      for (wavLen5=100*5; wavLen5<1000*5; wavLen5++)
      {
         double wavLen = (double)wavLen5 / 5.0;

         double absDipStr = 0.0;
         double ecdRotLen = 0.0;
         double ecdRotVel = 0.0;

         for (i=0; i<DIM; i++)
         {
            /* del in eV, GNorm in eV^-1, preFact dimensionless */

            double del    = eVxnm / wavLen - Eh2eV * totEneExc[i];
            double GNorm  = 1.0 / (fwhm * sqrt(PI));
            double preFac = GNorm * exp(-1.0 * del * del / (fwhm * fwhm)) *
                            (eVxnm / wavLen);

            absDipStr += totDipStr[i] * preFac * au2eps; /* M^-1 cm^-1 */
            ecdRotLen += totRotLen[i] * preFac * au2theta; /* deg cm^2 dmol^-1 */
            ecdRotVel += totRotVel[i] * preFac * au2theta; /* deg cm^2 dmol^-1 */
         }

          double asym_g_len = 0.0;
          double asym_g_vel = 0.0;
          if (absDipStr > 1.0e-12) 
          { 
             asym_g_len = ecdRotLen * theta2eps / absDipStr; 
             asym_g_vel = ecdRotVel * theta2eps / absDipStr; 
          }

         fprintf (file_txt, "%10.2f%20.10e%20.10e%20.10e\n", 
                  wavLen, ecdRotVel, ecdRotLen, absDipStr);
      }

      fprintf (file_dat, "#Wav.Len.  Dip.Str.  Rot.Vel.  Rot.Len.\n");
      for (i=0; i<DIM; i++)
      {
         fprintf (file_dat, "%20.6f%20.6f%20.6f%20.6f\n", 
                  Ehxnm / totEneExc[i], totDipStr[i], totRotVel[i], totRotLen[i]);
      }

      fclose(file_txt);
      fclose(file_dat);

      gsl_vector_free (eval);
      gsl_matrix_free (evec);
      gsl_matrix_free (data);
   }

   /* free the arrays */

   free (molStart);
   free (molEnd);

   /* end of main program */

   return 0;
}
