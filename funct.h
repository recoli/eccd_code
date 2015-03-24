/*****************************************************************************
 *  This file is part of the eccd_code program.                              *  
 *  Copyright (C) 2015 Xin Li <lixin.reco@gmail.com>                         *  
 *                                                                           *
 *  Filename:  funct.h                                                       *
 *  Function:  functions for memory handeling and matri solving              *
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

void* my_malloc(size_t bytes);

void my_eigen_symmv(gsl_matrix* data, int DIM,
                    gsl_vector* eval, gsl_matrix* evec);

void my_prog_die(char* txt);

void my_norm_vec(gsl_vector* x, int DIM);

void my_outer_prod(gsl_vector* x, gsl_vector* y, int DIM,
                   gsl_vector* z);
