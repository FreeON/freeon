/*
     This code is part of the MondoSCF suite of programs for linear scaling
     electronic structure theory and ab initio molecular dynamics.

     Copyright (2004). The Regents of the University of California. This
     material was produced under U.S. Government contract W-7405-ENG-36
     for Los Alamos National Laboratory, which is operated by the University
     of California for the U.S. Department of Energy. The U.S. Government has
     rights to use, reproduce, and distribute this software.  NEITHER THE
     GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
     OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.

     This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by the
     Free Software Foundation; either version 2 of the License, or (at your
     option) any later version. Accordingly, this program is distributed in
     the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE. See the GNU General Public License at www.gnu.org for details.

     While you may do as you like with this software, the GNU license requires
     that you clearly mark derivative software.  In addition, you are encouraged
     to return derivative works to the MondoSCF group for review, and possible
     disemination in future releases.
*/
/*------------------------------------------------------------------*/
/*    C wrappers for PHiPAC DGEMMs                                  */
/*    Author: Matt Challacombe                                      */
/*------------------------------------------------------------------*/

#include "config.h"

#if defined (PHIPAC)
/* C wrapper for the PHiPAC DGEMM routine mm_double_NN_1.  
   Computes C(M, N)=A(M, K).B(K, M)+beta*C(M, N) 
*/
extern void mm_double_NN_1(int, int, int, double*, double*, double*, int, int, int, double);
void dgemm_nn_(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{mm_double_NN_1(*N, *K, *M, B, A, C, *K, *M, *M, *beta);}
void dgemm_nn__(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{mm_double_NN_1(*N, *K, *M, B, A, C, *K, *M, *M, *beta);}
void dgemm_nn(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{mm_double_NN_1(*N, *K, *M, B, A, C, *K, *M, *M, *beta);}

/* C wrapper for the PHiPAC DGEMM routine mm_double_NN_c.  
   Computes C(M, N)=alpha*A(M, K).B(K, M)+beta*C(M, N) */

extern void mm_double_NN_c(int, int, int, double*, double*, double*, int, int, int, double, double);
void dgemm_nnc_(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{mm_double_NN_c(*N, *K, *M, B, A, C, *K, *M, *M, *alpha, *beta);}
void dgemm_nnc__(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{mm_double_NN_c(*N, *K, *M, B, A, C, *K, *M, *M, *alpha, *beta);}
void dgemm_nnc(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{mm_double_NN_c(*N, *K, *M, B, A, C, *K, *M, *M, *alpha, *beta);}

/*---------------------------------------------------------------------------------
   C wrapper for the PHiPAC DGEMM routine mm_double_TN_1.  
   Computes C(M, N)=AT(M, K).B(K, M)+beta*C(M, N) */

extern void mm_double_NT_1(int, int, int, double*, double*, double*, int, int, int, double);
void dgemm_tn_(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{mm_double_NT_1(*N, *K, *M, B, A, C, *K, *K, *M, *beta);}
void dgemm_tn__(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{mm_double_NT_1(*N, *K, *M, B, A, C, *K, *K, *M, *beta);}
void dgemm_tn(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{mm_double_NT_1(*N, *K, *M, B, A, C, *K, *K, *M, *beta);}


/* C wrapper for the PHiPAC DGEMM routine mm_double_TN_1.  
   Computes C(M, N)=alpha*AT(M, K).B(K, M)+beta*C(M, N) */

extern void mm_double_NT_c(int, int, int, double*, double*, double*, int, int, int, double, double);
void dgemm_tnc_(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{mm_double_NT_c(*N, *K, *M, B, A, C, *K, *K, *M, *alpha, *beta);}
void dgemm_tnc__(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{mm_double_NT_c(*N, *K, *M, B, A, C, *K, *K, *M, *alpha, *beta);}
void dgemm_tnc(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{mm_double_NT_c(*N, *K, *M, B, A, C, *K, *K, *M, *alpha, *beta);}

/*---------------------------------------------------------------------------------
   C wrapper for the PHiPAC DGEMM routine dgemm_NT
   Computes C(M, N)=A(M, K).BT(K, M)+beta*C(M, N) */

extern void mm_double_TN_1(int, int, int, double*, double*, double*, int, int, int, double);
void dgemm_nt_(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{mm_double_TN_1(*N, *K, *M, B, A, C, *M, *N, *M, *beta);}
void dgemm_nt__(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{mm_double_TN_1(*N, *K, *M, B, A, C, *M, *N, *M, *beta);}
void dgemm_nt(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{mm_double_TN_1(*N, *K, *M, B, A, C, *M, *N, *M, *beta);}


/* C wrapper for the PHiPAC DGEMM routine dgemm_NTc
   Computes C(M, N)=alpha*A(M, K).BT(K, M)+beta*C(M, N) */

extern void mm_double_TN_c(int, int, int, double*, double*, double*, int, int, int, double, double);
void dgemm_ntc_(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{mm_double_TN_c(*N, *K, *M, B, A, C, *M, *N, *M, *alpha, *beta);}
void dgemm_ntc__(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{mm_double_TN_c(*N, *K, *M, B, A, C, *M, *N, *M, *alpha, *beta);}
void dgemm_ntc(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{mm_double_TN_c(*N, *K, *M, B, A, C, *M, *N, *M, *alpha, *beta);}

#else

extern void dgemm_(char *, char *, int *, int *, int *, double *, double *,
  int *, double *, int *, double *, double *, int *);

/* Define a factor of one so that we can use "& one" later.
 */
double one = 1.0;

/* C wrapper for dgemm NN_1
 */
void dgemm_nn_(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{dgemm_("N", "N", M, N, K, & one, A, M, B, K, beta, C, M);}
void dgemm_nn__(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{dgemm_("N", "N", M, N, K, & one, A, M, B, K, beta, C, M);}
void dgemm_nn(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{dgemm_("N", "N", M, N, K, & one, A, M, B, K, beta, C, M);}

/* C wrapper for dgemm_NN_c
 */
void dgemm_nnc_(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{dgemm_("N", "N", M, N, K, alpha, A, M, B, K, beta, C, M);}
void dgemm_nnc__(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{dgemm_("N", "N", M, N, K, alpha, A, M, B, K, beta, C, M);}
void dgemm_nnc(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{dgemm_("N", "N", M, N, K, alpha, A, M, B, K, beta, C, M);}

/* C wrapper for dgemm_TN_1
 */
void dgemm_tn_(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{dgemm_("T", "N", M, N, K, & one, A, M, B, N, beta, C, M);}
void dgemm_tn__(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{dgemm_("T", "N", M, N, K, & one, A, M, B, N, beta, C, M);}
void dgemm_tn(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{dgemm_("T", "N", M, N, K, & one, A, M, B, N, beta, C, M);}


/* C wrapper for dgemm_TN_c
 */
void dgemm_tnc_(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{dgemm_("T", "N", M, N, K, alpha, A, M, B, N, beta, C, M);}
void dgemm_tnc__(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{dgemm_("T", "N", M, N, K, alpha, A, M, B, N, beta, C, M);}
void dgemm_tnc(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{dgemm_("T", "N", M, N, K, alpha, A, M, B, N, beta, C, M);}

/* C wrapper for dgemm_NT_1
 */
void dgemm_nt_(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{dgemm_("N", "T", M, N, K, & one, A, M, B, N, beta, C, M);}
void dgemm_nt__(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{dgemm_("N", "T", M, N, K, & one, A, M, B, N, beta, C, M);}
void dgemm_nt(int* M, int* K, int* N, double* beta, double* A, double* B, double* C)
{dgemm_("N", "T", M, N, K, & one, A, M, B, N, beta, C, M);}


/* C wrapper for dgemm_NT_c
 */
void dgemm_ntc_(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{dgemm_("N", "T", M, N, K, alpha, A, M, B, N, beta, C, M);}
void dgemm_ntc__(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{dgemm_("N", "T", M, N, K, alpha, A, M, B, N, beta, C, M);}
void dgemm_ntc(int* M, int* K, int* N, double* alpha, double* beta, double* A, double* B, double* C)
{dgemm_("N", "T", M, N, K, alpha, A, M, B, N, beta, C, M);}

#endif
