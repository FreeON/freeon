/*----------------------------------------------------------------------------------
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copywrite 2000, The University of California
!
------------------------------------------------------------------------------------*/
/* C wrapper for the PHiPAC DGEMM routine mm_double_NN_1.  
   Computes C(M,N)=A(M,K).B(K,M)+beta*C(M,N) */

extern void mm_double_NN_1(int,int,int,double*,double*,double*,int,int,int,double);
void dgemm_nn_(int* M, int* K,int* N,double* beta,double* A,double* B,double* C)
{mm_double_NN_1(*N,*K,*M,B,A,C,*K,*M,*M,*beta);}
void dgemm_nn__(int* M, int* K,int* N,double* beta,double* A,double* B,double* C)
{mm_double_NN_1(*N,*K,*M,B,A,C,*K,*M,*M,*beta);}

/* C wrapper for the PHiPAC DGEMM routine mm_double_NN_c.  
   Computes C(M,N)=alpha*A(M,K).B(K,M)+beta*C(M,N) */

extern void mm_double_NN_c(int,int,int,double*,double*,double*,int,int,int,double,double);
void dgemm_nnc_(int* M, int* K,int* N,double* alpha,double* beta,double* A,double* B,double* C)
{mm_double_NN_c(*N,*K,*M,B,A,C,*K,*M,*M,*alpha,*beta);}
void dgemm_nnc__(int* M, int* K,int* N,double* alpha,double* beta,double* A,double* B,double* C)
{mm_double_NN_c(*N,*K,*M,B,A,C,*K,*M,*M,*alpha,*beta);}

/*---------------------------------------------------------------------------------
   C wrapper for the PHiPAC DGEMM routine mm_double_TN_1.  
   Computes C(M,N)=AT(M,K).B(K,M)+beta*C(M,N) */

extern void mm_double_NT_1(int,int,int,double*,double*,double*,int,int,int,double);
void dgemm_tn_(int* M, int* K,int* N,double* beta,double* A,double* B,double* C)
{mm_double_NT_1(*N,*K,*M,B,A,C,*K,*K,*M,*beta);}
void dgemm_tn__(int* M, int* K,int* N,double* beta,double* A,double* B,double* C)
{mm_double_NT_1(*N,*K,*M,B,A,C,*K,*K,*M,*beta);}


/* C wrapper for the PHiPAC DGEMM routine mm_double_TN_1.  
   Computes C(M,N)=alpha*AT(M,K).B(K,M)+beta*C(M,N) */

extern void mm_double_NT_c(int,int,int,double*,double*,double*,int,int,int,double,double);
void dgemm_tnc_(int* M, int* K,int* N,double* alpha,double* beta,double* A,double* B,double* C)
{mm_double_NT_c(*N,*K,*M,B,A,C,*K,*K,*M,*alpha,*beta);}
void dgemm_tnc__(int* M, int* K,int* N,double* alpha,double* beta,double* A,double* B,double* C)
{mm_double_NT_c(*N,*K,*M,B,A,C,*K,*K,*M,*alpha,*beta);}

/*---------------------------------------------------------------------------------
   C wrapper for the PHiPAC DGEMM routine dgemm_NT
   Computes C(M,N)=A(M,K).BT(K,M)+beta*C(M,N) */

extern void mm_double_TN_1(int,int,int,double*,double*,double*,int,int,int,double);
void dgemm_nt_(int* M, int* K,int* N,double* beta,double* A,double* B,double* C)
{mm_double_TN_1(*N,*K,*M,B,A,C,*M,*N,*M,*beta);}
void dgemm_nt__(int* M, int* K,int* N,double* beta,double* A,double* B,double* C)
{mm_double_TN_1(*N,*K,*M,B,A,C,*M,*N,*M,*beta);}


/* C wrapper for the PHiPAC DGEMM routine dgemm_NTc
   Computes C(M,N)=alpha*A(M,K).BT(K,M)+beta*C(M,N) */

extern void mm_double_TN_c(int,int,int,double*,double*,double*,int,int,int,double,double);
void dgemm_ntc_(int* M, int* K,int* N,double* alpha,double* beta,double* A,double* B,double* C)
{mm_double_TN_c(*N,*K,*M,B,A,C,*M,*N,*M,*alpha,*beta);}
void dgemm_ntc__(int* M, int* K,int* N,double* alpha,double* beta,double* A,double* B,double* C)
{mm_double_TN_c(*N,*K,*M,B,A,C,*M,*N,*M,*alpha,*beta);}
