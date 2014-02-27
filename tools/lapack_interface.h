/** @file
 *
 * The lapack interface.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __LAPACK_INTERFACE_H
#define __LAPACK_INTERFACE_H

#include "config.h"

/** The dgemm() interface. We are assuming that the dgemm function is using a
 * Fortran style interface, hence all the pointers in the argument list.
 *
 * @param[in] TRANSA
 *           TRANSA is CHARACTER*1
 *            On entry, TRANSA specifies the form of op( A ) to be used in
 *            the matrix multiplication as follows:
 *
 *               TRANSA = 'N' or 'n',  op( A ) = A.
 *
 *               TRANSA = 'T' or 't',  op( A ) = A**T.
 *
 *               TRANSA = 'C' or 'c',  op( A ) = A**T.
 * @param[in]	TRANSB
 *           TRANSB is CHARACTER*1
 *            On entry, TRANSB specifies the form of op( B ) to be used in
 *            the matrix multiplication as follows:
 *
 *               TRANSB = 'N' or 'n',  op( B ) = B.
 *
 *               TRANSB = 'T' or 't',  op( B ) = B**T.
 *
 *               TRANSB = 'C' or 'c',  op( B ) = B**T.
 * @param[in]	M
 *           M is INTEGER
 *            On entry,  M  specifies  the number  of rows  of the  matrix
 *            op( A )  and of the  matrix  C.  M  must  be at least  zero.
 * @param[in]	N
 *           N is INTEGER
 *            On entry,  N  specifies the number  of columns of the matrix
 *            op( B ) and the number of columns of the matrix C. N must be
 *            at least zero.
 * @param[in]	K
 *           K is INTEGER
 *            On entry,  K  specifies  the number of columns of the matrix
 *            op( A ) and the number of rows of the matrix op( B ). K must
 *            be at least  zero.
 * @param[in]	ALPHA
 *           ALPHA is DOUBLE PRECISION.
 *            On entry, ALPHA specifies the scalar alpha.
 * @param[in]	A
 *           A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
 *            k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
 *            Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
 *            part of the array  A  must contain the matrix  A,  otherwise
 *            the leading  k by m  part of the array  A  must contain  the
 *            matrix A.
 * @param[in]	LDA
 *           LDA is INTEGER
 *            On entry, LDA specifies the first dimension of A as declared
 *            in the calling (sub) program. When  TRANSA = 'N' or 'n' then
 *            LDA must be at least  max( 1, m ), otherwise  LDA must be at
 *            least  max( 1, k ).
 * @param[in]	B
 *           B is DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
 *            n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
 *            Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
 *            part of the array  B  must contain the matrix  B,  otherwise
 *            the leading  n by k  part of the array  B  must contain  the
 *            matrix B.
 * @param[in]	LDB
 *           LDB is INTEGER
 *            On entry, LDB specifies the first dimension of B as declared
 *            in the calling (sub) program. When  TRANSB = 'N' or 'n' then
 *            LDB must be at least  max( 1, k ), otherwise  LDB must be at
 *            least  max( 1, n ).
 * @param[in]	BETA
 *           BETA is DOUBLE PRECISION.
 *            On entry,  BETA  specifies the scalar  beta.  When  BETA  is
 *            supplied as zero then C need not be set on input.
 * @param[in,out]	C
 *           C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).
 *            Before entry, the leading  m by n  part of the array  C must
 *            contain the matrix  C,  except when  beta  is zero, in which
 *            case C need not be set on entry.
 *            On exit, the array  C  is overwritten by the  m by n  matrix
 *            ( alpha*op( A )*op( B ) + beta*C ).
 * @param[in]	LDC
 *           LDC is INTEGER
 *            On entry, LDC specifies the first dimension of C as declared
 *            in  the  calling  (sub)  program.   LDC  must  be  at  least
 *            max( 1, m ).
 */
extern "C"
void F77_FUNC(dgemm, DGEMM) (const char *const TRANSA,
    const char *const TRANSB, const int *const M, const int *const N,
    const int *const K, const double *const ALPHA, const double *const A,
    const int *const LDA, const double *const B,
    const int *const LDB, const double *const BETA,
    double *const C, const int *const LDC);

/** DSYEV computes all eigenvalues and, optionally, eigenvectors of a
 *  real symmetric matrix A.
 *
 * @param [in]	JOBZ
 *           JOBZ is CHARACTER*1
 *           = 'N':  Compute eigenvalues only;
 *           = 'V':  Compute eigenvalues and eigenvectors.
 * @param [in]	UPLO
 *           UPLO is CHARACTER*1
 *           = 'U':  Upper triangle of A is stored;
 *           = 'L':  Lower triangle of A is stored.
 * @param [in]	N
 *           N is INTEGER
 *           The order of the matrix A.  N >= 0.
 * @param [in,out]	A
 *           A is DOUBLE PRECISION array, dimension (LDA, N)
 *           On entry, the symmetric matrix A.  If UPLO = 'U', the
 *           leading N-by-N upper triangular part of A contains the
 *           upper triangular part of the matrix A.  If UPLO = 'L',
 *           the leading N-by-N lower triangular part of A contains
 *           the lower triangular part of the matrix A.
 *           On exit, if JOBZ = 'V', then if INFO = 0, A contains the
 *           orthonormal eigenvectors of the matrix A.
 *           If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
 *           or the upper triangle (if UPLO='U') of A, including the
 *           diagonal, is destroyed.
 * @param [in]	LDA
 *           LDA is INTEGER
 *           The leading dimension of the array A.  LDA >= max(1,N).
 * @param [out]	W
 *           W is DOUBLE PRECISION array, dimension (N)
 *           If INFO = 0, the eigenvalues in ascending order.
 * @param [out]	WORK
 *           WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
 *           On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
 * @param [in]	LWORK
 *           LWORK is INTEGER
 *           The length of the array WORK.  LWORK >= max(1,3*N-1).
 *           For optimal efficiency, LWORK >= (NB+2)*N,
 *           where NB is the blocksize for DSYTRD returned by ILAENV.
 *           If LWORK = -1, then a workspace query is assumed; the routine
 *           only calculates the optimal size of the WORK array, returns
 *           this value as the first entry of the WORK array, and no error
 *           message related to LWORK is issued by XERBLA.
 * @param [out]	INFO
 *           INFO is INTEGER
 *           = 0:  successful exit
 *           < 0:  if INFO = -i, the i-th argument had an illegal value
 *           > 0:  if INFO = i, the algorithm failed to converge; i
 *                 off-diagonal elements of an intermediate tridiagonal
 *                 form did not converge to zero.
 */
extern "C"
void F77_FUNC(dsyev, DSYEV) (const char *const JOBZ, const char *const UPLO,
    const int *const N, double *const A, const int *const LDA,
    double *const W, double *const WORK, const int *const LWORK,
    int *const INFO);

#endif
