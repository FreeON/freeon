/* Provide the dgemm interface.
 *
 *  Arguments
 *  ==========
 *
 *  TRANSA - CHARACTER*1.
 *           On entry, TRANSA specifies the form of op( A ) to be used in
 *           the matrix multiplication as follows:
 *
 *              TRANSA = 'N' or 'n',  op( A ) = A.
 *
 *              TRANSA = 'T' or 't',  op( A ) = A'.
 *
 *              TRANSA = 'C' or 'c',  op( A ) = A'.
 *
 *           Unchanged on exit.
 *
 *  TRANSB - CHARACTER*1.
 *           On entry, TRANSB specifies the form of op( B ) to be used in
 *           the matrix multiplication as follows:
 *
 *              TRANSB = 'N' or 'n',  op( B ) = B.
 *
 *              TRANSB = 'T' or 't',  op( B ) = B'.
 *
 *              TRANSB = 'C' or 'c',  op( B ) = B'.
 *
 *           Unchanged on exit.
 *
 *  M      - INTEGER.
 *           On entry,  M  specifies  the number  of rows  of the  matrix
 *           op( A )  and of the  matrix  C.  M  must  be at least  zero.
 *           Unchanged on exit.
 *
 *  N      - INTEGER.
 *           On entry,  N  specifies the number  of columns of the matrix
 *           op( B ) and the number of columns of the matrix C. N must be
 *           at least zero.
 *           Unchanged on exit.
 *
 *  K      - INTEGER.
 *           On entry,  K  specifies  the number of columns of the matrix
 *           op( A ) and the number of rows of the matrix op( B ). K must
 *           be at least  zero.
 *           Unchanged on exit.
 *
 *  ALPHA  - DOUBLE PRECISION.
 *           On entry, ALPHA specifies the scalar alpha.
 *           Unchanged on exit.
 *
 *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
 *           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
 *           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
 *           part of the array  A  must contain the matrix  A,  otherwise
 *           the leading  k by m  part of the array  A  must contain  the
 *           matrix A.
 *           Unchanged on exit.
 *
 *  LDA    - INTEGER.
 *           On entry, LDA specifies the first dimension of A as declared
 *           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
 *           LDA must be at least  max( 1, m ), otherwise  LDA must be at
 *           least  max( 1, k ).
 *           Unchanged on exit.
 *
 *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
 *           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
 *           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
 *           part of the array  B  must contain the matrix  B,  otherwise
 *           the leading  n by k  part of the array  B  must contain  the
 *           matrix B.
 *           Unchanged on exit.
 *
 *  LDB    - INTEGER.
 *           On entry, LDB specifies the first dimension of B as declared
 *           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
 *           LDB must be at least  max( 1, k ), otherwise  LDB must be at
 *           least  max( 1, n ).
 *           Unchanged on exit.
 *
 *  BETA   - DOUBLE PRECISION.
 *           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
 *           supplied as zero then C need not be set on input.
 *           Unchanged on exit.
 *
 *  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
 *           Before entry, the leading  m by n  part of the array  C must
 *           contain the matrix  C,  except when  beta  is zero, in which
 *           case C need not be set on entry.
 *           On exit, the array  C  is overwritten by the  m by n  matrix
 *           ( alpha*op( A )*op( B ) + beta*C ).
 *
 *  LDC    - INTEGER.
 *           On entry, LDC specifies the first dimension of C as declared
 *           in  the  calling  (sub)  program.   LDC  must  be  at  least
 *           max( 1, m ).
 *           Unchanged on exit.
 */

#include "config.h"
#include "lal.h"

void
lal_dgemm (const char *transA, const char *transB, const int M, const int N,
    const int K, const double alpha, const matrix_t A, const int lda,
    const matrix_t B, const int ldb, const double beta, matrix_t C,
    const int ldc)
{
}
