/** @file
 *
 * Some macros for indexing.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __INDEX_H
#define __INDEX_H

/** The linear offset inside a dense square matrix block. Column-major order. */
#define BLOCK_INDEX(i, j, iLower, jLower, blocksize) ((i-iLower)+(j-jLower)*blocksize)

/** The linear offset inside a dense non-square matrix block. Column-major
 * order. */
#define BLOCK_INDEX_NONSQUARE(i, j, iLower, jLower, M, N) ((i-iLower)+(j-jLower)*M)

/** The linear offset into a 3D array. */
#define BLOCK_INDEX_3(i, j, k, N) (i+j*N+k*N*N)

/** The linear tree index. */
#define CHILD_INDEX(i, j) ((i << 1) | j)

#endif
