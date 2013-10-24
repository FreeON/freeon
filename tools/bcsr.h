/** @file
 *
 * The header file for the BCSR class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __BCSR_H
#define __BCSR_H

/** The BCSR Matrix. */
class BCSR
{
  private:

    /** The number of rows. */
    int M;

    /** The number of columns. */
    int N;

    /** The number of spin matrices. */
    int NSMat;

    /** The number of atoms. */
    int NAtoms;

    /** The number of non-zero elements. */
    int numberNonZero;

    /** The number of blocks. */
    int numberBlocks;

    /** The number of basis functions. */
    int numberBasisFunctions;

    /** The block sizes per atom. */
    int *blockSize;

    /** The offsets. */
    int *offset;

    /** The array of row indices. */
    int *rowPointer;

    /** The array of column indices. */
    int *columnPointer;

    /** The array of block indices. */
    int *blockPointer;

    /** The non-zero matrix elements. */
    double *matrix;

  public:

    BCSR (char *filename);
    ~BCSR (void);
    void getSpectralBounds (int method, double *minBound, double *maxBound);
    int getNumberNonZero (void);
    double getElement (int i);
    void toDense (int *M, int *N, double **ADense);
    void toStr (void);
    void put (char *filename);
    void toMM (char *filename);
};

#endif
