/** @file
 *
 * The implementation of the BCSR class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "bcsr.h"
#include "index.h"
#include "lapack_interface.h"
#include "logger.h"

#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

/** The constructor.
 *
 * @param filename The filename of the BCSR file.
 */
BCSR::BCSR (char *filename)
{
  FILE *fd = NULL;

  if((fd = fopen(filename, "r")) == NULL)
  {
    ABORT("error opening BCSR file \"%s\"\n", filename);
  }

  int result;

  if((result = fread(&NSMat, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&NAtoms, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&numberNonZero, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&numberBlocks, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&numberBasisFunctions, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&M, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }
  if((result = fread(&N, sizeof(int), 1, fd)) != 1) { ABORT("read error (read %d items)\n", result); }

  blockSize = new int[NAtoms];
  offset = new int[NAtoms];

  if((result = fread(blockSize, sizeof(int), NAtoms, fd)) != NAtoms) { ABORT("read error\n"); }
  if((result = fread(offset, sizeof(int), NAtoms, fd)) != NAtoms) { ABORT("read error\n"); }

  rowPointer = new int[NAtoms+1];
  columnPointer = new int[numberBlocks];
  blockPointer = new int[numberBlocks];
  matrix = new double[numberNonZero];

  double *dummy = new double[numberBlocks];

  if((result = fread(rowPointer, sizeof(int), NAtoms+1, fd)) != NAtoms+1) { ABORT("read error\n"); }
  if((result = fread(columnPointer, sizeof(int), numberBlocks, fd)) != numberBlocks) { ABORT("read error\n"); }
  if((result = fread(blockPointer, sizeof(int), numberBlocks, fd)) != numberBlocks) { ABORT("read error\n"); }
  if((result = fread(dummy, sizeof(double), numberBlocks, fd)) != numberBlocks) { ABORT("read error\n"); }
  if((result = fread(matrix, sizeof(double), numberNonZero, fd)) != numberNonZero) { ABORT("read error\n"); }

  delete[] dummy;

  /* Fortran counts starting at 1, C starts from 0. */
  for(int i = 0; i < NAtoms; i++)
  {
    offset[i]--;
  }
  for(int i = 0; i < NAtoms+1; i++)
  {
    rowPointer[i]--;
  }
  for(int i = 0; i < numberBlocks; i++)
  {
    columnPointer[i]--;
    blockPointer[i]--;
  }

  if((result = fclose(fd)) != 0)
  {
    ABORT("error closing BCSR file\n");
  }

  INFO("read BCSR matrix, size %dx%d, %d nonzeros, %1.4f percent nonzero elements\n",
      M, N, numberNonZero, 100*numberNonZero/(double) (N*N));
}

/** The destructor.
 */
BCSR::~BCSR (void)
{
  delete[] blockSize;
  delete[] offset;
  delete[] rowPointer;
  delete[] columnPointer;
  delete[] blockPointer;
  delete[] matrix;
}

/** Get the spectral bounds of the matrix by using the Gershgorin circle
 * theorem.
 *
 * Estimate spectral bounds via Gersgorin approximation, @f$ \left[
 * F_{min}-F_{max} \right] @f$.
 *
 * In detail:
 * @f[
 *   R_{i} = \sum_{j \neq i} \left| F_{ij} \right|
 * @f]
 * @f[
 *   F_{\mathrm{max}} = \max_{i} \left\{ F_{ii} + R_{i} \right\}
 * @f]
 * @f[
 *   F_{\mathrm{min}} = \min_{i} \left\{ F_{ii} - R_{i} \right\}
 * @f]
 *
 * @param method The method to use. method = 0 is Gershgorin; method = 1 is
 * full eigensolve.
 * @param minBound [out] The lower bound.
 * @param maxBound [out] The upper bound.
 */
void BCSR::getSpectralBounds (int method, double *minBound, double *maxBound)
{
  switch(method)
  {
    case 0:
      {
        assert(M == N);

        double *R = new double[M];
        double *diagonal = new double[M];

        memset(R, 0, sizeof(double)*M);

        for(int atom = 0; atom < NAtoms; atom++)
        {
          int MBlock = blockSize[atom];
          int rowOffset = offset[atom];

          for(int iRow = rowPointer[atom]; iRow < rowPointer[atom+1]; iRow++)
          {
            int NBlock = blockSize[columnPointer[iRow]];
            int columnOffset = offset[columnPointer[iRow]];

            for(int iSMat = 0; iSMat < NSMat; iSMat++)
            {
              int pointer = blockPointer[iRow]+iSMat*MBlock*NBlock;
              int i_block;
              int j_block;

              switch(iSMat+1)
              {
                case 1:
                  i_block = rowOffset;
                  j_block = columnOffset;
                  break;

                case 2:
                  i_block = rowOffset;
                  j_block = columnOffset+numberBasisFunctions;
                  break;

                case 3:
                  i_block = rowOffset+numberBasisFunctions;
                  j_block = columnOffset;
                  break;

                case 4:
                  i_block = rowOffset+numberBasisFunctions;
                  j_block = columnOffset+numberBasisFunctions;
                  break;

                default:
                  ABORT("error (NSMat = %d)\n", NSMat);
                  break;
              }

              for(int i = 0; i < MBlock; i++) {
                for(int j = 0; j < NBlock; j++)
                {
                  int index = pointer+i+j*MBlock;
                  if(i+i_block == j+j_block)
                  {
                    diagonal[i+i_block] = matrix[index];
                  }

                  else
                  {
                    R[i+i_block] += fabs(matrix[index]);
                  }
                }
              }
            }
          }
        }

        *minBound = diagonal[0];
        *maxBound = diagonal[0];

        for(int i = 0; i < M; i++)
        {
          if(*minBound > diagonal[i]-R[i])
          {
            *minBound = diagonal[i]-R[i];
          }

          if(*maxBound < diagonal[i]+R[i])
          {
            *maxBound = diagonal[i]+R[i];
          }
        }

        delete[] diagonal;
        delete[] R;
      }

      break;

    case 1:
#ifdef DSYEV
      {
        int M, N;
        double *ADense;

        toDense(&M, &N, &ADense);
        assert(M == N);

        double *eigenvalue = new double[N];
        int lwork = 3*N;
        double *work = new double[lwork];
        int info;

        DSYEV((char*) "N", (char*) "U", &N, ADense, &N, eigenvalue, work, &lwork, &info);

        if(info != 0)
        {
          ABORT("error in dsyev\n");
        }
        delete[] work;

        *minBound = eigenvalue[0];
        *maxBound = eigenvalue[0];

        for(int i = 0; i < N; i++)
        {
          if(*minBound > eigenvalue[i])
          {
            *minBound = eigenvalue[i];
          }

          if(*maxBound < eigenvalue[i])
          {
            *maxBound = eigenvalue[i];
          }
        }

        delete[] eigenvalue;
      }
#else
      ABORT("full eigensolve was not built because lapack was missing at compile time\n");
#endif
      break;

    default:
      ABORT("unknow method\n");
      break;
  }
}

/** Return the number of non-zero elements.
 *
 * @return The number of non-zero elements.
 */
int BCSR::getNumberNonZero (void)
{
  return numberNonZero;
}

/** Get a particular non-zero element.
 *
 * @param i The index of the non-zero element.
 *
 * @return The value of that element.
 */
double BCSR::getElement (int i)
{
  assert(i >= 0);
  assert(i < numberNonZero);

  return matrix[i];
}

/** Convert a BCSR matrix into a dense matrix.
 *
 * @param M [out] The number of rows.
 * @param N [out] The number of columns.
 * @param ADense [out] The dense matrix.
 */
void BCSR::toDense (int *M, int *N, double **ADense)
{
  *M = this->M;
  *N = this->N;

  *ADense = new double[this->M*this->N];
  memset(*ADense, 0, sizeof(double)*this->M*this->N);

  for(int atom = 0; atom < NAtoms; atom++)
  {
    int MBlock = blockSize[atom];
    int rowOffset = offset[atom];

    for(int iRow = rowPointer[atom]; iRow < rowPointer[atom+1]; iRow++)
    {
      int NBlock = blockSize[columnPointer[iRow]];
      int columnOffset = offset[columnPointer[iRow]];

      for(int iSMat = 0; iSMat < NSMat; iSMat++)
      {
        int pointer = blockPointer[iRow]+iSMat*MBlock*NBlock;
        int i_block;
        int j_block;

        switch(iSMat+1)
        {
          case 1:
            i_block = rowOffset;
            j_block = columnOffset;
            break;

          case 2:
            i_block = rowOffset;
            j_block = columnOffset+numberBasisFunctions;
            break;

          case 3:
            i_block = rowOffset+numberBasisFunctions;
            j_block = columnOffset;
            break;

          case 4:
            i_block = rowOffset+numberBasisFunctions;
            j_block = columnOffset+numberBasisFunctions;
            break;

          default:
            ABORT("error\n");
            break;
        }

        for(int i = 0; i < MBlock; i++) {
          for(int j = 0; j < NBlock; j++)
          {
            (*ADense)[BLOCK_INDEX_NONSQUARE(i+i_block, j+j_block, 0, 0, this->M, this-N)] =
              matrix[pointer+i+j*MBlock];
          }
        }
      }
    }
  }
}

/** Print some information on the BCSR matrix.
 */
void BCSR::toStr (const bool verbose)
{
  printf("BCSR: M             = %d\n", M);
  printf("BCSR: N             = %d\n", N);
  printf("BCSR: NBasF         = %d\n", numberBasisFunctions);
  printf("BCSR: NSMat         = %d\n", NSMat);
  printf("BCSR: NAtoms        = %d\n", NAtoms);
  if(verbose)
  {
    printf("BCSR: block sizes   = {");
    for(int i = 0; i < NAtoms; i++)
    {
      printf(" %d", blockSize[i]);
    }
    printf(" }\n");
    printf("BCSR: offset        = {");
    for(int i = 0; i < NAtoms; i++)
    {
      printf(" %d", offset[i]);
    }
    printf(" }\n");
    printf("BCSR: rowPointer    = {");
    for(int i = 0; i < NAtoms+1; i++)
    {
      printf(" %d", rowPointer[i]);
    }
    printf(" }\n");
    printf("BCSR: columnPointer = {");
    for(int i = 0; i < numberBlocks; i++)
    {
      printf(" %d", columnPointer[i]);
    }
    printf(" }\n");
    printf("BCSR: blockPointer  = {");
    for(int i = 0; i < numberBlocks; i++)
    {
      printf(" %d", blockPointer[i]);
    }
    printf(" }\n");
  }
  printf("BCSR: numberNonZero = %d\n", numberNonZero);
  printf("BCSR: numberBlocks  = %d\n", numberBlocks);
}

/** Put BCSR to file.
 *
 * @param filename The filename to write to.
 */
void BCSR::put (char *filename)
{
  ABORT("FIXME\n");
}

/** Write a BCSR matrix into MatrixMarket format.
 *
 * @param filename The filename to write to.
 *
 * @return A string that contains the matrix in MatrixMarket format.
 */
void BCSR::toMM (char *filename)
{
  int fd = open(filename, O_CREAT | O_EXCL | O_WRONLY, 00644);

  if(fd == -1)
  {
    if(errno == EEXIST)
    {
      ABORT("file \"%s\" already exists\n", filename);
    }

    else
    {
      ABORT("error accessing file: %s\n", strerror(errno));
    }
  }

  else
  {
    FILE *fstream = fdopen(fd, "w");

    fprintf(fstream, "%%%%MatrixMarket matrix coordinate double general\n");
    fprintf(fstream, "%% %d x %d --> %d elements\n", M, N, M*N);
    fprintf(fstream, "%d %d %d\n", M, N, numberNonZero);

    for(int atom = 0; atom < NAtoms; atom++)
    {
      int MBlock = blockSize[atom];
      int rowOffset = offset[atom];

      for(int iRow = rowPointer[atom]; iRow < rowPointer[atom+1]; iRow++)
      {
        int NBlock = blockSize[columnPointer[iRow]];
        int columnOffset = offset[columnPointer[iRow]];

        for(int iSMat = 0; iSMat < NSMat; iSMat++)
        {
          int pointer = blockPointer[iRow]+iSMat*MBlock*NBlock;
          int i_block;
          int j_block;

          switch(iSMat+1)
          {
            case 1:
              i_block = rowOffset;
              j_block = columnOffset;
              break;

            case 2:
              i_block = rowOffset;
              j_block = columnOffset+numberBasisFunctions;
              break;

            case 3:
              i_block = rowOffset+numberBasisFunctions;
              j_block = columnOffset;
              break;

            case 4:
              i_block = rowOffset+numberBasisFunctions;
              j_block = columnOffset+numberBasisFunctions;
              break;

            default:
              ABORT("error\n");
              break;
          }

          for(int i = 0; i < MBlock; i++) {
            for(int j = 0; j < NBlock; j++)
            {
              fprintf(fstream, "%d %d % e\n", 1+i+i_block, 1+j+j_block, Aij);
            }
          }
        }
      }
    }

    if(fclose(fstream) != 0)
    {
      ABORT("error closing file\n");
    }
  }
}
