C------------------------------------------------------------------------------
C    This code is part of the MondoSCF suite of programs for linear scaling
C    electronic structure theory and ab initio molecular dynamics.
C
C    Copyright (2004). The Regents of the University of California. This
C    material was produced under U.S. Government contract W-7405-ENG-36
C    for Los Alamos National Laboratory, which is operated by the University
C    of California for the U.S. Department of Energy. The U.S. Government has
C    rights to use, reproduce, and distribute this software.  NEITHER THE
C    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
C    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
C
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by the
C    Free Software Foundation; either version 2 of the License, or (at your
C    option) any later version. Accordingly, this program is distributed in
C    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
C    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
C    PURPOSE. See the GNU General Public License at www.gnu.org for details.
C
C    While you may do as you like with this software, the GNU license requires
C    that you clearly mark derivative software.  In addition, you are encouraged
C    to return derivative works to the MondoSCF group for review, and possible
C    disemination in future releases.
C------------------------------------------------------------------------------
      SUBROUTINE QSDis(X, Y1, N, KFLAG)
      IMPLICIT  REAL*8(A-H,O-Z), INTEGER(I-N)
C***BEGIN PROLOGUE  QuickSort1
C***PURPOSE  Sort a REAL array and optionally make the same interchanges in
C            an auxiliary INTEGER array.  The array may be sorted in increasing
C            or decreasing order.  A slightly modified QUICKSORT
C            algorithm is used.
C***LIBRARY   SLATEC
C***CATEGORY  N6A2B
C***TYPE      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
C***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   SSORT sorts array X and optionally makes the same interchanges in
C   array Y.  The array X may be sorted in increasing order or
C   decreasing order.  A slightly modified quicksort algorithm is used.
C
C   Description of Parameters
C      X - array of values to be sorted   (usually abscissas)
C      Y - array to be (optionally) carried along
C      N - number of values in array X to be sorted
C      KFLAG - control parameter
C            =  2  means sort X in increasing order and carry Y along.
C            =  1  means sort X in increasing order (ignoring Y)
C            = -1  means sort X in decreasing order (ignoring Y)
C            = -2  means sort X in decreasing order and carry Y along.
C
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   761101  DATE WRITTEN
C   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891009  Removed unreferenced statement labels.  (WRB)
C   891024  Changed category.  (WRB)
C   891024  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   901012  Declared all variables; changed X,Y to SX,SY. (M. McClain)
C   920501  Reformatted the REFERENCES section.  (DWL, WRB)
C   920519  Clarified error messages.  (DWL)
C   920801  Declarations section rebuilt and code restructured to use
C           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
C***END PROLOGUE  SSORT
C     .. Scalar Arguments ..
      INTEGER KFLAG, N
C     .. Array Arguments ..
      REAL*8 X(N)
      INTEGER Y1(3,N)
C     .. Local Scalars ..
      REAL*8 R, T, TTK
      INTEGER TY,TZ,TTY,TTZ
      INTEGER I, IJ, J, K, KK, L, M, NN
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
C     .. External Subroutines ..
c      EXTERNAL XERMSG
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
C***FIRST EXECUTABLE STATEMENT  SSORT
      NN = N
      IF (NN .LT. 1) THEN
         STOP 'NN < 1 in QuickSortDis'
      ENDIF
C
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
         write(*,*)
     +      'The sort control parameter, K, is not 2, 1, -1, or -2.'
         STOP 'Problem in QuickSortDis.'
      ENDIF
C
C     Alter array X to get decreasing order if needed
C
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            X(I) = -X(I)
   10    CONTINUE
      ENDIF
C
      IF (KK .EQ. 2) GO TO 100
C
C     Sort X and carry Y along
C
  100 M = 1
      I = 1
      J = NN
      R = 0.375E0
C
  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
C
  120 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = X(IJ)
      TY1 = Y1(1,IJ)
      TY2 = Y1(2,IJ)
      TY3 = Y1(3,IJ)


C
C     If first element of array is greater than T, interchange with T
C
      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
         Y1(1,IJ) = Y1(1,I)
         Y1(2,IJ) = Y1(2,I)
         Y1(3,IJ) = Y1(3,I)

         Y1(1,I) = TY1
         Y1(2,I) = TY2
         Y1(3,I) = TY3

         TY1 = Y1(1,IJ)
         TY2 = Y1(2,IJ)
         TY3 = Y1(3,IJ)

      ENDIF
      L = J
C
C     If last element of array is less than T, interchange with T
C
      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)

         Y1(1,IJ) = Y1(1,J)
         Y1(2,IJ) = Y1(2,J)
         Y1(3,IJ) = Y1(3,J)

         Y1(1,J) = TY1
         Y1(2,J) = TY2
         Y1(3,J) = TY3

         TY1 = Y1(1,IJ)
         TY2 = Y1(2,IJ)
         TY3 = Y1(3,IJ)

C
C        If first element of array is greater than T, interchange with T
C
         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)

            Y1(1,IJ) = Y1(1,I)
            Y1(2,IJ) = Y1(2,I)
            Y1(3,IJ) = Y1(3,I)

            Y1(1,I) = TY1
            Y1(2,I) = TY2
            Y1(3,I) = TY3

            TY1 = Y1(1,IJ)
            TY2 = Y1(2,IJ)
            TY3 = Y1(3,IJ)

         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
  130 L = L-1
      IF (X(L) .GT. T) GO TO 130
C
C     Find an element in the first half of the array which is greater
C     than T
C
  140 K = K+1
      IF (X(K) .LT. T) GO TO 140
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT

         TTY1 = Y1(1,L)
         TTY2 = Y1(2,L)
         TTY3 = Y1(3,L)

         Y1(1,L) = Y1(1,K)
         Y1(2,L) = Y1(2,K)
         Y1(3,L) = Y1(3,K)

         Y1(1,K) = TTY1
         Y1(2,K) = TTY2
         Y1(3,K) = TTY3

         GO TO 130
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160
C
C     Begin again on another portion of the unsorted array
C
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1
C
  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = X(I+1)

      TY1 = Y1(1,I+1)
      TY2 = Y1(2,I+1)
      TY3 = Y1(3,I+1)

      IF (X(I) .LE. T) GO TO 170
      K = I
C
  180 X(K+1) = X(K)
      Y1(1,K+1) = Y1(1,K)
      Y1(2,K+1) = Y1(2,K)
      Y1(3,K+1) = Y1(3,K)


      K = K-1
      IF (T .LT. X(K)) GO TO 180
      X(K+1) = T
      Y1(1,K+1) = TY1
      Y1(2,K+1) = TY2
      Y1(3,K+1) = TY3

      GO TO 170
C
C     Clean up
C
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            X(I) = -X(I)
  200    CONTINUE
      ENDIF
      RETURN
      END
