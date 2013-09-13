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
!-----------------------------------------------------------------
!     Quick sort of doubles, carrying allong doubles
!
      SUBROUTINE DblIntSort77(N,DX,DY,KFLAG)
      INTEGER KFLAG, N
      REAL*8 DX(N)
      REAL*8 R, T, TT
      INTEGER DY(N)
      INTEGER I, IJ, J, K, KK, L, M, NN, TY, TTY
      INTEGER IL(21), IU(21)
      INTRINSIC ABS, INT
C
C      KFLAG=-2
C
      NN = N
      IF (NN .LT. 1) THEN
         write(*,*)
     +      'the number of values to be sorted is not positive.',nn
         CALL Trap()
         RETURN
      ENDIF
C
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
         write(*,*)
     +      'The sort control parameter, K, is not 2, 1, -1, or -2.'
         CALL Trap()
      ENDIF
C
C     Alter array DX to get decreasing order if needed
C
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            DX(I) = -DX(I)
   10    CONTINUE
      ENDIF
C
  100 M = 1
      I = 1
      J = NN
      R = 0.375D0
C
  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437D0) THEN
         R = R+3.90625D-2
      ELSE
         R = R-0.21875D0
      ENDIF
C
  120 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = DX(IJ)
      TY = DY(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (DX(I) .GT. T) THEN
         DX(IJ) = DX(I)
         DX(I) = T
         T = DX(IJ)
         DY(IJ) = DY(I)
         DY(I) = TY
         TY = DY(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than T, interchange with T
C
      IF (DX(J) .LT. T) THEN
         DX(IJ) = DX(J)
         DX(J) = T
         T = DX(IJ)
         DY(IJ) = DY(J)
         DY(J) = TY
         TY = DY(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (DX(I) .GT. T) THEN
            DX(IJ) = DX(I)
            DX(I) = T
            T = DX(IJ)
            DY(IJ) = DY(I)
            DY(I) = TY
            TY = DY(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
  130 L = L-1
      IF (DX(L) .GT. T) GO TO 130
C
C     Find an element in the first half of the array which is greater
C     than T
C
  140 K = K+1
      IF (DX(K) .LT. T ) GO TO 140
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = DX(L)
         DX(L) = DX(K)
         DX(K) = TT
         TTY = DY(L)
         DY(L) = DY(K)
         DY(K) = TTY
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
      IF (I .EQ. J ) GO TO 150
      T = DX(I+1)
      TY = DY(I+1)
      IF (DX(I) .LE. T) GO TO 170
      K = I
C
  180 DX(K+1) = DX(K)
      DY(K+1) = DY(K)
      K = K-1
      IF (T .LT. DX(K)) GO TO 180
      DX(K+1) = T
      DY(K+1) = TY
      GO TO 170
C
C     Clean up
C
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            DX(I) = -DX(I)
  200    CONTINUE
      ENDIF
      RETURN
      END

!-----------------------------------------------------------------
!     Quick sort modified to handle INTEGER*8 IX
!
      SUBROUTINE I8SORT (IX,IY,N,KFLAG)
C***BEGIN PROLOGUE  ISORT
C***PURPOSE  Sort an array and optionally make the same interchanges in
C            an auxiliary array.  The array may be sorted in increasing
C            or decreasing order.  A slightly modified QUICKSORT
C            algorithm is used.
C***LIBRARY   SLATEC
C***CATEGORY  N6A2A
C***TYPE      INTEGER (SSORT-S, DSORT-D, ISORT-I)
C***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Kahaner, D. K., (NBS)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   ISORT sorts array IX and optionally makes the same interchanges in
C   array IY.  The array IX may be sorted in increasing order or
C   decreasing order.  A slightly modified quicksort algorithm is used.
C
C   Description of Parameters
C      IX - integer array of values to be sorted
C      IY - integer array to be (optionally) carried along
C      N  - number of values in integer array IX to be sorted
C      KFLAG - control parameter
C            =  2  means sort IX in increasing order and carry IY along.
C            =  1  means sort IX in increasing order (ignoring IY)
C            = -1  means sort IX in decreasing order (ignoring IY)
C            = -2  means sort IX in decreasing order and carry IY along.
C
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   761118  DATE WRITTEN
C   810801  Modified by David K. Kahaner.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891009  Removed unreferenced statement labels.  (WRB)
C   891009  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   901012  Declared all variables; changed X,Y to IX,IY. (M. McClain)
C   920501  Reformatted the REFERENCES section.  (DWL, WRB)
C   920519  Clarified error messages.  (DWL)
C   920801  Declarations section rebuilt and code restructured to use
C           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
C***END PROLOGUE  ISORT
C     .. Scalar Arguments ..
      INTEGER KFLAG, N
C     .. Array Arguments ..
      INTEGER*8 IX(N)
      INTEGER   IY(N)
C     .. Local Scalars ..
      REAL R
      INTEGER*8 T,TT
      INTEGER I, IJ, J, K, KK, L, M, NN, TTY, TY
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
C***FIRST EXECUTABLE STATEMENT  ISORT
      NN = N
      IF (NN .LT. 1) THEN
C        The number of values to be sorted is not positive
         WRITE(*,*)' Died from error 1 in I8Sort '
         CALL Trap()
      ENDIF
C
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
C        The sort control parameter, K, is not 2, 1, -1, or -2
         WRITE(*,*)' Died from error 2 in I8Sort '
         CALL Trap()
      ENDIF
C
C     Alter array IX to get decreasing order if needed
C
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            IX(I) = -IX(I)
   10    CONTINUE
      ENDIF
C
      IF (KK .EQ. 2) GO TO 100
C
C     Sort IX only
C
      M = 1
      I = 1
      J = NN
      R = 0.375E0
C
   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
C
   30 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = IX(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (IX(I) .GT. T) THEN
         IX(IJ) = IX(I)
         IX(I) = T
         T = IX(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than than T, interchange with T
C
      IF (IX(J) .LT. T) THEN
         IX(IJ) = IX(J)
         IX(J) = T
         T = IX(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (IX(I) .GT. T) THEN
            IX(IJ) = IX(I)
            IX(I) = T
            T = IX(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
   40 L = L-1
      IF (IX(L) .GT. T) GO TO 40
C
C     Find an element in the first half of the array which is greater
C     than T
C
   50 K = K+1
      IF (IX(K) .LT. T) GO TO 50
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = IX(L)
         IX(L) = IX(K)
         IX(K) = TT
         GO TO 40
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
      GO TO 70
C
C     Begin again on another portion of the unsorted array
C
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1
C
   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = IX(I+1)
      IF (IX(I) .LE. T) GO TO 80
      K = I
C
   90 IX(K+1) = IX(K)
      K = K-1
      IF (T .LT. IX(K)) GO TO 90
      IX(K+1) = T
      GO TO 80
C
C     Sort IX and carry IY along
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
      T = IX(IJ)
      TY = IY(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (IX(I) .GT. T) THEN
         IX(IJ) = IX(I)
         IX(I) = T
         T = IX(IJ)
         IY(IJ) = IY(I)
         IY(I) = TY
         TY = IY(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than T, interchange with T
C
      IF (IX(J) .LT. T) THEN
         IX(IJ) = IX(J)
         IX(J) = T
         T = IX(IJ)
         IY(IJ) = IY(J)
         IY(J) = TY
         TY = IY(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (IX(I) .GT. T) THEN
            IX(IJ) = IX(I)
            IX(I) = T
            T = IX(IJ)
            IY(IJ) = IY(I)
            IY(I) = TY
            TY = IY(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
  130 L = L-1
      IF (IX(L) .GT. T) GO TO 130
C
C     Find an element in the first half of the array which is greater
C     than T
C
  140 K = K+1
      IF (IX(K) .LT. T) GO TO 140
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = IX(L)
         IX(L) = IX(K)
         IX(K) = TT
         TTY = IY(L)
         IY(L) = IY(K)
         IY(K) = TTY
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
      T = IX(I+1)
      TY = IY(I+1)
      IF (IX(I) .LE. T) GO TO 170
      K = I
C
  180 IX(K+1) = IX(K)
      IY(K+1) = IY(K)
      K = K-1
      IF (T .LT. IX(K)) GO TO 180
      IX(K+1) = T
      IY(K+1) = TY
      GO TO 170
C
C     Clean up
C
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            IX(I) = -IX(I)
  200    CONTINUE
      ENDIF
      RETURN
      END

      SUBROUTINE IntIntSort77(N,IA,JA, KFLAG)

C       modified to work with integer arrays IA and JA. WMC 97
c
C     ABSTRACT
C         This routine sorts an integer array IA and makes the same
C         interchanges in the integer array JA and the double precision
C         array A.
C         The array IA may be sorted in increasing order or decreasing
C         order.  A slightly modified quicksort algorithm is used.
C
C     DESCRIPTION OF PARAMETERS
C        IA - Integer array of values to be sorted.
C        JA - Integer array to be carried along.
C         N - Number of values in integer array IA to be sorted.
C     KFLAG - Control parameter
C           = 1 means sort IA in INCREASING order.
C           =-1 means sort IA in DECREASING order.
c
      implicit integer (i-m)
c
      real*8 R
c
      INTEGER N,NN
      DIMENSION IA(N), JA(N)
c
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
c
      INTRINSIC ABS, INT

      NN = N
      IF (NN.LT.1) THEN
        write(*,*)' The number of values to be sorted was not positive.'
        CALL Trap()
      ENDIF
      IF( N.EQ.1 ) RETURN
      KK = ABS(KFLAG)
      IF ( KK.NE.1 ) THEN
        write(*,*)'The sort control parameter, K, was not 1 or -1.'
        CALL Trap()
      ENDIF
C
C     Alter array IA to get decreasing order if needed.
C
      IF( KFLAG.LT.1 ) THEN
         DO 20 I=1,NN
            IA(I) = -IA(I)
 20      CONTINUE
      ENDIF
C
C     Sort IA and carry JA and A along.
C     And now...Just a little black magic...
      M = 1
      I = 1
      J = NN
      R = .375D0
 210  IF( R.LE.0.5898437D0 ) THEN
         R = R + 3.90625D-2
      ELSE
         R = R-.21875D0
      ENDIF
 225  K = I
C
C     Select a central element of the array and save it in location
C     it, jt, at.
C
      IJ = I + INT ((J-I)*R)
      IT = IA(IJ)
      JT = JA(IJ)
C
C     If first element of array is greater than it, interchange with it.
C
      IF( IA(I).GT.IT ) THEN
         IA(IJ) = IA(I)
         IA(I)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(I)
         JA(I)  = JT
         JT     = JA(IJ)
      ENDIF
      L=J
C
C     If last element of array is less than it, swap with it.
C
      IF( IA(J).LT.IT ) THEN
         IA(IJ) = IA(J)
         IA(J)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(J)
         JA(J)  = JT
         JT     = JA(IJ)
C
C     If first element of array is greater than it, swap with it.
C
         IF ( IA(I).GT.IT ) THEN
            IA(IJ) = IA(I)
            IA(I)  = IT
            IT     = IA(IJ)
            JA(IJ) = JA(I)
            JA(I)  = JT
            JT     = JA(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is
C     smaller than it.
C
  240 L=L-1
      IF( IA(L).GT.IT ) GO TO 240
C
C     Find an element in the first half of the array which is
C     greater than it.
C
  245 K=K+1
      IF( IA(K).LT.IT ) GO TO 245
C
C     Interchange these elements.
C
      IF( K.LE.L ) THEN
         IIT   = IA(L)
         IA(L) = IA(K)
         IA(K) = IIT
         JJT   = JA(L)
         JA(L) = JA(K)
         JA(K) = JJT
         GOTO 240
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted.
C
      IF( L-I.GT.J-K ) THEN
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
      GO TO 260
C
C     Begin again on another portion of the unsorted array.
C
  255 M = M-1
      IF( M.EQ.0 ) GO TO 300
      I = IL(M)
      J = IU(M)
  260 IF( J-I.GE.1 ) GO TO 225
      IF( I.EQ.J ) GO TO 255
      IF( I.EQ.1 ) GO TO 210
      I = I-1
  265 I = I+1
      IF( I.EQ.J ) GO TO 255
      IT = IA(I+1)
      JT = JA(I+1)
      IF( IA(I).LE.IT ) GO TO 265
      K=I
  270 IA(K+1) = IA(K)
      JA(K+1) = JA(K)
      K = K-1
      IF( IT.LT.IA(K) ) GO TO 270
      IA(K+1) = IT
      JA(K+1) = JT
      GO TO 265
C
C     Clean up, if necessary.
C
  300 IF( KFLAG.LT.1 ) THEN
         DO 310 I=1,NN
            IA(I) = -IA(I)
 310     CONTINUE
      ENDIF
      RETURN
      END



      SUBROUTINE IntSort77(N,IA,KFLAG)

C       modified to work with only the integer array IA. WMC 99
c
C     ABSTRACT
C         This routine sorts an integer array IA and makes the same
C         interchanges in the integer array JA and the double precision
C         array A.
C         The array IA may be sorted in increasing order or decreasing
C         order.  A slightly modified quicksort algorithm is used.
C
C     DESCRIPTION OF PARAMETERS
C        IA - Integer array of values to be sorted.
C        JA - Integer array to be carried along.
C         N - Number of values in integer array IA to be sorted.
C     KFLAG - Control parameter
C           = 1 means sort IA in INCREASING order.
C           =-1 means sort IA in DECREASING order.
c
      implicit integer (i-m)
c
      real*8 R
c
      INTEGER N,NN
      DIMENSION IA(N)
c
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
c
      INTRINSIC ABS, INT

      NN = N
      IF (NN.LT.1) THEN
        write(*,*)' The number of values to be sorted was not positive.'
        CALL Trap()
      ENDIF
      IF( N.EQ.1 ) RETURN
      KK = ABS(KFLAG)
      IF ( KK.NE.1 ) THEN
        write(*,*)'The sort control parameter, K, was not 1 or -1.'
        CALL Trap()
      ENDIF
C
C     Alter array IA to get decreasing order if needed.
C
      IF( KFLAG.LT.1 ) THEN
         DO 20 I=1,NN
            IA(I) = -IA(I)
 20      CONTINUE
      ENDIF
C
C     Sort IA and carry JA and A along.
C     And now...Just a little black magic...
      M = 1
      I = 1
      J = NN
      R = .375D0
 210  IF( R.LE.0.5898437D0 ) THEN
         R = R + 3.90625D-2
      ELSE
         R = R-.21875D0
      ENDIF
 225  K = I
C
C     Select a central element of the array and save it in location
C     it, jt, at.
C
      IJ = I + INT ((J-I)*R)
      IT = IA(IJ)
C
C     If first element of array is greater than it, interchange with it.
C
      IF( IA(I).GT.IT ) THEN
         IA(IJ) = IA(I)
         IA(I)  = IT
         IT     = IA(IJ)
      ENDIF
      L=J
C
C     If last element of array is less than it, swap with it.
C
      IF( IA(J).LT.IT ) THEN
         IA(IJ) = IA(J)
         IA(J)  = IT
         IT     = IA(IJ)
C
C     If first element of array is greater than it, swap with it.
C
         IF ( IA(I).GT.IT ) THEN
            IA(IJ) = IA(I)
            IA(I)  = IT
            IT     = IA(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is
C     smaller than it.
C
  240 L=L-1
      IF( IA(L).GT.IT ) GO TO 240
C
C     Find an element in the first half of the array which is
C     greater than it.
C
  245 K=K+1
      IF( IA(K).LT.IT ) GO TO 245
C
C     Interchange these elements.
C
      IF( K.LE.L ) THEN
         IIT   = IA(L)
         IA(L) = IA(K)
         IA(K) = IIT
         GOTO 240
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted.
C
      IF( L-I.GT.J-K ) THEN
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
      GO TO 260
C
C     Begin again on another portion of the unsorted array.
C
  255 M = M-1
      IF( M.EQ.0 ) GO TO 300
      I = IL(M)
      J = IU(M)
  260 IF( J-I.GE.1 ) GO TO 225
      IF( I.EQ.J ) GO TO 255
      IF( I.EQ.1 ) GO TO 210
      I = I-1
  265 I = I+1
      IF( I.EQ.J ) GO TO 255
      IT = IA(I+1)
      IF( IA(I).LE.IT ) GO TO 265
      K=I
  270 IA(K+1) = IA(K)
      K = K-1
      IF( IT.LT.IA(K) ) GO TO 270
      IA(K+1) = IT
      GO TO 265
C
C     Clean up, if necessary.
C
  300 IF( KFLAG.LT.1 ) THEN
         DO 310 I=1,NN
            IA(I) = -IA(I)
 310     CONTINUE
      ENDIF
      RETURN
      END
