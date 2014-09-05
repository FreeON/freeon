!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
!    SORTING AND ORDERING (VIA SPACE FILLING CURVES)
!    Author: Matt Challacombe
MODULE Order
   USE DerivedTypes
   USE GlobalCharacters
   USE GlobalScalars
   USE GlobalObjects
   USE ProcessControl
   USE MondoLogger
   USE MemMan
   USE Opt_Trav_Band ; USE AnnealMap
   USE Parse
   IMPLICIT NONE
   INTERFACE Sort
      MODULE PROCEDURE Sort_DBL_INT,Sort_INT_INT, Sort_INT_VECT
   END INTERFACE
   INTERFACE Random
      MODULE PROCEDURE RANDOM_INT,RANDOM_DBL
   END INTERFACE
   INTERFACE
      SUBROUTINE DblIntSort77(N,X,Y,Ordr)
         USE DerivedTypes
         INTEGER                      :: N,Ordr
         REAL(DOUBLE),DIMENSION(N)    :: X
         INTEGER,     DIMENSION(N)    :: Y
      END SUBROUTINE
      SUBROUTINE SFCOrder77(N,R,Point,Key,Hilbert)
         USE DerivedTypes
         INTEGER                      :: N
         LOGICAL                      :: Hilbert
         REAL(DOUBLE), DIMENSION(3,N) :: R
         INTEGER(INT8),DIMENSION(N)   :: Key
         INTEGER,      DIMENSION(N)   :: Point
      END SUBROUTINE

!      FUNCTION Interleave(Ix,Iy,Iz)
!         IMPLICIT NONE
!         INTEGER,INTENT(IN) :: Ix,Iy,Iz
!         INTEGER, PARAMETER :: INT8=SELECTED_INT_KIND(18)
!         INTEGER(INT8)            :: Interleave
!       END FUNCTION Interleave

   END INTERFACE
   CONTAINS

       FUNCTION RANDOM_INT(Limits)
          INTEGER                :: RANDOM_INT
          INTEGER, DIMENSION(2)  :: Limits
          REAL(DOUBLE)           :: Delta
          REAL(DOUBLE), EXTERNAL :: Random

          Delta=DBLE(Limits(2)-Limits(1)+1)
          RANDOM_INT=Limits(1)+INT(Delta*Random())

       END FUNCTION RANDOM_INT

      FUNCTION RANDOM_DBL(Limits)
         REAL(DOUBLE)              :: RANDOM_DBL
         REAL(DOUBLE),DIMENSION(2) :: Limits
         REAL(DOUBLE)              :: Delta
         REAL(DOUBLE), EXTERNAL    :: Random

         Delta=Limits(2)-Limits(1)+0.0D0
         RANDOM_DBL=Limits(1)+Delta*Random()

      END FUNCTION RANDOM_DBL

! the following function is buggy
!     FUNCTION RANDOM_INT(Limits)
!        INTEGER               :: RANDOM_INT,Delta
!        INTEGER, SAVE         :: JRan=10408
!        INTEGER, DIMENSION(2) :: Limits
!        INTEGER, PARAMETER    :: Im=259200,Ia=7141,Ic=54773
!        JRan=MOD(JRan*Ia+Ic,Im)
!        Delta=Limits(2)-Limits(1)+1
!        RANDOM_INT=Limits(1)+(Delta*JRan)/Im
!        IF(RANDOM_INT>Limits(2).OR.RANDOM_INT<Limits(1)) &
!           CALL Halt(' Limits hosed in RANDOM_INT ')
!     END FUNCTION RANDOM_INT


!--------------------------------------------------------------
!    F90 wrapper for SFCOrder77, which circumvents the lack
!    of INTEGER(KIND=8) (INTEGER*8) support for cheazy
!    F90 compilers (pgf,nag...)
     SUBROUTINE SFCOrder(N,R,Point,SFC_KEY)
        INTEGER                       :: N
        INTEGER                       :: SFC_KEY
        TYPE(DBL_RNK2)                :: R
        TYPE(INT_VECT)                :: Point
!
        INTEGER(INT8),ALLOCATABLE, &
                         DIMENSION(:) :: IKey
        TYPE(DBL_VECT)                :: RKey
        INTEGER                       :: I
        INTEGER::N_Neighbor
        REAL(DOUBLE)::Gamma,Alpha,NLFac

        IF(SFC_KEY==SFC_RANDOM)THEN
           CALL MondoLog(DEBUG_NONE, "Order", "Random order")
           CALL New(RKey,N)
           DO I=1,N
              Point%I(I)=I
              CALL RANDOM_NUMBER(RKey%D(I))
           ENDDO
           CALL Sort_DBL_INT(RKey,Point,N)
           CALL Delete(RKey)
        ELSEIF(SFC_KEY==SFC_PEANO)THEN
           CALL MondoLog(DEBUG_NONE, "Order", "Peano order")
           ALLOCATE(IKey(N))
           CALL SFCOrder77(N,R%D,Point%I,IKey,.FALSE.)
           DEALLOCATE(IKey)
        ELSEIF(SFC_KEY==SFC_HILBERT)THEN
           CALL MondoLog(DEBUG_NONE, "Order", "Hilbert order")
           ALLOCATE(IKey(N))
           CALL SFCOrder77(N,R%D,Point%I,IKey,.TRUE.)
           DEALLOCATE(IKey)
        ELSEIF(SFC_KEY==SFC_TRAVEL)THEN
           CALL MondoLog(DEBUG_NONE, "Order", "Travel order")
          CALL OpenASCII(InpFile,Inp)
          IF(.NOT.OptDblQ(Inp,'Gamma',Gamma)) THEN
            WRITE(*,*) 'Cannot find Gamma in ',InpFile
            STOP
          ENDIF
          IF(.NOT.OptDblQ(Inp,'NLFac',NLFac)) THEN
            WRITE(*,*) 'Cannot find NLFac in ',InpFile
            STOP
          ENDIF
          IF(.NOT.OptDblQ(Inp,'Alpha',Alpha)) THEN
            WRITE(*,*) 'Cannot find Alpha in ',InpFile
            STOP
          ENDIF
          IF(.NOT.OptIntQ(Inp,'N_Neighbor',N_Neighbor)) THEN
            WRITE(*,*) 'Cannot find N_Neighbor in ',InpFile
            STOP
          ENDIF
          Close(Unit=Inp,STATUS='KEEP')
          CALL Anneal(R,N,Point,Gamma,NLFac,N_Neighbor,Alpha)
        ELSEIF(SFC_KEY==SFC_TABLETRAV)THEN
          CALL TableAnneal(R,N,Point)
        ELSE
           CALL MondoHalt(-100,'Bad SFC_Key in SFCOrder')
        ENDIF
     END SUBROUTINE SFCOrder

     SUBROUTINE Sort_DBL_INT(X,Y,N_O,Ordr_O)
        TYPE(DBL_VECT)     :: X
        TYPE(INT_VECT)     :: Y
        INTEGER,OPTIONAL   :: N_O,Ordr_O
        INTEGER                       :: N,Ordr
        Ordr=-2
        N=MIN(SIZE(X%D),SIZE(Y%I))
        IF(PRESENT(N_O))THEN
           IF(N_O>N)CALL Halt(' Dimensioning off in DblSort ')
           N=N_O
        ENDIF
        IF(PRESENT(Ordr_O))Ordr=Ordr_O
        CALL DblIntSort77(N,X%D,Y%I,Ordr)
    END SUBROUTINE Sort_DBL_INT


     SUBROUTINE Sort_INT_INT(X,Y,N_O,Ordr_O)
        TYPE(INT_VECT)     :: X
        TYPE(INT_VECT)     :: Y
        INTEGER,OPTIONAL   :: N_O,Ordr_O
        INTEGER            :: N,Ordr
        Ordr=-1
        N=MIN(SIZE(X%I),SIZE(Y%I))
        IF(PRESENT(N_O))THEN
           IF(N_O>N)CALL Halt(' Dimensioning off in DblSort ')
           N=N_O
        ENDIF
        IF(PRESENT(Ordr_O))Ordr=Ordr_O
        CALL IntIntSort77(N,X%I,Y%I,Ordr)
    END SUBROUTINE Sort_INT_INT


     SUBROUTINE Sort_INT_VECT(X,N_O,Ordr_O)
        TYPE(INT_VECT)     :: X
        INTEGER,OPTIONAL   :: N_O,Ordr_O
        INTEGER            :: N,Ordr
        Ordr=-1
        N=SIZE(X%I)
        IF(PRESENT(N_O))THEN
           IF(N_O>N)CALL Halt(' Dimensioning off in DblSort ')
           N=N_O
        ENDIF
        IF(PRESENT(Ordr_O))Ordr=Ordr_O
        CALL IntSort77(N,X%I,Ordr)
    END SUBROUTINE Sort_INT_VECT
  END MODULE Order

  SUBROUTINE SFCOrder77(N,R,Point,Key,Hilbert)
    Use DerivedTypes
    IMPLICIT NONE
    INTEGER,PARAMETER::BitNum=21 ! number of bits to represent an integer
    INTEGER::IMax,I,J,N,Ix,Iy,Iz,Point(N)
    INTEGER(INT8):: Key(N)
    INTEGER(INT8),EXTERNAL:: Interleave,HilbertKey
    REAL(DOUBLE)::R(3,N),RMin(3),Ratio,MaxDiff,Diff
    LOGICAL::Hilbert

    RMin(1)=1.D10; RMin(2)=1.D10; RMin(3)=1.D10

    DO J=1,N
       Point(J)=J
       RMin(1)=DMIN1(RMin(1),R(1,J))
       RMin(2)=DMIN1(RMin(2),R(2,J))
       RMin(3)=DMIN1(RMin(3),R(3,J))
    ENDDO

    MaxDiff = -1.0D0 ! any negative value will do
    DO J = 1,N
      DO I = 1, 3
        Diff = R(I,J) - RMin(I)
        IF( Diff > MaxDiff ) THEN
          MaxDiff = Diff
        ENDIF
      ENDDO
    ENDDO

    IF( MaxDiff .LE. 0) THEN
      WRITE(*,*) 'MaxDiff = ',MaxDiff
      STOP 'ERROR: MaxDiff must be positive!'
    ENDIF
    Ratio = (2**BitNum-1.0)/MaxDiff

    IMax = -Big_Int
    IF(Hilbert)THEN
       DO J=1,N
          Ix=DNINT((R(1,J)-RMin(1))*Ratio)
          IF(Ix > IMax) IMax = Ix
          Iy=DNINT((R(2,J)-RMin(2))*Ratio)
          IF(Iy > IMax) IMax = Iy
          Iz=DNINT((R(3,J)-RMin(3))*Ratio)
          IF(Iz > IMax) IMax = Iz
          Key(J)=Interleave(BitNum,Ix,Iy,Iz)
          Key(J)=HilbertKey(BitNum,Key(J))
       END DO

       IF( IMax /= (2**BitNum-1)) THEN
!        stop 'ERROR: numerical accuracy problem !!'
       END IF
    ELSE ! Peano option
       DO J=1,N
          Ix=DNINT((R(1,J)-RMin(1))*Ratio)
          Iy=DNINT((R(2,J)-RMin(2))*Ratio)
          Iz=DNINT((R(3,J)-RMin(3))*Ratio)
          Key(J)=Interleave(BitNum,Ix,Iy,Iz)
       END DO
    ENDIF
    CALL I8Sort(Key,Point,N,2)
  END SUBROUTINE SFCOrder77

!-----------------------------------------------------------------
!   Bit shuffle for a triple of integers
  FUNCTION Interleave(BitNum,Ix,Iy,Iz)
     USE DerivedTypes
     IMPLICIT NONE
     INTEGER I,K,Ix,Iy,Iz,BitNum,EndNum
     INTEGER(INT8) Interleave
     K=0
     Interleave=0
     EndNum = 3*(BitNum-1)
     DO I=0,EndNum,3
        IF(BTEST(Iz,K))Interleave=IBSET(Interleave,I  )
        IF(BTEST(Iy,K))Interleave=IBSET(Interleave,I+1)
        IF(BTEST(Ix,K))Interleave=IBSET(Interleave,I+2)
        K=K+1
     ENDDO
  END FUNCTION Interleave

!-----------------------------------------------------------------
! Hilbert ordering for three-dimensional data, (3*BitNum)-bit
! implementation (stored in 64-bit integer).
! Interleaved key based on Algortithm H2 of Faloutsos and Roseman, PODS 89.
  FUNCTION HilbertKey(BitNum,Key)
     USE DerivedTypes
     IMPLICIT NONE
     INTEGER I,J,K,L,BitNum,BeginNum
     INTEGER(INT8) Key
     INTEGER(INT8) HilbertKey
     INTEGER SubKey(21)
     INTEGER ToState(0:7,12)
     INTEGER ToBinry(0:7,12)
!--------------------------------------------------------------
!  State table, based on Fig 3 , T. Bially,
!  IEEE Trans. on Information Theory, 1969, pp. 658, vol. IT-15
!                                 0, 1, 2, 3, 4, 5, 6, 7
     DATA (ToState(I,1 ),I=0,7) / 9,12, 2, 1, 2, 3, 9, 1/
     DATA (ToBinry(I,1 ),I=0,7) / 0, 1, 3, 2, 7, 6, 4, 5/
     DATA (ToState(I,2 ),I=0,7) / 5, 3, 2, 2, 3, 5, 1, 4/
     DATA (ToBinry(I,2 ),I=0,7) / 4, 3, 5, 2, 7, 0, 6, 1/
     DATA (ToState(I,3 ),I=0,7) / 2, 3, 6, 3, 1, 7, 7, 1/
     DATA (ToBinry(I,3 ),I=0,7) / 6, 5, 1, 2, 7, 4, 0, 3/
     DATA (ToState(I,4 ),I=0,7) /10, 9, 4, 2, 5, 2, 4, 9/
     DATA (ToBinry(I,4 ),I=0,7) / 6, 7, 5, 4, 1, 0, 2, 3/
     DATA (ToState(I,5 ),I=0,7) / 5, 2, 5, 6, 8, 4, 4, 8/
     DATA (ToBinry(I,5 ),I=0,7) / 2, 1, 5, 6, 3, 0, 4, 7/
     DATA (ToState(I,6 ),I=0,7) / 6, 6, 5, 3, 7, 8, 3, 5/
     DATA (ToBinry(I,6 ),I=0,7) / 2, 5, 3, 4, 1, 6, 0, 7/
     DATA (ToState(I,7 ),I=0,7) / 6, 7,11,12,11, 7, 6, 3/
     DATA (ToBinry(I,7 ),I=0,7) / 4, 5, 7, 6, 3, 2, 0, 1/
     DATA (ToState(I,8 ),I=0,7) / 8, 6,10,11, 8,11, 5, 6/
     DATA (ToBinry(I,8 ),I=0,7) / 2, 3, 1, 0, 5, 4, 6, 7/
     DATA (ToState(I,9 ),I=0,7) /12,10, 1, 4,10,12, 9, 9/
     DATA (ToBinry(I,9 ),I=0,7) / 0, 7, 1, 6, 3, 4, 2, 5/
     DATA (ToState(I,10),I=0,7) / 8, 4, 4, 8,10, 9,10,11/
     DATA (ToBinry(I,10),I=0,7) / 4, 7, 3, 0, 5, 6, 2, 1/
     DATA (ToState(I,11),I=0,7) / 7, 8,12,10,11,11,10,12/
     DATA (ToBinry(I,11),I=0,7) / 6, 1, 7, 0, 5, 2, 4, 3/
     DATA (ToState(I,12),I=0,7) / 1, 7, 7, 1, 9,12,11,12/
     DATA (ToBinry(I,12),I=0,7) / 0, 3, 7, 4, 1, 2, 6, 5/

! Form 3 bit SubKeys from left to right (Step 3 of H2)
     DO I=1,BitNum
        SubKey(I)=0
     ENDDO
     K=1
     BeginNum=3*(BitNum-1)
     DO I=BeginNum,0,-3
        IF(BTEST(Key,I  ))SubKey(K)=IBSET(SubKey(K),0)
        IF(BTEST(Key,I+1))SubKey(K)=IBSET(SubKey(K),1)
        IF(BTEST(Key,I+2))SubKey(K)=IBSET(SubKey(K),2)
        IF(K.GT.BitNum)STOP ' gt BitNum '
        IF(I.LT.0)STOP ' lt 0 '
        K=K+1
     ENDDO
! Change each 3 bit SubKey according to the output key
! from the current state and move to a new state (Step 4 of H2)
     I=1
     DO K=1,BitNum
        L=SubKey(K)
        SubKey(K)=ToBinry(L,I)
        I        =ToState(L,I)
     ENDDO

! Reassemble key from sub keys, with SubKey(21)
! the rightmost 3 bits (Step 5 of H2)
     HilbertKey=0
     I=0
     DO K=BitNum,1,-1
        IF(BTEST(SubKey(K),0))HilbertKey=IBSET(HilbertKey,I  )
        IF(BTEST(SubKey(K),1))HilbertKey=IBSET(HilbertKey,I+1)
        IF(BTEST(SubKey(K),2))HilbertKey=IBSET(HilbertKey,I+2)
        I=I+3
     ENDDO
  END FUNCTION HilbertKey
