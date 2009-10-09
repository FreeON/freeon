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
MODULE CholFactor
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalObjects
   USE GlobalCharacters
   USE InOut
   USE MemMan
   USE PrettyPrint
   USE ParsingConstants
   USE MondoLogger

IMPLICIT NONE

CONTAINS
!
!---------------------------------------------------------------
!
   SUBROUTINE GetBorders(Border,IA,JA,NDim)
     INTEGER,DIMENSION(:) :: Border,IA,JA
     INTEGER :: NCart,I,J,K,L,NDim
     ! pass in ordered matrix A !
     DO I=1,NDim
       IF(IA(I+1)==IA(I)) THEN
         Border(I)=I
       ELSE
         Border(I)=JA(IA(I))
       ENDIF
     ENDDO
   END SUBROUTINE GetBorders
!
!---------------------------------------------------------------
!
   SUBROUTINE SymbCholComplete(Border,N,IU,JU)
     INTEGER,DIMENSION(:) :: Border,IU,JU
     REAL(DOUBLE),DIMENSION(1) :: AN
     TYPE(INT_VECT) :: IL,JL
     INTEGER :: I,J,K,L,N,NZ
     !
     NZ=SIZE(JU)
     CALL New(IL,N+1)
     CALL New(JL,NZ)
     NZ=0
     IL%I(1)=1
     DO I=1,N
       DO J=Border(I),I-1
         NZ=NZ+1
         JL%I(NZ)=J
       ENDDO
       IL%I(I+1)=NZ+1
     ENDDO
     CALL TransPose1x1(IL%I,JL%I,AN,N,N,IU,JU,AN,'symb')
     CALL Delete(IL)
     CALL Delete(JL)
   END SUBROUTINE SymbCholComplete
!
!---------------------------------------------------------------
!
   SUBROUTINE SymbCholInCompl(IA,JA,N,IU,JU)
     TYPE(INT_VECT)       :: IP
     INTEGER,DIMENSION(:) :: IA,JA,IU,JU
     INTEGER              :: I,J,K,L,N
     INTEGER              :: NM,NH,JP,JPI,JPP,MIN,IAA,IAB,JJ
     INTEGER              :: LAST,LH,IUA,IUB
     !
     CALL New(IP,N)
     NM=N-1
     NH=N+1
     IU(:)=0
     IP%I(:)=0
     10      CONTINUE
     JP=1
     DO 90 I=1,NM
       JPI=JP
       JPP=N+JP-1
       MIN=NH
       IAA=IA(I)
       IAB=IA(I+1)-1
       IF(IAB.LT.IAA) GO TO 30
       DO 20 J=IAA,IAB
         JJ=JA(J)
         JU(JP)=JJ
         JP=JP+1
         IF(JJ.LT.MIN) MIN=JJ
         IU(JJ)=I
       20 CONTINUE
       30 LAST=IP%I(I)
       IF(LAST.EQ.0) GO TO 60
       L=LAST
       40 L=IP%I(L)
       LH=L+1
       IUA=IU(L)
       IUB=IU(LH)-1
       IF(LH.EQ.I) IUB=JPI-1
       IU(I)=I
       DO 50 J=IUA,IUB
         JJ=JU(J)
         IF(IU(JJ).EQ.I) GO TO 50
         JU(JP)=JJ
         JP=JP+1
         IU(JJ)=I
         IF(JJ.LT.MIN) MIN=JJ
       50 CONTINUE
       IF(JP.EQ.JPP) GO TO 70
       IF(L.NE.LAST) GO TO 40
       60 IF(MIN.EQ.NH) GO TO 90
       70 L=IP%I(MIN)
       IF(L.EQ.0) GO TO 80
       IP%I(I)=IP%I(L)
       IP%I(L)=I
       GO TO 90
       80 IP%I(MIN)=I
       IP%I(I)=I
     90 IU(I)=JPI
     IU(N)=JP
     IU(NH)=JP
     !
     CALL Delete(IP)
   END SUBROUTINE SymbCholInCompl
!
!----------------------------------------------------------------
!
   SUBROUTINE NumCholFact(IA,JA,AN,AD,IU,JU,N,UN,DI,DoAbs_O)
     TYPE(INT_VECT)       :: IP,IUP
     INTEGER,DIMENSION(:) :: IA,JA,IU,JU
     REAL(DOUBLE),DIMENSION(:) :: AN,AD,DI,UN
     REAL(DOUBLE)              :: UM
     INTEGER              :: N,I,IH,IUA,IUB,IAA,IAB,J,LAST,LN,L
     INTEGER              :: IUC,IUD,JJ
     LOGICAL,OPTIONAL     :: DoAbs_O
     LOGICAL              :: DoAbs
     !
     IF(PRESENT(DoAbs_O)) THEN
       DoAbs=DoAbs_O
     ELSE
       DoAbs=.FALSE.
     ENDIF
     !
     CALL New(IP,N)
     CALL New(IUP,N)
     !
     IP%I=0
     DO 130 I=1,N
       IH=I+1
       IUA=IU(I)
       IUB=IU(IH)-1
       IF(IUB.LT.IUA) GO TO 40
       DO 20 J=IUA,IUB
     20 DI(JU(J))=Zero
       IAA=IA(I)
       IAB=IA(IH)-1
       IF(IAB.LT.IAA) GO TO 40
       DO 30 J=IAA,IAB
       30 DI(JA(J))=AN(J)
       40 DI(I)=AD(I)
       LAST=IP%I(I)
       IF(LAST.EQ.0) GO TO 90
       LN=IP%I(LAST)
       50 L=LN
       LN=IP%I(L)
       IUC=IUP%I(L)
       IUD=IU(L+1)-1
       UM=UN(IUC)*DI(L)
       DO 60 J=IUC,IUD
         JJ=JU(J)
       60 DI(JJ)=DI(JJ)-UN(J)*UM
       UN(IUC)=UM
       IUP%I(L)=IUC+1
       IF(IUC.EQ.IUD) GO TO 80
       J=JU(IUC+1)
       JJ=IP%I(J)
       IF(JJ.EQ.0) GO TO 70
       IP%I(L)=IP%I(JJ)
       IP%I(JJ)=L
       GO TO 80
       70 IP%I(J)=L
       IP%I(L)=L
       80 IF(L.NE.LAST) GO TO 50
       90 DI(I)=One/DI(I)  !!!! negative diagonal may occure here
       IF(IUB.LT.IUA) GO TO 120
       DO 100 J=IUA,IUB
       100 UN(J)=DI(JU(J))
       J=JU(IUA)
       JJ=IP%I(J)
       IF(JJ.EQ.0) GO TO 110
       IP%I(I)=IP%I(JJ)
       IP%I(JJ)=I
       GO TO 120
       110 IP%I(J)=I
       IP%I(I)=I
       120 IUP%I(I)=IUA
       130 CONTINUE
       !
     CALL Delete(IP)
     CALL Delete(IUP)
   END SUBROUTINE NumCholFact
!
!---------------------------------------------------------------------
!
   SUBROUTINE Perm1x1(Perm,IA,JA,AN,Symb_O)
     INTEGER,DIMENSION(:)      :: Perm,IA,JA
     REAL(DOUBLE),DIMENSION(:) :: AN
     TYPE(INT_VECT)            :: IAT1,JAT1
     TYPE(DBL_VECT)            :: ANT1
     INTEGER                   :: I,J,NZ,N,K,L
     CHARACTER(LEN=DCL)        :: Char
     LOGICAL,OPTIONAL          :: Symb_O
     !
     Char='full'
     IF(PRESENT(Symb_O)) THEN
       IF(Symb_O) THEN
         Char='Symb'
       ENDIF
     ENDIF
     !
     N=SIZE(IA)-1
     NZ=IA(N+1)-1
     DO I=1,N
       DO J=IA(I),IA(I+1)-1
         K=JA(J)
         JA(J)=Perm(K)
       ENDDO
     ENDDO
     CALL New(IAT1,N+1)
     CALL New(JAT1,NZ)
     CALL New(ANT1,NZ)
     CALL TransPose1x1(IA,JA,AN,N,N,IAT1%I,JAT1%I,ANT1%D,TRIM(Char))
     DO I=1,N
       DO J=IAT1%I(I),IAT1%I(I+1)-1
         K=JAT1%I(J)
         JAT1%I(J)=Perm(K)
       ENDDO
     ENDDO
     CALL TransPose1x1(IAT1%I,JAT1%I,ANT1%D,N,N,IA,JA,AN,TRIM(Char))
     CALL Delete(IAT1)
     CALL Delete(JAT1)
     CALL Delete(ANT1)
   END SUBROUTINE Perm1x1
!
!---------------------------------------------------------------------
!
   SUBROUTINE ITransPose1x1(IA,JA,AN,N,M,IAT,JAT,ANT,Char)
     INTEGER,DIMENSION(:) :: IA,JA,IAT,JAT
     INTEGER,DIMENSION(:) :: AN,ANT
     INTEGER :: N,M,I,J,K,L,MH,NH,IAB,IAA,JP
     CHARACTER(LEN=*) :: Char
     ! transpose an NxM matrix 'A'
     MH=M+1
     NH=N+1
     IAT(:)=0
     IAB=IA(NH)-1
     DO 20 I=1,IAB
       J=JA(I)+2
       IF(J.LE.MH) IAT(J)=IAT(J)+1
     20 CONTINUE
     IAT(1)=1
     IAT(2)=1
     IF(M.EQ.1) GO TO 40
     DO 30 I=3,MH
     30 IAT(I)=IAT(I)+IAT(I-1)
     40 DO 60 I=1,N
     IAA=IA(I)
     IAB=IA(I+1)-1
     IF(IAB.LT.IAA) GO TO 60
     IF(Char=='full') THEN
       DO 50 JP=IAA,IAB !!!! symbolic and numeric
         J=JA(JP)+1
         K=IAT(J)
         JAT(K)=I
         ANT(K)=AN(JP)
       50 IAT(J)=K+1
     ELSE
       DO 51 JP=IAA,IAB !!!! symbolic only
         J=JA(JP)+1
         K=IAT(J)
         JAT(K)=I
       51 IAT(J)=K+1
     ENDIF
     60 CONTINUE
     !
   END SUBROUTINE ITransPose1x1
!
!---------------------------------------------------------------------
!
   SUBROUTINE TransPose1x1(IA,JA,AN,N,M,IAT,JAT,ANT,Char)
     INTEGER,DIMENSION(:) :: IA,JA,IAT,JAT
     REAL(DOUBLE),DIMENSION(:) :: AN,ANT
     INTEGER :: N,M,I,J,K,L,MH,NH,IAB,IAA,JP
     CHARACTER(LEN=*) :: Char
     ! transpose an NxM matrix 'A'
     MH=M+1
     NH=N+1
     IAT(:)=0
     IAB=IA(NH)-1
     DO 20 I=1,IAB
       J=JA(I)+2
       IF(J.LE.MH) IAT(J)=IAT(J)+1
     20 CONTINUE
     IAT(1)=1
     IAT(2)=1
     IF(M.EQ.1) GO TO 40
     DO 30 I=3,MH
     30 IAT(I)=IAT(I)+IAT(I-1)
     40 DO 60 I=1,N
     IAA=IA(I)
     IAB=IA(I+1)-1
     IF(IAB.LT.IAA) GO TO 60
     IF(Char=='full') THEN
       DO 50 JP=IAA,IAB !!!! symbolic and numeric
         J=JA(JP)+1
         K=IAT(J)
         JAT(K)=I
         ANT(K)=AN(JP)
       50 IAT(J)=K+1
     ELSE
       DO 51 JP=IAA,IAB !!!! symbolic only
         J=JA(JP)+1
         K=IAT(J)
         JAT(K)=I
       51 IAT(J)=K+1
     ENDIF
     60 CONTINUE
     !
   END SUBROUTINE TransPose1x1
!
!----------------------------------------------------------------------
!
   SUBROUTINE SymbOrder(IA,JA,N,M)
     INTEGER,DIMENSION(:) :: IA,JA
     TYPE(INT_VECT) :: IAT1,JAT1,IAT2,JAT2
     INTEGER :: I,J,K,L,M,N,NZ
     REAL(DOUBLE)   :: Aux1(2),Aux2(2)
     NZ=IA(N+1)-1
     !
     CALL New(IAT1,M+1)
     CALL New(JAT1,NZ)
     CALL TransPose1x1(IA,JA,Aux1,N,M,IAT1%I,JAT1%I,Aux2,'symb')
     CALL New(IAT2,N+1)
     CALL New(JAT2,NZ)
     CALL TransPose1x1(IAT1%I,JAT1%I,Aux1,M,N,IAT2%I,JAT2%I,Aux2,'symb')
     IA(1:N+1)=IAT2%I
     JA(1:NZ)=JAT2%I
     CALL Delete(IAT1)
     CALL Delete(JAT1)
     CALL Delete(IAT2)
     CALL Delete(JAT2)
   END SUBROUTINE SymbOrder
!
!----------------------------------------------------------------------
!
   SUBROUTINE IShow1x1(GcSRowPt,GcSColPt,GcSMTrix,Char,N,M,N1_O,M1_O)
     INTEGER,DIMENSION(:)   :: GcSRowPt,GcSColPt
     INTEGER,DIMENSION(:)   :: GcSMTrix
     INTEGER                :: I,J,N,NDim,K,L,M
     INTEGER,OPTIONAL       :: N1_O,M1_O
     TYPE(INT_RNK2)         :: Aux
     CHARACTER(LEN=*)       :: Char
     !
     WRITE(6,*) Char
     CALL ISp1x1ToFull(GcSRowPt,GcSColPt,GcSMTrix,N,M,Aux)
     DO I=1,N
       WRITE(6,100) (Aux%I(I,J),J=1,M)
     ENDDO
     100 FORMAT(100I5)
     CALL Delete(Aux)
   END SUBROUTINE IShow1x1
!
!----------------------------------------------------------------------
!
   SUBROUTINE Show1x1(GcSRowPt,GcSColPt,GcSMTrix,Char,N,M,N1_O,M1_O)
     INTEGER,DIMENSION(:)   :: GcSRowPt,GcSColPt
     REAL(DOUBLE),DIMENSION(:) :: GcSMTrix
     INTEGER        :: I,J,N,NDim,K,L,M,II,JJ
     INTEGER,OPTIONAL       :: N1_O,M1_O
     TYPE(DBL_RNK2) :: Aux,Aux2
     CHARACTER(LEN=*) :: Char
     !
     CALL Sp1x1ToFull(GcSRowPt,GcSColPt,GcSMTrix,N,M,Aux)
     IF(PRESENT(N1_O).AND.PRESENT(M1_O)) THEN
       CALL New(Aux2,(/N-N1_O+1,M-M1_O+1/))
       II=0
       DO I=N1_O,N
         II=II+1
         JJ=0
         DO J=M1_O,M
           JJ=JJ+1
           Aux2%D(II,JJ)=Aux%D(I,J)
         ENDDO
       ENDDO
       CALL PPrint(Aux2,Char,Unit_O=6)
       CALL Delete(Aux2)
     ELSE
     ! CALL PPrint(Aux,Char,Unit_O=6)
       WRITE(*,*) TRIM(Char)
       DO I=1,N
         WRITE(*,123) (Aux%D(I,J),J=1,M)
       ENDDO
       123 FORMAT(20F6.3)
     ENDIF
     CALL Delete(Aux)
   END SUBROUTINE Show1x1
!
!------------------------------------------------------------------
!
   SUBROUTINE ISp1x1ToFull(IA,JA,AN,N,M,FullMat)
     INTEGER,DIMENSION(:)      :: IA,JA
     INTEGER,DIMENSION(:)      :: AN
     INTEGER                   :: I,J,N,L,M
     TYPE(INT_RNK2)            :: FullMat
     !
     CALL New(FullMat,(/N,M/))
     FullMat%I=Zero
     DO I=1,N
       DO J=IA(I),IA(I+1)-1
         L=JA(J)
         FullMat%I(I,L)=AN(J)
       ENDDO
     ENDDO
   END SUBROUTINE ISp1x1ToFull
!
!------------------------------------------------------------------
!
   SUBROUTINE Sp1x1ToFull(IA,JA,AN,N,M,FullMat)
     INTEGER,DIMENSION(:)   :: IA,JA
     REAL(DOUBLE),DIMENSION(:) :: AN
     INTEGER        :: I,J,N,L,M
     TYPE(DBL_RNK2) :: FullMat
     !
     CALL New(FullMat,(/N,M/))
     FullMat%D=Zero
     DO I=1,N
       DO J=IA(I),IA(I+1)-1
         L=JA(J)
         FullMat%D(I,L)=AN(J)
       ENDDO
     ENDDO
   END SUBROUTINE Sp1x1ToFull
!
!---------------------------------------------------------------------
!
   SUBROUTINE TopToSp1x1(Top,IA,JA)
     INTEGER,DIMENSION(:,:)      :: Top
     TYPE(INT_VECT)              :: IA,JA
     INTEGER                     :: I,J,Dim1,NZ
     !
     Dim1=SIZE(Top,1)
     CALL New(IA,Dim1+1)
     IA%I(1)=1
     DO I=1,Dim1
       IA%I(I+1)=IA%I(I)+Top(I,1)
     ENDDO
     NZ=IA%I(Dim1+1)-1
     CALL New(JA,NZ)
     JA%I=0
     NZ=0
     DO I=1,Dim1
       DO J=1,Top(I,1)
         NZ=NZ+1
         JA%I(NZ)=Top(I,J+1)
       ENDDO
     ENDDO
   END SUBROUTINE TopToSp1x1
!
!---------------------------------------------------------------------
!
   SUBROUTINE IntCToTop(IntCs,Top)
     INTEGER,DIMENSION(:,:) :: Top
     TYPE(INTC)             :: IntCs
     INTEGER                :: I,J
     !
     DO I=1,IntCs%N
       IF(IntCs%Def%C(I)(1:4)=='STRE') THEN
         Top(I,1)=2
       ELSE IF(IntCs%Def%C(I)(1:4)=='BEND') THEN
         Top(I,1)=3
       ELSE IF(IntCs%Def%C(I)(1:4)=='TORS') THEN
         Top(I,1)=4
       ELSE IF(IntCs%Def%C(I)(1:4)=='OUTP') THEN
         Top(I,1)=4
       ELSE IF(IntCs%Def%C(I)(1:4)=='LINB') THEN
         Top(I,1)=3
       ELSE IF(IntCs%Def%C(I)(1:4)=='CART') THEN
         Top(I,1)=1
       ENDIF
       DO J=1,Top(I,1)
         Top(I,J+1)=IntCs%Atoms%I(I,J)
       ENDDO
     ENDDO
   END SUBROUTINE IntCToTop
!
!---------------------------------------------------------------------
!
   SUBROUTINE Sp1x1ToTop(Top,IA,JA)
     TYPE(INT_RNK2)  :: Top
     TYPE(INT_VECT)  :: IA,JA
     INTEGER         :: Dim1,NZ,Max2,I,J
     !
     Dim1=SIZE(IA%I)-1
     Max2=0
     DO I=1,Dim1 ; Max2=MAX(Max2,IA%I(I+1)-IA%I(I)) ; ENDDO
     CALL New(Top,(/Dim1,Max2+1/))
     Top%I=0
     DO I=1,Dim1
       Top%I(I,1)=IA%I(I+1)-IA%I(I)
       NZ=IA%I(I)
       DO J=1,Top%I(I,1)
         Top%I(I,J+1)=JA%I(NZ+J-1)
       ENDDO
     ENDDO
   END SUBROUTINE Sp1x1ToTop
!
!---------------------------------------------------------------------
!
   SUBROUTINE FillInto1x1(Gc,GcSRowPt,GcSColPt,GcSDiag,GcSMTrix, &
       GcSNon0,NCart,Char)
     TYPE(BCSR) :: Gc
     TYPE(INT_VECT) :: GcSRowPt,GcSColPt,AuxI
     TYPE(DBL_VECT) :: GcSMTrix,GcSDiag
     TYPE(DBL_RNK2) :: Aux
     INTEGER :: NCart,I,J,K,L,M,N,JJ,GcSNon0,I1,I2,III,NN,LL,MM,NZ
     CHARACTER(LEN=*) :: Char
     REAL(DOUBLE)   :: Threshold,Sum
     !
     Threshold=1.D-8
     !
       GcSRowPt%I(1)=1
       GcSNon0=0
       GcSDiag%D=Zero
     DO I=1,Gc%Natms
       IF(OffS%I(I)>NCart) EXIT
       I1=Gc%RowPt%I(I)
       I2=Gc%RowPt%I(I+1)
       ! count number of nonzero columns at block-row I
        N=0
       DO J=I1,I2-1
         N=N+BSiz%I(Gc%ColPt%I(J))
       ENDDO
       ! fill non-zeros of row into Aux
       III=BSiz%I(I)
       CALL New(Aux,(/III,N/))
       CALL New(AuxI,N)
           NN=0
       DO J=I1,I2-1
         K=Gc%BlkPt%I(J)
         L=Gc%ColPt%I(J)
         DO JJ=1,BSiz%I(L)
           NN=NN+1
           AuxI%I(NN)=OffS%I(L)+JJ-1
           Aux%D(1:III,NN)=Gc%MTrix%D(K+(JJ-1)*III:)
         ENDDO
       ENDDO
       IF(Char=='lower') THEN
       ! fill aux matrix into 1x1 blocked sparse row, lower triangle only
         DO J=1,III
           LL=OffS%I(I)-1+J
           DO JJ=1,NN
             Sum=Aux%D(J,JJ)
             MM=AuxI%I(JJ)
             IF(MM==LL) GcSDiag%D(LL)=Sum
             IF(ABS(Sum)>Threshold) THEN
               IF(MM<LL) THEN
                 GcSNon0=GcSNon0+1
                 GcSColPt%I(GcSNon0)=MM
                 GcSMTrix%D(GcSNon0)=Sum
               ENDIF
             ENDIF
           ENDDO
           GcSRowPt%I(LL+1)=GcSNon0+1
         ENDDO
       ELSE IF(Char=='upper') THEN
       ! fill aux matrix into 1x1 blocked sparse row, upper triangle only
         DO J=1,III
           LL=OffS%I(I)-1+J
           DO JJ=1,NN
             Sum=Aux%D(J,JJ)
             MM=AuxI%I(JJ)
             IF(MM==LL) GcSDiag%D(LL)=Sum
             IF(ABS(Sum)>Threshold) THEN
               IF(MM>LL) THEN
                 GcSNon0=GcSNon0+1
                 GcSColPt%I(GcSNon0)=MM
                 GcSMTrix%D(GcSNon0)=Sum
               ENDIF
             ENDIF
           ENDDO
           GcSRowPt%I(LL+1)=GcSNon0+1
         ENDDO
       ELSE IF(Char=='full') THEN
         ! fill aux matrix into 1x1 blocked sparse row, full mtrx
         DO J=1,III
           LL=OffS%I(I)-1+J
           DO JJ=1,NN
             Sum=Aux%D(J,JJ)
             MM=AuxI%I(JJ)
             IF(MM==LL) GcSDiag%D(LL)=Sum
             IF(ABS(Sum)>Threshold.OR.MM==LL) THEN
                 GcSNon0=GcSNon0+1
                 GcSColPt%I(GcSNon0)=MM
                 GcSMTrix%D(GcSNon0)=Sum
             ENDIF
           ENDDO
           GcSRowPt%I(LL+1)=GcSNon0+1
         ENDDO
       ENDIF
       CALL Delete(AuxI)
       CALL Delete(Aux)
     ENDDO
   END SUBROUTINE FillInto1x1
!
!--------------------------------------------------------------------
!
   SUBROUTINE UpperTr(IA,JA,AN,IA2,JA2,AN2,AD2)
     INTEGER,DIMENSION(:) :: IA,JA,IA2,JA2
     REAL(DOUBLE),DIMENSION(:) :: AN,AN2,AD2
     INTEGER :: I,J,K,L,N,NZ
     N=SIZE(IA)-1
     NZ=0
     AD2=Zero
     IA2(1)=1
     DO I=1,N
       DO J=IA(I),IA(I+1)-1
         K=JA(J)
         IF(K==I) AD2(I)=AN(J)
         IF(K<=I) CYCLE
         NZ=NZ+1
         JA2(NZ)=K
         AN2(NZ)=AN(J)
       ENDDO
       IA2(I+1)=NZ+1
     ENDDO
   END SUBROUTINE UpperTr
!
!--------------------------------------------------------------------
!
   SUBROUTINE Scale1x1(Scale,IA,JA,AN,AD,NCart)
     INTEGER,DIMENSION(:) :: IA,JA
     REAL(DOUBLE),DIMENSION(:) :: AN,AD,Scale
     REAL(DOUBLE) :: FactI,FactK
     INTEGER :: NCart,I,J,K,L,N,M
     DO I=1,NCart
       FactI=Scale(I)
       AD(I)=FactI*AD(I)*FactI
       DO J=IA(I),IA(I+1)-1
         K=JA(J)
         FactK=Scale(K)
         AN(J)=FactI*FactK*AN(J)
       ENDDO
     ENDDO
   END SUBROUTINE Scale1x1
!
!--------------------------------------------------------------------
!
   SUBROUTINE CholFactSolve(IU,JU,UN,DI,B,N,X)
     INTEGER,DIMENSION(:) :: IU,JU
     REAL(DOUBLE),DIMENSION(:) :: UN,DI,B,X
     INTEGER :: N,I,J,K,NM,IUA,IUB
     REAL(DOUBLE) :: XX
     !
     ! Solve Ut*D*U X = B
     !
     NM=N-1
     X(1:N)=B(1:N)
     DO 40 K=1,NM
       IUA=IU(K)
       IUB=IU(K+1)-1
       XX=X(K)
       IF(IUB.LT.IUA) GO TO 30
       DO 20 I=IUA,IUB
       20 X(JU(I))=X(JU(I))-UN(I)*XX
       30 X(K)=XX*DI(K)
     40 CONTINUE
     X(N)=X(N)*DI(N)
     K=NM
     IF(K<=0) RETURN
     50 IUA=IU(K)
     IUB=IU(K+1)-1
     IF(IUB.LT.IUA) GO TO 70
     XX=X(K)
     DO 60 I=IUA,IUB
     60 XX=XX-UN(I)*X(JU(I))
     X(K)=XX
     70 K=K-1
     IF(K.GT.0) GO TO 50
   END SUBROUTINE CholFactSolve
!
!-----------------------------------------------------------------
!
   SUBROUTINE Symmetrize1x1(IA,JA,A_O)
     TYPE(INT_VECT)                     :: IA,JA
     TYPE(DBL_VECT),OPTIONAL            :: A_O
     REAL(DOUBLE)                       :: Aux1(2),Aux2(2)
     TYPE(INT_VECT)                     :: ITr,JTr
     TYPE(DBL_VECT)                     :: ATr
     TYPE(INT_VECT)                     :: IC,JC
     TYPE(DBL_VECT)                     :: CN
     INTEGER                            :: N,NZ
     !
     N=SIZE(IA%I)-1
     NZ=SIZE(JA%I)
     CALL New(ITr,N+1)
     CALL New(JTr,NZ)
     IF(PRESENT(A_O)) THEN
       CALL New(ATr,NZ)
       CALL TransPose1x1(IA%I,JA%I,A_O%D,N,N,ITr%I,JTr%I,ATr%D,'full')
       CALL AddMat_1x1(IA%I,JA%I,A_O%D,ITr%I,JTr%I,ATr%D,IC,JC,CN,N,N,SymbOnly_O=.FALSE.)
     ELSE
       CALL TransPose1x1(IA%I,JA%I,Aux1,N,N,ITr%I,JTr%I,Aux2,'symb')
       CALL AddMat_1x1(IA%I,JA%I,Aux1,ITr%I,JTr%I,Aux2,IC,JC,CN,N,N,SymbOnly_O=.TRUE.)
     ENDIF
     CALL Delete(IA)
     CALL Delete(JA)
     CALL New(IA,SIZE(IC%I))
     CALL New(JA,SIZE(JC%I))
     CALL SetEq(IA,IC)
     CALL SetEq(JA,JC)
     IF(PRESENT(A_O)) THEN
       CALL Delete(A_O)
       CALL New(A_O,SIZE(CN%D))
       CALL SetEq(A_O,CN)
     ENDIF
     CALL Delete(ITr)
     CALL Delete(JTr)
     CALL Delete(IC)
     CALL Delete(JC)
     IF(PRESENT(A_O)) THEN
       CALL Delete(A_O)
       CALL Delete(ATr)
       CALL Delete(CN)
     ENDIF
   END SUBROUTINE Symmetrize1x1
!
!-----------------------------------------------------------------
!
   SUBROUTINE Plot_1x1(IA,JA,Name,NatmsLoc)
      INTEGER,DIMENSION(:) :: IA,JA
      CHARACTER(LEN=*)     :: Name
      INTEGER              :: I,K
      INTEGER              :: NatmsLoc
      !
      CALL OpenASCII(TRIM(Name)//'_PlotFile_1',Plt,NewFile_O=.TRUE.)
      DO I=1,NAtmsLoc
        DO K=IA(I),IA(I+1)-1
          WRITE(Plt,1)JA(K),NAtmsLoc-I
        ENDDO
      ENDDO
      CLOSE(Plt)

      CALL OpenASCII(TRIM(Name)//'_GnuPlotMe',Plt,NewFile_O=.TRUE.)
      WRITE(Plt,2)
      WRITE(Plt,3)TRIM(Name)//'.eps'
      WRITE(Plt,4)DBLE(NAtmsLoc+1)
      WRITE(Plt,5)DBLE(NAtmsLoc)
      WRITE(Plt,6); WRITE(Plt,7); WRITE(Plt,8)
      WRITE(Plt,*)'set pointsize '//TRIM(FltToChar(50.D0/DBLE(NAtmsLoc)))
      WRITE(Plt,*)"plot '"//TRIM(Name)//"_PlotFile_1' using 1:2 notitle with points 1 "
      CLOSE(Plt)
     1   FORMAT(2(1x,I16))
     2   FORMAT('set term  postscript eps  "Times-Roman" 18')
     !!    2   FORMAT('set term  jpeg transparent')
     3   FORMAT('set output "',A,'"')
     4   FORMAT('set xrange [0:',F12.3,']')
     5   FORMAT('set yrange [-1:',F12.3,']')
     6   FORMAT('set size square ')
     7   FORMAT('set noxtics ')
     8   FORMAT('set noytics ')
     9   FORMAT('plot [0 : ',I12,' ] ',I12,', \\')
     10  FORMAT('                    ',I12,', \\')
   END SUBROUTINE Plot_1x1
!
!--------------------------------------------------------------------
!
   SUBROUTINE RCMOrder(Perm,Iperm,NDim,IA,JA)
     INTEGER,DIMENSION(:) :: Perm,Iperm,IA,JA
     INTEGER :: I,J,K,L,NDim
     TYPE(INT_VECT) :: XLs
     !
     CALL New(XLs,NDim+1)
       CALL GenRCM(NDim,IA,JA,Perm,Xls%I)
     CALL Delete(XLs)
     DO I=1,NDim
       IPerm(Perm(I))=I
     ENDDO
   END SUBROUTINE RCMOrder
!
!--------------------------------------------------------------------
!
   SUBROUTINE NoOrder(Perm,IPerm,NDim)
     INTEGER,DIMENSION(:) :: Perm,IPerm
     INTEGER              :: I,NDim
     DO I=1,NDim
       Perm(I)=I
       IPerm(I)=I
     ENDDO
   END SUBROUTINE NoOrder
!
!--------------------------------------------------------------------
!
   SUBROUTINE BMatrFact(B,NCart,IA,JA,AN,AD,&
       Perm,IPerm,GcScale,Char)
     TYPE(BMATR) :: B
     INTEGER     :: NCart,I,J,K,L,NIntC,NDim,NZ,JJ,KK,LL
     INTEGER     :: ColOff,From,To
     TYPE(INT_VECT) :: IA,JA,IA2,JA2
     TYPE(DBL_VECT) :: AN,AD,AN2
     TYPE(INT_VECT) :: Perm,IPerm
     TYPE(DBL_VECT) :: GcScale
     CHARACTER(LEN=*):: Char
     !
     NIntC=SIZE(B%IB%I,1)
     NDim=NIntC+NCart
     CALL New(IA2,NDim+1)
     CALL New(JA2,12*NIntC+NDim)
     CALL New(AN2,12*NIntC+NDim)
     CALL New(AD,NDim)
     !
     ColOff=NIntC
     IF(Char=='Gc') THEN
       ColOff=NIntC
       From=0
     ENDIF
     IF(Char=='Gi') THEN
       ColOff=0
       From=NCart
     ENDIF
       To  =From+NIntC
     !
     NZ=0
     IA2%I(1)=1
     DO I=1,From
       IA2%I(I+1)=NZ+1
     ENDDO
     DO I=1,NIntC
       DO J=1,4
         K=B%IB%I(I,J)
         IF(K==0) CYCLE
         LL=3*(J-1)
         KK=ColOff+3*(K-1)
         DO JJ=1,3
           NZ=NZ+1
           JA2%I(NZ)=KK+JJ
           AN2%D(NZ)=B%B%D(I,LL+JJ)
         ENDDO
       ENDDO
       IA2%I(From+I+1)=NZ+1
     ENDDO
     DO I=To+1,NDim
       IA2%I(I+1)=NZ+1
     ENDDO
     !
     CALL New(IA,NDim+1)
     CALL New(JA,NZ)
     CALL New(AN,NZ)
     IF(Char=='Gi') THEN
       CALL TransPose1x1(IA2%I,JA2%I,AN2%D,NDim,NDim,&
                         IA%I,JA%I,AN%D,'full')
     ELSE
       IA%I(1:NDim+1)=IA2%I(1:NDim+1)
       JA%I(1:NZ)=JA2%I(1:NZ)
       AN%D(1:NZ)=AN2%D(1:NZ)
     ENDIF
     CALL Delete(IA2)
     CALL Delete(JA2)
     CALL Delete(AN2)
     !
     ! Prepare factor with unit diagonals
     !
     IF(Char=='Gc') THEN
       AD%D(1:NintC)=1.D-6
       AD%D(NintC+1:NDim)=1.D-3
     ENDIF
     IF(Char=='Gi') THEN
       AD%D(1:NCart)=1.D-3
       AD%D(NCart+1:NDim)=1.D-6
     ENDIF
     !
     CALL New(Perm,NDim)
     CALL New(IPerm,NDim)
     CALL New(GcScale,NDim)
     DO I=1,NDim
       Perm%I(I)=I
       IPerm%I(I)=I
       GcScale%D(I)=One
     ENDDO
   END SUBROUTINE BMatrFact
!
!--------------------------------------------------------------------
!
   SUBROUTINE UInvX(IU,JU,UN,X,N,W)
     INTEGER,DIMENSION(:) :: IU,JU
     INTEGER              :: N,I,J,K,L,IUA,IUB,IP
     REAL(DOUBLE),DIMENSION(:) :: UN,X,W
     REAL(DOUBLE)              :: Z
     !
     W(1:N)=X(1:N)
     I=N-1
     20 IUA=IU(I)
     IUB=IU(I+1)-1
     IF(IUB.LT.IUA) GO TO 40
     Z=W(I)
     DO 30 IP=IUA,IUB
     30 Z=Z-UN(IP)*W(JU(IP))
     W(I)=Z
     40 I=I-1
     IF(I.GT.0) GO TO 20
   END SUBROUTINE UInvX
!
!--------------------------------------------------------------------
!
   SUBROUTINE UtInvX(IU,JU,UN,X,N,W)
     INTEGER,DIMENSION(:) :: IU,JU
     INTEGER              :: N,I,J,K,L,IUA,IUB,IP,NM
     REAL(DOUBLE),DIMENSION(:) :: UN,X,W
     REAL(DOUBLE)              :: Z
     !
     NM=N-1
     W(1:N)=X(1:N)
     DO 30 I=1,NM
       IUA=IU(I)
       IUB=IU(I+1)-1
       IF(IUB<IUA) GO TO 30
       Z=W(I)
       DO 20 IP=IUA,IUB
         J=JU(IP)
     20 W(J)=W(J)-Z*(UN(IP))
     30 CONTINUE
   END SUBROUTINE UtInvX
!
!--------------------------------------------------------------------
!
   SUBROUTINE ProdMatrVect(IA,JA,AN,B,N,C)
     INTEGER,DIMENSION(:)      :: IA,JA
     REAL(DOUBLE),DIMENSION(:) :: AN,B,C
     REAL(DOUBLE)              :: U
     INTEGER                   :: N,I,J,IAA,IAB,K
     ! C=A*B
     DO 20 I=1,N
       U=Zero
       IAA=IA(I)
       IAB=IA(I+1)-1
       IF(IAB<IAA) GO TO 20
       DO 10 K=IAA,IAB
       10 U=U+AN(K)*B(JA(K))
     20 C(I)=U
   END SUBROUTINE ProdMatrVect
!
!--------------------------------------------------------------------
!
   SUBROUTINE BackSubst(IU,JU,UN,UD,VectIn,VectOut)
     INTEGER,DIMENSION(:) :: IU,JU
     REAL(DOUBLE),DIMENSION(:) :: UN,UD,VectIn,VectOut
     INTEGER :: I,J,K,L,N
     REAL(DOUBLE) :: Sum,UDI
     !
     N=SIZE(IU)-1
     DO I=1,N
       UDI=UD(I)
       Sum=VectIn(I)
       DO J=IU(I),IU(I+1)-1
         Sum=Sum-UN(J)*VectIn(JU(J))
       ENDDO
       VectOut(I)=Sum/UDI
     ENDDO
   END SUBROUTINE BackSubst
!
!--------------------------------------------------------------------
!
   SUBROUTINE GetGc(NCart,ISpB,JSpB,ASpB,IGc,JGc,AGc,SymbOnly_O)
     TYPE(INT_VECT)  :: ISpB,JSpB
     TYPE(DBL_VECT)  :: ASpB
     TYPE(INT_VECT)  :: ISpBt,JSpBt
     TYPE(DBL_VECT)  :: ASpBt
     TYPE(INT_VECT)  :: IGc,JGc
     TYPE(DBL_VECT)  :: AGc
     TYPE(INT_VECT)  :: JGc2,IX
     TYPE(DBL_VECT)  :: AGc2,X
     INTEGER         :: I,J,K,L
     INTEGER         :: NZSpB,NZGc,NIntC,NCart
     LOGICAL,OPTIONAL:: SymbOnly_O
     !
     NIntC=SIZE(ISpB%I)-1
     NZSpB=ISpB%I(NIntC+1)-1
     !
     ! estimate NZGc
     !
     NZGc=0
     DO I=1,NIntC
       NZGc=NZGc+(ISpB%I(I+1)-ISpB%I(I))**2
     ENDDO
     !
     ! Transpose SpB
     !
     CALL New(ISpBt,NCart+1)
     CALL New(JSpBt,NZSpB)
     CALL New(ASpBt,NZSpB)
     CALL TransPose1x1(ISpB%I,JSpB%I,ASpB%D,NIntC,NCart, &
          ISpBt%I,JSpBt%I,ASpBt%D,'full')
     !
     CALL MatMul_1x1(ISpBt%I,JSpBt%I,ASpBt%D,ISpB%I,JSpB%I,ASpB%D, &
                     IGc,JGc,AGc,NCart,NIntC,NCart,SymbOnly_O=SymbOnly_O)
     IF(PRESENT(SymbOnly_O)) THEN
       IF(SymbOnly_O) THEN
         CALL Delete(ISpBt)
         CALL Delete(JSpBt)
         CALL Delete(ASpBt)
         RETURN
       ENDIF
     ELSE
       CALL ThreshMatr(IGc,JGc,AGc,1.D-7)
       CALL Delete(ISpBt)
       CALL Delete(JSpBt)
       CALL Delete(ASpBt)
     ENDIF
   END SUBROUTINE GetGc
!
!--------------------------------------------------------------------
!
   SUBROUTINE BtoSpB_1x1(B,ISpB,JSpB,ASpB)
     TYPE(BMATR)     :: B
     TYPE(INT_VECT)  :: ISpB,JSpB
     TYPE(DBL_VECT)  :: ASpB
     TYPE(INT_VECT)  :: JSpB2
     TYPE(DBL_VECT)  :: ASpB2
     INTEGER         :: I,J,K,L,JJ,KK,NCart
     INTEGER         :: NIntC,NZSpB
     !
     NIntC=SIZE(B%IB%I,1)
     NZSpB=(12+9)*NIntC
     CALL New(ISpB,NIntC+1)
     CALL New(JSpB2,NZSpB)
     CALL New(ASpB2,NZSpB)
     NZSpB=0
     ISpB%I(1)=1
     DO I=1,NIntC
       DO J=1,4
         K=B%IB%I(I,J)
         IF(K==0) EXIT
         JJ=3*(J-1)
         KK=3*(K-1)
         DO L=1,3
           NZSpB=NZSpB+1
           KK=KK+1
           JJ=JJ+1
           JSpB2%I(NZSpB)=KK
           ASpB2%D(NZSpB)=B%B%D(I,JJ)
         ENDDO
       ENDDO
       !
       ! Copy lattice part (B%BL)
       !
       IF(B%BLI%I(I)/=0) THEN
         NCart=B%BLI%I(I)
         DO J=1,9
           NZSpB=NZSpB+1
           JSpB2%I(NZSpB)=NCart+J
           ASpB2%D(NZSpB)=B%BL%D(I,J)
         ENDDO
       ENDIF
       ISpB%I(I+1)=NZSpB+1
     ENDDO
     CALL New(JSpB,NZSpB)
     CALL New(ASpB,NZSpB)
       JSpB%I(1:NZSpB)=JSpB2%I(1:NZSpB)
       ASpB%D(1:NZSpB)=ASpB2%D(1:NZSpB)
     CALL Delete(JSpB2)
     CALL Delete(ASpB2)
     !
     CALL ThreshMatr(ISpB,JSpB,ASpB,1.D-8)
   END SUBROUTINE BtoSpB_1x1
!
!--------------------------------------------------------------------
!
   SUBROUTINE IMatMul_1x1(IA,JA,AN,IB,JB,BN,IC,JC,CN,NP,NQ,NR,SymbOnly_O)
     INTEGER,DIMENSION(:) :: IA,JA,IB,JB
     INTEGER,DIMENSION(:) :: AN,BN
     TYPE(INT_VECT)       :: IC,JC
     TYPE(INT_VECT)       :: CN
     INTEGER              :: I,J,K,L,NP,NQ,NR,NZC
     LOGICAL,OPTIONAL     :: SymbOnly_O
     ! C=A*B
     ! A : NP x NQ matrix
     ! B : NQ x NR matrix
     ! C : NP x NR matrix
     !
     CALL MatMulSymbDriver(IA,JA,IB,JB,NP,NQ,NR,IC,JC)
     IF(PRESENT(SymbOnly_O)) THEN
       IF(SymbOnly_O) RETURN
     ENDIF
     NZC=IC%I(NP+1)-1
     CALL New(CN,NZC)
     CALL IMatMulNum(IA,JA,AN,IB,JB,BN,IC%I,JC%I,NP,NR,CN%I)
   END SUBROUTINE IMatMul_1x1
!
!--------------------------------------------------------------------
!
   SUBROUTINE MatMul_1x1(IA,JA,AN,IB,JB,BN,IC,JC,CN,NP,NQ,NR,SymbOnly_O)
     INTEGER,DIMENSION(:) :: IA,JA,IB,JB
     REAL(DOUBLE),DIMENSION(:) :: AN,BN
     TYPE(INT_VECT)   :: IC,JC
     TYPE(DBL_VECT)   :: CN
     INTEGER          :: I,J,K,L,NP,NQ,NR,NZC
     LOGICAL,OPTIONAL :: SymbOnly_O
     ! C=A*B
     ! A : NP x NQ matrix
     ! B : NQ x NR matrix
     ! C : NP x NR matrix
     !
     CALL MatMulSymbDriver(IA,JA,IB,JB,NP,NQ,NR,IC,JC)
     IF(PRESENT(SymbOnly_O)) THEN
       IF(SymbOnly_O) RETURN
     ENDIF
     NZC=IC%I(NP+1)-1
     CALL New(CN,NZC)
     CALL MatMulNum(IA,JA,AN,IB,JB,BN,IC%I,JC%I,NP,NR,CN%D)
   END SUBROUTINE MatMul_1x1
!
!--------------------------------------------------------------------
!
   SUBROUTINE MatMulSymbDriver(IA,JA,IB,JB,NP,NQ,NR,IC,JC)
     TYPE(INT_VECT)       :: IC,JC,JC2,ColSum
     INTEGER,DIMENSION(:) :: IA,JA,IB,JB
     INTEGER              :: NP,NQ,NR,NZC,I,J,K
     !
     ! first, calc. upper-limit for filing of matrix C
     CALL New(ColSum,NQ)
     ColSum%I=0
     DO I=1,NP
       DO J=IA(I),IA(I+1)-1
         K=JA(J)
         ColSum%I(K)=ColSum%I(K)+1
       ENDDO
     ENDDO
     NZC=0
     DO I=1,NQ
       NZC=NZC+ColSum%I(I)*(IB(I+1)-IB(I))
     ENDDO
     CALL Delete(ColSum)
     CALL New(IC,NP+1)
     CALL New(JC2,NZC)
     CALL MatMulSymb(IA,JA,IB,JB,NP,NQ,NR,IC%I,JC2%I)
     NZC=IC%I(NP+1)-1
     CALL New(JC,NZC)
     JC%I(1:NZC)=JC2%I(1:NZC)
     CALL Delete(JC2)
   END SUBROUTINE MatMulSymbDriver
!
!--------------------------------------------------------------------
!
   SUBROUTINE MatMulSymb(IA,JA,IB,JB,NP,NQ,NR,IC,JC)
     !
     ! A : NP x NQ matrix
     ! B : NQ x NR matrix
     ! C : NP x NR matrix
     ! IX: scratch
     !
     ! This routine carries out symbolic multiplication of C = A*B
     !
     INTEGER      :: NP,NQ,NR
     INTEGER      :: IA(NP+1),JA(*)
     INTEGER      :: IB(NQ+1),JB(*)
     INTEGER      :: IC(NP+1),JC(*)
     INTEGER      :: I,J,K,L
     INTEGER      :: IP,IAA,IAB,JP,IBA,IBB,KP
     TYPE(INT_VECT):: IX
     !
     CALL New(IX,NR)
     IP=1
     IX%I(:)=0
     DO 40 I=1,NP
     IC(I)=IP
     IAA=IA(I)
     IAB=IA(I+1)-1
     IF(IAB<IAA) GO TO 40
     DO 30 JP=IAA,IAB
     J=JA(JP)
     IBA=IB(J)
     IBB=IB(J+1)-1
     IF(IBB<IBA) GO TO 30
     DO 20 KP=IBA,IBB
     K=JB(KP)
     IF(IX%I(K)==I) GO TO 20
     JC(IP)=K
     IP=IP+1
     IX%I(K)=I
20   CONTINUE
30   CONTINUE
40   CONTINUE
     IC(NP+1)=IP
     CALL Delete(IX)
   END SUBROUTINE MatMulSymb
!
!--------------------------------------------------------------------
!
   SUBROUTINE IMatMulNum(IA,JA,AN,IB,JB,BN,IC,JC,NP,NR,CN)
     !
     ! Numerical multiplication of two matrices: C = A*B
     ! all numbers are integer
     !
     ! NP: the number of rows of the first and result matrices
     ! NR: the number of coulumns of the result matrix
     ! from matrix size data
     ! X: scratch array
     !
     INTEGER,DIMENSION(:)      :: IA,JA,IB,JB,IC,JC
     INTEGER,DIMENSION(:)      :: AN,BN,CN
     INTEGER                   :: I,J,K,L
     INTEGER                   :: NP,NR,ICA,ICB,IAA,IAB,JP,IBA,IBB,KP
     INTEGER                   :: A
     TYPE(INT_VECT)            :: X
     !
     CALL New(X,NR)
     DO 50 I=1,NP
     ICA=IC(I)
     ICB=IC(I+1)-1
     IF(ICB.LT.ICA) GO TO 50
     DO 10 J=ICA,ICB
10   X%I(JC(J))=Zero
     IAA=IA(I)
     IAB=IA(I+1)-1
     DO 30 JP=IAA,IAB
     J=JA(JP)
     A=AN(JP)
     IBA=IB(J)
     IBB=IB(J+1)-1
     IF(IBB<IBA) GO TO 30
     DO 20 KP=IBA,IBB
     K=JB(KP)
     X%I(K)=X%I(K)+A*BN(KP)
20   CONTINUE
30   CONTINUE
     DO 40 J=ICA,ICB
40   CN(J)=X%I(JC(J))
50   CONTINUE
     CALL Delete(X)
   END SUBROUTINE IMatMulNum
!
!--------------------------------------------------------------------
!
   SUBROUTINE MatMulNum(IA,JA,AN,IB,JB,BN,IC,JC,NP,NR,CN)
     !
     ! Numerical multiplication of two matrices: C = A*B
     !
     ! NP: the number of rows of the first and result matrices
     ! NR: the number of coulumns of the result matrix
     ! from matrix size data
     ! X: scratch array
     !
     INTEGER,DIMENSION(:)      :: IA,JA,IB,JB,IC,JC
     REAL(DOUBLE),DIMENSION(:) :: AN,BN,CN
     INTEGER                   :: I,J,K,L
     INTEGER                   :: NP,NR,ICA,ICB,IAA,IAB,JP,IBA,IBB,KP
     REAL(DOUBLE)              :: A
     TYPE(DBL_VECT)            :: X
     !
     CALL New(X,NR)
     DO 50 I=1,NP
     ICA=IC(I)
     ICB=IC(I+1)-1
     IF(ICB.LT.ICA) GO TO 50
     DO 10 J=ICA,ICB
10   X%D(JC(J))=Zero
     IAA=IA(I)
     IAB=IA(I+1)-1
     DO 30 JP=IAA,IAB
     J=JA(JP)
     A=AN(JP)
     IBA=IB(J)
     IBB=IB(J+1)-1
     IF(IBB<IBA) GO TO 30
     DO 20 KP=IBA,IBB
     K=JB(KP)
     X%D(K)=X%D(K)+A*BN(KP)
20   CONTINUE
30   CONTINUE
     DO 40 J=ICA,ICB
40   CN(J)=X%D(JC(J))
50   CONTINUE
     CALL Delete(X)
   END SUBROUTINE MatMulNum
!
!--------------------------------------------------------------------
!
   SUBROUTINE CholXMatrXChol(CholData,IB,JB,BN,N,IX,JX,XN)
     !
     ! X= (D^-1/2)*(U^-T)*B*(U^-1)*(D^-1/2) (E.g. for preconditioning)
     ! B,X: NxN matrix
     ! UT : NxN upper triangle, with unit diagonals (diags not stored)!
     !
     TYPE(Cholesky) :: CholData
     TYPE(INT_VECT) :: IB,JB,IX,JX,IX2,JX2
     TYPE(DBL_VECT) :: BN,XN,X2N
     INTEGER        :: N,NC,I,J
     TYPE(DBL_VECT) :: Diag12
     !
     CALL Halt('This routine is under development, effect scaling matrices CholData%GcScale%D should be added')
     NC=N
     CALL New(Diag12,N)
     DO I=1,N ; Diag12%D(I)=SQRT(ABS(CholData%ChDiag%D(I))) ; ENDDO
     !
     CALL Perm1x1(CholData%Perm%I,IB%I,JB%I,BN%D)
     !
     CALL MatrXUInv(CholData%ChRowPt%I,CholData%ChColPt%I, &
                    CholData%ChFact%D, &
                    IB%I,JB%I,BN%D,IX2,JX2,X2N,N,NC)
     !
     CALL MatrXD(IX2%I,JX2%I,X2N%D,Diag12%D,N,NC)
     !
     CALL UTInvXMatr(CholData%ChRowPt%I,CholData%ChColPt%I, &
                     CholData%ChFact%D, &
                     IX2%I,JX2%I,X2N%D,IX,JX,XN,N,NC)
     CALL DxMatr(Diag12%D,IX%I,JX%I,XN%D,N,NC)
     !
     CALL Delete(IX2) ; CALL Delete(JX2) ; CALL Delete(X2N)
     CALL Delete(Diag12)
     !
     CALL Perm1x1(CholData%IPerm%I,IX%I,JX%I,XN%D)
     !
   END SUBROUTINE CholXMatrXChol
!
!--------------------------------------------------------------------
!
   SUBROUTINE InvMatXMatr(CholData,IB,JB,BN,IX,JX,XN,N,NC,Thresh)
     !
     ! X=(U^-1* D^-1 * U^-T) B
     ! B,X: NxNC matrix
     ! U : NxN upper triangle, with unit diagonals (diags not stored)!
     !
     TYPE(Cholesky)            :: CholData
     INTEGER,DIMENSION(:)      :: IB,JB
     REAL(DOUBLE),DIMENSION(:) :: BN
     TYPE(INT_VECT)            :: IX,JX,IX1,JX1,IBc,JBc
     TYPE(DBL_VECT)            :: XN,XN1,BNc
     INTEGER                   :: N,NC,NZ
     REAL(DOUBLE)              :: Thresh
     !
     NZ=IB(N+1)-1
     CALL New(IBc,N+1)
     CALL New(JBc,NZ)
     CALL New(BNc,NZ)
     IBc%I(1:N+1)=IB(1:N+1) ; JBc%I(1:NZ)=JB(1:NZ) ; BNc%D(1:NZ)=BN(1:NZ)
     !
     CALL PermRow(IBc%I,JBc%I,BNc%D,CholData%Perm%I,N,NC)
     CALL ScaleRow(IBc%I,JBc%I,BNc%D,CholData%GcScale%D)
     CALL UTInvXMatr(CholData%ChRowPt%I,CholData%ChColPt%I, &
                     CholData%ChFact%D, &
                     IBc%I,JBc%I,BNc%D,IX1,JX1,XN1,N,NC)
     CALL Delete(IBc)
     CALL Delete(JBc)
     CALL Delete(BNc)
     !
     CALL DxMatr(CholData%ChDiag%D,IX1%I,JX1%I,XN1%D,N,NC)
     !
     CALL UInvXMatr(CholData%ChRowPt%I,CholData%ChColPt%I, &
                    CholData%ChFact%D, &
                    IX1%I,JX1%I,XN1%D,IX,JX,XN,N,NC)
     CALL Delete(IX1)
     CALL Delete(JX1)
     CALL Delete(XN1)
     CALL ScaleRow(IX%I,JX%I,XN%D,CholData%GcScale%D)
     CALL PermRow(IX%I,JX%I,XN%D,CholData%IPerm%I,N,NC)
     CALL ThreshMatr(IX,JX,XN,Thresh)
     !
   END SUBROUTINE InvMatXMatr
!
!--------------------------------------------------------------------
!
   SUBROUTINE MatrXUInv(IU,JU,UN,IB,JB,BN,IX,JX,XN,N,NC)
     !
     ! X=B * U^-1
     ! B,X: NCxN matrix
     ! UT : NxN upper triangle, with unit diagonals (diags not stored)!
     !
     INTEGER,DIMENSION(:)      :: IU,JU,IB,JB
     REAL(DOUBLE),DIMENSION(:) :: UN,BN
     TYPE(INT_VECT)            :: IX,JX,IBT,JBT,IXT,JXT
     TYPE(DBL_VECT)            :: XN,BNT,XNT
     INTEGER                   :: N,NC,NZU,NZB,NZX
     !
     NZB=SIZE(JB)
     CALL New(IBT,N+1)
     CALL New(JBT,NZB)
     CALL New(BNT,NZB)
     CALL TransPose1x1(IB,JB,BN,NC,N,IBT%I,JBT%I,BNT%D,'full')
     !
     CALL UTInvXMatr(IU,JU,UN,IBT%I,JBT%I,BNT%D,IXT,JXT,XNT,N,NC)
     CALL Delete(IBT)
     CALL Delete(JBT)
     CALL Delete(BNT)
     !
     NZX=SIZE(JXT%I)
     CALL New(IX,NC+1)
     CALL New(JX,NZX)
     CALL New(XN,NZX)
     CALL TransPose1x1(IXT%I,JXT%I,XNT%D,N,NC,IX%I,JX%I,XN%D,'full')
     !
     CALL Delete(IXT)
     CALL Delete(JXT)
     CALL Delete(XNT)
   END SUBROUTINE MatrXUInv
!
!--------------------------------------------------------------------
!
   SUBROUTINE UInvXMatr(IU,JU,UN,IB,JB,BN,IX,JX,XN,N,NC)
     INTEGER,DIMENSION(:)      :: IU,JU,IB,JB
     REAL(DOUBLE),DIMENSION(:) :: UN,BN
     TYPE(INT_VECT)            :: IX,JX,IX2,JX2,IUT,JUT
     TYPE(DBL_VECT)            :: XN,UNT
     INTEGER                   :: N,NC,NZU,NZB,NZX
     !
     ! X=UT^-1 * B
     ! B,X: NxNC matrix
     ! U : NxN upper triangle, with unit diagonals (diags not stored)!
     !
     NZX=((IB(N+1)-1)+(IU(N+1)-1))*4 !!! this is an estimate of the numb. of nonzeros in X
     CALL New(IX,N+1)
     CALL New(JX2,NZX)
     CALL UInvXMatrSymb(IU,JU,IB,JB,N,NC,IX%I,JX2%I)
     !
     NZX=IX%I(N+1)-1
     CALL New(JX,NZX)
     CALL New(XN,NZX)
     JX%I(1:NZX)=JX2%I(1:NZX)
     CALL Delete(JX2)
     !
     CALL UInvXMatrNum(IU,JU,UN,IB,JB,BN,IX%I,JX%I,N,NC,XN%D)
     !
   END SUBROUTINE UInvXMatr
!
!--------------------------------------------------------------------
!
   SUBROUTINE UTInvXMatr(IU,JU,UN,IB,JB,BN,IX,JX,XN,N,NC)
     INTEGER,DIMENSION(:)      :: IU,JU,IB,JB
     REAL(DOUBLE),DIMENSION(:) :: UN,BN
     TYPE(INT_VECT)            :: IX,JX,IX2,JX2,IUT,JUT
     TYPE(DBL_VECT)            :: XN,UNT
     INTEGER                   :: N,NC,NZU,NZB,NZX
     !
     ! X=UT^-1 * B
     ! B,X: NxNC matrix
     ! UT : NxN upper triangle, with unit diagonals (diags not stored)!
     !
     NZU=SIZE(JU)
     NZB=SIZE(JB)
     !
     CALL New(IUT,N+1)
     CALL New(JUT,NZU)
     CALL New(UNT,NZU)
     !
     CALL TransPose1x1(IU,JU,UN,N,N,IUT%I,JUT%I,UNT%D,'full')
     !
     NZX=((IB(N+1)-1)+(IU(N+1)-1))*4 !!! this is an estimate of the numb. of nonzeros in X
     CALL New(IX,N+1)
     CALL New(JX2,NZX)
     CALL UTInvXMatrSymb(IUT%I,JUT%I,IB,JB,N,NC,IX%I,JX2%I)
     !
     NZX=IX%I(N+1)-1
     CALL New(JX,NZX)
     CALL New(XN,NZX)
     JX%I(1:NZX)=JX2%I(1:NZX)
     CALL Delete(JX2)
     !
     CALL UTInvXMatrNum(IUT%I,JUT%I,UNT%D,IB,JB,BN,IX%I,JX%I,N,NC,XN%D)
     !
     CALL Delete(IUT)
     CALL Delete(JUT)
     CALL Delete(UNT)
     !
   END SUBROUTINE UTInvXMatr
!
!--------------------------------------------------------------------
!
   SUBROUTINE UTInvXMatrSymb(IUT,JUT,IB,JB,N,NC,IX,JX)
     INTEGER :: IUT(:),JUT(:),IB(:),JB(:),N,NC,IX(:),JX(:)
     INTEGER :: IBA,IBB,JP,I,JJ,J,IUA,IUB,IXA,IXB,KP,K
     TYPE(INT_VECT) :: IP
     !
     ! X=UT^-1 * B
     ! B,X: NxNC matrix
     ! UT : NxN lower triangle, with unit diagonals (diags not stored)!
     !
     CALL New(IP,NC)
     IP%I(:)=0
     JP=1
     IX(1)=1
     DO 70 I=1,N
       IBA=IB(I)
       IBB=IB(I+1)-1
       IF(IBB<IBA) GO TO 30
       DO 20 JJ=IBA,IBB
         J=JB(JJ)
         IP%I(J)=I
         JX(JP)=J
       20 JP=JP+1
       30 IUA=IUT(I)
       IUB=IUT(I+1)-1
       IF(IUB<IUA) GO TO 60
       DO 50 JJ=IUA,IUB
         J=JUT(JJ)
         IXA=IX(J)
         IXB=IX(J+1)-1
         IF(IXB<IXA) GO TO 50
         DO 40 KP=IXA,IXB
           K=JX(KP)
           IF(IP%I(K)==I) GO TO 40
           IP%I(K)=I
           JX(JP)=K
           JP=JP+1
         40 CONTINUE
       50 CONTINUE
       60 IX(I+1)=JP
     70 CONTINUE
     CALL Delete(IP)
   END SUBROUTINE UTInvXMatrSymb
!
!--------------------------------------------------------------------
!
   SUBROUTINE UTInvXMatrNum(IUT,JUT,UNT,IB,JB,BN,IX,JX,N,NC,XN)
     INTEGER        :: IUT(:),JUT(:),IB(:),JB(:),IX(:),JX(:),N,NC
     REAL(DOUBLE)   :: UNT(:),BN(:),XN(:)
     REAL(DOUBLE)   :: A
     TYPE(DBL_VECT) :: X
     INTEGER        :: I,IH,IXA,IXB,IP,IBA,IBB,IUA,IUB,JP,J,IXC,IXD,KP,K
     !
     ! X=UT^-1 * B
     ! B,X: NxNC matrix
     ! UT : NxN upper triangle, with unit diagonals (diags not stored)!
     !
     CALL New(X,NC)
     !
     DO 80 I=1,N
       IH=I+1
       IXA=IX(I)
       IXB=IX(IH)-1
       IF(IXB<IXA) GO TO 80
       DO 10 IP=IXA,IXB
       10 X%D(JX(IP))=Zero
       IBA=IB(I)
       IBB=IB(IH)-1
       IF(IBB<IBA) GO TO 30
       DO 20 IP=IBA,IBB
       20 X%D(JB(IP))=BN(IP)
       30 IUA=IUT(I)
       IUB=IUT(IH)-1
       IF(IUB<IUA) GO TO 60
       DO 50 JP=IUA,IUB
         J=JUT(JP)
         IXC=IX(J)
         IXD=IX(J+1)-1
         IF(IXD<IXC) GO TO 50
         A=UNT(JP)
         DO 40 KP=IXC,IXD
           K=JX(KP)
         40 X%D(K)=X%D(K)-A*XN(KP)
       50 CONTINUE
       60 CONTINUE
       DO 70 IP=IXA,IXB
       70 XN(IP)=X%D(JX(IP))
     80 CONTINUE
     CALL Delete(X)
   END SUBROUTINE UTInvXMatrNum
!
!--------------------------------------------------------------------
!
   SUBROUTINE UInvXMatrSymb(IU,JU,IB,JB,N,NC,IX,JX)
     INTEGER :: IU(:),JU(:),IB(:),JB(:),N,NC,IX(:),JX(:)
     INTEGER :: IBA,IBB,JP,I,JJ,J,IUA,IUB,IXA,IXB,KP,K,NZMax
     INTEGER :: IStart
     TYPE(INT_VECT) :: IP
     !
     ! X=U^-1 * B
     ! B,X: NxNC matrix
     ! U : NxN lower triangle, with unit diagonals (diags not stored)!
     !
     CALL New(IP,NC)
     NZMax=SIZE(JX)
     IP%I(:)=0
     JP=NZMax
     IX(N+1)=NZMax+1
     DO I=N,1,-1  ! reversing the order
       IBA=IB(I)   ! reverse order
       IBB=IB(I+1)-1
       IF(IBB<IBA) GO TO 30
       DO JJ=IBA,IBB
         J=JB(JJ)
         IP%I(J)=I  ! initialize elements of IP
         JX(JP)=J   ! all elements of B contribute to X
         JP=JP-1    ! count back for sure nonzeros of X
       ENDDO
       30 CONTINUE
       IUA=IU(I)      ! row of U starts here
       IUB=IU(I+1)-1  ! ends here
       IF(IUB<IUA) GO TO 60
       DO JJ=IUA,IUB
         J=JU(JJ)   ! in the Jth colmn of U there is a non-zero
         IXA=IX(J)  ! scan the Jth row of the solution (J<I for sure)
         IXB=IX(J+1)-1
         IF(IXB<IXA) CYCLE
         DO KP=IXA,IXB
           K=JX(KP)
           IF(IP%I(K)==I) CYCLE
           IP%I(K)=I
           JX(JP)=K
           JP=JP-1
         ENDDO
       ENDDO
       60 CONTINUE
       IX(I)=JP+1
     ENDDO
     CALL Delete(IP)
     !
     ! Now, pull back numbering, so that IX(1)=1
     !
     IStart=IX(1)
     J=IStart-1
     DO I=1,N+1
       IX(I)=IX(I)-J
     ENDDO
     JJ=0
     DO I=IStart,NZMax
       JJ=JJ+1
       JX(JJ)=JX(I)
     ENDDO
   END SUBROUTINE UInvXMatrSymb
!
!--------------------------------------------------------------------
!
   SUBROUTINE UInvXMatrNum(IU,JU,UN,IB,JB,BN,IX,JX,N,NC,XN)
     INTEGER        :: IU(:),JU(:),IB(:),JB(:),IX(:),JX(:),N,NC
     REAL(DOUBLE)   :: UN(:),BN(:),XN(:)
     REAL(DOUBLE)   :: A
     TYPE(DBL_VECT) :: X
     INTEGER        :: I,IH,IXA,IXB,IP,IBA,IBB,IUA,IUB,JP,J,IXC,IXD,KP,K
     !
     ! X=U^-1 * B
     ! B,X: NxNC matrix
     ! U : NxN upper triangle, with unit diagonals (diags not stored)!
     !
     CALL New(X,NC)
     !
     DO I=N,1,-1
       IH=I+1
       IXA=IX(I)
       IXB=IX(IH)-1
       IF(IXB<IXA) CYCLE
       DO IP=IXA,IXB
         X%D(JX(IP))=Zero ! zero expanded accumulator for actual row of X
       ENDDO
       IBA=IB(I)
       IBB=IB(IH)-1
       IF(IBB<IBA) GO TO 30
       DO IP=IBA,IBB
         X%D(JB(IP))=BN(IP) ! Put contribution of B into X
       ENDDO
       30 CONTINUE
       !
       IUA=IU(I)
       IUB=IU(IH)-1
       IF(IUB<IUA) GO TO 60
       DO JP=IUA,IUB
         J=JU(JP)
         IXC=IX(J)
         IXD=IX(J+1)-1
         IF(IXD<IXC) CYCLE
         A=UN(JP)
         DO KP=IXC,IXD  ! scan the Jth row of the solution
           K=JX(KP)
           X%D(K)=X%D(K)-A*XN(KP)
         ENDDO
       ENDDO
       60 CONTINUE
       !
       DO IP=IXA,IXB
         XN(IP)=X%D(JX(IP))
       ENDDO
     ENDDO
     CALL Delete(X)
   END SUBROUTINE UInvXMatrNum
!
!--------------------------------------------------------------------
!
   SUBROUTINE CHKGcInv(B,CholData,IUtr,JUtr,AUtr,IUtrT,JUtrT,AUtrT)
     TYPE(BMATR)    :: B
     TYPE(Cholesky) :: CholData
     TYPE(INT_VECT) :: ISpB,JSpB,IGc,JGc,IUtr,JUtr,IUtrT,JUtrT
     TYPE(DBL_VECT) :: ASpB,AGc,AUtr,AUtrT,Vect
     INTEGER        :: NCart,NIntC,NZUtr,I
     !
     NCart=SIZE(CholData%ChRowPt%I)-1
     NIntC=SIZE(B%IB%I,1)
     !
     CALL BtoSpB_1x1(B,ISpB,JSpB,ASpB)
     CALL PermCol(ISpB%I,JSpB%I,CholData%Perm%I)
     CALL ScaleCol(ISpB%I,JSpB%I,ASpB%D,CholData%GcScale%D)
     CALL MatrXUInv(CholData%ChRowPt%I,CholData%ChColPt%I, &
       CholData%ChFact%D,ISpB%I,JSpB%I,ASpB%D, &
       IUtr,JUtr,AUtr,NCart,NIntC)
     !
     CALL New(Vect,NCart)
     DO I=1,NCart ; Vect%D(I)=SQRT(CholData%ChDiag%D(I)) ; ENDDO
     CALL MatrXD(IUtr%I,JUtr%I,AUtr%D,Vect%D,NIntC,NCart)
     CALL Delete(Vect)
     !
     NZUtr=IUtr%I(NIntC+1)-1
     CALL New(IUtrT,NCart+1)
     CALL New(JUtrT,NZUtr)
     CALL New(AUtrT,NZUtr)
     CALL TransPose1x1(IUtr%I,JUtr%I,AUtr%D,NIntC,NCart, &
       IUtrT%I,JUtrT%I,AUtrT%D,'full')
     !
     CALL MatMul_1x1(IUtrT%I,JUtrT%I,AUtrT%D,IUtr%I,JUtr%I,AUtr%D, &
       IGc,JGc,AGc,NCart,NIntC,NCart)
     !
     CALL Show1x1(IGc%I,JGc%I,AGc%D,'identity? ',NCart,NCart)
     !
     CALL Delete(IGc)
     CALL Delete(JGc)
     CALL Delete(AGc)
!    CALL Delete(IUtrT)
!    CALL Delete(JUtrT)
!    CALL Delete(AUtrT)
!    CALL Delete(IUtr)
!    CALL Delete(JUtr)
!    CALL Delete(AUtr)
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
   END SUBROUTINE CHKGcInv
!
!--------------------------------------------------------------------
!
   SUBROUTINE DxMatr(D,IA,JA,AN,N,NC)
     INTEGER,DIMENSION(:)      :: IA,JA
     REAL(DOUBLE),DIMENSION(:) :: AN,D
     REAL(DOUBLE)              :: Sum
     INTEGER                   :: I,J,K,N,NC
     !
     ! A=D*A
     ! A : NxNC matrix
     ! D : NxN
     !
     DO I=1,N
       Sum=D(I)
       DO J=IA(I),IA(I+1)-1
         AN(J)=Sum*AN(J)
       ENDDO
     ENDDO
   END SUBROUTINE DxMatr
!
!--------------------------------------------------------------------
!
   SUBROUTINE MatrXD(IA,JA,AN,D,N,NC)
     INTEGER,DIMENSION(:)      :: IA,JA
     REAL(DOUBLE),DIMENSION(:) :: AN,D
     REAL(DOUBLE)              :: Sum
     INTEGER                   :: I,J,JJ,K,N,NC
     !
     ! A=A*D
     ! A : NCxN matrix
     ! D : NCxNC
     !
     DO I=1,N
       DO J=IA(I),IA(I+1)-1
         JJ=JA(J)
         AN(J)=D(JJ)*AN(J)
       ENDDO
     ENDDO
   END SUBROUTINE MatrXD
!
!----------------------------------------------------------------
!
   SUBROUTINE CholFactGi(ISpB,JSpB,ASpB,NCart,NIntC, &
                       CholData,Print2,Shift_O)
     TYPE(Cholesky)        :: CholData
     TYPE(INT_VECT)        :: ISpB,JSpB
     TYPE(DBL_VECT)        :: ASpB
     TYPE(INT_VECT)        :: ISpBt,JSpBt
     TYPE(DBL_VECT)        :: ASpBt
     INTEGER               :: NZ,NCart,NIntC
     REAL(DOUBLE),OPTIONAL :: Shift_O
     LOGICAL               :: Print2
     !
     NZ=SIZE(JSpB%I)
     CALL New(ISpBt,NCart+1)
     CALL New(JSpBt,NZ)
     CALL New(ASpBt,NZ)
     CALL TransPose1x1(ISpB%I,JSpB%I,ASpB%D,NIntC,NCart, &
                       ISpBt%I,JSpBt%I,ASpBt%D,'full')
     CALL CholFact(ISpBt,JSpBt,ASpBt,NIntC,NCart, &
                   CholData,Print2,Shift_O=Shift_O)
     CALL Delete(ISpBt)
     CALL Delete(JSpBt)
     CALL Delete(ASpBt)
   END SUBROUTINE CholFactGi
!
!----------------------------------------------------------------
!
   SUBROUTINE CholFact(ISpB,JSpB,ASpB,NCart,NIntC, &
                       CholData,Print,Shift_O,ILow_O)
     TYPE(Cholesky) :: CholData
     TYPE(INT_VECT) :: ISpB,JSpB
     TYPE(DBL_VECT) :: ASpB
     TYPE(INT_VECT) :: IGc,JGc
     TYPE(DBL_VECT) :: AGc
     INTEGER        :: NCart,NIntC
     REAL(DOUBLE)   :: SparsitySpB,SparsityGc,SparsityCh
     INTEGER        :: NZSpB,NZGc,NZCh
     LOGICAL        :: Print
     REAL(DOUBLE),OPTIONAL :: Shift_O
     INTEGER,OPTIONAL :: ILow_O
     !
     ! Compute Gc=Bt*B
     !
     NZSpB=ISpB%I(NIntC+1)-1
     CALL GetGc(NCart,ISpB,JSpB,ASpB,IGc,JGc,AGc)
     !
     IF(PRESENT(Shift_O)) THEN
       CALL SpectrShift_1x1(IGc,JGc,AGc,Shift_O)
     ENDIF
     IF(PRESENT(ILow_O)) THEN
      !CALL Plot_1x1(IGc%I,JGc%I,'Gc_1',NCart)
       CALL CompleteDiag_1x1(IGc,JGc,AGc,ILow_O)
      !CALL Plot_1x1(IGc%I,JGc%I,'Gc_2',NCart)
     ENDIF
     NZGc=IGc%I(NCart+1)-1
     !
     CALL TriangFact(IGc,JGc,AGc,CholData)
       NZCh=CholData%ChRowPt%I(NCart+1)-1

     IF(Print) THEN
       SparsitySpB=DBLE(NZSpB)/DBLE(NIntC*NCart)*100.D0
       SparsityGc=DBLE(NZGc)/DBLE(NCart*NCart)*100.D0
       SparsityCh=DBLE(NZCh)/DBLE(NCart*NCart-NCart)*200.D0

       CALL MondoLog(DEBUG_NONE, "CholFact", 'Sparsity of  B = '//TRIM(FltToMedmChar(SparsitySpB)) &
                    //" %, # of SpB%NNon0 = "//TRIM(IntToChar(NZSpB)))
       CALL MondoLog(DEBUG_NONE, "CholFact", 'Sparsity of Gc = '//TRIM(FltToMedmChar(SparsityGc))  &
                    //" %, # of  Gc%NNon0 = "//TRIM(IntToChar(NZGc)))
       CALL MondoLog(DEBUG_NONE, "CholFact", 'Sparsity of Ch = '//TRIM(FltToMedmChar(SparsityCh))  &
                    //" %, # of  Ch%NNon0 = "//TRIM(IntToChar(NZCh)))
     ENDIF
   END SUBROUTINE CholFact
!
!--------------------------------------------------------------------
!
   SUBROUTINE TestChFact(ChDiag,ChRowPt,ChColPt,ChFact)
     TYPE(INT_VECT) :: ChRowPt,ChColPt
     TYPE(DBL_VECT) :: ChDiag,ChFact
     INTEGER        :: I,J,K,L,NDim
     TYPE(DBL_RNK2) :: Aux,Aux2,Aux3
     REAL(DOUBLE)   :: Sum
     !
     NDim=SIZE(ChDiag%D)
     CALL New(Aux,(/NDim,NDim/))
     CALL New(Aux2,(/NDim,NDim/))
     CALL New(Aux3,(/NDim,NDim/))
     Aux%D=Zero
     Aux3%D=Zero
     DO I=1,NDim
         Aux3%D(I,I)=One/ChDiag%D(I)
         Aux%D(I,I)=One
       DO J=ChRowPt%I(I),ChRowPt%I(I+1)-1
         K=ChColPt%I(J)
         Aux%D(I,K)=ChFact%D(J)
       ENDDO
     ENDDO
     !
    !CALL PPrint(Aux,'Cholesky factor',unit_o=6)
     WRITE(*,*) 'Cholesky factors'
     DO I=1,NDim
       WRITE(*,123) (Aux%D(I,J),J=1,NDim)
     ENDDO
     CALL DGEMM_NNc(NDim,NDim,NDim,One,Zero,Aux3%D,Aux%D, &
                         Aux2%D)
WRITE(*,*) 'aux2'
DO I=1,NDim
  WRITE(*,123) (Aux2%D(I,J),J=1,NDim)
ENDDO
     CALL DGEMM_TNc(NDim,NDim,NDim,One,Zero,Aux%D,Aux2%D, &
                         Aux3%D)
    !CALL PPrint(Aux3,'check Cholesky factors',unit_o=6)
     WRITE(*,*) 'check Cholesky factors'
     DO I=1,NDim
       WRITE(*,123) (Aux3%D(I,J),J=1,NDim)
     ENDDO
     123 FORMAT(20F6.3)
     !
     CALL Delete(Aux)
     CALL Delete(Aux2)
     CALL Delete(Aux3)
   END SUBROUTINE TestChFact
!
!--------------------------------------------------------------------
!
   SUBROUTINE PermVect(VectIn,VectOut,Perm)
     INTEGER,DIMENSION(:) :: Perm
     REAL(DOUBLE),DIMENSION(:) :: VectIn,VectOut
     INTEGER                   :: I,J,NCart
     NCart=SIZE(Perm)
     DO I=1,NCart
       VectOut(Perm(I))=VectIn(I)
     ENDDO
   END SUBROUTINE PermVect
!
!--------------------------------------------------------------------
!
   SUBROUTINE PermCol(IA,JA,Perm)
     INTEGER,DIMENSION(:) :: IA,JA,Perm
     INTEGER              :: I,J,K,L,N
     !
     N=SIZE(IA)-1
     DO I=1,N
       DO J=IA(I),IA(I+1)-1
         JA(J)=Perm(JA(J))
       ENDDO
     ENDDO
   END SUBROUTINE PermCol
!
!--------------------------------------------------------------------
!
   SUBROUTINE PermRow(IA,JA,AN,Perm,N,M)
     INTEGER,DIMENSION(:)      :: IA,JA,Perm
     REAL(DOUBLE),DIMENSION(:) :: AN
     INTEGER                   :: I,J,K,L,N,NZ,M
     TYPE(INT_VECT)            :: IAT,JAT
     TYPE(DBL_VECT)            :: ANT
     !
     NZ=SIZE(JA)
     CALL New(IAT,M+1)
     CALL New(JAT,NZ)
     CALL New(ANT,NZ)
     !
     CALL TransPose1x1(IA,JA,AN,N,M,IAT%I,JAT%I,ANT%D,'full')
     CALL PermCol(IAT%I,JAT%I,Perm)
     CALL TransPose1x1(IAT%I,JAT%I,ANT%D,M,N,IA,JA,AN,'full')
     CALL Delete(IAT)
     CALL Delete(JAT)
     CALL Delete(ANT)
   END SUBROUTINE PermRow
!
!--------------------------------------------------------------------
!
   SUBROUTINE ScaleCol(IA,JA,AN,Scale)
     INTEGER,DIMENSION(:)      :: IA,JA
     REAL(DOUBLE),DIMENSION(:) :: AN,Scale
     INTEGER                   :: I,J,K,JJ,N
     !
     N=SIZE(IA)-1
     DO I=1,N
       DO J=IA(I),IA(I+1)-1
         JJ=JA(J)
         AN(J)=Scale(JJ)*AN(J)
       ENDDO
     ENDDO
   END SUBROUTINE ScaleCol
!
!--------------------------------------------------------------------
!
   SUBROUTINE ScaleRow(IA,JA,AN,Scale)
     INTEGER,DIMENSION(:)      :: IA,JA
     REAL(DOUBLE),DIMENSION(:) :: AN,Scale
     INTEGER                   :: I,J,K,JJ,N
     !
     N=SIZE(IA)-1
     DO I=1,N
       DO J=IA(I),IA(I+1)-1
         AN(J)=Scale(I)*AN(J)
       ENDDO
     ENDDO
   END SUBROUTINE ScaleRow
!
!--------------------------------------------------------------------
!
   SUBROUTINE ScaleVect(Vect,Scale)
     REAL(DOUBLE),DIMENSION(:) :: Vect,Scale
     INTEGER                   :: I,N
     N=SIZE(Scale)
     DO I=1,N
       Vect(I)=Vect(I)*Scale(I)
     ENDDO
   END SUBROUTINE ScaleVect
!
!--------------------------------------------------------------------
!
   subroutine genrcm ( n, xadj, iadj, perm, xls )
     integer n
     integer iadj(*)
     integer i
     integer iccsze
     integer mask(n)
     integer nlvl
     integer num
     integer perm(n)
     integer root
     integer xadj(n+1)
     integer xls(n+1)
     !
     !  MASK marks variables that have been numbered.
     !
     mask(1:n)=1
     num = 1
     do i = 1, n
       !
       !  For each masked connected component...
       !
       if ( mask(i) /= 0 ) then
         root = i
         !
         !  Find a pseudo-peripheral node root.
         !  Note that the level structure found by
         !  fnroot is stored starting at perm(num).
         !
            call fnroot (root,xadj,iadj,mask,nlvl,xls,perm(num),n)
         !
         !  RCM orders the component using ROOT as the starting node.
         !
         call rcm ( root, xadj, iadj, mask, perm(num), iccsze, xls, n )
         num = num + iccsze
         if ( num > n ) then
           return
         end if
       end if
     end do
   end SUBROUTINE GenRCM
!
!--------------------------------------------------------------------
!
   subroutine rcm ( root, xadj, iadj, mask, perm, iccsze, deg, n )
     integer n
     integer iadj(*)
     integer deg(n)
     integer fnbr
     integer i
     integer iccsze
     integer j
     integer jstop
     integer jstrt
     integer k
     integer l
     integer lbegin
     integer lnbr
     integer lperm
     integer lvlend
     integer mask(n)
     integer nbr
     integer node
     integer perm(n)
     integer root
     integer xadj(n+1)
     !
     !  Find the degrees of the nodes in the
     !  component specified by MASK and ROOT.
     !
     call degree ( root, xadj, iadj, mask, deg, iccsze, perm, n )
     mask(root) = 0
     if ( iccsze <= 1 ) then
       return
     end if
     lvlend = 0
     lnbr = 1
     !
     !  LBEGIN and LVLEND point to the beginning and
     !  the end of the current level respectively.
     !
     10 continue
     lbegin = lvlend + 1
     lvlend = lnbr
     do i = lbegin, lvlend
       !
       !  For each node in current level...
       !
       node = perm(i)
       jstrt = xadj(node)
       jstop = xadj(node+1)-1
       !
       !  Find the unnumbered neighbors of NODE.
       !
       !  FNBR and LNBR point to the first and last unnumbered neighbors
       !  of the current node in PERM.
       !
       fnbr = lnbr+1
       do j = jstrt, jstop
         nbr = iadj(j)
         if ( mask(nbr) /= 0 ) then
           lnbr = lnbr+1
           mask(nbr) = 0
           perm(lnbr) = nbr
         end if
       end do
       if ( fnbr >= lnbr ) then
         goto 60
       end if
       !
       !  Sort the neighbors of node in increasing order by degree.
       !  Linear insertion is used.
       !
       k = fnbr
       30   continue
       l = k
       k = k+1
       nbr = perm(k)
       40   continue
       if ( l > fnbr ) then
         lperm = perm(l)
         if ( deg(lperm) > deg(nbr) ) then
           perm(l+1) = lperm
           l = l-1
           goto 40
         end if
       end if
       perm(l+1) = nbr
       if ( k < lnbr ) then
         goto 30
       end if
       60   continue
     end do
     if ( lnbr > lvlend ) then
       goto 10
     end if
     !
     !  We now have the Cuthill-McKee ordering.  Reverse it.
     !
       call ivec_reverse ( iccsze, perm )
   end SUBROUTINE RCM
!
!--------------------------------------------------------------------
!
   subroutine ivec_reverse ( n, a )
     integer n
     integer a(n)
     integer i
     do i = 1, n/2
       call i_swap ( a(i), a(n+1-i) )
     end do
   end subroutine ivec_reverse
!
!--------------------------------------------------------------------
!
   subroutine i_swap ( i, j )
   !
     integer i
     integer j
     integer k
   !
     k = i
     i = j
     j = k
   end subroutine i_swap
!
!--------------------------------------------------------------------
!
   subroutine degree ( root, xadj, iadj, mask, deg, iccsze, ls, n )
     integer n
     integer iadj(*)
     integer deg(n)
     integer i
     integer iccsze
     integer ideg
     integer j
     integer jstop
     integer jstrt
     integer lbegin
     integer ls(n)
     integer lvlend
     integer lvsize
     integer mask(n)
     integer nbr
     integer node
     integer root
     integer xadj(n+1)
     !
     !  The array XADJ is used as a temporary marker to
     !  indicate which nodes have been considered so far.
     !
     ls(1) = root
     xadj(root) = -xadj(root)
     lvlend = 0
     iccsze = 1
     !
     !  LBEGIN is the pointer to the beginning of the current level, and
     !  LVLEND points to the end of this level.
     !
     10 continue
     lbegin = lvlend+1
     lvlend = iccsze
     !
     !  Find the degrees of nodes in the current level,
     !  and at the same time, generate the next level.
     !
     do i = lbegin, lvlend
       node = ls(i)
       jstrt = -xadj(node)
       jstop = abs ( xadj(node+1) ) - 1
       ideg = 0
       do j = jstrt, jstop
         nbr = iadj(j)
         if ( mask(nbr) /= 0 ) then
           ideg = ideg+1
           if ( xadj(nbr) >= 0 ) then
             xadj(nbr) = -xadj(nbr)
             iccsze = iccsze+1
             ls(iccsze) = nbr
           end if
         end if
       end do
       deg(node) = ideg
     end do
     !
     !  Compute the current level width.
     !
     lvsize = iccsze - lvlend
     !
     !  If the current level width is nonzero, generate another level.
     !
     if ( lvsize > 0 ) then
       goto 10
     end if
     !
     !  Reset XADJ to its correct sign and return.
     !
     do i = 1, iccsze
       node = ls(i)
       xadj(node) = -xadj(node)
     end do
   end subroutine degree
!
!--------------------------------------------------------------------
!
   subroutine fnroot ( root, xadj, iadj, mask, nlvl, xls, ls, n )
     integer n
     integer iadj(*)
     integer iccsze
     integer j
     integer jstrt
     integer k
     integer kstop
     integer kstrt
     integer ls(n)
     integer mask(n)
     integer mindeg
     integer nabor
     integer ndeg
     integer nlvl
     integer node
     integer nunlvl
     integer root
     integer xadj(n+1)
     integer xls(n+1)
     !
     !  Determine the level structure rooted at ROOT.
     !
     call rootls ( root, xadj, iadj, mask, nlvl, xls, ls, n )
     iccsze = xls(nlvl+1)-1
     if ( nlvl == 1 .or. nlvl == iccsze ) then
       return
     end if
     !
     !  Pick a node with minimum degree from the last level.
     !
     10 continue
     jstrt = xls(nlvl)
     mindeg = iccsze
     root = ls(jstrt)
     if ( iccsze > jstrt ) then
       do j = jstrt, iccsze
         node = ls(j)
         ndeg = 0
         kstrt = xadj(node)
         kstop = xadj(node+1)-1
         do k = kstrt, kstop
           nabor = iadj(k)
           if ( mask(nabor) > 0 ) then
             ndeg = ndeg+1
           end if
         end do
         if ( ndeg < mindeg ) then
           root = node
           mindeg = ndeg
         end if
       end do
     end if
     !
     !  Generate its rooted level structure.
     !
     call rootls ( root, xadj, iadj, mask, nunlvl, xls, ls, n )
     if ( nunlvl <= nlvl ) then
       return
     end if
     nlvl = nunlvl
     if ( nlvl < iccsze ) then
       goto 10
     end if
   end subroutine fnroot
!
!--------------------------------------------------------------------
!
   subroutine rootls ( root, xadj, iadj, mask, nlvl, xls, ls, n )
     integer n
     integer iadj(*)
     integer i
     integer iccsze
     integer j
     integer jstop
     integer jstrt
     integer lbegin
     integer ls(n)
     integer lvlend
     integer lvsize
     integer mask(n)
     integer nbr
     integer nlvl
     integer node
     integer root
     integer xadj(n+1)
     integer xls(n+1)
     !
     mask(root) = 0
     ls(1) = root
     nlvl = 0
     lvlend = 0
     iccsze = 1
     !
     !  LBEGIN is the pointer to the beginning of the current level, and
     !  LVLEND points to the end of this level.
     !
     10 continue
     lbegin = lvlend + 1
     lvlend = iccsze
     nlvl = nlvl + 1
     xls(nlvl) = lbegin
     !
     !  Generate the next level by finding all the
     !  masked neighbors of nodes
     !  in the current level.
     !
     do i = lbegin, lvlend
       node = ls(i)
       jstrt = xadj(node)
       jstop = xadj(node+1)-1
       do j = jstrt, jstop
         nbr = iadj(j)
         if ( mask(nbr) /= 0 ) then
           iccsze = iccsze + 1
           ls(iccsze) = nbr
           mask(nbr) = 0
         end if
       end do
     end do
     !
     !  Compute the current level width.
     !  If it is nonzero, generate the next level.
     !
     lvsize = iccsze-lvlend
     if ( lvsize > 0 ) then
       goto 10
     end if
     !
     !  Reset MASK to one for the nodes in the level structure.
     !
     xls(nlvl+1) = lvlend + 1
     do i = 1, iccsze
       node = ls(i)
       mask(node) = 1
     end do
   end subroutine rootls
!
!--------------------------------------------------------------------
!
   SUBROUTINE AddMat_1x1(IA,JA,AN,IB,JB,BN,IC,JC,CN,N,M,SymbOnly_O)
     ! addition of two NxM matrices
     INTEGER,DIMENSION(:)      :: IA,JA,IB,JB
     REAL(DOUBLE),DIMENSION(:) :: AN,BN
     INTEGER                   :: I,J,K,L,N,M,NZC
     TYPE(INT_VECT)            :: IC,JC,JC1
     TYPE(DBL_VECT)            :: CN
     LOGICAL,OPTIONAL          :: SymbOnly_O
     !
     CALL New(IC,N+1)
     NZC=IA(N+1)-1 + IB(N+1)-1
     CALL New(JC1,NZC)
     !
     CALL AddSymb(IA,JA,IB,JB,N,M,IC%I,JC1%I)
     !
     NZC=IC%I(N+1)-1
     CALL New(JC,NZC)
     JC%I(1:NZC)=JC1%I(1:NZC)
     CALL Delete(JC1)
     !
     IF(PRESENT(SymbOnly_O)) THEN
       IF(SymbOnly_O) RETURN
     ENDIF
     !
     CALL New(CN,NZC)
     CALL AddNum(iA,jA,AN,iB,jB,BN,N,M,iC%I,jC%I,CN%D)
   END SUBROUTINE AddMat_1x1
!
!--------------------------------------------------------------------
!
   SUBROUTINE AddSymb(IA,JA,IB,JB,N,M,IC,JC)
     ! symbolic addition of two NxM matrices
     INTEGER,DIMENSION(:) :: IA,JA,IB,JB,IC,JC
     TYPE(INT_VECT)       :: IX
     INTEGER              :: N,M,IP,I,J,IAA,IAB,JP,IBA,IBB
     !
     CALL New(IX,M)
     !
     IP=1
     IX%I=0
     DO 50 I=1,N
       IC(I)=IP
       IAA=IA(I)
       IAB=IA(I+1)-1
       IF(IAB<IAA) GO TO 30
       DO 20 JP=IAA,IAB
         J=JA(JP)
         JC(IP)=J
         IP=IP+1
       20 IX%I(J)=I
       30 IBA=IB(I)
       IBB=IB(I+1)-1
       IF(IBB<IBA) GO TO 50
       DO 40 JP=IBA,IBB
         J=JB(JP)
         IF(IX%I(J)==I) GO TO 40
         JC(IP)=J
         IP=IP+1
       40 CONTINUE
     50 CONTINUE
       IC(N+1)=IP
     !
     CALL Delete(IX)
   END SUBROUTINE AddSymb
!
!--------------------------------------------------------------------
!
   SUBROUTINE AddNum(iA,jA,sA,iB,jB,sB,N,M,iC,jC,sC)
     !numerical addition of two NxM, 1x1 blocked matrices
     INTEGER,DIMENSION(:)      :: iA,jA,iB,jB,iC,jC
     REAL(DOUBLE),DIMENSION(:) :: sA,sB,sC
     TYPE(DBL_VECT)            :: W
     INTEGER                   :: N,M,I,J,K
     !
     CALL New(W,M)
     DO 100 I=1,N
       DO 200 K=iC(I),iC(I+1)-1
         W%D(jC(K))=Zero
       200     CONTINUE
       DO 300 K=iA(I),iA(I+1)-1
         W%D(jA(K))=sA(K)
       300     CONTINUE
       DO 400 K=iB(I),iB(I+1)-1
         J=jB(K)
         W%D(J)=W%D(J)+sB(K)
       400     CONTINUE
       DO 500 K=iC(I),iC(I+1)-1
         sC(K)=W%D(jC(K))
       500     CONTINUE
     100  CONTINUE
     CALL Delete(W)
   END SUBROUTINE AddNum
!
!--------------------------------------------------------------------
!
   SUBROUTINE ThreshMatr(IGc,JGc,AGc,Thresh)
     TYPE(INT_VECT)            :: IGc,JGc,IGc2,JGc2
     TYPE(DBL_VECT)            :: AGc,AGc2
     INTEGER                   :: I,J,NRow,NZ
     REAL(DOUBLE)              :: Sum,Thresh
     !
     NRow=SIZE(IGc%I)-1
     NZ=SIZE(JGc%I)
     CALL New(IGc2,NRow+1)
     CALL New(JGc2,NZ)
     CALL New(AGc2,NZ)
     NZ=0
     IGc2%I(1)=1
     DO I=1,NRow
       DO J=IGc%I(I),IGc%I(I+1)-1
         Sum=AGc%D(J)
         IF(ABS(Sum)>Thresh) THEN
           NZ=NZ+1
           JGc2%I(NZ)=JGc%I(J)
           AGc2%D(NZ)=Sum
         ENDIF
       ENDDO
       IGc2%I(I+1)=NZ+1
     ENDDO
     !
     IGc%I=IGc2%I
     JGc%I(1:NZ)=JGc2%I(1:NZ)
     AGc%D(1:NZ)=AGc2%D(1:NZ)
     !
     CALL Delete(IGc2)
     CALL Delete(JGc2)
     CALL Delete(AGc2)
   END SUBROUTINE ThreshMatr
!
!--------------------------------------------------------------------
!
   SUBROUTINE SpectrShift_1x1(IGc,JGc,AGc,Shift)
     TYPE(INT_VECT)      :: IGc,JGc,IUnit,JUnit,IC,JC
     TYPE(DBL_VECT)      :: AGc,AUnit,CN
     REAL(DOUBLE)        :: Shift
     INTEGER             :: NDim,I,J,NZ
     !
     NDim=SIZE(IGc%I)-1
     CALL GetUnit_1x1(IUnit,JUnit,AUnit,NDim)
     AUnit%D=Shift*AUnit%D
     CALL AddMat_1x1(IGc%I,JGc%I,AGc%D,IUnit%I,JUnit%I,AUnit%D,IC,JC,CN,NDim,NDim)
     CALL Delete(IGc)
     CALL Delete(JGc)
     CALL Delete(AGc)
     NZ=IC%I(NDim+1)-1
     CALL New(IGc,NDim+1)
     CALL New(JGc,NZ)
     CALL New(AGc,NZ)
     IGc%I(1:NDim+1)=IC%I(1:NDim+1)
     JGc%I(1:NZ)=JC%I(1:NZ)
     AGc%D(1:NZ)=CN%D(1:NZ)
     !
     CALL Delete(IC)
     CALL Delete(JC)
     CALL Delete(CN)
     CALL Delete(IUnit)
     CALL Delete(JUnit)
     CALL Delete(AUnit)
   END SUBROUTINE SpectrShift_1x1
!
!--------------------------------------------------------------------
!
   SUBROUTINE GetUnit_1x1(IUnit,JUnit,AUnit,NDim)
     TYPE(INT_VECT)      :: IUnit,JUnit
     TYPE(DBL_VECT)      :: AUnit
     INTEGER             :: NDim,I,J
     !
     CALL New(IUnit,NDim+1)
     CALL New(JUnit,NDim)
     CALL New(AUnit,NDim)
     !
     IUnit%I(1)=1
     DO I=1,NDim
       IUnit%I(I+1)=I+1
       JUnit%I(I)=I
       AUnit%D(I)=One
     ENDDO
   END SUBROUTINE GetUnit_1x1
!
!--------------------------------------------------------------------
!
   SUBROUTINE CompleteDiag_1x1(IGc,JGc,AGc,ILow)
     TYPE(INT_VECT)      :: IGc,JGc,IUnit,JUnit,IC,JC
     TYPE(DBL_VECT)      :: AGc,AUnit,CN
     REAL(DOUBLE)        :: Shift
     INTEGER             :: NDim,I,J,NZ,K,ILow
     !
     NDim=SIZE(IGc%I)-1
     CALL GetUnit_1x1(IUnit,JUnit,AUnit,NDim)
     AUnit%D=Zero
     DO I=ILow+1,NDim
       AUnit%D(I)=One
     ENDDO
     !
     CALL AddMat_1x1(IGc%I,JGc%I,AGc%D,IUnit%I,JUnit%I,AUnit%D, &
                     IC,JC,CN,NDim,NDim)
     CALL Delete(IGc)
     CALL Delete(JGc)
     CALL Delete(AGc)
     NZ=IC%I(NDim+1)-1
     CALL New(IGc,NDim+1)
     CALL New(JGc,NZ)
     CALL New(AGc,NZ)
     IGc%I(1:NDim+1)=IC%I(1:NDim+1)
     JGc%I(1:NZ)=JC%I(1:NZ)
     AGc%D(1:NZ)=CN%D(1:NZ)
     !
     CALL Delete(IC)
     CALL Delete(JC)
     CALL Delete(CN)
     CALL Delete(IUnit)
     CALL Delete(JUnit)
     CALL Delete(AUnit)
   END SUBROUTINE CompleteDiag_1x1
!
!-------------------------------------------------------------------
!
   SUBROUTINE TriangFact(IGc,JGc,AGc,CholData)
     TYPE(Cholesky) :: CholData
     TYPE(INT_VECT) :: BorderLine
     TYPE(INT_VECT) :: ChColPt1
     TYPE(DBL_VECT) :: ChDiag,ChFact
     TYPE(INT_VECT) :: IGc,JGc
     TYPE(INT_VECT) :: IGcU,JGcU
     TYPE(DBL_VECT) :: DGc,AGc
     TYPE(DBL_VECT) :: DGcU,AGcU
     INTEGER        :: N,I,J,K,I1,I2,J1,J2,K1,K2,R,S,NZ
     INTEGER        :: NCart,GcSNon0,ChFillEst,NIntC
     INTEGER        :: NZSpB,NZGc,NZCh
     !
     NZGc=SIZE(JGc%I)
     NCart=SIZE(IGc%I)-1
     !
     ! Find RCM ordering
     !
     CALL New(CholData%Perm,NCart)
     CALL New(CholData%IPerm,NCart)
     CALL RCMOrder(CholData%IPerm%I,CholData%Perm%I,NCart,IGc%I,JGc%I)
     !
     ! Permute Gc
     !
     CALL Perm1x1(CholData%Perm%I,IGc%I,JGc%I,AGc%D)
     !
     ! Determine profile of Gc
     !
     CALL New(BorderLine,NCart)
     CALL GetBorders(BorderLine%I,IGc%I,JGc%I,NCart)
     !
     ! Reduce Gc to upper triangle
     !
     NZ=(IGc%I(NCart+1)-1-NCart)/2
     CALL New(IGcU,NCart+1)
     CALL New(JGcU,NZ)
     CALL New(AGcU,NZ)
     CALL New(DGcU ,NCart)
       CALL UpperTr(IGc%I,JGc%I,AGc%D,IGcU%I,JGcU%I,AGcU%D,DGcU%D)
     CALL Delete(IGc)
     CALL Delete(JGc)
     CALL Delete(AGc)
     !
     ! Scale matrix
     !
     CALL New(CholData%GcScale,NCart)
     DO I=1,NCart
       CholData%GcScale%D(I)=One/SQRT(DGcU%D(I))
     ENDDO
     CALL Scale1x1(CholData%GcScale%D,IGcU%I,JGcU%I,AGcU%D,DGcU%D,NCart)
     !
     ! Allocate storage for symbolic Cholesky factor
     !
     ChFillEst=NCart*500
     CALL New(CholData%ChRowPt,NCart+1)
     CALL New(ChColPt1,ChFillEst)
     !
     ! Determine symbolic structure of the Cholesky factor
     !
     CALL SymbCholComplete(BorderLine%I,NCart,&
               CholData%ChRowPt%I,ChColPt1%I)
     NZCh=CholData%ChRowPt%I(NCart+1)-1
     !
     CALL New(CholData%ChColPt,NZCh)
     CholData%ChColPt%I(1:NZCh)=ChColPt1%I(1:NZCh)
     CALL Delete(ChColPt1)
     !
     ! Numeric Cholesky factorization
     !
     CALL New(CholData%ChDiag,NCart)
     CALL New(CholData%ChFact,NZCh)
     CALL NumCholFact(IGcU%I,JGcU%I,AGcU%D, &
       DGcU%D,CholData%ChRowPt%I,CholData%ChColPt%I,NCart, &
       CholData%ChFact%D,CholData%ChDiag%D,.FALSE.)
     !
     CALL Delete(JGcU)
     CALL Delete(IGcU)
     CALL Delete(AGcU)
     CALL Delete(DGcU)
     CALL Delete(BorderLine)
   END SUBROUTINE TriangFact
!
!--------------------------------------------------------------------
!
END MODULE
