MODULE GetAtomPairMod
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE ShellPairStruct  
CONTAINS

  SUBROUTINE GetAtomPair2(AtmInfo,AtmPair,BS,PBC)
    TYPE(AtomInfo)                :: AtmInfo
    TYPE(AtomPr)  , DIMENSION(:)  :: AtmPair
    TYPE(BSET)                    :: BS
    REAL(DOUBLE)  , DIMENSION(3)  :: PBC
    INTEGER      :: CF12,CF1,CF2,I1,I2,II,JJ,IJ,iCell
    INTEGER      :: MinL1,MaxL1,Type1,MinL2,MaxL2,Type2
    INTEGER      :: StartL1,StopL1,StartL2,StopL2
    REAL(DOUBLE) :: Z1,Z2,Expt,InvExpt,R12,XiR12,RX,RY,RZ,CNT

    RX=PBC(1)
    RY=PBC(2)
    RZ=PBC(3)

    AtmInfo%Atm12X=AtmInfo%Atm1X-AtmInfo%Atm2X
    AtmInfo%Atm12Y=AtmInfo%Atm1Y-AtmInfo%Atm2Y
    AtmInfo%Atm12Z=AtmInfo%Atm1Z-AtmInfo%Atm2Z

    R12=(AtmInfo%Atm12X-RX)**2+(AtmInfo%Atm12Y-RY)**2+(AtmInfo%Atm12Z-RZ)**2
    CF12=0
    DO CF1=1,BS%NCFnc%I(AtmInfo%K1)
       MinL1=BS%ASymm%I(1,CF1,AtmInfo%K1)
       MaxL1=BS%ASymm%I(2,CF1,AtmInfo%K1)
       Type1=MaxL1*(MaxL1+1)/2+MinL1+1
       StartL1=BS%LStrt%I(CF1,AtmInfo%K1)
       StopL1=BS%LStop%I(CF1,AtmInfo%K1)
       DO CF2=1,BS%NCFnc%I(AtmInfo%K2)
          CF12=CF12+1
          MinL2=BS%ASymm%I(1,CF2,AtmInfo%K2)
          MaxL2=BS%ASymm%I(2,CF2,AtmInfo%K2)
          Type2=MaxL2*(MaxL2+1)/2+MinL2+1
             AtmPair(CF12)%SP%IntType=Type1*100+Type2
             StartL2=BS%LStrt%I(CF2,AtmInfo%K2)
             StopL2=BS%LStop%I(CF2,AtmInfo%K2)
             II=0
             ! Assume primitives are ordered (exponents in decressing order).
             DO I1=BS%NPFnc%I(CF1,AtmInfo%K1),1,-1
                Z1=BS%Expnt%D(I1,CF1,AtmInfo%K1)
                JJ=0
                DO I2=BS%NPFnc%I(CF2,AtmInfo%K2),1,-1
                   Z2=BS%Expnt%D(I2,CF2,AtmInfo%K2)
                   Expt=Z1+Z2
                   InvExpt=1.0d0/Expt
                   XiR12=Z2*Z1*InvExpt*R12
                   JJ=JJ+1
                   IJ=JJ+II
                   AtmPair(CF12)%SP%Cst(1,IJ)=Expt
                   AtmPair(CF12)%SP%Cst(2,IJ)=(Z1*AtmInfo%Atm1X+Z2*(AtmInfo%Atm2X+RX))*InvExpt
                   AtmPair(CF12)%SP%Cst(3,IJ)=(Z1*AtmInfo%Atm1Y+Z2*(AtmInfo%Atm2Y+RY))*InvExpt
                   AtmPair(CF12)%SP%Cst(4,IJ)=(Z1*AtmInfo%Atm1Z+Z2*(AtmInfo%Atm2Z+RZ))*InvExpt
                   AtmPair(CF12)%SP%Cst(5,IJ)=5.914967172796D0*EXP(-XiR12)*InvExpt*  &
                                             BS%CCoef%D(StopL1,I1,CF1,AtmInfo%K1)*   &
                                             BS%CCoef%D(StopL2,I2,CF2,AtmInfo%K2)
                   IF((Type1.NE.2.AND.Type2==2).OR.(Type2.NE.2.AND.Type1==2))THEN
                      Cnt=BS%CCoef%D(StopL1,I1,CF1,AtmInfo%K1)*BS%CCoef%D(StopL2,I2,CF2,AtmInfo%K2)
                      AtmPair(CF12)%SP%Cst(6,IJ)=BS%CCoef%D(StartL1,  I1,CF1,AtmInfo%K1)*BS%CCoef%D(StartL2,I2,CF2,AtmInfo%K2)/Cnt
                      AtmPair(CF12)%SP%Cst(7,IJ)=BIG_DBL
                      AtmPair(CF12)%SP%Cst(8,IJ)=BIG_DBL
                   ELSEIF(Type1==2.AND.Type2==2)THEN
                      Cnt=BS%CCoef%D(StopL1,I1,CF1,AtmInfo%K1)*BS%CCoef%D(StopL2,I2,CF2,AtmInfo%K2)
                      AtmPair(CF12)%SP%Cst(6,IJ)=BS%CCoef%D(StartL1,  I1,CF1,AtmInfo%K1)*BS%CCoef%D(StartL2,  I2,CF2,AtmInfo%K2)/Cnt
                      AtmPair(CF12)%SP%Cst(7,IJ)=BS%CCoef%D(StartL1+1,I1,CF1,AtmInfo%K1)*BS%CCoef%D(StartL2,  I2,CF2,AtmInfo%K2)/Cnt
                      AtmPair(CF12)%SP%Cst(8,IJ)=BS%CCoef%D(StartL1  ,I1,CF1,AtmInfo%K1)*BS%CCoef%D(StartL2+1,I2,CF2,AtmInfo%K2)/Cnt
                   ELSE
                      AtmPair(CF12)%SP%Cst(6,IJ)=BIG_DBL
                      AtmPair(CF12)%SP%Cst(7,IJ)=BIG_DBL
                      AtmPair(CF12)%SP%Cst(8,IJ)=BIG_DBL
                   ENDIF
                ENDDO
                ! We can skipp out the loop because we did not get any significant primitives.
                !                   IF(JJ.EQ.0) EXIT
                II=II+JJ
                !
             ENDDO
!             AtmPair(CF12)%SP%L=II
!             AtmPair(CF12)%SP%AtmInfo%Atm1X=AtmInfo%Atm1X
!             AtmPair(CF12)%SP%AtmInfo%Atm1Y=AtmInfo%Atm1Y
!             AtmPair(CF12)%SP%AtmInfo%Atm1Z=AtmInfo%Atm1Z
!             AtmPair(CF12)%SP%AtmInfo%Atm2X=AtmInfo%Atm2X+RX
!             AtmPair(CF12)%SP%AtmInfo%Atm2Y=AtmInfo%Atm2Y+RY
!             AtmPair(CF12)%SP%AtmInfo%Atm2Z=AtmInfo%Atm2Z+RZ




       AtmPair(CF12)%SP%L=II
       !
       ! We reorder the atomic positions if Type2 > Type1.
       ! Needed for the integral evaluations.
       IF(Type1.GE.Type2) THEN
          AtmPair(CF12)%SP%AtmInfo%Atm1X=AtmInfo%Atm1X
          AtmPair(CF12)%SP%AtmInfo%Atm1Y=AtmInfo%Atm1Y
          AtmPair(CF12)%SP%AtmInfo%Atm1Z=AtmInfo%Atm1Z
          !
          ! AtmInfo must be related to the atoms in the working cell ONLY.
          ! Then we add the PBC's to have the right atomic position.
          AtmPair(CF12)%SP%AtmInfo%Atm2X=AtmInfo%Atm2X+RX
          AtmPair(CF12)%SP%AtmInfo%Atm2Y=AtmInfo%Atm2Y+RY
          AtmPair(CF12)%SP%AtmInfo%Atm2Z=AtmInfo%Atm2Z+RZ
       ELSE
          !
          ! AtmInfo must be related to the atoms in the working cell ONLY.
          ! Then we add the PBC's to have the right atomic position.
          AtmPair(CF12)%SP%AtmInfo%Atm1X=AtmInfo%Atm2X+RX
          AtmPair(CF12)%SP%AtmInfo%Atm1Y=AtmInfo%Atm2Y+RY
          AtmPair(CF12)%SP%AtmInfo%Atm1Z=AtmInfo%Atm2Z+RZ
          AtmPair(CF12)%SP%AtmInfo%Atm2X=AtmInfo%Atm1X
          AtmPair(CF12)%SP%AtmInfo%Atm2Y=AtmInfo%Atm1Y
          AtmPair(CF12)%SP%AtmInfo%Atm2Z=AtmInfo%Atm1Z
       ENDIF



          ENDDO
       ENDDO
       IF(CF12.GT.SIZE(AtmPair)) STOP 'Increase the size of -AtmPair-'
     END SUBROUTINE GetAtomPair2

     SUBROUTINE GetAdrB(I,J,Ind,A,E_O)
!----------------------------------------------------------------------
!H Finds the CSR address corresponding to the dense matrix index pair.
!H Assumes that the collumn index array is in order (binary search)
!H Inputs:
!H   I     = row index of dense matrix
!H   J     = collumn index of dense matrix
!H   IType = 0 halt if element not found
!H         = 1 return index 0 if not found
!H   A     = A sparse matrix
!H   E_O   = 0 call Halt if I cant find the address
!H         = 1 return Ind=0 if I cant find the address
!H Outputs:
!H   Ind   = sparse matrix index corresponding to (I,J)
!----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN)          :: I,J
    INTEGER,INTENT(OUT)         :: Ind
    TYPE(BCSR),INTENT(IN)       :: A
    INTEGER,OPTIONAL,INTENT(IN) :: E_O
    INTEGER                     :: iMid,iBgn,iEnd,Err
    iBgn=A%RowPt%I(I)
    iEnd=A%RowPt%I(I+1)-1
    Ind=0
    Err=0
    IF (PRESENT(E_O)) Err=E_O
    SELECT CASE (Err)
    CASE (0)
1000   CONTINUE
       iMid=(iEnd+iBgn)/2
       IF (A%ColPt%I(iMid)==J) GOTO 2000
       IF (iBgn>=iEnd) THEN
          WRITE(*,*) "I=",I," J=",J
          CALL Halt(' Matrix element not found in ONX:GetAdrB')
       ENDIF
       IF (A%ColPt%I(iMid)>J) THEN
          iEnd=iMid-1
       ELSE
          iBgn=iMid+1
       ENDIF
       GOTO 1000
2000   CONTINUE
       Ind=iMid
    CASE (1)
3000   CONTINUE
       iMid=(iEnd+iBgn)/2
       IF (A%ColPt%I(iMid)==J) GOTO 4000
       IF (iBgn>=iEnd) RETURN
       IF (A%ColPt%I(iMid)>J) THEN
          iEnd=iMid-1
       ELSE
          iBgn=iMid+1
       ENDIF
       GOTO 3000
4000   CONTINUE
       Ind=iMid
    CASE DEFAULT
       CALL Halt(' Illegal switch in ONX:GetAdrB')
    END SELECT
  END SUBROUTINE GetAdrB


  SUBROUTINE GetColIdx(At,D,ColIdx)
!H---------------------------------------------------------------------------------
!H SUBROUTINE GetColIdx(At,D,ColIdx)
!H
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    INTEGER       , INTENT(IN ) :: At
    TYPE(BCSR)    , INTENT(IN ) :: D
    TYPE(INT_VECT), INTENT(OUT) :: ColIdx
    !-------------------------------------------------------------------
    INTEGER                     :: Ci
    !-------------------------------------------------------------------
    !
    CALL INT_VECT_EQ_INT_SCLR(NAtoms,ColIdx%I(1),0)
    !
    DO Ci=D%RowPt%I(At),D%RowPt%I(At+1)-1
       ColIdx%I(D%ColPt%I(Ci))=Ci
    ENDDO
    !
  END SUBROUTINE GetColIdx

END MODULE GetAtomPairMod
   
   
