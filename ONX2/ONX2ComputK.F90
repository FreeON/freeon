MODULE ONX2ComputK
!H=================================================================================
!H MODULE ONX2ComputK
!H This MODULE contains:
!H  PUBLIC:
!H  o SUB 
!H
!H  PRIVATE:
!H  o SUB 
!H
!H  OPTIONS:
!H  DEBUGING: Use -DONX2_DBUG to print some stuff.
!H  INFO    : Use -DONX2_INFO to print some stuff.
!H
!H  Comments:
!H
!H---------------------------------------------------------------------------------
  !
!#define ONX2_DBUG
#ifdef ONX2_DBUG
#define ONX2_INFO
#endif
  !
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE ONX2DataType
  USE InvExp
  !
  IMPLICIT NONE
  !PRIVATE
  !
!--------------------------------------------------------------------------------- 
! PUBLIC DECLARATIONS
!--------------------------------------------------------------------------------- 
!  PUBLIC  :: ComputK
  !
!--------------------------------------------------------------------------------- 
! PRIVATE DECLARATIONS
!--------------------------------------------------------------------------------- 


  !
  !/scratch1/valeryw/MondoExec/ONX2 h2_ultrasmall_G1_12650159 FockBuild 1 1 1 0 1 1
  !/scratch1/valeryw/MondoExec/ONX2 h2_ultrasmall_G2_12650502 FockBuild 2 1 1 1 1 1
  !
  !/scratch1/valeryw/MondoExec/ONX2 h2_small_G1_14204734 FockBuild 1 1 1 0 1 1
  !/scratch1/valeryw/MondoExec/ONX2 h2_small_G2_14204833 FockBuild 1 1 1 0 1 1
  !/scratch1/valeryw/MondoExec/ONX2 h2_small_G3_14204863 FockBuild 1 1 1 0 1 1
  !/scratch1/valeryw/MondoExec/ONX2 h2_small_G4_14207257 FockBuild 1 1 1 0 1 1
  !
  !                                                                                            ONX2   ONX
  !/scratch1/valeryw/MondoExec/ONX2 h2_2349837 FockBuild 1 1 1 0 1 1
  !/scratch1/valeryw/MondoExec/ONX2 h2_2347151 FockBuild 1 1 1 0 1 1
  !/scratch1/valeryw/MondoExec/ONX2 h2_small_perio_15408107 FockBuild 1 1 1 0 1 1   STO-6G
  !/scratch1/valeryw/MondoExec/ONX2 h2_small_15405733 FockBuild 1 1 1 0 1 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !periodic
  !/scratch1/valeryw/MondoExec/ONX2 h2_per_5253061 FockBuild 1 1 1 0 1 1         TZ           12.58   18.73
  !/scratch1/valeryw/MondoExec/ONX2 h2_per_5253091 FockBuild 1 1 1 0 1 1         STO-6G
  !/scratch1/valeryw/MondoExec/ONX2 h2_per_big_5253209 FockBuild 1 1 1 0 1 1     TZ         1742.07 2573.15
  !
  !/scratch1/valeryw/MondoExec/ONX2 h2_small_P1_11022061 FockBuild 2 1 1 1 1 1   TZ
  !/scratch1/valeryw/MondoExec/ONX2 h2_small_P2_11022228 FockBuild 1 1 1 0 1 1   STO-6G
  !/scratch1/valeryw/MondoExec/ONX2 h2_small_P3_11023886 FockBuild 1 1 1 0 1 1   STO-2G
  !/scratch1/valeryw/MondoExec/ONX2 h2_small_P4_11026235 FockBuild 1 1 1 0 1 1   STO-6G H9
  !/scratch1/valeryw/MondoExec/ONX2 h2_small_P5_11028291 FockBuild 1 1 1 0 1 1   STO-6G H32
  !
  !/scratch1/valeryw/MondoExec/ONX2 h2_big_P1_11031381 FockBuild 1 1 1 0 1 1     STO-2G       
  !/scratch1/valeryw/MondoExec/ONX2 h2_big_P2_11031496 FockBuild 1 1 1 0 1 1     STO-3G      942.13  882.70
  !/scratch1/valeryw/MondoExec/ONX2 h2_big_P3_11031710 FockBuild 1 1 1 0 1 1     STO-6G      
  !/scratch1/valeryw/MondoExec/ONX2 h2_big_P4_1139814 FockBuild 1 1 1 0 1 1     6-31G**     4408.99 5510.54 
  !
  !
  !/scratch1/valeryw/MondoExec/ONX2 h2_big_3D_P1_11031971 FockBuild 1 1 1 0 1 1  STO-6G      639.53  687.63
  !/scratch1/valeryw/MondoExec/ONX2 h2_big_3D_P2_5800007 FockBuild 1 1 1 0 1 1   STO-3G      146.32  195.34
  !/scratch1/valeryw/MondoExec/ONX2 h2_big_3D_P3_5800526 FockBuild 1 1 1 0 1 1   STO-3G      813.33 1057.54
  !/scratch1/valeryw/MondoExec/ONX h2_big_3D_P4_1139794 GuessEqCore 0 1 1 0 1 1  6-31G**    1133.20 1971.64
  !
  !/scratch1/valeryw/MondoExec/ONX2 h2_big_1_1652013 DensitySuperposition 0 1 1 0 1 1 STO-6G 2098.97 1961.22
CONTAINS
  !
  !
  SUBROUTINE ComputK(D,Kx,ListC,ListD,GM,BS,CS_OUT)
!H---------------------------------------------------------------------------------
!H SUBROUTINE ComputK(D,Kx,ListC,ListD,GM,BS,CS_OUT)
!H
!H---------------------------------------------------------------------------------
    !
    USE ONXGet, ONLY: GetAdrB
    !
    IMPLICIT NONE
    !
    !
    ! Kab = Dcd*(ac(R)|bd(R'))
    !
    !-------------------------------------------------------------------
    TYPE(BCSR)                 , INTENT(IN   ) :: D
    TYPE(BCSR)                 , INTENT(INOUT) :: Kx
    TYPE(CList2) , DIMENSION(:), POINTER       :: ListC,ListD
    TYPE(CRDS)                 , INTENT(IN   ) :: GM
    TYPE(BSET)                 , INTENT(IN   ) :: BS
    TYPE(CellSet)              , INTENT(IN   ) :: CS_OUT
    !-------------------------------------------------------------------
    TYPE(ANode2), POINTER      :: AtAListTmp,AtAList,AtBListTmp,AtBList
    TYPE(AtomInfo)             :: ACAtmInfo,BDAtmInfo
    INTEGER                    :: AtA,AtB,AtC,AtD,KA,KB,KC,KD,CFA,CFB,CFC,CFD
    INTEGER                    :: ci,iPtrD,iPtrK,NBFC,NBFD,NBFA,NBFB
    INTEGER                    :: CFAC,CFBD
    INTEGER                    :: Off,Ind
    REAL(DOUBLE)               :: Dcd
    REAL(DOUBLE), DIMENSION(12000) :: C,CC  ! this should be declarated somewhere
    INTEGER      :: LocNInt
    REAL(DOUBLE) :: NInts
    !-------------------------------------------------------------------
    TYPE(AtomPr), DIMENSION(50) :: ACAtmPair,BDAtmPair  ! this should be declarated somewhere
    TYPE(OffSt) :: OffSet
    !-------------------------------------------------------------------
    REAL(DOUBLE), EXTERNAL :: DGetAbsMax
    !-------------------------------------------------------------------
    INTEGER :: ACR,BDR
    real(DOUBLE) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
    integer :: i,iA,iB
    real(DOUBLE) :: time1,time2,SumInt
    integer :: MaxCont,IntType

    integer :: iint,iprint,isize
    !
    REAL(DOUBLE), PARAMETER :: ThresholdTwoE=-1.0D-15
    !-------------------------------------------------------------------
    !
    NULLIFY(AtAListTmp,AtAList,AtBListTmp,AtBList)
    !
    !Simple check Simple check Simple check Simple check
    isize=0
    do i=1,natoms
       isize=MAX(isize,BS%BfKnd%I(GM%AtTyp%I(i)))
       IF((BS%BfKnd%I(GM%AtTyp%I(i)))**4.GT.SIZE(C)) THEN
          write(*,*) 'size',(BS%BfKnd%I(GM%AtTyp%I(i)))**4
          STOP 'Incrase the size of C'
       ENDIF
    enddo

    write(*,*) 'size C=',isize**4

    if(CS_OUT%NCells.GT.SIZE(ACAtmPair)) then
       write(*,*) 'size(ACAtmPair)',size(ACAtmPair),'.LT.',CS_OUT%NCells
       STOP 'Incrase the size of ACAtmPair and BDAtmPair'
    endif
    !Simple check Simple check Simple check Simple check


    !
    SumInt=0.0d0
    MaxCont=0
    iint=0
    !
    LocNInt=0
    NInts=0.0d0
    !
    DO AtC=1,NAtoms ! Run over AtC.
       KC=GM%AtTyp%I(AtC)
       NBFC=BS%BfKnd%I(KC)
       ACAtmInfo%Atm2X=GM%Carts%D(1,AtC)
       ACAtmInfo%Atm2Y=GM%Carts%D(2,AtC)
       ACAtmInfo%Atm2Z=GM%Carts%D(3,AtC)
       ACAtmInfo%K2=KC
       !
       ! Get AtA List.
       AtAListTmp=>ListC(AtC)%GoList
       !
       DO ci=D%RowPt%I(AtC),D%RowPt%I(AtC+1)-1 ! Run over AtD
          AtD = D%ColPt%I(ci)
          iPtrD= D%BlkPt%I(ci)
          KD=GM%AtTyp%I(AtD)
          NBFD=BS%BfKnd%I(KD)
          BDAtmInfo%Atm2X=GM%Carts%D(1,AtD)
          BDAtmInfo%Atm2Y=GM%Carts%D(2,AtD)
          BDAtmInfo%Atm2Z=GM%Carts%D(3,AtD)
          BDAtmInfo%K2=KD
          !
          ! Get max of the block density matrix.
          Dcd=DGetAbsMax(NBFC*NBFD,D%MTrix%D(iPtrD))
#ifdef ONX2_DBUG
          WRITE(*,*) 'Max(Dcd)=',Dcd
#endif
          !
          ! Get AtB List.
          AtBListTmp=>ListD(AtD)%GoList
          !
          AtAList=>AtAListTmp
          !
          DO ! Run over AtA
             AtA=AtAList%Atom
             KA=GM%AtTyp%I(AtA)
             NBFA=BS%BfKnd%I(KA)
             !
             !Skip out Dcd*sqrt(ac(R)|ac(R))*sqrt(b_1d(0)|b_1d(0))<Thresh
             !write(*,*) '-----'
!!$             IF(CS_OUT%NCells.EQ.1) THEN
!!$                !write(*,'(A,2I4,E22.15,A,2I4,E22.15,E22.15)') '(ac|ac)',AtA,AtC,AtAList%SqrtInt(1),&
!!$                !     &                                ' (bd|bd)',1,AtD,AtBListTmp%SqrtInt(1),&
!!$                !     &                                 AtAList%SqrtInt(1)*AtBListTmp%SqrtInt(1)
!!$                ACAtmInfo%NCell=1
!!$                !IF(AtAList%SqrtInt(1)*Dcd*AtBListTmp%SqrtInt(1).LT.ThresholdTwoE) THEN
!!$                IF(AtAList%SqrtInt(1)*AtBListTmp%SqrtInt(1).LT.ThresholdTwoE) THEN
!!$                   write(*,'(A,E22.15,A,E22.15)') 'We skip 1', &
!!$                        !& AtAList%SqrtInt(1)*Dcd*AtBListTmp%SqrtInt(1), &
!!$                        & AtAList%SqrtInt(1)*AtBListTmp%SqrtInt(1), &
!!$                        & '.LT.',ThresholdTwoE  
!!$                   EXIT
!!$                ENDIF
!!$             ELSE
             !ACAtmInfo%NCell=GetNonNeglCell(AtAList,AtBListTmp%SqrtInt(1),Thresholds%TwoE)
             ACAtmInfo%NCell=GetNonNeglCell(AtAList,AtBListTmp%SqrtInt(1),ThresholdTwoE)
             IF(ACAtmInfo%NCell.EQ.0) THEN
                !write(*,'(A,E22.15,A,E22.15)') 'We skip 1', &
                !     !& AtAList%SqrtInt(1)*Dcd*AtBListTmp%SqrtInt(1), &
                !     & AtAList%SqrtInt(1)*AtBListTmp%SqrtInt(1), &
                !     & '.LT.',ThresholdTwoE  
                EXIT
             ENDIF
!!$             ENDIF
             !
             ACAtmInfo%Atm1X=GM%Carts%D(1,AtA)
             ACAtmInfo%Atm1Y=GM%Carts%D(2,AtA)
             ACAtmInfo%Atm1Z=GM%Carts%D(3,AtA)
             ACAtmInfo%K1=KA
             !
             ACAtmInfo%Atm12X=ACAtmInfo%Atm1X-ACAtmInfo%Atm2X
             ACAtmInfo%Atm12Y=ACAtmInfo%Atm1Y-ACAtmInfo%Atm2Y
             ACAtmInfo%Atm12Z=ACAtmInfo%Atm1Z-ACAtmInfo%Atm2Z
             !
             ! Get atom pair for BD.
             CALL GetAtomPair2N(ACAtmInfo,AtAList,ACAtmPair,BS,CS_OUT)
             !
             AtBList=>AtBListTmp
             !
             DO ! Run over AtB
                AtB=AtBList%Atom 
                !
                !if(AtB.eq.1.and.AtA.eq.2) then
                IF(AtB.LE.AtA) THEN ! Symmetry of the K matrix
                   !
                   !Skip out |Dcd|*sqrt(ac(0)|ac(0))*sqrt(bd(R')|bd(R'))<Thresh
                   !write(*,*) '-----'
!!$                   IF(CS_OUT%NCells.EQ.1) THEN
!!$                      !write(*,'(A,2I4,E22.15,A,2I4,E22.15,E22.15)') &
!!$                      !&                                ' (ac|ac)',AtA,AtC,AtAList%SqrtInt(1),&
!!$                      !&                                ' (bd|bd)',AtB,AtD,AtBListTmp%SqrtInt(1),&
!!$                      !&                                 AtAList%SqrtInt(1)*AtBListTmp%SqrtInt(1)
!!$                      !IFasaa(AtAList%SqrtInt(1)*Dcd*AtBListTmp%SqrtInt(1).LT.ThresholdTwoE) THEN
!!$                      IF(AtAList%SqrtInt(1)*AtBList%SqrtInt(1).LT.ThresholdTwoE) THEN
!!$                         write(*,'(A,E22.15,A,E22.15)') 'We skip 2', &
!!$                              !& AtAList%SqrtInt(1)*Dcd*AtBListTmp%SqrtInt(1), &
!!$                              & AtAList%SqrtInt(1)*AtBList%SqrtInt(1), &
!!$                              & '.LT.',ThresholdTwoE  
!!$                         EXIT
!!$                      ENDIF
!!$                   ELSE
                   !BDAtmInfo%NCell=GetNonNeglCell(AtBList,AtAList%SqrtInt(1),Thresholds%TwoE)
                   BDAtmInfo%NCell=GetNonNeglCell(AtBList,AtAList%SqrtInt(1),ThresholdTwoE)
                   IF(BDAtmInfo%NCell.EQ.0) THEN
!                         write(*,'(A,E22.15,A,E22.15)') 'We skip 2', &
!                              !& AtAList%SqrtInt(1)*Dcd*AtBListTmp%SqrtInt(1), &
!                              & AtAList%SqrtInt(1)*AtBList%SqrtInt(1), &
!                              & '.LT.',ThresholdTwoE  
                      EXIT
                   ENDIF
!!$                   ENDIF
                   !
                   KB=GM%AtTyp%I(AtB)
                   NBFB=BS%BfKnd%I(KB)
                   !
                   BDAtmInfo%Atm1X=GM%Carts%D(1,AtB)
                   BDAtmInfo%Atm1Y=GM%Carts%D(2,AtB)
                   BDAtmInfo%Atm1Z=GM%Carts%D(3,AtB)
                   BDAtmInfo%K1=KB
                   !
                   BDAtmInfo%Atm12X=BDAtmInfo%Atm1X-BDAtmInfo%Atm2X
                   BDAtmInfo%Atm12Y=BDAtmInfo%Atm1Y-BDAtmInfo%Atm2Y
                   BDAtmInfo%Atm12Z=BDAtmInfo%Atm1Z-BDAtmInfo%Atm2Z
                   !
                   ! Get address for Kx.
                   CALL GetAdrB(AtA,AtB,Ind,Kx,0)
                   iPtrK = Kx%BlkPt%I(Ind)
                   !
                   ! Get atom pair for BD.
                   CALL GetAtomPair2N(BDAtmInfo,AtBList,BDAtmPair,BS,CS_OUT)
                   !
                   !CALL DBL_VECT_EQ_DBL_SCLR(NBFA*NBFB*NBFC*NBFD,C(1),BIG_DBL)
                   CALL DBL_VECT_EQ_DBL_SCLR(NBFA*NBFB*NBFC*NBFD,C(1),0.0d0)
                   !C=BIG_DBL !TO REMOVE
                   !
                   CFAC=0
!!!!!!!!!!!!!!!!!!!!
                   DO ACR=1,ACAtmInfo%NCell
!!!!!!!!!!!!!!!!!!!!
                   OffSet%A=1
                   DO CFA=1,BS%NCFnc%I(KA) ! Run over blkfunc on A
                      OffSet%C=1
                      DO CFC=1,BS%NCFnc%I(KC) ! Run over blkfunc on C
                         CFAC=CFAC+1
                         MaxCont=MAX(MaxCont,ACAtmPair(CFAC)%SP%L)
                         CFBD=0
!!!!!!!!!!!!!!!!!!!!!!!!!!
                         DO BDR=1,BDAtmInfo%NCell
!!!!!!!!!!!!!!!!!!!!!!!!!!
                         OffSet%B=1
                         DO CFB=1,BS%NCFnc%I(KB) ! Run over blkfunc on B
                            OffSet%D=1
                            DO CFD=1,BS%NCFnc%I(KD) ! Run over blkfunc on D
                               CFBD=CFBD+1
                               !
                               ! Compute integral type.
                               IntType=ACAtmPair(CFAC)%SP%IntType*10000+BDAtmPair(CFBD)%SP%IntType
                               !
                               ! The integral interface.
                               INCLUDE 'ERIInterface.Inc'
                               !
                               NInts=NInts+DBLE(LocNInt)
                               OffSet%D=OffSet%D+BS%LStop%I(CFD,KD)-BS%LStrt%I(CFD,KD)+1
                            ENDDO ! End blkfunc on B and D
                            OffSet%B=OffSet%B+BS%LStop%I(CFB,KB)-BS%LStrt%I(CFB,KB)+1
                         ENDDO ! End blkfunc on A and C
!!!!!!!!!!!!!!!!!!!!!!!!!!
                         ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!
                         OffSet%C=OffSet%C+BS%LStop%I(CFC,KC)-BS%LStrt%I(CFC,KC)+1
                      ENDDO
                      OffSet%A=OffSet%A+BS%LStop%I(CFA,KA)-BS%LStrt%I(CFA,KA)+1
                   ENDDO
!!!!!!!!!!!!!!!!!!!!
                   ENDDO
!!!!!!!!!!!!!!!!!!!!
                   !CALL PrintMatrix(C(1),NBFA*NBFB,NBFC*NBFD,2,TEXT_O='Int matrix')
                   !CALL PrintMatrix(C(1),1,NBFA*NBFB*NBFC*NBFD,3,TEXT_O='Int matrix')
                   !CALL PrintMatrix(D%MTrix%D(iPtrD),NBFC,NBFD,2,TEXT_O='D matrix')
#ifdef ONX2_DBUG
                   WRITE(*,'(A,E22.15,4I4)') ' MaxInt=',MAXVAL(C(1:NBFA*NBFB*NBFC*NBFD)),AtA,AtC,AtB,AtD
#endif
                   !
                   ! Digest the block of integral.
                   CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,-1.0d0,C(1), &
                        &     NBFA*NBFB,D%MTrix%D(iPtrD),1,1.0d0, &
                        &     Kx%MTrix%D(IPtrK),1)

                   !CALL PrintMatrix(Kx%MTrix%D(IPtrK),NBFA,NBFB,2,TEXT_O='Int matrix')
                   !CALL Print_BCSR(Kx,'Kx',Unit_O=6)
                   !
                ENDIF
                !
!             endif
                !
                IF(.NOT.ASSOCIATED(AtBList%AtmNext)) EXIT
                AtBList=>AtBList%AtmNext
             ENDDO ! End AtB
             !
             IF(.NOT.ASSOCIATED(AtAList%AtmNext)) EXIT
             AtAList=>AtAList%AtmNext
          ENDDO ! End AtA
          !
       ENDDO ! End AtD
       !
    ENDDO ! End AtC
    !
!#ifdef ONX2_INFO
    WRITE(*,*) '-------------------------------------'
    WRITE(*,*) 'ComputK Statistic.'
    WRITE(*,'(A,F22.1)') ' Nbr ERI  =',NInts
    WRITE(*,'(A,I4)') ' Max Prim =',INT(SQRT(DBLE(MaxCont)))
    WRITE(*,*) '-------------------------------------'
!#endif
    !
  END SUBROUTINE ComputK
  !
  !
  INTEGER FUNCTION GetNonNeglCell(List,DFac,Trsh)
!H---------------------------------------------------------------------------------
!H INTEGER FUNCTION GetNonNeglCell(List,DFac,Trsh)
!H
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(ANode2), POINTER    :: List
    REAL(DOUBLE), INTENT(IN) :: DFac,Trsh
    !-------------------------------------------------------------------
    INTEGER                  :: I
    !-------------------------------------------------------------------
    !
    GetNonNeglCell=0
    !
    DO I=1,List%NCell
       IF(List%SqrtInt(I)*DFac.LE.Trsh) EXIT
       GetNonNeglCell=GetNonNeglCell+1
    ENDDO
    !
  END FUNCTION GetNonNeglCell
  !
  !
  SUBROUTINE GetAtomPair2N(AtmInfo,List,AtmPair,BS,CS_OUT)
!H---------------------------------------------------------------------------------
!H SUBROUTINE GetAtomPair2N(AtmInfo,List,AtmPair,BS,CS_OUT)
!H
!H---------------------------------------------------------------------------------
    USE Thresholding
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(AtomInfo)              , INTENT(IN   ) :: AtmInfo
    TYPE(ANode2  ), POINTER                     :: List
    TYPE(AtomPr  ), DIMENSION(:), INTENT(INOUT) :: AtmPair
    TYPE(BSET)                  , INTENT(IN   ) :: BS
    TYPE(CellSet)               , INTENT(IN   ) :: CS_OUT
    !-------------------------------------------------------------------
    INTEGER      :: CF12,CF1,CF2,I1,I2,II,JJ,IJ,iCell,Cell
    INTEGER      :: MinL1,MaxL1,Type1,MinL2,MaxL2,Type2
    INTEGER      :: StopL1,StartL1,StopL2,StartL2
    REAL(DOUBLE) :: Z1,Z2,Expt,InvExpt,R12,XiR12,RX,RY,RZ
    !-------------------------------------------------------------------
    !
    ! dbuging
    AtmPair(:)%SP%IntType=BIG_INT
    !
    
    CF12=0
    DO iCell=1,AtmInfo%NCell
       Cell=List%CellIdx(iCell)
       RX=CS_OUT%CellCarts%D(1,Cell)
       RY=CS_OUT%CellCarts%D(2,Cell)
       RZ=CS_OUT%CellCarts%D(3,Cell)
       !
       ! AtmInfo must be related to the atoms in the working cell ONLY. 
       ! Then we add the PBC's to have the right interatomic distance.
       R12=(AtmInfo%Atm12X-RX)**2+(AtmInfo%Atm12Y-RY)**2+(AtmInfo%Atm12Z-RZ)**2
       !
       DO CF1=1,BS%NCFnc%I(AtmInfo%K1)
          MinL1=BS%ASymm%I(1,CF1,AtmInfo%K1)
          MaxL1=BS%ASymm%I(2,CF1,AtmInfo%K1)
          Type1=MaxL1*(MaxL1+1)/2+MinL1+1
          StartL1=BS%LStrt%I(CF1,AtmInfo%K1)
          StopL1=BS%LStop%I(CF1,AtmInfo%K1)

          if(Type1==2) stop 'SP shell not yet supported.'

          DO CF2=1,BS%NCFnc%I(AtmInfo%K2)
             CF12=CF12+1
             !
             MinL2=BS%ASymm%I(1,CF2,AtmInfo%K2)
             MaxL2=BS%ASymm%I(2,CF2,AtmInfo%K2)
             Type2=MaxL2*(MaxL2+1)/2+MinL2+1
             StartL2=BS%LStrt%I(CF2,AtmInfo%K2)
             StopL2=BS%LStop%I(CF2,AtmInfo%K2)
             !
             AtmPair(CF12)%SP%IntType=Type1*100+Type2
             !
             II=0
             !
             !old          DO iCell=1,AtmInfo%NCell
             !old             Cell=List%CellIdx(iCell)
             !old             RX=CS_OUT%CellCarts%D(1,Cell)
             !old             RY=CS_OUT%CellCarts%D(2,Cell)
             !old             RZ=CS_OUT%CellCarts%D(3,Cell)
             !old             !
             !old             ! AtmInfo must be related to the atoms in the working cell ONLY. 
             !old             ! Then we add the PBC's to have the right interatomic distance.
             !old             R12=(AtmInfo%Atm12X-RX)**2+(AtmInfo%Atm12Y-RY)**2+(AtmInfo%Atm12Z-RZ)**2
             !
             ! We assume the primitives are ordered (exponants in decressing order).
             DO I1=BS%NPFnc%I(CF1,AtmInfo%K1),1,-1
                Z1=BS%Expnt%D(I1,CF1,AtmInfo%K1)
                JJ=0
                !
                ! We assume the primitives are ordered (exponants in decressing order).
                DO I2=BS%NPFnc%I(CF2,AtmInfo%K2),1,-1
                   Z2=BS%Expnt%D(I2,CF2,AtmInfo%K2)
                   Expt=Z1+Z2
                   InvExpt=1.0d0/Expt
                   XiR12=Z2*Z1*InvExpt*R12
                   IF(XiR12<PrimPairDistanceThreshold) THEN
                      JJ=JJ+1
                      IJ=JJ+II
                      AtmPair(CF12)%SP%Cst(1,IJ)=Expt
                      !
                      ! AtmInfo must be related to the atoms in the working cell ONLY. 
                      ! Then we add the PBC's to have the right atomic position.
                      AtmPair(CF12)%SP%Cst(2,IJ)=(Z1*AtmInfo%Atm1X+Z2*(AtmInfo%Atm2X+RX))*InvExpt
                      AtmPair(CF12)%SP%Cst(3,IJ)=(Z1*AtmInfo%Atm1Y+Z2*(AtmInfo%Atm2Y+RY))*InvExpt
                      AtmPair(CF12)%SP%Cst(4,IJ)=(Z1*AtmInfo%Atm1Z+Z2*(AtmInfo%Atm2Z+RZ))*InvExpt
                      !AtmPair(CF12)%SP%Cst(5,IJ)=5.914967172796D0*EXP(-XiR12)*InvExpt* &
                      AtmPair(CF12)%SP%Cst(5,IJ)=5.914967172796D0*EXPInv(XiR12)*InvExpt* &
                           &                     BS%CCoef%D(StopL1,I1,CF1,AtmInfo%K1)* &
                           &                     BS%CCoef%D(StopL2,I2,CF2,AtmInfo%K2)
                      !
                      ! TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
                      ! Here I need to add the correction factor for SP shell.
                      ! TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
                   ELSE
                      ! We can skipp out the loop because the primitives are ordered.
                      EXIT
                   ENDIF
                ENDDO
                ! We can skipp out the loop if we did not get any significant primitives.
                IF(JJ.EQ.0) EXIT
                II=II+JJ
             ENDDO
             !
             !old          ENDDO ! NCell
             !
             AtmPair(CF12)%SP%L=II
             !AtmPair(CF12)%SP%L(iCell)=II
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
             !
#ifdef ONX2_DBUG
             WRITE(*,'(3(A,I3),A,I5)') 'Type1',Type1,' Type2',Type2, &
                  &                    ' IntType',AtmPair(CF12)%SP%IntType
#endif
             !
          ENDDO
       ENDDO
    ENDDO ! NCell
    !
  END SUBROUTINE GetAtomPair2N
  !
  !
  SUBROUTINE GetAtomPair4(AtmInfo,AtmPair,BS,PBC)
!H---------------------------------------------------------------------------------
!H SUBROUTINE GetAtomPair4(AtmInfo,AtmPair,BS,PBC)
!H
!H---------------------------------------------------------------------------------
    USE Thresholding
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(AtomInfo)              , INTENT(IN   ) :: AtmInfo
    TYPE(AtomPr)  , DIMENSION(:), INTENT(INOUT) :: AtmPair
    TYPE(BSET)                  , INTENT(IN   ) :: BS
    REAL(DOUBLE)  , DIMENSION(3), INTENT(IN   ) :: PBC
    !-------------------------------------------------------------------
    INTEGER      :: CF12,CF1,CF2,I1,I2,II,JJ,IJ,iCell
    INTEGER      :: MinL1,MaxL1,Type1,MinL2,MaxL2,Type2
    INTEGER      :: StartL1,StopL1,StartL2,StopL2
    REAL(DOUBLE) :: Z1,Z2,Expt,InvExpt,R12,XiR12,RX,RY,RZ
    !-------------------------------------------------------------------
    !
    RX=PBC(1)
    RY=PBC(2)
    RZ=PBC(3)
    !
    ! AtmInfo must be related to the atoms in the working cell ONLY. 
    ! Then we add the PBC's to have the right interatomic distance.
    R12=(AtmInfo%Atm12X-RX)**2+(AtmInfo%Atm12Y-RY)**2+(AtmInfo%Atm12Z-RZ)**2
    !
    CF12=0
    DO CF1=1,BS%NCFnc%I(AtmInfo%K1)
       MinL1=BS%ASymm%I(1,CF1,AtmInfo%K1)
       MaxL1=BS%ASymm%I(2,CF1,AtmInfo%K1)
       Type1=MaxL1*(MaxL1+1)/2+MinL1+1
       StartL1=BS%LStrt%I(CF1,AtmInfo%K1)
       StopL1=BS%LStop%I(CF1,AtmInfo%K1)

       if(Type1==2) stop 'SP shell not yet supported.'
       DO CF2=1,BS%NCFnc%I(AtmInfo%K2)
          CF12=CF12+1
          !
          MinL2=BS%ASymm%I(1,CF2,AtmInfo%K2)
          MaxL2=BS%ASymm%I(2,CF2,AtmInfo%K2)
          Type2=MaxL2*(MaxL2+1)/2+MinL2+1
          AtmPair(CF12)%SP%IntType=Type1*100+Type2
          StartL2=BS%LStrt%I(CF2,AtmInfo%K2)
          StopL2=BS%LStop%I(CF2,AtmInfo%K2)
          !
          II=0
          !
          ! We assume the primitives are ordered (exponants in decressing order).
          DO I1=BS%NPFnc%I(CF1,AtmInfo%K1),1,-1
             Z1=BS%Expnt%D(I1,CF1,AtmInfo%K1)
             JJ=0
             !
             !We assume the primitives are ordered (exponants in decressing order).
             DO I2=BS%NPFnc%I(CF2,AtmInfo%K2),1,-1
                Z2=BS%Expnt%D(I2,CF2,AtmInfo%K2)
                Expt=Z1+Z2
                InvExpt=1.0d0/Expt
                XiR12=Z2*Z1*InvExpt*R12
                IF(XiR12<PrimPairDistanceThreshold) THEN
                   JJ=JJ+1
                   IJ=JJ+II
                   AtmPair(CF12)%SP%Cst(1,IJ)=Expt
                   !
                   ! AtmInfo must be related to the atoms in the working cell ONLY. 
                   ! Then we add the PBC's to have the right atomic position.
                   AtmPair(CF12)%SP%Cst(2,IJ)=(Z1*AtmInfo%Atm1X+Z2*(AtmInfo%Atm2X+RX))*InvExpt
                   AtmPair(CF12)%SP%Cst(3,IJ)=(Z1*AtmInfo%Atm1Y+Z2*(AtmInfo%Atm2Y+RY))*InvExpt
                   AtmPair(CF12)%SP%Cst(4,IJ)=(Z1*AtmInfo%Atm1Z+Z2*(AtmInfo%Atm2Z+RZ))*InvExpt
                   !AtmPair(CF12)%SP%Cst(5,IJ)=5.914967172796D0*EXP(-XiR12)* &
                   AtmPair(CF12)%SP%Cst(5,IJ)=5.914967172796D0*EXPInv(XiR12)*InvExpt* &
                        &                     BS%CCoef%D(StopL1,I1,CF1,AtmInfo%K1)* &
                        &                     BS%CCoef%D(StopL2,I2,CF2,AtmInfo%K2)
                   !
                   ! TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
                   ! Here I need to add the correction factor for SP shell.
                   ! TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
                ELSE
                   ! We can skipp out the loop because the primitives are ordered.
                   EXIT
                ENDIF
             ENDDO
             ! We can skipp out the loop because we did not get any significant primitives.
             IF(JJ.EQ.0) EXIT
             II=II+JJ
             !
          ENDDO
          !
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
          !
#ifdef ONX2_DBUG
          write(*,*) 'Printing from GetAtomPair4'
          write(*,'(3(A,I3),A,I5)') 'Type1',Type1,' Type2',Type2, &
               &                    ' IntType',AtmPair(CF12)%SP%IntType
#endif
          !
       ENDDO
    ENDDO
    !
  END SUBROUTINE GetAtomPair4
  !
  !
  SUBROUTINE MakeList(List,GM,BS,CS_OUT)
!H---------------------------------------------------------------------------------
!H SUBROUTINE MakeList(List,GM,BS,CS_OUT)
!H
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(CList2) , DIMENSION(:), POINTER :: List
    TYPE(CRDS)   , INTENT(IN)            :: GM
    TYPE(BSET)   , INTENT(IN)            :: BS
    TYPE(CellSet), INTENT(IN)            :: CS_OUT
    !-------------------------------------------------------------------
    TYPE(ANode2) , POINTER               :: AtAList,AtAListTmp,NodeA
    TYPE(AtomInfo)                       :: ACAtmInfo
    INTEGER                              :: AtA,AtC,KA,KC,CFA,CFC,iCell,CFAC
    INTEGER                              :: NCell,I,IntType,LocNInt
    REAL(DOUBLE)                         :: RInt,AC2,NInts
    !-------------------------------------------------------------------
    TYPE(AtomPr) , DIMENSION(50)         :: ACAtmPair ! this should be declared somewhere
    REAL(DOUBLE) , DIMENSION(50)         :: RIntCell  ! this should be declared somewhere
    REAL(DOUBLE) , DIMENSION(12000)      :: C         ! this should be declared somewhere
    INTEGER      , DIMENSION(50)         :: IndxCell  ! this should be declared somewhere
    !-------------------------------------------------------------------
    REAL(DOUBLE) , EXTERNAL              :: DGetAbsMax
    !-------------------------------------------------------------------
    REAL(DOUBLE), PARAMETER :: ThresholdIntegral=-1.0D-15
    REAL(DOUBLE), PARAMETER :: ThresholdDistance=1.0D+99
    !-------------------------------------------------------------------
    !
    integer :: isize
    !Simple check
    isize=0
    do i=1,natoms
       isize=MAX(isize,BS%BfKnd%I(GM%AtTyp%I(i)))
       IF((BS%BfKnd%I(GM%AtTyp%I(i)))**4.GT.SIZE(C)) THEN
          write(*,*) 'size',(BS%BfKnd%I(GM%AtTyp%I(i)))**4
          STOP 'Incrase the size of C'
       ENDIF
    enddo
    write(*,*) 'size C=',isize**4
    !
    !
    NULLIFY(AtAList,AtAListTmp,NodeA)
    NInts=0.0d0
    !
    DO AtC=1,NAtoms ! Run over AtC
       !
       KC=GM%AtTyp%I(AtC)
       !
       ACAtmInfo%Atm2X=GM%Carts%D(1,AtC)
       ACAtmInfo%Atm2Y=GM%Carts%D(2,AtC)
       ACAtmInfo%Atm2Z=GM%Carts%D(3,AtC)
       ACAtmInfo%K2=KC
       !
       DO AtA=1,NAtoms ! Run over AtA
          !
          KA=GM%AtTyp%I(AtA)
          !
          ACAtmInfo%Atm1X=GM%Carts%D(1,AtA)
          ACAtmInfo%Atm1Y=GM%Carts%D(2,AtA)
          ACAtmInfo%Atm1Z=GM%Carts%D(3,AtA)
          ACAtmInfo%K1=KA
          !
          ACAtmInfo%Atm12X=ACAtmInfo%Atm1X-ACAtmInfo%Atm2X
          ACAtmInfo%Atm12Y=ACAtmInfo%Atm1Y-ACAtmInfo%Atm2Y
          ACAtmInfo%Atm12Z=ACAtmInfo%Atm1Z-ACAtmInfo%Atm2Z
          !
          ! Check the interatomic distance in the working box.
          AC2=(ACAtmInfo%Atm12X)**2+(ACAtmInfo%Atm12Y)**2+(ACAtmInfo%Atm12Z)**2
          !
          ! Cycle if needed.
          IF(AC2.GT.ThresholdDistance) CYCLE
          !
          ! Initialize some cell variables.
          NCell=0
          !
          DO iCell=1,CS_OUT%NCells ! Run over R
             !
             ! Check the interatomic distance between boxes.
             AC2 =  (ACAtmInfo%Atm12X-CS_OUT%CellCarts%D(1,iCell))**2+ &
                  & (ACAtmInfo%Atm12Y-CS_OUT%CellCarts%D(2,iCell))**2+ &
                  & (ACAtmInfo%Atm12Z-CS_OUT%CellCarts%D(3,iCell))**2
             !
             ! Cycle if needed.
             IF(AC2.GT.ThresholdDistance) CYCLE
             !
             ! Get the atom pair.
             CALL GetAtomPair4(ACAtmInfo,ACAtmPair,BS,CS_OUT%CellCarts%D(1,iCell))
             !
             ! Initialize some cell variables.
             RInt=0.0d0
             !
             !C=BIG_DBL !TO REMOVE
             !
             DO CFAC=1,BS%NCFnc%I(KA)*BS%NCFnc%I(KC) ! Run over blkfunc on A,C
                !
                ! Compute integral type.
                IntType=ACAtmPair(CFAC)%SP%IntType
                !
                ! The integral interface.
                INCLUDE 'ERIListInterface.Inc'
                !
                RInt=MAX(RInt,DGetAbsMax(LocNInt,C(1)))
#ifdef ONX2_DBUG
                WRITE(*,'(2(A,E22.15),2(A,I6))') 'RInt',RInt,' RIntLocal',DGetAbsMax(LocNInt,C(1)), &
                     &    ' LocNInt',LocNInt,' IntType',IntType
#endif
                !
                NInts=NInts+DBLE(LocNInt)
             ENDDO ! End over blkfunc on A,C
             !
#ifdef ONX2_DBUG
             WRITE(*,'(A,E22.15,4I4)') ' MaxInt =',RInt,AtA,AtC,AtA,AtC
#endif
             RInt=DSQRT(RInt)
             !
             ! Keep the cell if needed.
             IF(RInt.GT.ThresholdIntegral) THEN
                NCell=NCell+1
                RIntCell(NCell)=RInt
                IndxCell(NCell)=iCell
                !write(*,'(A,I3,A,I3,A,I3,A,E25.15,A,I3)') &
                !    & 'AtA',AtA,' AtC',AtC,' NCell',NCell,' RInt',RInt,' iCell',iCell
             ENDIF
             !
          ENDDO ! End R
          !
          ! Check if no cell.
          IF(NCell.EQ.0) CYCLE
          !
          ! Allocate a new node.
          ALLOCATE(NodeA)
          NULLIFY(NodeA%AtmNext)
          NodeA%Atom=AtA
          NodeA%NCell=NCell
          !
          ! Order the Cell list.
          IF(CS_OUT%NCells.GT.0) CALL QuickSortDis(RIntCell(1),IndxCell(1),NCell,-2)
          !
          ! Allocate the Cell array.
          ALLOCATE(NodeA%CellIdx(NCell))
          ALLOCATE(NodeA%SqrtInt(NCell))
          !
          ! Copy the arrays.
          DO I=1,NCell
             NodeA%SqrtInt(I)=RIntCell(I)
             NodeA%CellIdx(I)=IndxCell(I)
          ENDDO
          !
          ! Insert the new node (at the right place).
          CALL InsertNode(List(AtC)%GoList,NodeA)
          !
       ENDDO ! End AtA
       !
    ENDDO ! End AtC
    !
!#ifdef ONX2_INF
    WRITE(*,*) '-------------------------------------'
    WRITE(*,*) 'MakeList Statistic.'
    WRITE(*,'(A,F22.1)') ' Nbr ERI=',NInts
    WRITE(*,*) '-------------------------------------'
!#endif
    !
  END SUBROUTINE MakeList
  !
  !
  SUBROUTINE InsertNode(BegPtr,NewPtr)
!H---------------------------------------------------------------------------------
!H SUBROUTINE InsertNode(BegPtr,NewPtr)
!H
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(ANode2), POINTER :: BegPtr,NewPtr
    !-------------------------------------------------------------------
    TYPE(ANode2), POINTER :: TmpPtr
    !-------------------------------------------------------------------
    !
    NULLIFY(TmpPtr)
    !
    ! Check if Beg is null, then insert node and exit.
    IF(.NOT.ASSOCIATED(BegPtr)) THEN
       BegPtr=>NewPtr
       NULLIFY(NewPtr)
       RETURN
    ENDIF
    !
    ! Check if the first element in the list is smaller
    ! than the ones in the new node.
    IF(BegPtr%SqrtInt(1).LE.NewPtr%SqrtInt(1)) THEN
       NewPtr%AtmNext=>BegPtr
       BegPtr=>NewPtr       
       NULLIFY(NewPtr)
       RETURN
    ENDIF
    !
    ! General search.
    TmpPtr=>BegPtr
    DO
       ! We are at the end.
       IF(.NOT.ASSOCIATED(TmpPtr%AtmNext)) THEN
          TmpPtr%AtmNext=>NewPtr
          NULLIFY(NewPtr)
          EXIT
       ENDIF
       !
       IF(TmpPtr%AtmNext%SqrtInt(1).LE.NewPtr%SqrtInt(1)) THEN
          NewPtr%AtmNext=>TmpPtr%AtmNext
          TmpPtr%AtmNext=>NewPtr
          NULLIFY(NewPtr)
          EXIT
       ENDIF
       TmpPtr=>TmpPtr%AtmNext
       !
    ENDDO
    !
  END SUBROUTINE InsertNode
  !
  !
  SUBROUTINE AllocList(List,LSize)
!H---------------------------------------------------------------------------------
!H SUBROUTINE AllocList(List,LSize)
!H
!H---------------------------------------------------------------------------------
    !
    TYPE(CList2), DIMENSION(:), POINTER    :: List
    INTEGER                   , INTENT(IN) :: LSize
    !-------------------------------------------------------------------
    INTEGER                                :: I
    !-------------------------------------------------------------------
    !
    ALLOCATE(List(LSize))
    !
    DO I=1,LSize
       NULLIFY(List(I)%GoList)
    ENDDO
    !
  END SUBROUTINE AllocList
  !
  !
  SUBROUTINE DeAllocList(List)
!H---------------------------------------------------------------------------------
!H SUBROUTINE DeAllocList(List)
!H
!H---------------------------------------------------------------------------------
    !
    TYPE(CList2), DIMENSION(:), POINTER :: List
    !-------------------------------------------------------------------
    TYPE(ANode2)              , POINTER :: ListA,ListATmp
    INTEGER                             :: LSize,I
    !-------------------------------------------------------------------
    !
    NULLIFY(ListA,ListATmp)
    !
    LSize=SIZE(List)
    !
    DO I=1,LSize
       !
       ListA=>List(I)%GoList
       !
       DO
          IF(.NOT.ASSOCIATED(ListA)) EXIT
          !
          IF(.NOT.ASSOCIATED(ListA%CellIdx).OR. &
               & .NOT.ASSOCIATED(ListA%SqrtInt)) STOP 'Array not allocate for the list.'
          !
          IF(ASSOCIATED(ListA%CellIdx)) DEALLOCATE(ListA%CellIdx)
          IF(ASSOCIATED(ListA%SqrtInt)) DEALLOCATE(ListA%SqrtInt)
          !
          ListATmp=>ListA
          ListA=>ListA%AtmNext
          !
          DEALLOCATE(ListATmp)
          !
       ENDDO
       NULLIFY(List(I)%GoList)
       !
    ENDDO
    !
    DEALLOCATE(List)
    !
  END SUBROUTINE DeAllocList
  !
  !
  SUBROUTINE PrintList2(List)
!H---------------------------------------------------------------------------------
!H SUBROUTINE PrintList2(List)
!H
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(CList2), DIMENSION(:), POINTER :: List
    !-------------------------------------------------------------------
    TYPE(ANode2)              , POINTER :: ListA,ListATmp
    INTEGER                             :: LSize,I,J
    !-------------------------------------------------------------------
    !
    NULLIFY(ListA,ListATmp)
    !
    LSize=SIZE(List)
    !
    DO I=1,LSize
       !
       ListA=>List(I)%GoList
       !
       DO
          IF(.NOT.ASSOCIATED(ListA)) EXIT
          !
          WRITE(*,'(A,I4,A,I4)') 'AtC',I,', AtA',ListA%Atom
          IF(.NOT.ASSOCIATED(ListA%CellIdx).OR. &
               & .NOT.ASSOCIATED(ListA%SqrtInt)) STOP 'Array not allocate for the list.'
          DO J=1,SIZE(ListA%CellIdx)
             WRITE(*,'(A,I4,A,E22.15)') 'CellIdx=',ListA%CellIdx(J),' SqrtInt=',ListA%SqrtInt(J)
          ENDDO
          !
          ListA=>ListA%AtmNext
          !
       ENDDO
       !
    ENDDO
    !
  END SUBROUTINE PrintList2
  !
  !
  SUBROUTINE PrintMatrix(V,M,N,IOpt,IOut_O,SHFTM_O,SHFTN_O,TEXT_O)
    IMPLICIT NONE
!    ++++ PRINT OUT A RECTANGULAR MATRIX +++++
!    
!    V    : MATRIX MxN
!    M    : NUMBER OF ROW
!    N    : NUMBER OF COLUMN
!
!    IOpt = 0; MAX COLUMN = 10 
!    IOpt = 1; MAX COLUMN =  7
!    IOpt = 2; MAX COLUMN =  5
!    IOpt = 3; MAX COLUMN =  3 
!
    INTEGER                                  :: IMax,IMin,NMAX,I,J
    INTEGER                                  :: IOut,SHFTM,SHFTN
    INTEGER                     , INTENT(IN) :: M,N,IOpt
    INTEGER         , OPTIONAL  , INTENT(IN) :: IOut_O,SHFTM_O,SHFTN_O
    REAL(DOUBLE)                             :: V(M,N)
    CHARACTER(LEN=*), OPTIONAL  , INTENT(IN) :: TEXT_O
!
    IOut = 6
    IF(PRESENT(IOut_O)) IOut = IOut_O
!
    SHFTM = 0
    SHFTN = 0
    IF(PRESENT(SHFTM_O)) SHFTM = SHFTM_O
    IF(PRESENT(SHFTN_O)) SHFTN = SHFTN_O
!
    WRITE(IOut,100)
!
    SELECT CASE(IOPT)
    CASE(0); NMAX = 10
    CASE(1); NMAX =  7
    CASE(2); NMAX =  5
    CASE(3); NMAX =  3
    CASE DEFAULT
       ! WRITE ERROR MESSAGE OR PUT SOMETHING THERE
       RETURN
    END SELECT
!
    IF(N.EQ.0) RETURN
! 
    IF(PRESENT(TEXT_O)) WRITE(IOUT,*) TEXT_O
    IMAX = 0
!    
    DO WHILE (IMAX.LT.M)
       IMIN = IMAX+1
       IMAX = IMAX+NMAX
       IF(IMAX .GT. M) IMAX = M
       SELECT CASE(IOPT)
       CASE(0); WRITE(IOUT,1000) (I+SHFTN,I=IMIN,IMAX)
       CASE(1); WRITE(IOUT,2000) (I+SHFTN,I=IMIN,IMAX)
       CASE(2); WRITE(IOUT,3000) (I+SHFTN,I=IMIN,IMAX)
       CASE(3); WRITE(IOUT,4000) (I+SHFTN,I=IMIN,IMAX)
       END SELECT
       DO J = 1,N
          SELECT CASE(IOPT)
          CASE(0); WRITE(IOUT,1100) J+SHFTM,(V(J,I),I=IMIN,IMAX)
          CASE(1); WRITE(IOUT,2100) J+SHFTM,(V(J,I),I=IMIN,IMAX)
          CASE(2); WRITE(IOUT,3100) J+SHFTM,(V(J,I),I=IMIN,IMAX)
          CASE(3); WRITE(IOUT,4100) J+SHFTM,(V(J,I),I=IMIN,IMAX)
          END SELECT
       ENDDO
       WRITE (IOUT,100)
    ENDDO
    !
    WRITE (IOUT,100)
    !    
100 FORMAT('')
1000 FORMAT(6X,10(4X,I3,4X))
1100 FORMAT(I5,1X,10F11.5  )
2000 FORMAT(6X,7(6X,I3,6X ))
2100 FORMAT(I5,1X,7F15.10  )
3000 FORMAT(6X,7(7X,I3,6X ))
3100 FORMAT(I5,1X,7E16.8   )
4000 FORMAT(6X,7(7X,I3,6X ))
4100 FORMAT(I5,1X,7E16.8   )
  END SUBROUTINE PrintMatrix
  !
END MODULE ONX2ComputK

