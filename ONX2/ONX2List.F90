MODULE ONX2List
!H=================================================================================
!H MODULE ONX2List
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
  USE ONXParameters
  !
#ifdef PARALLEL
  USE MondoMPI
  USE FastMatrices
#endif
  !
  IMPLICIT NONE
  !PRIVATE
  !
!--------------------------------------------------------------------------------- 
! PUBLIC DECLARATIONS
!--------------------------------------------------------------------------------- 
  PUBLIC  :: MakeList
  PUBLIC  :: MakeGList
  !
!--------------------------------------------------------------------------------- 
! PRIVATE DECLARATIONS
!--------------------------------------------------------------------------------- 
  PRIVATE :: GetAtomPair
  PRIVATE :: GetAtomPairG_
  !
CONTAINS
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
    INTEGER                              :: NCell,I,IntType,LocNInt,NBFA,NBFC
    REAL(DOUBLE)                         :: RInt,AC2,NInts
    !-------------------------------------------------------------------
    REAL(DOUBLE) , DIMENSION(CS_OUT%NCells) :: RIntCell
    INTEGER      , DIMENSION(CS_OUT%NCells) :: IndxCell
    !-------------------------------------------------------------------
    TYPE(AtomPr) , DIMENSION(100)        :: ACAtmPair ! this should be declared somewhere
    REAL(DOUBLE) , DIMENSION(38416)      :: C         ! this should be declared somewhere
    !-------------------------------------------------------------------
    REAL(DOUBLE) , EXTERNAL              :: DGetAbsMax
    !-------------------------------------------------------------------
    REAL(DOUBLE), PARAMETER :: ThresholdIntegral=-1.0D-15
    !REAL(DOUBLE), PARAMETER :: ThresholdDistance=1.0D+99
    !REAL(DOUBLE), PARAMETER :: ThresholdIntegral=1.0D-12
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
#ifdef PARALLEL
  IF(MyID.EQ.ROOT) &
#endif
!    write(*,*) 'size C=',isize**4
    !
    !
    NULLIFY(AtAList,AtAListTmp,NodeA)
    NInts=0.0d0
    !
!#ifdef PARALLEL
    ! TODO TODO TODO TODO TODO TODO TODO TODO
!#else
    DO AtC=1,NAtoms ! Run over AtC
!#endif
       !
       KC=GM%AtTyp%I(AtC)
       NBFC=BS%BfKnd%I(KC)
       ACAtmInfo%Atm2X=GM%Carts%D(1,AtC)
       ACAtmInfo%Atm2Y=GM%Carts%D(2,AtC)
       ACAtmInfo%Atm2Z=GM%Carts%D(3,AtC)
       ACAtmInfo%K2=KC
       !
       DO AtA=1,NAtoms ! Run over AtA
          !
          KA=GM%AtTyp%I(AtA)
          NBFA=BS%BfKnd%I(KA)
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
          ! Set the range for range of exchange.
          DisRange=MAX(DisRange,SQRT(AC2)*1.01D0)
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
             CALL GetAtomPair(ACAtmInfo,ACAtmPair,BS,CS_OUT%CellCarts%D(1,iCell))
             !
             ! Initialize some cell variables.
             RInt=0.0d0
             !
             DO CFAC=1,BS%NCFnc%I(KA)*BS%NCFnc%I(KC) ! Run over blkfunc on A,C
                !
                ! Compute integral type.
                IntType=ACAtmPair(CFAC)%SP%IntType
                !
                CALL DBL_VECT_EQ_DBL_SCLR(NBFA*NBFA*NBFC*NBFC,C(1),0.0d0) !I need less zeroing!
                !
                ! The integral interface.
                INCLUDE 'ERIListInterface.Inc'
                !
                RInt=MAX(RInt,DGetAbsMax(LocNInt,C(1)))
                !
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
#ifdef ONX2_INFO
#ifdef PARALLEL
  IF(MyID.EQ.ROOT) THEN
#endif
    WRITE(*,*) '-------------------------------------'
    WRITE(*,*) 'MakeList Statistic.'
    WRITE(*,'(A,F22.1)') ' Nbr ERI=',NInts
    WRITE(*,*) '-------------------------------------'
#ifdef PARALLEL
  ENDIF
#endif
#endif
    !
  END SUBROUTINE MakeList
  !
  !
  SUBROUTINE MakeGList(List,GM,BS,CS_OUT)
!H---------------------------------------------------------------------------------
!H SUBROUTINE MakeGList(List,GM,BS,CS_OUT)
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
    INTEGER                              :: NCell,I,IntType,LocNInt,NBFA,NBFC
    REAL(DOUBLE)                         :: RInt,AC2,NInts
    !-------------------------------------------------------------------
    REAL(DOUBLE) , DIMENSION(CS_OUT%NCells) :: RIntCell
    INTEGER      , DIMENSION(CS_OUT%NCells) :: IndxCell
    !-------------------------------------------------------------------
    TYPE(AtomPr) , DIMENSION(   MaxShelPerAtmBlk**2*CS_OUT%NCells) :: ACAtmPair
    REAL(DOUBLE) , DIMENSION(12*MaxFuncPerAtmBlk**4) :: C
    !-------------------------------------------------------------------
    REAL(DOUBLE) , EXTERNAL              :: DGetAbsMax
    !-------------------------------------------------------------------
    REAL(DOUBLE), PARAMETER :: ThresholdIntegral=-1.0D-15
    !REAL(DOUBLE), PARAMETER :: ThresholdDistance=1.0D+99
    !REAL(DOUBLE), PARAMETER :: ThresholdIntegral=1.0D-12
    REAL(DOUBLE), PARAMETER :: ThresholdDistance=1.0D+99
    !-------------------------------------------------------------------
    !
    integer :: isize,NIntBlk
    !Simple check
    isize=0
    do i=1,natoms
       isize=MAX(isize,BS%BfKnd%I(GM%AtTyp%I(i)))
       IF((BS%BfKnd%I(GM%AtTyp%I(i)))**4.GT.SIZE(C)) THEN
          write(*,*) 'size',(BS%BfKnd%I(GM%AtTyp%I(i)))**4
          STOP 'Incrase the size of C'
       ENDIF
    enddo
!#ifdef PARALLEL
!  IF(MyID.EQ.ROOT) &
!#endif
    !write(*,*) 'size C=',isize**4
!    write(*,*) 'In MakeList 1'
    !
    !
    NULLIFY(AtAList,AtAListTmp,NodeA)
    NInts=0.0d0
    !
!#ifdef PARALLEL
    ! TODO TODO TODO TODO TODO TODO TODO TODO
!#else
    DO AtC=1,NAtoms ! Run over AtC
!#endif
       !
       KC=GM%AtTyp%I(AtC)
       NBFC=BS%BfKnd%I(KC)
       !
       ACAtmInfo%Atm2X=GM%Carts%D(1,AtC)
       ACAtmInfo%Atm2Y=GM%Carts%D(2,AtC)
       ACAtmInfo%Atm2Z=GM%Carts%D(3,AtC)
       ACAtmInfo%K2=KC

!    write(*,*) 'In MakeList 2'

       !
       DO AtA=1,NAtoms ! Run over AtA
          !
          KA=GM%AtTyp%I(AtA)
          NBFA=BS%BfKnd%I(KA)
          !
          ACAtmInfo%Atm1X=GM%Carts%D(1,AtA)
          ACAtmInfo%Atm1Y=GM%Carts%D(2,AtA)
          ACAtmInfo%Atm1Z=GM%Carts%D(3,AtA)
          ACAtmInfo%K1=KA

!    write(*,*) 'In MakeList 3'

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

!    write(*,*) 'In MakeList 4'

          ! Set the range for range of exchange.
          DisRange=MAX(DisRange,SQRT(AC2)*1.01D0)
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

!    write(*,*) 'In MakeList 5'

             !
             ! Cycle if needed.
             IF(AC2.GT.ThresholdDistance) CYCLE
             !
             ! Get the atom pair.
             !CALL GetAtomPair(ACAtmInfo,ACAtmPair,BS,CS_OUT%CellCarts%D(1,iCell))
             CALL GetAtomPairG_(ACAtmInfo,ACAtmPair,BS,CS_OUT%CellCarts%D(1,iCell))
             !
             ! Initialize some cell variables.
             RInt=0.0d0
             !

!    write(*,*) 'In MakeList 6'

             CFAC=0
             !DO CFAC=1,BS%NCFnc%I(KA)*BS%NCFnc%I(KC)
             DO CFA=1,BS%NCFnc%I(KA) ! Run over blkfunc on A
             DO CFC=1,BS%NCFnc%I(KC) ! Run over blkfunc on C
                !
                CFAC=CFAC+1
                !
                ! Compute integral type.
                IntType=ACAtmPair(CFAC)%SP%IntType
                !
                LocNInt=(BS%LStop%I(CFA,KA)-BS%LStrt%I(CFA,KA)+1)**2* &
                     &  (BS%LStop%I(CFC,KC)-BS%LStrt%I(CFC,KC)+1)**2 *12
                !
                !CALL DBL_VECT_EQ_DBL_SCLR(LocNInt,C(1),0.0d0) !I need less zeroing *12!
                CALL DBL_VECT_EQ_DBL_SCLR(LocNInt,C(1),0.0d0) !I need less zeroing *12!
                !
                ! The integral interface.
                !write(*,*) 'In'
                INCLUDE 'ERIListInterface.Inc'
                !INCLUDE 'DERIListInterface.Inc'
                !write(*,*) 'Out'
                !

                !write(*,*) 'C',C(1:LocNInt)

                RInt=MAX(RInt,DGetAbsMax(LocNInt,C(1)))
                !
#ifdef GONX2_DBUG
                WRITE(*,'(2(A,E22.15),4(A,I6))') 'RInt',RInt,' RIntLocal',DGetAbsMax(LocNInt,C(1)), &
                     & ' LocNInt',LocNInt,' IntType',IntType,' CFC',CFC,' CFA',CFA
#endif
                !
                NInts=NInts+DBLE(LocNInt)
             ENDDO ! End over blkfunc on A
             ENDDO ! End over blkfunc on C
             !
             !stop 'in onx2list'
             !if(atc==1.and.ata==1) RInt=(0.123456789d0)**2
             !if(atc==2.and.ata==2) RInt=(0.123456789d0)**2
             !if(atc==3.and.ata==3) RInt=(0.123456789d0)**2

#ifdef GONX2_DBUG
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
#ifdef GONX2_INFO
#ifdef PARALLEL
  IF(MyID.EQ.ROOT) THEN
#endif
    WRITE(*,*) '-------------------------------------'
    WRITE(*,*) 'MakeList Statistic.'
    WRITE(*,'(A,F22.1)') ' Nbr ERI=',NInts
    WRITE(*,*) '-------------------------------------'
#ifdef PARALLEL
  ENDIF
#endif
#endif
    !
!  write(*,*) 'End of MakeGList'
  END SUBROUTINE MakeGList
  !
  !
  SUBROUTINE GetAtomPair(AtmInfo,AtmPair,BS,PBC)
!H---------------------------------------------------------------------------------
!H SUBROUTINE GetAtomPair(AtmInfo,AtmPair,BS,PBC)
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
          write(*,*) 'Printing from GetAtomPair'
          write(*,'(3(A,I3),A,I5)') 'Type1',Type1,' Type2',Type2, &
               &                    ' IntType',AtmPair(CF12)%SP%IntType
#endif
          !
       ENDDO
    ENDDO
    !
    IF(CF12.GT.SIZE(AtmPair)) STOP 'Increase the size of -AtmPair-'
    !
  END SUBROUTINE GetAtomPair
  !
  !
  SUBROUTINE GetAtomPairG_(AtmInfo,AtmPair,BS,PBC)
!H---------------------------------------------------------------------------------
!H SUBROUTINE GetAtomPairG(AtmInfo,AtmPair,BS,PBC)
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
    INTEGER      :: ISwitch1,ISwitch2
    REAL(DOUBLE) :: Z1,Z2,Expt,InvExpt,R12,XiR12,RX,RY,RZ
    LOGICAL      :: Switch
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
          ! .NOT.(Do we need to switch the centers?)
          Switch=Type1.GE.Type2
          !
          IF(Switch) THEN
             ISwitch1=6
             ISwitch2=7
          ELSE
             ISwitch1=7
             ISwitch2=6
          ENDIF
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
                   AtmPair(CF12)%SP%Cst(ISwitch1,IJ)=2.0d0*Z1!->ISwitch1=6,7
                   AtmPair(CF12)%SP%Cst(ISwitch2,IJ)=2.0d0*Z2!->ISwitch2=7,6
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
          IF(Switch) THEN
             !old IF(Type1.GE.Type2) THEN
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
          write(*,*) 'Printing from GetAtomPair'
          write(*,'(3(A,I3),A,I5)') 'Type1',Type1,' Type2',Type2, &
               &                    ' IntType',AtmPair(CF12)%SP%IntType
#endif
          !
       ENDDO
    ENDDO
    !
    IF(CF12.GT.SIZE(AtmPair)) STOP 'Increase the size of -AtmPair-'
    !
  END SUBROUTINE GetAtomPairG_
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
  SUBROUTINE PrintList(List)
!H---------------------------------------------------------------------------------
!H SUBROUTINE PrintList(List)
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
  END SUBROUTINE PrintList
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
       DO J = 1,M
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
END MODULE ONX2List

