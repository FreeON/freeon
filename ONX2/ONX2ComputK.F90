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
#ifndef PARALLEL
#undef ONX2_PARALLEL
#endif
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
#ifdef ONX2_PARALLEL
  USE MondoMPI
  USE FastMatrices
#endif
  !
  IMPLICIT NONE
  PRIVATE
  !
!--------------------------------------------------------------------------------- 
! PUBLIC DECLARATIONS
!--------------------------------------------------------------------------------- 
  PUBLIC  :: ComputK
  !
!--------------------------------------------------------------------------------- 
! PRIVATE DECLARATIONS
!--------------------------------------------------------------------------------- 
  PRIVATE :: GetNonNeglCell
  PRIVATE :: GetAtomPair2N
  !
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
#ifdef ONX2_PARALLEL
    TYPE(fastmat)              , POINTER       :: D
    TYPE(fastmat)              , POINTER       :: Kx
    TYPE(FASTMAT)              , POINTER       :: P
    TYPE(SRST   )              , POINTER       :: U

    TYPE(FASTMAT)              , POINTER       :: Q
    TYPE(SRST   )              , POINTER       :: V
#else
    TYPE(BCSR)                 , INTENT(INout) :: D
    TYPE(BCSR)                 , INTENT(INOUT) :: Kx
#endif
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
    !-------------------------------------------------------------------
    REAL(DOUBLE), DIMENSION(MaxFuncPerAtmBlk**4              ) :: C
    TYPE(AtomPr), DIMENSION(MaxShelPerAtmBlk**2*CS_OUT%NCells) :: ACAtmPair,BDAtmPair
    !-------------------------------------------------------------------
    INTEGER                    :: LocNInt
    REAL(DOUBLE)               :: NInts
    !-------------------------------------------------------------------
    TYPE(ONX2OffSt) :: OffSet
    !-------------------------------------------------------------------
    REAL(DOUBLE), EXTERNAL :: DGetAbsMax
    !-------------------------------------------------------------------
    REAL(DOUBLE) :: TmBeg,TmEnd

    INTEGER :: ACR,BDR
    real(DOUBLE) :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
    integer :: i,iA,iB
    real(DOUBLE) :: time1,time2,SumInt
    integer :: MaxCont,IntType

    integer :: iint,iprint,isize
    !
    REAL(DOUBLE), PARAMETER :: ThresholdTwoE=-1.0D-15
    !REAL(DOUBLE), PARAMETER :: ThresholdTwoE=1.0D-12
    !-------------------------------------------------------------------
    !
    NULLIFY(AtAListTmp,AtAList,AtBListTmp,AtBList)
!!$    do iint=1,size(D%MTrix%D)
!!$    D%MTrix%D(iint)=dble(iint)
!!$    enddo
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
    !write(*,*) 'size C=',isize**4
    if(CS_OUT%NCells.GT.SIZE(ACAtmPair)) then
       write(*,*) 'size(ACAtmPair)',size(ACAtmPair),'.LT.',CS_OUT%NCells
       STOP 'Incrase the size of ACAtmPair and BDAtmPair'
    endif
    !write(*,*) 'CS_OUT%NCells=',CS_OUT%NCells
    tmp1=0.0d0
    tmp4=0.0d0
    !Simple check Simple check Simple check Simple check
    !
    SumInt=0.0d0
    MaxCont=0
    iint=0
    !
    LocNInt=0
    NInts=0.0d0
    !
#ifdef ONX2_PARALLEL
    P => D%Next                              ! Loop over atom C
    DO                               
       IF(.NOT.ASSOCIATED(P)) EXIT   
       AtC = P%Row
#else
    DO AtC=1,NAtoms ! Run over AtC.
#endif
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
#ifdef ONX2_PARALLEL
       U => P%RowRoot ! Run over AtD
       DO                                
          IF(.NOT.ASSOCIATED(U)) EXIT    
          IF(U%L.NE.U%R) THEN            
             U => U%Next                 
             CYCLE                       
          ENDIF
          AtD = U%L                      
          ! Set Time.                    
          TmBeg = MPI_WTIME()            
#else
       DO ci=D%RowPt%I(AtC),D%RowPt%I(AtC+1)-1 ! Run over AtD
          AtD = D%ColPt%I(ci)
          iPtrD= D%BlkPt%I(ci)
#endif
          KD=GM%AtTyp%I(AtD)
          NBFD=BS%BfKnd%I(KD)
          BDAtmInfo%Atm2X=GM%Carts%D(1,AtD)
          BDAtmInfo%Atm2Y=GM%Carts%D(2,AtD)
          BDAtmInfo%Atm2Z=GM%Carts%D(3,AtD)
          BDAtmInfo%K2=KD
          !
          ! Get max of the block density matrix.
#ifdef ONX2_PARALLEL
          Dcd=DGetAbsMax(NBFC*NBFD,U%MTrix(1,1))
#else
          Dcd=DGetAbsMax(NBFC*NBFD,D%MTrix%D(iPtrD))
#endif
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
             ACAtmInfo%NCell=GetNonNeglCell(AtAList,AtBListTmp%SqrtInt(1)*Dcd,ThresholdTwoE)
             IF(ACAtmInfo%NCell.EQ.0) THEN
                !write(*,'(A,E22.15,A,E22.15)') 'We skip 1', &
                !     !& AtAList%SqrtInt(1)*Dcd*AtBListTmp%SqrtInt(1), &
                !     & AtAList%SqrtInt(1)*AtBListTmp%SqrtInt(1), &
                !     & '.LT.',ThresholdTwoE  
                EXIT
             ENDIF
!!$             ENDIF
             !
             ! Find the row in Kx.
#ifdef ONX2_PARALLEL
             Q => FindFastMatRow_1(Kx,AtA)
#endif
             !
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
             call cpu_time(tmp2)
             CALL GetAtomPair2N(ACAtmInfo,AtAList,ACAtmPair,BS,CS_OUT)
             call cpu_time(tmp3)
             tmp1=tmp1+tmp3-tmp2
             !
             !
             AtBList=>AtBListTmp
             !
             DO ! Run over AtB
                AtB=AtBList%Atom 
                !
                !if(AtB.eq.1.and.AtA.eq.2) then
!                if(AtB.eq.1.and.AtD.eq.1.and.AtA.eq.2.and.AtC.eq.2) then
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
                   BDAtmInfo%NCell=GetNonNeglCell(AtBList,AtAList%SqrtInt(1)*Dcd,ThresholdTwoE)
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
                   ! Get atom pair for BD.
                   call cpu_time(tmp5)
                   CALL GetAtomPair2N(BDAtmInfo,AtBList,BDAtmPair,BS,CS_OUT)
                   call cpu_time(tmp6)
                   tmp4=tmp4+tmp6-tmp5
                   !
                   !CALL DBL_VECT_EQ_DBL_SCLR(NBFA*NBFB*NBFC*NBFD,C(1),BIG_DBL)
                   CALL DBL_VECT_EQ_DBL_SCLR(NBFA*NBFB*NBFC*NBFD,C(1),0.0d0)
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
                            ENDDO ! End blkfunc on D
                            OffSet%B=OffSet%B+BS%LStop%I(CFB,KB)-BS%LStrt%I(CFB,KB)+1
                         ENDDO ! End blkfunc on B
!!!!!!!!!!!!!!!!!!!!!!!!!!
                         ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!
                         OffSet%C=OffSet%C+BS%LStop%I(CFC,KC)-BS%LStrt%I(CFC,KC)+1
                      ENDDO ! End blkfunc on C
                      OffSet%A=OffSet%A+BS%LStop%I(CFA,KA)-BS%LStrt%I(CFA,KA)+1
                   ENDDO ! End blkfunc on A
!!!!!!!!!!!!!!!!!!!!
                   ENDDO
!!!!!!!!!!!!!!!!!!!!
!                   write(*,'(4I4)') AtA,AtB,AtC,AtD
#ifdef ONX2_DBUG
                   WRITE(*,'(A,E22.15,4I4)') ' MaxInt=',MAXVAL(C(1:NBFA*NBFB*NBFC*NBFD)),AtA,AtC,AtB,AtD
#endif
                   !
                   ! Get address for Kx and digest the block of integral.
#ifdef ONX2_PARALLEL
                   V => InsertSRSTNode(Q%RowRoot,AtB)
                   IF(.NOT.ASSOCIATED(V%MTrix)) THEN
                      ALLOCATE(V%MTrix(NBFA,NBFB),STAT=MemStatus)
                      !CALL IncMem(MemStatus,0,NBFA*NBFB,'AddFASTMATBlok')
                      CALL DBL_VECT_EQ_DBL_SCLR(NBFA*NBFB,V%MTrix(1,1),0.0d0)
                   ENDIF
                   !
                   CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,-1.0d0,C(1), &
                        &     NBFA*NBFB,U%MTrix(1,1),1,1.0d0, &
                        &     V%MTrix(1,1),1)
#else
                   CALL GetAdrB(AtA,AtB,Ind,Kx,0)
                   iPtrK = Kx%BlkPt%I(Ind)
                   !
                   CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,-1.0d0,C(1), &
                        &     NBFA*NBFB,D%MTrix%D(iPtrD),1,1.0d0, &
                        &     Kx%MTrix%D(IPtrK),1)
#endif
                   !
                   !CALL PrintMatrix(Kx%MTrix%D(IPtrK),NBFA,NBFB,2,TEXT_O='Int matrix')
                   !CALL Print_BCSR(Kx,'Kx',Unit_O=6)
                   !
                ENDIF
                !
!                endif
                !
                IF(.NOT.ASSOCIATED(AtBList%AtmNext)) EXIT
                AtBList=>AtBList%AtmNext
             ENDDO ! End AtB
             !
             IF(.NOT.ASSOCIATED(AtAList%AtmNext)) EXIT
             AtAList=>AtAList%AtmNext
          ENDDO ! End AtA
          !
#ifdef ONX2_PARALLEL
          ! Set Time.
          TmEnd = MPI_WTIME()
          !Add Time.
          U%Part = U%Part+TmEnd-TmBeg
          U => U%Next
#endif
       ENDDO ! End AtD
       !
#ifdef ONX2_PARALLEL
       P => P%Next
#endif
    ENDDO ! End AtC
    !
#ifdef ONX2_INFO
#ifdef ONX2_PARALLEL
    IF(MyID.EQ.ROOT) THEN
#endif
       WRITE(*,*) '-------------------------------------'
       WRITE(*,*) 'ComputK Statistic.'
       WRITE(*,'(A,F22.1)') ' Nbr ERI  =',NInts
       WRITE(*,'(A,I4)') ' Max Prim =',INT(SQRT(DBLE(MaxCont)))
       WRITE(*,*) '-------------------------------------'
#ifdef ONX2_PARALLEL
    ENDIF
#endif
#endif
  !
  END SUBROUTINE ComputK
  !
  !
  REAL(DOUBLE) FUNCTION GetAbsD(D,NC,ND,StrtC,StopC,StrtD,StopD)
    INTEGER :: NC,ND,StrtC,StopC,StrtD,StopD
    REAL(DOUBLE) :: D(NC,ND)
    INTEGER :: I,J
    GetAbsD=0.0d0
    DO I=StrtC,StopC
       DO J=StrtD,StopD
          GetAbsD=GetAbsD+ABS(D(I,J))
       ENDDO
    ENDDO
  END FUNCTION GetAbsD
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
#ifdef ONX2_DBUG
    AtmPair(:)%SP%IntType=BIG_INT
#endif
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
END MODULE ONX2ComputK

