MODULE HONXComputD2K
!H=================================================================================
!H MODULE HONXComputD2K
!H This MODULE contains:
!H  PUBLIC:
!H  o SUB 
!H
!H  PRIVATE:
!H  o SUB 
!H
!H  OPTIONS:
!H  DEBUGING: Use -DGONX2_DBUG to print some stuff.
!H  INFO    : Use -DGONX2_INFO to print some stuff.
!H
!H  Comments:
!H
!H---------------------------------------------------------------------------------
  !
!!$#ifndef PARALLEL
!!$#undef ONX2_PARALLEL
!!$#endif
  !
#ifdef GONX2_DBUG
#define GONX2_INFO
#endif
  !
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE Thresholding
  USE ONX2DataType
  USE InvExp
  USE ONXParameters
  !
!!$#ifdef ONX2_PARALLEL
!!$  USE MondoMPI
!!$  USE FastMatrices
!!$#endif
  !
  USE LinAlg
  USE ONX2List
  !
  IMPLICIT NONE
  !PRIVATE
  !
!--------------------------------------------------------------------------------- 
! PUBLIC DECLARATIONS
!--------------------------------------------------------------------------------- 
  PUBLIC  :: ComputDK2
  !
!--------------------------------------------------------------------------------- 
! PRIVATE DECLARATIONS
!--------------------------------------------------------------------------------- 
  PRIVATE :: GetAtomPairG
  !
CONTAINS
  !
  !
!!$#ifdef ONX2_PARALLEL
!!$  SUBROUTINE ComputDK(DFMcd,DFMab,GradX,ListC,ListD,GM,BS,CS_OUT)
!!$#else
  SUBROUTINE ComputDK(D,HessX,ListC,ListD,GM,BS,CS_OUT)
!!$#endif
!H---------------------------------------------------------------------------------
!H SUBROUTINE ComputDK(D,HessX,ListC,ListD,GM,BS,CS_OUT)
!H
!H---------------------------------------------------------------------------------
    !
    USE ONXGet, ONLY: GetAdrB
    !
    IMPLICIT NONE
    !
    !
    ! HessX = Dcd*(ac(R)|bd(R'))''*Dab
    !
    !-------------------------------------------------------------------
!!$#ifdef ONX2_PARALLEL
!!$    TYPE(FASTMAT)              , POINTER       :: DFMcd,DFMab
!!$    TYPE(FASTMAT)              , POINTER       :: P,Q
!!$    TYPE(SRST   )              , POINTER       :: U,V
!!$#else
    TYPE(BCSR)                 , INTENT(INout) :: D
!!$#endif
    TYPE(BCSR)                 , INTENT(INOUT) :: HessX
    TYPE(CList2) , DIMENSION(:), POINTER       :: ListC,ListD
    TYPE(CRDS)                 , INTENT(IN   ) :: GM
    TYPE(BSET)                 , INTENT(IN   ) :: BS
    TYPE(CellSet)              , INTENT(IN   ) :: CS_OUT
    !-------------------------------------------------------------------
    TYPE(ANode2), POINTER      :: AtAListTmp,AtAList,AtBListTmp,AtBList
    TYPE(AtomInfo)             :: ACAtmInfo,BDAtmInfo
    INTEGER                    :: AtA,AtB,AtC,AtD,KA,KB,KC,KD,CFA,CFB,CFC,CFD
    INTEGER                    :: ci,iPtrD,iPtrD2,iPtrK,NBFC,NBFD,NBFA,NBFB
    INTEGER                    :: CFAC,CFBD
    INTEGER                    :: NCFncA,NCFncB,NCFncC,NCFncD
    INTEGER                    :: Off,Ind,LocNInt,IntType
    INTEGER                    :: ACR,BDR,IXYZ,NIntBlk,Indx
!!$#ifdef ONX2_PARALLEL
!!$    REAL(DOUBLE)               :: TmBeg,TmEnd
!!$#endif
    REAL(DOUBLE)               :: TmpGradA,TmpGradC,TmpGradB,TmpGradD
    REAL(DOUBLE)               :: Dcd,Dab,NInts
    !-------------------------------------------------------------------
    REAL(DOUBLE) , DIMENSION(12*MaxFuncPerAtmBlk**4) :: C  ! MUST CHANGE THE ARRAY SIZE.
    REAL(DOUBLE) , DIMENSION(   MaxFuncPerAtmBlk**2) :: Work
    TYPE(AtomPrG), DIMENSION(   MaxShelPerAtmBlk**2*CS_OUT%NCells) :: ACAtmPair,BDAtmPair
    REAL(DOUBLE), DIMENSION(MaxShelPerAtmBlk**2) :: DMcd,DMab
    !REAL(DOUBLE), DIMENSION(100) :: DMcd,DMab
    !-------------------------------------------------------------------
    TYPE(ONX2OffSt)            :: OffSet
    TYPE(INT_VECT )            :: BColIdx
    !-------------------------------------------------------------------
    REAL(DOUBLE), EXTERNAL     :: DGetAbsMax
    REAL(DOUBLE), EXTERNAL     :: DDOT
    !-------------------------------------------------------------------
    !
    integer :: isize,i
    !
    !
    ! Initialize.
    NULLIFY(AtAListTmp,AtAList,AtBListTmp,AtBList)
!!$#ifdef ONX2_PARALLEL
!!$    NULLIFY(P,Q,U,V)
!!$#endif
    !
    !Simple check Simple check Simple check Simple check
    isize=0
    do i=1,natoms
       isize=MAX(isize,BS%BfKnd%I(GM%AtTyp%I(i)))
       IF(12*(BS%BfKnd%I(GM%AtTyp%I(i)))**4.GT.SIZE(C)) THEN
          write(*,*) 'size',12*BS%BfKnd%I(GM%AtTyp%I(i))**4
          STOP 'Incrase the size of C'
       ENDIF
    enddo
    !write(*,*) 'size C=',12*isize**4
    if(CS_OUT%NCells.GT.SIZE(ACAtmPair)) then
       write(*,*) 'size(ACAtmPair)',size(ACAtmPair),'.LT.',CS_OUT%NCells
       STOP 'Incrase the size of ACAtmPair and BDAtmPair'
    endif
    !write(*,*) 'CS_OUT%NCells=',CS_OUT%NCells
    !Simple check Simple check Simple check Simple check
    !call Print_BCSR(D,'Density matrix',Unit_O=6)
    !
    CALL New(BColIdx,Natoms)
    !
    LocNInt=0
    NInts=0.0D0
    !
!!$#ifdef ONX2_PARALLEL
!!$    P => DFMcd%Next ! Run over AtC.
!!$    DO                               
!!$       IF(.NOT.ASSOCIATED(P)) EXIT   
!!$       AtC = P%Row                   
!!$       !write(*,*) 'AtC=',AtC,'MyID',MyID
!!$#else
    DO AtC=1,NAtoms ! Run over AtC.
!!$#endif
       !
       KC=GM%AtTyp%I(AtC)
       NBFC=BS%BfKnd%I(KC)
       NCFncC=BS%NCFnc%I(KC)
       !
       ACAtmInfo%Atm2X=GM%Carts%D(1,AtC)
       ACAtmInfo%Atm2Y=GM%Carts%D(2,AtC)
       ACAtmInfo%Atm2Z=GM%Carts%D(3,AtC)
       ACAtmInfo%K2=KC
       !
       ! Get AtA List.
       AtAListTmp=>ListC(AtC)%GoList
       !
!!$#ifdef ONX2_PARALLEL
!!$       U => P%RowRoot ! Run over AtD
!!$       DO
!!$          IF(.NOT.ASSOCIATED(U)) EXIT
!!$          IF(U%L.NE.U%R) THEN
!!$             U => U%Next
!!$             CYCLE
!!$          ENDIF
!!$          AtD = U%L
!!$          !write(*,*) 'AtC=',AtC,'AtD=',AtD,'MyID',MyID
!!$          ! Set Time.
!!$          TmBeg = MPI_WTIME()
!!$#else
       DO ci=D%RowPt%I(AtC),D%RowPt%I(AtC+1)-1 ! Run over AtD
          AtD=D%ColPt%I(ci)
          iPtrD=D%BlkPt%I(ci)
!!$#endif
          !
          KD=GM%AtTyp%I(AtD)
          NBFD=BS%BfKnd%I(KD)
          NCFncD=BS%NCFnc%I(KD)
          !
          BDAtmInfo%Atm2X=GM%Carts%D(1,AtD)
          BDAtmInfo%Atm2Y=GM%Carts%D(2,AtD)
          BDAtmInfo%Atm2Z=GM%Carts%D(3,AtD)
          BDAtmInfo%K2=KD
          !
          ! Get max of the block density matrix.
!!$#ifdef ONX2_PARALLEL
!!$          Dcd=DGetAbsMax(NBFC*NBFD,U%MTrix(1,1))
!!$          !write(*,*) 'Dcd',Dcd,'AtC=',AtC,'AtD=',AtD,'MyID',MyID
!!$#else
          Dcd=DGetAbsMax(NBFC*NBFD,D%MTrix%D(iPtrD))
!!$#endif
          !
!!$#ifdef ONX2_PARALLEL
!!$          CALL GetAbsDenBlk(U%MTrix(1,1),NBFC,NBFD,DMcd(1),    &
!!$               &            BS%NCFnc%I(KC),BS%NCFnc%I(KD),     &
!!$               &            BS%LStrt%I(1,KC),BS%LStop%I(1,KC), &
!!$               &            BS%LStrt%I(1,KD),BS%LStop%I(1,KD)  )
!!$#else
          CALL GetAbsDenBlk(D%MTrix%D(iPtrD),NBFC,NBFD,DMcd(1),&
               &            BS%NCFnc%I(KC),BS%NCFnc%I(KD),     &
               &            BS%LStrt%I(1,KC),BS%LStop%I(1,KC), &
               &            BS%LStrt%I(1,KD),BS%LStop%I(1,KD)  )
!!$#endif
          !
#ifdef GONX2_DBUG
          WRITE(*,*) 'Max(Dcd)=',Dcd
#endif
          !
          ! Get AtB List.
          AtBListTmp=>ListD(AtD)%GoList
          !
          AtAList=>AtAListTmp
          !
          RnOvA: DO ! Run over AtA
             AtA=AtAList%Atom
!!$#ifdef ONX2_PARALLEL
!!$             !if(myid==0)write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
!!$             !if(myid==0)write(*,*)'AtA in list',AtA,ASSOCIATED(Q)
!!$             !if(myid==1)write(*,*)'Q%Row<AtA',Q%Row,'<',AtA
!!$             ! May change that to FindFASTMATRow
!!$             Q => DFMab%Next
!!$             DO
!!$                IF(.NOT.ASSOCIATED(Q)) EXIT RnOvA
!!$                IF(Q%Row.EQ.AtA) EXIT
!!$                IF(Q%Row.GT.AtA) THEN
!!$                   IF(.NOT.ASSOCIATED(AtAList%AtmNext)) EXIT RnOvA
!!$                   AtAList=>AtAList%AtmNext
!!$                   CYCLE RnOvA
!!$                ENDIF
!!$                Q => Q%Next
!!$             ENDDO
!!$             !if(myid==0)write(*,*) 'We find AtC=',AtC,'AtD=',AtD,'AtA',AtA,'MyID',MyID
!!$#endif
             !
             KA=GM%AtTyp%I(AtA)
             NBFA=BS%BfKnd%I(KA)
             NCFncA=BS%NCFnc%I(KA)
             !
             ACAtmInfo%NCell=GetNonNeglCell(AtAList,AtBListTmp%SqrtInt(1)*Dcd,Thresholds%TwoE &
#ifdef GTRESH
               & )
#else
               & *(-1d0))
#endif
             !write(*,*) 'ACAtmInfo%NCell',ACAtmInfo%NCell
             IF(ACAtmInfo%NCell.EQ.0) EXIT RnOvA
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
             CALL GetAtomPairG(ACAtmInfo,AtAList,ACAtmPair,BS,CS_OUT)
             !
             AtBList=>AtBListTmp
             !
#ifndef ONX2_PARALLEL
             CALL GetColIdx(AtA,D,BColIdx)
#endif
             !
             !#ifdef ONX2_PARALLEL   !I may not need
             !V => Q%RowRoot         !I may not need
             !#endif                 !I may not need
             !
             RnOvB: DO ! Run over AtB
                AtB=AtBList%Atom 
                KB=GM%AtTyp%I(AtB)
                NBFB=BS%BfKnd%I(KB)
                NCFncB=BS%NCFnc%I(KB)
                !
!!$#ifdef ONX2_PARALLEL
!!$                !if(myid==1)write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++'
!!$                ! May change that to FindFASTMATCol
!!$                V => Q%RowRoot
!!$                DO
!!$                   IF(.NOT.ASSOCIATED(V)) EXIT RnOvB !use the binary tree would be better.
!!$                   IF(V%L.NE.V%R) THEN
!!$                      V => V%Next
!!$                      CYCLE
!!$                   ENDIF
!!$                   !
!!$                   IF(V%L.EQ.AtB) EXIT
!!$                   IF(V%L.GT.AtB) THEN
!!$                      IF(.NOT.ASSOCIATED(AtBList%AtmNext)) EXIT RnOvB
!!$                      AtBList=>AtBList%AtmNext
!!$                      CYCLE RnOvB
!!$                   ENDIF
!!$                   V => V%Next
!!$                ENDDO
!!$                !if(myid==1)
!!$                !write(*,*) 'We find AtC=',AtC,'AtD=',AtD,'AtA',AtA,'AtB',AtB,'MyID',MyID
!!$#else
                Ind=BColIdx%I(AtB)
                IF(Ind.GT.0) THEN ! Skip out if no density matrix element.
                   iPtrD2 = D%BlkPt%I(Ind)
!!$#endif
                   !
                   ! Get max of the block density matrix.
!!$#ifdef ONX2_PARALLEL
!!$                   Dab=DGetAbsMax(NBFA*NBFB,V%MTrix(1,1))
!!$#else
                   Dab=DGetAbsMax(NBFA*NBFB,D%MTrix%D(iPtrD2))
!!$#endif
                   !
                   BDAtmInfo%NCell=GetNonNeglCell(AtBList,AtAList%SqrtInt(1)*Dab*Dcd*Half,Thresholds%TwoE &
#ifdef GTRESH
                   & )
#else
                   & *(-1d0))
#endif
                   IF(BDAtmInfo%NCell.EQ.0) EXIT RnOvB
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
!!$#ifdef ONX2_PARALLEL
!!$                   CALL GetAbsDenBlk(V%MTrix(1,1),NBFA,NBFB,DMab(1),     &
!!$                        &            BS%NCFnc%I(KA),BS%NCFnc%I(KB),      &
!!$                        &            BS%LStrt%I(1,KA),BS%LStop%I(1,KA),  &
!!$                        &            BS%LStrt%I(1,KB),BS%LStop%I(1,KB)   )
!!$#else
                   CALL GetAbsDenBlk(D%MTrix%D(iPtrD2),NBFA,NBFB,DMab(1),&
                        &            BS%NCFnc%I(KA),BS%NCFnc%I(KB),      &
                        &            BS%LStrt%I(1,KA),BS%LStop%I(1,KA),  &
                        &            BS%LStrt%I(1,KB),BS%LStop%I(1,KB)   )
!!$#endif
                   !
                   ! Get atom pair for BD.
                   CALL GetAtomPairG(BDAtmInfo,AtBList,BDAtmPair,BS,CS_OUT)
                   !
                   NIntBlk=NBFA*NBFB*NBFC*NBFD
                   !
                   CALL DBL_VECT_EQ_DBL_SCLR(12*NIntBlk,C(1),0.0D0)
                   !
                   CFAC=0
!!!!!!!!!!!!!!!!!!!!
                   DO ACR=1,ACAtmInfo%NCell
!!!!!!!!!!!!!!!!!!!!
                   OffSet%A=1
                   DO CFA=1,NCFncA ! Run over blkfunc on A
                      OffSet%C=1
                      DO CFC=1,NCFncC ! Run over blkfunc on C
                         CFAC=CFAC+1
                         CFBD=0                         
!!!!!!!!!!!!!!!!!!!!!!!!!!
                         DO BDR=1,BDAtmInfo%NCell
!!!!!!!!!!!!!!!!!!!!!!!!!!
                         OffSet%B=1
                         DO CFB=1,NCFncB ! Run over blkfunc on B
                            OffSet%D=1
                            DO CFD=1,NCFncD ! Run over blkfunc on D
                               CFBD=CFBD+1
                               !
                               !if(DMcd((CFC-1)*NCFncD+CFD)>Dcd) stop
                               !if(DMab((CFA-1)*NCFncB+CFB)>Dab) stop
                               !
#ifdef GTRESH
                               IF(DMcd((CFC-1)*NCFncD+CFD)*DMab((CFA-1)*NCFncB+CFB)*Half* &
                                    & AtAList%SqrtInt(1)*AtBList%SqrtInt(1)>Thresholds%TwoE) THEN
#endif
                                  !
                                  ! Compute integral type.
                                  IntType=ACAtmPair(CFAC)%SP%IntType*10000+BDAtmPair(CFBD)%SP%IntType
                                  !
                                  ! The integral interface.
                                  INCLUDE 'D2ERIInterface.Inc'
                                  !
                                  NInts=NInts+DBLE(LocNInt)
#ifdef GTRESH
                               ENDIF
#endif
                               !
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
!!!!!!!!!!!!!!!!!!!!!!!!!!
                   ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!
                   !
                   ! Compute the exchange-hessian.
                   DO IXYZ=1,3
                      !
                      !-----------------------------------------------------
                      ! AtA
                      Indx=(IXYZ-1)*NIntBlk+1
                      !
!!$#ifdef ONX2_PARALLEL
!!$                      CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,C(Indx), &
!!$                           &     NBFA*NBFB,U%MTrix(1,1),1,0.0d0,    &
!!$                           &     Work(1),1)
!!$                      TmpGradA=-DDOT(NBFA*NBFB,V%MTrix(1,1),1,Work(1),1)
!!$#else
                      CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,C(Indx), &
                           &     NBFA*NBFB,D%MTrix%D(iPtrD),1,0.0d0,    &
                           &     Work(1),1)
                      TmpGradA=-DDOT(NBFA*NBFB,D%MTrix%D(iPtrD2),1,Work(1),1)
!!$#endif
                      GradX%D(IXYZ,AtA)=GradX%D(IXYZ,AtA)+TmpGradA
                      !
                      !-----------------------------------------------------
                      ! AtC
                      Indx=(3+IXYZ-1)*NIntBlk+1
!!$#ifdef ONX2_PARALLEL
!!$                      CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,C(Indx), &
!!$                           &     NBFA*NBFB,U%MTrix(1,1),1,0.0d0,    &
!!$                           &     Work(1),1)
!!$                      TmpGradC=-DDOT(NBFA*NBFB,V%MTrix(1,1),1,Work(1),1)
!!$#else
                      CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,C(Indx), &
                           &     NBFA*NBFB,D%MTrix%D(iPtrD),1,0.0d0,    &
                           &     Work(1),1)
                      TmpGradC=-DDOT(NBFA*NBFB,D%MTrix%D(iPtrD2),1,Work(1),1)
!!$#endif
                      GradX%D(IXYZ,AtC)=GradX%D(IXYZ,AtC)+TmpGradC
                      !
                      !-----------------------------------------------------
                      ! AtB
                      Indx=(6+IXYZ-1)*NIntBlk+1
!!$#ifdef ONX2_PARALLEL
!!$                      CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,C(Indx), &
!!$                           &     NBFA*NBFB,U%MTrix(1,1),1,0.0d0,    &
!!$                           &     Work(1),1)
!!$                      TmpGradB=-DDOT(NBFA*NBFB,V%MTrix(1,1),1,Work(1),1)
!!$#else
                      CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,C(Indx), &
                           &     NBFA*NBFB,D%MTrix%D(iPtrD),1,0.0d0,    &
                           &     Work(1),1)
                      TmpGradB=-DDOT(NBFA*NBFB,D%MTrix%D(iPtrD2),1,Work(1),1)
!!$#endif
                      GradX%D(IXYZ,AtB)=GradX%D(IXYZ,AtB)+TmpGradB
                      !
                      !-----------------------------------------------------
                      ! AtD
                      GradX%D(IXYZ,AtD)=GradX%D(IXYZ,AtD)-(TmpGradA+TmpGradC+TmpGradB)
                      !
                   ENDDO
                   !
#ifndef ONX2_PARALLEL
                ENDIF ! Skip out if no density matrix element.
#endif
                !
                IF(.NOT.ASSOCIATED(AtBList%AtmNext)) EXIT RnOvB
                AtBList=>AtBList%AtmNext
                !
             ENDDO RnOvB ! End AtB
             !
             IF(.NOT.ASSOCIATED(AtAList%AtmNext)) EXIT RnOvA
             AtAList=>AtAList%AtmNext
             !
          ENDDO RnOvA ! End AtA
          !
!!$#ifdef ONX2_PARALLEL
!!$          ! Set Time.
!!$          TmEnd = MPI_WTIME()
!!$          !Add Time.
!!$          U%Part = U%Part+TmEnd-TmBeg
!!$          U => U%Next
!!$#endif
          !
       ENDDO ! End AtD
       !
!!$#ifdef ONX2_PARALLEL
!!$       P => P%Next
!!$#endif
       !
    ENDDO ! End AtC
    !
    CALL Delete(BColIdx)
    !
    WRITE(*,100) NInts,12D0*CS_OUT%NCells**2*DBLE(MaxNon0-1)**2, &
         &       NInts/(12D0*CS_OUT%NCells**2*DBLE(MaxNon0-1)**2)*100D0
100 FORMAT(' NInts = ',E8.2,' NIntTot = ',E8.2,' Ratio = ',E8.2,'%')
    !
  END SUBROUTINE ComputDK
  !
  !
  SUBROUTINE GetAtomPairG(AtmInfo,List,AtmPair,BS,CS_OUT)
!H---------------------------------------------------------------------------------
!H SUBROUTINE GetAtomPairG(AtmInfo,List,AtmPair,BS,CS_OUT)
!H
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(AtomInfo)              , INTENT(IN   ) :: AtmInfo
    TYPE(ANode2  ), POINTER                     :: List
    TYPE(AtomPrG ), DIMENSION(:), INTENT(INOUT) :: AtmPair
    TYPE(BSET)                  , INTENT(IN   ) :: BS
    TYPE(CellSet)               , INTENT(IN   ) :: CS_OUT
    !-------------------------------------------------------------------
    INTEGER      :: CF12,CF1,CF2,I1,I2,II,JJ,IJ,iCell,Cell
    INTEGER      :: MinL1,MaxL1,Type1,MinL2,MaxL2,Type2
    INTEGER      :: StopL1,StartL1,StopL2,StartL2
    INTEGER      :: ISwitch1,ISwitch2
    REAL(DOUBLE) :: Z1,Z2,Expt,InvExpt,R12,XiR12,RX,RY,RZ
    LOGICAL      :: Switch
    !-------------------------------------------------------------------
    !
    ! dbuging
#ifdef GONX2_DBUG
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
                      AtmPair(CF12)%SP%Cst(5,IJ)=5.914967172796D0*EXP(-XiR12)*InvExpt* &
                      !AtmPair(CF12)%SP%Cst(5,IJ)=5.914967172796D0*EXPInv(XiR12)*InvExpt* &
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
                ! We can skipp out the loop if we did not get any significant primitives.
                IF(JJ.EQ.0) EXIT
                II=II+JJ
             ENDDO
             !
             AtmPair(CF12)%SP%L=II
             !
             ! We reorder the atomic positions if Type2 > Type1.
             ! Needed for the integral evaluations.
             IF(Switch) THEN
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
#ifdef GONX2_DBUG
             WRITE(*,'(3(A,I3),A,I5)') 'Type1',Type1,' Type2',Type2, &
                  &                    ' IntType',AtmPair(CF12)%SP%IntType
#endif
             !
          ENDDO
       ENDDO
    ENDDO ! NCell
    !
  END SUBROUTINE GetAtomPairG
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
  !
  !
END MODULE GONX2ComputDK











