MODULE GONX2ComputDK
!H=================================================================================
!H MODULE GONX2ComputDK
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
#ifndef PARALLEL
#undef ONX2_PARALLEL
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
#ifdef ONX2_PARALLEL
  USE MondoMPI
  USE FastMatrices
#endif
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
  PUBLIC  :: ComputDK
  !
!--------------------------------------------------------------------------------- 
! PRIVATE DECLARATIONS
!--------------------------------------------------------------------------------- 
  PRIVATE :: GetAtomPairG
  !
CONTAINS
  !
  !
#ifdef ONX2_PARALLEL
  SUBROUTINE ComputDK(DFMcd,DFMab,GradX,BoxX,ListC,ListD,OffArr,GMc,BSc,CS_OUT)
#else
  SUBROUTINE ComputDK(D,GradX,BoxX,ListC,ListD,OffArr,GMc,BSc,CS_OUT)
#endif
!H---------------------------------------------------------------------------------
!H SUBROUTINE ComputDK(D,GradX,ListC,ListD,GMc,BSc,CS_OUT)
!H
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !
    ! GradX = Dcd*(ac(R)|bd(R'))'*Dab
    !
    !-------------------------------------------------------------------
#ifdef ONX2_PARALLEL
    TYPE(FASTMAT)              , POINTER :: DFMcd,DFMab
    TYPE(FASTMAT)              , POINTER :: P,Q
    TYPE(SRST   )              , POINTER :: U,V
#else
    TYPE(BCSR)                           :: D
#endif
    TYPE(DBL_RNK2)                       :: GradX,BoxX
    TYPE(INT_RNK2)                       :: OffArr
    TYPE(CList) , DIMENSION(:) , POINTER :: ListC,ListD
    TYPE(CRDS)                           :: GMc
    TYPE(BSET)                           :: BSc
    TYPE(CellSet)                        :: CS_OUT
    !-------------------------------------------------------------------
    TYPE(ANode), POINTER       :: AtAListTmp,AtAList,AtBListTmp,AtBList
    TYPE(AtomInfo)             :: ACAtmInfo,BDAtmInfo
    INTEGER                    :: AtA,AtB,AtC,AtD,KA,KB,KC,KD,CFA,CFB,CFC,CFD
    INTEGER                    :: ci,iPtrD,iPtrD2,NBFC,NBFD,NBFA,NBFB
    INTEGER                    :: iErr
    INTEGER                    :: NCFncA,NCFncB,NCFncC,NCFncD
    INTEGER                    :: Ind,LocNInt,IntType,Split
    INTEGER                    :: IXYZ,JXYZ,NIntBlk,Indx,iFAC,iFBD
    INTEGER                    :: OT,OAL,OBL,LDAL,LDBL,GOAL,GOBL
    INTEGER                    :: OA,OB,OC,OD,LDA,LDB,LDC,LDD,GOA,GOB,GOC,GOD
#ifdef ONX2_PARALLEL
    REAL(DOUBLE)               :: TmBeg,TmEnd
#endif
    REAL(DOUBLE)               :: TmpGradA,TmpGradC,TmpGradB
    REAL(DOUBLE)               :: Dcd,Dab,NInts,NIntsTot
    !-------------------------------------------------------------------
    REAL(DOUBLE) , DIMENSION(MaxFuncPerAtmBlk**2) :: Work
    REAL(DOUBLE) , DIMENSION(MaxShelPerAtmBlk**2) :: DMcd,DMab
    REAL(DOUBLE) , DIMENSION(12*MaxFuncPerAtmBlk**4) :: C  !,C_
    !STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS
    REAL(DOUBLE) , DIMENSION( 9*MaxFuncPerAtmBlk**4) :: CC
    !STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS
    TYPE(AtomPrG), DIMENSION(:), ALLOCATABLE :: ACAtmPair,BDAtmPair !size=MaxShelPerAtmBlk**2*CS_OUT%NCells
    !-------------------------------------------------------------------
    TYPE(ONX2OffSt)            :: OffSet
    TYPE(INT_VECT )            :: BColIdx
    !-------------------------------------------------------------------
    REAL(DOUBLE), EXTERNAL     :: DGetAbsMax
    REAL(DOUBLE), EXTERNAL     :: DDOT
    REAL(DOUBLE), EXTERNAL     :: MondoTimer
    !-------------------------------------------------------------------
    !
    integer :: isize,i
    logical :: DoStress
    DoStress=.FALSE.
    !
    ! Initialize.
    NULLIFY(AtAListTmp,AtAList,AtBListTmp,AtBList)
#ifdef ONX2_PARALLEL
    NULLIFY(P,Q,U,V)
#endif
    !
    ! Allocate arrays.
    ALLOCATE(ACAtmPair(MaxShelPerAtmBlk**2*CS_OUT%NCells), &
         &   BDAtmPair(MaxShelPerAtmBlk**2*CS_OUT%NCells),STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In ComputDK: Allocation problem.')
    !
    !
    !Simple check Simple check Simple check Simple check
    isize=0
    do i=1,natoms
       isize=MAX(isize,BSc%BfKnd%I(GMc%AtTyp%I(i)))
       IF(12*(BSc%BfKnd%I(GMc%AtTyp%I(i)))**4.GT.SIZE(C)) THEN
          write(*,*) 'size',12*BSc%BfKnd%I(GMc%AtTyp%I(i))**4
          STOP 'Incrase the size of C'
       ENDIF
    enddo
    !write(*,*) 'size C=',12*isize**4
    if(CS_OUT%NCells.GT.SIZE(ACAtmPair)) then
       write(*,*) 'size(ACAtmPair)',size(ACAtmPair),'.LT.',CS_OUT%NCells
       STOP 'Incrase the size of ACAtmPair and BDAtmPair'
    endif
    !write(*,*) 'CS_OUT%NCells=',CS_OUT%NCells
    !write(*,*) 'CS_OUT%CellCarts=',CS_OUT%CellCarts%D(2,1:CS_OUT%NCells)
    !Simple check Simple check Simple check Simple check
    !
    CALL New(BColIdx,NAtoms)
    !
    LocNInt=0
    NInts=0.0D0
    NIntsTot=0.0D0
    !
#ifdef ONX2_PARALLEL
    P => DFMcd%Next ! Run over AtC.
    DO                               
       IF(.NOT.ASSOCIATED(P)) EXIT   
       AtC = P%Row                   
       !write(*,*) 'AtC=',AtC,'MyID',MyID
#else
    DO AtC=1,NAtoms ! Run over AtC.
#endif
       !
       KC=GMc%AtTyp%I(AtC)
       NBFC=BSc%BfKnd%I(KC)
       NCFncC=BSc%NCFnc%I(KC)
       !
       ACAtmInfo%Atm2X=GMc%Carts%D(1,AtC)
       ACAtmInfo%Atm2Y=GMc%Carts%D(2,AtC)
       ACAtmInfo%Atm2Z=GMc%Carts%D(3,AtC)
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
          !write(*,*) 'AtC=',AtC,'AtD=',AtD,'MyID',MyID
          ! Set Time.
          TmBeg = MondoTimer()
#else
       DO ci=D%RowPt%I(AtC),D%RowPt%I(AtC+1)-1 ! Run over AtD
          AtD=D%ColPt%I(ci)
          iPtrD=D%BlkPt%I(ci)
#endif
          !
          KD=GMc%AtTyp%I(AtD)
          NBFD=BSc%BfKnd%I(KD)
          NCFncD=BSc%NCFnc%I(KD)
          !
          BDAtmInfo%Atm2X=GMc%Carts%D(1,AtD)
          BDAtmInfo%Atm2Y=GMc%Carts%D(2,AtD)
          BDAtmInfo%Atm2Z=GMc%Carts%D(3,AtD)
          BDAtmInfo%K2=KD
          !
          ! Get max of the block density matrix.
#ifdef ONX2_PARALLEL
          Dcd=DGetAbsMax(NBFC*NBFD,U%MTrix(1,1))
          !write(*,*) 'Dcd',Dcd,'AtC=',AtC,'AtD=',AtD,'MyID',MyID
#else
          Dcd=DGetAbsMax(NBFC*NBFD,D%MTrix%D(iPtrD))
#endif
          !
#ifdef ONX2_PARALLEL
          CALL GetAbsDenBlk(U%MTrix(1,1),NBFC,NBFD,DMcd(1),    &
               &            BSc%NCFnc%I(KC),BSc%NCFnc%I(KD),     &
               &            BSc%LStrt%I(1,KC),BSc%LStop%I(1,KC), &
               &            BSc%LStrt%I(1,KD),BSc%LStop%I(1,KD)  )
#else
          CALL GetAbsDenBlk(D%MTrix%D(iPtrD),NBFC,NBFD,DMcd(1),&
               &            BSc%NCFnc%I(KC),BSc%NCFnc%I(KD),     &
               &            BSc%LStrt%I(1,KC),BSc%LStop%I(1,KC), &
               &            BSc%LStrt%I(1,KD),BSc%LStop%I(1,KD)  )
#endif
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
#ifdef ONX2_PARALLEL
             !
             Q => DFMab%Next
             DO
                IF(.NOT.ASSOCIATED(Q)) EXIT RnOvA
                IF(Q%Row.EQ.AtA) EXIT
                IF(Q%Row.GT.AtA) THEN
                   IF(.NOT.ASSOCIATED(AtAList%AtmNext)) EXIT RnOvA
                   AtAList=>AtAList%AtmNext
                   CYCLE RnOvA
                ENDIF
                Q => Q%Next
             ENDDO
             !if(myid==0)write(*,*) 'We find AtC=',AtC,'AtD=',AtD,'AtA',AtA,'MyID',MyID
#endif
             !
             KA=GMc%AtTyp%I(AtA)
             NBFA=BSc%BfKnd%I(KA)
             NCFncA=BSc%NCFnc%I(KA)
             !
             ACAtmInfo%NFPair=GetNonNFPair(AtAList,AtBListTmp%RInt(1)*Dcd,Thresholds%TwoE &
#ifdef GTRESH
             & )
#else
             & *(-1d0))
#endif
             IF(ACAtmInfo%NFPair.EQ.0) EXIT RnOvA
             !
             ACAtmInfo%Atm1X=GMc%Carts%D(1,AtA)
             ACAtmInfo%Atm1Y=GMc%Carts%D(2,AtA)
             ACAtmInfo%Atm1Z=GMc%Carts%D(3,AtA)
             ACAtmInfo%K1=KA
             !
             ACAtmInfo%Atm12X=ACAtmInfo%Atm1X-ACAtmInfo%Atm2X
             ACAtmInfo%Atm12Y=ACAtmInfo%Atm1Y-ACAtmInfo%Atm2Y
             ACAtmInfo%Atm12Z=ACAtmInfo%Atm1Z-ACAtmInfo%Atm2Z
             !
             ! Get atom pair for BD.
             CALL GetAtomPairG(ACAtmInfo,AtAList,ACAtmPair,BSc,CS_OUT)
             !
             AtBList=>AtBListTmp
             !
#ifndef ONX2_PARALLEL
             CALL GetColIdx(AtA,D,BColIdx)
#endif
             !
             RnOvB: DO ! Run over AtB
                AtB=AtBList%Atom 
                KB=GMc%AtTyp%I(AtB)
                NBFB=BSc%BfKnd%I(KB)
                NCFncB=BSc%NCFnc%I(KB)
                !
#ifdef ONX2_PARALLEL
                !---------------------------------------------------------------------
                ! V => Q%RowRoot
                ! DO
                !    IF(.NOT.ASSOCIATED(V)) EXIT RnOvB
                !    IF(V%L.NE.V%R) THEN
                !       V => V%Next
                !       CYCLE
                !    ENDIF
                !    !
                !    IF(V%R.EQ.AtB) EXIT
                !    IF(V%L.GT.AtB) THEN
                !       IF(.NOT.ASSOCIATED(AtBList%AtmNext)) EXIT RnOvB
                !       AtBList=>AtBList%AtmNext
                !       CYCLE RnOvB
                !    ENDIF
                !    V => V%Next
                ! ENDDO
                !---------------------------------------------------------------------
                ! V => Q%RowRoot
                ! DO
                !    IF(.NOT.ASSOCIATED(V)) THEN
                !       IF(.NOT.ASSOCIATED(AtBList%AtmNext)) EXIT RnOvB
                !       AtBList=>AtBList%AtmNext
                !       CYCLE RnOvB
                !    ENDIF
                !    IF(V%L.NE.V%R) THEN
                !       V => V%Next
                !       CYCLE
                !    ENDIF
                !    !
                !    IF(V%R.EQ.AtB) EXIT
                !    IF(V%L.GT.AtB) THEN
                !       IF(.NOT.ASSOCIATED(AtBList%AtmNext)) EXIT RnOvB
                !       AtBList=>AtBList%AtmNext
                !       CYCLE RnOvB
                !    ENDIF
                !    V => V%Next
                ! ENDDO
                !---------------------------------------------------------------------
                V=>Q%RowRoot
                DO
                   IF(.NOT.ASSOCIATED(V)) THEN
                      IF(.NOT.ASSOCIATED(AtBList%AtmNext)) EXIT RnOvB
                      AtBList=>AtBList%AtmNext
                      CYCLE RnOvB
                   ENDIF
                   IF(AtB.EQ.V%L.AND.AtB.EQ.V%R) EXIT
                   Split=IntervalSplit(V%L,V%R)
                   IF(AtB.GE.V%L.AND.AtB.LE.Split) THEN
                      V=>V%Left
                   ELSE
                      V=>V%Right
                   ENDIF
                ENDDO
                !write(*,*) 'We find AtC=',AtC,'AtD=',AtD,'AtA',AtA,'AtB',AtB,'MyID',MyID
#else
                Ind=BColIdx%I(AtB)
                IF(Ind.GT.0) THEN ! Skip out if no density matrix element.
                   iPtrD2 = D%BlkPt%I(Ind)
#endif
                   !
                   ! Get max of the block density matrix.
#ifdef ONX2_PARALLEL
                   Dab=DGetAbsMax(NBFA*NBFB,V%MTrix(1,1))
#else
                   Dab=DGetAbsMax(NBFA*NBFB,D%MTrix%D(iPtrD2))
#endif
                   !
                   BDAtmInfo%NFPair=GetNonNFPair(AtBList,AtAList%RInt(1)*Dab*Dcd*Half,Thresholds%TwoE &
#ifdef GTRESH
                   & )
#else
                   & *(-1d0))
#endif
                   IF(BDAtmInfo%NFPair.EQ.0) EXIT RnOvB
                   !
                   BDAtmInfo%Atm1X=GMc%Carts%D(1,AtB)
                   BDAtmInfo%Atm1Y=GMc%Carts%D(2,AtB)
                   BDAtmInfo%Atm1Z=GMc%Carts%D(3,AtB)
                   BDAtmInfo%K1=KB
                   !
                   BDAtmInfo%Atm12X=BDAtmInfo%Atm1X-BDAtmInfo%Atm2X
                   BDAtmInfo%Atm12Y=BDAtmInfo%Atm1Y-BDAtmInfo%Atm2Y
                   BDAtmInfo%Atm12Z=BDAtmInfo%Atm1Z-BDAtmInfo%Atm2Z
                   !
#ifdef ONX2_PARALLEL
                   CALL GetAbsDenBlk(V%MTrix(1,1),NBFA,NBFB,DMab(1),     &
                        &            BSc%NCFnc%I(KA),BSc%NCFnc%I(KB),      &
                        &            BSc%LStrt%I(1,KA),BSc%LStop%I(1,KA),  &
                        &            BSc%LStrt%I(1,KB),BSc%LStop%I(1,KB)   )
#else
                   CALL GetAbsDenBlk(D%MTrix%D(iPtrD2),NBFA,NBFB,DMab(1),&
                        &            BSc%NCFnc%I(KA),BSc%NCFnc%I(KB),      &
                        &            BSc%LStrt%I(1,KA),BSc%LStop%I(1,KA),  &
                        &            BSc%LStrt%I(1,KB),BSc%LStop%I(1,KB)   )
#endif
                   !
                   ! Get atom pair for BD.
                   CALL GetAtomPairG(BDAtmInfo,AtBList,BDAtmPair,BSc,CS_OUT)
                   !
                   NIntBlk=NBFA*NBFB*NBFC*NBFD
                   !
                   CALL DBL_VECT_EQ_DBL_SCLR(12*NIntBlk,C(1),0.0D0)
                   !STRESS STRESS STRESS STRESS STRESS STRESS STRESS
                   IF(DoStress) CALL DBL_VECT_EQ_DBL_SCLR(9*NIntBlk,CC(1),0.0D0)
                   !STRESS STRESS STRESS STRESS STRESS STRESS STRESS
                   DO iFAC=1,ACAtmInfo%NFPair
#ifdef GTRESH
                      IF(Dcd*Dab*Half*AtAList%RInt(iFAC)*AtBList%RInt(1).LT.Thresholds%TwoE) EXIT
#endif
                      CFA=AtAList%Indx(1,iFAC)
                      CFC=AtAList%Indx(2,iFAC)
                      !
                      ! BraSwitch
                      IF(ACAtmPair(iFAC)%SP%Switch) THEN
                         OAL=OffArr%I(CFC,KC)-1;LDAL=NBFA*NBFB;GOAL=4; !A
                         OBL=OffArr%I(CFA,KA)  ;LDBL=1        ;GOBL=1; !C
                      ELSE
                         OAL=OffArr%I(CFA,KA)  ;LDAL=1        ;GOAL=1; !A
                         OBL=OffArr%I(CFC,KC)-1;LDBL=NBFA*NBFB;GOBL=4; !C
                      ENDIF
                      !
                      DO iFBD=1,BDAtmInfo%NFPair 
                         CFB=AtBList%Indx(1,iFBD)
                         CFD=AtBList%Indx(2,iFBD)
#ifdef GTRESH
                         IF(DMcd((CFC-1)*NCFncD+CFD)*DMab((CFA-1)*NCFncB+CFB)*Half* &
                              & AtAList%RInt(iFAC)*AtBList%RInt(iFBD)>Thresholds%TwoE) THEN
#endif
                            !
                            ! KetSwitch
                            IF(BDAtmPair(iFBD)%SP%Switch) THEN
                               OC=OffArr%I(CFD,KD)-1;LDC=NBFA*NBFB*NBFC;GOC=10; !B
                               OD=OffArr%I(CFB,KB)-1;LDD=NBFA          ;GOD=7 ; !D
                            ELSE
                               OC=OffArr%I(CFB,KB)-1;LDC=NBFA          ;GOC=7 ; !B
                               OD=OffArr%I(CFD,KD)-1;LDD=NBFA*NBFB*NBFC;GOD=10; !D
                            ENDIF
                            !
                            ! BraKetSwitch
                            IF(ACAtmPair(iFAC)%SP%IntType.LT.BDAtmPair(iFBD)%SP%IntType) THEN
                               OA =OC ;OC =OAL ; 
                               LDA=LDC;LDC=LDAL;
                               GOA=GOC;GOC=GOAL;
                               !
                               OB =OD ;OD =OBL ;
                               LDB=LDD;LDD=LDBL;
                               GOB=GOD;GOD=GOBL;
                            ELSE
                               OA =OAL ;OB =OBL ;
                               LDA=LDAL;LDB=LDBL;
                               GOA=GOAL;GOB=GOBL;
                            ENDIF
                            !
                            ! Compute integral type.
                            IntType=ACAtmPair(iFAC)%SP%IntType*10000+BDAtmPair(iFBD)%SP%IntType
                            !
                            ! The integral interface.
                            INCLUDE 'dERIInterfaceB.Inc'
                            !
                            NInts=NInts+DBLE(LocNInt)
#ifdef GTRESH
                         ENDIF
#endif
                         !
                      ENDDO
                   ENDDO
                   !
                   !STRESS STRESS STRESS STRESS STRESS STRESS
                   IF(DoStress) THEN
                      DO IXYZ=1,3
                         DO JXYZ=1,3
                            Indx=3*NIntBlk*(IXYZ-1)+NIntBlk*(JXYZ-1)+1
                            ! Compute Stress Components.
#ifdef ONX2_PARALLEL
                            CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,CC(Indx), &
                                 &     NBFA*NBFB,U%MTrix(1,1),1,0.0d0,   &
                                 &     Work(1),1)
                            BoxX%D(IXYZ,JXYZ)=BoxX%D(IXYZ,JXYZ) &
                                 & -DDOT(NBFA*NBFB,V%MTrix(1,1),1,Work(1),1)
#else
                            CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,CC(Indx), &
                                 &     NBFA*NBFB,D%MTrix%D(iPtrD),1,0.0d0,   &
                                 &     Work(1),1)
                            BoxX%D(IXYZ,JXYZ)=BoxX%D(IXYZ,JXYZ) &
                                 & -DDOT(NBFA*NBFB,D%MTrix%D(iPtrD2),1,Work(1),1)
#endif
                            !
                         ENDDO
                      ENDDO
                   ENDIF
                   !STRESS STRESS STRESS STRESS STRESS STRESS
                   !
                   ! Compute the exchange-forces.
                   DO IXYZ=1,3
                      !
                      !-----------------------------------------------------
                      ! AtA
                      Indx=(IXYZ-1)*NIntBlk+1
                      !
#ifdef ONX2_PARALLEL
                      CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,C(Indx), &
                           &     NBFA*NBFB,U%MTrix(1,1),1,0.0d0,    &
                           &     Work(1),1)
                      TmpGradA=-DDOT(NBFA*NBFB,V%MTrix(1,1),1,Work(1),1)
#else
                      CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,C(Indx), &
                           &     NBFA*NBFB,D%MTrix%D(iPtrD),1,0.0d0,    &
                           &     Work(1),1)
                      TmpGradA=-DDOT(NBFA*NBFB,D%MTrix%D(iPtrD2),1,Work(1),1)
#endif
                      GradX%D(IXYZ,AtA)=GradX%D(IXYZ,AtA)+TmpGradA
                      !
                      !-----------------------------------------------------
                      ! AtC
                      Indx=(3+IXYZ-1)*NIntBlk+1
#ifdef ONX2_PARALLEL
                      CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,C(Indx), &
                           &     NBFA*NBFB,U%MTrix(1,1),1,0.0d0,    &
                           &     Work(1),1)
                      TmpGradC=-DDOT(NBFA*NBFB,V%MTrix(1,1),1,Work(1),1)
#else
                      CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,C(Indx), &
                           &     NBFA*NBFB,D%MTrix%D(iPtrD),1,0.0d0,    &
                           &     Work(1),1)
                      TmpGradC=-DDOT(NBFA*NBFB,D%MTrix%D(iPtrD2),1,Work(1),1)
#endif
                      GradX%D(IXYZ,AtC)=GradX%D(IXYZ,AtC)+TmpGradC
                      !
                      !-----------------------------------------------------
                      ! AtB
                      Indx=(6+IXYZ-1)*NIntBlk+1
#ifdef ONX2_PARALLEL
                      CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,C(Indx), &
                           &     NBFA*NBFB,U%MTrix(1,1),1,0.0d0,    &
                           &     Work(1),1)
                      TmpGradB=-DDOT(NBFA*NBFB,V%MTrix(1,1),1,Work(1),1)
#else
                      CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,C(Indx), &
                           &     NBFA*NBFB,D%MTrix%D(iPtrD),1,0.0d0,    &
                           &     Work(1),1)
                      TmpGradB=-DDOT(NBFA*NBFB,D%MTrix%D(iPtrD2),1,Work(1),1)
#endif
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
#ifdef ONX2_PARALLEL
          ! Set Time.
          TmEnd = MondoTimer()
          !Add Time.
          U%Part = U%Part+TmEnd-TmBeg
          U => U%Next
#endif
          !
       ENDDO ! End AtD
       !
#ifdef ONX2_PARALLEL
       P => P%Next
#endif
       !
    ENDDO ! End AtC
    !
    CALL Delete(BColIdx)
    !
    ! DeAllocate arrays.
    DEALLOCATE(ACAtmPair,BDAtmPair,STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In ComputDK: Deallocation problem.')
    !
    !
    WRITE(*,100) NInts,12D0*CS_OUT%NCells**2*DBLE(NBasF)**4, &
         &       NInts/(12D0*CS_OUT%NCells**2*DBLE(NBasF)**4)*100D0
100 FORMAT(' NInts = ',E8.2,' NIntTot = ',E8.2,' Ratio = ',E8.2,'%')
    !
  END SUBROUTINE ComputDK
  !
  !
  SUBROUTINE GetAtomPairG(AtmInfo,List,AtmPair,BSc,CS_OUT)
!H---------------------------------------------------------------------------------
!H SUBROUTINE GetAtomPairG(AtmInfo,List,AtmPair,BSc,CS_OUT)
!H
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(AtomInfo)               :: AtmInfo
    TYPE(ANode)   , POINTER      :: List
    TYPE(AtomPrG) , DIMENSION(:) :: AtmPair
    TYPE(BSET)                   :: BSc
    TYPE(CellSet)                :: CS_OUT
    !-------------------------------------------------------------------
    INTEGER      :: CF1,CF2,I1,I2,II,JJ,IJ,Cell
    INTEGER      :: MinL1,MaxL1,Type1,MinL2,MaxL2,Type2
    INTEGER      :: StopL1,StartL1,StopL2,StartL2
    INTEGER      :: ISwitch1,ISwitch2,iNFPair
    REAL(DOUBLE) :: Z1,Z2,Expt,InvExpt,R12,XiR12,RX,RY,RZ
    LOGICAL      :: Switch
    !-------------------------------------------------------------------
    !
    ! dbuging
#ifdef GONX2_DBUG
    AtmPair(:)%SP%IntType=BIG_INT
#endif
    !
    DO iNFPair=1,AtmInfo%NFPair
       CF1 =List%Indx(1,iNFPair) !A,B
       CF2 =List%Indx(2,iNFPair) !C,D
       Cell=List%Indx(3,iNFPair)
       RX=CS_OUT%CellCarts%D(1,Cell)
       RY=CS_OUT%CellCarts%D(2,Cell)
       RZ=CS_OUT%CellCarts%D(3,Cell)
       !
       ! AtmInfo must be related to the atoms in the working cell ONLY.
       ! Then we add the PBC's to have the right interatomic distance.
       R12=(AtmInfo%Atm12X-RX)**2+(AtmInfo%Atm12Y-RY)**2+(AtmInfo%Atm12Z-RZ)**2
       !
       MinL1=BSc%ASymm%I(1,CF1,AtmInfo%K1)
       MaxL1=BSc%ASymm%I(2,CF1,AtmInfo%K1)
       Type1=MaxL1*(MaxL1+1)/2+MinL1+1
       StartL1=BSc%LStrt%I(CF1,AtmInfo%K1)
       StopL1=BSc%LStop%I(CF1,AtmInfo%K1)
       !
       if(Type1==2) stop 'SP shell not yet supported.'
       !
       MinL2=BSc%ASymm%I(1,CF2,AtmInfo%K2)
       MaxL2=BSc%ASymm%I(2,CF2,AtmInfo%K2)
       Type2=MaxL2*(MaxL2+1)/2+MinL2+1
       StartL2=BSc%LStrt%I(CF2,AtmInfo%K2)
       StopL2=BSc%LStop%I(CF2,AtmInfo%K2)
       !
       if(Type2==2) stop 'SP shell not yet supported.'
       !
       !>new
       Switch=Type1.LT.Type2
       AtmPair(iNFPair)%SP%Switch=Switch
       IF(Switch) THEN
          ISwitch1=10
          ISwitch2=9
          AtmPair(iNFPair)%SP%IntType=Type2*100+Type1
       ELSE
          ISwitch1=9
          ISwitch2=10
          AtmPair(iNFPair)%SP%IntType=Type1*100+Type2
       ENDIF
       !<<<new
       !oAtmPair(iNFPair)%SP%IntType=Type1*100+Type2
       !o!
       !o! .NOT.(Do we need to switch the centers?)
       !oSwitch=Type1.GE.Type2
       !o!
       !oIF(Switch) THEN
       !o   ISwitch1=9
       !o   ISwitch2=10
       !o!          ISwitch1=6
       !o!          ISwitch2=7
       !oELSE
       !o!          ISwitch1=7
       !o!          ISwitch2=6
       !o   ISwitch1=10
       !o   ISwitch2=9
       !oENDIF
       !
       II=0
       !
       ! We assume the primitives are ordered (exponants in decressing order).
       DO I1=BSc%NPFnc%I(CF1,AtmInfo%K1),1,-1
          Z1=BSc%Expnt%D(I1,CF1,AtmInfo%K1)
          JJ=0
          !
          ! We assume the primitives are ordered (exponants in decressing order).
          DO I2=BSc%NPFnc%I(CF2,AtmInfo%K2),1,-1
             Z2=BSc%Expnt%D(I2,CF2,AtmInfo%K2)
             Expt=Z1+Z2
             InvExpt=1.0d0/Expt
             XiR12=Z2*Z1*InvExpt*R12
             IF(XiR12<PrimPairDistanceThreshold) THEN
                JJ=JJ+1
                IJ=JJ+II
                AtmPair(iNFPair)%SP%Cst(1,IJ)=Expt
                !
                ! AtmInfo must be related to the atoms in the working cell ONLY.
                ! Then we add the PBC's to have the right atomic position.
                AtmPair(iNFPair)%SP%Cst(2,IJ)=(Z1*AtmInfo%Atm1X+Z2*(AtmInfo%Atm2X+RX))*InvExpt
                AtmPair(iNFPair)%SP%Cst(3,IJ)=(Z1*AtmInfo%Atm1Y+Z2*(AtmInfo%Atm2Y+RY))*InvExpt
                AtmPair(iNFPair)%SP%Cst(4,IJ)=(Z1*AtmInfo%Atm1Z+Z2*(AtmInfo%Atm2Z+RZ))*InvExpt
                AtmPair(iNFPair)%SP%Cst(5,IJ)=5.914967172796D0*EXP(-XiR12)*InvExpt* &
                     !AtmPair(CF12)%SP%Cst(5,IJ)=5.914967172796D0*EXPInv(XiR12)*InvExpt* &
                     &                     BSc%CCoef%D(StopL1,I1,CF1,AtmInfo%K1)* &
                     &                     BSc%CCoef%D(StopL2,I2,CF2,AtmInfo%K2)
                !
                AtmPair(iNFPair)%SP%Cst(ISwitch1,IJ)=2.0d0*Z1!->ISwitch1=6,7
                AtmPair(iNFPair)%SP%Cst(ISwitch2,IJ)=2.0d0*Z2!->ISwitch2=7,6
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
       AtmPair(iNFPair)%SP%L=II
       !
       ! We reorder the atomic positions if Type2 > Type1.
       ! Needed for the integral evaluations.
       IF(Type1.GE.Type2) THEN
       !oIF(Switch) THEN
          AtmPair(iNFPair)%SP%AtmInfo%Atm1X=AtmInfo%Atm1X
          AtmPair(iNFPair)%SP%AtmInfo%Atm1Y=AtmInfo%Atm1Y
          AtmPair(iNFPair)%SP%AtmInfo%Atm1Z=AtmInfo%Atm1Z
          !
          ! AtmInfo must be related to the atoms in the working cell ONLY.
          ! Then we add the PBC's to have the right atomic position.
          AtmPair(iNFPair)%SP%AtmInfo%Atm2X=AtmInfo%Atm2X+RX
          AtmPair(iNFPair)%SP%AtmInfo%Atm2Y=AtmInfo%Atm2Y+RY
          AtmPair(iNFPair)%SP%AtmInfo%Atm2Z=AtmInfo%Atm2Z+RZ
       ELSE
          !
          ! AtmInfo must be related to the atoms in the working cell ONLY.
          ! Then we add the PBC's to have the right atomic position.
          AtmPair(iNFPair)%SP%AtmInfo%Atm1X=AtmInfo%Atm2X+RX
          AtmPair(iNFPair)%SP%AtmInfo%Atm1Y=AtmInfo%Atm2Y+RY
          AtmPair(iNFPair)%SP%AtmInfo%Atm1Z=AtmInfo%Atm2Z+RZ
          AtmPair(iNFPair)%SP%AtmInfo%Atm2X=AtmInfo%Atm1X
          AtmPair(iNFPair)%SP%AtmInfo%Atm2Y=AtmInfo%Atm1Y
          AtmPair(iNFPair)%SP%AtmInfo%Atm2Z=AtmInfo%Atm1Z
       ENDIF
       !
#ifdef GONX2_DBUG
       WRITE(*,'(3(A,I3),A,I5)') 'Type1',Type1,' Type2',Type2, &
            &                    ' IntType',AtmPair(iNFPair)%SP%IntType
#endif
    ENDDO
    !
  END SUBROUTINE GetAtomPairG
  !
  !
  INTEGER FUNCTION GetNonNFPair(List,DFac,Trsh)
!H---------------------------------------------------------------------------------
!H INTEGER FUNCTION GetNonNFPair(List,DFac,Trsh)
!H
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(ANode ), POINTER    :: List
    REAL(DOUBLE), INTENT(IN) :: DFac,Trsh
    !-------------------------------------------------------------------
    INTEGER                  :: I
    !-------------------------------------------------------------------
    !
    GetNonNFPair=0
    !
    DO I=1,List%NFPair
       IF(List%RInt(I)*DFac.LE.Trsh) EXIT
       GetNonNFPair=GetNonNFPair+1
    ENDDO
    !
  END FUNCTION GetNonNFPair
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
    INTEGER, INTENT(IN) :: At
    TYPE(BCSR)          :: D
    TYPE(INT_VECT)      :: ColIdx
    !-------------------------------------------------------------------
    INTEGER             :: Ci
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











