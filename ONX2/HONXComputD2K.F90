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
!H  DEBUGING: Use -DHONX_DBUG to print some stuff.
!H  INFO    : Use -DHONX_INFO to print some stuff.
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
  PUBLIC  :: ComputD2K
  !
!--------------------------------------------------------------------------------- 
! PRIVATE DECLARATIONS
!--------------------------------------------------------------------------------- 
  PRIVATE :: GetAtomPairH
  !
CONTAINS
  !
  !
#ifdef ONX2_PARALLEL
  SUBROUTINE ComputD2K(DFMcd,DFMab,HessX,MixHessX,LatHessX,DoLHess,ListC,ListD, &
       &               OffArr,GMc,BSc,CS_OUT)
#else
  SUBROUTINE ComputD2K(D          ,HessX,MixHessX,LatHessX,DoLHess,ListC,ListD, &
       &               OffArr,GMc,BSc,CS_OUT)
#endif
!H---------------------------------------------------------------------------------
!H SUBROUTINE ComputD2K(D,HessX,MixHessX,LatHessX,DoLHess,ListC,ListD,OffArr,GMc,BSc,CS_OUT)
!H
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !
    ! HessX = Dcd*(ac(R)|bd(R'))''*Dab
    !
    !-------------------------------------------------------------------
#ifdef ONX2_PARALLEL
    TYPE(FASTMAT)              , POINTER :: DFMcd,DFMab
    TYPE(FASTMAT)              , POINTER :: P,Q
    TYPE(SRST   )              , POINTER :: U,V
#else
    TYPE(BCSR)                           :: D
#endif
    TYPE(DBL_RNK2)                       :: HessX,MixHessX,LatHessX
    TYPE(INT_RNK2)                       :: OffArr
    TYPE(CList) , DIMENSION(:) , POINTER :: ListC,ListD
    TYPE(CRDS)                           :: GMc
    TYPE(BSET)                           :: BSc
    TYPE(CellSet)                        :: CS_OUT
    LOGICAL, INTENT(IN)                  :: DoLHess
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
    REAL(DOUBLE)               :: TmpGradA,TmpGradC,TmpGradB,TmpHess(12,12)!TmpHess(78)
    REAL(DOUBLE)               :: Dcd,Dab,NInts,NIntsTot
    !-------------------------------------------------------------------
    REAL(DOUBLE) , DIMENSION(MaxFuncPerAtmBlk**2) :: Work
    REAL(DOUBLE) , DIMENSION(MaxShelPerAtmBlk**2) :: DMcd,DMab
    REAL(DOUBLE) , DIMENSION(78*MaxFuncPerAtmBlk**4) :: C
    !LATTICE HESSIAN LATTICE HESSIAN LATTICE HESSIAN LATTICE HESSIAN 
    REAL(DOUBLE) , DIMENSION( 9*MaxFuncPerAtmBlk**4) :: CC
    !LATTICE HESSIAN LATTICE HESSIAN LATTICE HESSIAN LATTICE HESSIAN 
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
    IF(iErr.NE.0) CALL Halt('In ComputD2K: Allocation problem.')
    !
    !
    !Simple check Simple check Simple check Simple check
    isize=0
    do i=1,natoms
       isize=MAX(isize,BSc%BfKnd%I(GMc%AtTyp%I(i)))
       IF(78*(BSc%BfKnd%I(GMc%AtTyp%I(i)))**4.GT.SIZE(C)) THEN
          write(*,*) 'size',78*BSc%BfKnd%I(GMc%AtTyp%I(i))**4
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
    !Check for wrong InvBoxSh and BoxShape
    IF(    ABS(GMc%PBC%InvBoxSh%D(2,1)).GT.1D-15.OR.ABS(GMc%PBC%InvBoxSh%D(3,1)).GT.1D-15.OR.&
         & ABS(GMc%PBC%InvBoxSh%D(3,2)).GT.1D-15.OR.ABS(GMc%PBC%BoxShape%D(2,1)).GT.1D-15.OR.&
         & ABS(GMc%PBC%BoxShape%D(3,1)).GT.1D-15.OR.ABS(GMc%PBC%BoxShape%D(3,2)).GT.1D-15) THEN
       WRITE(*,*) 'The following guys MUST be ZERO!'
       WRITE(*,*) 'InvBoxSh:',GMc%PBC%InvBoxSh%D(2,1),GMc%PBC%InvBoxSh%D(3,1),GMc%PBC%InvBoxSh%D(3,2)
       WRITE(*,*) 'BoxShape:',GMc%PBC%BoxShape%D(2,1),GMc%PBC%BoxShape%D(3,1),GMc%PBC%BoxShape%D(3,2)
       STOP 'STOP in HONXComputD2K.F90'
    ENDIF
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
#ifdef HONX_DBUG
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
             CALL GetAtomPairH(ACAtmInfo,AtAList,ACAtmPair,BSc,CS_OUT)
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
                   CALL GetAtomPairH(BDAtmInfo,AtBList,BDAtmPair,BSc,CS_OUT)
                   !
                   NIntBlk=NBFA*NBFB*NBFC*NBFD
                   !
                   CALL DBL_VECT_EQ_DBL_SCLR(78*NIntBlk,C(1),0.0D0)
                   !LATTICE HESSIAN LATTICE HESSIAN LATTICE HESSIAN LATTICE HESSIAN 
                   IF(DoLHess) CALL DBL_VECT_EQ_DBL_SCLR(9*NIntBlk,CC(1),0.0D0)
                   !LATTICE HESSIAN LATTICE HESSIAN LATTICE HESSIAN LATTICE HESSIAN 
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
                            IF(IntType.EQ.1010101) THEN
                               CALL hIntB1010101(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                                    BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                                    ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                                    OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NIntBlk, &
                                    GMc%PBC,C(1))
                               write(*,*) C(1)
                            ELSE
                               write(*,*) 'Doesn''t support that integral type:',IntType
                               STOP 'In ComputD2K'
                            ENDIF
                            !HINCLUDE 'hERIInterfaceB.Inc'
                            !
                            NInts=NInts+DBLE(LocNInt)
#ifdef GTRESH
                         ENDIF
#endif
                         !
                      ENDDO
                   ENDDO
                   !
                   !LATTICE HESSIAN LATTICE HESSIAN LATTICE HESSIAN LATTICE HESSIAN 
!!$                   IF(DoLHess) THEN
!!$                      DO JXYZ=1,3
!!$                         DO IXYZ=1,3
!!$                            IF(GMc%PBC%AutoW%I(IXYZ).EQ.1.AND.GMc%PBC%AutoW%I(JXYZ).EQ.1) THEN
!!$                               Indx=3*NIntBlk*(JXYZ-1)+NIntBlk*(IXYZ-1)+1
!!$                               ! Compute Stress Components.
!!$#ifdef ONX2_PARALLEL
!!$                               CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,CC(Indx), &
!!$                                    &     NBFA*NBFB,U%MTrix(1,1),1,0.0d0,   &
!!$                                    &     Work(1),1)
!!$                               LatHessX%D(IXYZ,JXYZ)=LatHessX%D(IXYZ,JXYZ) &
!!$                                    & -DDOT(NBFA*NBFB,V%MTrix(1,1),1,Work(1),1)
!!$#else
!!$                               !write(*,*) 'CC(Indx)',CC(Indx)
!!$                               CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,CC(Indx), &
!!$                                    &     NBFA*NBFB,D%MTrix%D(iPtrD),1,0.0d0,   &
!!$                                    &     Work(1),1)
!!$                               LatHessX%D(IXYZ,JXYZ)=LatHessX%D(IXYZ,JXYZ) &
!!$                                    & -DDOT(NBFA*NBFB,D%MTrix%D(iPtrD2),1,Work(1),1)
!!$#endif
!!$                            !
!!$                            ENDIF
!!$                         ENDDO
!!$                      ENDDO
!!$                   ENDIF
                   !LATTICE HESSIAN LATTICE HESSIAN LATTICE HESSIAN LATTICE HESSIAN 
                   !
!!$                   ! Compute the exchange-Hessian.
!!$                   DO jAtm=1,3
!!$                   DO JXYZ=1,3
!!$                      DO iAtm=jAtm,3
!!$                         IMin=1
!!$                         IF(jAtm.EQ.iAtm)IMin=JXYZ
!!$                         DO IXYZ=1,3
!!$
!!$                         !TODO TODO TODO TODO TODO TODO 
!!$                         Indx=1!TODO TODO TODO TODO TODO
!!$                         !TODO TODO TODO TODO TODO TODO  
!!$                         ii=3*(iAtm-1)+IXYZ
!!$                         jj=3*(jAtm-1)+JXYZ
!!$                         !ij=ii+(2*12-jj)*(jj-1)/2
!!$#ifdef ONX2_PARALLEL
!!$                         CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,C(Indx), &
!!$                              &     NBFA*NBFB,U%MTrix(1,1),1,0.0d0,    &
!!$                              &     Work(1),1)
!!$                         TmpHess(ii,jj)=-DDOT(NBFA*NBFB,V%MTrix(1,1),1,Work(1),1)
!!$#else
!!$                         CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,C(Indx), &
!!$                              &     NBFA*NBFB,D%MTrix%D(iPtrD),1,0.0d0,    &
!!$                              &     Work(1),1)
!!$                         TmpHess(ii,jj)=-DDOT(NBFA*NBFB,D%MTrix%D(iPtrD2),1,Work(1),1)
!!$#endif
!!$                         TmpHess(jj,ii)=TmpHess(ii,jj)
!!$                      ENDDO
!!$                      ENDDO
!!$                   ENDDO
!!$                   ENDDO
!!$                   
!!$                   APt=3*(AtA-1)
!!$                   CPt=3*(AtC-1)
!!$                   BPt=3*(AtB-1)
!!$                   DPt=3*(AtD-1)
!!$
!!$                   DO J=1,3
!!$                   DO I=1,3
!!$                      !
!!$                      iAtm=0;jAtm=0
!!$                      ii=3*(iAtm-1)+I;jj=3*(jAtm-1)+J
!!$                      HessX%D(APt+I,APt+J)=HessX%D(APt+I,APt+J)+TmpHess(ii,jj)
!!$                      !
!!$                      iAtm=1;jAtm=0
!!$                      ii=3*(iAtm-1)+I;jj=3*(jAtm-1)+J
!!$                      HessX%D(CPt+I,APt+J)=HessX%D(CPt+I,APt+J)+TmpHess(ii,jj)
!!$                      !
!!$                      iAtm=2;jAtm=0
!!$                      ii=3*(iAtm-1)+I;jj=3*(jAtm-1)+J
!!$                      HessX%D(BPt+I,APt+J)=HessX%D(Bpt+I,Apt+J)+TmpHess(ii,jj)
!!$                      !
!!$                      iAtm=3;jAtm=0
!!$                      ii=3*(iAtm-1)+I;jj=3*(jAtm-1)+J
!!$                      HessX%D(DPt+I,APt+J)=HessX%D(Dpt+I,Apt+J)+TmpHess(ii,jj)
!!$                      !
!!$                      
!!$                      !....
!!$                   ENDDO
!!$                   ENDDO



                   DO IXYZ=1,3
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
                      HessX%D(IXYZ,AtA)=HessX%D(IXYZ,AtA)+TmpGradA
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
                      HessX%D(IXYZ,AtC)=HessX%D(IXYZ,AtC)+TmpGradC
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
                      HessX%D(IXYZ,AtB)=HessX%D(IXYZ,AtB)+TmpGradB
                      !
                      !-----------------------------------------------------
                      ! AtD
                      HessX%D(IXYZ,AtD)=HessX%D(IXYZ,AtD)-(TmpGradA+TmpGradC+TmpGradB)
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
    IF(iErr.NE.0) CALL Halt('In ComputD2K: Deallocation problem.')
    !
    !
    WRITE(*,100) NInts,12D0*CS_OUT%NCells**2*DBLE(NBasF)**4, &
         &       NInts/(12D0*CS_OUT%NCells**2*DBLE(NBasF)**4)*100D0
100 FORMAT(' NInts = ',E8.2,' NIntTot = ',E8.2,' Ratio = ',E8.2,'%')
    !
  END SUBROUTINE ComputD2K
  !
  !
! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (s s|s s) 
! ---------------------------------------------------------- 
  SUBROUTINE hIntB1010101(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, &
       OA,LDA,OB,LDB,OC,LDC,OD,LDD,GOA,GOB,GOC,GOD,NINT,PBC,HESSIAN)
    USE DerivedTypes
    USE VScratchB
    USE GlobalScalars
    USE ShellPairStruct
    USE GammaF2
    IMPLICIT REAL(DOUBLE) (W)
    INTEGER        :: LBra,LKet,NINT,CDOffSet
    REAL(DOUBLE)   :: PrmBufB(10,LBra),PrmBufK(10,LKet)
    TYPE(SmallAtomInfo) :: ACInfo,BDInfo
    TYPE(PBCInfo) :: PBC
    REAL(DOUBLE)  :: HESSIAN(NINT,12)
    REAL(DOUBLE)  :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz
    REAL(DOUBLE)  :: Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz
    REAL(DOUBLE)  :: PQx,PQy,PQz,FPQx,FPQy,FPQz
    REAL(DOUBLE)  :: Zeta,Eta,Omega,Up,Uq,Upq
    REAL(DOUBLE)  :: T,ET,TwoT,InvT,SqInvT
    REAL(DOUBLE)  :: Alpha,Beta,Gamma
    REAL(DOUBLE), DIMENSION(4) :: HRRTmp 
    REAL(DOUBLE), DIMENSION(1,1,1) :: HRR 
    REAL(DOUBLE), DIMENSION(4,1,1) :: HRRA,HRRB 
    REAL(DOUBLE), DIMENSION(1,4,1) :: HRRC 
    real(double) HRRAA(10,1,1),HRRBB(10,1,1),HRRCC(1,10,1),HRRAB(10,1,1),HRRAC(4,4,1),HRRBC(4,4,1)
    real(double)  :: AlpAlp,BetBet,GamGam,AlpBet,AlpGam,BetGam
    REAL(DOUBLE)  :: VRR(4,4,0:1)
    REAL(DOUBLE)  :: TOm,PQJ(3),FP(9)
    INTEGER       :: OffSet,OA,LDA,GOA,OB,LDB,GOB,OC,LDC,GOC,OD,LDD,GOD,I,J,K,L,IJ,i1,i2
    EXTERNAL InitDbl
    HRR  =BIG_DBL
    HRRA =BIG_DBL
    HRRB =BIG_DBL
    HRRC =BIG_DBL
    HRRAA=BIG_DBL
    HRRBB=BIG_DBL
    HRRCC=BIG_DBL
    HRRAB=BIG_DBL
    HRRAC=BIG_DBL
    HRRBC=BIG_DBL
    !
    CALL InitDbl(1*1,HRR(1,1,1))
    CALL InitDbl(4*1,HRRA(1,1,1))
    CALL InitDbl(4*1,HRRB(1,1,1))
    CALL InitDbl(1*4,HRRC(1,1,1))
    !
    CALL InitDbl(10,HRRAA(1,1,1))
    CALL InitDbl(10,HRRBB(1,1,1))
    CALL InitDbl(10,HRRCC(1,1,1))
    CALL InitDbl(10,HRRAB(1,1,1))
    CALL InitDbl(16,HRRAC(1,1,1))
    CALL InitDbl(16,HRRBC(1,1,1))
    !
    Ax=ACInfo%Atm1X
    Ay=ACInfo%Atm1Y
    Az=ACInfo%Atm1Z
    Bx=ACInfo%Atm2X
    By=ACInfo%Atm2Y
    Bz=ACInfo%Atm2Z
    Cx=BDInfo%Atm1X
    Cy=BDInfo%Atm1Y
    Cz=BDInfo%Atm1Z
    Dx=BDInfo%Atm2X
    Dy=BDInfo%Atm2Y
    Dz=BDInfo%Atm2Z
    ABx=Ax-Bx
    ABy=Ay-By
    ABz=Az-Bz
    CDx=Cx-Dx
    CDy=Cy-Dy
    CDz=Cz-Dz
    !
    DO J=1,LKet ! K^2 VRR |N0) loop 
       Eta=PrmBufK(1,J)
       Qx=PrmBufK(2,J)
       Qy=PrmBufK(3,J)
       Qz=PrmBufK(4,J)
       Uq=PrmBufK(5,J)
       Gamma =PrmBufK(9,J)
       QCx=Qx-Cx
       QCy=Qy-Cy
       QCz=Qz-Cz
       DO K=1,LBra ! K^2 VRR (M0| loop 
          Zeta=PrmBufB(1,K)
          Px=PrmBufB(2,K)
          Py=PrmBufB(3,K)
          Pz=PrmBufB(4,K)
          Up=PrmBufB(5,K)
          Alpha =PrmBufB(9,K)
          Beta  =PrmBufB(10,K)
          r1xZpE=One/(Zeta+Eta)
          Upq=SQRT(r1xZpE)*Up*Uq
          HfxZpE=Half/(Zeta+Eta)
          r1x2E=Half/Eta
          r1x2Z=Half/Zeta
          ExZpE=Eta*r1xZpE
          ZxZpE=Zeta*r1xZpE
          Omega=Eta*Zeta*r1xZpE
          PAx=Px-Ax
          PAy=Py-Ay
          PAz=Pz-Az
          PQx=Px-Qx
          PQy=Py-Qy
          PQz=Pz-Qz
          !
          WPx = -Eta*PQx*r1xZpE
          WPy = -Eta*PQy*r1xZpE
          WPz = -Eta*PQz*r1xZpE
          WQx = Zeta*PQx*r1xZpE
          WQy = Zeta*PQy*r1xZpE
          WQz = Zeta*PQz*r1xZpE
          T=Omega*(PQx*PQx+PQy*PQy+PQz*PQz)
          IF(T<Gamma_Switch)THEN
             L=AINT(T*Gamma_Grid)
             ET=EXP(-T)
             TwoT=Two*T
             W2=(F2_0(L)+T*(F2_1(L)+T*(F2_2(L)+T*(F2_3(L)+T*F2_4(L)))))
             W1=+3.333333333333333D-01*(TwoT*W2+ET)
             W0=TwoT*W1+ET
             VRR(1,1,0)=Upq*W0
             VRR(1,1,1)=Upq*W1
             VRR(1,1,2)=Upq*W2
          ELSE
             InvT=One/T
             SqInvT=DSQRT(InvT)
             VRR(1,1,0)=+8.862269254527580D-01*Upq*SqInvT
             SqInvT=SqInvT*InvT
             VRR(1,1,1)=+4.431134627263790D-01*Upq*SqInvT
             SqInvT=SqInvT*InvT
             VRR(1,1,2)=+6.646701940895685D-01*Upq*SqInvT
          ENDIF
          ! Generating (p0|s0)^(1)
          VRR(2,1,1)=PAx*VRR(1,1,1)+WPx*VRR(1,1,2) 
          VRR(3,1,1)=PAy*VRR(1,1,1)+WPy*VRR(1,1,2) 
          VRR(4,1,1)=PAz*VRR(1,1,1)+WPz*VRR(1,1,2) 
          ! Generating (p0|s0)^(0)
          VRR(2,1,0)=PAx*VRR(1,1,0)+WPx*VRR(1,1,1) 
          VRR(3,1,0)=PAy*VRR(1,1,0)+WPy*VRR(1,1,1) 
          VRR(4,1,0)=PAz*VRR(1,1,0)+WPz*VRR(1,1,1) 
          ! Generating (d0|s0)^(0)
          VRR( 5,1,0)=PAx*VRR(2,1,0)+r1x2Z*(VRR(1,1,0)-ExZpE*VRR(1,1,1))+WPx*VRR(2,1,1)
          VRR( 6,1,0)=PAx*VRR(3,1,0)+WPx*VRR(3,1,1)
          VRR( 7,1,0)=PAy*VRR(3,1,0)+r1x2Z*(VRR(1,1,0)-ExZpE*VRR(1,1,1))+WPy*VRR(3,1,1)
          VRR( 8,1,0)=PAx*VRR(4,1,0)+WPx*VRR(4,1,1)
          VRR( 9,1,0)=PAy*VRR(4,1,0)+WPy*VRR(4,1,1)
          VRR(10,1,0)=PAz*VRR(4,1,0)+r1x2Z*(VRR(1,1,0)-ExZpE*VRR(1,1,1))+WPz*VRR(4,1,1)
          ! Generating (s0|p0)^(1)
          VRR(1,2,1)=QCx*VRR(1,1,1)+WQx*VRR(1,1,2)
          VRR(1,3,1)=QCy*VRR(1,1,1)+WQy*VRR(1,1,2)
          VRR(1,4,1)=QCz*VRR(1,1,1)+WQz*VRR(1,1,2)
          ! Generating (s0|p0)^(0)
          VRR(1,2,0)=QCx*VRR(1,1,0)+WQx*VRR(1,1,1)
          VRR(1,3,0)=QCy*VRR(1,1,0)+WQy*VRR(1,1,1)
          VRR(1,4,0)=QCz*VRR(1,1,0)+WQz*VRR(1,1,1)
          ! Generating (s0|d0)^(0)
          VRR(1, 5,0)=r1x2E*VRR(1,1,0)+QCx*VRR(1,2,0)-r1x2E*ZxZpE*VRR(1,1,1)+WQx*VRR(1,2,1)
          VRR(1, 6,0)=QCx*VRR(1,3,0)+WQx*VRR(1,3,1)
          VRR(1, 7,0)=r1x2E*VRR(1,1,0)+QCy*VRR(1,3,0)-r1x2E*ZxZpE*VRR(1,1,1)+WQy*VRR(1,3,1)
          VRR(1, 8,0)=QCx*VRR(1,4,0)+WQx*VRR(1,4,1)
          VRR(1, 9,0)=QCy*VRR(1,4,0)+WQy*VRR(1,4,1)
          VRR(1,10,0)=r1x2E*VRR(1,1,0)+QCz*VRR(1,4,0)-r1x2E*ZxZpE*VRR(1,1,1)+WQz*VRR(1,4,1)
          ! Generating (p0|p0)^(0)
          VRR(2,2,0)=QCx*VRR(2,1,0)+HfxZpE*VRR(1,1,1)+WQx*VRR(2,1,1)
          VRR(2,3,0)=QCy*VRR(2,1,0)+WQy*VRR(2,1,1)
          VRR(2,4,0)=QCz*VRR(2,1,0)+WQz*VRR(2,1,1)
          VRR(3,2,0)=QCx*VRR(3,1,0)+WQx*VRR(3,1,1)
          VRR(3,3,0)=QCy*VRR(3,1,0)+HfxZpE*VRR(1,1,1)+WQy*VRR(3,1,1)
          VRR(3,4,0)=QCz*VRR(3,1,0)+WQz*VRR(3,1,1)
          VRR(4,2,0)=QCx*VRR(4,1,0)+WQx*VRR(4,1,1)
          VRR(4,3,0)=QCy*VRR(4,1,0)+WQy*VRR(4,1,1)
          VRR(4,4,0)=QCz*VRR(4,1,0)+HfxZpE*VRR(1,1,1)+WQz*VRR(4,1,1)
          ! Contracting ... 
          AlpAlp=Alpha*Alpha
          BetBet=Beta *Beta
          GamGam=Gamma*Gamma
          AlpBet=Alpha*Beta
          AlpGam=Alpha*Gamma
          BetGam=Beta *Gamma
          HRR(1,1,1)=HRR(1,1,1)+VRR(1,1,0)
          DO I1=1,10
             HRRAA(I1,1,1)=HRRAA(I1,1,1)+AlpAlp*VRR(I1,1,0)
             HRRBB(I1,1,1)=HRRBB(I1,1,1)+BetBet*VRR(I1,1,0)
             HRRAB(I1,1,1)=HRRAB(I1,1,1)+AlpBet*VRR(I1,1,0)
             HRRCC(1,I1,1)=HRRCC(1,I1,1)+GamGam*VRR(1,I1,0)
          ENDDO
          DO I2=1,4
             HRRA(I2,1,1)=HRRA(I2,1,1)+Alpha*VRR(I2,1,0)
             HRRB(I2,1,1)=HRRB(I2,1,1)+Beta *VRR(I2,1,0)
             HRRC(1,I2,1)=HRRC(1,I2,1)+Gamma*VRR(1,I2,0)
             DO I1=1,4
                HRRAC(I1,I2,1)=HRRAC(I1,I2,1)+AlpGam*VRR(I1,I2,0)
                HRRBC(I1,I2,1)=HRRBC(I1,I2,1)+BetGam*VRR(I1,I2,0)
             ENDDO
          ENDDO
          !
       ENDDO ! (M0| loop
    ENDDO ! |N0) loop

    



    !TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
    !Nedd to collect integrals in HESSIAN(...)
    !TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
    
  END SUBROUTINE hIntB1010101
  !
  !
  SUBROUTINE GetAtomPairH(AtmInfo,List,AtmPair,BSc,CS_OUT)
!H---------------------------------------------------------------------------------
!H SUBROUTINE GetAtomPairH(AtmInfo,List,AtmPair,BSc,CS_OUT)
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
#ifdef HONX_DBUG
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
#ifdef HONX_DBUG
       WRITE(*,'(3(A,I3),A,I5)') 'Type1',Type1,' Type2',Type2, &
            &                    ' IntType',AtmPair(iNFPair)%SP%IntType
#endif
    ENDDO
    !
  END SUBROUTINE GetAtomPairH
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
END MODULE HONXComputD2K











