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
#undef DIAG_QCTC

MODULE ParallelQCTC
  USE Globals
  USE MondoMPI
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalObjects
  USE ERIGlobals
  USE InOut
  USE PoleTree
  IMPLICIT NONE

#ifdef PARALLEL
  TYPE(BBox) :: EBB_GlobalBox  ! GlobalBoundingBox
  TYPE(BBox) :: RootBox ! Root bounding box
  REAL(DOUBLE) :: BP(6)
  REAL(DOUBLE),PARAMETER :: EP = 1.0D-8
  TYPE(DBL_RNK2) :: LC,RC
  REAL(DOUBLE) :: RootBDist
  TYPE(INT_VECT) :: BegNQInd,EndNQInd
  TYPE(INT_VECT) :: BegPrInd,EndPrInd
  TYPE(INT_VECT) :: BegAtInd,EndAtInd
  INTEGER :: BranchTier
  TYPE PointerArray
    TYPE(PoleNode),POINTER :: Ptr
  END TYPE PointerArray
  TYPE(PointerArray),ALLOCATABLE :: PA1(:),PA2(:)
  INTEGER :: RunPointer
  INTEGER :: NodesVisit,NodesVisit1,IntNum,DblNum
  INTEGER :: GIntNum,GDblNum

  INTEGER,ALLOCATABLE :: IntArr(:)
  REAL(DOUBLE),ALLOCATABLE :: DblArr(:)
  INTEGER :: IntIndex,DblIndex

  INTEGER,ALLOCATABLE :: RecIntArr(:)
  REAL(DOUBLE),ALLOCATABLE :: RecDblArr(:)

  INTEGER :: TotPrCount
  TYPE(DBL_RNK2) ::  PosTimePair
  INTEGER :: GlobalCount,GQLineLoc
CONTAINS
!------------------------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------------------------
  SUBROUTINE GetDistrRho(Name,Args,SCFCycle)
    TYPE(ARGMT) :: Args
    INTEGER :: NNDist,TotNNDist,NDist,SCFCycle,IOS,I,LMNLen,ENExpt,TotENExpt,&
               ENDist,TotENDist,ENCoef,TotENCoef,NewNExpt,TotNDistm,&
               NNCoef,TotNNCoef,TotNCoef,NewNDist,NewNCoef,TotNDist,IErr,Sumi

    REAL(DOUBLE) :: SumECoef,TotSumECoef
    CHARACTER(LEN=*) :: Name
    CHARACTER(Len=DCL) :: Filename
    TYPE(HGRho) :: Tmp
    TYPE(INT_VECT) :: CSize,disp

    Filename = TrixFile('Rho'//IntToChar(MyID),Args,SCFCycle)
    Open(Unit=Seq,File=TRIM(Filename),Status='old',form='unformatted',access='sequential')

    Read(Unit=Seq,Err=202,IOSTAT=IOS) Tmp%NSDen,Tmp%NExpt,Tmp%NDist,Tmp%NCoef!<<<SPIN
    CALL New(Tmp%NQ  ,Tmp%NExpt)
    CALL New(Tmp%OffQ,Tmp%NExpt)
    CALL New(Tmp%OffR,Tmp%NExpt)
    CALL New(Tmp%Lndx,Tmp%NExpt)
    CALL New(Tmp%Expt,Tmp%NExpt)
    CALL New(Tmp%Qx,  Tmp%NDist)
    CALL New(Tmp%Qy,  Tmp%NDist)
    CALL New(Tmp%Qz,  Tmp%NDist)
    CALL New(Tmp%Co,  Tmp%NCoef)
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Tmp%NQ%I  (i),i=1,Tmp%NExpt)
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Tmp%OffQ%I(i),i=1,Tmp%NExpt)
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Tmp%OffR%I(i),i=1,Tmp%NExpt)
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Tmp%Lndx%I(i),i=1,Tmp%NExpt)
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Tmp%Expt%D(i),i=1,Tmp%NExpt)
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Tmp%Qx%D  (i),i=1,Tmp%NDist)
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Tmp%Qy%D  (i),i=1,Tmp%NDist)
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Tmp%Qz%D  (i),i=1,Tmp%NDist)
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Tmp%Co%D  (i),i=1,Tmp%NCoef)
    Close(UNIT=Seq,STATUS='KEEP')

#ifdef DIAG_QCTC
IF(MyID == 0) THEN
   write(*,*) 'Tmp%NQ%I(Tmp%NExpt) =', Tmp%NQ%I(Tmp%NExpt)
   write(*,*) 'Tmp%Lndx%I(Tmp%NExpt) =', Tmp%Lndx%I(Tmp%NExpt)
   LMNLen = LHGTF(Tmp%Lndx%I(Tmp%NExpt))
   write(*,*) 'LMNLen =', LMNLen
ENDIF
#endif

    ENExpt = Tmp%NExpt-1
    ENDist = Sum(Tmp%NQ%I(1:ENExpt))
    NNDist = Tmp%NQ%I(Tmp%NExpt)
    NNCoef = LHGTF(Tmp%Lndx%I(Tmp%NExpt))*Tmp%NQ%I(Tmp%NExpt)
    ENCoef = Tmp%NCoef-NNCoef
    NDist = Sum(Tmp%NQ%I(1:Tmp%NExpt))
    IF(NDist /= Tmp%NDist) THEN
      STOP 'ERR: something is wrong!'
    ENDIF

    SumECoef = 0.0D0
    DO I = 1, ENCoef
      SumECoef = SumECoef + Tmp%Co%D(I)
    ENDDO

    TotENExpt= AllReduce(ENExpt)
    TotENDist = AllReduce(ENDist)
    TotENCoef = AllReduce(ENCoef)
    TotNNCoef = AllReduce(NNCoef)
    TotNDist = AllReduce(NDist)
    TotNCoef = AllReduce(Tmp%NCoef)
    TotNNDist = AllReduce(NNDist)

    TotSumECoef = AllReduce(SumECoef)
#ifdef DIAG_QCTC
IF(MyID ==0 ) THEN
   write(*,*) 'TotENExpt = ',TotENExpt
   write(*,*) 'TotENDist = ',TotENDist
   write(*,*) 'TotENCoef = ',TotENCoef
   write(*,*) 'TotNNCoef = ',TotNNCoef
   write(*,*) 'TotNCoef = ',TotNCoef
   write(*,*) 'TotSumECoef = ',TotSumECoef
   write(*,*) 'TotNDist = ',TotNDist
   write(*,*) 'TotNNDist = ',TotNNDist
ENDIF
#endif

    !! electronic + 1
    NewNExpt = TotENExpt + 1
    NewNDist = TotNDist
    NewNCoef = TotNCoef
#ifdef DIAG_QCTC
IF(MyID == 0) THEN
    write(*,*) 'NewNExpt,NewNDist,NewNCoef=',NewNExpt,NewNDist,NewNCoef
ENDIF
#endif
    CALL New_HGRho(Rho,(/NewNExpt,NewNDist,NewNCoef,Tmp%NSDen/))!<<<SPIN

    CALL New(CSize,NPrc-1,M_O=0)
    CALL New(disp,NPrc-1,M_O=0)
    CALL MPI_AllGather(ENExpt,1,MPI_INTEGER,CSize%I(0),1,MPI_INTEGER,MONDO_COMM,IErr)
#ifdef DIAG_QCTC
IF(MyID == 1) THEN
  sumi = 0
  DO i = 0, NPrc-1
    sumi  = sumi + Csize%I(i)
  ENDDO
  WRITE(*,*) 'sumi = ', sumi
ENDIF
#endif
    disp%I(0) = 0
    do i = 1, nprc-1
      disp%i(i) = disp%i(i-1)+csize%i(i-1)
    enddo
    CALL MPI_AllGatherV(Tmp%NQ%I(1),ENExpt,MPI_INTEGER,Rho%NQ%I(1),Csize%I(0),disp%I(0),MPI_INTEGER, &
                        MONDO_COMM,IErr)
    CALL MPI_AllGatherV(Tmp%Lndx%I(1),ENExpt,MPI_INTEGER,Rho%Lndx%I(1),Csize%I(0),disp%I(0),MPI_INTEGER, &
                        MONDO_COMM,IErr)
    CALL MPI_AllGatherV(Tmp%Expt%d(1),ENExpt,MPI_DOUBLE_PRECISION,Rho%Expt%d(1),Csize%I(0),disp%I(0), &
                        MPI_DOUBLE_PRECISION,MONDO_COMM,IErr)
    Rho%NQ%i(NewNExpt) = TotNNDist
    Rho%Lndx%i(NewNExpt) = Tmp%Lndx%i(Tmp%NExpt)
    Rho%Expt%D(NewNExpt) = Tmp%expt%d(Tmp%NExpt)
    if(Rho%Lndx%i(NewNExpt) /= 0) THEN
      stop 'must be zero!!!'
    endif
    if(Rho%expt%d(newnexpt) /= nuclearexpnt) then
      stop 'must be nuclearexpt!'
    endif

if(myid == 1) then
  do i = 1, totenexpt-1
    ! write(*,*) 'totrho%nq%i(i=',i,') = ',totrho%nq%i(i)
    ! write(*,*) 'totrho%lndx%i(i=',i,') = ',totrho%lndx%i(i)
    ! write(*,*) 'totrho%expt%d(i=',i,') = ',totrho%expt%d(i)
  enddo
endif

    CALL MPI_AllGather(ENDist,1,MPI_INTEGER,CSize%I(0),1,MPI_INTEGER,MONDO_COMM,IErr)
    disp%I(0) = 0
    do i = 1, nprc-1
      disp%i(i) = disp%i(i-1)+csize%i(i-1)
    enddo
    call mpi_allgatherv(Tmp%Qx%d(1),ENDist,MPI_DOUBLE_PRECISION,Rho%Qx%d(1), &
                        CSize%I(0),disp%i(0),MPI_DOUBLE_PRECISION,MONDO_COMM,IErr)
    call mpi_allgatherv(Tmp%Qy%d(1),ENDist,MPI_DOUBLE_PRECISION,Rho%Qy%d(1), &
                        CSize%I(0),disp%i(0),MPI_DOUBLE_PRECISION,MONDO_COMM,IErr)
    call mpi_allgatherv(Tmp%Qz%d(1),ENDist,MPI_DOUBLE_PRECISION,Rho%Qz%d(1), &
                       CSize%I(0),disp%i(0),MPI_DOUBLE_PRECISION,MONDO_COMM,IErr)

    CALL MPI_AllGather(NNDist,1,MPI_INTEGER,CSize%I(0),1,MPI_INTEGER,MONDO_COMM,IErr)

    disp%I(0) = 0
    do i = 1, nprc-1
      disp%i(i) = disp%i(i-1)+csize%i(i-1)
    enddo
    call mpi_allgatherv(tmp%Qx%d(ENDist+1),NNDist,MPI_DOUBLE_PRECISION,Rho%Qx%d(TotENDist+1), &
         & CSize%I(0),disp%i(0),MPI_DOUBLE_PRECISION,MONDO_COMM,IErr)
    call mpi_allgatherv(tmp%Qy%d(ENDist+1),NNDist,MPI_DOUBLE_PRECISION,Rho%Qy%d(TotENDist+1), &
         & CSize%I(0),disp%i(0),MPI_DOUBLE_PRECISION,MONDO_COMM,IErr)
    call mpi_allgatherv(tmp%Qz%d(ENDist+1),NNDist,MPI_DOUBLE_PRECISION,Rho%Qz%d(TotENDist+1), &
         & CSize%I(0),disp%i(0),MPI_DOUBLE_PRECISION,MONDO_COMM,IErr)


    ! now Co
    CALL MPI_AllGather(ENCoef,1,MPI_INTEGER,CSize%I(0),1,MPI_INTEGER,MONDO_COMM,IErr)
    disp%I(0) = 0
    do i = 1, nprc-1
      disp%i(i) = disp%i(i-1)+csize%i(i-1)
    enddo
    call mpi_allgatherv(Tmp%Co%d(1),ENCoef,MPI_DOUBLE_PRECISION,Rho%Co%d(1), &
         & CSize%I(0),disp%i(0),MPI_DOUBLE_PRECISION,MONDO_COMM,IErr)
    CALL MPI_AllGather(NNCoef,1,MPI_INTEGER,CSize%I(0),1,MPI_INTEGER,MONDO_COMM,IErr)
    disp%I(0) = 0
    do i = 1, nprc-1
      disp%i(i) = disp%i(i-1)+csize%i(i-1)
    enddo
    call mpi_allgatherv(Tmp%Co%d(ENCoef+1),NNCoef,MPI_DOUBLE_PRECISION, &
         & Rho%Co%d(TotENCoef+1),CSize%I(0),disp%i(0),MPI_DOUBLE_PRECISION,MONDO_COMM,IErr)

    ! finally OffR and OffQ
    Rho%OffQ%i(1) = 0
    Rho%OffR%i(1) = 0
    do i = 2, NewNExpt
      Rho%OffQ%i(i) = Rho%OffQ%i(i-1)+Rho%NQ%i(i-1)
      Rho%OffR%i(i) = Rho%OffR%i(i-1)+Rho%NQ%i(i-1)*LHGTF(Rho%Lndx%I(i-1))
    ! write(*,*) 'i = ', i, ', Rho%OffQ%i(i) = ', Rho%OffQ%i(i), ', Rho%OffR%i(i) = ', Rho%OffR%i(i)
    enddo
    RETURN
 202 CALL Halt('Died in GetRho, IOSTAT='//TRIM(IntToChar(IOS)))
  END SUBROUTINE GetDistrRho
!------------------------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------------------------
  SUBROUTINE NukE_ENPart(GM_Loc)
    TYPE(CRDS),  INTENT(IN)  :: GM_Loc
    INTEGER :: NProc,NAtoms
    CALL New(BegAtInd,NPrc-1,0)
    CALL New(EndAtInd,NPrc-1,0)
    NProc = NPrc
    NAtoms = GM_Loc%NAtms
    CALL ENPart(NProc,NAtoms,BegAtInd,EndAtInd)
  END SUBROUTINE NukE_ENPart
!------------------------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------------------------
  SUBROUTINE ParaRhoToPoleTree
    INTEGER :: IErr
    ALLOCATE(PA1(0:NPrc-1),PA2(0:NPrc-1),STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In ParaRhoToPoleTree: Allocation problem.')
    CALL MakeTree1
    CALL BuildLocalTrees
    CALL ConcatTrees
    DO CurrentTier = BranchTier-1, 0, -1
      CALL MakePoleTree(PR1)
    ENDDO
    PR1%Ell = SPEll
    PoleRoot => PR1
!
!!$    IF(MyID==0) THEN
!!$       WRITE(*,*) 'BranchTier = ',BranchTier
!!$       CALL CheckNodes(PoleRoot)
!!$    ENDIF
!!$    IF(.TRUE.) STOP
!
  END SUBROUTINE ParaRhoToPoleTree
!------------------------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------------------------
  SUBROUTINE ConcatTrees
    RunPointer = -1
    CALL SetLink(PR1)
  END SUBROUTINE ConcatTrees
!------------------------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE SetLink(P)
    TYPE(PoleNode),POINTER :: P,Left,Right
    IF(P%Box%Tier == BranchTier-1) THEN
      DEALLOCATE(P%Descend%Travrse)
      DEALLOCATE(P%Descend)
      RunPointer = RunPointer + 1
      P%Descend => PA2(RunPointer)%Ptr

      RunPointer = RunPointer + 1
      P%Descend%Travrse => PA2(RunPointer)%Ptr
    ELSE
      Left => P%Descend
      Right => P%Descend%Travrse
      CALL SetLink(Left)
      CALL SetLink(Right)
    ENDIF
  END SUBROUTINE SetLink
!------------------------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------------------------
  SUBROUTINE BuildLocalTrees
    TYPE(PoleNode),POINTER :: TP
    INTEGER :: I,J,IErr,SendTo,RecvFr,ActIntRec,ActDblRec,IntSize,DblSize
    INTEGER,ALLOCATABLE :: SF(:),Dest(:),NodesNumArr(:),IntNumArr(:),DblNumArr(:)
    INTEGER,DIMENSION(MPI_STATUS_SIZE) :: IntStatus
    INTEGER,DIMENSION(MPI_STATUS_SIZE) :: DblStatus
!
    MaxTier = BranchTier
    CALL NewPoleNode(PA2(MyID)%Ptr,BranchTier)
    TP => PA2(MyID)%Ptr
    TP%Bdex = PA1(MyID)%Ptr%Bdex
    TP%Edex = PA1(MyID)%Ptr%Edex
    TP%NQ   = TP%Edex - TP%Bdex + 1
    CALL SplitPole(TP)
    DO CurrentTier = MaxTier,BranchTier, -1
      CALL MakePoleTree(TP)
    ENDDO
!
    IntNum = 0
    DblNum = 0
    NodesVisit = 0
    CALL FindSize(TP)
!
    ALLOCATE(IntArr(IntNum),STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In ParaQ 100: Allocation problem.')
    ALLOCATE(DblArr(DblNum),STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In ParaQ 200: Allocation problem.')
!
    CALL PackTree(TP)
!
    IF(IntNum /= IntIndex) THEN
      STOP 'ERR: Missing integer during packing!'
    ENDIF
    IF(DblNum /= DblIndex) THEN
      STOP 'ERR: Missing doubles during packing!'
    ENDIF
    IF(NodesVisit /= NodesVisit1) THEN
      STOP 'ERR: Missing nodes ?'
    ENDIF
!
    ALLOCATE(NodesNumArr(0:NPrc-1),IntNumArr(0:NPrc-1),DblNumArr(0:NPrc-1),STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In ParaQ 400: Allocation problem.')
!
    CALL MPI_AllGather(NodesVisit,1,MPI_INTEGER,NodesNumArr(0),1,MPI_INTEGER,MONDO_COMM,IErr)
    CALL MPI_AllGather(IntNum,1,MPI_INTEGER,IntNumArr(0),1,MPI_INTEGER,MONDO_COMM,IErr)
    CALL MPI_AllGather(DblNum,1,MPI_INTEGER,DblNumArr(0),1,MPI_INTEGER,MONDO_COMM,IErr)
    CALL MPI_AllReduce(IntNum,GIntNum,1,MPI_INTEGER,MPI_MAX,MONDO_COMM,IErr)
    CALL MPI_AllReduce(DblNum,GDblNum,1,MPI_INTEGER,MPI_MAX,MONDO_COMM,IErr)
!
#ifdef DIAG_QCTC
    IF(MyID == 0) THEN
       CALL OpenASCII(OutFile,Out)
       WRITE(*,'(A,F10.5)') 'int space(MB) for tree pack = ',GIntNum*4.0D-6
       WRITE(*,'(A,F10.5)') 'dbl space(MB) for tree pack = ',GDblNum*8.0D-6
       WRITE(Out,'(A,F10.5)') 'int space(MB) for tree pack = ',GIntNum*4.0D-6
       WRITE(Out,'(A,F10.5)') 'dbl space(MB) for tree pack = ',GDblNum*8.0D-6
       CLOSE(Out,STATUS='KEEP')
    ENDIF
#endif
!
!!$    IF(MyID == 1) THEN
!!$       WRITE(*,*) 'NodesNumArr= ',NodesNumArr(:)
!!$       WRITE(*,*) 'IntNumArr  = ',IntNumArr(:)
!!$       WRITE(*,*) 'DblNumArr  = ',DblNumArr(:)
!!$    ENDIF
!
    ALLOCATE(RecIntArr(GIntNum),STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In ParaQ 700: Allocation problem.')
    ALLOCATE(RecDblArr(GDblNum),STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In ParaQ 800: Allocation problem.')

    !VWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVW
    DO I=NPrc-1,0,-1
       IntSize=IntNumArr(I)
       DblSize=DblNumArr(I)
       IF(MyID .EQ. I) THEN
!         Copy arrays
          CALL INT_VECT_EQ_INT_VECT(IntSize,RecIntArr(1),IntArr(1))
          CALL DBL_VECT_EQ_DBL_VECT(DblSize,RecDblArr(1),DblArr(1))
       ENDIF
!      BCast the arrays
       CALL MPI_BCAST(RecIntArr(1),IntSize,MPI_INTEGER         ,I,MONDO_COMM,IErr)
       CALL MPI_BCAST(RecDblArr(1),DblSize,MPI_DOUBLE_PRECISION,I,MONDO_COMM,IErr)
!      Build tree
       IF(MyID .NE. I) THEN
          CALL CopyTree(I)
          IF(IntIndex /= IntNumArr(I)) STOP 'ERR: should be the same!'
          IF(DblIndex /= DblNumArr(I)) STOP 'ERR: should be the same!'
          IF(NodesVisit /= NodesNumArr(I)) STOP 'ERR: should be the same!'
       ENDIF
    ENDDO
    !VWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVWVW
!!$    DO I = 1, NPrc-1
!!$      DO J = 0, NPrc-1
!!$        SendTo = MODULO(J+I,NPrc)
!!$        Dest(J) = SendTo
!!$        SF(J) = 0
!!$      ENDDO
!!$      DO J = 0, NPrc-1
!!$        IF(SF(J) == 0 .AND. SF(Dest(J)) == 0) THEN
!!$          SF(J) = 1
!!$          SF(Dest(J)) = 2
!!$        ENDIF
!!$      ENDDO
!!$      SendTo = MODULO(MyID+I,NPrc)
!!$      RecvFr = MODULO(MyID-I,NPrc)
!!$
!!$      IF(SF(MyID) == 1) THEN
!!$         CALL MPI_Send(IntArr(1),IntNum,MPI_INTEGER,SendTo         ,MyID,MONDO_COMM,IErr)
!!$         CALL MPI_Send(DblArr(1),DblNum,MPI_DOUBLE_PRECISION,SendTo,MyID,MONDO_COMM,IErr)
!!$
!!$         CALL MPI_Recv(RecIntArr(1),IntNumArr(RecvFr),MPI_INTEGER, &
!!$              & RecvFr,RecvFr,MONDO_COMM,IntStatus,IErr)
!!$         CALL MPI_Recv(RecDblArr(1),DblNumArr(RecvFr),MPI_DOUBLE_PRECISION, &
!!$              & RecvFr,RecvFr,MONDO_COMM,DblStatus,IErr)
!!$!orig         CALL MPI_Recv(RecIntArr(1),GIntNum,MPI_INTEGER,RecvFr         ,RecvFr, &
!!$!orig              & MONDO_COMM,IntStatus,IErr)
!!$!orig         CALL MPI_Recv(RecDblArr(1),GDblNum,MPI_DOUBLE_PRECISION,RecvFr,RecvFr, &
!!$!orig              & MONDO_COMM,DblStatus,IErr)
!!$
!!$        CALL MPI_Get_Count(IntStatus,MPI_INTEGER,ActIntRec,IErr)
!!$        CALL MPI_Get_Count(DblStatus,MPI_DOUBLE_PRECISION,ActDblRec,IErr)
!!$        IF(ActIntRec /= IntNumArr(RecvFr)) STOP 'ERR: 1: Missing int while receiving ?'
!!$        IF(ActDblRec /= DblNumArr(RecvFr)) STOP 'ERR: 1: Missing dbl while receiving ?'
!!$
!!$        CALL CopyTree(RecvFr)
!!$
!!$        IF(IntIndex /= IntNumArr(RecvFr)) STOP 'ERR: should be the same!'
!!$        IF(DblIndex /= DblNumArr(RecvFr)) STOP 'ERR: should be the same!'
!!$        IF(NodesVisit /= NodesNumArr(RecvFr)) STOP 'ERR: should be the same!'
!!$
!!$      ELSE
!!$
!!$         CALL MPI_Recv(RecIntArr(1),IntNumArr(RecvFr),MPI_INTEGER, &
!!$              & RecvFr,RecvFr,MONDO_COMM,IntStatus,IErr)
!!$         CALL MPI_Recv(RecDblArr(1),DblNumArr(RecvFr),MPI_DOUBLE_PRECISION, &
!!$              & RecvFr,RecvFr,MONDO_COMM,DblStatus,IErr)
!!$!orig         CALL MPI_Recv(RecIntArr(1),GIntNum,MPI_INTEGER,RecvFr         ,RecvFr, &
!!$!orig              & MONDO_COMM,IntStatus,IErr)
!!$!orig         CALL MPI_Recv(RecDblArr(1),GDblNum,MPI_DOUBLE_PRECISION,RecvFr,RecvFr, &
!!$!orig              & MONDO_COMM,DblStatus,IErr)
!!$
!!$        CALL MPI_Get_Count(IntStatus,MPI_INTEGER,ActIntRec,IErr)
!!$        CALL MPI_Get_Count(DblStatus,MPI_DOUBLE_PRECISION,ActDblRec,IErr)
!!$        IF(ActIntRec /= IntNumArr(RecvFr)) STOP 'ERR: 1: Missing int while receiving ?'
!!$        IF(ActDblRec /= DblNumArr(RecvFr)) STOP 'ERR: 1: Missing dbl while receiving ?'
!!$
!!$        CALL CopyTree(RecvFr)
!!$
!!$        IF(IntIndex /= IntNumArr(RecvFr)) STOP 'ERR: should be the same!'
!!$        IF(DblIndex /= DblNumArr(RecvFr)) STOP 'ERR: should be the same!'
!!$        IF(NodesVisit /= NodesNumArr(RecvFr)) STOP 'ERR: should be the same!'
!!$
!!$        CALL MPI_Send(IntArr(1),IntNum,MPI_INTEGER,SendTo         ,MyID,MONDO_COMM,IErr)
!!$        CALL MPI_Send(DblArr(1),DblNum,MPI_DOUBLE_PRECISION,SendTo,MyID,MONDO_COMM,IErr)
!!$
!!$     ENDIF
!!$  ENDDO
!!
    DEALLOCATE(RecIntArr,STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In ParaQ 900: DeAllocation problem.')
    DEALLOCATE(RecDblArr)
    IF(iErr.NE.0) CALL Halt('In ParaQ 1000: DeAllocation problem.')
    DEALLOCATE(IntArr)
    IF(iErr.NE.0) CALL Halt('In ParaQ 1100: DeAllocation problem.')
    DEALLOCATE(DblArr)
    IF(iErr.NE.0) CALL Halt('In ParaQ 1200: DeAllocation problem.')
!
  END SUBROUTINE BuildLocalTrees
!------------------------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------------------------
  SUBROUTINE CopyTree(I)
    INTEGER :: I
!
    NodesVisit = 0
    IntIndex = 0
    DblIndex = 0
    CALL RecurCopy(PA2(I)%Ptr)
!
  END SUBROUTINE CopyTree
!------------------------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE RecurCopy(P)
    TYPE(PoleNode),POINTER :: P
    INTEGER :: SizeofSP,SizeofCo
    integer :: ierr
!
    NodesVisit = NodesVisit + 1
!   Allocate the Node
    CALL NewPoleNode(P,0)
!   Unpack the Logical Leaf
    IntIndex = IntIndex+1
    IF(RecIntArr(IntIndex) == 100)  THEN
       P%Leaf=.TRUE.
    ELSEIF(RecIntArr(IntIndex) == 200) THEN
       P%Leaf=.FALSE.
    ELSE
       CALL Halt("Error: Integer Array not unpacked")
    ENDIF
!   Unpack the Integers
    IntIndex = IntIndex+1;P%Bdex  = RecIntArr(IntIndex)
    IntIndex = IntIndex+1;P%Edex  = RecIntArr(IntIndex)
    IntIndex = IntIndex+1;P%NQ    = RecIntArr(IntIndex)
    IntIndex = IntIndex+1;P%Ell   = RecIntArr(IntIndex)
    IntIndex = IntIndex+1;P%EllCD = RecIntArr(IntIndex)
!   Unpack the Doubles
    DblIndex = DblIndex+1;P%Zeta     = RecDblArr(DblIndex)
    DblIndex = DblIndex+1;P%Strength = RecDblArr(DblIndex)
    DblIndex = DblIndex+1;P%DMax2    = RecDblArr(DblIndex)
    DblIndex = DblIndex+1;P%WCoef    = RecDblArr(DblIndex)
!   Unpack the Box
    IntIndex = IntIndex+1;P%Box%Tier   = RecIntArr(IntIndex)
    IntIndex = IntIndex+1;P%Box%Number = RecIntArr(IntIndex)
    DblIndex = DblIndex+3;P%Box%Center(1:3)   = RecDblArr(DblIndex-2:DblIndex)
    DblIndex = DblIndex+3;P%Box%Half(1:3)     = RecDblArr(DblIndex-2:DblIndex)
    DblIndex = DblIndex+3;P%Box%BndBox(1:3,1) = RecDblArr(DblIndex-2:DblIndex)
    DblIndex = DblIndex+3;P%Box%BndBox(1:3,2) = RecDblArr(DblIndex-2:DblIndex)
!   Unpack and Allocate the SP's
    SizeofSP = LSP(P%Ell)+1
!
    ALLOCATE(P%S(0:SizeofSP-1),STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In RecurCopy 200: Allocation problem.')
    ALLOCATE(P%C(0:SizeofSP-1),STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In RecurCopy 300: Allocation problem.')
!
    DblIndex = DblIndex+SizeofSP;P%C(0:SizeofSP-1) = RecDblArr(DblIndex-SizeofSP+1:DblIndex)
    DblIndex = DblIndex+SizeofSP;P%S(0:SizeofSP-1) = RecDblArr(DblIndex-SizeofSP+1:DblIndex)
!   Unpack and Allocate the HG's
    IF(P%Leaf) THEN
       SizeofCo = LHGTF(P%Ell)
!
       ALLOCATE(P%Co(1:SizeofCo),STAT=iErr)
       IF(iErr.NE.0) CALL Halt('In RecurCopy 400: Allocation problem.')
!
       DblIndex = DblIndex+SizeofCo;P%Co(1:SizeofCo) = RecDblArr(DblIndex-SizeofCo+1:DblIndex)
       NULLIFY(P%Descend)
       NULLIFY(P%Travrse)
    ELSE
      CALL RecurCopy(P%Descend)
      CALL RecurCopy(P%Descend%Travrse)
    ENDIF
  END SUBROUTINE RecurCopy
!------------------------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------------------------
  SUBROUTINE PackTree(P)
    TYPE(PoleNode),POINTER :: P
    NodesVisit1 = 0
    IntIndex = 0
    DblIndex = 0
    CALL RecurPack(P)
  END SUBROUTINE PackTree
!------------------------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE RecurPack(P)
    TYPE(PoleNode),POINTER :: P,Left,Right
    INTEGER :: SizeofCo,SizeofSP
!
    NodesVisit1 = NodesVisit1 + 1
!    WRITE(*,*) NodesVisit1,IntIndex,DblIndex
!   Store the logical Leaf
    IntIndex = IntIndex+1
    IF(      P%Leaf) IntArr(IntIndex) = 100
    IF(.NOT. P%Leaf) IntArr(IntIndex) = 200
!   Store the Integers
    IntIndex = IntIndex+1;IntArr(IntIndex) = P%Bdex
    IntIndex = IntIndex+1;IntArr(IntIndex) = P%Edex
    IntIndex = IntIndex+1;IntArr(IntIndex) = P%NQ
    IntIndex = IntIndex+1;IntArr(IntIndex) = P%Ell
    IntIndex = IntIndex+1;IntArr(IntIndex) = P%EllCD
!   Store the Doubles
    DblIndex = DblIndex+1;DblArr(DblIndex) = P%Zeta
    DblIndex = DblIndex+1;DblArr(DblIndex) = P%Strength
    DblIndex = DblIndex+1;DblArr(DblIndex) = P%DMax2
    DblIndex = DblIndex+1;DblArr(DblIndex) = P%WCoef
!   Store the Box
    IntIndex = IntIndex+1;IntArr(IntIndex) = P%Box%Tier
    IntIndex = IntIndex+1;IntArr(IntIndex) = P%Box%Number
    DblIndex = DblIndex+3;DblArr(DblIndex-2:DblIndex) = P%Box%Center(1:3)
    DblIndex = DblIndex+3;DblArr(DblIndex-2:DblIndex) = P%Box%Half(1:3)
    DblIndex = DblIndex+3;DblArr(DblIndex-2:DblIndex) = P%Box%BndBox(1:3,1)
    DblIndex = DblIndex+3;DblArr(DblIndex-2:DblIndex) = P%Box%BndBox(1:3,2)
!   Store the SP's
    SizeofSP = LSP(P%Ell)+1
    DblIndex = DblIndex+SizeofSP;DblArr(DblIndex-SizeofSP+1:DblIndex) = P%C(0:SizeofSP-1)
    DblIndex = DblIndex+SizeofSP;DblArr(DblIndex-SizeofSP+1:DblIndex) = P%S(0:SizeofSP-1)
!   Store the HG's
    IF(P%Leaf) THEN
       SizeofCo = LHGTF(P%Ell)
       DblIndex = DblIndex+SizeofCo;DblArr(DblIndex-SizeofCo+1:DblIndex) = P%Co(1:SizeofCo)
    ELSE
!     Recur Down
      Left => P%Descend
      Right => P%Descend%Travrse
      CALL RecurPack(Left)
      CALL RecurPack(Right)
    ENDIF

  END SUBROUTINE RecurPack
!------------------------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE FindSize(P)
    TYPE(PoleNode),POINTER :: P,Left,Right
    INTEGER :: SizeofCo,SizeofSP
!
    NodesVisit = NodesVisit + 1
    IntNum   = IntNum + 8
    DblNum   = DblNum + 16
    SizeofSP = LSP(P%Ell)+1
    DblNum   = DblNum + SizeofSP + SizeofSP
    IF(P%Leaf) THEN
       SizeofCo = LHGTF(P%Ell)
       DblNum   = DblNum + SizeofCo
    ELSE
!     Recur Down
      Left => P%Descend
      Right => P%Descend%Travrse
      CALL FindSize(Left)
      CALL FindSize(Right)
    ENDIF
  END SUBROUTINE FindSize
!------------------------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------------------------
  SUBROUTINE MakeTree1
    PoleNodes = 0
    MaxTier = 0
    RunPointer = -1
    CALL NewPoleNode(PR1,0)
    PR1%Bdex = 1
    PR1%Edex = Rho%NDist
    PR1%NQ = Rho%NDist
    BranchTier = NINT(LOG(NPrc*1.0D0)/LOG(2.0D0))
    IF(2.0D0**BranchTier /= NPrc) THEN
      STOP 'ERR: Check NPrc!'
    ENDIF
    CALL SplitPole1(PR1)
  END SUBROUTINE MakeTree1
!------------------------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE SplitPole1(P)
    TYPE(PoleNode),POINTER :: P, Left,Right
    IF(P%NQ == 1) THEN
      STOP  'Err: Nq must be greater than 1!'
    ELSE
      IF(P%Box%Tier == BranchTier) THEN
        RunPointer = RunPointer + 1
        PA1(RunPointer)%Ptr => P
      ELSE
        CALL NewPoleNode(P%Descend,P%Box%Tier+1)
        CALL NewPoleNode(P%Descend%Travrse,P%Box%Tier+1)
        Left => P%Descend
        Right => P%Descend%Travrse
        CALL SplitPoleBox(P,Left,Right)
        CALL SplitPole1(Left)
        CALL SplitPole1(Right)
        CALL NewSPArrays(P)
      ENDIF
    ENDIF
  END SUBROUTINE SplitPole1
!------------------------------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------------------------------
  SUBROUTINE CheckInBox(GMLoc)
    REAL(DOUBLE) :: Dist
    INTEGER :: I,NQ,iq,iadd,zq,TotNQ
    REAL(DOUBLE),DIMENSION(3) :: PQ
    TYPE(CRDS) :: GMLoc
    LOGICAL :: Fail

    RootBDist = ZERO
    TotNQ = 0
    DO zq = 1, Rho%NExpt
      NQ = Rho%NQ%I(zq)
      DO iq = 1, NQ
        TotNQ = TotNQ + 1
        iadd = Rho%OffQ%I(zq) + iq
        PQ(1) = Rho%Qx%D(iadd)-GMLoc%PBC%CellCenter%D(1)
        PQ(2) = Rho%Qy%D(iadd)-GMLoc%PBC%CellCenter%D(2)
        PQ(3) = Rho%Qz%D(iadd)-GMLoc%PBC%CellCenter%D(3)

        Dist = SQRT(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
        IF(Dist > RootBDist) THEN
          RootBDist = Dist
        ENDIF
      ENDDO

    ENDDO
#ifdef DIAG_QCTC
    WRITE(*,*) 'EN partition: TotNQ = ',TotNQ, ', RootBDist = ',RootBDist
#endif
    CALL ENPart(NPrc,TotNQ,BegNQInd,EndNQInd)
  END SUBROUTINE CheckInBox
!------------------------------------------------------------------------------------------------------
! count number of primitives
!------------------------------------------------------------------------------------------------------
  SUBROUTINE EBBCountJBlock(Count,Pair,PoleRoot)
    TYPE(AtomPair)  :: Pair
    TYPE(PoleNode), POINTER  :: PoleRoot
    INTEGER :: KA,KB,CFA,CFB,PFA,PFB,Count,J
    TYPE(BBox) :: NodeBox
    REAL(DOUBLE) :: ZA,ZB

    KA = Pair%KA
    KB = Pair%KB
    Prim%AB2=Pair%AB2
    Prim%A = Pair%A
    Prim%B = Pair%B
    Count = 0
    DO CFA = 1, BS%NCFnc%I(KA)
      DO CFB = 1, BS%NCFnc%I(KB)
        DO PFA = 1, BS%NPFnc%I(CFA,KA)
          DO PFB = 1, BS%NPFnc%I(CFB,KB)
            Prim%ZA = BS%Expnt%D(PFA,CFA,KA)
            Prim%ZB = BS%Expnt%D(PFB,CFB,KB)
            Prim%Zeta = Prim%ZA + Prim%ZB
            Prim%Xi = Prim%ZA*Prim%ZB/Prim%Zeta
            IF(TestPrimPair(Prim%Xi,Prim%AB2)) THEN
              ZA = Prim%ZA
              ZB = Prim%ZB
              Prim%P=(ZA*Prim%A+ZB*Prim%B)/Prim%Zeta
              NodeBox%BndBox(1:3,1) = Prim%P(1:3)
              NodeBox%BndBox(1:3,2) = Prim%P(1:3)
              NodeBox = ExpandBox(NodeBox,EP)
              DO J = 1, 3
                RootBox%BndBox(J,1) = MIN(RootBox%BndBox(J,1),NodeBox%BndBox(J,1))
                RootBox%BndBox(J,2) = MAX(RootBox%BndBox(J,2),NodeBox%BndBox(J,2))
              ENDDO
              Count = Count + 1
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE EBBCountJBlock



  ! count number of primitives in the box
  SUBROUTINE BoxEBBPairCount(Count,Pair,PoleRoot)
    TYPE(AtomPair)  :: Pair
    TYPE(PoleNode), POINTER  :: PoleRoot
    INTEGER :: KA,KB,CFA,CFB,PFA,PFB,Count,J
    TYPE(BBox) :: NodeBox
    REAL(DOUBLE) :: ZA,ZB

    KA = Pair%KA
    KB = Pair%KB
    Prim%AB2=Pair%AB2
    Prim%A = Pair%A
    Prim%B = Pair%B
    Count = 0
    DO CFA = 1, BS%NCFnc%I(KA)
      DO CFB = 1, BS%NCFnc%I(KB)
        DO PFA = 1, BS%NPFnc%I(CFA,KA)
          DO PFB = 1, BS%NPFnc%I(CFB,KB)
            Prim%ZA = BS%Expnt%D(PFA,CFA,KA)
            Prim%ZB = BS%Expnt%D(PFB,CFB,KB)
            Prim%Zeta = Prim%ZA + Prim%ZB
            Prim%Xi = Prim%ZA*Prim%ZB/Prim%Zeta
            IF(TestPrimPair(Prim%Xi,Prim%AB2)) THEN
              ZA = Prim%ZA
              ZB = Prim%ZB
              Prim%P=(ZA*Prim%A+ZB*Prim%B)/Prim%Zeta
              IF(Prim%P(1) >= LC%D(1,MyID+1) .AND. Prim%P(1) < RC%D(1,MyID+1) &
                 .AND. Prim%P(2) >= LC%D(2,MyID+1) .AND. Prim%P(2) < RC%D(2,MyID+1) &
                 .AND. Prim%P(3) >= LC%D(3,MyID+1) .AND. Prim%P(3) < RC%D(3,MyID+1)) THEN
                Count = Count + 1
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE BoxEBBPairCount


  ! equaltime partition
  SUBROUTINE ET_Part
    INTEGER :: NIter,I,J,IErr,NV,NS,RI,PI,CI,RootNum,Dir
    TYPE(BBox) :: NodeBox,LBox
    INTEGER :: P2(0:30)
    TYPE(DBL_RNK2) :: RepLC,RepRC
    TYPE(DBL_VECT) :: VTm,Tau
    REAL(DOUBLE),PARAMETER :: TauB=1.0D-5,TauE=1.0D-4
    REAL(DOUBLE) :: NewCost,TVol,OVol,NVol,Diff,TotCost,T1,T2,T3,LD,MD,PCost,x0,x1,x2,&
                    Imbalance,MaxCost,WasteCost,f0,f1,f2,LocalT,TotLT
    TYPE(INT_VECT) :: XYZ
    TYPE(DBL_VECT) :: Pos

    DO I = 1, TotPrCount
      NodeBox%BndBox(1:3,1) = (/PosTimePair%D(1,I),PosTimePair%D(2,I),PosTimePair%D(3,I)/)
      NodeBox%BndBox(1:3,2) = (/PosTimePair%D(1,I),PosTimePair%D(2,I),PosTimePair%D(3,I)/)
      NodeBox = ExpandBox(NodeBox,EP)
      IF(I == 1) THEN
        LBox%BndBox(1:3,1:2) = NodeBox%BndBox(1:3,1:2)
      ELSE
        DO J = 1, 3
          LBox%BndBox(J,1) = MIN(LBox%BndBox(J,1),NodeBox%BndBox(J,1))
          LBox%BndBox(J,2) = MAX(LBox%BndBox(J,2),NodeBox%BndBox(J,2))
        ENDDO
      ENDIF
    ENDDO
    CALL MPI_AllReduce(LBox%BndBox(1,1),EBB_GlobalBox%BndBox(1,1),&
         3,MPI_DOUBLE_PRECISION,MPI_MIN,MONDO_COMM,IErr)
    CALL MPI_AllReduce(LBox%BndBox(1,2),EBB_GlobalBox%BndBox(1,2),&
         3,MPI_DOUBLE_PRECISION,MPI_MAX,MONDO_COMM,IErr)
    IF(MyID == 0) THEN
      DO I = 1, 3
        Diff = ABS(EBB_GlobalBox%BndBox(I,1)-RootBox%BndBox(I,1))
        IF(Diff > 1.0D-8) THEN
          WRITE(*,*) 'I = ',I, ', Diff = ',Diff
          WRITE(*,*) 'EBB_GlobalBox%BndBox(I,1) = ',EBB_GlobalBox%BndBox(I,1)
          WRITE(*,*) 'RootBox%BndBox(I,1) = ',RootBox%BndBox(I,1)
          STOP 'ERR: Boxes not consistent!'
        ENDIF
      ENDDO
    ENDIF

    P2(0) = 1
    DO I = 1, 30
      P2(I) = P2(I-1)*2
    ENDDO

    NV = NPrc
    NS = NINT(LOG(NV*1.0D0)/LOG(2.0D0))
    IF(P2(NS) /= NV) THEN
      STOP 'ERR: In ET_Part, P2 problem!'
    ENDIF

    CALL New(RepLC,(/3,NV/))
    CALL New(RepRC,(/3,NV/))
    CALL New(VTm,NV)


    TotCost = Reduce(ETimer(2))
    IF(MyID == 0) THEN
      RepLC%D(1:3,1) = EBB_GlobalBox%BndBox(1:3,1)
      RepRC%D(1:3,1) = EBB_GlobalBox%BndBox(1:3,2)
      VTm%D(1) = TotCost
      CALL New(Tau,NS)
      IF(NS == 1) THEN
        Tau%D(1) = TauB
      ELSE
        T1 = LOG(TauB)
        T2 = LOG(TauE)
        DO I = 1, NS
          T3 = T1 + (I-1.0D0)*(T2-T1)/(NS-1.0D0)
          Tau%D(I) = EXP(T3)
        ENDDO
      ENDIF
      RootNum = NV-1
      CALL New(XYZ,RootNum)
      CALL New(Pos,RootNum)
      RI = 0
      DO I = 1, NS
        DO PI = 1, P2(I-1)
          CI = PI + P2(I-1)
          MD = 0.0D0
          DO J = 1, 3
            LD = RepRC%D(J,PI)-RepLC%D(J,PI)
            IF(LD <= 0) THEN
              STOP 'ERR: LD not positive!'
            ENDIF
            IF(LD > MD) THEN
              MD = LD
              Dir = J
            ENDIF
          ENDDO
          RI = RI + 1
          XYZ%I(RI) = Dir
          RepLC%D(1:3,CI) = RepLC%D(1:3,PI)
          RepRC%D(1:3,CI) = RepRC%D(1:3,PI)
          PCost = VTm%D(PI)
          x0 = RepLC%D(Dir,PI)
          x1 = RepRC%D(Dir,PI)
          f0 = -0.5D0
          f1 = 0.5D0
          NIter = 0
          BP(1:3) = RepLC%D(1:3,PI)
          BP(4:6) = RepRC%D(1:3,PI)
          DO
            NIter = NIter + 1
            IF(NIter > 100) THEN
              STOP 'ERR: Too many iterations in ParaQ'
            ENDIF
            x2 = (x0+x1)*0.50D0
            RepRC%D(Dir,PI) = x2
            BP(3+Dir) = x2
            CALL MPI_BCast(BP(1),6,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
            LocalT = CalTime()
            TotLT = Reduce(LocalT)
            f2 = TotLT*1.0D0/(PCost*1.0D0)-0.50D0
            IF(ABS(f2) < Tau%D(I)) THEN
              RepLC%D(Dir,CI) = x2
              VTm%D(PI) = TotLT
              VTm%D(CI) = PCost-TotLT
              POS%D(RI) = x2
              EXIT
            ELSE
              IF(f2*f0 < 0) THEN
                x1 = x2; f1 = f2
              ELSE
                x0 = x2; f0 = f2
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      BP(:) = 0.0
      CALL MPI_BCast(BP(1),6,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)

      TVol = 1.0D0
      DO I = 1, 3
        LD = EBB_GlobalBox%BndBox(I,2)-EBB_GlobalBox%BndBox(I,1)
        IF(LD <= 0) THEN
          STOP 'ERR: LD Not positive!'
        ENDIF
        TVol = TVol*LD
      ENDDO
      OVol = TVol
      NVol = ZERO
      DO I = 1, NV
        TVol = 1.0D0
        DO J = 1, 3
          LD = RepRC%D(J,I)-RepLC%D(J,I)
          IF(LD <= 0) THEN
            STOP 'ERR: LD Not positive1!'
          ENDIF
          TVol = TVol*LD
        ENDDO
        NVol = NVol + TVol
      ENDDO
      Diff = ABS(NVol-OVol)
      IF(Diff > 1.0D-8) THEN
         WRITE(*,*) 'NVol = ',NVol
         WRITE(*,*) 'OVol = ',OVol
         STOP 'ERR: Diff in Vol!'
      ENDIF
      NewCost = ZERO
      DO I = 1, NV
        NewCost = NewCost + VTm%D(I)
      ENDDO
      Diff = ABS(NewCost - TotCost)
      IF(Diff > 1.0D-8) THEN
        WRITE(*,*) 'TotCost = ',TotCost
        WRITE(*,*) 'NewCost = ',NewCost
        STOP 'ERR: Diff in TotCost!'
      ENDIF
      MaxCost = -10000.00
      DO I = 1, NV
        IF(VTm%D(I) > MaxCost) THEN
          MaxCost = VTm%D(I)
        ENDIF
      ENDDO
      WasteCost = Zero
      DO I = 1, NV
        WasteCost = WasteCost + (MaxCost-VTm%D(I))
      ENDDO
      Imbalance = WasteCost/(NV*MaxCost)
#ifdef DIAG_QCTC
      write(*,*) 'Predicted Imbalance is ', Imbalance
#endif
      CALL Put(NPrc,'QLineLoc')
      CALL Put(XYZ,'QETDir')
      CALL Put(POS,'QETRoot')
      CALL Delete(XYZ)
      CALL Delete(Pos)
      CALL Delete(Tau)
    ELSE
      DO
        CALL MPI_BCast(BP(1),6,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
        IF(SUM(ABS(BP(:))) < 1.0D-14) THEN
          EXIT
        ENDIF
        LocalT = CalTime()
        TotLT = Reduce(LocalT)
      ENDDO

    ENDIF
    CALL Delete(RepLC)
    CALL Delete(RepRC)
    CALL Delete(VTm)

  END SUBROUTINE ET_Part

  FUNCTION CalTime()
    REAL(DOUBLE) :: CubeV,OverlapV,P(3),Q(3),Ovlap(3),CalTime,LE(1:3),UE(3)
    INTEGER :: I,KQ,J
    TYPE(BBox) :: NodeBox


    CalTime = 0.0D0
    DO J = 1, TotPrCount
      NodeBox%BndBox(1:3,1) = PosTimePair%D(1:3,J)
      NodeBox%BndBox(1:3,2) = PosTimePair%D(1:3,J)
      NodeBox = ExpandBox(NodeBox,EP)
      LE(1:3) = NodeBox%BndBox(1:3,1)
      UE(1:3) = NodeBox%BndBox(1:3,2)
      P(1:3) = BP(1:3)
      Q(1:3) = BP(4:6)
      IF(UE(1) <= P(1) .OR. LE(1) >= Q(1) .OR. &
       UE(2) <= P(2) .OR. LE(2) >= Q(2) .OR. &
       UE(3) <= P(3) .OR. LE(3) >= Q(3)) THEN
      ELSE
        DO I = 1, 3
          IF(LE(I) < P(I)) THEN
            IF(UE(I) <= Q(I)) THEN
              Ovlap(I) = UE(I) - P(I)
            ELSE
              Ovlap(I) = Q(I) - P(I)
            ENDIF
          ELSE
            IF(UE(I) <= Q(I)) THEN
               Ovlap(I) = UE(I) - LE(I)
            ELSE
               Ovlap(I) = Q(I) - LE(I)
            ENDIF
          ENDIF
        ENDDO
        OverlapV = Ovlap(1)*Ovlap(2)*Ovlap(3)
        IF(OverlapV <= 0.0D0) STOP 'ERROR: OverlapV is negative!'
        CubeV = (UE(1)-LE(1))*(UE(2)-LE(2))*(UE(3)-LE(3))
        CalTime = CalTime + PosTimePair%D(4,J)*OverlapV/CubeV
      ENDIF
    ENDDO
  END FUNCTION CalTime


  ! initialization needed for ET
  SUBROUTINE EqualTimeSetUp()
    INTEGER :: MaxPrNum,Ind,ReducedTotPrCount,IErr,I,J,NV,NS,AtA,AtB,NC,PrCount,PI,CI,Dir
    REAL(DOUBLE) :: LD,MD,x2,TVol,NewV,OldV,Diff
    TYPE(AtomPair) :: Pair
    REAL(DOUBLE),DIMENSION(3) :: B
    INTEGER :: P2(0:30)
    TYPE(INT_VECT) :: XYZ
    TYPE(DBL_VECT) :: POS
    TYPE(INT_VECT) :: TotPrNumArr


    CALL NEW(LC,(/3,NPrc/))
    CALL NEW(RC,(/3,NPrc/))
    P2(0) = 1
    DO I = 1, 30
      P2(I) = P2(I-1)*2
    ENDDO
    NV = NPrc
    NS = NINT(LOG(NV*1.0D0)/LOG(2.0D0))
    IF(P2(NS) /= NV) THEN
      STOP 'ERR: In ET_Part, P2 problem!'
    ENDIF
    CALL Get(GQLineLoc,'QLineLoc')

    CALL MPI_Bcast(GQLineLoc,1,MPI_INTEGER,0,MONDO_COMM,IErr)
    IF(GQLineLoc == 0) THEN
      CALL New(BegPrInd,NPrc-1,0)
      CALL New(EndPrInd,NPrc-1,0)
    ENDIF

    CALL New(XYZ,NPrc-1)
    CALL New(POS,NPrc-1)
    CALL Get(XYZ,'QETDir')
    CALL Get(Pos,'QETRoot')

    IF(MyID == 0) THEN
      TotPrCount = 0
      RootBox%BndBox(1:3,1) = 1.0D+99
      RootBox%BndBox(1:3,2) = -1.0D+99
      DO AtA = 1, NAtoms
        DO AtB = 1, NAtoms
          IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
            B = Pair%B
            DO NC = 1, CS_OUT%NCells
              Pair%B = B + CS_OUT%CellCarts%D(:,NC)
              Pair%AB2 = (Pair%A(1)-Pair%B(1))**2+(Pair%A(2)-Pair%B(2))**2+(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                CALL EBBCountJBlock(PrCount,Pair,PoleRoot)
                TotPrCount = TotPrCount + PrCount
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      GlobalCount = TotPrCount
#ifdef DIAG_QCTC
      WRITE(*,*) 'Root: TotPrCount = ',TotPrCount
#endif

      LC%D(1:3,1) = RootBox%BndBox(1:3,1)
      RC%D(1:3,1) = RootBox%BndBox(1:3,2)


#ifdef DIAG_QCTC
      WRITE(*,*) 'GQLineLoc = ',GQLineLoc
#endif
      IF(GQLineLoc < 0) THEN
        STOP 'ERR: GQLineLoc is negative!'
      ELSE IF(GQLineLoc == 0) THEN

        CALL ENPart(NPrc,GlobalCount,BegPrInd,EndPrInd)

        DO I = 1, NS
          DO PI = 1, P2(I-1)
            CI = P2(I-1) + PI
            LC%D(1:3,CI) = LC%D(1:3,PI)
            RC%D(1:3,CI) = RC%D(1:3,PI)
            MD = 0.0D0
            DO J = 1, 3
              LD = RC%D(J,PI)-LC%D(J,PI)
              IF(LD <= 0) STOP 'ERR: Negative length!'
              IF(LD > MD) THEN
                MD = LD
                Dir = J
              ENDIF
            ENDDO
            x2 = (RC%D(Dir,PI)+LC%D(Dir,PI))*0.50D0
            RC%D(Dir,PI) = x2
            LC%D(Dir,CI) = x2
          ENDDO
        ENDDO

      ELSE
        Ind = 0
        DO I = 1, NS
          DO PI = 1, P2(I-1)
            CI = P2(I-1) + PI
            LC%D(1:3,CI) = LC%D(1:3,PI)
            RC%D(1:3,CI) = RC%D(1:3,PI)
            Ind = Ind + 1
            Dir = XYZ%I(Ind)
            x2 = POS%D(Ind)
            IF(Dir < 1 .OR. Dir > 3) STOP 'ERR: Wrong Dir'
            RC%D(Dir,PI) = x2
            LC%D(Dir,CI) = x2
          ENDDO
        ENDDO
      ENDIF
      TVol = 1.0D0
      DO I = 1, 3
        LD = RootBox%BndBox(I,2)-RootBox%BndBox(I,1)
        TVol = LD*TVol
      ENDDO
      OldV = TVol
      NewV = ZERO
      DO I = 1, NV
        TVol = 1.0D0
        DO J = 1, 3
          LD = RC%D(J,I)-LC%D(J,I)
          IF(LD <= 0) THEN
            STOP 'ERR: 1:  LD not positive!'
          ENDIF
          TVol = LD*TVol
        ENDDO
        NewV = NewV + TVol
      ENDDO
      Diff = ABS(NewV-OldV)
      IF(Diff > 1.0D-8) THEN
        STOP 'ERR: Vol not conserved!'
      ENDIF

    ENDIF
    CALL MPI_Bcast(GlobalCount,1,MPI_INTEGER,0,MONDO_COMM,IErr)

    IF(GQLineLoc == 0) THEN
      CALL MPI_Bcast(BegPrInd%I(0),NPrc,MPI_INTEGER,0,MONDO_COMM,IErr)
      CALL MPI_Bcast(EndPrInd%I(0),NPrc,MPI_INTEGER,0,MONDO_COMM,IErr)
      TotPrCount = EndPrInd%I(MyID) - BegPrInd%I(MyID) + 1
    ELSE
      CALL MPI_Bcast(LC%D(1,1),3*NV,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
      CALL MPI_Bcast(RC%D(1,1),3*NV,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)

      TotPrCount = 0
      DO AtA = 1, NAtoms
        DO AtB = 1, NAtoms
          IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
            B = Pair%B
            DO NC = 1, CS_OUT%NCells
              Pair%B = B + CS_OUT%CellCarts%D(:,NC)
              Pair%AB2 = (Pair%A(1)-Pair%B(1))**2+(Pair%A(2)-Pair%B(2))**2+(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                CALL BoxEBBPairCount(PrCount,Pair,PoleRoot)
                TotPrCount = TotPrCount + PrCount
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    CALL New(TotPrNumArr,NPrc-1,0)
    CALL MPI_Gather(TotPrCount,1,MPI_INTEGER,TotPrNumArr%I(0),1,MPI_INTEGER,0,MONDO_COMM,IErr)
    ReducedTotPrCount = Reduce(TotPrCount)
    IF(MyID == 0) THEN
      MaxPrNum = -1
      DO I = 0, NPrc-1
        MaxPrNum = MAX(TotPrNumArr%I(I),MaxPrNum)
      ENDDO
#ifdef DIAG_QCTC
      CALL OpenASCII(OutFile,Out)
      WRITE(*,'(A,F10.5)') 'Max space(MB) needed for ET cost array = ',MaxPrNum*32.0D-6
      WRITE(Out,'(A,F10.5)') 'Max space(MB) needed for ET cost array = ',MaxPrNum*32.0D-6
      CLOSE(Out,STATUS='KEEP')
#endif
      IF(ReducedTotPrCount /= GlobalCount) THEN
        STOP 'ERR: Losing Primitive pairs!'
      ENDIF
    ENDIF
    CALL Delete(TotPrNumArr)
    CALL New(PosTimePair,(/4,TotPrCount/))
  END SUBROUTINE EqualTimeSetUp

  SUBROUTINE ENPart(N,TotN,BegInd,EndInd)
    INTEGER :: N,I,TotN,Num,CheckTotN
    TYPE(INT_VECT) :: BegInd,EndInd

    IF(N == 0) STOP 'ERR: N in ENPart is zero!'
    DO I = 1, N
      EndInd%I(I-1) = NINT((TotN*1.0D0/(N*1.0D0))*(I*1.0D0))
    ENDDO
    IF(EndInd%I(N-1) /= TotN) STOP 'ERR: integer problem!'
    BegInd%I(0) = 1
    DO I = 1, N-1
      BegInd%I(I) = EndInd%I(I-1)+1
    ENDDO
    CheckTotN = 0
    DO I = 0, N-1
      Num = EndInd%I(I)-BegInd%I(I)+1
      IF(Num <= 0) STOP 'ERR: too little nq to split ? '
      CheckTotN = CheckTotN + Num
    ENDDO
    if(CheckTotN /= TotN) THEN
      STOP 'ERR: missing nq ?'
    endif
  END SUBROUTINE ENPart

#endif
END MODULE ParallelQCTC
