PROGRAM GONX2
  !
#ifndef PARALLEL
#undef ONX2_PARALLEL
#endif
  !
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE Functionals
  !
  USE ONX2DataType
  USE ONXCtrSclg   , ONLY: TrnMatBlk
  USE ONX2List     , ONLY: AllocList,DeAllocList,MakeGList,PrintList
  USE GONX2ComputDK, ONLY: ComputDK

  !
#ifdef ONX2_PARALLEL
  USE MondoMPI
  USE FastMatrices
  USE ONXGet       , ONLY: Get_Essential_RowCol,GetOffArr,Set_DFASTMAT_EQ_DBCSR2
  USE PartDrv      , ONLY: PDrv_Initialize,PDrv_Finalize
#else
  USE ONXGet       , ONLY: GetOffArr
#endif
  !
  IMPLICIT NONE
  !
#ifdef ONX2_PARALLEL
  TYPE(FASTMAT), POINTER         :: DFMcd,DFMab,KxFM
#else
  TYPE(BCSR)                     :: D
#endif
  TYPE(BSET)                     :: BSc
  TYPE(CRDS)                     :: GMc
  TYPE(ARGMT)                    :: Args
!--------------------------------------------------------------------------------
#ifdef ONX2_PARALLEL
  TYPE(DBL_RNK2)                 :: GradTmp
  TYPE(INT_VECT)                 :: APt,BPt,CPt,DPt
#endif
  TYPE(DBL_RNK2)                 :: GradX,GradAux,BoxX
  TYPE(DBL_VECT)                 :: GTmp
  REAL(DOUBLE)                   :: KScale
!--------------------------------------------------------------------------------
#ifdef ONX2_PARALLEL
  TYPE(DBL_VECT)                 :: TmGxArr,TmMLArr,TmTMArr,TmALArr,TmDLArr
  INTEGER                        :: CMin,CMax,DMin,DMax,IErr
  INTEGER                        :: ANbr,BNbr,CNbr,DNbr
#endif
  REAL(DOUBLE)                   :: Time1,Time2
  REAL(DOUBLE)                   :: TmTM,TmML,TmGx,TmAL,TmDL
  CHARACTER(LEN=*),PARAMETER     :: Prog='GONX2'
  REAL(DOUBLE),EXTERNAL          :: MondoTimer
!--------------------------------------------------------------------------------
  TYPE(INT_RNK2) :: OffArr
#ifdef ONX2_PARALLEL
  TYPE(CList), DIMENSION(:), POINTER :: ListC,ListD
#else
  TYPE(CList), DIMENSION(:), POINTER :: ListC
#endif
!--------------------------------------------------------------------------------
  !
#ifdef ONX2_PARALLEL
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
#else
  CALL StartUp(Args,Prog)
#endif
  !
  CALL Get(BSiz ,'atsiz',Tag_O=PrvBase)
  CALL Get(OffS ,'atoff',Tag_O=PrvBase)
  CALL Get(NBasF,'nbasf',Tag_O=PrvBase)
  !
  CALL Get(BSc,Tag_O=CurBase)
  CALL Get(GMc,Tag_O=CurGeom)
  !
  !------------------------------------------------
  ! Initialization and allocations.
#ifdef ONX2_PARALLEL
  NULLIFY(ListC,ListD)
#else
  NULLIFY(ListC)
#endif
  !
  TmTM=Zero;TmML=Zero;TmGx=Zero;TmAL=Zero;TmDL=Zero
  CALL New(GradX,(/3,NAtoms/))
  CALL DBL_VECT_EQ_DBL_SCLR(3*NAtoms,GradX%D(1,1),0.0D0)
  !
  CALL New(GradAux,(/3,NAtoms/))
  CALL DBL_VECT_EQ_DBL_SCLR(3*NAtoms,GradAux%D(1,1),0.0D0)
  !
  CALL New(OffArr,(/BSc%NCtrt,BSc%NKind/))
  !
  CALL New(BoxX,(/3,3/))
  BoxX%D=0.0D0
  !
  CALL GetBufferSize(GMc,BSc,GMc,BSc)
  !
  CALL GetOffArr(OffArr,BSc)
  !
  !------------------------------------------------
  ! Get denstiy matrix.
  !
#ifdef ONX2_PARALLEL
  IF(MyID.EQ.ROOT) &
#endif
  WRITE(*,*) '-------- We are in GONX2 --------'
  !
  SELECT CASE(SCFActn)
  CASE('ForceEvaluation')
#ifdef ONX2_PARALLEL
     CALL PDrv_Initialize(DFMcd,TrixFile('D',Args,0),'GONXPart',Args,WhenPartS_O='GEO')
#else
     CALL Get(D,TrixFile('D',Args,0))
#endif
  CASE DEFAULT
     CALL Halt('GONX2: Do not recognize this action <'//TRIM(SCFActn)//'>.')
  END SELECT
  !
  !------------------------------------------------
  ! Normalize the density matrix
  !
#ifdef ONX2_PARALLEL
  Time1 = MondoTimer()
  CALL TrnMatBlk(BSc,GMc,DFMcd)
  Time2 = MondoTimer()
#else
  CALL CPU_TIME(Time1)
  CALL TrnMatBlk(BSc,GMc,D)
  CALL CPU_TIME(Time2)
#endif
  TmTM = Time2-Time1
  !
  !------------------------------------------------
  ! Allocate the list(s) and get the buffer sizes.
  !WRITE(*,*) 'allocate List'
#ifdef ONX2_PARALLEL
  Time1 = MondoTimer()
  CALL Get_Essential_RowCol(DFMcd,CPt,CNbr,CMin,CMax,DPt,DNbr,DMin,DMax)
  !write(*,*) 'CMin',CMin,'CMax',CMax,'MyID',MyID
  !write(*,*) 'DMin',DMin,'DMax',DMax,'MyID',MyID
  !
  CALL AllocList(ListC,CMin,CMax)
  CALL AllocList(ListD,DMin,DMax)
  Time2 = MondoTimer()
#else
  CALL CPU_TIME(Time1)
  CALL AllocList(ListC,1,NAtoms)
  CALL CPU_TIME(Time2)
#endif
  TmAL = Time2-Time1
  !WRITE(*,*) 'allocate List: ok',Time2-Time1
  !
  !------------------------------------------------
  ! Make the distribution list(s).
  !WRITE(*,*) 'make List'
#ifdef ONX2_PARALLEL
  Time1 = MondoTimer()
  CALL MakeGList(ListC,GMc,BSc,CS_OUT,CPt,CNbr,APt,ANbr)
  CALL MPI_Barrier(MONDO_COMM,IErr)
  CALL MakeGList(ListD,GMc,BSc,CS_OUT,DPt,DNbr,BPt,BNbr)
  CALL MPI_Barrier(MONDO_COMM,IErr)
  Time2 = MondoTimer()
#else
  CALL CPU_TIME(Time1)
  CALL MakeGList(ListC,GMc,BSc,CS_OUT)
  CALL CPU_TIME(Time2)
#endif
  TmML = Time2-Time1
  !WRITE(*,*) 'make List: ok',Time2-Time1
  !
  !------------------------------------------------
  ! Print list.
  !
#ifdef ONX2_DBUG
  !WRITE(*,*) 'Print List'
#ifdef ONX2_PARALLEL
  CALL PrintList(ListC)
  CALL PrintList(ListD)
#else
  CALL PrintList(ListC)
#endif
  !WRITE(*,*) 'Print List:ok'
#endif
  !
  !------------------------------------------------
  ! Get the second density matrix and normalization.
  !
#ifdef ONX2_PARALLEL
  CALL GetDab(DFMab,APt,ANbr,BPt,BNbr,Args)
  Time1 = MondoTimer()
  CALL TrnMatBlk(BSc,GMc,DFMab)
  Time2 = MondoTimer()
  TmTM = TmTM+Time2-Time1
#endif
  !
  !------------------------------------------------
  ! Compute Exchange Forces.
  !WRITE(*,*) 'DKx'
#ifdef ONX2_PARALLEL
  Time1 = MondoTimer()
  CALL ComputDK(DFMcd,DFMab,GradX,ListC,ListD,OffArr,GMc,BSc,CS_OUT)
  Time2 = MondoTimer()
#else
  CALL CPU_TIME(Time1)
  CALL ComputDK(D,GradX,BoxX,ListC,ListC,OffArr,GMc,BSc,CS_OUT)
  CALL CPU_TIME(Time2)
#endif
  TmGx = Time2-Time1
  !WRITE(*,*) 'DKx:ok',Time2-Time1
  !
  !------------------------------------------------
  ! Free up some space. Deallocate the list(s).
  !WRITE(*,*) 'deallocate List'
#ifdef ONX2_PARALLEL
  Time1 = MondoTimer()
  CALL DeAllocList(ListC)
  CALL DeAllocList(ListD)
  Time2 = MondoTimer()
#else
  CALL CPU_TIME(Time1)
  CALL DeAllocList(ListC)
  CALL CPU_TIME(Time2)
#endif
  TmDL = Time2-Time1
  !WRITE(*,*) 'deallocate List:ok',Time2-Time1
  !
  !------------------------------------------------
  ! Redistribute partition informations.
#ifdef ONX2_PARALLEL
  CALL PDrv_Finalize(DFMcd,CollectInPar_O=.TRUE.)
  CALL Delete_FastMat1(DFMcd)
  CALL Delete_FastMat1(DFMab)
#else
  CALL Delete(D)
#endif
  !
  !------------------------------------------------
  ! Add Exchange Gradient and save.
  !
  CALL Get(GradAux,'gradients',Tag_O=CurGeom)
  KScale=ExactXScale(ModelChem)
#ifdef ONX2_PARALLEL
  CALL New(GradTmp,(/3,NAtoms/))
  CALL DBL_VECT_EQ_DBL_SCLR(3*NAtoms,GradTmp%D(1,1),0.0d0)
  CALL MPI_REDUCE(GradX%D(1,1),GradTmp%D(1,1),3*NAtoms,MPI_DOUBLE_PRECISION, &
       &          MPI_SUM,ROOT,MONDO_COMM,IErr)
  ! Print out.
  CALL New(GTmp,3*NAtoms)
  CALL CartRNK2ToCartRNK1(GTmp%D,GradTmp%D)
  CALL PChkSum(GTmp,'dKx/dR['//TRIM(CurGeom)//']',Proc_O=Prog)  
  CALL Delete(GTmp)
  !
  CALL DAXPY(3*NAtoms,KScale,GradTmp%D(1,1),1,GradAux%D(1,1),1)
  CALL Delete(GradTmp)
#else
  ! Print out.
  CALL New(GTmp,3*NAtoms)
  CALL CartRNK2ToCartRNK1(GTmp%D,GradX%D)
  CALL PChkSum(GTmp,'dKx/dR['//TRIM(CurGeom)//']',Proc_O=Prog)  
  CALL Delete(GTmp)
  !
  CALL DAXPY(3*NAtoms,KScale,GradX%D(1,1),1,GradAux%D(1,1),1)
  !
#endif
  !
#ifdef 0
!StressStressStressStressStressStressStressStress
!  write(*,*) 'CS_OUT%NCells=',CS_OUT%NCells
!  write(*,*) 'LatFrc(1,1)',GMc%PBC%LatFrc%D(1,1)
!  write(*,*) 'XBox',BoxX%D(1,1)
!  write(*,*) 'Tot Stress',GMc%PBC%LatFrc%D(1,1)+BoxX%D(1,1)
!StressStressStressStressStressStressStressStress
  write(*,*) 'Grad Kx'
  do i=1,natoms
     write(*,100) i,GradX%D(:,i)
  enddo

  !write(*,*) 'Grad before'
  !do i=1,natoms
  !   write(*,100) i,GradAux%D(:,i)-GradX%D(:,i)
  !enddo

  !write(*,*) 'Grad Tot'
  !do i=1,natoms
  !   write(*,100) i,GradAux%D(:,i)
  !enddo

100 format(I4,2X,3E24.16)
#endif
  !
  CALL Put(GradAux,'gradients',Tag_O=CurGeom)
  !
  !------------------------------------------------
  ! Timing.
  !
!#ifdef GONX2_INFO
#ifdef ONX2_PARALLEL
  !
  ! End Total Timing
  CALL New(TmGxArr,NPrc)
  CALL New(TmMLArr,NPrc)
  CALL New(TmTMArr,NPrc)
  CALL New(TmALArr,NPrc)
  CALL New(TmDLArr,NPrc)
  !
  CALL MPI_Gather(TmGx,1,MPI_DOUBLE_PRECISION,TmGxArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmML,1,MPI_DOUBLE_PRECISION,TmMLArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmTM,1,MPI_DOUBLE_PRECISION,TmTMArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmAL,1,MPI_DOUBLE_PRECISION,TmALArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmDL,1,MPI_DOUBLE_PRECISION,TmDLArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  !
  IF(MyID.EQ.ROOT) THEN
     !
     ! Imbalance stuff.
     CALL PImbalance(TmGxArr ,NPrc,Prog_O='ComputeGK')
     !CALL PImbalance(TmKTArr,NPrc,Prog_O='GONX'     )
     !
     WRITE(*,1001) SUM(TmALArr%D )/DBLE(NPrc),MINVAL(TmALArr%D ),MAXVAL(TmALArr%D )
     WRITE(*,1002) SUM(TmMLArr%D )/DBLE(NPrc),MINVAL(TmMLArr%D ),MAXVAL(TmMLArr%D )
     WRITE(*,1003) SUM(TmTMArr%D )/DBLE(NPrc),MINVAL(TmTMArr%D ),MAXVAL(TmTMArr%D )
     WRITE(*,1004) SUM(TmGxArr%D )/DBLE(NPrc),MINVAL(TmGxArr%D ),MAXVAL(TmGxArr%D )
     WRITE(*,1005) SUM(TmDLArr%D )/DBLE(NPrc),MINVAL(TmDLArr%D ),MAXVAL(TmDLArr%D )
     !WRITE(*,1006) SUM(NERIsArr%D)           ,MINVAL(NERIsArr%D),MAXVAL(NERIsArr%D)
     !
1001 FORMAT(' GONX: Ave TmAL = ',F15.2,', Min TmAL = ',F15.2,', Max TmAL = ',F15.2)
1002 FORMAT(' GONX: Ave TmML = ',F15.2,', Min TmML = ',F15.2,', Max TmML = ',F15.2)
1003 FORMAT(' GONX: Ave TmTM = ',F15.2,', Min TmTM = ',F15.2,', Max TmTM = ',F15.2)
1004 FORMAT(' GONX: Ave TmGx = ',F15.2,', Min TmGx = ',F15.2,', Max TmGx = ',F15.2)
1005 FORMAT(' GONX: Ave TmDL = ',F15.2,', Min TmDL = ',F15.2,', Max TmDL = ',F15.2)
     !1006 FORMAT(' ONX: Tot ERI  = ',F15.2,', Min ERI  = ',F15.2,', Max ERI  = ',F15.2)
  ENDIF
  !
  CALL Delete(TmALArr)
  CALL Delete(TmMLArr)
  CALL Delete(TmTMArr)
  CALL Delete(TmGxArr)
  CALL Delete(TmDLArr)
  !
#else
  !
  WRITE(*,1001) TmAL
  WRITE(*,1002) TmML
  WRITE(*,1003) TmTM
  WRITE(*,1004) TmGx
  WRITE(*,1005) TmDL
  !WRITE(*,1006) ....
  !
1001 FORMAT(' GONX: Ave TmAL = ',F15.2)
1002 FORMAT(' GONX: Ave TmML = ',F15.2)
1003 FORMAT(' GONX: Ave TmTM = ',F15.2)
1004 FORMAT(' GONX: Ave TmGx = ',F15.2)
1005 FORMAT(' GONX: Ave TmDL = ',F15.2)
  !1006 FORMAT(' ONX: Tot ERI  = ',F15.2)
  !
!#endif
#endif
  !
  !------------------------------------------------
  ! Deallocate Arrays.
  !
  CALL Delete(GradX)
  CALL Delete(GradAux)
  CALL Delete(BoxX)
  CALL Delete(OffArr)
  !
  CALL Delete(BSc)
  CALL Delete(GMc)
  !
  CALL ShutDown(Prog)
  !
#ifdef ONX2_PARALLEL
CONTAINS
  SUBROUTINE GetDab(DFMab,APt,ANbr,BPt,BNbr,Args)
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(FASTMAT), POINTER    :: DFMab
    TYPE(INT_VECT)            :: APt,BPt
    INTEGER      , INTENT(IN) :: ANbr,BNBr
    TYPE(ARGMT)               :: Args
    !-------------------------------------------------------------------
    TYPE(BCSR )     :: A
    TYPE(DBCSR)     :: B
    TYPE(INT_VECT)  :: ANBrArr,BNBrArr,APtArr,BPtArr,ADisp,BDisp
    INTEGER         :: iPrc,ANbrTot,BNbrTot
    INTEGER         :: I,J,JG,M,MN,MN1,P,NAtms,NBlks,NNon0,ICol,IRow
    LOGICAL         :: ReAllocate
    !-------------------------------------------------------------------
    INTEGER, EXTERNAL :: IBinSrch
    !-------------------------------------------------------------------
    INTEGER :: TotBlk
    !
    !Send ANbr,BNbr to ROOT
    CALL New(ANBrArr,NPrc)
    CALL New(BNBrArr,NPrc)
    CALL Gather(ANBr,ANBrArr)
    CALL Gather(BNBr,BNBrArr)
    !
    IF(MyID==0) THEN
       !write(*,*) 'ANBrArr',ANBrArr%I
       !write(*,*) 'BNBrArr',BNBrArr%I
    ENDIF
    !
    ! Compute some indicies.
    IF(MyID.EQ.ROOT) THEN
       CALL New(ADisp,NPrc+1)
       CALL New(BDisp,NPrc+1)
       ANbrTot=ANBrArr%I(1)
       BNbrTot=BNBrArr%I(1)
       ADisp%I(1)=0
       BDisp%I(1)=0
       DO iPrc=2,NPrc
          ADisp%I(iPrc)=ANBrArr%I(iPrc-1)+ADisp%I(iPrc-1)
          BDisp%I(iPrc)=BNBrArr%I(iPrc-1)+BDisp%I(iPrc-1)
          ANbrTot=ANbrTot+ANBrArr%I(iPrc)
          BNbrTot=BNbrTot+BNBrArr%I(iPrc)
       ENDDO
       ADisp%I(NPrc+1)=ANbrTot
       BDisp%I(NPrc+1)=BNbrTot
       !write(*,*) 'ANbrTot',ANbrTot
       !write(*,*) 'ADisp',ADisp%I
       !write(*,*) 'BNbrTot',BNbrTot
       !write(*,*) 'BDisp',BDisp%I
       CALL New(APtArr,ANbrTot)
       CALL New(BPtArr,BNbrTot)
    ENDIF
    !
    !Alloc the arrays A,B on ROOT to get the atom lists
    CALL Gather(APt,APtArr,ANBr,ANBrArr,ADisp)
    CALL Gather(BPt,BPtArr,BNBr,BNBrArr,BDisp)
    !
    !IF(MyId==ROOT) write(*,*) 'APtArr',APtArr%I
    !IF(MyId==ROOT) write(*,*) 'BPtArr',BPtArr%I
    !
    ! Get the B matrix from disc.
    CALL Get_BCSR(A,TrixFile('D',Args,0))
    !
    !CALL PChkSum(A,'Density matrix on root',Unit_O=6)  
    !
    CALL New(B)
    !
!------------------------------------------------
!        Distribute to each processor
!
    DO iPrc=NPrc-1,0,-1
       IF(MyId==ROOT)THEN
          B%NAtms=0
          B%NBlks=1
          B%NNon0=1
          B%RowPt%I(1)=1
          !
          !vw beg and end must be set somewhere else.
          IRow=1
          DO I=1,NAtoms !Row min, Row max.
             ! Is it the right row?
             IF(IBinSrch(APtArr%I(ADisp%I(iPrc+1)+1),I,ADisp%I(iPrc+2)-ADisp%I(iPrc+1)).LT.1) CYCLE
             !if(iPrc==1) write(*,*) 'I',I
             ICol=1
             M=BSiz%I(I)
             DO J=A%RowPt%I(I),A%RowPt%I(I+1)-1
                JG=A%ColPt%I(J)
                ! Is it the right col?
                IF(IBinSrch(BPtArr%I(BDisp%I(iPrc+1)+1),JG,BDisp%I(iPrc+2)-BDisp%I(iPrc+1)).LT.1) CYCLE
                !if(iPrc==1) write(*,*) 'JG',JG
                B%ColPt%I(B%NBlks)=JG
                B%BlkPt%I(B%NBlks)=B%NNon0
                B%NBlks=B%NBlks+1
                MN=M*BSiz%I(JG);MN1=MN-1
                P=A%BlkPt%I(J)         
                B%MTrix%D(B%NNon0:B%NNon0+MN1)=A%MTrix%D(P:P+MN1) 
                B%NNon0=B%NNon0+MN
             ENDDO
             B%NAtms=B%NAtms+1
             B%RowPt%I(B%NAtms+1)=B%NBlks
          ENDDO
          B%NBlks=B%NBlks-1
          B%NNon0=B%NNon0-1
! MINUS 
          IF(iPrc/=ROOT)THEN
             CALL Send(B%NAtms,iPrc,1)               
             CALL Send(B%NBlks,iPrc,2)               
             CALL Send(B%NNon0,iPrc,3)
             CALL Send(B%RowPt,B%NAtms+1,iPrc,4)               
             CALL Send(B%ColPt,B%NBlks,iPrc,5)               
             CALL Send(B%BlkPt,B%NBlks,iPrc,6)               
             CALL Send(B%MTrix,B%NNon0,iPrc,7)
          ENDIF
       ELSEIF(MyId==iPrc)THEN
          CALL Recv(B%NAtms,ROOT,1)               
          CALL Recv(B%NBlks,ROOT,2)               
          CALL Recv(B%NNon0,ROOT,3)
          ReAllocate=(SIZE(B%RowPt%I)<B%NAtms+1).OR. &
               (SIZE(B%ColPt%I)<B%NBlks)  .OR. &
               (SIZE(B%BlkPt%I)<B%NBlks)  .OR. &
               (SIZE(B%MTrix%D)<B%NNon0)
          IF(ReAllocate)THEN     
             NAtms=B%NAtms; NBlks=B%NBlks; NNon0=B%NNon0
             CALL Delete(B)     
             CALL New(B,(/NAtms,NBlks,NNon0/))
          ENDIF
          CALL Recv(B%RowPt,B%NAtms+1,ROOT,4)               
          CALL Recv(B%ColPt,B%NBlks,ROOT,5)               
          CALL Recv(B%BlkPt,B%NBlks,ROOT,6)               
          CALL Recv(B%MTrix,B%NNon0,ROOT,7)               
       ENDIF
    ENDDO
    !------------------------------------------------
    CALL BCast(A%NBlks)
    IF(MyId==ROOT)THEN
       CALL SetEq(B%GRwPt,A%RowPt,NAtoms+1)
       CALL SetEq(B%GClPt,A%ColPt,A%NBlks)
    ENDIF
    CALL BCast(B%GRwPt,NAtoms+1)
    CALL BCast(B%GClPt,A%NBlks)
    B%Node=MyId
    B%GUpDate=STATUS_TRUE
    !------------------------------------------------
    !
    !IF(MyID==1)write(*,*) 'B MyID=',MyID,B%MTrix%D
    !IF(MyID==1)write(*,*) 'B MyID=',MyID,B%RowPt%I
    !IF(MyID==1)write(*,*) 'B MyID=',MyID,B%ColPt%I
    !
    IF(MyID.EQ.ROOT) TotBlk=A%NBlks
    CALL BCast(TotBlk)
    !
    WRITE(*,100) DBLE(B%NBlks)/DBLE(TotBlk)*100d0,DBLE(B%NBlks)/(DBLE(NAtoms)**2)*100d0,MyID
100 FORMAT(' Remaining block, Rel=',F8.3,'% Abs=',F8.3,'% MyID=',I4)
    !
    ! Copy the DBCSR to a FastMat.
    CALL New_FASTMAT(DFMab,0,(/0,0/))
    CALL Set_DFASTMAT_EQ_DBCSR2(DFMab,B,APt,ANBr)
    !CALL PChkSum_FASTMAT2(DFMab,'Density matrix ab',Unit_O=6)
    !
    ! Delete arrays.
    IF(MyID.EQ.ROOT) THEN
       CALL Delete(APtArr)
       CALL Delete(BPtArr)
       CALL Delete(A)
    ENDIF
    !
    CALL Delete(B)
    CALL Delete(ANBrArr)
    CALL Delete(BNBrArr)
    !
  END SUBROUTINE GetDab
  !
  !
!!$  SUBROUTINE Set_DFASTMAT_EQ_DBCSR2(B,A,RowPt,RowNbr)
!!$    TYPE(FASTMAT),POINTER :: B,C
!!$    TYPE(SRST),POINTER    :: S
!!$    TYPE(DBCSR)           :: A
!!$    TYPE(INT_VECT)        :: RowPt
!!$    INTEGER, INTENT(IN)   :: RowNbr
!!$    INTEGER               :: I,IRow,J,JP,P,N,M
!!$    !---------------------------------------------------------------------
!!$    !
!!$    ! Set some pointers.
!!$    NULLIFY(C,S)
!!$    !
!!$    ! Check for prior allocation
!!$    IF(ASSOCIATED(B))THEN
!!$       CALL Delete_FASTMAT1(B)
!!$       ! Begin with a new header Node     
!!$       CALL New_FASTMAT(B,0,(/0,0/))
!!$    ENDIF
!!$    !
!!$    DO I = 1,RowNbr
!!$       IRow=RowPt%I(I)
!!$       M = BSiz%I(IRow)
!!$       IF(A%RowPt%I(I+1)-A%RowPt%I(I)>1) THEN
!!$          ! Set current row link
!!$          C => FindFastMatRow_1(B,IRow)
!!$          DO JP = A%RowPt%I(I),A%RowPt%I(I+1)-1
!!$             J = A%ColPt%I(JP)
!!$             P = A%BlkPt%I(JP)
!!$             N = BSiz%I(J)
!!$             ! Add bloks to this sparse row search tree
!!$             CALL AddFASTMATBlok(C,IRow,J,RESHAPE(A%MTrix%D(P:P+M*N-1),(/M,N/)))
!!$          ENDDO
!!$       ENDIF
!!$    ENDDO
!!$    !
!!$    CALL FlattenAllRows(B)
!!$    !
!!$  END SUBROUTINE Set_DFASTMAT_EQ_DBCSR2
  !
  SUBROUTINE PChkSum_FASTMAT2(A,Name,Unit_O,Proc_O,ChkInPar_O)
!H---------------------------------------------------------------------------------
!H SUBROUTINE  PChkSum_FASTMAT(A,Name,Unit_O,Proc_O,ChkInPar_O)
!H
!H---------------------------------------------------------------------------------
    TYPE(FASTMAT)    , POINTER              :: A
    CHARACTER(LEN=*) , INTENT(IN)           :: Name
    CHARACTER(LEN=*) , INTENT(IN), OPTIONAL :: Proc_O
    INTEGER          , INTENT(IN), OPTIONAL :: Unit_O
    LOGICAL          , INTENT(IN), OPTIONAL :: ChkInPar_O
    !-------------------------------------------------------------------
    TYPE(FASTMAT)    , POINTER              :: R
    TYPE(SRST   )    , POINTER              :: U
    INTEGER                                 :: I,PU,J,M,N
    REAL(DOUBLE)                            :: Chk
    LOGICAL                                 :: InPara
    CHARACTER(LEN=DEFAULT_CHR_LEN)       :: ChkStr
    !-------------------------------------------------------------------
    REAL(DOUBLE), EXTERNAL               :: DDOT
    !-------------------------------------------------------------------
    !
    ! IF(PrintFlags%Key/=DEBUG_MAXIMUM.AND. &
    !  PrintFlags%Chk/=DEBUG_CHKSUMS)RETURN
    !
    NULLIFY(R,U)
    Chk=Zero
    !
    ! Flatten A.
    CALL FlattenAllRows(A)
    !
    R => A%Next
    DO
       IF(.NOT.ASSOCIATED(R)) EXIT
       I = R%Row
       M = BSiz%I(I)
       U => R%RowRoot
       DO
          IF(.NOT.ASSOCIATED(U)) EXIT
          IF(U%L.EQ.U%R) THEN
             IF(ASSOCIATED(U%MTrix)) THEN
                J = U%L
                N = BSiz%I(J)
                Chk=Chk+DDOT(M*N,U%MTrix(1,1),1,U%MTrix(1,1),1)
             ENDIF
          ENDIF
          U => U%Next
       ENDDO
       R => R%Next
    ENDDO
    Chk=SQRT(Chk) 
    ChkStr=CheckSumString(Chk,Name,Proc_O)
    ! Write check string
    ! PU=OpenPU(Unit_O=Unit_O)
    ! WRITE(PU,'(1x,A)')TRIM(ChkStr)
    WRITE(*,'(1x,A,I3)')TRIM(ChkStr),MyID
    ! CALL ClosePU(PU)
    !
  END SUBROUTINE PChkSum_FASTMAT2

#endif

END PROGRAM GONX2

