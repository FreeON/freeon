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
MODULE PartDrv
!H=================================================================================
!H MODULE PartDrv
!H This MODULE contains:
!H  PUBLIC:
!H  o SUB PDrv_Initialize
!H  o SUB PDrv_Finalize
!H
!H  PRIVATE:
!H  o SUB PDrv_Part_1D
!H  o SUB PDrv_Part_2D
!H  o SUB PDrv_Part_Level1
!H  o SUB PDrv_Part_Level2
!H  o SUB PDrv_Load_Level2
!H  o SUB PDrv_Gather_Part
!H  o SUB Set_Universal_BegEnd
!H  o SUB Save_Universal_BegEnd
!H  o SUB Load_Universal_BegEnd
!H  o SUB Get_RowDist
!H  o SUB Get_ColDist
!H  o SUB Get_RowColDist
!H  o SUB Get_Sum_Row
!H  o SUB Get_Dist_Prc_2D
!H  o SUB Set_DBCSR_EQ_BCSR_Part
!H  o SUB Get_DBCSR_Part
!H  o FUN Get_Optimal_Bound
!H  o FUN BinSrch
!H  o FUN Probe
!H
!H  OPTIONS:
!H  DEBUGING: Use -DPARTDRV_DBUG to print some stuff.
!H  INFO    : Use -DPARTDRV_INFO to print some stuff.
!H
!H  Comments:
!H
!H---------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
! vw comments:
! o clean up the code.
! o Should replace the format in PDrv_Gather_Part by CSR.
! o Should replace the format in Set_FASTMAT_EQ_SHITY -> Set_FASTMAT_EQ_CSR.
!
! TODO TODO TODO TODO TODO TODO TODO TODO
! x There are problems now with NAtoms ~ NPrc with 1D part., should
!   do some tests with that, if <= 0 then use Beg and End, + WARNING.
!----------------------------------------------------------------------------------
  !
#ifndef PARALLEL
#undef ONX2_PARALLEL
#endif
  !
#ifdef ONX2_PARALLEL
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE MondoMPI
  USE FastMatrices
  USE InOut
  !
  IMPLICIT NONE
  PRIVATE
  !
!---------------------------------------------------------------------------------
! PUBLIC DECLARATIONS
!---------------------------------------------------------------------------------
  PUBLIC  :: PDrv_Initialize
  PUBLIC  :: PDrv_Finalize
  !
!---------------------------------------------------------------------------------
! PRIVATE DECLARATIONS
!---------------------------------------------------------------------------------
  PRIVATE :: PDrv_Part_1D
  PRIVATE :: PDrv_Part_2D
  PRIVATE :: PDrv_Part_Level1
  PRIVATE :: PDrv_Part_Level2
  PRIVATE :: PDrv_Load_Level2
  PRIVATE :: PDrv_Gather_Part
  PRIVATE :: Set_Universal_BegEnd
  PRIVATE :: Save_Universal_BegEnd
  PRIVATE :: Load_Universal_BegEnd
  PRIVATE :: Get_RowDist
  PRIVATE :: Get_ColDist
  PRIVATE :: Get_RowColDist
  PRIVATE :: Get_Sum_Row
  PRIVATE :: Get_DistPrc_2D
  PRIVATE :: Set_DBCSR_EQ_BCSR_Part
  PRIVATE :: Get_DBCSR_Part
  PRIVATE :: Get_Optimal_Bound
  PRIVATE :: BinSrch
  PRIVATE :: Probe
  !
!---------------------------------------------------------------------------------
! SEMI-GLOBAL VARIABLES
!---------------------------------------------------------------------------------
  LOGICAL                       , PRIVATE, SAVE :: IsFirst      = .TRUE.
  LOGICAL                       , PRIVATE, SAVE :: PartAtLevel1 = .TRUE.
  TYPE(INT_VECT)                , PRIVATE, SAVE :: Beg_Save    ! Temporary buffer for Beg
  TYPE(INT_VECT)                , PRIVATE, SAVE :: End_Save    ! Temporary buffer for End
  CHARACTER(LEN=DEFAULT_CHR_LEN), PRIVATE, SAVE :: ParType
  CHARACTER(LEN=DEFAULT_CHR_LEN), PRIVATE, SAVE :: PartName
  !
!---------------------------------------------------------------------------------
! SEMI-GLOBAL PARAMETER
!---------------------------------------------------------------------------------
  CHARACTER(LEN=12), PARAMETER, PRIVATE :: PART_NO_PART  = 'Part_No_Part'
  CHARACTER(LEN=13), PARAMETER, PRIVATE :: PART_NZERO_1D = 'Part_NZero_1D'
  CHARACTER(LEN=13), PARAMETER, PRIVATE :: PART_NZERO_2D = 'Part_NZero_2D'
  CHARACTER(LEN=12), PARAMETER, PRIVATE :: PART_TIME_1D  = 'Part_Time_1D'
  CHARACTER(LEN=12), PARAMETER, PRIVATE :: PART_TIME_2D  = 'Part_Time_2D'
  !
CONTAINS
  !
  SUBROUTINE PDrv_Initialize(AFastMat,Name,PartN,Args,PartS_O,FirstPartS_O,WhenPartS_O,PFix_O,CheckPoint_O)
!H---------------------------------------------------------------------------------
!H SUBROUTINE PDrv_Initialize(AFastMat,Name,PartN,Args,PartS_O,FirstPartS_O,
!H                            WhenPartS_O,PFix_O,CheckPoint_O)
!H
!H---------------------------------------------------------------------------------
    USE Parse
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(FASTMAT)             , POINTER    :: AFastMat
    TYPE(ARGMT  )             , INTENT(IN) :: Args
    CHARACTER(LEN=*)          , INTENT(IN) :: Name,PartN
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: PFix_O
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: FirstPartS_O,PartS_O,WhenPartS_O
    LOGICAL,          OPTIONAL, INTENT(IN) :: CheckPoint_O
    !-------------------------------------------------------------------
    TYPE(DBCSR  )                          :: A
    CHARACTER(LEN=DEFAULT_CHR_LEN)         :: FirstPartS,PartS,WhenPartS
    INTEGER                                :: Cycl,Basis,Geom,iONXPartExist,iGONXPartExist
    INTEGER :: iErr
    !-------------------------------------------------------------------
    !
    ! Initialize some variables.
    Cycl  = Args%I%I(1)
    Basis = Args%I%I(2)
    Geom  = Args%I%I(3)
    !
    ! Open Input File
    CALL OpenASCII(InpFile,Inp)
    !
    ! When do we need to do the partition?
    WhenPartS='SCF'
    IF(PRESENT(WhenPartS_O)) WhenPartS=WhenPartS_O
    !
    ! Initialize semi-global variables.
    IsFirst = .TRUE.
    SELECT CASE(WhenPartS)
    CASE('SCF')
       !IF(OptKeyQ(Inp,'Guess','Core')) THEN
       !   !We use a Core Guess.
       !   IsFirst = Cycl.LE.1.AND.Basis.LE.1!.AND.Geom.LE.1
       !ELSE
       !   !We do not use a Core Guess.
       !   IsFirst = Cycl.LE.0.AND.Basis.LE.1!.AND.Geom.LE.1
       !ENDIF
       iONXPartExist=0
       CALL Get(iONXPartExist,'ONXPartExist')
       IsFirst=iONXPartExist.LE.0
       !
       ! It is no more the first relative iteration.
       IF(IsFirst) THEN
          iONXPartExist=NPrc
          CALL Put(iONXPartExist,'ONXPartExist')
       ENDIF
    CASE('GEO')
       iGONXPartExist=0
       CALL Get(iGONXPartExist,'GONXPartExist')
       IsFirst=iGONXPartExist.LE.0
       !
       ! It is no more the first relative iteration.
       IF(IsFirst) THEN
          iGONXPartExist=NPrc
          CALL Put(iGONXPartExist,'GONXPartExist')
       ENDIF
    CASE DEFAULT
       CALL Halt('PartDrv: Does not regonize when to do the Partition <'//TRIM(WhenPartS)//'>.')
    END SELECT
    !
#ifdef PARTDRV_DBUG
    IF(MYID.EQ.ROOT) THEN
       WRITE(*,*) 'Cycl',Cycl
       WRITE(*,*) 'Basis',Basis
       WRITE(*,*) 'Geom',Geom
       WRITE(*,*) 'IsFirst',IsFirst
       WRITE(*,*) 'Do we use Core Guess? ',OptKeyQ(Inp,'Guess','Core')! OptKeyQ(Inp,GUESS_OPTION,GUESS_CORE)
    ENDIF
#endif
    !
    !Close the Input file.
    CLOSE(Inp)
    !
    IF(PRESENT(FirstPartS_O)) THEN
       FirstPartS = FirstPartS_O
       IF(TRIM(FirstPartS).NE.PART_NZERO_1D.OR.TRIM(FirstPartS).NE.PART_NZERO_2D) &
            & STOP 'Partition in PDrv_Initialize not supported.'
    ELSE
       FirstPartS = PART_NZERO_2D!PART_NZERO_1D !PART_NZERO_2D!
    ENDIF
    IF(PRESENT(PartS_O)) THEN
       PartS = PartS_O
    ELSE
       PartS = PART_TIME_2D!PART_TIME_1D !PART_TIME_2D!
    ENDIF
    !
    PartName = PartN
    !
    ! Quick return if needed.
    IF(TRIM(PartS).EQ.PART_NO_PART.OR.TRIM(FirstPartS).EQ.PART_NO_PART) RETURN
    !
    ! Set variable for part. at level 1 only.
    PartAtLevel1 = TRIM(PartS).EQ.PART_NZERO_1D.OR.TRIM(PartS).EQ.PART_NZERO_2D!.OR.IsFirst
    !
    ! Select the partition.
    ParType = PartS
    IF(IsFirst) ParType = FirstPartS
    !
#ifdef PARTDRV_INFO
    IF(MyID.EQ.ROOT) THEN
       WRITE(*,*) '________________________________________________________'
       WRITE(*,*) ''
       WRITE(*,*) '  Is First Cycle',IsFirst,'. Is Partition at Level1 only', &
            &        PartAtLevel1
       WRITE(*,*) '________________________________________________________'
    ENDIF
#endif
    !
    ! Save the current Beg and End.
    CALL Save_Universal_BegEnd()
    !
    ! Get the partition and send the partitioned matrix to the Proc's.
    CALL Get_DBCSR_Part(A,Name,ParType,PFix_O)
    !
    !CALL CheckSum_DBCSR2(A,'A')
    ! Copy the DBCSR to a FastMat.
    CALL New_FASTMAT(AFastMat,0,(/0,0/),NSMat_O=A%NSMat)
    CALL Set_DFASTMAT_EQ_DBCSR(AFastMat,A)
    !CALL PChkSum_FASTMAT2(AFastMat,'AFastMat')
    CALL Delete(A)
    !
    ! Select the next partition if first.
    IF(IsFirst) ParType = PartS
    !
    ! if needed...
    IF(.NOT.PartAtLevel1) CALL SetFASTMATPart(AFastMat)
    !
  END SUBROUTINE PDrv_Initialize
  !
  !
  SUBROUTINE PDrv_Finalize(AFastMat,CollectInPar_O)
!H---------------------------------------------------------------------------------
!H SUBROUTINE PDrv_Finalize(AFastMat,CollectInPar_O)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(FASTMAT), POINTER              :: AFastMat
    LOGICAL      , OPTIONAL, INTENT(IN) :: CollectInPar_O
    !-------------------------------------------------------------------
    TYPE(FASTMAT), POINTER              :: BFastMat
    LOGICAL                             :: CollectInPar
    !-------------------------------------------------------------------
    NULLIFY(BFastMat)
    !
    ! Test if we need to do a partition at level 1.
    IF(.NOT.PartAtLevel1) THEN
       !
       IF(PRESENT(CollectInPar_O)) THEN
          CollectInPar = CollectInPar_O
       ELSE
          CollectInPar = .TRUE.
       ENDIF
       !
       !If we needed to the collect from Proc's
       IF(CollectInPar) THEN
          CALL PDrv_Gather_Part(AFastMat,BFastMat)
       ELSE
          BFastMat => AFastMat
       ENDIF
       !
       ! Now we have the fastmat partition on root
       ! Do partitioning at level 2.
       CALL PDrv_Part_Level2(BFastMat,ParType)
       !
       ! Delete tempory FastMat.
       IF(CollectInPar) THEN
          CALL Delete_FASTMAT1(BFastMat)
       ELSE
          NULLIFY(BFastMat)
       ENDIF
       !
    ENDIF
    !
    ! Send the old Beg and End.
    CALL Load_Universal_BegEnd()
    !
    ! BroadCast Beg and End.
    CALL BCast(Beg)
    CALL BCast(End)
    !
  END SUBROUTINE PDrv_Finalize
  !-------------------------------------------------------------------------------
  !     Get a DBCSR matrix
  SUBROUTINE Get_DBCSR_Part(A,Name,ParType,PFix_O,CheckPoint_O)
    IMPLICIT NONE
    TYPE(DBCSR)                               :: A
    LOGICAL,          OPTIONAL, INTENT(IN   ) :: CheckPoint_O
    CHARACTER(LEN=*),           INTENT(IN   ) :: Name,ParType
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN   ) :: PFix_O
    TYPE(BCSR)                                :: B
#ifdef ONX2_PARALLEL
    TYPE(FASTMAT), POINTER                    :: BFastMat
    TYPE(INT_VECT)                            :: PrcDist
    LOGICAL                                   :: InParTemp
    INTEGER                                   :: NBlks,NNon0
    integer :: TotNBlks,TotNNon0
#endif
    !-------------------------------------------------------------------------------
    NULLIFY(BFastMat)
#ifdef ONX2_PARALLEL
    IF(PRESENT(Checkpoint_O))THEN
       InParTemp=InParallel
       ! We must turn off the parallel broadcast at this point
       ! so that we only gather from HDF to the ROOT node
       InParallel=.FALSE.
    ENDIF
#endif
    !
    ! Get the B matrix from disc.
    CALL Get_BCSR(B,Name,PFix_O,CheckPoint_O)
    !
#ifdef ONX2_PARALLEL
    ! Do partition on Root.
    IF(MyID==0) THEN
       NBlks = B%NBlks
       NNon0 = B%NNon0
       !
       ! Allocate distribution array.
       CALL New(PrcDist,NPrc*4,0)
       CALL SetEq(PrcDist,0)
       !
       !
       IF(IsFirst.OR.PartAtLevel1) THEN
          ! Copy BCSR to FastMat.
          CALL Set_FASTMAT_EQ_BCSR(BFastMat,B)
          ! Do partitioning at level1.
          CALL PDrv_Part_Level1(BFastMat,PrcDist,ParType)!,'Part_NZero_2D')! ParType)
          ! Delete tempory FastMat.
          CAll Delete_FASTMAT1(BFastMat)
          !
       ELSE
          ! Load partitioning at level2.
          CALL PDrv_Load_Level2(PrcDist)
          !
       ENDIF
       !
       ! Set the new partition pointer.
       CALL Set_Universal_BegEnd(PrcDist)
       !
       !
    ENDIF
    IF(PRESENT(Checkpoint_O)) THEN
       InParallel=InParTemp
    ENDIF
    !
    ! BroadCast the new Beg and End Array.
    CALL BCast(Beg)
    CALL BCast(End)
#endif
    !
    ! Send the BCSR -> DBCSR
    !old CALL SetEq(A,B)
    CALL Set_DBCSR_EQ_BCSR_Part(A,B,PrcDist)
    A%Node = MyId
    !
    ! Checking if DM is well distributed.
    TotNBlks=Reduce(A%NBlks)
    TotNNon0=Reduce(A%NNon0)
    IF(MyID.EQ.ROOT) THEN
       IF(TotNBlks.NE.B%NBlks) THEN
          WRITE(*,*) 'TotNBlks=',TotNBlks,' B%NBlks=',B%NBlks
          STOP 'Err1: Not consitent distribution of the DM'
       ENDIF
       IF(TotNNon0.NE.B%NNon0) THEN
          WRITE(*,*) 'TotNNon0=',TotNNon0,' B%NNon0=',B%NNon0
          STOP 'Err2: Not consitent distribution of the DM'
       ENDIF
    ENDIF
    !
    ! Delete BCSR on Root.
    CALL Delete(B)

    ! Delete tempory array.
    IF(MyID.EQ.ROOT) CALL Delete(PrcDist)
    !
!    ! Transform DBCSR into FastMat (for sub Get_FASTMAT_Part).
!#ifdef PARALLEL_ONX
!    !CALL New_FASTMAT(DFastMat,0,(/0,0/))
!    !CALL Set_DFASTMAT_EQ_DBCSR(DFastMat,D)
!    !CALL Delete(D)
!#endif
  END SUBROUTINE Get_DBCSR_Part
  !
  !
!============================================================================
!     Scatter a serial BCSR matrix from ROOT to a distributed BCSR matrix
!============================================================================
  SUBROUTINE Set_DBCSR_EQ_BCSR_Part(B,A,PrcDist)
    TYPE(BCSR)                  :: A
    TYPE(DBCSR)                 :: B
    INTEGER                     :: I,Id,J,JG,M,MN,MN1,P,  &
         &                         NAtms,NBlks,NNon0!,DUM(3)
    LOGICAL                     :: ReAllocate
    TYPE(INT_VECT)              :: PrcDist
    !TYPE(INT_VECT) :: V
!-----------------------------------------------------------------------
!        Allocate if required
!
    CALL BCast(A%NSMat) !<<< SPIN
    IF(.NOT.AllocQ(B%Alloc))CALL New(B,NSMat_O=A%NSMat) !<<< SPIN
!------------------------------------------------
!        Distribute to each processor
!
    DO Id=NPrc-1,0,-1
       IF(MyId==ROOT)THEN
          B%NAtms=0
          B%NBlks=1
          B%NNon0=1
          B%RowPt%I(1)=1
          !
          !vw beg and end must be set somewhere else.
          DO I=PrcDist%I(Id*4),PrcDist%I(Id*4+1) !Row min, Row max.
             M=BSiz%I(I)
             DO J=A%RowPt%I(I),A%RowPt%I(I+1)-1
                JG=A%ColPt%I(J)
                IF(JG.LT.PrcDist%I(Id*4+2)) CYCLE !Col min.
                IF(JG.GT.PrcDist%I(Id*4+3)) CYCLE !Col max.
                B%ColPt%I(B%NBlks)=JG
                B%BlkPt%I(B%NBlks)=B%NNon0
                B%NBlks=B%NBlks+1
                MN=M*BSiz%I(JG)*A%NSMat !<<< SPIN
                MN1=MN-1
                P=A%BlkPt%I(J)
                !CALL DBL_VECT_EQ_DBL_VECT(MN1,B%MTrix%D(B%NNon0),A%MTrix%D(P))
                B%MTrix%D(B%NNon0:B%NNon0+MN1)=A%MTrix%D(P:P+MN1)
                B%NNon0=B%NNon0+MN
             ENDDO
             B%NAtms=B%NAtms+1
             B%RowPt%I(B%NAtms+1)=B%NBlks
          ENDDO
          B%NBlks=B%NBlks-1
          B%NNon0=B%NNon0-1
! MINUS
          IF(Id/=ROOT)THEN
             !DUM(1)=B%NAtms
             !DUM(2)=B%NBlks
             !DUM(3)=B%NNon0
             !CALL MPI_SEND(DUM(1),3,MPI_INTEGER,Id,1,MONDO_COMM,iErr)
             !N=B%NAtms+1+2*B%NBlks
             !CALL INT_VECT_EQ_INT_VECT(B%NAtms+1,V%I(1)                ,B%RowPt%I(1))
             !CALL INT_VECT_EQ_INT_VECT(B%NBlks  ,V%I(B%NAtms+2)        ,B%ColPt%I(1))
             !CALL INT_VECT_EQ_INT_VECT(B%NBlks  ,V%I(B%NAtms+2+B%NBlks),B%BlkPt%I(1))
             !CALL MPI_SEND(V%I(1),N,MPI_INTEGER,Id,2,MONDO_COMM,iErr)
             !CALL MPI_SEND(B%MTrix%D(1),B%NNon0,MPI_DOUBLE_PRECISION,Id,7,MONDO_COMM,iErr)
             CALL Send(B%NAtms,Id,1)
             CALL Send(B%NBlks,Id,2)
             CALL Send(B%NNon0,Id,3)
             CALL Send(B%RowPt,B%NAtms+1,Id,4)
             CALL Send(B%ColPt,B%NBlks,Id,5)
             CALL Send(B%BlkPt,B%NBlks,Id,6)
             CALL Send(B%MTrix,B%NNon0,Id,7)
          ENDIF
       ELSEIF(MyId==Id)THEN
          !CALL MPI_RECV(DUM(1),3,MPI_INTEGER,ROOT,1,MONDO_COMM,MPI_STATUS_IGNORE,iErr)
          !B%NAtms=DUM(1)
          !B%NBlks=DUM(2)
          !B%NNon0=DUM(3)
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
          !N=B%NAtms+1+2*B%NBlks
          !CALL MPI_RECV(V%I(1),N,MPI_INTEGER,ROOT,2,MONDO_COMM,MPI_STATUS_IGNORE,iErr)
          !CALL INT_VECT_EQ_INT_VECT(B%NAtms+1,B%RowPt%I(1),V%I(1))
          !CALL INT_VECT_EQ_INT_VECT(B%NBlks  ,B%ColPt%I(1),V%I(B%NAtms+2))
          !CALL INT_VECT_EQ_INT_VECT(B%NBlks  ,B%BlkPt%I(1),V%I(B%NAtms+2+B%NBlks))
          !CALL MPI_RECV(B%MTrix%D(1),B%NNon0,MPI_DOUBLE_PRECISION, &
          !         &    ROOT,7,MONDO_COMM,MPI_STATUS_IGNORE,iErr)
          CALL Recv(B%RowPt,B%NAtms+1,ROOT,4)
          CALL Recv(B%ColPt,B%NBlks,ROOT,5)
          CALL Recv(B%BlkPt,B%NBlks,ROOT,6)
          CALL Recv(B%MTrix,B%NNon0,ROOT,7)
       ENDIF
    ENDDO
!-------------------------------------------------------
    CALL BCast(A%NBlks)
    IF(MyId==ROOT)THEN
       CALL SetEq(B%GRwPt,A%RowPt,NAtoms+1)
       CALL SetEq(B%GClPt,A%ColPt,A%NBlks)
    ENDIF
    CALL BCast(B%GRwPt,NAtoms+1)
    CALL BCast(B%GClPt,A%NBlks)
    B%Node=MyId
    B%GUpDate=STATUS_TRUE
    !
  END SUBROUTINE Set_DBCSR_EQ_BCSR_Part
  !
  !
  SUBROUTINE PDrv_Part_Level1(A,PrcDist,ParType)
!H---------------------------------------------------------------------------------
!H SUBROUTINE PDrv_Part_Level1(A,PrcDist,ParType)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(FASTMAT )  , POINTER       :: A
    TYPE(INT_VECT)                  :: PrcDist
    CHARACTER(LEN=*), INTENT(IN   ) :: ParType
    !-------------------------------------------------------------------
    !
    IF(MyID.NE.ROOT) RETURN
    !
#ifdef PARTDRV_INFO
    WRITE(*,*) '___________________________________________'
    WRITE(*,*) ''
    WRITE(*,*) '  We are doing a ',TRIM(ParType),' partition.'
    WRITE(*,*) '___________________________________________'
#endif
    !
    ! First level partitioning.
    SELECT CASE(TRIM(ParType))
    CASE('Part_NZero_1D'); CALL PDrv_Part_1D(A,'Non0',PrcDist)
    CASE('Part_NZero_2D'); CALL PDrv_Part_2D(A,'Non0',PrcDist)
    CASE DEFAULT
       CALL Halt(' Partition Scheme NOT supported in PDrv_Part, STOP. ')
    END SELECT
    !
  END SUBROUTINE PDrv_Part_Level1
  !
  !
  SUBROUTINE PDrv_Part_Level2(A,ParType)
!H---------------------------------------------------------------------------------
!H SUBROUTINE PDrv_Part_Level2(A,ParType)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(FASTMAT)   , POINTER    :: A
    CHARACTER(LEN=*), INTENT(IN) :: ParType
    !-------------------------------------------------------------------
    TYPE(INT_VECT)               :: PrcDist
    !-------------------------------------------------------------------
    !
    ! To be sure.
    IF(MyID.NE.ROOT) RETURN
    !
#ifdef PARTDRV_INFO
    WRITE(*,*) '__________________________________________'
    WRITE(*,*) ''
    WRITE(*,*) '  We are doing a ',TRIM(ParType),' partition.'
    WRITE(*,*) '__________________________________________'
#endif
    !
    !
    CALL New(PrcDist,NPrc*4,0)
    CALL SetEq(PrcDist,0)
    !
    ! Second level partitioning.
    SELECT CASE(TRIM(ParType))
    CASE('Part_Time_1D'); CALL PDrv_Part_1D(A,'Time',PrcDist)
    CASE('Part_Time_2D'); CALL PDrv_Part_2D(A,'Time',PrcDist)
    CASE DEFAULT
       CALL Halt(' Partition Scheme NOT supported in PDrv_Part, STOP. ')
    END SELECT
    !
    ! Delete disribution array.
    CALL Delete(PrcDist)
    !
  END SUBROUTINE PDrv_Part_Level2
  !
  !
  SUBROUTINE PDrv_Load_Level2(PrcDist)
!H---------------------------------------------------------------------------------
!H SUBROUTINE PDrv_Load_Level2(PrcDist)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(INT_VECT)                :: PrcDist
    !-------------------------------------------------------------------
    INTEGER                       :: iPrc
    LOGICAL                       :: InParTemp
    !-------------------------------------------------------------------
    !
    IF(.NOT.AllocQ(PrcDist%Alloc)) CALL New(PrcDist,NPrc*4,0)
    CALL SetEq(PrcDist,0)
    !
    ! We must turn off the parallel broadcast.
    InParTemp  = InParallel
    InParallel = .FALSE.
    !
    ! Get the old 1D-level2 distribution.
    CALL Get(PrcDist,TRIM(PartName))
    !CALL Get(PrcDist,'ONXPart')
    InParallel = InParTemp
    !
    !-TODO- Do an easy check here.
    !
    !
  END SUBROUTINE PDrv_Load_Level2
  !
  !
  SUBROUTINE PDrv_Gather_Part(A,B)
!H---------------------------------------------------------------------------------
!H SUBROUTINE PDrv_Gather_Part(A,B)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(FASTMAT ), POINTER :: A,B
    !-------------------------------------------------------------------
    TYPE(FASTMAT ), POINTER :: P
    TYPE(SRST    ), POINTER :: U
    TYPE(INT_VECT)          :: RowTot,ColTot,RowLoc,ColLoc,VNbrNode,Disp
    TYPE(DBL_VECT)          :: ParTot,ParLoc
    INTEGER                 :: NbrNode,NbrNodeTot,iPrc,I
    INTEGER                 :: Col,Row
    !-------------------------------------------------------------------
    !
    IF(.NOT.ASSOCIATED(A)) CALL Halt(' A is null in PDrv_Collect, STOP. ')
    !
    NULLIFY(P)
    NULLIFY(U)
    !
    P => A%Next
    NbrNode = 0
    ! Get number of nodes in the FASTMAT. !Can be done in another way.
    DO
       IF(.NOT.ASSOCIATED(P)) EXIT
       U => P%RowRoot
       DO
          IF(.NOT.ASSOCIATED(U)) EXIT
          IF(U%L == U%R) THEN
             NbrNode = NbrNode+1
          ENDIF
          U => U%Next
       ENDDO
       P => P%Next
    ENDDO
    !
    ! Create array on ROOT.
    IF(MyID.EQ.ROOT) THEN
       CALL New(VNbrNode,NPrc-1,M_O=0)
       VNbrNode%I(:) = 0
    ENDIF
    !
    ! Allocate local buffers.
    CALL New(RowLoc,NbrNode)
    CALL New(ColLoc,NbrNode)
    CALL New(ParLoc,NbrNode)
    !
    ! Fill Row, Col, Partition array.
    ! A must be associated!
    I = 0
    P => A%Next
    DO
       IF(.NOT.ASSOCIATED(P)) EXIT
       U => P%RowRoot
       Row = P%Row
       DO
          IF(.NOT.ASSOCIATED(U)) EXIT
          IF(U%L == U%R) THEN
             I = I+1
             Col = U%R
             !
             ! Collect local info.
             RowLoc%I(I) = Row
             ColLoc%I(I) = Col
             ParLoc%D(I) = U%Part
          ENDIF
          U => U%Next
       ENDDO
       P => P%Next
    ENDDO
    !
    ! Collect the number of nodes per Proc.
    CALL Gather(NbrNode,VNbrNode)
    !
    ! Compute some indicies.
    IF(MyID.EQ.ROOT) THEN
       CALL New(Disp,NPrc-1,M_O=0)
       NbrNodeTot = VNbrNode%I(0)
       Disp%I(0) = 0
       DO iPrc = 1,NPrc-1
          Disp%I(iPrc) = VNbrNode%I(iPrc-1)+Disp%I(iPrc-1)
          NbrNodeTot = NbrNodeTot+VNbrNode%I(iPrc)
       ENDDO
    ENDIF
    !
    ! Allocate Global buffer on Root (Can be done in a better way!)
    ! (MUST BE CHANGED!, should use st like DBCSR).
    IF(MyID.EQ.ROOT) THEN
       CALL New(RowTot,NbrNodeTot)
       CALL New(ColTot,NbrNodeTot)
       CALL New(ParTot,NbrNodeTot)
    ENDIF
    !
    ! Print some stuff if needed.
#ifdef PARTDRV_DBUG
    IF(MYID==0) write(*,*) 'Root ','NbrNodeTot',NbrNodeTot
    IF(MYID==0) write(*,*) 'Disp',Disp%I
    IF(MYID==0) write(*,*) 'V#Node',VNbrNode%I
#endif
    !
    ! Gather the arrays on Root (Can be done in a better way!) (MUST BE CHANGED!).
    CALL Gather(RowLoc,RowTot,NbrNode,VNbrNode,Disp)
    CALL Gather(ColLoc,ColTot,NbrNode,VNbrNode,Disp)
    CALL Gather(ParLoc,ParTot,NbrNode,VNbrNode,Disp)
    !
    ! Print some stuff if needed.
#ifdef PARTDRV_DBUG
    IF(MYID==0) WRITE(*,*) RowTot%I
#endif
    !
    ! Delete local buffers (MUST BE CHANGED!).
    CALL Delete(RowLoc)
    CALL Delete(ColLoc)
    CALL Delete(ParLoc)
    !
    ! Transform Shity->FastMat (MUST BE CHANGED!)
    IF(MyID.EQ.ROOT) CALL Set_FASTMAT_EQ_SHITY(B,RowTot,ColTot,ParTot,NbrNodeTot)
    !
    ! Delete Root buffers.
    IF(MyID.EQ.ROOT) THEN
       CALL Delete(RowTot)
       CALL Delete(ColTot)
       CALL Delete(ParTot)
    ENDIF
    !
  END SUBROUTINE PDrv_Gather_Part
  !
  !
  SUBROUTINE Set_FASTMAT_EQ_SHITY(A,VRow,VCol,VPar,N)
!H---------------------------------------------------------------------------------
!H
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(FASTMAT ), POINTER    :: A
    TYPE(INT_VECT)             :: VRow,VCol
    TYPE(DBL_VECT)             :: VPar
    INTEGER       , INTENT(IN) :: N
    !-------------------------------------------------------------------
    TYPE(FASTMAT ), POINTER    :: C
    TYPE(SRST    ), POINTER    :: P
    INTEGER                    :: Row,Col,I
    REAL(DOUBLE)               :: Par
    !-------------------------------------------------------------------
    !
    NULLIFY(C,P)
    ! Some checks.
    IF(ASSOCIATED(A))THEN
       CALL Delete_FASTMAT1(A)
    ENDIF
    CALL New_FASTMAT(A,0,(/0,0/))
    !
    ! Copy vector stuff in the FASTMAT.
    DO I = 1,N
       Row = VRow%I(I)
       Col = VCol%I(I)
       Par = VPar%D(I)
       ! Set current row link.
       C => FindFastMatRow_1(A,Row)
       ! Set current col link.
       P => InsertSRSTNode(C%RowRoot,Col)
       ! Put value.
       P%Part = Par
    ENDDO
    !
    ! Build the linked list.
    CALL FlattenAllRows(A)
    !
  END SUBROUTINE Set_FASTMAT_EQ_SHITY
  !
  !
  SUBROUTINE PDrv_Part_1D(A,ParType,PrcDist)
!H---------------------------------------------------------------------------------
!H SUBROUTINE PDrv_Part_1D(A,ParType)
!H
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(FASTMAT)   , POINTER       :: A
    CHARACTER(LEN=*), INTENT(IN   ) :: ParType
    TYPE(INT_VECT)                  :: PrcDist
    !-------------------------------------------------------------------
    INTEGER                         :: IPrc
    TYPE(DBL_VECT)                  :: DeltaT
    !-------------------------------------------------------------------
    !
    ! Create array.
    CALL New(DeltaT,NPrc,0)
    CALL SetEq(DeltaT,Zero)
    !
    ! Compute row distribution.
    CALL Get_RowDist(A,PrcDist,ParType,DeltaT_O=DeltaT)
    !
    ! Simple check.
    ! ToDo....
    !
    ! Set columns partition.
    DO iPrc = 0,NPrc-1
       PrcDist%I(iPrc*4+2)=1
       PrcDist%I(iPrc*4+3)=NAtoms
    ENDDO
    !
    SELECT CASE(TRIM(ParType))
    CASE('Non0')
       !
#ifdef PARTDRV_INFO
       ! Print the new partition.
       DO iPrc = 0,NPrc-1
          WRITE(*,1000) iPrc,INT(DeltaT%D(iPrc)), &
               &        PrcDist%I(iPrc*4  ),PrcDist%I(iPrc*4+1), &
               &        PrcDist%I(iPrc*4+2),PrcDist%I(iPrc*4+3)
       ENDDO
       WRITE(*,1010) INT(SUM(DeltaT%D(0:NPrc-1)))
       WRITE(*,1011) INT(MINVAL(DeltaT%D(0:NPrc-1))),INT(MAXVAL(DeltaT%D(0:NPrc-1)))
       WRITE(*,1012) INT(SUM(DeltaT%D(0:NPrc-1))/DBLE(NPrc)), &
            &        INT(SQRT(DOT_PRODUCT(DeltaT%D(0:NPrc-1)-(SUM(DeltaT%D(0:NPrc-1))/DBLE(NPrc)), &
            &                 DeltaT%D(0:NPrc-1)-(SUM(DeltaT%D(0:NPrc-1))/DBLE(NPrc)))/DBLE(NPrc)))
#endif
       !
    CASE('Time')
       ! Save new partition on the disc.
       CALL Put(PrcDist,TRIM(PartName))
       !CALL Put(PrcDist,'ONXPart')
       !
#ifdef PARTDRV_INFO
       ! Print the new partition.
       DO IPrc = 0,NPrc-1
          WRITE(*,1100) iPrc,DeltaT%D(iPrc), &
               &        PrcDist%I(iPrc*4  ),PrcDist%I(iPrc*4+1), &
               &        PrcDist%I(iPrc*4+2),PrcDist%I(iPrc*4+3)
       ENDDO
       WRITE(*,1110) SUM(DeltaT%D)
       WRITE(*,1111) MINVAL(DeltaT%D(0:NPrc-1)),MAXVAL(DeltaT%D(0:NPrc-1))
       WRITE(*,1112) SUM(DeltaT%D(0:NPrc-1))/DBLE(NPrc), &
            &        SQRT(DOT_PRODUCT(DeltaT%D(0:NPrc-1)-(SUM(DeltaT%D(0:NPrc-1))/DBLE(NPrc)), &
            &             DeltaT%D(0:NPrc-1)-(SUM(DeltaT%D(0:NPrc-1))/DBLE(NPrc)))/DBLE(NPrc))
#endif
       !
    CASE DEFAULT
       CALL Halt(' Unknown Partition in PDrv_Part_1D. ')
    END SELECT
    !
    ! Delete array.
    CALL Delete(DeltaT)
    !
    ! Format declaretions.
1000 FORMAT(' IPrc = ',I3,', Non0 = ',  I10,', NewBeg = ',I4,', NewEnd = ',I4, &
          & ', ColBeg = ',I4,', ColEnd = ',I4)
1010 FORMAT(' SumNon0 = ',  I10)
1011 FORMAT(' MinNon0 = ',  I10,' MaxNon0 = ',  I10)
1012 FORMAT(' AveNon0 = ',  I10,' VarNon0 = ',  I10)
1100 FORMAT(' IPrc = ',I3,', Time = ',E10.4,', RowBeg = ',I4,', RowEnd = ',I4, &
          & ', ColBeg = ',I4,', ColEnd = ',I4)
1110 FORMAT(' SumTime = ',E10.4)
1111 FORMAT(' MinTime = ',E10.4,' MaxTime = ',E10.4)
1112 FORMAT(' AveTime = ',E10.4,' VarTime = ',E10.4)
    !
  END SUBROUTINE PDrv_Part_1D
  !
  !
  SUBROUTINE Get_RowDist(A,PrcDist,ParType,W_O,NumPrc_O,Copy_O,DeltaT_O)
!H---------------------------------------------------------------------------------
!H SUBROUTINE Get_RowDist(A,PrcDist,ParType,W_O,NumPrc_O,Copy_O,DeltaT_O)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(FASTMAT )          , POINTER       :: A
    TYPE(INT_VECT)                          :: PrcDist
    INTEGER       , OPTIONAL, INTENT(IN   ) :: NumPrc_O
    CHARACTER(LEN=*)        , INTENT(IN   ) :: ParType
    TYPE(DBL_VECT), OPTIONAL                :: W_O
    TYPE(INT_VECT), OPTIONAL                :: Copy_O
    TYPE(DBL_VECT), OPTIONAL                :: DeltaT_O
    !-------------------------------------------------------------------
    TYPE(DBL_VECT)                          :: T
    TYPE(DBL_VECT)                          :: W
    TYPE(INT_VECT)                          :: Copy
    REAL(DOUBLE  )                          :: Bs,Wt,OpB
    INTEGER                                 :: p,pT,sp,sp_old,NumPrc,iCopy
    !-------------------------------------------------------------------
    real(double) :: tOpB
    integer :: NbrZero,iLoop
    !
    CALL New(T,NAtoms,0)
    CALL SetEq(T,Zero)
    !
    ! Set number of proc's.
    NumPrc = NPrc
    IF(PRESENT(NumPrc_O)) NumPrc = NumPrc_O
    !
    ! Set weight.
    CALL New(W,NumPrc+1)
    CALL SetEq(W,One)
    IF(PRESENT(W_O)) CALL SetEq(W,W_O)
    !
    ! Set shift.
    CALL New(Copy,NumPrc)
    CALL SetEq(Copy,1)
    IF(PRESENT(Copy_O)) CALL SetEq(Copy,Copy_O)
    !
    ! Sum the number of Non-Zero elements in each row.
    CALL Get_Sum_Row(A,T,ParType)
    !
    ! Do the pre-Ali-partitioning on the rows.
    OpB = Get_Optimal_Bound(T,NumPrc_O)
    !
    !
    !------------------------------------------------>>>
    tOpB = OpB
    iLoop = 1
    DO
    !<<<------------------------------------------------
    !
       Bs = OpB*W%D(1)
       pT = 0
       sp_old = 1
       DO p = 0,NumPrc-1
          sp = BinSrch(T,Bs,NAtoms)
          IF(p.EQ.NumPrc-1) sp = NAtoms
          DO iCopy = 1,Copy%I(p+1)
             PrcDist%I(pT*4  ) = sp_old
             PrcDist%I(pT*4+1) = sp
             pT = pT+1
          ENDDO
          Bs = T%D(sp)+OpB*W%D(p+2)
          sp_old = sp+1
       ENDDO
    !
    !------------------------------------------------>>>
       NbrZero = 0
       DO p = 0,NPrc-1
          IF(PrcDist%I(p*4).GT.PrcDist%I(p*4+1)) NbrZero = NbrZero+1
       ENDDO
       IF(NbrZero.EQ.0) EXIT
#ifdef PARTDRV_INFO
       WRITE(*,*) 'There is(are)',NbrZero,' zero(s).'
#endif
       !
       ! Reduce the optimal bottleneck.
       OpB = tOpB*(One-DBLE(iLoop)*0.0025D+00)
       iLoop = iLoop+1
       !
       ! Exit if to many loops.
       IF(iLoop.GT.20) EXIT
    ENDDO
    !<<<------------------------------------------------
    !
    !
    ! If needed copy.
    IF(PRESENT(DeltaT_O)) THEN
       DO p = 0,NumPrc-1
          DeltaT_O%D(p) = T%D(PrcDist%I(p*4+1))-T%D(PrcDist%I(p*4)-1)
       ENDDO
    ENDIF
    !
#ifdef PARTDRV_DBUG
    ! Print if needed.
    WRITE(*,1000) OpB
    pT = 0
    DO p = 0,NumPrc-1
       WRITE(*,1100) p,T%D(PrcDist%I(pT*4+1))-T%D(PrcDist%I(pT*4)-1), &
            &        PrcDist%I(pT*4),PrcDist%I(pT*4+1)
       pT = pT+Copy%I(p+1)
    ENDDO
    WRITE(*,*) ''
    !
    ! Format declaretions.
1000 FORMAT(' OpB  = ',E12.6)
1100 FORMAT(' Pat = ',I3,', Time = ',E12.6,', New Beg = ',I4,', New End = ',I4)
#endif
    !
    ! Delete local array.
    CALL Delete(T)
    CALL Delete(W)
    CALL Delete(Copy)
    !
  END SUBROUTINE Get_RowDist
  !
  !
  SUBROUTINE Get_ColDist(A,ParType,NbrPrt,ColP,PrcDist,NumPrc_O,DeltaT_O)
!H---------------------------------------------------------------------------------
!H SUBROUTINE Get_ColDist(A,ParType,NbrPrt,ColP,PrcDist,NumPrc_O,DeltaT_O)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(FASTMAT ), POINTER                 :: A
    TYPE(INT_VECT)                          :: ColP
    TYPE(INT_VECT)                          :: PrcDist
    TYPE(DBL_VECT)               , OPTIONAL :: DeltaT_O
    CHARACTER(LEN=*),INTENT(IN  )           :: ParType
    INTEGER       , INTENT(IN   )           :: NbrPrt
    INTEGER       , INTENT(IN   ), OPTIONAL :: NumPrc_O
    !-------------------------------------------------------------------
    TYPE(FASTMAT ), POINTER                 :: P
    TYPE(SRST    ), POINTER                 :: U
    TYPE(DBL_VECT)                          :: T
    INTEGER                                 :: iPrt,Col,Row,RowMin
    INTEGER                                 :: RowMax,M,N,NumPrc,AtA
    INTEGER                                 :: iPrc,sp,sp_old,iPrcT,iPrcT1
    INTEGER                                 :: IType
    REAL(DOUBLE)                            :: OpB,WMin,Bs
    INTEGER       , PARAMETER               :: TYPE_NON0 = 1
    INTEGER       , PARAMETER               :: TYPE_TIME = 2
    !-------------------------------------------------------------------
    !
    NULLIFY(P,U)
    ! Set number of proc's.
    NumPrc = NPrc
    IF(PRESENT(NumPrc_O)) NumPrc = NumPrc_O
    !
    ! Select partition type.
    SELECT CASE(TRIM(ParType))
    CASE('Non0'); IType = TYPE_NON0
    CASE('Time'); IType = TYPE_TIME
    CASE DEFAULT; STOP
    END SELECT
    !
    ! Allocate local array.
    CALL New(T,NAtoms,0)
    !
    ! Set some variables.
    iPrcT = 0
    iPrcT1 = 0
    P => A%Next
    DO iPrt = 1,NbrPrt
       !
       ! Get the min and max row for the give partition.
       RowMin = PrcDist%I(iPrcT1*4)
       RowMax = PrcDist%I(iPrcT1*4+1)
#ifdef PARTDRV_DBUG
       WRITE(*,*) 'RowMin',RowMin,'RowMax',RowMax
#endif
       !
       ! Initialize.
       CALL SetEq(T,Zero)
       !
       ! Let's go.
       DO
          IF(.NOT.ASSOCIATED(P)) EXIT
          Row = P%Row
          IF(Row.LT.RowMin) STOP 'Problem in Get_ColDist.'
          IF(Row.GT.RowMax) EXIT
          M = BSiz%I(Row)
          U => P%RowRoot
          DO
             IF(.NOT.ASSOCIATED(U)) EXIT
             IF(U%L.EQ.U%R) THEN
                Col = U%L
                N = BSiz%I(Col)
                !
                ! Select the wanted partition.
                SELECT CASE(IType)
                CASE(TYPE_NON0); T%D(Col) = T%D(Col)+DBLE(M*N)
                CASE(TYPE_TIME); T%D(Col) = T%D(Col)+U%Part
                END SELECT
             ENDIF
             U => U%Next
          ENDDO
          P => P%Next
       ENDDO
       !
       ! Max and Min elements per row.
       ! WMaxI = MAXVAL(TI%I)
       WMin = MINVAL(T%D)
       !
       ! Simple check.
       IF(WMin.LT.0) &
            & CALL Halt(' A Non0 element in ..... is smaller that 0! ')
#ifdef PARTDRV_DBUG
       WRITE(*,*) 'SUM IS:',SUM(T%D)
#endif
       !
       ! Prefix sum.
       DO AtA = 1,NAtoms
          T%D(AtA) = T%D(AtA)+T%D(AtA-1)
       ENDDO
       !
       ! Do the pre-Ali-partitioning on the rows.
       OpB = Get_Optimal_Bound(T,ColP%I(iPrt))
       !
       ! Find the partition.
       sp_old = 1
       Bs = OpB
       DO iPrc = 1,ColP%I(iPrt)
          sp = BinSrch(T,Bs,NAtoms)
          IF(iPrc.EQ.ColP%I(iPrt)) sp = NAtoms
          PrcDist%I(iPrcT1*4+2) = sp_old
          PrcDist%I(iPrcT1*4+3) = sp
          Bs = T%D(sp)+OpB
          sp_old = sp+1
          iPrcT1=iPrcT1+1
       ENDDO
#ifdef PARTDRV_DBUG
       WRITE(*,*) 'PrcDist',PrcDist%I
       WRITE(*,*) 'OpB',OpB,iPrt
#endif
       !
       ! If needed copy.
       IF(PRESENT(DeltaT_O)) THEN
          DO iPrc = 1,ColP%I(iPrt)
             DeltaT_O%D(iPrcT) = T%D(PrcDist%I(iPrcT*4+3))-T%D(PrcDist%I(iPrcT*4+2)-1)
             iPrcT = iPrcT+1
          ENDDO
       ENDIF
       !
    ENDDO
    !
    ! Delete the local array.
    CALL Delete(T)
    !
  END SUBROUTINE Get_ColDist
  !
  !
  SUBROUTINE Get_RowColDist(A,PrcDist,ParType,W_O,NumPrc_O,DeltaT_O)
!H---------------------------------------------------------------------------------
!H SUBROUTINE Get_RowColDist(A,PrcDist,ParType,W_O,NumPrc_O,DeltaT_O)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(FASTMAT )          , POINTER       :: A
    TYPE(INT_VECT)                          :: PrcDist
    TYPE(DBL_VECT), OPTIONAL                :: W_O
    TYPE(DBL_VECT), OPTIONAL                :: DeltaT_O
    INTEGER       , OPTIONAL, INTENT(IN   ) :: NumPrc_O
    CHARACTER(LEN=*)        , INTENT(IN   ) :: ParType
    !-------------------------------------------------------------------
    TYPE(INT_VECT)               :: ColP
    TYPE(DBL_VECT)               :: Weight
    INTEGER                      :: RowP,MaxCol,ij,i,j
    !-------------------------------------------------------------------
    !
    ! Get the # of proc's per row and the weights.
    CALL Get_DistPrc_2D(RowP,ColP,Weight)
    !
    ! Get the # max of proc per row.
    MaxCol = MAXVAL(ColP%I)
    !
    CALL SetEq(PrcDist,0)
    !
    ! Get the partial row distribution.
    CALL Get_RowDist(A,PrcDist,ParType,NumPrc_O=RowP,Copy_O=ColP,W_O=Weight)
    !
    ! Get the full row-col partition.
    CALL Get_ColDist(A,ParType,RowP,ColP,PrcDist,DeltaT_O=DeltaT_O)
    !
#ifdef PARTDRV_DBUG
    DO I = 1,RowP
       WRITE(*,*) 'RowP',RowP,'ColP',ColP%I(i)
    ENDDO
#endif
    !
    ! Delete local array.
    CALL Delete(ColP  )
    CALL Delete(Weight)
    !
  END SUBROUTINE Get_RowColDist
  !
  !
  SUBROUTINE PDrv_Part_2D(A,ParType,PrcDist)
!H---------------------------------------------------------------------------------
!H SUBROUTINE PDrv_Part_2D(A,ParType,PrcDist)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(FASTMAT)   , POINTER       :: A
    CHARACTER(LEN=*), INTENT(IN   ) :: ParType
    TYPE(INT_VECT)                  :: PrcDist
    !-------------------------------------------------------------------
    INTEGER                      :: Row,Col,M,N,AtA,IPrc,OpBI
    INTEGER                      :: WMaxI,WMinI,WtI,LBI,UBI,midBI,BI
    TYPE(DBL_VECT)               :: DeltaT
    !-------------------------------------------------------------------
    !
    ! Initialize array.
    CALL New(DeltaT,NPrc,0)
    CALL SetEq(DeltaT,Zero)
    !
    ! Get the distribution.
    CALL Get_RowColDist(A,PrcDist,ParType,DeltaT_O=DeltaT)
    !
    ! Simple check.
    ! ToDo....
    !

    !
    ! Select the desired parition.
    SELECT CASE(ParType)
    CASE('Non0')
       !
#ifdef PARTDRV_INFO
       ! Print if needed.
       DO iPrc = 0,NPrc-1
          WRITE(*,1000) iPrc,INT(DeltaT%D(iPrc)), &
               &        PrcDist%I(iPrc*4  ),PrcDist%I(iPrc*4+1), &
               &        PrcDist%I(iPrc*4+2),PrcDist%I(iPrc*4+3)
       ENDDO
       WRITE(*,1010) INT(SUM(DeltaT%D(0:NPrc-1)))
       WRITE(*,1011) INT(MINVAL(DeltaT%D(0:NPrc-1))),INT(MAXVAL(DeltaT%D(0:NPrc-1)))
       WRITE(*,1012) INT(SUM(DeltaT%D(0:NPrc-1))/DBLE(NPrc)), &
            &        INT(SQRT(DOT_PRODUCT(DeltaT%D(0:NPrc-1)-(SUM(DeltaT%D(0:NPrc-1))/DBLE(NPrc)), &
            &                 DeltaT%D(0:NPrc-1)-(SUM(DeltaT%D(0:NPrc-1))/DBLE(NPrc)))/DBLE(NPrc)))
#endif
       !
    CASE('Time')
       !
       ! Save new partition on the disc.
       CALL Put(PrcDist,TRIM(PartName))
       !CALL Put(PrcDist,'ONXPart')
       !
#ifdef PARTDRV_INFO
       ! Print if needed.
       DO IPrc = 0,NPrc-1
          WRITE(*,1100) iPrc,DeltaT%D(iPrc), &
               &        PrcDist%I(iPrc*4  ),PrcDist%I(iPrc*4+1), &
               &        PrcDist%I(iPrc*4+2),PrcDist%I(iPrc*4+3)
       ENDDO
       WRITE(*,1110) SUM(DeltaT%D(0:NPrc-1))
       WRITE(*,1111) MINVAL(DeltaT%D(0:NPrc-1)),MAXVAL(DeltaT%D(0:NPrc-1))
       WRITE(*,1112) SUM(DeltaT%D(0:NPrc-1))/DBLE(NPrc), &
            &        SQRT(DOT_PRODUCT(DeltaT%D(0:NPrc-1)-(SUM(DeltaT%D(0:NPrc-1))/DBLE(NPrc)), &
            &             DeltaT%D(0:NPrc-1)-(SUM(DeltaT%D(0:NPrc-1))/DBLE(NPrc)))/DBLE(NPrc))
#endif
       !
    CASE DEFAULT
       CALL Halt(' Unknown Partition in PDrv_Part_2D. ')
    END SELECT
    !
    ! Delete array.
    CALL Delete(DeltaT)
    !
    ! Format declarations.
1000 FORMAT(' IPrc = ',I3,', Non0 = ',  I10,', RowBeg = ',I4,', RowEnd = ',I4, &
          & ', ColBeg = ',I4,', ColEnd = ',I4)
1010 FORMAT(' SumNon0 = ',  I10)
1011 FORMAT(' MinNon0 = ',  I10,' MaxNon0 = ',  I10)
1012 FORMAT(' AveNon0 = ',  I10,' SdvNon0 = ',  I10)
1100 FORMAT(' IPrc = ',I3,', Time = ',E10.4,', RowBeg = ',I4,', RowEnd = ',I4, &
          & ', ColBeg = ',I4,', ColEnd = ',I4)
1110 FORMAT(' SumTime = ',E10.4)
1111 FORMAT(' MinTime = ',E10.4,' MaxTime = ',E10.4)
1112 FORMAT(' AveTime = ',E10.4,' SdvTime = ',E10.4)
    !
  END SUBROUTINE PDrv_Part_2D
  !
  !
  SUBROUTINE Get_Sum_Row(A,T,ParType)
!H---------------------------------------------------------------------------------
!H SUBROUTINE Get_Sum_Row(A,T,ParType)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(FASTMAT )  , POINTER     :: A
    TYPE(DBL_VECT)                :: T
    CHARACTER(LEN=*), INTENT(IN ) :: ParType
    !-------------------------------------------------------------------
    TYPE(FASTMAT )  , POINTER     :: P
    TYPE(SRST    )  , POINTER     :: U
    INTEGER                       :: Row,Col,AtA,M,N,IType
    REAL(DOUBLE)                  :: WMax,WMin,Wt,LB,UB
    INTEGER         , PARAMETER   :: TYPE_NON0=1
    INTEGER         , PARAMETER   :: TYPE_TIME=2
    !-------------------------------------------------------------------
    !
    NULLIFY(P,U)
    !
    SELECT CASE(TRIM(ParType))
    CASE('Non0'); IType=TYPE_NON0
    CASE('Time'); IType=TYPE_TIME
    CASE DEFAULT; STOP
    END SELECT
    !
    ! Create local array.
    ! do some check here
    CALL SetEq(T,Zero)
    !
    ! Sum the time in each row.
    P => A%Next
    DO
       IF(.NOT.ASSOCIATED(P)) EXIT
       Row = P%Row
       M = BSiz%I(Row)
       U => P%RowRoot
       DO
          IF(.NOT.ASSOCIATED(U)) EXIT
          IF(U%L.EQ.U%R) THEN
             Col = U%L
             N = BSiz%I(Col)
             SELECT CASE(IType)
             CASE(TYPE_NON0); T%D(Row) = T%D(Row)+DBLE(M*N)
             CASE(TYPE_TIME); T%D(Row) = T%D(Row)+U%Part
             END SELECT
          ENDIF
          U => U%Next
       ENDDO
       P => P%Next
    ENDDO
    !
    ! Max elements per row.
    WMax = MAXVAL(T%D)
    WMin = MINVAL(T%D)
    !
    ! Simple check.
    IF(WMin.LT.Zero) &
         & CALL Halt(' The Non0 or Time in Get_Sum_Row is smaller that zero! ')
    !
    ! Prefix sum.
    DO AtA = 1,NAtoms
       T%D(AtA) = T%D(AtA)+T%D(AtA-1)
    ENDDO
    !
  END SUBROUTINE Get_Sum_Row
  !
  !
  FUNCTION Get_Optimal_Bound(T,NumPrc_O) RESULT(OpB)
!H---------------------------------------------------------------------------------
!H FUNCTION Get_Optimal_Bound(T,NumPrc_O) RESULT(OpB)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(DBL_VECT)                :: T
    INTEGER, OPTIONAL, INTENT(IN) :: NumPrc_O
    REAL(DOUBLE)                  :: OpB
    !-------------------------------------------------------------------
    INTEGER                       :: NumPrc
    REAL(DOUBLE)                  :: midB,UB,LB,Wt,WMax,Tol
    !-------------------------------------------------------------------
    !
    ! Set number of proc's
    NumPrc = NPrc
    IF(PRESENT(NumPrc_O)) NumPrc = NumPrc_O
    !
    ! Max elements per row.
    WMax = MAXVAL(T%D)
    !
    ! Set some variables.
    Wt = T%D(NAtoms)
    LB = Wt/DBLE(NumPrc)!DBLE(NPrc)
    UB = LB+WMax
    !
    ! Set Tol.
    ! This is a test, Tol seems to need to be about 1% of LBD.
    Tol = 0.01D+00*LB
    !
    ! Do the pre-Ali-partitioning.
    DO
       midB = (UB+LB)/Two
       IF(Probe(T,midB,Wt,NumPrc)) THEN
          UB = midB
       ELSE
          LB = midB+One
       ENDIF
       IF(UB-LB.LE.Tol) EXIT
    ENDDO
    !
    ! Now UB is optimal.
    OpB = UB
    !
  END FUNCTION Get_Optimal_Bound
  !
  !
  SUBROUTINE Get_DistPrc_2D(RowP,ColP,Weight)
!H---------------------------------------------------------------------------------
!H SUBROUTINE Get_Dist_Prc_2D(RowP,ColP,Weight)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    INTEGER       , INTENT(OUT) :: RowP
    TYPE(INT_VECT)              :: ColP
    TYPE(DBL_VECT)              :: Weight
    !-------------------------------------------------------------------
    INTEGER                     :: ColTot,iRow
    !-------------------------------------------------------------------
    !
    ! Get the number of Row.
    RowP = CEILING(SQRT(DBLE(NPrc)))
    !
    ! Allocate some some arrays.
    IF(AllocQ(ColP%Alloc  )) CALL Delete(ColP  )
    IF(AllocQ(Weight%Alloc)) CALL Delete(Weight)
    CALL New(ColP  ,RowP  )
    CALL New(Weight,RowP+1)
    !
    CALL SetEq(ColP,0)
    CALL SetEq(Weight,Zero)
    !
    ! Set the the number of proc's for the first row.
    ColP%I(1)   = NINT(SQRT(DBLE(NPrc)))
    !
    ! Set the first weight.
    Weight%D(1) = DBLE(ColP%I(1)*RowP)/DBLE(NPrc)
    !
    ! Set the partial number of proc's.
    ColTot = ColP%I(1)
    DO iRow = 2,RowP
       IF(MOD(iRow,2).EQ.0) THEN
          ColP%I(iRow) = FLOOR(  DBLE(NPrc-ColTot)/DBLE(RowP+1-iRow))
       ELSE
          ColP%I(iRow) = CEILING(DBLE(NPrc-ColTot)/DBLE(RowP+1-iRow))
       ENDIF
       Weight%D(iRow) = DBLE(ColP%I(iRow)*RowP)/DBLE(NPrc)
       ColTot = ColTot+ColP%I(iRow)
    ENDDO
    !
  END SUBROUTINE Get_DistPrc_2D
  !
  !
  PURE FUNCTION BinSrch(DVec,DVal,NDim) RESULT(Idx)
!H---------------------------------------------------------------------------------
!H PURE FUNCTION BinSrch(DVec,DVal,NDim) RESULT(Idx)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    REAL(DOUBLE)  , INTENT(IN) :: DVal
    INTEGER       , INTENT(IN) :: NDim
    TYPE(DBL_VECT), INTENT(IN) :: DVec
    !-------------------------------------------------------------------
    INTEGER                    :: Idx,IMin,IMax,IMid
    LOGICAL                    :: IsEqual
    !-------------------------------------------------------------------
    IMin = 0
    IMax = NDim+1
    IsEqual = .FALSE.
    DO WHILE(IMax-IMin.GT.1)
       IMid = (IMax+IMin)/2
       IF(DVal.GT.DVec%D(IMid)) THEN
          IMin = IMid
       ELSEIF(DVal.LT.DVec%D(IMid)) THEN
          IMax = IMid
       ELSE
          IsEqual = .TRUE.
          EXIT
       ENDIF
    ENDDO
    !
    IF(IMax.GT.NDim) IMax = NDim
    IF(IMin.LT.   1) IMin = 1
    !
    IF(IsEqual) THEN
       Idx = IMid
    ELSEIF(ABS(DVec%D(IMax)-DVal).LT.ABS(DVec%D(IMin)-DVal)) THEN
       Idx = IMax
    ELSE
       Idx = IMin
    ENDIF
  END FUNCTION BinSrch
  !
  !
  LOGICAL FUNCTION Probe(T,B,Wt,NumPrc)
!H---------------------------------------------------------------------------------
!H LOGICAL FUNCTION Probe(T,B,Wt,NumPrc)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(DBL_VECT)             :: T
    REAL(DOUBLE  ), INTENT(IN) :: B,Wt
    INTEGER       , INTENT(IN) :: NumPrc
    !-------------------------------------------------------------------
    REAL(DOUBLE  )             :: Bs
    INTEGER                    :: p,sp
    !-------------------------------------------------------------------
    Bs = B
    p = 1
    DO WHILE(p.LE.NumPrc.AND.Bs.LT.Wt)   !vw changed that.
       sp = BinSrch(T,Bs,NAtoms)
       Bs = T%D(sp)+B
       p = p+1
    ENDDO
    IF(Bs.LT.Wt) THEN
       Probe = .FALSE.
    ELSE
       Probe = .TRUE.
    ENDIF
  END FUNCTION Probe
  !
  !
  SUBROUTINE Set_Universal_BegEnd(PrcDist)
!H---------------------------------------------------------------------------------
!H SUBROUTINE Set_Universal_BegEnd(PrcDist)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(INT_VECT)             :: PrcDist
    !-------------------------------------------------------------------
    INTEGER                    :: iPrc
    !-------------------------------------------------------------------
    IF(MyID.EQ.ROOT) THEN
#ifdef PARTDRV_DBUG
       WRITE(*,*) 'Before change: Beg%I',Beg%I
       WRITE(*,*) 'Before change: End%I',End%I
#endif
       DO iPrc = 0,NPrc-1
          Beg%I(iPrc) = PrcDist%I(iPrc*4  )
          End%I(iPrc) = PrcDist%I(iPrc*4+1)
       ENDDO
#ifdef PARTDRV_DBUG
       WRITE(*,*) 'After change: Beg%I',Beg%I
       WRITE(*,*) 'After change: End%I',End%I
#endif
    ENDIF
  END SUBROUTINE Set_Universal_BegEnd
  !
  !
  SUBROUTINE Save_Universal_BegEnd()
!H---------------------------------------------------------------------------------
!H SUBROUTINE Save_Universal_BegEnd()
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    INTEGER :: iPrc
    !-------------------------------------------------------------------
    ! Save the current Beg and End.
    IF(MyID.EQ.ROOT) THEN
       IF(AllocQ(Beg_Save%Alloc)) CALL Delete(Beg_Save)
       IF(AllocQ(End_Save%Alloc)) CALL Delete(End_Save)
       CALL New(Beg_Save,NPrc-1,0)
       CALL New(End_Save,NPrc-1,0)
       DO iPrc = 0,NPrc-1
          Beg_Save%I(iPrc) = Beg%I(iPrc)
          End_Save%I(iPrc) = End%I(iPrc)
       ENDDO
    ENDIF
  END SUBROUTINE Save_Universal_BegEnd
  !
  !
  SUBROUTINE Load_Universal_BegEnd()
!H---------------------------------------------------------------------------------
!H SUBROUTINE Load_Universal_BegEnd()
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    INTEGER :: iPrc
    !-------------------------------------------------------------------
    IF(MyID.EQ.ROOT) THEN
       IF(.NOT.AllocQ(Beg_Save%Alloc)) CALL Delete(Beg_Save)
       IF(.NOT.AllocQ(End_Save%Alloc)) CALL Delete(End_Save)
       DO iPrc = 0,NPrc-1
          Beg%I(iPrc) = Beg_Save%I(iPrc)
          End%I(iPrc) = End_Save%I(iPrc)
       ENDDO
       CALL Delete(Beg_Save)
       CALL Delete(End_Save)
    ENDIF
  END SUBROUTINE Load_Universal_BegEnd
  !
  !
#endif
  !
#ifdef BLABLA1000
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
    INTEGER                                 :: I,PU,J,M,N,iErr
    REAL(DOUBLE)                            :: Chk
    LOGICAL                                 :: InPara
    CHARACTER(LEN=DEFAULT_CHR_LEN)       :: ChkStr
    REAL(DOUBLE), EXTERNAL               :: DDOT
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
       if(myid==9)write(*,*) 'Row=',I
       M = BSiz%I(I)
       U => R%RowRoot
       DO
          IF(.NOT.ASSOCIATED(U)) EXIT
          IF(U%L.EQ.U%R) THEN
             IF(ASSOCIATED(U%MTrix)) THEN
                J = U%L
                if(myid==9)write(*,*) 'Col=',J
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
    DO I=0,NPrc
       IF(MyID.EQ.I) then
          WRITE(*,'(1X,A,1X,I4)')TRIM(ChkStr),I
       ENDIF
       DO J=1,10;CALL MPI_BARRIER(MONDO_COMM,iErr);ENDDO
    ENDDO
  END SUBROUTINE PChkSum_FASTMAT2
  SUBROUTINE CheckSum_DBCSR2(A,Name,Proc_O,Unit_O)
    TYPE(DBCSR)               :: A
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proc_O
    INTEGER,         OPTIONAL,INTENT(IN) :: Unit_O
    REAL(DOUBLE)              :: Chk
    REAL(DOUBLE), EXTERNAL    :: DDot
    REAL(DOUBLE)              :: DotPrd
    CHARACTER(LEN=*)          :: Name
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: ChkStr
    INTEGER :: I,PU,iErr,J
    Chk=Zero
    DO I=1,A%NNon0
       Chk=Chk+A%MTrix%D(I)*A%Mtrix%D(I)
    ENDDO
    Chk=SQRT(Chk)
    ChkStr=CheckSumString(Chk,Name,Proc_O)
    DO I=0,NPrc
       IF(MyID.EQ.I) then
          WRITE(*,'(1X,A,1X,I4)')TRIM(ChkStr),I
       ENDIF
       DO J=1,10;CALL MPI_BARRIER(MONDO_COMM,iErr);ENDDO
    ENDDO
  END SUBROUTINE CheckSum_DBCSR2
#endif
END MODULE PartDrv
