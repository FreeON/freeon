MODULE ParseParallel
  USE InOut
#ifdef NAG
   USE F90_UNIX
#endif
  USE ParallelKeys
  USE ControlStructures
CONTAINS
  !============================================================================
  !
  !============================================================================
  SUBROUTINE LoadParallel(N,O,G,B,M)
    TYPE(FileNames)    :: N
    TYPE(Options)      :: O
    TYPE(Geometries)   :: G
    TYPE(BasisSets)    :: B
    TYPE(Parallel)     :: M
    INTEGER            :: I,J
    CHARACTER(LEN=DCL) :: MFile
#ifdef PARALLEL
    !-------------------------------------------------------------------------!  
    CALL OpenASCII(N%IFile,Inp)
#if !defined(MPI2)
    ! We are doing MPI-1, where mpirun is spawned underneath MondoSCF
    ! Note that this does not work (so far anyway) on AIX with POE
    !
    ! Obtain MPI invocation: mpirun, mpiexe, prun, etc.
    IF(.NOT.OptCharQ(Inp,MPI_INVOCATION,M%Invoking))  & 
         CALL MondoHalt(PRSE_ERROR,' MPI invocation not found. Check for option ' &
         //TRIM(MPI_INVOCATION))
    ! Obtain flag for the number of processors; -n, -np, -proc etc
    IF(.NOT.OptCharQ(Inp,MPI_PROCESSOR_FLAG,M%ProcFlag))M%ProcFlag="-np"
    ! Obtain the total number of processors to employ (hard uper limit)
    IF(.NOT.OptIntQ(Inp,MPI_PROCESSOR_NUMBER,M%NProc))THEN
       M%NProc=2
       CALL Warn(' Parallel code defaulting to 2 processors ')
    ENDIF
    ! Obtain the flag for the machine file; -machinefile, etc 
    IF(.NOT.OptCharQ(Inp,MPI_MACHINE_FLAG,M%MachFlag))M%MachFlag=" "
    ! Obtain the machine file  
    IF(.NOT.OptCharQ(Inp,MPI_MACHINE_FILE,M%MachFile))M%MachFile=" "
    ! Can be an env variable like $PBSNODES or $MONDO_NODE_FILE
    IF(SCAN(M%MachFile,'$')/=0)THEN
       MFile=M%MachFile
       CALL GETENV(MFile,M%MachFile)
    ENDIF
#else
    ! We are doing MPI-2, where MondoSCF must be invoked with mpirun.
    ! Note that this does not work (so far anyway) on True64 with LSF, since
    ! COMPAQ MPI is just MPICH.  
    M%NProc=MSize()    
#endif
    ! Parse for the number of processors in the spatial dimension
    IF(.NOT.OptIntQ(Inp,MPI_SPATIAL_PROC,M%NSpace))THEN
       M%NSpace=M%NProc
       CALL Warn(' # of spatial proc not specified, defaulting to NSpace=NProc='//TRIM(IntToChar(M%NSpace)))
    ENDIF
    CLOSE(UNIT=Inp,STATUS='KEEP')    
    ! Set up parallelism in space and time    
    CALL SpaceTimeSetUp(M%NProc,M%NSpace,G%Clones,M%Clumps,M%Clump)
    ! Space for parallel sparse matrix indecies
    ALLOCATE(M%Beg(1:G%Clones,1:B%NBSets))
    ALLOCATE(M%End(1:G%Clones,1:B%NBSets))
    ALLOCATE(M%GLO(1:G%Clones,1:B%NBSets))
    DO I=1,G%Clones
       DO J=1,B%NBSets
          CALL GreedyDBCSRPartition(B%BSets(I,J)%NAtms,M%NSpace,B%BSiz(I,J),M%Beg(I,J),M%End(I,J),M%GLO(I,J))
#ifdef FULL_ON_FRONT_END_DEBUG
          WRITE(*,*)I,J,' Beg = ',M%Beg(I,J)%I
          WRITE(*,*)I,J,' End = ',M%End(I,J)%I
#endif
       ENDDO
    ENDDO
    ! Just punt for now on setting DBCSR limits 
    M%MxAtsNode=MIN(B%MxAts,CEILING(Two*DBLE(B%MxAts)/DBLE(M%NSpace)))
    M%MxBlkNode=MIN(B%MxBlk,CEILING(Two*DBLE(B%MxBlk)/DBLE(M%NSpace)))
    M%MxN0sNode=MIN(B%MxN0s,CEILING(Two*DBLE(B%MxN0s)/DBLE(M%NSpace)))
#else
    CALL SpaceTimeSetUp(1,1,G%Clones,M%Clumps,M%Clump)
#endif
  END SUBROUTINE LoadParallel
  !============================================================================
  ! BREAK PARALLELSIM IN SPACE AND TIME INTO GROUPS OF CLONES 
  !============================================================================
  SUBROUTINE SpaceTimeSetUp(NProc,NSPace,Clones,Clumps,Clump)
    TYPE(INT_RNK2) :: Clump
    INTEGER        :: ClonesPerClump,LeftOvers,NProc,NSpace,Clones,Clumps,iCLUMP
    !--------------------------------------------------------------------------!
    ! Number of clones in a group (clump)
    ClonesPerClump=NProc/NSpace
    ! Number of groups
    Clumps=Clones/ClonesPerClump
    ! The remainder of NSpace proc jobs that dont fit evenly into a group
    LeftOvers=MOD(Clones*NSpace,NProc)/NSpace
    ! Include leftovers as a reduced clump
    Clumps=Clumps+LeftOvers
    ! Array to hold the grouping info
    CALL New(Clump,(/3,Clumps/))
    DO iCLUMP=1,Clumps
       ! Number of clones in a clump
       Clump%I(1,iCLUMP)=ClonesPerClump
       ! Clone offset of this clump
       Clump%I(2,iCLUMP)=(iCLUMP-1)*ClonesPerClump+1
       ! Number of proc to use per clump
       Clump%I(3,iCLUMP)=NSpace*ClonesPerClump
    ENDDO
    IF(LeftOvers>0)THEN
       ! And here we add in the leftovers
       Clump%I(1,Clumps)=LeftOvers
       Clump%I(2,Clumps)=Clumps*ClonesPerClump+1
       Clump%I(3,Clumps)=LeftOvers*NSpace
    ENDIF
!#ifdef FULL_ON_FRONT_END_DEBUG
    WRITE(*,*)' Clumps         = ',Clumps
    WRITE(*,*)' ClonesPerClump = ',ClonesPerClump
    WRITE(*,*)' LeftOvers      = ',LeftOvers
    DO iCLUMP=1,Clumps
       WRITE(*,*)'CLUMP(',iCLUMP,')= ',Clump%I(:,iCLUMP)
    ENDDO
!#endif
  END SUBROUTINE SpaceTimeSetUp
  !============================================================================
  ! GREEDY LOOK AHEAD DOMAIN DECOMPOSITION TO PARTITION DBCSR MATRICES
  !============================================================================      
  SUBROUTINE GreedyDBCSRPartition(Atoms,NSpace,BlokSize,RowBeg,RowEnd,GLOffSet)
    TYPE(INT_VECT)                      :: BlokSize,RowBeg,RowEnd,GLOffSet
    INTEGER                             :: Atoms,NSpace,NAtsAv,I,K,IPrc,ISet
    INTEGER, PARAMETER                  :: Mns=1,Pls=2
    INTEGER, DIMENSION(Mns:Pls)         :: Beg3,End2
    INTEGER, DIMENSION(Mns:Pls,Mns:Pls) :: End3,Dv
    !-------------------------------------------------------------------------!
    ! Allocate domain limits
    ! Greedy, look ahead algorithm for decomposition
    CALL New(RowBeg,NSpace-1,0)
    CALL New(RowEnd,NSpace-1,0)
    CALL New(GLOffSet,NSpace-1,0)
    IF(Atoms<NSpace)THEN
       ! Silly hack to allow more processors than atoms ...
       RowBeg%I(:)=0
       RowEnd%I(:)=0
       DO IPrc=0,Atoms-1
          RowBeg%I(IPrc)=IPrc+1
          RowEnd%I(IPrc)=IPrc+1
       ENDDO
    ELSE
       RowBeg%I(0)=1
       RowEnd%I(NSpace-1)=Atoms
       DO IPrc=0,NSpace-2
          ! Compute running average to section
          NAtsAv=DBLE(SUM(BlokSize%I(RowBeg%I(IPrc):Atoms)))/DBLE(NSpace-IPrc)
          ! Forcast RowBeg and RowEnd for (de/inc)rements of (+/- 1) 
          DO K=RowBeg%I(IPrc),Atoms
             IF(SUM(BlokSize%I(RowBeg%I(IPrc):K))>=NAtsAv)THEN
                End2(Pls)=K
                EXIT            
             ENDIF
          ENDDO
          End2(Mns)=End2(Pls)-1
          Beg3(Pls)=End2(Pls)+1
          Beg3(Mns)=End2(Mns)+1
          End3(Pls,Pls)=Atoms            
          DO K=Beg3(Pls),Atoms   
             IF(SUM(BlokSize%I(Beg3(Pls):K))>=NAtsAv)THEN
                End3(Pls,Pls)=K
                EXIT            
             ENDIF
          ENDDO
          End3(Pls,Mns)=End3(Pls,Pls)-1
          DO K=Beg3(Mns),Atoms    
             IF(SUM(BlokSize%I(Beg3(Mns):K))>=NAtsAv)THEN
                End3(Mns,Pls)=K
                EXIT            
             ENDIF
          ENDDO
          End3(Mns,Mns)=End3(Mns,Pls)-1
          ! These are deviations from the running average for each choice
          ! of a place to section and the possilbe sectioning in the next iteration
          Dv(Pls,Pls)= &
               ABS(NAtsAv-SUM(BlokSize%I(RowBeg%I(IPrc):End2(Pls))))   &
              +ABS(NAtsAv-SUM(BlokSize%I(Beg3(Pls):End3(Pls,Pls)))) 
          Dv(Pls,Mns)= &
               ABS(NAtsAv-SUM(BlokSize%I(RowBeg%I(IPrc):End2(Pls))))   &
              +ABS(NAtsAv-SUM(BlokSize%I(Beg3(Pls):End3(Pls,Mns)))) 
          Dv(Mns,Pls)= &
               ABS(NAtsAv-SUM(BlokSize%I(RowBeg%I(IPrc):End2(Mns))))   &
              +ABS(NAtsAv-SUM(BlokSize%I(Beg3(Mns):End3(Mns,Pls)))) 
          Dv(Mns,Mns)= &
               ABS(NAtsAv-SUM(BlokSize%I(RowBeg%I(IPrc):End2(Mns))))   &
              +ABS(NAtsAv-SUM(BlokSize%I(Beg3(Mns):End3(Mns,Mns)))) 
          ! Pick the best section based on minimizing the I and I+1 deviation from the average
          IF(MIN(Dv(Pls,Pls),Dv(Pls,Mns))<MIN(Dv(Mns,Pls),Dv(Mns,Mns)))THEN
             RowEnd%I(IPrc)=End2(Pls)
             RowBeg%I(IPrc+1)=RowEnd%I(IPrc)+1
          ELSE
             RowEnd%I(IPrc)=MAX(RowBeg%I(IPrc),End2(Mns))
             RowBeg%I(IPrc+1)=RowEnd%I(IPrc)+1
          ENDIF
       ENDDO
    ENDIF
    ! Calculate global-local DBCSR off-sets
    DO I=0,NSpace-1
       GLOffSet%I(I)=RowEnd%I(I)-RowBeg%I(I)+1
    ENDDO
    DO I=1,NSpace-1
       GLOffSet%I(I)=GLOffSet%I(I)+GLOffSet%I(I-1)
    ENDDO
    DO I=NSpace-1,1,-1
       GLOffSet%I(I)=GLOffSet%I(I-1)
    ENDDO
    GLOffSet%I(0)=0
  END SUBROUTINE GreedyDBCSRPartition

END MODULE ParseParallel
