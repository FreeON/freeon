MODULE PunchHDF
  USE InOut
  USE MemMan
  USE CellSets
  USE Indexing
  USE AtomPairs
  USE PrettyPrint
  USE GlobalScalars
  USE ControlStructures
  USE InCoords
  IMPLICIT NONE
CONTAINS
  !==============================================================================
  !
  !==============================================================================
  SUBROUTINE InitArchive(N)
    TYPE(FileNames)  :: N
    !---------------------------------------------------------------------------!
    HDF_CurrentID=InitHDF(N%HFile)
    HDF_CurrentID=OpenHDF(N%HFile)
    CALL Put(N%RFile,'oldinfo')
    CALL Put(N%IFile,'inPutfile')
    CALL Put(N%HFile,'infofile')
    CALL Put(N%LFile,'logfile')
    CALL Put(N%OFile,'outPutfile')
    CALL CloseHDF(HDF_CurrentID)
  END SUBROUTINE InitArchive
  !==============================================================================
  !
  !==============================================================================
  SUBROUTINE MPIsArchive(N,NSpace,Clump)
    TYPE(FileNames)      :: N
    INTEGER,DIMENSION(2) :: Clump
    INTEGER              :: NSpace
    TYPE(INT_VECT)   :: ST
    !---------------------------------------------------------------------------!
    HDF_CurrentID=OpenHDF(N%HFile)
#ifdef PARALLEL    
    CALL New(ST,3)
    ST%I=(/NSpace,Clump(1),Clump(2)/)
    CALL Put(ST,'SpaceTime')
    CALL Delete(ST)
#else
    CALL Put(Clump(2),'SpaceTime')
#endif
    CALL CloseHDF(HDF_CurrentID)
  END SUBROUTINE MPIsArchive

  SUBROUTINE StateArchive(N,S)
    TYPE(FileNames) :: N
    TYPE(State)     :: S
    HDF_CurrentID=OpenHDF(N%HFile)
    CALL Put(S%Current,'current_state')
    CALL Put(S%Previous,'previous_state')
    CALL Put(S%Action,'action_state')
    CALL Put(S%SubAction,'subaction_state')
    CALL CloseHDF(HDF_CurrentID)
  END SUBROUTINE StateArchive
  !==============================================================================
  ! THIS IS WHERE ALL SORTS OF MISC DATA SPACE IS INITIALIZED IN THE HDF5 FILE
  ! IN ORDER TO AVOID CHANGING THE DATA SPACE WHEN THE HDF FILE HAS BEEN OPENED
  ! SIMULTANEOUSLY BY EACH CLONE
  !==============================================================================
  SUBROUTINE InitClones(N,P,B,G)
    TYPE(FileNames)  :: N
    TYPE(Parallel)   :: P
    TYPE(Geometries) :: G    
    TYPE(BasisSets)  :: B
    TYPE(DBL_VECT)   :: DoubleVect
    TYPE(BCSR)       :: DM
    CHARACTER(LEN=DCL) :: FailedProgram
    INTEGER          :: iGEO,iCLONE,iBAS,HDFFileID,I,MaxEll,MaxAtoms,MaxBloks,MaxNon0s
    CHARACTER(LEN=2) :: chGEO
    TYPE(INT_VECT)   :: ETDirArr
    TYPE(DBL_VECT)   :: ETRootArr
    !---------------------------------------------------------------------------!
    chGEO=IntToChar(iGEO)
    HDFFileID=OpenHDF(N%HFile)
    HDF_CurrentID=HDFFileID
    ! Find max dimensions for BCSR matrices
    MaxAtoms=0
    MaxBloks=0
    MaxNon0s=0
    DO iBAS=1,B%NBSets
       MaxAtoms=MAX(MaxAtoms,B%MxAts(iBAS))
       MaxBloks=MAX(MaxBloks,B%MxBlk(iBAS))
       MaxNon0s=MAX(MaxNon0s,B%MxN0s(iBAS))
    ENDDO
    CALL New(DM,(/MaxAtoms,MaxBloks,MaxNon0s/))
    !
    CALL Put(G%Clones,'clones')
    DO iCLONE=1,G%Clones
       HDF_CurrentID=InitHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
       ! Create data space for grouped objects to preserve structure of HDF5 
       ! to avoid multiple clones simultaneously creating new data
       CALL New(DoubleVect,3)
       DoubleVect%D=BIG_DBL
       CALL Put(DoubleVect,'dipole')
       CALL Delete(DoubleVect)
       CALL New(DoubleVect,6)
       DoubleVect%D=BIG_DBL
       CALL Put(DoubleVect,'quadrupole')
       CALL Delete(DoubleVect)
       CALL Put(BIG_DBL,'e_nucleartotal')
       CALL Put(BIG_DBL,'exc')
       CALL Put(BIG_DBL,'homolumogap')
       CALL Put(BIG_DBL,'e_electronictotal')
       CALL Put(BIG_DBL,'etot')
       CALL Put(BIG_DBL,'dmax')
       CALL Put(BIG_DBL,'diiserr')
       CALL Put(.TRUE.,'programfailed')
#ifdef PARALLEL
       CALL Put(0,'LineLocExist')
       CALL New(ETDirArr,P%NSpace-1)
       CALL New(ETRootArr,P%NSpace-1)
       ETDirArr%I(:) = 10
       ETRootArr%D(:) = 10.0
       CALL Put(ETDirArr,'ETDirArr')
       CALL Put(ETRootArr,'ETRootArr')
       CALL Delete(ETDirArr)
       CALL Delete(ETRootArr)
#endif
       DO I=1,LEN(FailedProgram)
          FailedProgram(I:I)='X'
       ENDDO
       CALL Put(FailedProgram,'failedprogram')
       MaxEll=G%Clone(iCLONE)%PBC%PFFMaxEll
       CALL Put(MaxEll,'MaxEll')
       CALL New(DoubleVect,LSP(2*MaxEll),0)
       DoubleVect%D=BIG_DBL
       CALL Put(DoubleVect,'PFFTensorC')
       CALL Put(DoubleVect,'PFFTensorS')
       CALL Delete(DoubleVect) 
       CALL Put(DM,'CurrentDM',CheckPoint_O=.TRUE.)
       CALL CloseHDFGroup(HDF_CurrentID)
    ENDDO
    CALL Delete(DM)
    CALL CloseHDF(HDFFileID)
  END SUBROUTINE InitClones
!==============================================================================
!
!==============================================================================
  SUBROUTINE GeomArchive(cBAS,cGEO,N,B,G)
    TYPE(FileNames)  :: N
    TYPE(BasisSets)  :: B
    TYPE(Geometries) :: G    
    TYPE(CellSet)    :: CS
    INTEGER          :: cBAS,cGEO,iCLONE,HDFFileID
    CHARACTER(LEN=2) :: chGEO
    !---------------------------------------------------------------------------!
    chGEO=IntToChar(cGEO)
    HDFFileID=OpenHDF(N%HFile)
    DO iCLONE=1,G%Clones
       G%Clone(iCLONE)%Confg=cGEO
!      Set the correct PBC cell set list
       CALL SetLatticeVectors(G%Clone(iCLONE),CS,B%AtomPairThresh(iCLONE,cBAS))
!      Make sure everything is wrapped correctly
       CALL WrapAtoms(G%Clone(iCLONE))
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
!      Put the geometry to this group ...
       CALL Put(G%Clone(iCLONE),chGEO)
!      ... and the corresponding lattice vectors which for now ONLY DEPEND
!      ON THE BASIS SET.  THIS WILL HAVE TO BE CHANGED WHEN WE ADD BOX FORCES! 
       CALL Put(CS,Tag_O=IntToChar(cBAS))
!!$       CALL Put(CS,'CS_OUT',Tag_O=IntToChar(cBAS))
!!$       CALL Put(CS,'CS_IN' ,Tag_O=IntToChar(cBAS))
       ! Close this clones group
       CALL CloseHDFGroup(HDF_CurrentID)
       ! And free memory for the the lattice vectors 
       CALL Delete(CS%CellCarts)
    ENDDO
    CALL CloseHDF(HDFFileID)
  END SUBROUTINE GeomArchive
!
!---------------------------------------------------------------------
!
   SUBROUTINE GDIISArch(Nams,iCLONE,XYZ_O,Vect_O,Tag_O)
     TYPE(FileNames)                      :: Nams
     INTEGER                              :: iCLONE
     INTEGER                              :: IGeom,NCart
     INTEGER                              :: NatmsLoc,J,I
     REAL(DOUBLE),DIMENSION(:,:),OPTIONAL :: XYZ_O 
     REAL(DOUBLE),DIMENSION(:),OPTIONAL   :: Vect_O
     CHARACTER(LEN=*),OPTIONAL            :: Tag_O
     TYPE(DBL_RNK2)                       :: XYZ 
     TYPE(DBL_VECT)                       :: AuxVect
     TYPE(DBL_VECT)                       :: Vect
     INTEGER                              :: GDIISMemory
     INTEGER                              :: SRMemory
     INTEGER                              :: RefMemory
     INTEGER                              :: CartGradMemory
     INTEGER                              :: IntGradMemory
     INTEGER                              :: HDFFileID
     !
     HDFFileID=OpenHDF(Nams%HFile)
     HDF_CurrentID= &
       OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
     !
     ! Increment GDIISMemory
     !
     IF(PRESENT(Tag_O).AND.TRIM(Tag_O)=='SR') THEN
       CALL Get(SrMemory,'SRMemory')
       SRMemory=SRMemory+1
       IGeom=SRMemory
       CALL Put(SRMemory,'SRMemory')
     ENDIF
     !
     IF(PRESENT(Tag_O).AND.TRIM(Tag_O)=='Ref') THEN
       CALL Get(RefMemory,'RefMemory')
       RefMemory=RefMemory+1
       IGeom=RefMemory
       CALL Put(RefMemory,'RefMemory')
     ENDIF
     !
     IF(PRESENT(Tag_O).AND.TRIM(Tag_O)=='CartGrad') THEN
       CALL Get(CartGradMemory,'CartGradMemory')
       CartGradMemory=CartGradMemory+1
       IGeom=CartGradMemory
       CALL Put(CartGradMemory,'CartGradMemory')
     ENDIF
     !
     IF(PRESENT(XYZ_O)) THEN
       NatmsLoc=SIZE(XYZ_O,2)
       NCart=3*NatmsLoc
       !
       ! Put Geometry of simple relaxation set into HDF, for GDIIS.
       !
       CALL New(AuxVect,NCart)
         CALL CartRNK2ToCartRNK1(AuxVect%D,XYZ_O)
         CALL Put(AuxVect,TRIM(Tag_O)//TRIM(IntToChar(IGeom)))
       CALL Delete(AuxVect)
       !
     ELSE IF(PRESENT(Vect_O)) THEN
       !
       NCart=SIZE(Vect_O)
       CALL New(Vect,NCart)
       Vect%D=Vect_O
         CALL Put(Vect,TRIM(Tag_O)//TRIM(IntToChar(IGeom)))
       CALL Delete(Vect)
     ELSE
       CALL Halt('No input coordinates in GDIISArch')
     ENDIF
     !
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(HDFFileID)
   END SUBROUTINE GDIISArch
!
!---------------------------------------------------------------------
!
!==============================================================================
! SET UP SUMMATION OF LATTICE VECTORS, ACCOUNTING FOR CELL SIZE AND SHAPE
!==============================================================================
  SUBROUTINE SetLatticeVectors(G,C,AtomPairThresh,Rad_O)
    TYPE(CRDS)                :: G
    TYPE(CellSet)             :: C
    REAL(DOUBLE), OPTIONAL    :: Rad_O
    REAL(DOUBLE)              :: AtomPairThresh,Radius
    INTEGER                   :: IL
!------------------------------------------------------------------------------
!   First we create a CellSet with a Distance of MaxBoxDim+SQRT(AtomPairThresh)
!   Some of these boxes will be to far away, next we determine how close the 
!   closest Faces of each cell are, if they are within SQRT(AtomPairThresh), 
!   then we keep them
!------------------------------------------------------------------------------
    IF(PRESENT(Rad_O)) THEN
       Radius = Rad_O
       CALL New_CellSet_Sphere(C,G%PBC%AutoW,G%PBC%BoxShape,Radius)   
    ELSE
       IF(G%PBC%PFFOvRide) THEN
          IL = G%PBC%PFFMaxLay
          CALL New_CellSet_Cube(C,G%PBC%AutoW,G%PBC%BoxShape,(/IL,IL,IL/))
       ELSE
!          Radius = (One+1.D-14)*MaxAtomDist(G)+SQRT(AtomPairThresh)
          Radius = (One+1.D-14)*MaxBoxDim(G)+SQRT(AtomPairThresh)
          CALL New_CellSet_Sphere(C,G%PBC%AutoW,G%PBC%BoxShape,Radius)
       ENDIF
    ENDIF
!
    CALL Sort_CellSet(C)
    C%Radius = SQRT(C%CellCarts%D(1,1)**2+C%CellCarts%D(2,1)**2+C%CellCarts%D(3,1)**2)
  END SUBROUTINE SetLatticeVectors
!==============================================================================
!
!==============================================================================
  SUBROUTINE BSetArchive(cBAS,N,O,G,B,M)    
    TYPE(FileNames)  :: N
    TYPE(Options)    :: O
    TYPE(Geometries) :: G
    TYPE(BasisSets)  :: B
    TYPE(Parallel)   :: M
    INTEGER          :: cBAS,iCLONE,HDFFileID
    CHARACTER(LEN=2) :: chBAS
    !---------------------------------------------------------------------------!
    chBAS=IntToChar(cBAS)
    HDFFileID=OpenHDF(N%HFile)
    HDF_CurrentID=HDFFileID
    CALL Put(B%MxAts(cBAS),'maxatms',chBAS)
    CALL Put(B%MxBlk(cBAS),'maxblks',chBAS)
    CALL Put(B%MxN0s(cBAS),'maxnon0',chBAS)
    CALL Put(O%Models(cBAS),'ModelChemistry',chBAS)
#ifdef PARALLEL
    CALL Put(M%MxAtsNode(cBAS),'maxatmsnode',Tag_O=chBAS)
    CALL Put(M%MxBlkNode(cBAS),'maxblksnode',Tag_O=chBAS)
    CALL Put(M%MxN0sNode(cBAS),'maxnon0node',Tag_O=chBAS)
#endif
    DO iCLONE=1,G%Clones
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
       CALL Put(B%NExpt(cBAS),'nexpt',chBAS)
       CALL Put(O%Thresholds(cBAS),chBAS)
       CALL Put(B%BSets(iCLONE,cBAS),Tag_O=chBAS)
       CALL Put(B%BSiz(iCLONE,cBAS),'atsiz',Tag_O=chBAS)
       CALL Put(B%OffS(iCLONE,cBAS),'atoff',Tag_O=chBAS)
       CALL Put(B%DExpt(iCLONE,cBAS),'dexpt',Tag_O=chBAS)
       CALL Put(B%Lndex(iCLONE,cBAS),'lndex',Tag_O=chBAS)
#ifdef PARALLEL
       CALL Put(M%Beg(iCLONE,cBAS),'beg',Tag_O=chBAS)
       CALL Put(M%End(iCLONE,cBAS),'end',Tag_O=chBAS)
       CALL Put(M%GLO(iCLONE,cBAS),'dbcsroffsets',Tag_O=chBAS)
#endif
       CALL CloseHDFGroup(HDF_CurrentID)
    ENDDO
    CALL CloseHDF(HDFFileID)
  END SUBROUTINE BSetArchive

END MODULE PunchHDF
