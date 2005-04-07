MODULE PunchHDF
  USE InOut
  USE MemMan
  USE CellSets
  USE Indexing
  USE AtomPairs
  USE PrettyPrint
  USE GlobalScalars
  USE ControlStructures
  USE OptionKeys
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

    ! Not sure if this will cause a bug...
!
!    CALL Put(S%Action,'action_state')
!    CALL Put(S%SubAction,'subaction_state')
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
    CHARACTER(LEN=DCL) :: chGEO
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
       CALL Put(.FALSE.,'archivedensity')
!      MD Stuff
       CALL Put(.FALSE.,'DoingMD')
       CALL Put(0,'DMPOrder')
#ifdef PARALLEL
       CALL Put(0,'LineLocExist')
       CALL New(ETDirArr,P%NSpace-1)
       CALL New(ETRootArr,P%NSpace-1)
       ETDirArr%I(:) = 10
       ETRootArr%D(:) = 10.0
       CALL Put(ETDirArr,'ETDirArr')
       CALL Put(ETRootArr,'ETRootArr')
       CALL Put(0,'QLineLoc')
       CALL Put(ETDirArr,'QETDir')
       CALL Put(ETRootArr,'QETRoot')
       CALL Delete(ETDirArr)
       CALL Delete(ETRootArr)
       CALL Put(0,'GONXPartExist')
       CALL Put(0,'ONXPartExist')
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

    TYPE(FileNames)    :: N
    TYPE(BasisSets)    :: B
    TYPE(Geometries)   :: G    
    TYPE(CellSet)      :: CS_IN,CS_OUT
    INTEGER            :: cBAS,cGEO,iCLONE,HDFFileID,I
    CHARACTER(LEN=DCL) :: chGEO
!---------------------------------------------------------------------------!
    chGEO=IntToChar(cGEO)
    HDFFileID=OpenHDF(N%HFile)
    DO iCLONE=1,G%Clones
       G%Clone(iCLONE)%Confg=cGEO
!      Set the correct PBC cell set list
       CALL SetLatticeVectors(G%Clone(iCLONE),B%AtomPairThresh(iCLONE,cBAS),CS_IN,CS_OUT)
!
!      Make sure everything is wrapped correctly
       G%Clone(iCLONE)%Carts%D = G%Clone(iCLONE)%Carts%D      
       CALL CalFracCarts(G%Clone(iCLONE))
       CALL WrapAtoms(G%Clone(iCLONE))
       G%Clone(iCLONE)%Carts%D = G%Clone(iCLONE)%Carts%D  
!
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
!      If we have ECPs, temporarily reset this geometries nuclear charges
       IF(B%BSets(iCLONE,cBAS)%HasECPs) &
          CALL SetAtomCharges(G%Clone(iCLONE),B%BSets(iCLONE,cBAS))
!      Put the geometry to this group ...
       CALL Put(G%Clone(iCLONE),chGEO)
!      Store the CellSets
       CALL Put(CS_IN ,'CS_IN' ,Tag_O=IntToChar(cBAS))
       CALL Put(CS_OUT,'CS_OUT',Tag_O=IntToChar(cBAS))
!      Close this clones group
       CALL CloseHDFGroup(HDF_CurrentID)
!      And free memory for the the lattice vectors 
       CALL Delete(CS_IN%CellCarts)
       CALL Delete(CS_OUT%CellCarts)
!      And unset the nuclear charges in the case of ECPs
       IF(B%BSets(iCLONE,cBAS)%HasECPs) &
          CALL UnSetAtomCharges(G%Clone(iCLONE),B%BSets(iCLONE,cBAS))
    ENDDO
    CALL CloseHDF(HDFFileID)
  END SUBROUTINE GeomArchive
!==============================================================================
!
!==============================================================================
  SUBROUTINE GeomReArchive(N,O,G)
    TYPE(FileNames)  :: N
    TYPE(Options)    :: O
    TYPE(Geometries) :: G    
    TYPE(CellSet)    :: CS
    INTEGER          :: cBAS,cGEO,iCLONE,HDFFileID
    INTEGER          :: LastGeo,LastBas,iGEO,NCart,NatmsLoc,iGEOStart
    INTEGER          :: iGEOMem
    CHARACTER(LEN=DCL) :: chGEO
    TYPE(CRDS)       :: GMLoc,GMAux
    !
    !-----------------------------------------------------------------!
    ! it is supposed that the number of clones did not change at restart
    !
    iGEOMem=11
    LastGeo=O%RestartState%I(3)
    LastBas=O%RestartState%I(2)
    IGEOStart=MAX(LastGeo-iGEOMem,2)
    !
    NatmsLoc=G%Clone(1)%Natms
    NCart=3*NatmsLoc
    !
    DO iCLONE=1,G%Clones
      ! Reachive very first geometry, as it serves 
      ! as a reference for PBC optimizations
      CALL RearchIGEO(N,O,G,iCLONE,1,LastGeo) 
      ! Reachive GMLoc-s 
      DO iGEO=iGEOStart,LastGeo
        CALL RearchIGEO(N,O,G,iCLONE,iGEO,LastGeo)
      ENDDO
!      ! Rearchive CS 
!      HDFFileID=OpenHDF(N%RFile)
!      HDF_CurrentID=OpenHDFGroup(HDFFileID, &
!                    "Clone #"//TRIM(IntToChar(iCLONE)))
!        CALL Get(CS,Tag_O=IntToChar(LastBAS))
!      CALL CloseHDFGroup(HDF_CurrentID)
!      CALL CloseHDF(HDFFileID)
!      HDFFileID=OpenHDF(N%HFile)
!      HDF_CurrentID=OpenHDFGroup(HDFFileID, &
!                    "Clone #"//TRIM(IntToChar(iCLONE)))
!        CALL Put(CS,Tag_O=IntToChar(LastBAS))
!      CALL CloseHDFGroup(HDF_CurrentID)
!      CALL CloseHDF(HDFFileID)
!      CALL Delete(CS)
    ENDDO
  END SUBROUTINE GeomReArchive
!
!---------------------------------------------------------------
!
   SUBROUTINE RearchIGEO(N,O,G,iCLONE,iGEO,LastGeo)
     TYPE(FileNames)  :: N
     TYPE(Options)    :: O
     TYPE(Geometries) :: G    
     INTEGER          :: iCLONE,HDFFileID
     INTEGER          :: iGEO,NCart,NatmsLoc,LastGeo
     CHARACTER(LEN=DCL) :: chGEO
     TYPE(CRDS)       :: GMLoc,GMAux
     !
     NatmsLoc=G%Clone(1)%Natms
     NCart=3*NatmsLoc
     !
     HDFFileID=OpenHDF(N%RFile)
     HDF_CurrentID=OpenHDFGroup(HDFFileID, &
                   "Clone #"//TRIM(IntToChar(iCLONE)))
       chGEO=IntToChar(iGeo)
       CALL Get(GMLoc,chGEO)
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(HDFFileID)
     !
     HDFFileID=OpenHDF(N%HFile)
     HDF_CurrentID=OpenHDFGroup(HDFFileID, &
                   "Clone #"//TRIM(IntToChar(iCLONE)))
       chGEO=IntToChar(iGEO)
       CALL Put(GMLoc,chGEO)
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(HDFFileID)
     CALL Delete(GMLoc)
   END SUBROUTINE RearchIGEO
!
!---------------------------------------------------------------
!
!==============================================================================
! SET UP SUMMATION OF LATTICE VECTORS, ACCOUNTING FOR CELL SIZE AND SHAPE
!==============================================================================
  SUBROUTINE SetLatticeVectors(G,AtomPairThresh,CS_IN,CS_OUT,Rad_O)
    TYPE(CRDS)                :: G
    TYPE(CellSet)             :: CS_IN,CS_OUT
    REAL(DOUBLE), OPTIONAL    :: Rad_O
    REAL(DOUBLE)              :: AtomPairThresh,Radius_in,Radius_out
    INTEGER                   :: IL
    INTEGER                   :: MaxCell
!------------------------------------------------------------------------------
!   CS_OUT
    MaxCell=1
    IF(G%PBC%Dimen==0) THEN
       MaxCell=1
    ELSEIF(G%PBC%Dimen==1) THEN
       MaxCell=50
    ELSEIF(G%PBC%Dimen==2) THEN
       MaxCell=500
    ELSEIF(G%PBC%Dimen==3) THEN
       MaxCell=5000
    ENDIF 
!
    IF(PRESENT(Rad_O)) THEN
       Radius_out = 1.0D0*Rad_O
       Radius_in  = 1.5D0*Rad_O
    ELSE
       Radius_out = (One+1.D-14)*MaxBoxDim(G)+SQRT(AtomPairThresh)
       Radius_in  = (One+1.D-14)*MaxBoxDim(G)+2.0D0*SQRT(AtomPairThresh)
    ENDIF
!   Determine PFFOverde
    IF(G%PBC%PFFOvRide) THEN
       IL = G%PBC%PFFMaxLay
       CALL New_CellSet_Cube(CS_OUT,G%PBC%AutoW%I,G%PBC%BoxShape%D,(/IL,IL,IL/)      ,MaxCell_O=MaxCell)
       CALL New_CellSet_Cube(CS_IN ,G%PBC%AutoW%I,G%PBC%BoxShape%D,(/2*IL,2*IL,2*IL/),MaxCell_O=MaxCell)
    ELSE
       CALL New_CellSet_Sphere(CS_OUT,G%PBC%AutoW%I,G%PBC%BoxShape%D,Radius_out ,MaxCell_O=MaxCell)
       CALL New_CellSet_Sphere(CS_IN ,G%PBC%AutoW%I,G%PBC%BoxShape%D,Radius_in  ,MaxCell_O=MaxCell)
    ENDIF
!
    CALL Sort_CellSet(CS_IN)
    CS_IN%Radius  = SQRT(CS_IN%CellCarts%D(1,1)**2 +CS_IN%CellCarts%D(2,1)**2 +CS_IN%CellCarts%D(3,1)**2)
    CALL Sort_CellSet(CS_OUT)
    CS_OUT%Radius = SQRT(CS_OUT%CellCarts%D(1,1)**2+CS_OUT%CellCarts%D(2,1)**2+CS_OUT%CellCarts%D(3,1)**2)
!
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
!---------------------------------------------------------------------
!
  SUBROUTINE InitState(S,O)
    TYPE(State)    :: S
    TYPE(Options)  :: O
    INTEGER        :: LastGEO,LastBAS
    !
    IF(O%Guess==GUESS_EQ_RESTART) THEN
      LastBAS=O%RestartState%I(2)
      LastGEO=O%RestartState%I(3)
    ELSE
      LastBAS=1
      LastGEO=1
    ENDIF
    S%Previous%I=(/0,LastBAS,LastGEO/)
    S%Current%I=(/1,LastBAS,LastGEO/)
  END SUBROUTINE InitState
!
!-----------------------------------------------------------------
!
  SUBROUTINE InitGlobal(C)
    TYPE(Controls) :: C
    !
    ! Initialize the HDF archival file
    CALL InitArchive(C%Nams)
    ! Initialize HDF files for clones
    CALL InitClones(C%Nams,C%MPIs,C%Sets,C%Geos)
    ! Do rearchivation if requested
    IF(C%Opts%Guess==GUESS_EQ_RESTART) THEN
      CALL GeomReArchive(C%Nams,C%Opts,C%Geos)
    ENDIF
    ! Set initial state vector
    CALL InitState(C%Stat,C%Opts)
  END SUBROUTINE InitGlobal

  SUBROUTINE SetAtomCharges(G,B)
    TYPE(CRDS) :: G
    TYPE(BSET) :: B
    INTEGER    :: I
    DO I=1,G%NAtms
       G%AtNum%D(I)=G%AtNum%D(I)-B%NCoreEl%D(G%AtTyp%I(I))
    ENDDO
    G%NElec=0
    DO I=1,G%NAtms
       IF(G%AtNum%D(i)<105.D0) THEN
          G%NElec=G%NElec+G%AtNum%D(I)
       ENDIF
    ENDDO
    G%NElec=G%NElec-G%TotCh
  END SUBROUTINE SetAtomCharges

  SUBROUTINE UnSetAtomCharges(G,B)
    TYPE(CRDS) :: G
    TYPE(BSET) :: B
    INTEGER    :: I
    DO I=1,G%NAtms
       G%AtNum%D(I)=G%AtNum%D(I)+B%NCoreEl%D(G%AtTyp%I(I))
    ENDDO
    G%NElec=0
    DO I=1,G%NAtms
       IF(G%AtNum%D(i)<105.D0) THEN
          G%NElec=G%NElec+G%AtNum%D(I)
       ENDIF
    ENDDO
    G%NElec=G%NElec-G%TotCh
  END SUBROUTINE UnSetAtomCharges
!
!-------------------------------------------------------------------
!
   SUBROUTINE GetRefXYZ(HFileIn,RefXYZ,iCLONE)
     TYPE(FileNames)      :: Nams
     CHARACTER(LEN=*)     :: HFileIn
     INTEGER              :: HDFFileID,iCLONE
     TYPE(DBL_RNK2)       :: RefXYZ
     !
     HDFFileID=OpenHDF(HFileIn)
     HDF_CurrentID= &
       OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
     CALL Get(RefXYZ,'Abcartesians',Tag_O=TRIM(IntToChar(1)))
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(HDFFileID)
   END SUBROUTINE GetRefXYZ
!
!-------------------------------------------------------------------
!
   SUBROUTINE CollectPast(RefXYZ,SRStruct,RefStruct,RefGrad,SRDispl, &
                          HFileIn,NatmsLoc,Mem,iGEO,iCLONE)
     TYPE(DBL_RNK2)   :: SRStruct,RefStruct,RefGrad,SRDispl
     TYPE(DBL_RNK2)   :: Aux
     TYPE(PBCInfo)    :: PBC
     TYPE(DBL_VECT)   :: Vect
     INTEGER          :: NatmsLoc,Mem,iGEO,iCLONE,NCartS,NatmsS
     CHARACTER(LEN=*) :: HFileIn
     INTEGER          :: HDFFileID,ICount,I,J,IStart,NCart,IGeom,K,L
     REAL(DOUBLE),DIMENSION(:,:) :: RefXYZ
     !
     NCart=3*NatmsLoc
     NatmsS=NatmsLoc-3
     NCartS=3*NatmsS
     IStart=iGEO-Mem+1
     !
     ! Get GDIIS memory of Cartesian coords and grads
     !
     CALL New(SRStruct,(/NCart,Mem/))
     CALL New(RefStruct,(/NCart,Mem/))
     CALL New(RefGrad,(/NCart,Mem/))
     CALL New(SRDispl,(/NCart,Mem/))
     !
     CALL New(Vect,NCartS)
     CALL New(Aux,(/3,NatmsS/))
     CALL New(PBC)
     !
     HDFFileID=OpenHDF(HFileIn)
     HDF_CurrentID= &
       OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
     !
     DO IGeom=IStart,iGEO
       ICount=IGeom-IStart+1
       CALL Get(PBC,Tag_O=TRIM(IntToChar(IGeom)))
       !!!!
!       WRITE(*,*) 'WARNING! Hardwired cleaning of lattice forces in CollectPast'
!       WRITE(Out,*) 'WARNING! Hardwired cleaning of lattice forces in CollectPast'
!         PBC%LatFrc%D(2:3,1)=Zero
!         PBC%LatFrc%D(3,2)=Zero
!         CALL Put(PBC,Tag_O=TRIM(IntToChar(IGeom)))
!       !!!!
       L=NCartS
       DO K=1,3
         DO J=1,3
           L=L+1
           RefStruct%D(L,ICount)=PBC%BoxShape%D(J,K)
             RefGrad%D(L,ICount)=  PBC%LatFrc%D(J,K)
         ENDDO
       ENDDO
       CALL Get(Aux,'Displ',Tag_O=TRIM(IntToChar(IGeom)))
       CALL CartRNK2ToCartRNK1(Vect%D,Aux%D)
       DO J=1,NCartS ; SRStruct%D(J,ICount)=Vect%D(J) ; ENDDO
       CALL Get(Aux,'Abcartesians',Tag_O=TRIM(IntToChar(IGeom)))
       CALL ConvertToXYZRef(Aux%D,RefXYZ,PBC%Dimen, &
                            BoxShape_O=PBC%BoxShape%D)
       CALL CartRNK2ToCartRNK1(Vect%D,Aux%D)
       DO J=1,NCartS ; RefStruct%D(J,ICount)=Vect%D(J) ; ENDDO
       CALL Get(Aux,'gradients',Tag_O=TRIM(IntToChar(IGeom)))
       CALL CartRNK2ToCartRNK1(Vect%D,Aux%D)
       DO J=1,NCartS ; RefGrad%D(J,ICount)=Vect%D(J) ; ENDDO
     ENDDO
     CALL Delete(Aux)
     CALL Delete(Vect)
     CALL Delete(PBC)
       SRDispl%D=SRStruct%D-RefStruct%D
     !
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(HDFFileID)
   END SUBROUTINE CollectPast
! 
!-------------------------------------------------------------------
!
END MODULE PunchHDF
