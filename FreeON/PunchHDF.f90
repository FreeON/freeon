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

MODULE PunchHDF
  USE InOut
  USE MemMan
  USE PBC
  USE Indexing
  USE AtomPairs
  USE PrettyPrint
  USE GlobalScalars
  USE ControlStructures
  USE OptionKeys

  IMPLICIT NONE

CONTAINS

  SUBROUTINE InitArchive(N)
    TYPE(FileNames)  :: N

    CALL MondoLog(DEBUG_MAXIMUM, "InitArchive", "initializing hdf file")

    HDF_CurrentID=InitHDF(N%HFile)
    HDF_CurrentID=OpenHDF(N%HFile)
    CALL Put(N%RFile,'oldinfo')
    CALL Put(N%IFile,'inPutfile')
    CALL Put(N%HFile,'infofile')
    CALL Put(N%LFile,'logfile')
    CALL Put(N%OFile,'outPutfile')
    CALL Put(RecycleHDF, "RecycleHDF")
    CALL CloseHDF(HDF_CurrentID)
  END SUBROUTINE InitArchive

  SUBROUTINE MPIsArchive(N,NSpace,Clump)
    TYPE(FileNames)      :: N
    INTEGER,DIMENSION(2) :: Clump
    INTEGER              :: NSpace
    TYPE(INT_VECT)       :: ST

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

  !==============================================================================
  ! THIS IS WHERE ALL SORTS OF MISC DATA SPACE IS INITIALIZED IN THE HDF5 FILE
  ! IN ORDER TO AVOID CHANGING THE DATA SPACE WHEN THE HDF FILE HAS BEEN OPENED
  ! SIMULTANEOUSLY BY EACH CLONE
  !==============================================================================
  SUBROUTINE InitClones(N,P,B,G,D)
    TYPE(FileNames)    :: N
    TYPE(Parallel)     :: P
    TYPE(Geometries)   :: G
    TYPE(BasisSets)    :: B
    TYPE(Dynamics)     :: D
    TYPE(DBL_VECT)     :: DoubleVect
    TYPE(DBL_RNK2)     :: DoubleRnk2
    TYPE(BCSR)         :: DM
    CHARACTER(LEN=DCL) :: FailedProgram
    INTEGER            :: iGEO,iCLONE,iBAS,HDFFileID,I,MaxEll,MaxAtoms,MaxBloks,MaxNon0s,J
    CHARACTER(LEN=DCL) :: chGEO
    TYPE(INT_VECT)     :: ETDirArr,IntVect
    TYPE(DBL_VECT)     :: ETRootArr
    TYPE(DBL_RNK2)     :: DblMat
    TYPE(CMPoles)      :: MP

    CALL MondoLog(DEBUG_MAXIMUM, "InitClones", "initializing hdf file")

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
#ifdef PARALLEL
    !If parallel and we have clones, we need to put the D matrices in the HDF.
    IF(G%Clones.GT.1) CALL New(DM,(/MaxAtoms,MaxBloks,MaxNon0s/))
#endif
    CALL Put(G%Clones,'clones')

    ! If we are doing NEB, then put the endpoints to HDF as well
    IF(SIZE(G%Clone)==G%Clones+2)THEN
      ! All thats going into the endpoints is the geometry
      HDF_CurrentID=InitHDFGroup(HDFFileID,"Clone #0")
      CALL Put(G%Clone(0))
      CALL CloseHDFGroup(HDF_CurrentID)
      HDF_CurrentID=InitHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(G%Clones+1)))
      CALL Put(G%Clone(G%Clones+1))
      CALL CloseHDFGroup(HDF_CurrentID)
    ENDIF
    !
    DO iCLONE=1,G%Clones
      HDF_CurrentID=InitHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
      ! Create data space for grouped objects to preserve structure of HDF5
      ! to avoid multiple clones simultaneously creating new data

      ! BEGIN STUPID MULTIPOLE REDUNDANCY
      CALL New(MP)
      MP%DPole%D=BIG_DBL
      MP%QPole%D=BIG_DBL
      CALL Put(MP)
      CALL Delete(MP)
      CALL New(DoubleVect,3)
      DoubleVect%D=BIG_DBL
      CALL Put(DoubleVect,'dipole')
      CALL Delete(DoubleVect)
      CALL New(DoubleVect,6)
      DoubleVect%D=BIG_DBL
      CALL Put(DoubleVect,'quadrupole')
      CALL Delete(DoubleVect)
      ! END STUPID REDUNDANCY

      CALL Put(BIG_DBL,'e_nucleartotal')
      CALL Put(BIG_DBL,'exc')
      CALL Put(BIG_DBL,'homolumogap')
      CALL Put(-1.0D0, "Entropy")
      CALL Put(BIG_DBL,'e_electronictotal')
      CALL Put(BIG_DBL,'etot')
      CALL Put(BIG_DBL,'dmax')
      CALL Put(BIG_DBL,'diiserr')
      CALL New(IntVect,2)
      IntVect%I=0
      CALL Put(IntVect,'diisinfo')
      CALL Delete(IntVect)
      CALL New(DblMat,(/DIIS_MAX_MATRIX_SIZE,DIIS_MAX_MATRIX_SIZE/))
      DblMat%D=Zero
      CALL Put(DblMat,'diismtrix')
      CALL Delete(DblMat)
      CALL Put(.TRUE.,'programfailed')
      CALL Put(.FALSE.,'archivedensity')

      ! MD and MC Stuff
      CALL Put(1 ,'MDIter')
      CALL Put(BIG_DBL,'MDTime')
      CALL Put(D%MDalpha, 'MDalpha')
      CALL Put(D%MDDampStep, 'MDDampStep')
      CALL Put(.FALSE.,'DoingMD')
      CALL Put(.FALSE.,'DoingHybridMC')
      CALL Put("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",'MDGeuss')
      CALL Put(BIG_DBL,'MCEtot0')
      CALL Put(BIG_DBL,'MCTemp0')
      CALL New(DoubleRnk2,(/3,MaxAtoms/))
      DoubleRnk2%D=BIG_DBL
      CALL Put(DoubleRnk2,'MCCarts0')
      CALL Delete(DoubleRnk2)
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
      CALL Put(MaxEll,'MaxPFFFEll')
      CALL New(DoubleVect,LSP(2*MaxEll),0)
      DoubleVect%D=BIG_DBL
      CALL Put(DoubleVect,'PFFTensorC')
      CALL Put(DoubleVect,'PFFTensorS')
      CALL Delete(DoubleVect)
#ifdef PARALLEL
      !If parallel and we have clones, we need to put the D matrices in the HDF.
      IF(G%Clones.GT.1) CALL Put(DM,'CurrentDM',CheckPoint_O=.TRUE.)
#endif
      CALL CloseHDFGroup(HDF_CurrentID)
    ENDDO
#ifdef PARALLEL
    !If parallel and we have clones, we need to put the D matrices in the HDF.
    IF(G%Clones.GT.1) CALL Delete(DM)
#endif
    CALL CloseHDF(HDFFileID)
  END SUBROUTINE InitClones

  SUBROUTINE GeomArchive(cBAS,cGEO,N,O,B,G)
    TYPE(Options)      :: O
    TYPE(FileNames)    :: N
    TYPE(BasisSets)    :: B
    TYPE(Geometries)   :: G
    TYPE(CellSet)      :: CS_IN,CS_OUT
    INTEGER            :: cBAS,cGEO,iCLONE,HDFFileID,I,NK,CF,PF
    CHARACTER(LEN=DCL) :: chGEO
    REAL(DOUBLE)       :: MinExpt

    chGEO=IntToChar(cGEO)

    CALL MondoLog(DEBUG_MAXIMUM, "GeomArchive", "Archiving data in HDF file "//TRIM(N%HFile))
    HDFFileID=OpenHDF(N%HFile)

    ! Set Confg in case we have more clones than 1:G%Clones, as is the case in
    ! NEB.
    DO iCLONE = LBOUND(G%Clone, 1), 1
      G%Clone(iCLONE)%Confg = cGEO
    ENDDO

    DO iCLONE = G%Clones+1, UBOUND(G%Clone, 1)
      G%Clone(iCLONE)%Confg = cGEO
    ENDDO

    DO iCLONE=1, G%Clones
      G%Clone(iCLONE)%Confg=cGEO
      CALL MkGeomPeriodic(G%Clone(iCLONE),B%BSets(iCLONE,cBAS), &
                          O%Thresholds(cBAS)%Dist,O%Thresholds(cBAS)%TwoE)
      HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
      ! If we have ECPs, temporarily reset this geometries nuclear charges
      IF(B%BSets(iCLONE,cBAS)%HasECPs) THEN
        CALL SetAtomCharges(G%Clone(iCLONE),B%BSets(iCLONE,cBAS))
      ENDIF
      ! Put the geometry to this group ...
      CALL Put(G%Clone(iCLONE),chGEO)
      ! Close this clones group
      CALL CloseHDFGroup(HDF_CurrentID)
      ! And unset the nuclear charges in the case of ECPs
      IF(B%BSets(iCLONE,cBAS)%HasECPs) THEN
        CALL UnSetAtomCharges(G%Clone(iCLONE),B%BSets(iCLONE,cBAS))
      ENDIF
    ENDDO
    CALL CloseHDF(HDFFileID)
    !CALL PPrint(G%Clone(1),Unit_O=6)
    CALL MondoLog(DEBUG_MAXIMUM, "GeomArchive", "done archiving")
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
    !-----------------------------------------------------------------!
    ! it is supposed that the number of clones did not change at restart
    !
    iGEOMem=101
    LastGeo=O%RestartState%I(3)
    LastBas=O%RestartState%I(2)
    IGEOStart=MAX(LastGeo-iGEOMem,1)
    !
    DO iCLONE=1,G%Clones
      ! Reachive very first geometry, as it serves
      ! as a reference for PBC optimizations
      CALL ReArchIGEO(N,O,G,iCLONE,1,LastGeo)
      ! Reachive GMLoc-s
      DO iGEO=iGEOStart,LastGeo
        CALL ReArchIGEO(N,O,G,iCLONE,iGEO,LastGeo)
      ENDDO
    ENDDO
  END SUBROUTINE GeomReArchive
  !
  !---------------------------------------------------------------
  !
  SUBROUTINE ReArchIGEO(N,O,G,iCLONE,iGEO,LastGeo)
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
  END SUBROUTINE ReArchIGEO

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
    CALL Put(O%NSMat(cBAS),'NSMat',chBAS)
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
    IF(O%Guess==GUESS_EQ_RESTART.OR.O%Guess==GUESS_EQ_NUGUESS) THEN
      LastBAS=O%RestartState%I(2)
      LastGEO=O%RestartState%I(3)
    ELSE
      LastBAS=1
      LastGEO=1
    ENDIF
    S%Previous%I=(/0,LastBAS,LastGEO/)
    S%Current%I=(/1,LastBAS,LastGEO/)
  END SUBROUTINE InitState
  !-----------------------------------------------------------------
  !
  SUBROUTINE StateArchive(N,G,S,Init_O)
    TYPE(FileNames) :: N
    TYPE(Geometries):: G
    TYPE(State)     :: S
    LOGICAL,OPTIONAL:: Init_O
    INTEGER         :: HDFFileID,iCLONE,J
    HDFFileID=OpenHDF(N%HFile)
    HDF_CurrentID=HDFFileID
    CALL Put(S%Current,'current_state')
    CALL Put(S%Previous,'previous_state')
    IF(PRESENT(Init_O))THEN
      IF(Init_O)THEN
        DO iCLONE=1,G%Clones
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
          CALL Put(BIG_DBL,'E_NuclearTotal',Stats_O=S%Current%I)
          CALL Put(BIG_DBL,'Exc',Stats_O=S%Current%I)
          CALL Put(BIG_DBL,'Etot',Stats_O=S%Current%I)
          CALL CloseHDFGroup(HDF_CurrentID)
        ENDDO
      ENDIF
    ENDIF
    CALL CloseHDF(HDFFileID)
  END SUBROUTINE StateArchive
  !-----------------------------------------------------------------
  !
  SUBROUTINE InitGlobal(C)
    TYPE(Controls) :: C

    CALL MondoLog(DEBUG_MAXIMUM, "InitGlobal", "initializing...")

    ! Initialize the HDF archival file
    CALL InitArchive(C%Nams)
    ! Initialize HDF files for clones
    CALL InitClones(C%Nams, C%MPIs, C%Sets, C%Geos, C%Dyns)
    ! Do rearchivation if requested
    IF(C%Opts%Guess==GUESS_EQ_RESTART.OR.C%Opts%Guess==GUESS_EQ_NUGUESS)THEN
      CALL GeomReArchive(C%Nams,C%Opts,C%Geos)
    ENDIF
    ! Set initial state vector
    CALL InitState(C%Stat,C%Opts)
  END SUBROUTINE InitGlobal

  SUBROUTINE SetAtomCharges(G,B)
    TYPE(CRDS) :: G
    TYPE(BSET) :: B
    INTEGER    :: I,NUnPEl
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
    ! We need to correct number of alpha/beta electron.
    NUnPEl=G%Multp-1
    G%NAlph=(G%NElec+NUnPEl)/2
    G%NBeta=(G%NElec-NUnPEl)/2
    !write(*,*) 'SetAtomCharges: NElec',G%NElec,' NAlph',G%NAlph,' NBeta',G%NBeta
  END SUBROUTINE SetAtomCharges

  SUBROUTINE UnSetAtomCharges(G,B)
    TYPE(CRDS) :: G
    TYPE(BSET) :: B
    INTEGER    :: I,NUnPEl
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
    ! We need to correct number of alpha/beta electron.
    NUnPEl=G%Multp-1
    G%NAlph=(G%NElec+NUnPEl)/2
    G%NBeta=(G%NElec-NUnPEl)/2
    !write(*,*) 'UnSetAtomCharges: NElec',G%NElec,' NAlph',G%NAlph,' NBeta',G%NBeta
  END SUBROUTINE UnSetAtomCharges
  !
  SUBROUTINE ResetThresholds(cBAS,N,O,G,ReScale)
    TYPE(FileNames)  :: N
    TYPE(Options)    :: O
    TYPE(Geometries) :: G
    TYPE(TOLS)       :: LocalTols
    REAL(DOUBLE)     :: ReScale
    INTEGER          :: cBAS,iCLONE,HDFFileID
    CHARACTER(LEN=2) :: chBAS
    !---------------------------------------------------------------------------!
    LocalTols=O%Thresholds(cBAS)
    LocalTols%Dist=LocalTols%Dist*ReScale
    LocalTols%TwoE=LocalTols%TwoE*ReScale
    !    LocalTols%Trix=LocalTols%Trix*ReScale
    !---------------------------------------------------------------------------!
    chBAS=IntToChar(cBAS)
    HDFFileID=OpenHDF(N%HFile)
    HDF_CurrentID=HDFFileID
    DO iCLONE=1,G%Clones
      HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
      CALL Put(LocalTols,chBAS)
      CALL CloseHDFGroup(HDF_CurrentID)
    ENDDO
    CALL CloseHDF(HDFFileID)
  END SUBROUTINE ResetThresholds
  !-------------------------------------------------------------------
  !
  SUBROUTINE GetRefXYZ(HFileIn,RefXYZ,iCLONE)
    TYPE(FileNames)      :: Nams
    CHARACTER(LEN=*)     :: HFileIn
    INTEGER              :: HDFFileID,iCLONE
    TYPE(DBL_RNK2)       :: RefXYZ
    !
    HDFFileID=OpenHDF(HFileIn)
    HDF_CurrentID= OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
    CALL Get(RefXYZ,'cartesians',Tag_O=TRIM(IntToChar(1)))
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
    HDF_CurrentID= OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
    !
    DO IGeom=IStart,iGEO
      ICount=IGeom-IStart+1
      CALL Get(PBC,Tag_O=TRIM(IntToChar(IGeom)))
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
      DO J=1,NCartS
        SRStruct%D(J,ICount)=Vect%D(J)
      ENDDO
      CALL Get(Aux,'cartesians',Tag_O=TRIM(IntToChar(IGeom)))
      CALL ConvertToXYZRef(Aux%D,RefXYZ,PBC%Dimen, &
           BoxShape_O=PBC%BoxShape%D)
      CALL CartRNK2ToCartRNK1(Vect%D,Aux%D)
      DO J=1,NCartS
        RefStruct%D(J,ICount)=Vect%D(J)
      ENDDO
      CALL Get(Aux,'gradients',Tag_O=TRIM(IntToChar(IGeom)))
      CALL CartRNK2ToCartRNK1(Vect%D,Aux%D)
      DO J=1,NCartS
        RefGrad%D(J,ICount)=Vect%D(J)
      ENDDO
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
