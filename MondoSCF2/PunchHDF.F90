MODULE PunchHDF
  USE InOut
  USE MemMan
  USE PrettyPrint
  USE ControlStructures
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
  SUBROUTINE MPIsArchive(N,G,M)
    TYPE(FileNames)  :: N
    TYPE(Geometries) :: G
    TYPE(Parallel)   :: M
#ifdef PARALLEL
    TYPE(INT_VECT)   :: ST
#endif
    !---------------------------------------------------------------------------!
#ifdef PARALLEL
    HDF_CurrentID=OpenHDF(N%HFile)
    CALL New(ST,2)
    ST%I=(/M%NSpace,G%Clones/)
    CALL Put(ST,'SpaceTime')
    CALL Delete(ST)
    CALL CloseHDF(HDF_CurrentID)
#endif
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
  !
  !==============================================================================
  SUBROUTINE InitClones(N,G)
    TYPE(FileNames)  :: N
    TYPE(Geometries) :: G    
    INTEGER          :: iGEO,iCLONE,HDFFileID
    CHARACTER(LEN=2) :: chGEO
    !---------------------------------------------------------------------------!
    chGEO=IntToChar(iGEO)
    HDFFileID=OpenHDF(N%HFile)
    DO iCLONE=1,G%Clones
       HDF_CurrentID=InitHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
       CALL CloseHDFGroup(HDF_CurrentID)
    ENDDO
    CALL CloseHDF(HDFFileID)
  END SUBROUTINE InitClones

  SUBROUTINE GeomArchive(iGEO,N,G)
    TYPE(FileNames)  :: N
    TYPE(Geometries) :: G    
    INTEGER          :: iGEO,iCLONE,HDFFileID
    CHARACTER(LEN=2) :: chGEO
    !---------------------------------------------------------------------------!
    chGEO=IntToChar(iGEO)
    HDFFileID=OpenHDF(N%HFile)
    DO iCLONE=1,G%Clones
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
!       CALL Print_CRDS(G%Clone(iCLONE))
       CALL Put(G%Clone(iCLONE),chGEO)
       CALL CloseHDFGroup(HDF_CurrentID)
    ENDDO
    CALL CloseHDF(HDFFileID)
  END SUBROUTINE GeomArchive
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
       CALL Put(M%Beg(iCLONE,cBAS),'beg',Tag_O=chBAS)
       CALL Put(M%End(iCLONE,cBAS),'end',Tag_O=chBAS)
       CALL Put(M%GLO(iCLONE,cBAS),'dbcsroffsets',Tag_O=chBAS)
       CALL Put(B%BSets(iCLONE,cBAS),Tag_O=chBAS)
       CALL Put(B%BSiz(iCLONE,cBAS),'atsiz',Tag_O=chBAS)
       CALL Put(B%OffS(iCLONE,cBAS),'atoff',Tag_O=chBAS)
       CALL Put(B%DExpt(iCLONE,cBAS),'dexpt',Tag_O=chBAS)
       CALL Put(B%Lndex(iCLONE,cBAS),'lndex',Tag_O=chBAS)
       CALL CloseHDFGroup(HDF_CurrentID)
    ENDDO
    CALL CloseHDF(HDFFileID)
  END SUBROUTINE BSetArchive

END MODULE PunchHDF
