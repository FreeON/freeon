MODULE Mechanics
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE GlobalObjects
  USE InOut
  IMPLICIT NONE
#ifdef MMech
  !  Global parameter for QM/MM
  LOGICAL,DIMENSION(2) :: MechFlag
CONTAINS
  SUBROUTINE InitMMech()
    ! Open HDF file
    CALL OpenHDF(TRIM(InfFile))
    CALL GET(MechFlag(1),'Ctrl_Mechanics1') 
    CALL GET(MechFlag(2),'Ctrl_Mechanics2') 
    CALL CloseHDF()    
  END SUBROUTINE InitMMech

  FUNCTION HasMM() 
    LOGICAL :: HasMM
    HasMM=MechFlag(1)
  END FUNCTION HasMM
  !
  FUNCTION HasQM() 
    LOGICAL :: HasQM
    HasQM=MechFlag(2)
  END FUNCTION HasQM
  !
  FUNCTION MMOnly() 
    LOGICAL :: MMOnly
    MMOnly=MechFlag(1).AND..NOT.MechFlag(2)
  END FUNCTION MMOnly
  !
  FUNCTION QMOnly() 
    LOGICAL :: QMOnly
    QMOnly=MechFlag(2).AND..NOT.MechFlag(1)
  END FUNCTION QMOnly
#endif
END MODULE Mechanics
