PROGRAM VisDX
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE AtomPairs
  USE RhoUtil
  USE PotUtil 
!
  CHARACTER(LEN=5),PARAMETER  :: Prog='VisDX'
!
  REAL(DOUBLE)                :: Del1
  REAL(DOUBLE),DIMENSION(3)   :: Origin1
  INTEGER                     :: Nx1,Ny1,Nz1
!-------------------------------------------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog)
!---------------------------------------------------------------------------
! Get the geometry
#ifdef MMech
  IF(HasMM()) THEN
    CALL Get(GM,Tag_O='GM_MM'//CurGeom)
  ELSE
    CALL Get(GM,Tag_O=CurGeom)
  ENDIF
#else
  CALL Get(GM,Tag_O=CurGeom)
#endif
  CALL PPrint(GM,TrixFile('xyz',Args,PWD_O=.TRUE.),Geo,'XYZ')
#ifdef PERIODIC
! Get the Outer Cell Set
  CALL Get_CellSet(CS_OUT,'CS_OUT'//CurBase//CurGeom)
#endif
!
! Create the Surfaces
!
  CALL RhoCubed(Args,Del1,Origin1,Nx1,Ny1,Nz1)
  CALL PotCubed(Args,Del1,Origin1,Nx1,Ny1,Nz1)
!
  CALL Delete(Args)
! Shutdown 
  CALL ShutDown(Prog)
END PROGRAM VisDX
