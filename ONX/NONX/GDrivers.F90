SUBROUTINE GDrivers(TBra,TKet,GDrv)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE ONXParameters
  USE ONXMemory
  IMPLICIT NONE
  TYPE(GradD),INTENT(INOUT) :: GDrv
  INTEGER,INTENT(IN)        :: TBra,TKet

  CALL GDLoader(GDrv%GDrv1%I,GDrv%GDrv2%I,GDrv%GDrv3%I,GDrv%GDrv4%I, &
                GDrv%GDrv5%I,TBra,TKet,GDrv%LG1,GDrv%LG2,GDrv%LG3,   &
                GDrv%LG4,GDrv%LG5,GDrv%LenDa,GDrv%LenDc,GDrv%LenDb,  &
                GDrv%LenDn)

END SUBROUTINE GDrivers
