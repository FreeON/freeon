MODULE ParsePeriodic
  USE ControlStructures
  USE PeriodicKeys
CONTAINS
  SUBROUTINE LoadPeriodic(N,G,P)
    TYPE(FileNames)  :: N
    TYPE(Geometries) :: G
    TYPE(Periodics)  :: P
#ifdef PERIODIC

#endif
  END SUBROUTINE LoadPeriodic
END MODULE ParsePeriodic

