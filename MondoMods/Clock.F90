MODULE Clock
   USE DerivedTypes
#ifdef PARALLEL
   USE MondoMPI
#endif
   IMPLICIT NONE
   INTERFACE 
      FUNCTION CPU_Seconds()
         USE DerivedTypes
         REAL(DOUBLE) CPU_Seconds
      END FUNCTION CPU_Seconds

      FUNCTION Wall_Seconds()
         USE DerivedTypes
         REAL(DOUBLE) Wall_Seconds
      END FUNCTION Wall_Seconds
   END INTERFACE
   CONTAINS
      FUNCTION CPUSec()
         REAL(DOUBLE) CPUSec
#ifdef PARALLEL
         CPUSec=0.0D0
#else
         CPUSec=CPU_Seconds()
#endif
      END FUNCTION CPUSec
      FUNCTION WallSec()
         REAL(DOUBLE) WallSec
         WallSec=Wall_Seconds()
      END FUNCTION WallSec
END MODULE
