MODULE MPIInclude
   IMPLICIT NONE   
#ifdef PARALLEL
   INCLUDE "mpif.h"
#endif
END MODULE
