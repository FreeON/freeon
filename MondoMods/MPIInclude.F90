!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 1999, The University of California
!
!
!    MPI INCLUDE 
!
MODULE MPIInclude
   IMPLICIT NONE   
#ifdef PARALLEL
   INCLUDE 'MONDO_MPI_INCLUDE.Inc'
#endif
END MODULE
