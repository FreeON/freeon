!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 1999, The University of California
!
MODULE GlobalObjects
   USE DerivedTypes
   IMPLICIT NONE
!-------------------------------------------------
!  Global time/flop accumulator
!
   TYPE(TIME), SAVE :: PerfMon
!-------------------------------------------------  
!  Global accumulator of memory statistics 
!
   TYPE(MEMS), SAVE :: MemStats
!-------------------------------------------------  
!  Global debug keys
!
   TYPE(DEBG), SAVE :: PrintFlags
!-------------------------------------------------  
!  Universal BCSR indices
!
   TYPE(INT_VECT), SAVE :: BSiz,OffS
#ifdef PARALLEL
!-------------------------------------------------  
!  Universal || domains
!
   TYPE(INT_VECT), SAVE :: Beg,End
!-----------------------------------------------------------  
!  Off set index for local RowPt -> global GRwPt conversion
!
   TYPE(INT_VECT), SAVE :: OffSt
!-----------------------------------------------------------  
!  Scheduling vector
!
!   TYPE(INT_VECT), SAVE :: ShredSched
#endif
END MODULE
