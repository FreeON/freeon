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
      SUBROUTINE MWait(Sec)
         REAL(DOUBLE) :: Sec,Start
         Start=WallSec()
         WRITE(*,*)' WAITING '
         DO 
           IF(WallSec()-Start>Sec)EXIT
         ENDDO
         WRITE(*,*)' WAITED FOR ',Sec,' Seconds '
      END SUBROUTINE MWait             
!
      FUNCTION CPUSec()
         REAL(DOUBLE) CPUSec
#ifdef PARALLEL
         CPUSec=0.0D0
#else
!         CALL CPU_TIME(CPUSec)
         CPUSec=CPU_Seconds()
#endif
      END FUNCTION CPUSec
      FUNCTION WallSec()
         REAL(DOUBLE) WallSec
         WallSec=Wall_Seconds()
      END FUNCTION WallSec


!REAL  FUNCTION SECOND()
! Calculate the CPU time since the last call to second
! Usage:  t0 = second()
!         your operations
!         tc = second() - t0

!REAL  :: T1
!REAL  :: TARRAY( 2 )
!REAL, EXTERNAL  :: ETIME

! ETIME returns uer time in element 1, system time in element 2
!T1 = ETIME( TARRAY )
!SECOND = TARRAY( 1 )

!RETURN
!END FUNCTION SECOND
END MODULE
