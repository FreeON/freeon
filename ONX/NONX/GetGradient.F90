  SUBROUTINE GetGradient(N,GD,SW,W)
  USE DerivedTypes
  IMPLICIT NONE
  TYPE(GradD)               :: GD
  INTEGER                   :: N,I,K
  REAL(DOUBLE)              :: SW(N,*)      ! Can we fix this? (yes)  
  REAL(DOUBLE)              :: W(N,*)       ! (maybe)
  REAL(DOUBLE)              :: rI4
  INTEGER                   :: I1,I2,I3,I4
  DO K=1,GD%LG4
    I1=GD%GDrv4%I(1,K)
    I2=GD%GDrv4%I(2,K)
    I3=GD%GDrv4%I(3,K)
    I4=GD%GDrv4%I(4,K)
    IF (I3.EQ.0) THEN
      DO I=1,N
        SW(I,I1)=W(I,I2)
      END DO
    ELSE IF (I4.EQ.1) THEN
      DO I=1,N
        SW(I,I1)=W(I,I2)-W(I,I3)
      END DO
    ELSE
      rI4=DBLE(I4)
      DO I=1,N
        SW(I,I1)=W(I,I2)-rI4*W(I,I3)
      END DO
    END IF
  END DO
  END SUBROUTINE GetGradient
