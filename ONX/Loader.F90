SUBROUTINE Loader(A,N,X)
  USE DerivedTypes
  IMPLICIT NONE
  INTEGER                                 :: N
  REAL(DOUBLE),DIMENSION(N),INTENT(INOUT) :: A
  REAL(DOUBLE)                            :: X

  write(76,*) "N=",N," X=",X
  A(1:N)=X

END SUBROUTINE Loader


