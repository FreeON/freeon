MODULE MatrixUtil
  IMPLICIT NONE

  PUBLIC
  INTEGER,PARAMETER :: q1=SELECTED_REAL_KIND(6,30)         ! single precision
  INTEGER,PARAMETER :: q2=SELECTED_REAL_KIND(15,305)       ! double precision
  REAL(q2),PARAMETER :: pi=3.141592653589793238462643_q2

  CONTAINS
!-----------------------------------------------------------------------------------!

  SUBROUTINE matrix_multip(A,C,B,l,m,n)
    INTEGER,INTENT(IN) :: l,m,n
    REAL(q2),INTENT(IN),DIMENSION(l,m) :: A
    REAL(q2),INTENT(IN),DIMENSION(m,n) :: C
    REAL(q2),INTENT(OUT),DIMENSION(l,n) :: B

    INTEGER :: j,k

    WRITE(*,*) A
    WRITE(*,*) C
    WRITE(*,*) B
    pause

    DO j=1,n
      B(:,j)=0.0_q2
      DO k=1,m
        B(:,j)=B(:,j)+c(k,j)*A(:,k)
      END DO
    END DO

  RETURN
  END SUBROUTINE matrix_multip

!-----------------------------------------------------------------------------------!

  SUBROUTINE matrix_vector(A,v,b,n,m)
    INTEGER,INTENT(IN) :: n,m
    REAL(q2),INTENT(IN),DIMENSION(n,m) :: A
    REAL(q2),INTENT(IN),DIMENSION(m) :: v
    REAL(q2),INTENT(OUT),DIMENSION(n) :: b

    INTEGER :: i

    b=0.0_q2
    DO i=1,m
      b=b+v(i)*A(:,i)
    END DO

  RETURN
  END SUBROUTINE matrix_vector

!-----------------------------------------------------------------------------------!

  SUBROUTINE transform_matrix(A,B,n,m)
    INTEGER,INTENT(IN) :: n,m
    REAL(q2),INTENT(IN),DIMENSION(n,m) :: A
    REAL(q2),INTENT(INOUT),DIMENSION(m,n) :: B

    INTEGER :: i,j

    DO i=1,n
      DO j=1,m
        B(j,i)=A(i,j)
      END DO
    END DO

  RETURN
  END SUBROUTINE

!-----------------------------------------------------------------------------------!

END MODULE MatrixUtil
