MODULE DGEMMS
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE PrettyPrint
!        INTERFACE 
!           SUBROUTINE DGEMM_NN(M,K,N,Beta,A,B,C)
!              IMPLICIT NONE
!              INTEGER,                  INTENT(IN)    :: M,K,N
!              INTEGER,PARAMETER                       :: DOUBLE=KIND(0.D0)  
!              REAL(DOUBLE),             INTENT(IN)    :: Beta
!              REAL(DOUBLE),DIMENSION(:,:),INTENT(IN)    :: A,B
!              REAL(DOUBLE),DIMENSION(:,:),INTENT(INOUT) :: C
!           END SUBROUTINE DGEMM_NN
!        END INTERFACE
  CONTAINS
   SUBROUTINE DGEMM__NN(L,M,N,Beta,A,B,C,Print_O)
      INTEGER                      :: I,J,K,L,M,N
      REAL(DOUBLE), DIMENSION(L,M) :: A 
      REAL(DOUBLE), DIMENSION(M,N) :: B 
      REAL(DOUBLE), DIMENSION(L,N),INTENT(INOUT) :: C 
      REAL(DOUBLE)                 :: Beta,CIJ
      LOGICAL,OPTIONAL :: Print_O
      C=C+Beta*MATMUL(A,B)
      RETURN
   END SUBROUTINE DGEMM__NN

   SUBROUTINE DGEMM__TN(L,M,N,Beta,A,B,C)
      INTEGER :: I,J,K,L,M,N
!     [MxN]<-[MxL].[LxN]
      REAL(DOUBLE), DIMENSION(L,M) :: A 
      REAL(DOUBLE), DIMENSION(L,N) :: B 
      REAL(DOUBLE), DIMENSION(M,N) :: C 
      REAL(DOUBLE)                 :: Beta,CIJ
      C=C+Beta*MATMUL(TRANSPOSE(A),B)
      RETURN
   END SUBROUTINE DGEMM__TN

   SUBROUTINE DGEMM__NT(L,M,N,Beta,A,B,C)
      INTEGER :: I,J,K,L,M,N
!     [LxN]<-[LxM].[MxN]
      REAL(DOUBLE), DIMENSION(L,N) :: A 
      REAL(DOUBLE), DIMENSION(M,N) :: B 
      REAL(DOUBLE), DIMENSION(L,M) :: C 
      REAL(DOUBLE)                 :: Beta,CIJ
      C=C+Beta*MATMUL(A,TRANSPOSE(B))
      RETURN
   END SUBROUTINE DGEMM__NT

   SUBROUTINE DGEMM__TT(L,M,N,Beta,A,B,C)
      INTEGER :: I,J,K,L,M,N
!     [NxM]<-[NxL].[LxM]
      REAL(DOUBLE), DIMENSION(L,N) :: A 

      REAL(DOUBLE), DIMENSION(M,L) :: B 
      REAL(DOUBLE), DIMENSION(N,M) :: C 
      REAL(DOUBLE)                 :: Beta,CIJ
      C=C+Beta*MATMUL(TRANSPOSE(A),TRANSPOSE(B))
      RETURN
   END SUBROUTINE DGEMM__TT

END MODULE
