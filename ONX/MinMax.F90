MODULE MinMax
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  IMPLICIT NONE
  
  INTERFACE GetMax
    MODULE PROCEDURE GetMax_INT_VECT, GetMax_DBL_VECT
  END INTERFACE

  INTERFACE GetMin
    MODULE PROCEDURE GetMin_INT_VECT, GetMin_DBL_VECT
  END INTERFACE

  CONTAINS

  FUNCTION GetMax_INT_VECT(N,A,M_O)
    INTEGER, INTENT(IN)           :: N
    TYPE(INT_VECT), INTENT(IN)    :: A
    INTEGER, OPTIONAL, INTENT(IN) :: M_O
    INTEGER                       :: GetMax_INT_VECT,I,J,M
    M=0; J=1; IF(PRESENT(M_O)) J=M_O
    DO I=J,N
      M=MAX(M,A%I(I))
    END DO
    GetMax_INT_VECT=M
  END FUNCTION GetMax_INT_VECT

  FUNCTION GetMax_DBL_VECT(N,A,M_O)
    INTEGER, INTENT(IN)           :: N
    TYPE(DBL_VECT), INTENT(IN)    :: A
    INTEGER, OPTIONAL, INTENT(IN) :: M_O
    REAL(DOUBLE)                  :: GetMax_DBL_VECT,M
    INTEGER                       :: I,J
    M=0.0D0; J=1; IF(PRESENT(M_O)) J=M_O
    DO I=J,N
      M=MAX(M,A%D(I))
    END DO
    GetMax_DBL_VECT=M
  END FUNCTION GetMax_DBL_VECT

  FUNCTION GetMin_INT_VECT(N,A,M_O)
    INTEGER, INTENT(IN)           :: N
    TYPE(INT_VECT), INTENT(IN)    :: A
    INTEGER, OPTIONAL, INTENT(IN) :: M_O
    INTEGER                       :: GetMin_INT_VECT,I,J,M
    M=0; J=1; IF(PRESENT(M_O)) J=M_O
    DO I=J,N
      M=MIN(M,A%I(I))
    END DO
    GetMin_INT_VECT=M
  END FUNCTION GetMin_INT_VECT

  FUNCTION GetMin_DBL_VECT(N,A,M_O)
    INTEGER, INTENT(IN)           :: N
    TYPE(DBL_VECT), INTENT(IN)    :: A
    INTEGER, OPTIONAL, INTENT(IN) :: M_O
    REAL(DOUBLE)                  :: GetMin_DBL_VECT,M
    INTEGER                       :: I,J
    M=0.0D0; J=1; IF(PRESENT(M_O)) J=M_O
    DO I=J,N
      M=MIN(M,A%D(I))
    END DO
    GetMin_DBL_VECT=M
  END FUNCTION GetMin_DBL_VECT


END MODULE
