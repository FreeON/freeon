MODULE Stats
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalObjects
   IMPLICIT NONE
  
   INTERFACE AbsMax 
     MODULE PROCEDURE AbsMax_DBL_VECT, AbsMax_INT_VECT
   END INTERFACE 

   CONTAINS

   FUNCTION AbsMax_DBL_VECT(N,A)
     IMPLICIT NONE
     INTEGER,INTENT(IN)         :: N
     TYPE(DBL_VECT), INTENT(IN) :: A
     REAL(DOUBLE)               :: AbsMax_DBL_VECT
     INTEGER                    :: I
     REAL(DOUBLE)               :: X
     X=0.0D0
     DO I=1,N
       X=MAX(X,ABS(A%D(I)))
     ENDDO
     AbsMax_DBL_VECT=X
   END FUNCTION AbsMax_DBL_VECT

   FUNCTION AbsMax_INT_VECT(N,A)
     IMPLICIT NONE
     INTEGER,INTENT(IN)         :: N
     TYPE(INT_VECT), INTENT(IN) :: A
     REAL(DOUBLE)               :: AbsMax_INT_VECT
     INTEGER                    :: I,J
     J=0
     DO I=1,N
       J=MAX(J,ABS(A%I(I)))
     ENDDO
     AbsMax_INT_VECT=J
   END FUNCTION AbsMax_INT_VECT

END MODULE
