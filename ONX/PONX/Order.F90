MODULE Order
   USE DerivedTypes
   USE GlobalCharacters
   USE GlobalScalars
   USE GlobalObjects
   USE ProcessControl
   USE MemMan
   IMPLICIT NONE
   INTERFACE Sort
      MODULE PROCEDURE Sort_DBL_INT,Sort_INT_INT, Sort_INT_VECT
   END INTERFACE

   CONTAINS

     SUBROUTINE Sort_DBL_INT(X,Y,N_O,Ordr_O)
        TYPE(DBL_VECT), INTENT(INOUT) :: X
        TYPE(INT_VECT), INTENT(INOUT) :: Y
        INTEGER,OPTIONAL,INTENT(IN)   :: N_O,Ordr_O 
        TYPE(DBCSR)                   :: A
        INTEGER                       :: N,Ordr 
        Ordr=-2
        N=MIN(SIZE(X%D),SIZE(Y%I))
        IF(PRESENT(N_O))THEN
           IF(N_O>N)CALL Halt(' Dimensioning off in DblSort ')
           N=N_O
        ENDIF
        IF(PRESENT(Ordr_O))Ordr=Ordr_O
        CALL DblIntSort77(N,X%D,Y%I,Ordr)
    END SUBROUTINE Sort_DBL_INT


     SUBROUTINE Sort_INT_INT(X,Y,N_O,Ordr_O)
        TYPE(INT_VECT), INTENT(INOUT) :: X
        TYPE(INT_VECT), INTENT(INOUT) :: Y
        INTEGER,OPTIONAL,INTENT(IN)   :: N_O,Ordr_O 
        INTEGER                       :: N,Ordr 
        Ordr=-1
        N=MIN(SIZE(X%I),SIZE(Y%I))
        IF(PRESENT(N_O))THEN
           IF(N_O>N)CALL Halt(' Dimensioning off in DblSort ')
           N=N_O
        ENDIF
        IF(PRESENT(Ordr_O))Ordr=Ordr_O
        CALL IntIntSort77(N,X%I,Y%I,Ordr)
    END SUBROUTINE Sort_INT_INT


     SUBROUTINE Sort_INT_VECT(X,N_O,Ordr_O)
        TYPE(INT_VECT), INTENT(INOUT) :: X
        INTEGER,OPTIONAL,INTENT(IN)   :: N_O,Ordr_O 
        INTEGER                       :: N,Ordr 
        Ordr=-1
        N=SIZE(X%I)
        IF(PRESENT(N_O))THEN
           IF(N_O>N)CALL Halt(' Dimensioning off in DblSort ')
           N=N_O
        ENDIF
        IF(PRESENT(Ordr_O))Ordr=Ordr_O
        CALL IntSort77(N,X%I,Ordr)
    END SUBROUTINE Sort_INT_VECT


END MODULE
