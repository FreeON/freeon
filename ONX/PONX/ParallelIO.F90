!  CALL PWrite(PrmBuf%D,LenF%PrmL,'Prm')
MODULE ParallelIO
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE Parse
  USE ONXParameters
#ifdef PARALLEL
  USE MondoMPI
#endif
   IMPLICIT NONE

   INTERFACE PWrite
     MODULE PROCEDURE PWrite_DBL_VECT, PWrite_DBL_RNK2, &
                      PWrite_INT_VECT
   END INTERFACE

   INTERFACE PRead
     MODULE PROCEDURE PRead_DBL_VECT,  &
                      PRead_INT_VECT
   END INTERFACE

   CONTAINS

   SUBROUTINE PWrite_DBL_VECT(A,Len,Name)
      REAL(DOUBLE),       INTENT(IN) :: A(*)
      INTEGER,            INTENT(IN) :: Len
      CHARACTER(LEN=*),   INTENT(IN) :: Name
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: FName
      INTEGER                        :: PID
      CALL GetEnv('MONDO_SCRATCH',MONDO_SCRATCH)
      CALL GetID(PID)
      FName=TRIM(MONDO_SCRATCH)//'/'//TRIM(IntToChar(PID))//'_'// &
            TRIM(IntToChar(MyID))//'.'//TRIM(Name)
      OPEN(UNIT=70,FILE=FName,FORM="UNFORMATTED")
      WRITE(70) A(1:Len)
      CLOSE(70)
   END SUBROUTINE PWrite_DBL_VECT

   SUBROUTINE PWrite_DBL_RNK2(A,L1,L2,Name)
      TYPE(DBL_RNK2),     INTENT(IN) :: A
      INTEGER,            INTENT(IN) :: L1,L2
      CHARACTER(LEN=*),   INTENT(IN) :: Name
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: FName
      INTEGER                        :: PID
      CALL GetEnv('MONDO_SCRATCH',MONDO_SCRATCH)
      CALL GetID(PID)
      FName=TRIM(MONDO_SCRATCH)//'/'//TRIM(IntToChar(PID))//'_'// &
            TRIM(IntToChar(MyID))//'.'//TRIM(Name)
      OPEN(UNIT=70,FILE=FName,FORM="UNFORMATTED")
!  WRITE(70) A(1:Len)
      CLOSE(70)
   END SUBROUTINE PWrite_DBL_RNK2


   SUBROUTINE PWrite_INT_VECT(A,Len,Name)
      INTEGER,            INTENT(IN) :: A(*)
      INTEGER,            INTENT(IN) :: Len
      CHARACTER(LEN=*),   INTENT(IN) :: Name
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: FName
      INTEGER                        :: PID
      CALL GetEnv('MONDO_SCRATCH',MONDO_SCRATCH)
      CALL GetID(PID)
      FName=TRIM(MONDO_SCRATCH)//'/'//TRIM(IntToChar(PID))//'_'// &
            TRIM(IntToChar(MyID))//'.'//TRIM(Name)
      OPEN(UNIT=70,FILE=FName,FORM="UNFORMATTED")
      WRITE(70) A(1:Len)
      CLOSE(70)
   END SUBROUTINE PWrite_INT_VECT

   SUBROUTINE PRead_DBL_VECT(A,Len,Name)
      REAL(DOUBLE),      INTENT(OUT) :: A(*)
      INTEGER,            INTENT(IN) :: Len
      CHARACTER(LEN=*),   INTENT(IN) :: Name
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: FName
      INTEGER                        :: PID
      CALL GetEnv('MONDO_SCRATCH',MONDO_SCRATCH)
      CALL GetID(PID)
      FName=TRIM(MONDO_SCRATCH)//'/'//TRIM(IntToChar(PID))//'_'// &
            TRIM(IntToChar(MyID))//'.'//TRIM(Name)
      OPEN(UNIT=70,FILE=FName,FORM="UNFORMATTED")
      READ(70) A(1:Len)
      CLOSE(70)
   END SUBROUTINE PRead_DBL_VECT

   SUBROUTINE PRead_INT_VECT(A,Len,Name)
      INTEGER,           INTENT(OUT) :: A(*)
      INTEGER,            INTENT(IN) :: Len
      CHARACTER(LEN=*),   INTENT(IN) :: Name
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: FName
      INTEGER                        :: PID
      CALL GetEnv('MONDO_SCRATCH',MONDO_SCRATCH)
      CALL GetID(PID)
      FName=TRIM(MONDO_SCRATCH)//'/'//TRIM(IntToChar(PID))//'_'// &
            TRIM(IntToChar(MyID))//'.'//TRIM(Name)
      OPEN(UNIT=70,FILE=FName,FORM="UNFORMATTED")
      READ(70) A(1:Len)
      CLOSE(70)
   END SUBROUTINE PRead_INT_VECT

   SUBROUTINE PErase(Name)
      CHARACTER(LEN=*),   INTENT(IN) :: Name
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: FName
      INTEGER                        :: PID
      CALL GetEnv('MONDO_SCRATCH',MONDO_SCRATCH)
      CALL GetID(PID)
      FName=TRIM(MONDO_SCRATCH)//'/'//TRIM(IntToChar(PID))//'_'// &
            TRIM(IntToChar(MyID))//'.'//TRIM(Name)
      OPEN(UNIT=70,FILE=FName,FORM="UNFORMATTED")
      CLOSE(70,STATUS="DELETE")
   END SUBROUTINE PErase


END MODULE

