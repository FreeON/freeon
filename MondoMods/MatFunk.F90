MODULE MatFunk
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE PrettyPrint
!
   INTEGER,SAVE   :: BIGBLOK,BLKLWORK
   TYPE(DBL_VECT) :: BLKWORK,BLKVALS
   TYPE(DBL_RNK2) :: BLKVECT,BLKMAT1,BLKMAT2
   CONTAINS
!
   SUBROUTINE SetDSYEVWork(BIGBLOK_O)
      INTEGER,OPTIONAL,INTENT(IN) :: BIGBLOK_O
      IF(PRESENT(BIGBLOK_O))THEN
         BIGBLOK=BIGBLOK_O
      ELSE
         BIGBLOK=MaxBlkSize
      ENDIF
      BLKLWORK=MAX(1,3*BIGBLOK)
      CALL New(BLKWORK,BLKLWORK)
      CALL New(BLKVALS,BIGBLOK)
      CALL New(BLKVECT,(/BIGBLOK,BIGBLOK/)) 
      CALL New(BLKMAT1,(/BIGBLOK,BIGBLOK/))
      CALL New(BLKMAT2,(/BIGBLOK,BIGBLOK/)) 
   END SUBROUTINE SetDSYEVWork
!
   SUBROUTINE UnSetDSYEVWork()
      CALL Delete(BLKWORK)
      CALL Delete(BLKVALS)
      CALL Delete(BLKVECT) 
      CALL Delete(BLKMAT1) 
      CALL Delete(BLKMAT2) 
   END SUBROUTINE UnSetDSYEVWork
!
   SUBROUTINE FunkOnSqMat(N,Funk,A,FOnA,PrintCond_O)
      INTEGER,                    INTENT(IN)    :: N
      REAL(DOUBLE),DIMENSION(N,N),INTENT(IN)    :: A
      REAL(DOUBLE),DIMENSION(N,N),INTENT(OUT)   :: FOnA
      LOGICAL, OPTIONAL                         :: PrintCond_O
      LOGICAL                                   :: PrintCond
      INTEGER                                   :: I,INFO      
      REAL(DOUBLE), PARAMETER                   :: EigenThreshold=1.D-8
      REAL(DOUBLE)                              :: CondA
      REAL(DOUBLE), EXTERNAL                    :: Funk
      CHARACTER(LEN=DEFAULT_CHR_LEN)            :: String
!----------------------------------------------------------------------------
      IF(PRESENT(PrintCond_O))THEN
         PrintCond=PrintCond_O
      ELSE
         PrintCond=.FALSE.
      ENDIF
!
      BLKVECT%D(1:N,1:N)=A(1:N,1:N)
!
      CALL DSYEV('V','U',N,BLKVECT%D,BIGBLOK,BLKVALS%D,BLKWORK%D,BLKLWORK,INFO)
      IF(INFO/=SUCCEED) &
         CALL Halt('DSYEV hosed in FunkOnSqMat. INFO='//TRIM(IntToChar(INFO)))  
!
      IF(PrintCond)THEN
         CondA=BLKVALS%D(N)/BLKVALS%D(1)
         String="Cond(A) = "//TRIM(DblToShrtChar(CondA))
         WRITE(*,*)TRIM(String)
      ENDIF
!
!     Apply the function to eigenvalues of the matrix, projecting out "zeros"
!
      BLKMAT1%D=Zero
      DO I=1,N
         IF(ABS(BLKVALS%D(I))<EigenThreshold)THEN
            BLKMAT1%D(I,I)=Zero
         ELSE
            BLKMAT1%D(I,I)=Funk(BLKVALS%D(I))
         ENDIF
      ENDDO     
!
!     Reconstruct the functionalized matrix
!
      CALL DGEMM_NNc(N,N,N,One,Zero,BLKVECT%D(1:N,1:N),BLKMAT1%D(1:N,1:N),BLKMAT2%D(1:N,1:N))
      CALL DGEMM_NTc(N,N,N,One,Zero,BLKMAT2%D(1:N,1:N),BLKVECT%D(1:N,1:N),     FOnA(1:N,1:N))
   END SUBROUTINE FunkOnSqMat
!
!
   FUNCTION Inverse(A) RESULT(InvA)
      REAL(DOUBLE), INTENT(IN)  :: A
      REAL(DOUBLE) :: InvA
      InvA=One/A
   END FUNCTION Inverse
!
   FUNCTION SqRoot(A) RESULT(SqRtA)
      REAL(DOUBLE), INTENT(IN)  :: A
      REAL(DOUBLE) :: SqRtA
      SqRtA=SQRT(A)
   END FUNCTION SqRoot
!
END MODULE MatFunk
