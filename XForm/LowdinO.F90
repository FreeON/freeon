PROGRAM LowdinO
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE SetXYZ
  USE LinAlg
  IMPLICIT NONE
#ifdef DSYEVD  
  INTERFACE DSYEVD
     SUBROUTINE DSYEVD(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,IWORK,LIWORK,INFO)
         USE GlobalScalars
         CHARACTER(LEN=1), INTENT(IN)    :: JOBZ, UPLO
         INTEGER,          INTENT(IN)    :: LDA, LIWORK, LWORK, N
         INTEGER,          INTENT(OUT)   :: INFO
         INTEGER,          INTENT(OUT)   :: IWORK(*)
         REAL(DOUBLE),     INTENT(INOUT) :: A(LDA,*)
         REAL(DOUBLE),     INTENT(OUT)   :: W(*)
         REAL(DOUBLE),     INTENT(OUT)   :: WORK(*)
     END SUBROUTINE DSYEVD
  END INTERFACE
#else
  INTERFACE DSYEV
     SUBROUTINE DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
         USE GlobalScalars
         CHARACTER(LEN=1), INTENT(IN)    :: JOBZ, UPLO
         INTEGER,          INTENT(IN)    :: LDA,  LWORK, N
         INTEGER,          INTENT(OUT)   :: INFO
         REAL(DOUBLE),     INTENT(INOUT) :: A(LDA,*)
         REAL(DOUBLE),     INTENT(OUT)   :: W(*)
         REAL(DOUBLE),     INTENT(OUT)   :: WORK(*)
     END SUBROUTINE DSYEV
  END INTERFACE
#endif
  TYPE(BCSR)                     :: S,X
  TYPE(DBL_RNK2)                 :: Vectors,Tmp1,Tmp2
  TYPE(DBL_VECT)                 :: Values,Work
  TYPE(INT_VECT)                 :: IWork
  TYPE(ARGMT)                    :: Args
  INTEGER                        :: I,J,K,LgN,LWORK,LIWORK,Info,Status
  CHARACTER(LEN=7),PARAMETER     :: Prog='LowdinO'
  REAL(DOUBLE)                   :: Chk
!--------------------------------------------------------------------
  CALL StartUp(Args,Prog)
  CALL New(Values,NBasF)
  CALL New(Vectors,(/NBasF,NBasF/))
  CALL Get(S,TrixFile('S',Args))
  CALL SetEq(Vectors,S)
  CALL Delete(S)
!--------------------------------------------------------------------
!
!
#ifdef DSYEVD
  DO K=4,10000
     IF(2**K>=NBasF)THEN
        LgN=K
        EXIT
     ENDIF   
  ENDDO
  LWORK=2*(1+5*NBasF+2*NBasF*LgN+3*NBasF**2)
  LIWORK=2*(2+5*NBasF)
  CALL New(Work,LWork)
  CALL New(IWork,LIWork)
  CALL DSYEVD('V','U',NBasF,Vectors%D,NBasF,Values%D, &
              Work%D,LWORK,IWork%I,LIWORK,Info)
  IF(Info/=SUCCEED) &
  CALL Halt('DSYEVD flaked in LowdinO. INFO='//TRIM(IntToChar(Info)))
  CALL Delete(Work)
  CALL Delete(IWork)
#else
  LWORK=MAX(1,3*NBasF)
  CALL New(Work,LWork)
  CALL DSYEV('V','U',NBasF,Vectors%D,NBasF,Values%D,Work%D,LWORK,Info)
  IF(Info/=SUCCEED) &
  CALL Halt('DSYEV flaked in LowdinO. INFO='//TRIM(IntToChar(Info)))
  CALL Delete(Work)
#endif
!--------------------------------------------------------------------
  CALL New(Tmp2,(/NBasF,NBasF/))
  CALL New(Tmp1,(/NBasF,NBasF/))
!--------------------------------------------------------------------
  Tmp1%D=Zero
  DO I=1,NBasF
     IF(Values%D(I)>1.D-12)THEN
        Tmp1%D(I,I)=One/SQRT(Values%D(I))
     ELSE
        Tmp1%D(I,I)=Zero
     ENDIF
  ENDDO     
  CALL DGEMM('N','N',NBasF,NBasF,NBasF,One,Vectors%D, &
             NBasF,Tmp1%D,NBasF,Zero,Tmp2%D,NBasF)
  CALL DGEMM('N','T',NBasF,NBasF,NBasF,One,Tmp2%D,    &
             NBasF,Vectors%D,NBasF,Zero,Tmp1%D,NBasF)
!--------------------------------------------------------------------
  CALL Delete(Values)
  CALL Delete(Vectors)
  CALL Delete(Tmp2)
!--------------------------------------------------------------------
  CALL SetEq(S,Tmp1)    
  CALL Delete(Tmp1)
  CALL Filter(X,S)      
!--------------------------------------------------------------------
  CALL Put(X,TrixFile('X',Args)) 
  CALL PChkSum(X,'X',Prog)
  CALL PPrint(X,'X')
  CALL Plot(X,'X')
!--------------------------------------------------------------------
  CALL Delete(S)
  CALL Delete(X)
  CALL ShutDown(Prog)
!
END PROGRAM LowdinO


#ifdef BROKEN_CODE

!
!!$  DO I=1,NBasF
!!$     DO J=1,NBasF
!!$        SUM = Zero
!!$        DO K=1,NBasF
!!$           IF(ABS(Values%D(K)) .GT. 1.D-10) THEN
!!$              SUM = SUM + Vectors%D(I,K)*Vectors%D(J,K)/SQRT(Values%D(K))
!!$           ENDIF
!!$        ENDDO
!!$        Tmp1%D(I,J) = SUM
!!$     ENDDO
!!$  ENDDO
!------------------------------------------------------------------------
!
!
  DO I=1,NBasF
     DO J=1,NBasF
        IF(ABS(Values%D(J)) .GT. 1.D-8) THEN
           Tmp1%D(I,J) = Vectors%D(I,J)/SQRT(Values%D(J))
           Tmp2%D(I,J) = Vectors%D(J,I)
        ENDIF
     ENDDO
  ENDDO
  CALL DGEMM('N','N',NBasF,NBasF,NBasF,One,Tmp1%D,NBasF,Tmp2%D,NBasF,Zero,Vectors%D,NBasF)
!
!
! Test for Inverse
!    
  CALL Get(S,TrixFile('S',Args))
  CALL SetEq(Tmp1,S)
  CALL Delete(S)
!
  DO I=1,NBasF
     DO J=1,NBasF
        SUM = Zero
        DO K=1,NBasF
           SUM = SUM + Vectors%D(I,K)*Vectors%D(K,J)
        ENDDO
        Tmp2%D(I,J) = SUM
     ENDDO
  ENDDO
!           
  Error = Zero
  DO I=1,NBasF
     DO J=1,NBasF
        SUM = Zero
        DO K=1,NBasF
           SUM = SUM + Tmp1%D(I,K)*Tmp2%D(K,J)
        ENDDO
        IF(I==J) THEN
           Error = Error + ABS((One-SUM))
        ELSE
           Error = Error + ABS(SUM)
        ENDIF
     ENDDO
  ENDDO
  WRITE(*,*) 'Num of Basis Fun   = ',NBasF
  WRITE(*,*) 'Lowest  Eigenvalue = ',Values%D(1)
  WRITE(*,*) 'Highest Eigenvalue = ',Values%D(NBasF)
  WRITE(*,*) 'Condition Number   = ',Values%D(NBasF)/Values%D(1)
  WRITE(*,*) 'Error in Inverse   = ',Error

#endif
