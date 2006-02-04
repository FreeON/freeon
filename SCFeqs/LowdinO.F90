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
  USE DenMatMethods,ONLY:MDiag_DSYEVD
  IMPLICIT NONE
  TYPE(BCSR)                     :: S,X
  TYPE(DBL_RNK2)                 :: Vectors,Tmp1,Tmp2
  TYPE(DBL_VECT)                 :: Values
  TYPE(ARGMT)                    :: Args
  INTEGER                        :: I,J,K,LgN,Info,Status,ISmall
  REAL(DOUBLE)                   :: Chk,CondS,OverlapEThresh,Scale
  CHARACTER(LEN=7),PARAMETER     :: Prog='LowdinO'
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
  EXTERNAL :: DGEMM_NT
!--------------------------------------------------------------------
  CALL StartUp(Args,Prog)
  ! Check for thresholds overide 
  CALL OpenASCII(InpFile,Inp)         
  IF(OptDblQ(Inp,Prog,OverlapEThresh))THEN
     Mssg=TRIM(ProcessName(Prog))//' OverlapEigenThreshold  = '  &
          //TRIM(DblToShrtChar(OverlapEThresh))
     CALL OpenASCII(OutFile,Out)         
     WRITE(Out,*)TRIM(Mssg)
     CLOSE(Out)
  ELSE
     OverlapEThresh=1.D-4
  ENDIF
  CLOSE(Inp)
  CALL New(Values,NBasF)
  CALL New(Vectors,(/NBasF,NBasF/))
  CALL Get(S,TrixFile('S',Args))
  CALL SetEq(Vectors,S)
  CALL Delete(S)
!--------------------------------------------------------------------
! Diag S
!
  CALL MDiag_DSYEVD(Vectors,NBasF,Values,0)
!--------------------------------------------------------------------
!
!
  IF(Values%D(1)<Zero)CALL Halt(' S matrix is not pos def')
  CondS=Values%D(NBasF)/Values%D(1)
  IF(CondS>1D4)CALL Warn('Illconditioning detected in MakeS: Cond(S)='//TRIM(DblToShrtChar(CondS)))
!--------------------------------------------------------------------
!
!
  CALL New(Tmp2,(/NBasF,NBasF/))
  CALL New(Tmp1,(/NBasF,NBasF/))
  !
  ISmall=0
  CALL DCOPY(NBasF**2,Vectors%D(1,1),1,Tmp2%D(1,1),1)
  DO I=1,NBasF
     IF(Values%D(I)>OverlapEThresh)THEN
        Scale=1D0/SQRT(Values%D(I))
        CALL DSCAL(NBasF,Scale,Tmp2%D(1,I),1)
     ELSE
        CALL DSCAL(NBasF,  0D0,Tmp2%D(1,I),1)
        ISmall=ISmall+1
     ENDIF
  ENDDO
  !
  CALL DGEMM_NT(NBasF,NBasF,NBasF,0D0,Tmp2%D(1,1),Vectors%D(1,1),Tmp1%D(1,1))
  !CALL DGEMM('N','T',NBasF,NBasF,NBasF,1D0,Tmp2%D(1,1), &
  !     &     NBasF,Vectors%D(1,1),NBasF,0D0,Tmp1%D(1,1),NBasF)
  !
  IF(ISmall.NE.0) CALL Warn('Removed '//TRIM(IntToChar(ISmall))//' eigenvalue(s) smaller than ' &
       & //TRIM(DblToShrtChar(OverlapEThresh)))
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
