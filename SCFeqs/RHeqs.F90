!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
PROGRAM RHEqs
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
#ifdef NAG
   USE F90_UNIX_ENV
#endif
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
  TYPE(BCSR)                     :: sP,sF,sX,sTmp1,sTmp2
  TYPE(DBL_RNK2)                 :: X,F,MO,P
  TYPE(DBL_VECT)                 :: EigenV,Work
  TYPE(INT_VECT)                 :: IWork
  TYPE(ARGMT)                    :: Args
  REAL(DOUBLE)                   :: CJK,HOMO,LUMO,dt
  INTEGER                        :: I,J,K,LgN,LWORK,LIWORK,Info
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FMatrix,PMatrix,XFile
  CHARACTER(LEN=5),PARAMETER     :: Prog='RHEqs'
  LOGICAL                        :: Present
#ifdef SCFDEVELOP
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: DevFile
  INTEGER                        :: TmpFmt
  INTEGER                        :: TmpMat
#endif
!--------------------------------------------------------------------
!
!
  CALL StartUp(Args,Prog)
!--------------------------------------------------------------------
!
!
  CALL New(F,(/NBasF,NBasF/))
  FMatrix=TrixFile('F_DIIS',Args,0)
  INQUIRE(FILE=FMatrix,EXIST=Present)
  IF(Present)THEN
     CALL Get(sF,FMatrix)
  ELSE
     CALL Get(sF,TrixFile('OrthoF',Args,0))    ! the orthogonalized Fock matrix
  ENDIF
  CALL SetEq(F,sF)
  CALL Delete(sF) 
!----------------------------------
!
!
  CALL New(EigenV,NBasF)
  CALL SetEq(EigenV,Zero)
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
  CALL DSYEVD('V','U',NBasF,F%D,NBasF,EigenV%D, &
              Work%D,LWORK,IWork%I,LIWORK,Info)
  IF(Info/=SUCCEED) &
  CALL Halt('DSYEVD flaked in RHEqs. INFO='//TRIM(IntToChar(Info)))
  CALL Delete(Work)
  CALL Delete(IWork)
#else
  LWORK=MAX(1,3*NBasF)
  CALL New(Work,LWork)
  CALL DSYEV('V','U',NBasF,F%D,NBasF,EigenV%D,Work%D,LWORK,Info)
  IF(Info/=SUCCEED) &
  CALL Halt('DSYEV flaked in RHEqs. INFO='//TRIM(IntToChar(Info)))
  CALL Delete(Work)
#endif
  HOMO=-1.D10
  LUMO= 1.D10
  DO I=1,NBasF
     IF(EigenV%D(I)<Zero)HOMO=MAX(HOMO,EigenV%D(I))
     IF(EigenV%D(I)>Zero)LUMO=MIN(LUMO,EigenV%D(I))
  ENDDO 
  IF(PrintFlags%Key>=DEBUG_MEDIUM)THEN 
     Mssg=ProcessName(Prog)//'HOMO = '//TRIM(DblToMedmChar(HOMO)) &
                          //', LUMO = '//TRIM(DblToMedmChar(LUMO))
     CALL OpenASCII(OutFile,Out)
     WRITE(Out,*)TRIM(Mssg)
     CLOSE(Out)
  ENDIF
  CALL Put(LUMO-HOMO,'HomoLumoGap')
  CALL Delete(EigenV)
!--------------------------------------------------------------
!
  CALL New(P,(/NBasF,NBasF/))
  CALL DBL_VECT_EQ_DBL_SCLR(NBasF*NBasF,P%D,Zero)
  DO K=1,Nel/2                             ! Closed shell only
     DO J=1,NBasF
        CJK=F%D(J,K)
        DO I=1,NBasF
           P%D(I,J)=P%D(I,J)+F%D(I,K)*CJK  ! P_{ij}=Sum^{N_El/2}_k MO_ik MO_kj
        ENDDO
    ENDDO
  ENDDO    
  CALL SetEq(sX,P)          !  sX=P
  CALL New(sP)              
  CALL Filter(sP,sX)        !  sP=Filter[sX]
  CALL Put(sP,'CurrentOrthoD',CheckPoint_O=.TRUE.)   
  CALL Put(sP,TrixFile('OrthoD',Args,1))
  CALL PChkSum(sP,'OrthoP['//TRIM(NxtCycl)//']',Prog)
  CALL PPrint(sP,'OrthoP['//TRIM(NxtCycl)//']')
  CALL Plot(sP,'OrthoP['//TRIM(NxtCycl)//']')
!---------------------------------------------------
!
!
  XFile=TrixFile('X',Args)
  INQUIRE(FILE=XFile,EXIST=Present)
  IF(Present)THEN
     CALL Get(sX,XFile)                  ! X =S^{-1/2}
     CALL Multiply(sX,sP,sTmp1)          ! T1=S^{-1/2}.P_Orthog
     CALL Multiply(sTmp1,sX,sP)          ! T1=S^{-1/2}.P_Orthog.S^{-1/2}
  ELSE
     CALL Get(sX,TRIM(SCFName)//'.Z')    ! X=Z=L^{-1}
     CALL Multiply(sX,sP,sTmp1)          ! T1=Z.P_Orthog
     CALL Get(sX,TRIM(SCFName)//'.ZT')   ! X =Z^t=L^{-t}           
     CALL Multiply(sTmp1,sX,sP)          ! F=Z.P_Orthog.Z^t
  ENDIF
  CALL Filter(sTmp1,sP)                  ! T1 =P_AO=Filter[Z.P_Orthog.Z^t]
!
  CALL Put(sTmp1,TrixFile('D',Args,1))
  CALL PChkSum(sTmp1,'P['//TRIM(NxtCycl)//']',Prog)
  CALL PPrint(sTmp1,'P['//TRIM(NxtCycl)//']')
  CALL Plot(sTmp1,'P['//TRIM(NxtCycl)//']')
#ifdef SCFDEVELOP
!----------------------------------------------------
! Read in and output intermediates for the development
! of new methods for solving the SCF equations
!
  CALL GetEnv('MONDO_HOME',DevFile)
  DevFile=TRIM(DevFile)//'/SCFeqs/MMA/SCFDev/scf.in'
  CALL OpenASCII(DevFile,35,NewFile_O=.TRUE.)
  PrintFlags%Fmt=DEBUG_MMASTYLE
  PrintFlags%Mat=DEBUG_MATRICES
  WRITE(35,*)' (* All matrices are in a non-orthogonal rep *)'
  WRITE(35,*)' NBasF  = ',NBasF,';'
  WRITE(35,*)' NEl    = ',NEl,';'
  WRITE(35,*)' NAtoms = ',NAtoms,';'  
  CALL Get(sX,TrixFile('F',Args,0))
  CALL PPrint(sX,'F',FileName_O=DevFile,Unit_O=35)   
  CALL Get(sX,TrixFile('S',Args))
  CALL PPrint(sX,'S',FileName_O=DevFile,Unit_O=35)   
  IF(Present)THEN
     CALL Get(sX,XFile)                  
     CALL PPrint(sX,'X',FileName_O=DevFile,Unit_O=35)   
  ELSE
     CALL Get(sX,TrixFile('Z',Args))
     CALL PPrint(sX,'Z',FileName_O=DevFile,Unit_O=35)   
     CALL Get(sX,TrixFile('ZT',Args))
     CALL PPrint(sX,'ZT',FileName_O=DevFile,Unit_O=35)   
  ENDIF
  CALL get(sX,TrixFile('D',Args,1))
  CALL PPrint(sX,'P',FileName_O=DevFile,Unit_O=35)   

  CLOSE(UNIT=35,STATUS='KEEP')
#endif
!----------------------------------------------
  CALL Delete(sX)
  CALL Delete(sP)
  CALL ShutDown(Prog)
END PROGRAM ! RHEqs


