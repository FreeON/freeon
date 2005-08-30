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
  TYPE(BCSR)                     :: sP,sF,sX,sTmp1,sTmp2
  TYPE(DBL_RNK2)                 :: X,F,MO,P
  TYPE(DBL_VECT)                 :: EigenV,Work
  TYPE(INT_VECT)                 :: IWork
  TYPE(ARGMT)                    :: Args
  REAL(DOUBLE)                   :: CJK,HOMO,LUMO,dt
  REAL(DOUBLE)                   :: Mu,Entrop,Ek,Z,A1,A2,H1,H3,H4,Fk,Sigma,Dum
  INTEGER                        :: I,J,K,LgN,LWORK,LIWORK,Info
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FMatrix,PMatrix,XFile,smearing
  CHARACTER(LEN=5),PARAMETER     :: Prog='RHEqs'
  LOGICAL                        :: Present,DensityArchive
!--------------------------------------------------------------------
!
!
  CALL StartUp(Args,Prog,Serial_O=.TRUE.)
!--------------------------------------------------------------------
! Initialize and parse some variables.
!
  CALL OpenASCII(InpFile,Inp)
  Smearing='NoSmearing'
  Sigma=0.002D0
  IF(OptKeyQ(Inp,'Smearing','MP')) Smearing='Methfessel-Paxton'
  IF(OptDblQ(Inp,'SmearingValue',Dum)) Sigma=Dum
  CLOSE(Inp)
!--------------------------------------------------------------------
!
!
  CALL New(F,(/NBasF,NBasF/))
  CALL New(sF)
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
  LWORK=MAX(1,3*NBasF+10)
  CALL New(Work,LWork)
! 
  CALL DSYEV('V','U',NBasF,F%D,NBasF,EigenV%D,Work%D,LWORK,Info)
  IF(Info/=SUCCEED)CALL Halt('DSYEV flaked in RHEqs. INFO='//TRIM(IntToChar(Info)))
!
  HOMO=EigenV%D(NEl/2)
  LUMO=EigenV%D(NEl/2+1)
  Mu=(HOMO+LUMO)*0.5D0
!
  IF(PrintFlags%Key>=DEBUG_MEDIUM)THEN 
     Mssg=ProcessName(Prog)//'HOMO = '//TRIM(DblToMedmChar(HOMO)) &
                         //', LUMO = '//TRIM(DblToMedmChar(LUMO))
     CALL OpenASCII(OutFile,Out)
     CALL PrintProtectL(Out)
     WRITE(Out,*)TRIM(Mssg)
     CALL PrintProtectR(Out)
     CLOSE(Out)
  ENDIF
  CALL Put(HOMO-LUMO,'HomoLumoGap')
!
  CALL Delete(Work)
!--------------------------------------------------------------
! Make a new closed shell, orthogonal density matrix
!
  CALL New(P,(/NBasF,NBasF/))
  P%D=Zero
!  CALL DBL_VECT_EQ_DBL_SCLR(NBasF*NBasF,P%D,Zero)
!
  SELECT CASE(Smearing)
  CASE('Methfessel-Paxton')
     Entrop=0D0
     DO K=1,NBasF
        ! PRB 40, 3616, 1989.
        ! Second order smearing N=2.
        Ek=EigenV%D(K)
        Z=(Ek-Mu)/Sigma
        A1=-1D0/( 4D0*SqrtPi)
        A2= 1D0/(32D0*SqrtPi)
        H1=2D0*Z
        H3=Z*(8D0*Z**2-12D0)
        H4=Z**2*(16D0*Z**2-48D0)+12D0
        ! Fractional occupation.
        Fk=0.5D0*(1D0-ERF(Z))+EXP(-Z**2)*(A1*H1+A2*H3)
        ! Entropic correction to the energy (Comp. Mat. Sci. 6, 15, 1996).
        Entrop=Entrop+Sigma*0.5D0*A2*H4*EXP(-Z**2)
        !write(*,*) Fk,Entrop
        DO J=1,NBasF
           CJK=F%D(J,K)*Fk
           DO I=1,NBasF
              P%D(I,J)=P%D(I,J)+F%D(I,K)*CJK
           ENDDO
        ENDDO
     ENDDO
     IF(PrintFlags%Key>=DEBUG_MEDIUM)THEN 
        Mssg=ProcessName(Prog)//'Sigma = '//TRIM(DblToShrtChar(Sigma)) &
             &               //', Entropic correction per atom = ' &
             &               //TRIM(DblToShrtChar(Entrop/DBLE(NAtoms)))
        CALL OpenASCII(OutFile,Out)
        CALL PrintProtectL(Out)
        WRITE(Out,*)TRIM(Mssg)
        CALL PrintProtectR(Out)
        CLOSE(Out)
     ENDIF
  CASE('NoSmearing')
     DO K=1,Nel/2                             ! Closed shell only
        DO J=1,NBasF
           CJK=F%D(J,K)
           DO I=1,NBasF
              P%D(I,J)=P%D(I,J)+F%D(I,K)*CJK  ! P_{ij}=Sum^{N_El/2}_k MO_ik MO_kj
           ENDDO
        ENDDO
     ENDDO
  CASE DEFAULT
     CALL Halt('RHEqs: Doesn''t regonize this Smearing <'//TRIM(Smearing)//'>')
  END SELECT
!
  CALL Delete(EigenV)
! 
  CALL SetEq(sX,P)          !  sX=P
  CALL New(sP)              
  CALL Filter(sP,sX)        !  sP=Filter[sX]
  CALL Put(sP,TrixFile('OrthoD',Args,1))
  CALL PChkSum(sP,'OrthoP['//TRIM(NxtCycl)//']',Prog)
  CALL PPrint(sP,'OrthoP['//TRIM(NxtCycl)//']')
  CALL Plot(sP,'OrthoP['//TRIM(NxtCycl)//']')
  ! Transform to non-orthogonal rep
  XFile=TrixFile('X',Args)
  INQUIRE(FILE=XFile,EXIST=Present)
  IF(Present)THEN
     CALL Get(sX,XFile)                  ! X =S^{-1/2}
     CALL Multiply(sX,sP,sTmp1)          ! T1=S^{-1/2}.P_Orthog
     CALL Multiply(sTmp1,sX,sP)          ! T1=S^{-1/2}.P_Orthog.S^{-1/2}
  ELSE
     CALL Get(sX,TrixFile('Z',Args))     ! X=Z=L^{-1}
     CALL Multiply(sX,sP,sTmp1)          ! T1=Z.P_Orthog
     CALL Get(sX,TrixFile('ZT',Args))    ! X =Z^t=L^{-t}           
     CALL Multiply(sTmp1,sX,sP)          ! F=Z.P_Orthog.Z^t
  ENDIF
  CALL Filter(sTmp1,sP)                  ! T1 =P_AO=Filter[Z.P_Orthog.Z^t]
  ! Archive the AO-DM
  CALL Get(DensityArchive,'ArchiveDensity')
  IF(DensityArchive) &
     CALL Put(sTmp1,'CurrentDM',CheckPoint_O=.TRUE.)
  CALL Put(sTmp1,TrixFile('D',Args,1))
  CALL PChkSum(sTmp1,'P['//TRIM(NxtCycl)//']',Prog)
  CALL PPrint(sTmp1,'P['//TRIM(NxtCycl)//']')
  CALL Plot(sTmp1,'P['//TRIM(NxtCycl)//']')
!----------------------------------------------
  CALL Delete(sX)
  CALL Delete(sP)
  CALL Delete(sTmp1)
  CALL ShutDown(Prog)
END PROGRAM  RHEqs


