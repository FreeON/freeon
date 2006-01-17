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
  TYPE(BCSR)                     :: sP,sF,sX,sTmp1,sTmp2
  TYPE(DBL_RNK2)                 :: X,F,MO,P    ,test1,test2,test3
  TYPE(DBL_VECT)                 :: EigenV,Work
  TYPE(INT_VECT)                 :: IWork
  TYPE(ARGMT)                    :: Args
  REAL(DOUBLE)                   :: CJK,HOMO,LUMO,dt
  REAL(DOUBLE)                   :: Mu,Entrop,Ek,Z,A1,A2,H1,H3,H4,Fk,Sigma,Dum
  INTEGER                        :: I,J,K,LgN,LWORK,LIWORK,Info,NRow,NCol,NSMat
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
  CALL New(sF)
  FMatrix=TrixFile('F_DIIS',Args,0)
  INQUIRE(FILE=FMatrix,EXIST=Present)
  IF(Present)THEN
     CALL Get(sF,FMatrix)
  ELSE
     CALL Get(sF,TrixFile('OrthoF',Args,0))    ! the orthogonalized Fock matrix
  ENDIF
  !
  IF(MyID.EQ.0) THEN
     NSMat=sF%NSMat
     SELECT CASE(NSMat)
     CASE(1);NRow=  NBasF;NCol=  NBasF
     CASE(2);NRow=  NBasF;NCol=2*NBasF
     CASE(4);NRow=2*NBasF;NCol=2*NBasF
     CASE DEFAULT;CALL Halt(' RHeqs: sF%NSMat doesn''t have an expected value! ')
     END SELECT
  ENDIF
  !
#ifdef PARALLEL
  CALL BCast(NRow)
  CALL BCast(NCol)
  CALL BCast(NSMat)
#endif
  !
  CALL New(F,(/NRow,NCol/))
  CALL SetEq(F,sF)
  CALL Delete(sF) 

write(*,*) 'RHeqs: NAlph',NAlph
write(*,*) 'RHeqs: NBeta',NBeta

!----------------------------------
!
!
  CALL New(EigenV,NCol)
  CALL SetEq(EigenV,Zero)
  LWORK=MAX(1,3*NRow+10)
  CALL New(Work,LWork)
! 
  SELECT CASE(NSMat)
  CASE(1)
     ! We just have one matrix.
     CALL DSYEV('V','U',NBasF,F%D(1,1),NBasF,EigenV%D(1),Work%D(1),LWORK,Info)
     IF(Info/=SUCCEED)CALL Halt('DSYEV flaked in RHEqs. INFO='//TRIM(IntToChar(Info)))
     !
     HOMO=EigenV%D(NEl/2)
     LUMO=EigenV%D(NEl/2+1)
  CASE(2)
     ! We have a block diagonal matrix, we diag each blocks separately.
     ! Alpha Spin
     CALL DSYEV('V','U',NBasF,F%D(1,1),NBasF,EigenV%D(1),Work%D(1),LWORK,Info)
     IF(Info/=SUCCEED)CALL Halt('DSYEV flaked in RHEqs. INFO='//TRIM(IntToChar(Info)))
     ! Beta Spin
     CALL DSYEV('V','U',NBasF,F%D(1,NBasF+1),NBasF,EigenV%D(NBasF+1),Work%D(1),LWORK,Info)
     IF(Info/=SUCCEED)CALL Halt('DSYEV flaked in RHEqs. INFO='//TRIM(IntToChar(Info)))
     !
     IF(EigenV%D(NAlph)-EigenV%D(NAlph+1).LT.EigenV%D(NBasF+NBeta)-EigenV%D(NBasF+NBeta+1)) THEN
        HOMO=EigenV%D(NAlph  )
        LUMO=EigenV%D(NAlph+1)
     ELSE
        HOMO=EigenV%D(NBasF+NBeta  )
        LUMO=EigenV%D(NBasF+NBeta+1)
     ENDIF
  CASE(4)
     ! We just have one matrix.
     CALL DSYEV('V','U',2*NBasF,F%D(1,1),NRow,EigenV%D(1),Work%D(1),LWORK,Info)
     IF(Info/=SUCCEED)CALL Halt('DSYEV flaked in RHEqs. INFO='//TRIM(IntToChar(Info)))
     !
     HOMO=EigenV%D(NEl)
     LUMO=EigenV%D(NEl+1)
  CASE DEFAULT;CALL Halt(' RHeqs: NSMat doesn''t have an expected value! ')
  END SELECT
!
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
  CALL New(P,(/NRow,NCol/))
  CALL DBL_VECT_EQ_DBL_SCLR(NRow*NCol,P%D(1,1),Zero)
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
     !
     SELECT CASE(NSMat)
     CASE(1)
        ! We have one density matrix to build.
        CALL BuildP0(Nel/2,NBasF,F%D(1,1),P%D(1,1))! Restricted 
     CASE(2)
        ! We have two density matrices to build.
        CALL BuildP0(NAlph,NBasF,F%D(1,1)      ,P%D(1,1)      )! Unrestricted Alpha
        CALL BuildP0(NBeta,NBasF,F%D(1,NBasF+1),P%D(1,NBasF+1))! Unrestricted Beta
     CASE(4)
        ! We have one density matrix to build.
        CALL BuildP0(Nel,2*NBasF,F%D(1,1),P%D(1,1))! Generalized 
        ! Need to recompute NAlph and NBeta at this level.
     CASE DEFAULT;CALL Halt(' RHeqs: NSMat doesn''t have an expected value! ')
     END SELECT
     !
  CASE DEFAULT;CALL Halt('RHEqs: Doesn''t regonize this Smearing <'//TRIM(Smearing)//'>')
  END SELECT
!
  CALL Delete(EigenV)
!
  CALL SetEq(sX,P,nsmat_o=nsmat)          !  sX=P
  CALL New(sP,nsmat_o=nsmat)              
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
  !
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
CONTAINS
  SUBROUTINE BuildP0(Nel,N,C,P)
    INTEGER      :: Nel,N
    REAL(DOUBLE) :: C(N,*),P(N,N)
    INTEGER      :: I,J,K
    REAL(DOUBLE) :: CJK
    DO K=1,Nel
       DO J=1,N
          CJK=C(J,K)
          DO I=1,N
             P(I,J)=P(I,J)+C(I,K)*CJK  ! P_{ij}=Sum^{N_El}_k MO_ik MO_kj
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE BuildP0
END PROGRAM  RHEqs


