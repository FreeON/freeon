PROGRAM FockNGrueven
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
  USE Functionals
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(DBCSR)                    :: F,V,T,J,K,X,Tmp1
#else
  TYPE(BCSR)                     :: F,V,T,J,K,X,Tmp1 
#endif
  TYPE(ARGMT)                    :: Args
  REAL(DOUBLE)                   :: KScale
  INTEGER                        :: ISCF
  CHARACTER(LEN=2)               :: Cycl
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: XFile,Mssg
  CHARACTER(LEN=12),PARAMETER    :: Prog='FockNGrueven'
  LOGICAL                        :: Present,ExchangeShift
!------------------------------------------------------------------ 
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  ISCF=Args%I%I(1)
  Cycl=IntToChar(ISCF)

!
  CALL New(F)
  CALL New(Tmp1)
! Start with the Coulomb matrix
  CALL Get(J,TrixFile('J',Args,0))                      ! J=J_Coulomb[Rho_Total(Nuc+El)]
! Then add the kinetic energy to get a preliminary Fock matrix
  CALL Get(T,TrixFile('T',Args))                        ! T=T_{Kinetic}
  CALL Add(J,T,F)                                       ! F=J+T    
  CALL Delete(T)
  CALL Delete(J)
  CALL New(K)
  IF(HasHF(ModelChem).AND.HasDFT(ModelChem))THEN
!    Add in Hartree-Fock exact exchange 
     CALL Get(K,TrixFile('K',Args,0))                      ! K=K_hf
     KScale=ExactXScale(ModelChem)
     CALL Multiply(K,KScale)                               ! K=KScale*K_hf
     CALL Add(F,K,Tmp1)                                    ! F=J+T+K
!    Now add in exchange-correlation
     CALL Get(K,TrixFile('Kxc',Args,0))                    ! K=K_xc
     CALL Add(K,Tmp1,F)                                    ! F=J+T+KShift*K_hf+K_xc
  ELSEIF(HasHF(ModelChem))THEN
!    Add in Hartree-Fock exact exchange 
     CALL Get(K,TrixFile('K',Args,0))                      ! K=K_hf
     CALL Add(F,K,Tmp1)                                    ! F=J+T+K_hf
     CALL SetEq(F,Tmp1)
  ELSEIF(HasDFT(ModelChem))THEN
!    Add in exchange-correlation
     CALL Get(K,TrixFile('Kxc',Args,0))                    ! K=K_{xc}
     CALL Add(K,F,Tmp1)                                    ! F=J+T+K_xc
     CALL SetEq(F,Tmp1)
  ENDIF
  CALL Delete(K)
!
  CALL Put(F,TrixFile('F',Args,0))
  CALL PChkSum(F,'F['//TRIM(Cycl)//']',Prog)
  CALL PPrint( F,'F['//TRIM(Cycl)//']')
  CALL Plot(   F,'F['//TRIM(Cycl)//']')

! Now transform to an orthogonal representation
  XFile=TrixFile('X',Args)
#ifdef PARALLEL
  IF(MyId==ROOT)THEN
#endif
     INQUIRE(FILE=XFile,EXIST=Present)
#ifdef PARALLEL
  ENDIF
  CALL BCast(Present)
#endif
  IF(Present)THEN
     CALL Get(X,XFile)                ! X =S^{-1/2}
     CALL Multiply(X,F,Tmp1)          ! T1=S^{-1/2}.F
     CALL Multiply(Tmp1,X,F)          ! T1=S^{-1/2}.F.S^{-1/2}
  ELSE
     CALL Get(X,TrixFile('ZT',Args))  ! X =Z^t=L^{-t}           
     CALL Multiply(X,F,Tmp1)          ! T1=Z^t.F_AO
     CALL Get(X,TrixFile('Z',Args))   ! X=Z=L^{-1}
     CALL Multiply(Tmp1,X,F)          ! F=Z^t.F_AO.Z
  ENDIF
  CALL Filter(Tmp1,F)                 ! T1 =F_Orthog=Filter[Z^t.F_AO.Z]
!
  CALL Put(Tmp1,TrixFile('OrthoF',Args,0)) 
  CALL PChkSum(Tmp1,'OrthoF['//TRIM(SCFCycl)//']',Prog)
  CALL PPrint(Tmp1,'OrthoF['//TRIM(SCFCycl)//']')
  CALL Plot(Tmp1,'OrthoF_'//TRIM(SCFCycl))
! Tidy up
  CALL Delete(F)
  CALL Delete(Tmp1)
  CALL ShutDown(Prog)
END PROGRAM FockNGrueven
