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
  CALL StartUp(Args,Prog)
  ISCF=Args%I%I(1)
  Cycl=IntToChar(ISCF)
!
  CALL New(F)
  CALL New(Tmp1)
  IF(Args%C%C(2)=='Core')THEN
     CALL Get(V,TrixFile('V',Args))                        ! V=V_{Nuc-El}
     CALL Get(T,TrixFile('T',Args))                        ! T=T_{Kinetic}
     CALL Add(V,T,F)                                       ! F=F_CoreH=V+T
     CALL Delete(V)
     CALL Delete(T)
  ELSEIF(Args%C%C(2)=='InkFok')THEN
!     CALL Get(J,TRIM(SCFName)//'_Cyc'//TRIM(Cycl)//'.DeltaJ') ! J=J_Coulomb[/\P]
!     CALL Get(K,TRIM(SCFName)//'_Cyc'//TRIM(Cycl)//'.DeltaK') ! K=K_Exchange[/\P]
     CALL Add(J,K,F)                                          ! F=/\F[/\P]
     CALL Delete(J)
     CALL Delete(K)
  ELSE
     CALL Get(J,TrixFile('J',Args,0))                      ! J=J_Coulomb[Rho_Total(Nuc+El)]
     CALL Get(T,TrixFile('T',Args))                        ! T=T_{Kinetic}
     CALL Add(J,T,F)                                       ! F=J+T    
     CALL Delete(T)
     CALL Delete(J)
     CALL New(K)
     IF(HasHF(ModelChem).AND.HasDFT(ModelChem))THEN
        CALL Get(K,TrixFile('K',Args,0))                      ! K=K_hf
        KScale=ExactXScale(ModelChem)
        CALL Multiply(K,KScale)                               ! K=KScale*K_hf
        CALL Add(F,K,Tmp1)                                    ! F=J+T+K
        CALL Get(K,TrixFile('Kxc',Args,0))                    ! K=K_xc
        CALL Add(K,Tmp1,F)                                    ! F=J+T+KShift*K_hf+K_xc
     ELSEIF(HasHF(ModelChem))THEN
        CALL Get(K,TrixFile('K',Args,0))                      ! K=K_hf
        CALL Add(F,K,Tmp1)                                    ! F=J+T+K_hf
        CALL SetEq(F,Tmp1)
     ELSEIF(HasDFT(ModelChem))THEN
        CALL Get(K,TrixFile('Kxc',Args,0))                    ! K=K_{xc}
        CALL Add(K,F,Tmp1)                                    ! F=J+T+K_xc
        CALL SetEq(F,Tmp1)
     ENDIF
     CALL Delete(K)
  ENDIF
  IF(Args%C%C(2)=='InkFok')THEN
!     CALL Put(F,TRIM(SCFName)//'_Cyc'//TRIM(Cycl)//'.DeltaF', &
!              BlksName_O='nfi'//TRIM(Cycl),                   &
!              Non0Name_O='nfm'//TRIM(Cycl))
!     CALL PChkSum(F,'DeltaF',Prog)
!     CALL PPrint( F,'DeltaF')
!     CALL Plot(   F,'DeltaF')
  ELSE
     CALL Put(F,TrixFile('F',Args,0),       &
              BlksName_O='nfi'//TRIM(Cycl), &
              Non0Name_O='nfm'//TRIM(Cycl))
     CALL PChkSum(F,'F['//TRIM(Cycl)//']',Prog)
     CALL PPrint( F,'F['//TRIM(Cycl)//']')
     CALL Plot(   F,'F['//TRIM(Cycl)//']')
  ENDIF
!----------------------------------------------
!
  XFile=TrixFile('X',Args)
  INQUIRE(FILE=XFile,EXIST=Present)
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
!--------------------------------------------------------------------
! 
!
  IF(Args%C%C(2)=='InkFok')THEN
     CALL Put(Tmp1,TrixFile('OrthoDeltaF',Args,0)) 
     CALL PChkSum(Tmp1,'OrthoDeltaF',Prog)
     CALL PPrint(Tmp1,'OrthoDeltaF')
     CALL Plot(Tmp1,'OrthoDeltaF')
  ELSE
     CALL Put(Tmp1,TrixFile('OrthoF',Args,0)) 
     CALL PChkSum(Tmp1,'OrthoF['//TRIM(SCFCycl)//']',Prog)
     CALL PPrint(Tmp1,'OrthoF['//TRIM(SCFCycl)//']')
     CALL Plot(Tmp1,'OrthoF_'//TRIM(SCFCycl))
  ENDIF
!--------------------------------
! Tidy up
!
  CALL Delete(F)
  CALL Delete(Tmp1)
  CALL ShutDown(Prog)
END PROGRAM FockNGrueven
