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
  USE McMurchie
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
  TYPE(CRDS)                     :: GM
  TYPE(BSET)                     :: BS
  REAL(DOUBLE)                   :: KScale,E_NukeNuke
  INTEGER                        :: ISCF
  CHARACTER(LEN=2)               :: Cycl
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: XFile,DevFile,Mssg
  CHARACTER(LEN=12),PARAMETER    :: Prog='FockNGrueven'
  LOGICAL                        :: Present,ExchangeShift
!------------------------------------------------------------------ 
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  ISCF=Args%I%I(1)
  Cycl=IntToChar(ISCF)
  CALL New(F)
  CALL New(Tmp1)
! Start with the Coulomb matrix
  CALL Get(J,TrixFile('J',Args,0))                      ! J=J_Coulomb[Rho_Total(Nuc+El)]
! Then add the kinetic energy to get a preliminary Fock matrix
  CALL Get(T,TrixFile('T',Args))                        ! T=T_{Kinetic}
  CALL Add(J,T,F)                                       ! F=J+T    
  CALL Delete(T)
  CALL Delete(J)
  IF(SCFActn/='GuessEqCore')THEN
     CALL New(K)
     IF(HasHF(ModelChem).AND.HasDFT(ModelChem))THEN
        ! Add in Hartree-Fock exact exchange 
        CALL Get(K,TrixFile('K',Args,0))                      ! K=K_hf
        KScale=ExactXScale(ModelChem)
        CALL Multiply(K,KScale)                               ! K=KScale*K_hf
        CALL Add(F,K,Tmp1)                                    ! F=J+T+K
        ! Now add in exchange-correlation
        CALL Get(K,TrixFile('Kxc',Args,0))                    ! K=K_xc
        CALL Add(K,Tmp1,F)                                    ! F=J+T+KShift*K_hf+K_xc
     ELSEIF(HasHF(ModelChem))THEN
        !  Add in Hartree-Fock exact exchange 
        CALL Get(K,TrixFile('K',Args,0))                      ! K=K_hf
        CALL Add(F,K,Tmp1)                                    ! F=J+T+K_hf
        CALL SetEq(F,Tmp1)
     ELSEIF(HasDFT(ModelChem))THEN
        !  Add in exchange-correlation
        CALL Get(K,TrixFile('Kxc',Args,0))                    ! K=K_{xc}
        CALL Add(K,F,Tmp1)                                    ! F=J+T+K_xc
        CALL SetEq(F,Tmp1)
     ENDIF
     CALL Delete(K)
#ifdef !PARALLEL
  ELSE
     IF(PrintFlags%Int==DEBUG_INTEGRAL)THEN
        ! This is a full dump run, with output of all matrices
        ! and integrals into the development directory MMA/SCFDev     
        CALL GetEnv('MONDO_HOME',DevFile)
        DevFile=TRIM(DevFile)//'/SCFeqs/MMA/SCFDev/scf.in'
        WRITE(*,*)' DUMP TAKEN IN ',DevFile
        CALL OpenASCII(DevFile,Tmp,NewFile_O=.TRUE.)
        ! Get basis set and geometry
        CALL Get(BS,Tag_O="1")
        CALL Get(GM,Tag_O="1")
        CALL Get(E_NukeNuke,'E_NuclearTotal')
        WRITE(Tmp,*)' (* All matrices are in an orthogonal rep *)'
        WRITE(Tmp,*)' NBasF  = ',NBasF,';'
        WRITE(Tmp,*)' NEl    = ',NEl,';'
        WRITE(Tmp,*)' ENukeNuke = ',E_NukeNuke,';'
        PrintFlags%Fmt=DEBUG_MMASTYLE
        PrintFlags%Mat=DEBUG_MATRICES
        ! Get and pprint the ao-orthog transformation matrix
        XFile=TrixFile('X',Args)
        INQUIRE(FILE=XFile,EXIST=Present)
        IF(Present)THEN
           CALL Get(X,XFile)                ! X =S^{-1/2}
           CALL PPrint(X,'X',FileName_O=DevFile,Unit_O=Tmp)
        ELSE
           CALL Get(X,TrixFile('ZT',Args))  ! X =Z^t=L^{-t}           
           CALL PPrint(X,'ZT',FileName_O=DevFile,Unit_O=Tmp)
           CALL Get(X,TrixFile('Z',Args))   ! X =Z=L^{-1}           
           CALL PPrint(X,'Z',FileName_O=DevFile,Unit_O=Tmp)
        ENDIF
        ! Print the core Hamiltonian
        CALL PPrint(F,'h',FileName_O=DevFile,Unit_O=Tmp)
        ! Print the two-electron integrals
        CALL PrintIntegrals(BS%NASym,2*BS%NASym,4*BS%NASym, &
                         BS%LMNLen,BS%LMNLen**4,BS,GM,3,Tmp) 
        CLOSE(UNIT=Tmp,STATUS='KEEP')
        STOP ! Done for good and all
     ENDIF
#endif
  ENDIF
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
