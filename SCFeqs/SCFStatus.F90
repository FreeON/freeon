PROGRAM SCFStatus
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalCharacters
   USE ProcessControl
   USE InOut
   USE PrettyPrint
   USE MemMan
   USE Parse
   USE Macros
   USE LinAlg
   USE Functionals
#ifdef PARALLEL
   USE MondoMPI
#endif
   IMPLICIT NONE
   TYPE(ARGMT)                     :: Args
#ifdef PARALLEL
   TYPE(DBCSR)                     :: P,Tmp1,Tmp2,Tmp3
#else
   TYPE(BCSR)                      :: P,Tmp1,Tmp2,Tmp3
#endif
   REAL(DOUBLE)                    :: E_el_tot,E_nuc_tot,E_es_tot,E_ECPs,KinE,ExchE,Exc, &
                                      Gap,Etot,DMax,Virial,DIISErr
#ifdef MMech
   REAL(DOUBLE)                    :: EBOND,EANGLE,ETorsion,ELJ,EOutOfPlane,MM_COUL,MM_ENERGY
   REAL(DOUBLE)                    :: E_C_EXCL,E_LJ_EXCL
#endif
   LOGICAL                         :: HasECPs
   CHARACTER(LEN=5*DEFAULT_CHR_LEN):: SCFMessage
   CHARACTER(LEN=9),PARAMETER      :: Prog='SCFStatus'  
   CHARACTER(LEN=2)                :: CurClone
!---------------------------------------------------------------------------------------
!  Macro the start up
   CALL StartUp(Args,Prog,Serial_O=.FALSE.)
!  Allocate some matrices
   CALL New(P)
   CALL New(Tmp1)
   CALL New(Tmp2)
!---------------------------------------------
!  Get the density matrix
   IF(SCFActn=='BasisSetSwitch' .OR. SCFActn=="RestartBasisSwitch") THEN
      ! If switching the density matrix or using a previous one from 
      ! restart use i+1 density matrix--its all that is available
      CALL Get(P,TrixFile('D',Args,1))
   ELSE
      CALL Get(P,TrixFile('D',Args,0))
   ENDIF
!---------------------------------------------
!  COMPUTE SOME EXPECTATION VALUES
!
!  KinE=<T>=Tr{P.T}
   CALL Get(Tmp1,TrixFile('T',Args))
#ifdef PARALLEL
   CALL Multiply(P,Tmp1,Tmp2)
   KinE=Two*Trace(Tmp2)     
#else
   KinE=Two*Trace(P,Tmp1)    
#endif 
   CALL Get(HasECPs,'hasecps',Tag_O=CurBase)
   IF(HasECPs)THEN
      ! Get the pseudopotential matrix U 
      CALL Get(Tmp1,TrixFile('U',Args))                        ! Tmp1=U_{ECP}
#ifdef PARALLEL
      CALL Multiply(P,Tmp1,Tmp2)
      E_ECPs=Two*Trace(Tmp2)     
#else
      E_ECPs=Two*Trace(P,Tmp1)    
#endif
   ELSE
      E_ECPs=Zero
   ENDIF
   ! E_el_tot=<Vee+Vne>=Tr{P.(Vee+Vne)}
   CALL Get(Tmp1,TrixFile('J',Args,0))
#ifdef PARALLEL
   CALL Multiply(P,Tmp1,Tmp2)
   E_el_tot=Trace(Tmp2)     
#else
   E_el_tot=Trace(P,Tmp1)    
#endif
   ! Total electrostatic energy icluding ECPs
   E_el_tot=E_el_tot+E_ECPs
   CALL Put(E_el_tot,'E_ElectronicTotal')
   ExchE=Zero
   Exc=Zero
   IF(SCFActn/="GuessEqCore")THEN
      ! ExchE=<Kx>=Tr{P.K}
      IF(HasHF(ModelChem))THEN
         CALL Get(Tmp1,TrixFile('K',Args,0))
#ifdef PARALLEL
         CALL Multiply(P,Tmp1,Tmp2)
         ExchE=ExactXScale(ModelChem)*Trace(Tmp2)     
#else
         ExchE=ExactXScale(ModelChem)*Trace(P,Tmp1)    
#endif
      ENDIF     
      !  Get the exchange correlation energy
      IF(HasDFT(ModelChem)) &
         CALL Get(Exc,'Exc',StatsToChar(Current))
   ENDIF
!  Get E_nuc_tot =<Vnn+Vne> 
   CALL Get(E_Nuc_Tot,'E_NuclearTotal',StatsToChar(Current))
   ! Total electrostaic energy
   E_es_tot=E_el_tot+E_nuc_tot
   ! Total SCF energy
   Etot=KinE+E_es_tot+Exc+ExchE
!   WRITE(*,*)' KinE = ',KinE
!   WRITE(*,*)' Elec = ',E_es_tot
!   WRITE(*,*)' Exc  = ',Exc
!   WRITE(*,*)' Exch = ',ExchE
!   WRITE(*,*)' ETot = ',ETot
   CALL Put(Etot,'Etot')
   CALL Put(Etot,'Etot',StatsToChar(Current))
!  The Virial
   Virial=E_es_tot/KinE
!--------------------------------------------------------
!  Find the largest block of the delta density matrix
!  Allows for checking between extrapolated or projected DMs
   IF(SCFActn=='BasisSetSwitch'.OR. SCFActn=="RestartBasisSwitch" .OR. SCFActn=='NumForceEvaluation')THEN
      DMax=Max(P)
   ELSE
      CALL Get(Tmp1,TrixFile('D',Args,0))
      CALL Get(Tmp2,TrixFile('D',Args,1))
      CALL Multiply(Tmp1,-One)
      CALL Add(Tmp1,Tmp2,P)
      DMax=Max(P)
   ENDIF
   CALL Put(DMax,'DMax')
!  IO for the delta density matrix
   IF(SCFActn=='InkFok')THEN
      CALL Put(P,TrixFile('DeltaD',Args,1))
      CALL PChkSum(P,'DeltaP['//TRIM(NxtCycl)//']',Prog)
      CALL PPrint( P,'DeltaP['//TRIM(NxtCycl)//']')
      CALL Plot(   P,'DeltaP_'//TRIM(NxtCycl))
   ENDIF

!-----------------------------------------------------------
!  Get DIIS Err and HOMO-LUMO Gap
   IF(Current(1)>=1.AND.SCFActn/='NumForceEvaluation') THEN  
      CALL Get(DIISErr,'diiserr')
   ELSE
      DIISErr=Zero
   ENDIF
   CALL Get(Gap,'HomoLumoGap')
!-----------------------------------------------------------  
!  PRINT STATISTICS
!  
   IF(NClones>1)THEN
      CurClone=IntToChar(MyClone)
   ELSE
      CurClone=""
   ENDIF
   SCFMessage=""
   IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
!     Fancy output
      SCFMessage=RTRN//'= = = = = SCFCycle #'//TRIM(SCFCycl)                &
                                 //', Basis #'//TRIM(CurBase)               &
                              //', Geometry #'//TRIM(CurGeom)              
      IF(NClones>1)THEN
         SCFMessage=TRIM(SCFMessage)//', Clone #'//TRIM(CurClone)              
      ENDIF
      SCFMessage=TRIM(SCFMessage)//' = = = = ='//RTRN//RTRN                    
!     Add in gap if RH
      IF(Gap/=Zero)                                                         &
         SCFMessage=TRIM(SCFMessage)                                        &
               //'       Gap     = '//TRIM(DblToShrtChar(-Gap))//RTRN         
!     Add in DIIS error 
      IF(DIISErr/=Zero)                                                     &
         SCFMessage=TRIM(SCFMessage)                                        &
               //'       DIISErr = '//TRIM(DblToShrtChar(DIISErr))//RTRN
      SCFMessage=TRIM(SCFMessage)                                           &
               //'       MaxDelD = '//TRIM(DblToShrtChar(DMax))//RTRN       &
               //'       <T>     = '//TRIM(DblToMedmChar(KinE))//RTRN       &
               //'       <V>     = '//TRIM(DblToMedmChar( E_es_tot))//RTRN         
!     Add in HF energy
      IF(ExchE/=Zero)                                                       &
         SCFMessage=TRIM(SCFMessage)                                        &
               //'       <HF>    = '//TRIM(DblToMedmChar(ExchE))//RTRN         
!     Add in DFT energy
      IF(Exc/=Zero)                                                         &
         SCFMessage=TRIM(SCFMessage)                                        &
               //'       <DFT>   = '//TRIM(DblToMedmChar(Exc))//RTRN         
!     Last but not least, total SCF energy
      SCFMessage=TRIM(SCFMessage)                                           &
               //'       <SCF>   = '//TRIM(FltToMedmChar(ETot))//RTRN
!     Add in MM energies
#ifdef MMech
      IF(HasMM()) THEN
         SCFMessage=TRIM(SCFMessage)                                        &
               //'     <EBOND>   = '//TRIM(DblToMedmChar(EBOND))//RTRN         
         SCFMessage=TRIM(SCFMessage)                                        &
               //'    <EANGLE>   = '//TRIM(DblToMedmChar(EANGLE))//RTRN         
         SCFMessage=TRIM(SCFMessage)                                        &
               //' <ETorsion>   = '//TRIM(DblToMedmChar(ETorsion))//RTRN         
         SCFMessage=TRIM(SCFMessage)                                        &
               //' <EOutOfPlane>   = '//TRIM(DblToMedmChar(EOutOfPlane))//RTRN         
         SCFMessage=TRIM(SCFMessage)                                        &
               //'<E_LENNARD_JONES>= '//TRIM(DblToMedmChar(ELJ))//RTRN         
         SCFMessage=TRIM(SCFMessage)                                        &
               //'   <MM_COUL>   = '//TRIM(DblToMedmChar(MM_COUL))//RTRN         
         SCFMessage=TRIM(SCFMessage)                                        &
               //' <MM_ENERGY>   = '//TRIM(FltToMedmChar(MM_ENERGY))//RTRN         
         SCFMessage=TRIM(SCFMessage)                                        &
               //'<TOTAL ENERGY> = '//TRIM(FltToMedmChar(Etot+MM_ENERGY))//RTRN         
      ENDIF
#endif
   ELSEIF(PrintFlags%Key>=DEBUG_NONE)THEN
      IF(NClones>1)THEN
         SCFMessage=ProcessName(Prog,'['//TRIM(SCFCycl)//','  &
                                        //TRIM(CurBase)//','  &
                                        //TRIM(CurGeom)//','//TRIM(CurClone)//']')
      ELSE
         SCFMessage=ProcessName(Prog,'['//TRIM(SCFCycl)//','  &
                                        //TRIM(CurBase)//','  &
                                        //TRIM(CurGeom)//']')
      ENDIF
      IF(SCFActn=='BasisSetSwitch')THEN
         SCFMessage=TRIM(SCFMessage)//' Basis set switch ... '       &
                                    //' MxD = '//TRIM(DblToShrtChar(DMax)) 

!      ELSEIF(SCFActn=='Restart')THEN
!         SCFMessage=TRIM(SCFMessage)//' Restarting ... '       &
!                                    //' MxD = '//TRIM(DblToShrtChar(DMax)) 
      ELSE
#ifdef MMech
IF(HasMM()) THEN
         SCFMessage=TRIM(SCFMessage)//' <SCF> = '//TRIM(FltToMedmChar(ETot)) &
                                    //' <MM_ENERGY> = '//TRIM(FltToMedmChar(MM_ENERGY)) &
                                    //' <TOTAL ENERGY> = '//TRIM(FltToMedmChar(Etot+MM_ENERGY)) &
                                    //', dD = '//TRIM(DblToShrtChar(DMax))
ELSE
#endif
         SCFMessage=TRIM(SCFMessage)//' <SCF> = '//TRIM(FltToMedmChar(ETot)) &
                                    //', dD = '//TRIM(DblToShrtChar(DMax))
#ifdef MMech
ENDIF
#endif
      ENDIF
!     Add in DIIS error
      IF(DIISErr/=Zero)                                                     &
         SCFMessage=TRIM(SCFMessage)                                        &
               //', DIIS = '//TRIM(DblToShrtChar(DIISErr))
   ENDIF
#ifdef PARALLEL
  IF(MyId==ROOT)THEN
#endif
     IF(SCFActn/='Silent')THEN
        CALL OpenASCII(OutFile,Out)
        WRITE(*,* )TRIM(SCFMessage)
        WRITE(Out,* )TRIM(SCFMessage)
        CLOSE(Out)
     ENDIF
#ifdef PARALLEL
  ENDIF
#endif
  CALL Delete(P)
  CALL Delete(Tmp1)
  CALL Delete(Tmp2)
  CALL ShutDown(Prog)
END PROGRAM SCFStatus


