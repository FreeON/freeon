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
   TYPE(DBCSR)                     :: P,Tmp1,Tmp2
#else
   TYPE(BCSR)                      :: P,Tmp1,Tmp2
#endif
   REAL(DOUBLE)                    :: E_el_tot,E_nuc_tot,E_es_tot,KinE,ExchE,Exc, &
                                      Gap,Etot,DMax,Virial,DIISErr
   CHARACTER(LEN=5*DEFAULT_CHR_LEN):: SCFMessage
   CHARACTER(LEN=9),PARAMETER      :: Prog='SCFStatus'  
!---------------------------------------------------------------------------------------
!  Macro the start up
   CALL StartUp(Args,Prog)
!  Allocate some matrices
   CALL New(P)
   CALL New(Tmp1)
   CALL New(Tmp2)
!---------------------------------------------
!  Get the density matrix 
   IF(SCFActn=='BasisSetSwitch')THEN 
!     if switching the density matrix, 
!     this is all that is available ...    
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
!  E_el_tot=<Vee+Vne>=Tr{P.(Vee+Vne)}
   CALL Get(Tmp1,TrixFile('J',Args,0))
#ifdef PARALLEL
   CALL Multiply(P,Tmp1,Tmp2)
   E_el_tot=Trace(Tmp2)     
#else
   E_el_tot=Trace(P,Tmp1)    
#endif
!  ExchE=<Kx>=Tr{P.K}
   IF(HasHF(ModelChem))THEN
      CALL Get(Tmp1,TrixFile('K',Args,0))
#ifdef PARALLEL
      CALL Multiply(P,Tmp1,Tmp2)
      ExchE=ExactXScale(ModelChem)*Trace(Tmp2)     
#else
      ExchE=ExactXScale(ModelChem)*Trace(P,Tmp1)    
#endif
   ELSE
      ExchE=Zero
   ENDIF
!  Get E_nuc_tot =<Vnn+Vne> 
   CALL Get(E_nuc_tot,'enn+ene',Tag_O=SCFCycl)
!  Get the exchange correlation energy
   IF(HasDFT(ModelChem))THEN      
      CALL Get(Exc,'Exc',Tag_O=SCFCycl)
   ELSE
      Exc=Zero
   ENDIF
   E_es_tot = E_el_tot+E_nuc_tot
   Etot=KinE+E_es_tot+Exc+ExchE
   CALL Put(Etot,'Etot',StatsToChar(Current))
!  The Virial
   Virial=E_es_tot/KinE
!--------------------------------------------------------
!  Find the largest block of the delta density matrix
!  Allows for checking between extrapolated or projected DMs
   DMax=1.D10
   IF(SCFActn/='BasisSetSwitch')THEN
      CALL Get(Tmp1,TrixFile('D',Args,0))
      CALL Get(Tmp2,TrixFile('D',Args,1))
      CALL Multiply(Tmp1,-One)
      CALL Add(Tmp1,Tmp2,P)
      DMax=Max(P)
   ENDIF
   CALL Put(DMax,'DMax',StatsToChar(Current))
!  IO for the delta density matrix
   IF(SCFActn=='InkFok')THEN
      CALL Put(P,TrixFile('DeltaD',Args,1))
      CALL PChkSum(P,'DeltaP['//TRIM(NxtCycl)//']',Prog)
      CALL PPrint( P,'DeltaP['//TRIM(NxtCycl)//']')
      CALL Plot(   P,'DeltaP_'//TRIM(NxtCycl))
   ENDIF
   CALL Delete(P)
   CALL Delete(Tmp1)
   CALL Delete(Tmp2)
!-----------------------------------------------------------
!  Get DIIS Err and HOMO-LUMO Gap
   IF(Current(1)>=1)THEN
      CALL Get(DIISErr,'diiserr',StatsToChar(Current))
   ELSE
      DIISErr=Zero
   ENDIF
   CALL Get(Gap,'HomoLumoGap')
!-----------------------------------------------------------  
!  PRINT STATISTICS
!   
   SCFMessage=""
   IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
!     Fancy output
      SCFMessage=RTRN//'= = = = = SCFCycle #'//TRIM(SCFCycl)                          &
                                 //', Basis #'//TRIM(CurBase)               &
                              //', Geometry #'//TRIM(CurGeom)               &
               //' = = = = ='//RTRN//RTRN                    
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
   ELSEIF(PrintFlags%Key>=DEBUG_NONE)THEN
      SCFMessage=ProcessName(Prog,'['//TRIM(SCFCycl)//','  &
                                     //TRIM(CurBase)//','  &
                                     //TRIM(CurGeom)//']') &
                //'<SCF> = '//TRIM(FltToMedmChar(ETot))     &
                //', dD = '//TRIM(DblToShrtChar(DMax))
!     Add in DIIS error
      IF(DIISErr/=Zero)                                                     &
         SCFMessage=TRIM(SCFMessage)                                        &
               //', DIIS = '//TRIM(DblToShrtChar(DIISErr))
   ENDIF
#ifdef PARALLEL
  IF(MyId==ROOT)THEN
#endif
     WRITE(*,* )TRIM(SCFMessage)
     IF(PrintFlags%Key>=DEBUG_MEDIUM)THEN
        CALL OpenASCII(OutFile,Out)
        WRITE(Out,* )TRIM(SCFMessage)
        CLOSE(Out)
     ENDIF
#ifdef PARALLEL
  ENDIF
#endif
  CALL ShutDown(Prog)
END PROGRAM SCFStatus


