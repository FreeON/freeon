PROGRAM ResponseStatus
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
   CHARACTER(LEN=14),PARAMETER     :: Prog='ResponseStatus'  
!---------------------------------------------------------------------------------------
!  Macro the start up
   CALL StartUp(Args,Prog,Serial_O=.FALSE.)
!  Allocate some matrices
   CALL New(P)
   CALL New(Tmp1)
   CALL New(Tmp2)
!---------------------------------------------
!  Get the density matrix 
   IF(SCFActn=='BasisSetSwitch'.OR.SCFActn=='Restart')THEN 
      ! If switching the density matrix or using a previous one from 
      ! restart use i+1 density matrix--its all that is available
      CALL Get(P,TrixFile('D',Args,1))
   ELSE
      CALL Get(P,TrixFile('D',Args,0))
   ENDIF
!---------------------------------------------
!  COMPUTE SOME EXPECTATION VALUES
!
   CALL Put(Etot,'Etot')
!--------------------------------------------------------
!  Find the largest block of the delta density matrix
!  Allows for checking between extrapolated or projected DMs
   CALL Get(Tmp1,TrixFile('D',Args,0))
   CALL Get(Tmp2,TrixFile('D',Args,1))
   CALL Multiply(Tmp1,-One)
   CALL Add(Tmp1,Tmp2,P)
   DMax=Max(P)

   CALL Put(DMax,'DMax')


   CALL Delete(P)
   CALL Delete(Tmp1)
   CALL Delete(Tmp2)
!-----------------------------------------------------------
!  Get DDIIS Err
   IF(Current(1)>=1.AND.SCFActn/='NumForceEvaluation')THEN
      CALL Get(DIISErr,'diiserr')
   ELSE
      DIISErr=Zero
   ENDIF
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
!     Add in DIIS error
      IF(DIISErr/=Zero)                                                     &
         SCFMessage=TRIM(SCFMessage)                                        &
               //', DIIS = '//TRIM(DblToShrtChar(DIISErr))
   ENDIF
   CALL OpenASCII(OutFile,Out)
   WRITE(*,* )TRIM(SCFMessage)
   WRITE(Out,* )TRIM(SCFMessage)
   CLOSE(Out)

   CALL ShutDown(Prog)
 END PROGRAM ResponseStatus


