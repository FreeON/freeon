!    FAST O(N lg N) COMPUTATION OF THE COULOMB MATRIX
!    Authors:  Matt Challacombe and CJ Tymczak
!==============================================================================
PROGRAM QCTC
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE Globals
  USE AtomPairs
  USE BraBloks
  USE QCTCThresholds
  USE PoleTree
#ifdef PERIODIC
  USE PBCFarField
#endif
  USE JGen
  USE NuklarE
  TYPE(BCSR)                 :: J,T1,T2
  REAL(DOUBLE)               :: E_Nuc_Tot
  TYPE(TIME)                 :: TimeMakeJ
  CHARACTER(LEN=4),PARAMETER :: Prog='QCTC'
#ifdef MMech
  TYPE(CRDS)                 :: GM_MM
  REAL(DOUBLE)               :: MM_COUL,E_C_EXCL,CONVF  
  INTEGER                    :: UOUT !!!!
#endif
!-------------------------------------------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog)
! Get basis set and geometry
#ifdef MMech
  IF(HasQM()) THEN
#endif
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
! Allocations 
  CALL NewBraBlok(BS)
#ifdef MMech
  ENDIF
  IF(HasMM()) THEN
    CALL Get(GM_MM,Tag_O='GM_MM'//CurGeom)
  ENDIF
#endif
! Get multipoles and density
  IF(SCFActn=='InkFok')THEN
     CALL Get(Rho,'DeltaRho',Args,0)
     CALL Get(RhoPoles,'Delta'//TRIM(SCFCycl))
  ELSE  
     CALL Get(Rho,'Rho',Args,0)
     CALL Get(RhoPoles,SCFCycl)
  ENDIF
! Set thresholds local to QCTC (for PAC and MAC)
  CALL SetLocalThresholds(Thresholds%TwoE)
! Initialize the auxiliary density arrays
  CALL InitRhoAux
! Setup global arrays for computation of multipole tensors
  CALL MultipoleSetUp(FFEll2)
! Build the global PoleTree representation of the total density
  CALL RhoToPoleTree
#ifdef PERIODIC
! Calculate the Number of Cells
  CALL SetCellNumber(GM)
  CALL PPrint(CS_OUT,'outer sum',Prog)
! Set the electrostatic background 
  CALL PBCFarFieldSetUp(PoleRoot)
#endif
! Delete the auxiliary density arrays
  CALL DeleteRhoAux
! Delete the Density
  CALL Delete(Rho)
#ifdef MMech
  IF(HasQM()) THEN
#endif
     ! Allocate J
     CALL New(J)
     ! Compute the Coulomb matrix J in O(N Lg N)
     CALL Elapsed_Time(TimeMakeJ,'Init')
     CALL MakeJ(J)
     CALL Elapsed_TIME(TimeMakeJ,'Accum')
     IF(SCFActn=='InkFok')THEN
        !    Add in correction if incremental J build
        CALL New(T1)
        CALL New(T2)
        CALL Get(T1,TrixFile('J',Args,-1))
        CALL Add(T1,J,T2)
        CALL Filter(T1,T2)
        CALL Delete(T2)
     ELSE
        CALL Filter(T1,J)
     ENDIF
     ! Put J to disk
     CALL Put(T1,TrixFile('J',Args,0))
     ! Compute the nuclear-total electrostatic energy in O(N Lg N)
     IF(SCFActn=='InkFok')THEN
        CALL Get(E_Nuc_Tot,'E_NuclearTotal',Tag_O=PrvCycl)
        E_Nuc_Tot=E_Nuc_Tot+NukE(GM)
     ELSE     
        E_Nuc_Tot=NukE(GM)
     ENDIF
     CALL Put(E_Nuc_Tot,'E_NuclearTotal',Tag_O=SCFCycl)
#ifdef MMech
  ENDIF !!!!  QM calculations
  IF(HasMM()) THEN
!
     MM_COUL = NukE(GM_MM)
     CALL Put(MM_COUL,'MM_COUL',Tag_O=CurGeom)
!
     IF(HasQM()) THEN
!
       CALL Get(E_C_EXCL,'E_C_EXCL',Tag_O=CurGeom)
!
     CALL OpenASCII(OutFile,UOut)
!
     CONVF=1000.D0*JtoHartree/C_Avogadro
!
! Print energies 
!
     write(uout,*) 'Energies in KJ/mol'        
     write(uout,*) 'E_MM_Coulomb    TOTAL= ',MM_COUL/CONVF
     write(uout,*) 'E_MM_Coulomb EXCLUDED= ',E_C_EXCL/CONVF
     write(uout,*) 'E_MM_Coulomb         = ',(MM_COUL-E_C_EXCL)/CONVF
!
     write(uout,*) 'Energies in atomic unit'        
     write(uout,*) 'E_MM_Coulomb    TOTAL= ',MM_COUL
     write(uout,*) 'E_MM_Coulomb EXCLUDED= ',E_C_EXCL
     write(uout,*) 'E_MM_Coulomb         = ',MM_COUL-E_C_EXCL
!
     write(*,*) 'E_MM_Coulomb    TOTAL= ',MM_COUL
!
     CLOSE(UNIT=UOut,STATUS='KEEP')
!
     ELSE
       write(*,*) 'E_MM_Coulomb    TOTAL= ',MM_COUL
     ENDIF
!
  ENDIF
#endif


!---------------------------------------------------------------
! Printing
#ifdef MMech
  IF(HasQM()) THEN
#endif
    CALL PChkSum(T1,'J['//TRIM(SCFCycl)//']',Prog)
    CALL PPrint( T1,'J['//TRIM(SCFCycl)//']')
    CALL Plot(   T1,'J['//TRIM(SCFCycl)//']')
#ifdef MMech
  ENDIF
#endif
!  WRITE(*,*) 'NukE[',TRIM(SCFCycl),'] = ',E_Nuc_Tot
#ifdef PERIODIC
! Print Periodic Info
  CALL Print_Periodic()
#endif
! Tidy up
#ifdef MMech
  IF(HasQM()) THEN
#endif
    CALL Delete(J)
    CALL Delete(T1)
    CALL Delete(BS)
    CALL Delete(GM)
#ifdef MMech
  ENDIF
  IF(HasMM()) THEN
    CALL Delete(GM_MM)
  ENDIF
#endif
  CALL Delete(Args)
  CALL Delete(RhoPoles)
! didn't count flops, any accumulation is residual
! from matrix routines
  PerfMon%FLOP=Zero 
! Shutdown 
  CALL ShutDown(Prog)
END PROGRAM QCTC
