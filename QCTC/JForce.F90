!    FAST O(N Lg N) COMPUTATION OF GRADIENTS OF THE COULOMB ENERGY 
!    WRT TO NUCLEAR COORDINATES
!    Author: Matt Challacombe
!==============================================================================
PROGRAM JForce
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
  USE PoleTree
#ifdef PERIODIC
  USE PBCFarField
#endif
  USE BlokTrPdJ
  USE BlokTrWdS
  USE NuklarE
#ifdef PARALLEL
  USE MondoMPI
  TYPE(DBCSR)                  :: P
#else
  TYPE(BCSR)                   :: P
#endif
  TYPE(AtomPair)               :: Pair
  TYPE(DBL_VECT)               :: Frc,JFrc
  INTEGER                      :: AtA,AtB,A1,A2,MA,NB,MN1,JP,Q
  REAL(DOUBLE)                 :: JFrcChk
  CHARACTER(LEN=6),PARAMETER   :: Prog='JForce'
#ifdef PERIODIC 
  INTEGER                      :: NC
  REAL(DOUBLE),DIMENSION(3)    :: B,F_nlm,nlm
  REAL(DOUBLE),DIMENSION(3,3)  :: LatFrc_J
#endif
#ifdef MMech
  INTEGER                      :: NatomsLoc
  TYPE(DBL_VECT)               :: MMJFrc
  TYPE(INT_VECT)               :: GlobalQMNum
  TYPE(INT_VECT)               :: AtmMark    
  TYPE(CRDS)                   :: GMLocMM,GMLoc
#else
  TYPE(CRDS)                   :: GMLoc
#endif
!-------------------------------------------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog)
#ifdef MMech
  IF(HasQM()) THEN
!    Get basis set and geometry
     CALL Get(BS,Tag_O=CurBase)
     CALL Get(GMLoc,Tag_O=CurGeom)
!    Allocations 
     CALL New(JFrc,3*GMLoc%Natms)
     JFrc%D(:)=Zero
     CALL NewBraBlok(BS,Gradients_O=.TRUE.)
     CALL Get(P,TrixFile('D',Args,1))
  ENDIF
  IF(HasMM()) THEN
    CALL Get(GMLocMM,Tag_O='GM_MM'//CurGeom)
    CALL New(AtmMark,GMLocMM%NAtms)
    CALL Get(AtmMark,'AtmMark')
    CALL New(MMJFrc,3*GMLocMM%NAtms)
    MMJFrc%D(:)=Zero
    CALL NewBraBlok(Gradients_O=.TRUE.)
  ENDIF
  CALL Get(Rho,'Rho',Args,1)
  IF(MMOnly()) THEN
     CALL Get(RhoPoles,CurGeom)
  ELSE
     CALL Get(RhoPoles,NxtCycl)
  ENDIF
#else
! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GMLoc,Tag_O=CurGeom)
! Allocations 
  CALL New(JFrc,3*GMLoc%Natms)
  JFrc%D(:)=Zero
  CALL NewBraBlok(BS,Gradients_O=.TRUE.)
  CALL Get(P,TrixFile('D',Args,1))
  CALL Get(Rho,'Rho',Args,1)
  CALL Get(RhoPoles,NxtCycl)
#endif   
! Set thresholds local to JForce (for PAC and MAC)
  CALL SetLocalThresholds(Thresholds%TwoE)
! Setup global arrays for computation of multipole tensors
  CALL InitRhoAux
! Setup global arrays for computation of multipole tensors
  CALL MultipoleSetUp(FFEll2)
! Build the global PoleTree representation of the total density
  CALL RhoToPoleTree
! Calculate the Number of Cells
! and Set the electrostatic background
#ifdef MMech 
#ifdef PERIODIC
  IF(HasMM()) THEN
     CALL SetCellNumber(GMLocMM)
     CALL PBCFarFieldSetUp(PoleRoot,GMLocMM)
  ELSE
     CALL SetCellNumber(GMLoc)
     CALL PBCFarFieldSetUp(PoleRoot,GMLoc)
  ENDIF
#endif
#else
#ifdef PERIODIC
  CALL SetCellNumber(GMLoc)
  CALL PBCFarFieldSetUp(PoleRoot,GMLoc)
  CALL PPrint(CS_OUT,'outer sum',Prog)
#endif
#endif
! Delete the auxiliary density arrays
  CALL DeleteRhoAux
! Delete the Density
  CALL Delete(Rho)
!--------------------------------------------------------------------------------
! Compute the Coulomb contribution to the force in O(N Lg N)
!--------------------------------------------------------------------------------
#ifdef MMech
  IF(HasQM()) THEN
#endif
  DO AtA=1,GMLoc%Natms
     MA=BSiz%I(AtA)
     A1=3*(AtA-1)+1
     A2=3*AtA
     JFrc%D(A1:A2)= dNukE(GMLoc,AtA)
     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
        AtB=P%ColPt%I(JP)
        IF(SetAtomPair(GMLoc,BS,AtA,AtB,Pair))THEN
           Q=P%BlkPt%I(JP)
           NB=BSiz%I(AtB)
           MN1=MA*NB-1
#ifdef PERIODIC
           B=Pair%B
           DO NC=1,CS_OUT%NCells
              Pair%B=B+CS_OUT%CellCarts%D(:,NC)
              Pair%AB2=(Pair%A(1)-Pair%B(1))**2 &
                      +(Pair%A(2)-Pair%B(2))**2 &
                      +(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair))THEN
                 F_nlm(1:3)    = TrPdJ(Pair,P%MTrix%D(Q:Q+MN1),GMLoc)
                 JFrc%D(A1:A2) = JFrc%D(A1:A2) + F_nlm(1:3)
              ENDIF
           ENDDO
#else
           JFrc%D(A1:A2)=JFrc%D(A1:A2)+TrPdJ(Pair,P%MTrix%D(Q:Q+MN1),GMLoc)
#endif
        ENDIF
     ENDDO
  ENDDO
! Closed shell...
  JFrc%D=Two*JFrc%D
#ifdef MMech
  ENDIF
#endif
!--------------------------------------------------------------------------------
! Print The JForce
!  CALL Print_Force(GMLoc,JFrc,'QM dJ/dR in au ')
!  JFrc%D(:)=JFrc%D(:)/KJPerMolPerAngstToHPerBohr
!  CALL Print_Force(GMLoc,JFrc,'QM dJ/dR in KJ/mol/A ')
!  JFrc%D(:)=JFrc%D(:)*KJPerMolPerAngstToHPerBohr
#ifdef MMech
  IF(HasMM()) THEN
     DO AtA=1,GMLocMM%NAtms
        IF(AtmMark%I(AtA)==0) THEN
           A1=3*(AtA-1)+1
           A2=3*AtA
           MMJFrc%D(A1:A2)=dNukE(GMLocMM,AtA)
        ENDIF
     ENDDO
!    'closed shell'
     MMJFrc%D=Two*MMJFrc%D
!    checks, prints
     CALL Print_Force(GMLocMM,MMJFrc,'MM dJ/dR in au ')
     MMJFrc%D=MMJFrc%D/KJPerMolPerAngstToHPerBohr
     CALL Print_Force(GMLocMM,MMJFrc,'MM dJ/dR in KJ/mol/A ')
     MMJFrc%D=MMJFrc%D*KJPerMolPerAngstToHPerBohr
  ENDIF
#endif 
!--------------------------------------------------------------------------------
#ifdef MMech
! Sum in contribution to total force
  IF(HasMM()) THEN
    NatmsLoc=GMLocMM%Natms
  ELSE IF(HasQM()) THEN
    NatmsLoc=GMLoc%Natms
  ELSE
    CALL Halt('No mechanics defined in JForce') 
  ENDIF
  CALL New(Frc,3*NatmsLoc)
  CALL Get(Frc,'GradE',Tag_O=CurGeom)
  CALL New(GlobalQMNum,NatmsLoc)
  IF(HasQM().AND.HasMM()) THEN
     CALL Get(GlobalQMNum,'GlobalQMNum')
  ELSE
     DO I=1,GMLoc%Natms; GlobalQMNum%I(I)=I ; ENDDO
     IF(HasQM()) THEN
        Do Ata=1,GMLoc%Natms 
           I=GlobalQMNum%I(Ata)
           IF(I/=0) THEN !!! don't save gradients on link atoms
              I1=3*(I-1)+1
              I2=3*I
              A1=3*(AtA-1)+1
              A2=3*AtA
              Frc%D(I1:I2)=Frc%D(I1:I2)+JFrc%D(A1:A2)
           ENDIF
        ENDDO
     ENDIF
     CALL PChkSum(Frc,'Frc bef dJ/dR added',Proc_O=Prog)  
     IF(HasMM()) THEN
        Frc%D=Frc%D+MMJFrc%D
        CALL PChkSum(MMJFrc,'MMJFrc bef dJ/dR added',Proc_O=Prog)  
     ENDIF
     CALL PChkSum(Frc,'Frc after dJ/dR added',Proc_O=Prog)  
     CALL Put(Frc,'GradE',Tag_O=CurGeom)
     CALL Delete(Frc)
     CALL Delete(GlobalQMNum)
! Tidy up
!--------------------------------------------------------------------------------
     IF(HasQM()) THEN
        CALL Delete(BS)
        CALL Delete(GMLoc)
        CALL Delete(JFrc)
        CALL Delete(P)
     ENDIF
     IF(HasMM()) THEN
        CALL Delete(GMLocMM)
        CALL Delete(MMJFrc)
        CALL Delete(AtmMark)
     ENDIF
  ENDIF
#else
  CALL New(Frc,3*GMLoc%Natms)
  CALL Get(Frc,'GradE',Tag_O=CurGeom)
  Frc%D=Frc%D+JFrc%D
  CALL PChkSum(Frc,'Frc after dJ/dR added',Proc_O=Prog)  
  CALL Put(Frc,'GradE',Tag_O=CurGeom)
  CALL Delete(Frc)
  CALL Delete(BS)
  CALL Delete(GMLoc)
  CALL Delete(JFrc)
#endif
  CALL Delete(RhoPoles)
  CALL DeleteBraBlok(Gradients_O=.TRUE.)
! didn't count flops, any accumulation is residual
! from matrix routines
  PerfMon%FLOP=Zero 
  CALL ShutDown(Prog)
END PROGRAM JForce
