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
  USE PBCFarField
  USE BlokTrPdJ
  USE BlokTrWdS
  USE NuklarE
  USE SetXYZ 
  IMPLICIT NONE

#ifdef PARALLEL
  TYPE(BCSR)                   :: P
  TYPE(DBL_VECT)               :: TotJFrc
  INTEGER                      :: IErr,TotFrcComp,MyAtomNum,TotAtomNum
  REAL(DOUBLE)                 :: JFrcBegTm,JFrcEndTm,JFrcTm
  TYPE(DBL_VECT)               :: TmJFrcArr
#else
  TYPE(BCSR)                   :: P
#endif
  TYPE(AtomPair)               :: Pair
  TYPE(DBL_VECT)               :: Frc,JFrc
  INTEGER                      :: AtA,AtB,A1,A2,MA,NB,MN1,JP,Q
  REAL(DOUBLE)                 :: JFrcChk
  CHARACTER(LEN=6),PARAMETER   :: Prog='JForce'
  INTEGER                      :: NC,I,J
  REAL(DOUBLE),DIMENSION(3)    :: B,nlm
  REAL(DOUBLE),DIMENSION(15)   :: F_nlm
  TYPE(DBL_RNK2)               :: LatFrc_J
#ifdef MMech
  INTEGER                      :: NatmsLoc,I,I1,I2
  TYPE(DBL_VECT)               :: MMJFrc
  TYPE(INT_VECT)               :: GlobalQMNum
  TYPE(INT_VECT)               :: AtmMark    
  TYPE(CRDS)                   :: GMLocMM,GMLoc
  TYPE(DBL_RNK2)               :: GradAux
#else
  TYPE(CRDS)                   :: GMLoc
#endif
  REAL(DOUBLE),EXTERNAL        :: MondoTimer
!-------------------------------------------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
!
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
  IF(MMOnly()) THEN
     CALL Get(Rho,'Rho',Args,Current(1))
     CALL Get(RhoPoles,CurGeom)
  ELSE
     CALL Get(Rho,'Rho',Args,1)
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
#ifdef PARALLEL
  CALL New(P,OnAll_O=.TRUE.)
#endif
  CALL Get(P,TrixFile('D',Args,1),BCast_O=.TRUE.)

#ifdef PARALLEL
  CALL GetDistrRho('Rho',Args,1)
#else
  CALL Get(Rho,'Rho',Args,1,Bcast_O=.TRUE.)  
#endif
  CALL Get(RhoPoles)
#endif   
! Set thresholds local to JForce (for PAC and MAC)
  CALL SetLocalThresholds(Thresholds%TwoE)
! Setup global arrays for computation of multipole tensors
  CALL InitRhoAux
! Setup global arrays for computation of multipole tensors
  CALL MultipoleSetUp()
! Build the global PoleTree representation of the total density
#ifdef PARALLEL
  CALL ParaRhoToPoleTree
#else
  CALL RhoToPoleTree
#endif
#ifdef MMech 
! Set the electrostatic background
  IF(HasMM()) THEN
     CALL PBCFarFieldSetUp(PoleRoot,GMLocMM)
  ELSE
     CALL PBCFarFieldSetUp(PoleRoot,GMLoc)
  ENDIF
#else
! Set the electrostatic background
  CALL PBCFarFieldSetUp(PoleRoot,GMLoc)
#endif
! Delete the auxiliary density arrays
  CALL DeleteRhoAux
! Delete the Density
  CALL Delete(Rho)
!--------------------------------------------------------------------------------
! Compute the Coulomb contribution to the force in O(N Lg N)
!--------------------------------------------------------------------------------
!
  CALL New(LatFrc_J,(/3,3/))
  LatFrc_J%D = Zero
!
#ifdef MMech
  IF(HasQM()) THEN
#endif
#ifdef PARALLEL
  MyAtomNum = End%I(MyId)-Beg%I(MyId)+1
  CALL MPI_AllReduce(MyAtomNum,TotAtomNum,1,MPI_INTEGER,MPI_SUM,MONDO_COMM,IErr)
  IF(TotAtomNum /= GMLoc%Natms) THEN
    WRITE(*,*) 'TotAtomNum = ',TotAtomNum
    WRITE(*,*) 'GMLoc%Natms = ',GMLoc%Natms
    STOP 'TotAtomNum not equal to GMLoc%Natms in JForce!'
  ENDIF
  JFrcBegTm = MondoTimer()
  DO AtA=Beg%I(MyId),End%I(MyId)
#else
  DO AtA=1,GMLoc%Natms
#endif
     MA=BSiz%I(AtA)
     A1=3*(AtA-1)+1
     A2=3*AtA
     F_nlm = dNukE(GMLoc,AtA)
     JFrc%D(A1:A2)= Two*F_nlm(1:3)
!    Store Inner Nuc Lattice Forces
     LatFrc_J%D(1:3,1) = LatFrc_J%D(1:3,1) + F_nlm(7:9)
     LatFrc_J%D(1:3,2) = LatFrc_J%D(1:3,2) + F_nlm(10:12)
     LatFrc_J%D(1:3,3) = LatFrc_J%D(1:3,3) + F_nlm(13:15)
!    Outer Nuc Lattice Forces
     nlm        = AtomToFrac(GMLoc,GMLoc%Carts%D(:,AtA))
     LatFrc_J%D = LatFrc_J%D + Two*LaticeForce(GMLoc,nlm,F_nlm(1:3))
!    Start AtB Loop
     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1 
        AtB=P%ColPt%I(JP)
        IF(SetAtomPair(GMLoc,BS,AtA,AtB,Pair)) THEN 
           Q=P%BlkPt%I(JP)
           NB=BSiz%I(AtB)
           MN1=MA*NB-1
           B=Pair%B
           DO NC=1,CS_OUT%NCells
              Pair%B=B+CS_OUT%CellCarts%D(:,NC)
              Pair%AB2=(Pair%A(1)-Pair%B(1))**2 &
                      +(Pair%A(2)-Pair%B(2))**2 &
                      +(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair))THEN
                 F_nlm = TrPdJ(Pair,P%MTrix%D(Q:Q+MN1),GMLoc)
                 IF(Pair%SameAtom) THEN
                    JFrc%D(A1:A2) = JFrc%D(A1:A2) +  Four*F_nlm(4:6)
                 ELSE
                    JFrc%D(A1:A2) = JFrc%D(A1:A2) + Eight*F_nlm(1:3)
                 ENDIF
!                Store Inner J Lattice Forces
                 LatFrc_J%D(1:3,1) = LatFrc_J%D(1:3,1) + Two*F_nlm(7:9)
                 LatFrc_J%D(1:3,2) = LatFrc_J%D(1:3,2) + Two*F_nlm(10:12)
                 LatFrc_J%D(1:3,3) = LatFrc_J%D(1:3,3) + Two*F_nlm(13:15)
!                Outer Lattice J Forces
                 nlm        = AtomToFrac(GMLoc,Pair%A) 
                 LatFrc_J%D = LatFrc_J%D+Four*LaticeForce(GMLoc,nlm,F_nlm(1:3))
                 nlm        = AtomToFrac(GMLoc,Pair%B) 
                 LatFrc_J%D = LatFrc_J%D+Four*LaticeForce(GMLoc,nlm,(F_nlm(4:6)-F_nlm(1:3)))
              ENDIF
           ENDDO
!
        ENDIF
     ENDDO
  ENDDO
! Dipole Correction
  IF(GMLoc%PBC%Dimen>0) THEN
     DO I=1,3
        LatFrc_J%D(I,I) = LatFrc_J%D(I,I)-E_DP/GMLoc%PBC%BoxShape%D(I,I)
     ENDDO
  ENDIF
   LatFrc_J%D = Zero
#ifdef PARALLEL
  JFrcEndTm = MondoTimer()
  JFrcTm = JFrcEndTm-JFrcBegTm
#endif
#ifdef MMech
  ENDIF
#endif
!
! Print The Forces and Lattice Forces
!
!!$  WRITE(*,*) 'JForce'
!!$  DO AtA=1,NAtoms
!!$     A1=3*(AtA-1)+1
!!$     A2=3*AtA
!!$     WRITE(*,*) JFrc%D(A1:A2)
!!$  ENDDO
!!$  WRITE(*,*) 'LatFrc_J'
!!$  DO I=1,3
!!$     WRITE(*,*) (LatFrc_J%D(I,J),J=1,3) 
!!$  ENDDO
!--------------------------------------------------------------------------------
! Print The JForce
!  CALL Print_Force(GMLoc,JFrc,'dJ/dR in au ')
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
     IF(PrintFlags%GeOp==DEBUG_GEOP) THEN
       CALL Print_Force(GMLocMM,MMJFrc,'MM dJ/dR in au ')
       MMJFrc%D=MMJFrc%D/KJPerMolPerAngstToHPerBohr
       CALL Print_Force(GMLocMM,MMJFrc,'MM dJ/dR in KJ/mol/A ')
       MMJFrc%D=MMJFrc%D*KJPerMolPerAngstToHPerBohr
     ENDIF
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
  CALL New(GradAux,(/3,NatmsLoc/))
  CALL Get(GradAux,'gradients',Tag_O=CurGeom)
  CALL CartRNK2ToCartRNK1(Frc%D,GradAux%D)
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
     CALL PChkSum(Frc ,'Frc bef dJ/dR added',Proc_O=Prog)  
     IF(HasMM()) THEN
        Frc%D=Frc%D+MMJFrc%D
        CALL PChkSum(MMJFrc,'MMJFrc bef dJ/dR added',Proc_O=Prog)  
     ENDIF
     CALL PChkSum(Frc,'Frc after dJ/dR added',Proc_O=Prog)  
     CALL CartRNK1ToCartRNK2(Frc%D,GradAux%D)
     CALL Put(GradAux,'gradients',Tag_O=CurGeom)
     CALL Delete(GradAux)
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
#ifdef PARALLEL
  TotFrcComp = 3*GMLoc%Natms
  CALL New(TotJFrc,TotFrcComp)
  CALL MPI_Reduce(JFrc%D(1),TotJFrc%D(1),TotFrcComp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MONDO_COMM,IErr)
  IF(MyID == 0) THEN
    JFrc%D(1:TotFrcComp) = TotJFrc%D(1:TotFrcComp)
  ENDIF
  CALL Delete(TotJFrc)
#endif
  CALL PChkSum(JFrc,'dJ/dR',Proc_O=Prog)  
! Sum in contribution to total force
  DO AtA=1,NAtoms
     A1=3*(AtA-1)+1
     A2=3*AtA
     GMLoc%Gradients%D(1:3,AtA) = GMLoc%Gradients%D(1:3,AtA)+JFrc%D(A1:A2)
  ENDDO
#ifdef NLATTFORCE
  WRITE(*,*) 'JForce: Not Putting Lattice Force to Disk'
#else
  GMLoc%PBC%LatFrc%D = GMLoc%PBC%LatFrc%D+LatFrc_J%D
#endif
  CALL Put(GMLoc,Tag_O=CurGeom)
!
  CALL Delete(BS)
  CALL Delete(GMLoc)
  CALL Delete(LatFrc_J)
  CALL Delete(JFrc)
#endif
  CALL Delete(RhoPoles)
  CALL DeleteBraBlok(Gradients_O=.TRUE.)
#ifdef PARALLEL
  CALL New(TmJFrcArr,NPrc)
  CALL MPI_Gather(JFrcTm,1,MPI_DOUBLE_PRECISION,TmJFrcArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  IF(MyID == ROOT) THEN
     ! Output needs a lot of work, and should not go to STDOUT
     ! Also, these statistics were written long ago by      
     ! Elapsed_TIME(T,Init_O,Proc_O) in PrettyPrint.  Why create
     ! another routine to do this????
     ! CALL PImbalance(TmJFrcArr,NPrc,Prog_O='JFrc')
  ENDIF
  CALL Delete(TmJFrcArr)
#endif
! didn't count flops, any accumulation is residual from matrix routines
  PerfMon%FLOP=Zero 
  CALL ShutDown(Prog)
END PROGRAM JForce
