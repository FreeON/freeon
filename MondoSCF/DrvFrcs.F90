!    DRIVER ROUTINES FOR DYNAMICS, OPTIMIZATION AND TS METHODS
!    Author:  Matt Challacombe and C. J. Tymczak
!------------------------------------------------------------------------------
MODULE DrvFrcs
  USE DerivedTypes
  USE GlobalScalars
#ifdef NAG
  USE F90_UNIX_PROC
#endif
  USE PrettyPrint
  USE SCFLocals
  USE Overlay
  USE ParsingKeys
  USE Functionals
  USE DrvSCFs
  USE ProcessControl    
  USE InOut
  USE AtomPairs
  USE LinALg
  IMPLICIT NONE 
  CONTAINS
!========================================================================================
!
!========================================================================================
    SUBROUTINE Forces(Ctrl,FileName_O,Unit_O)
      TYPE(SCFControls)           :: Ctrl
      CHARACTER(LEN=DCL),OPTIONAL :: FileName_O
      INTEGER,OPTIONAL            :: Unit_O
      INTEGER                     :: ICyc,IGeo,IBas,AtA,A1,PU
      LOGICAL                     :: Print
      TYPE(INT_VECT)              :: AtNum
      TYPE(DBL_VECT)              :: Frc
      INTEGER            :: Modl
!----------------------------------------------------------------------------------------
      Modl=Ctrl%Model(CBas)
!     Calculate analytic gradients piece by piece
      CALL LogSCF(Current,' Evaluating forces ')
      CtrlVect=SetCtrlVect(Ctrl,'ForceEvaluation')
!     The non-orthongonal response    
      CALL Invoke('SForce',CtrlVect)
!     Kinetic energy piece
      CALL Invoke('TForce',CtrlVect)
!     Coulomb part
      CALL Invoke('JForce',CtrlVect)
!     Exact Hartree-Fock exchange component
      IF(HasHF(Modl))CALL Invoke('XForce', CtrlVect)       
!     DFT exchange corrleation term
      IF(HasDFT(Modl))CALL Invoke('XCForce',CtrlVect)       
!---------------------------------------------------------------------------
!     If asking for just one force evaluation, we print forces to outfile
      IF(Ctrl%Grad==GRAD_ONE_FORCE)THEN
!        Get the current forces and atomic numbers 
         CALL OpenHDF(InfFile)
         CALL New(Frc,3*NAtoms)
         CALL New(AtNum,NAtoms)
         CALL Get(AtNum,'atomicnumbers',CurGeom)
         CALL Get(Frc,'GradE',CurGeom)
         CALL CloseHDF()
!        Print the gradients
         PU=OpenPU(FileName_O,Unit_O)
         WRITE(PU,33)
         WRITE(PU,*)'     Nuclear gradients for geometry #'//TRIM(CurGeom)//Rtrn
         WRITE(PU,34)
         WRITE(PU,33) 
         DO AtA = 1,NAtoms
            A1 = 3*(AtA-1)+1
            WRITE(PU,35) AtA,AtNum%I(AtA),Frc%D(A1:A1+2)
         ENDDO
         WRITE(PU,33) 
         CALL ClosePU(PU)
!        Tidy
         CALL Delete(AtNum)
         CALL Delete(Frc)
      ENDIF
   33 FORMAT(54('-'))
   34 FORMAT('  Atom      Z                Forces (au) ')
   35 FORMAT(' ',I4,'     ',I3,3(2X,F11.8))
!---------------------------------------------------
    END SUBROUTINE Forces
!========================================================================================
!
!========================================================================================
    SUBROUTINE GeOp(Ctrl)
       TYPE(SCFControls)         :: Ctrl
       INTEGER                   :: JGeo,KS
       REAL(DOUBLE)              :: RMSGrad,MaxGrad
       REAL(DOUBLE),DIMENSION(3) :: PrevState,CrntState
!-------------------------------------------------------
       CALL SetGlobalCtrlIndecies(Ctrl)           
       IF(CBas==1)THEN
!         Initialize file for geometry updates     
          CALL OpenASCII(GeoFile,Geo,NewFile_O=.TRUE.)
          CLOSE(UNIT=Geo,STATUS='KEEP')
       ENDIF
!      Compute starting energy
       CALL OneSCF(Ctrl)
!      Evaluate starting forces
       CALL Forces(Ctrl)
!      First Hessian update
       CALL Invoke('BFGSHess',CtrlVect)
       PrevState=Ctrl%Current
       Ctrl%Previous=Ctrl%Current
       Ctrl%Current=(/0,CBas,CGeo+1/)
       CrntState=Ctrl%Current
       CALL SetGlobalCtrlIndecies(Ctrl)           
!      Set starting stepsize to unity
       CALL OpenHDF(InfFile)
       CALL Put(One,'StepSize')
       CALL CloseHDF()
!      Take a step
       CtrlVect=SetCtrlVect(Ctrl,'QuNew')
       CALL Invoke('NewStep',CtrlVect)
!      Loop over geometries
       JGeo=CGeo
       DO WHILE(JGeo<=Ctrl%NGeom)
!         Evaluate an energy
          CALL OneSCF(Ctrl,Sum_O=.FALSE.)
!         Set to last succesful (minimizing) state
          Ctrl%Previous=PrevState
!         Evaluate forces 
          CALL Forces(Ctrl)
          KS=KeepStep(Ctrl)
!         Decide if we will keep this step, and check convergence
          IF(KS==1)THEN
!            Update Hessian
             CALL Invoke('BFGSHess',CtrlVect)
!            Write SCF statistics for minimizing geometry
             CALL SCFSummry(Ctrl)
!            Update status
             JGeo=JGeo+1
             Ctrl%Previous=Ctrl%Current
             PrevState=Ctrl%Current
             Ctrl%Current=(/0,CBas,JGeo/)
             CrntState=Ctrl%Current 
             CALL SetGlobalCtrlIndecies(Ctrl)           
          ELSEIF(KS==0)THEN
             Ctrl%Current=CrntState
             Ctrl%Previous=PrevState
             CALL SetGlobalCtrlIndecies(Ctrl)    
          ELSEIF(KS==-1)THEN
             RETURN
          ENDIF
!         Take another step 
          CtrlVect=SetCtrlVect(Ctrl,'QuNew')
          CALL Invoke('NewStep',CtrlVect)
       ENDDO
       CLOSE(Geo)
   END SUBROUTINE GeOp
!==============================================================================
!     Simple, very naive backtracking 
!==============================================================================
      FUNCTION KeepStep(Ctrl)
         INTEGER            :: KeepStep
         TYPE(SCFControls)  :: Ctrl
         TYPE(CRDS)         :: GM
         REAL(DOUBLE)       :: StepSz,E0,E1,RelErrE,EUncert,RMSGrad,MAXGrad
         LOGICAL            :: RelCnvrgd,RMSCnvrgd,MAXCnvrgd
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
!---------------------------------------------------------------------------
!        Open the InfFile
         CALL OpenHDF(Ctrl%Info)
!        Get the current step size
         CALL Get(StepSz,'StepSize')
!        Get the previous and current total energy
         CALL Get(E0,'Etot',Tag_O=StatsToChar(Ctrl%Previous))
         CALL Get(E1,'Etot',Tag_O=StatsToChar(Ctrl%Current))
         CALL Get(RMSGrad,'RMSGrad',CurGeom)
         CALL Get(MaxGrad,'MaxGrad',CurGeom)
!         WRITE(*,*)'E'//TRIM(StatsToChar(Ctrl%Previous))//' = '//TRIM(DblToChar(E0))
!         WRITE(*,*)'E'//TRIM(StatsToChar(Ctrl%Current))//' = '//TRIM(DblToChar(E1))
!        Check for decreasing energies and need for backtracking
         RelErrE=(E1-E0)/ABS(E0)
!        Check for convergence
         RelCnvrgd=RelErrE<1D1*ETol(Ctrl%AccL(CBas))
         RMSCnvrgd=RMSGrad<1D3*ETol(Ctrl%AccL(CBas))
         MAXCnvrgd=MAXGrad<1D3*ETol(Ctrl%AccL(CBas))
         IF(RelCnvrgd.AND.RMSCnvrgd.AND.MAXCnvrgd)THEN
            KeepStep=-1
            Mssg='Converged: dE = '//TRIM(DblToShrtChar(RelErrE)) &
              //', RMS = '//TRIM(DblToShrtChar(RMSGrad)) &
              //', MAX = '//TRIM(DblToShrtChar(MAXGrad)) &
              //', Step = '//TRIM(DblToShrtChar(StepSz))
         ELSEIF(RelErrE<Zero)THEN
!           If energy is going down, keep taking unit steps
            StepSz=One 
!           Dont back track
            KeepStep=1
            Mssg='GeomOpt: dE = '//TRIM(DblToShrtChar(RelErrE)) &
               //', RMS = '//TRIM(DblToShrtChar(RMSGrad)) &
               //', MAX = '//TRIM(DblToShrtChar(MAXGrad)) &
               //', Step = '//TRIM(DblToShrtChar(StepSz))
         ELSE
!           Decrease step by half (could get fancy here and try interpolation...)
            StepSz=StepSz*Half
            IF(StepSz<1.D-3) CALL MondoHalt(DRIV_ERROR,                    &
                                            ' Reached max resolution = '   &
                                           //TRIM(DblToShrtChar(RelErrE))  &
                                           //' in minimization of energy')
!           Mark for backtracking
            KeepStep=0
            Mssg='Baktrak: dE = '//TRIM(DblToShrtChar(RelErrE)) &
               //', RMS = '//TRIM(DblToShrtChar(RMSGrad)) &
               //', MAX = '//TRIM(DblToShrtChar(MAXGrad)) &
               //', Step = '//TRIM(DblToShrtChar(StepSz))
         ENDIF
         CALL Put(StepSz,'StepSize')
         WRITE(*,*)TRIM(Mssg)             
         IF(PrintFlags%Key>DEBUG_MINIMUM)THEN
            CALL OpenASCII(OutFile,Out)
            WRITE(Out,*)TRIM(Mssg)             
            CLOSE(Out)
         ENDIF
         IF(ABS(KeepStep)==1)THEN
            CALL Get(GM,CurGeom)
            GM%Confg=CGeo
            CALL PPrint(GM,GeoFile,Geo,'XYZ')
         ENDIF
         CALL CloseHDF()
      END FUNCTION KeepStep
!
#ifdef DEVELOPMENT

!========================================================================================
!
!========================================================================================
    SUBROUTINE NextMDGeometry(Ctrl)
      TYPE(SCFControls)            :: Ctrl
      INTEGER                      :: ICyc,IBas,IGeo,AtA,A1,I,J,K
      TYPE(CRDS)                   :: GM,GMOLD,GMNEW
      TYPE(DBL_VECT)               :: Frc
      REAL(DOUBLE)                 :: DELT,DELT2,D1MASS,D2MASS,E_MD_KIN,E_TOT_EL,VSCALE
      REAL(DOUBLE)                 :: MSS,X,Y,Z,VX,VY,VZ,AngM_X,AngM_Y,AngM_Z
      REAL(DOUBLE)                 :: RX,RY,RZ,MTOT,Volume
      REAL(DOUBLE),DIMENSION(3,3)  :: PTensor
      CHARACTER(LEN=50)            :: MDFile
!
      ICyc=Ctrl%Current(1)
      IBas=Ctrl%Current(2)
      IGeo=Ctrl%Current(3)
      CtrlVect=SetCtrlVect(Ctrl,'MD')
!
      DELT   = Ctrl%MDControls(1)
      VSCALE = Ctrl%MDControls(2)
!
      DELT2 = DELT*DELT
      CALL New(Frc,3*NAtoms)
!
      IF(IGeo==1) THEN
         CALL OpenHDF(InfFile)
         CALL Get(GM     ,Tag_O=IntToChar(IGeo))  
         CALL Get(GMNEW  ,Tag_O=IntToChar(IGeo+1))
         CALL Get(Frc    ,'GradE' ,Tag_O=IntToChar(IGeo))
         CALL CloseHDF()
         DO AtA=1,NAtoms
            A1  = 3*(AtA-1)+1
            D1MASS =  DELT/(GM%AtMss%D(AtA))
            D2MASS = DELT2/(GM%AtMss%D(AtA))
!
            GM%Vects%D(1,AtA)=VSCALE*GM%Vects%D(1,AtA)       
            GM%Vects%D(2,AtA)=VSCALE*GM%Vects%D(2,AtA)
            GM%Vects%D(3,AtA)=VSCALE*GM%Vects%D(3,AtA)
!
            GMNEW%Carts%D(1,AtA)=GM%Carts%D(1,AtA)+DELT*GM%Vects%D(1,AtA)-Half*D2MASS*Frc%D(A1)           
            GMNEW%Carts%D(2,AtA)=GM%Carts%D(2,AtA)+DELT*GM%Vects%D(2,AtA)-Half*D2MASS*Frc%D(A1+1)
            GMNEW%Carts%D(3,AtA)=GM%Carts%D(3,AtA)+DELT*GM%Vects%D(3,AtA)-Half*D2MASS*Frc%D(A1+2)
!
            GMNEW%Vects%D(1,AtA)=GM%Vects%D(1,AtA)-Half*D1MASS*Frc%D(A1)            
            GMNEW%Vects%D(2,AtA)=GM%Vects%D(2,AtA)-Half*D1MASS*Frc%D(A1+1)
            GMNEW%Vects%D(3,AtA)=GM%Vects%D(3,AtA)-Half*D1MASS*Frc%D(A1+2)
#ifdef PERIODIC
            CALL AtomCyclic(GMNEW,GMNEW%Carts%D(:,AtA))
            CALL FracCyclic(GMNEW,GMNEW%BoxCarts%D(:,AtA))
#endif
         ENDDO
!
      ELSEIF(IGeo .LT. Ctrl%NGeom) THEN
         CALL OpenHDF(InfFile)
         CALL Get(GM     ,Tag_O=IntToChar(IGeo))
         CALL Get(GMNEW  ,Tag_O=IntToChar(IGeo+1))
         CALL Get(Frc    ,'GradE' ,Tag_O=IntToChar(IGeo))
         CALL CloseHDF()
         DO AtA=1,NAtoms
            A1  = 3*(AtA-1)+1
            D1MASS =  DELT/GM%AtMss%D(AtA)
            D2MASS = DELT2/GM%AtMss%D(AtA)
!
            GM%Vects%D(1,AtA)=GM%Vects%D(1,AtA)-Half*D1MASS*Frc%D(A1)            
            GM%Vects%D(2,AtA)=GM%Vects%D(2,AtA)-Half*D1MASS*Frc%D(A1+1)
            GM%Vects%D(3,AtA)=GM%Vects%D(3,AtA)-Half*D1MASS*Frc%D(A1+2)
!
            GM%Vects%D(1,AtA)=VSCALE*GM%Vects%D(1,AtA)       
            GM%Vects%D(2,AtA)=VSCALE*GM%Vects%D(2,AtA)
            GM%Vects%D(3,AtA)=VSCALE*GM%Vects%D(3,AtA)
!
            GMNEW%Carts%D(1,AtA)=GM%Carts%D(1,AtA)+DELT*GM%Vects%D(1,AtA)-Half*D2MASS*Frc%D(A1)           
            GMNEW%Carts%D(2,AtA)=GM%Carts%D(2,AtA)+DELT*GM%Vects%D(2,AtA)-Half*D2MASS*Frc%D(A1+1)
            GMNEW%Carts%D(3,AtA)=GM%Carts%D(3,AtA)+DELT*GM%Vects%D(3,AtA)-Half*D2MASS*Frc%D(A1+2)
!
            GMNEW%Vects%D(1,AtA)=GM%Vects%D(1,AtA)-Half*D1MASS*Frc%D(A1)            
            GMNEW%Vects%D(2,AtA)=GM%Vects%D(2,AtA)-Half*D1MASS*Frc%D(A1+1)
            GMNEW%Vects%D(3,AtA)=GM%Vects%D(3,AtA)-Half*D1MASS*Frc%D(A1+2)
#ifdef PERIODIC
            CALL AtomCyclic(GMNEW,GMNEW%Carts%D(:,AtA))
            CALL FracCyclic(GMNEW,GMNEW%BoxCarts%D(:,AtA))
#endif          
         ENDDO
      ELSEIF(IGeo .EQ. Ctrl%NGeom) THEN
         CALL OpenHDF(InfFile)
         CALL Get(GM     ,Tag_O=IntToChar(IGeo))
         CALL Get(Frc    ,'GradE' ,Tag_O=IntToChar(IGeo))
         CALL CloseHDF()
         DO AtA=1,NAtoms
            A1  = 3*(AtA-1)+1
            D1MASS =  DELT/GM%AtMss%D(AtA)
            D2MASS = DELT2/GM%AtMss%D(AtA)
!
            GM%Vects%D(1,AtA)=GM%Vects%D(1,AtA)-Half*D1MASS*Frc%D(A1)            
            GM%Vects%D(2,AtA)=GM%Vects%D(2,AtA)-Half*D1MASS*Frc%D(A1+1)
            GM%Vects%D(3,AtA)=GM%Vects%D(3,AtA)-Half*D1MASS*Frc%D(A1+2)
!
            GM%Vects%D(1,AtA)=VSCALE*GM%Vects%D(1,AtA)       
            GM%Vects%D(2,AtA)=VSCALE*GM%Vects%D(2,AtA)
            GM%Vects%D(3,AtA)=VSCALE*GM%Vects%D(3,AtA)
#ifdef PERIODIC
            CALL AtomCyclic(GMNEW,GMNEW%Carts%D(:,AtA))
            CALL FracCyclic(GMNEW,GMNEW%BoxCarts%D(:,AtA))
#endif    
         ENDDO
      ENDIF 
!
      CALL OpenHDF(InfFile)
      CALL Put(GM     ,Tag_O=IntToChar(IGeo))
      CALL Put(GMNEW  ,Tag_O=IntToChar(IGeo+1))
      CALL CloseHDF()
      CALL Delete(Frc)

!
!     Calculate Center of Mass
!
      RX   = Zero
      RY   = Zero
      RZ   = Zero
      MTOT = Zero
      DO AtA=1,NAtoms
         RX   = RX   + GM%AtMss%D(AtA)*GM%Carts%D(1,AtA)
         RY   = RY   + GM%AtMss%D(AtA)*GM%Carts%D(2,AtA)
         RZ   = RZ   + GM%AtMss%D(AtA)*GM%Carts%D(3,AtA)
         MTOT = MTOT + GM%AtMss%D(AtA)
      ENDDO
      RX = RX/MTOT
      RY = RY/MTOT
      RZ = RZ/MTOT
!
!     Calculate Angular Momentum and Kinectic Energy
!
      E_MD_KIN = Zero
      AngM_X = Zero
      AngM_Y = Zero
      AngM_Z = Zero
      DO AtA=1,NAtoms
         X   = GM%Carts%D(1,AtA)-RX
         Y   = GM%Carts%D(2,AtA)-RY
         Z   = GM%Carts%D(3,AtA)-RZ
         VX  = GM%Vects%D(1,AtA)
         VY  = GM%Vects%D(2,AtA)
         VZ  = GM%Vects%D(3,AtA)
         MSS = GM%AtMss%D(AtA)
         E_MD_KIN = E_MD_KIN+Half*MSS*(VX*VX+VY*VY+VZ*VZ)
         AngM_X = AngM_X+MSS*(Y*VZ-Z*VY)
         AngM_Y = AngM_Y+MSS*(Z*VX-X*VZ)
         AngM_Z = AngM_Z+MSS*(X*VY-Y*VX)
      ENDDO
#ifdef PERIODIC
!
!     Calculate the Presure tensor (Need Work of 1 & 2 dim)
!
      PTensor = Zero 
!
      Volume = One
      DO I=1,3;IF(GM%AutoW(I)) Volume = Volume*GM%BoxShape%D(I,I);ENDDO


      DO AtA=1,NAtoms
         A1  = 3*(AtA-1)+1
         MSS = GM%AtMss%D(AtA)
         DO I=1,3
            DO J=1,3     
               IF(GM%AutoW(I) .AND. GM%AutoW(J) ) THEN
                  PTensor(I,J) = PTensor(I,J) + MSS*GM%Vects%D(I,AtA)*GM%Vects%D(J,AtA) &
                                              + GM%Carts%D(I,AtA)*Frc%D(A1+J-1)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      PTensor = PTensor/Volume
!
!     Store the Pressure Tensor
!
      MDFile   = 'MDPTensor_'// TRIM(SCF_NAME)// "_" // TRIM(PROCESS_ID) // '.out'
      OPEN(UNIT=30,FILE=MDFile,POSITION='APPEND',STATUS='UNKNOWN')
      WRITE(30,11) IGeo
      WRITE(30,12) PTensor(1,1),PTensor(1,2),PTensor(1,3)
      WRITE(30,12) PTensor(2,1),PTensor(2,2),PTensor(2,3)
      WRITE(30,12) PTensor(3,1),PTensor(3,2),PTensor(3,3)
#endif
!
!     Store the Energies and Angular Momentum
!
      CALL OpenHDF(InfFile)
      CALL Get(E_TOT_EL,'E_total',Tag_O=IntToChar(ICyc))
      CALL CloseHDF()
!
      MDFile   = 'MDEnergy_'// TRIM(SCF_NAME)// "_" // TRIM(PROCESS_ID) // '.out'
      OPEN(UNIT=30,FILE=MDFile,POSITION='APPEND',STATUS='UNKNOWN')
      WRITE(30,10) IGeo,E_MD_KIN,E_TOT_EL,E_MD_KIN+E_TOT_EL,AngM_X,AngM_Y,AngM_Z
      CLOSE(30)
10    FORMAT(2X,I6,2X,F22.16,2X,F22.16,2X,F22.16,4X,F14.8,2X,F14.8,2X,F14.8)
11    FORMAT(2X,I6)
12    FORMAT(2X,F22.16,2X,F22.16,2X,F22.16)
!
    END SUBROUTINE NextMDGeometry


!========================================================================================
!
!========================================================================================
    SUBROUTINE Force_Num_1(Ctrl)
      TYPE(SCFControls)                :: Ctrl
      TYPE(CRDS)                       :: GM,GMOLD      
      TYPE(BCSR)                       :: S,P,F,T1,T2
      INTEGER                          :: ICyc,IGeo,IBas
      INTEGER                          :: AtA,IX,II,A1,A2,IS
      INTEGER, DIMENSION(3)            :: Begin
      REAL(DOUBLE)                     :: DDelta
      REAL(DOUBLE),DIMENSION(2)        :: ES,ET,EJ,EN,EX,EXC 
      REAL(DOUBLE),DIMENSION(NAtoms,3) :: F_S,F_T,F_J,F_N,F_XC,F_X,F_TOT
      TYPE(DBL_VECT)                   :: Frc_Num
      CHARACTER(LEN=50)                :: ForceFile
!
      ICyc=Ctrl%Current(1)
      IBas=Ctrl%Current(2)
      IGeo=Ctrl%Current(3)
      CtrlVect=SetCtrlVect(Ctrl,'Force')
      Begin = (/ICyc,IBas,IGeo/)
!
      CALL OpenHDF(InfFile)
      CALL Get(GM   ,Tag_O=IntToChar(IGeo))
      CALL Get(GMOLD,Tag_O=IntToChar(IGeo))
      CALL CloseHDF()
!
      F_T   = Zero
      F_J   = Zero
      F_N   = Zero
      F_XC  = Zero
      F_TOT = Zero
!
      DDelta = 1.D-6
      DO AtA = 1,1!Natoms
         DO IX = 1,1
            DO II=1,2
!              Change the Geometry
               IF(II==1) THEN
                  GM%Carts%D(IX,AtA) = GMOLD%Carts%D(IX,AtA)+DDelta
               ELSEIF(II==2) THEN
                  GM%Carts%D(IX,AtA) = GMOLD%Carts%D(IX,AtA)-DDelta
               ENDIF
               CALL OpenHDF(InfFile)
               CALL Put(GM,Tag_O=IntToChar(IGeo))
               CALL CloseHDF()
!======================================================================================
!              SCF
! 
               DO IS = 1,3
                  Ctrl%Current(1) = Ctrl%Current(1)+1
                  IF(IS .EQ. 1) CALL OneEMats(Ctrl)
                  CALL SCFCycle(Ctrl,Begin)
                  IF(ConvergedQ(Ctrl))GOTO 999
999               CONTINUE 
               ENDDO
!======================================================================================
!              GET ENERGIES
!
               CALL OpenHDF(InfFile)
               CALL Get(ET(II), 'Kin'     ,Tag_O=IntToChar(Ctrl%Current(1)))
               CALL Get(EJ(II), 'ene+eee' ,Tag_O=IntToChar(Ctrl%Current(1)))
               CALL Get(EN(II), 'enn+ene' ,Tag_O=IntToChar(Ctrl%Current(1)))
               CALL Get(EXC(II),'Exc'     ,Tag_O=IntToChar(Ctrl%Current(1)))
               CALL CloseHDF() 
            ENDDO
            F_T(AtA,IX)  = (ET(1)-ET(2))/(Two*DDelta)
            F_J(AtA,IX)  = (EJ(1)-EJ(2))/(Two*DDelta)
            F_N(AtA,IX)  = (EN(1)-EN(2))/(Two*DDelta)
            F_XC(AtA,IX) = (EXC(1)-EXC(2))/(Two*DDelta)
            F_TOT(AtA,IX) =  F_T(AtA,IX)+F_J(AtA,IX)+F_N(AtA,IX)+F_XC(AtA,IX)
         ENDDO
      ENDDO      
!
      CALL New(Frc_Num,3*NAtoms)
      DO AtA = 1,NAtoms
         A1=3*(AtA-1)+1
         Frc_Num%D(A1)   = F_TOT(AtA,1)
         Frc_Num%D(A1+1) = F_TOT(AtA,2)
         Frc_Num%D(A1+2) = F_TOT(AtA,3)
      ENDDO
!
      ForceFile   = 'Force_Num_'// TRIM(SCF_NAME)// "_" // TRIM(PROCESS_ID) // '.out'
      OPEN(UNIT=30,FILE=ForceFile,POSITION='APPEND',STATUS='UNKNOWN')
      WRITE(30,12) IGeo
      DO AtA = 1,NAtoms
         A1 = 3*(AtA-1)+1
         WRITE(30,11) AtA,GM%Carts%D(:,AtA),Frc_Num%D(A1:A1+2)
      ENDDO
      CLOSE(30)
      CALL Delete(Frc_Num)
!
11    FORMAT('ATOM#',I6,2X,D20.12,2X,D20.12,2X,D20.12,2X,D20.12,2X,D20.12,2X,D20.12)
12    FORMAT('Geometry #',I6)
!
    END SUBROUTINE Force_Num_1
!========================================================================================
!
!========================================================================================
    SUBROUTINE Force_Num_2(Ctrl)
      TYPE(SCFControls)                :: Ctrl
      TYPE(CRDS)                       :: GM,GMOLD
      TYPE(BCSR)                       :: S,P,F,T1,T2
      INTEGER                          :: ICyc,IGeo,IBas
      INTEGER                          :: AtA,IX,II,A1,A2,IS
      REAL(DOUBLE)                     :: DDelta
      REAL(DOUBLE),DIMENSION(2)        :: ES,ET,EJ,EN,EX,EXC
      REAL(DOUBLE),DIMENSION(NAtoms,3) :: F_S,F_T,F_J,F_N,F_X,F_XC,F_TOT
      TYPE(DBL_VECT)                   :: Frc_Num
!
      ICyc=Ctrl%Current(1)
      IBas=Ctrl%Current(2)
      IGeo=Ctrl%Current(3)
      CtrlVect=SetCtrlVect(Ctrl,'Direct')
!
      CALL OpenHDF(InfFile)
      CALL Get(GM   ,Tag_O=IntToChar(IGeo))
      CALL Get(GMOLD,Tag_O=IntToChar(IGeo))
      CALL CloseHDF()
!
!     Load globals 
!
      CALL OpenHDF(InfFile)
      CALL Get(MaxAtms,  'maxatms',   Tag_O=IntToChar(IBas))
      CALL Get(MaxBlks,  'maxblks',   Tag_O=IntToChar(IBas))
      CALL Get(MaxNon0,  'maxnon0',   Tag_O=IntToChar(IBas))
      CALL Get(NBasF,    'nbasf',     Tag_O=IntToChar(IBas))
      CALL New(BSiz,NAtoms)
      CALL New(OffS,NAtoms)
      CALL Get(BSiz,'atsiz',          Tag_O=IntToChar(IBas))
      CALL Get(OffS,'atoff',          Tag_O=IntToChar(IBas))
      MaxBlkSize=0
      DO II=1,NAtoms; MaxBlkSize=MAX(MaxBlkSize,BSiz%I(II)); ENDDO
      CALL SetThresholds(IntToChar(IBas))
      CALL CloseHDF()
!
!     Get Matices
!
      CALL New(S)
      CALL New(F)
      CALL New(P)
      CALL New(T1)
      CALL New(T2)
!      
      CALL OpenHDF(InfFile)
      CALL Get(ModelChem,'ModelChemistry',Tag_O=IntToChar(IBas))
!
      CALL Get(P,TrixFile('D',OffSet_O=0,Name_O=CtrlVect(1),Stats_O=Ctrl%Current(1:3)))
      CALL Get(F,TrixFile('F',OffSet_O=0,Name_O=CtrlVect(1),Stats_O=Ctrl%Current(1:3)))
      CALL Get(S,TrixFile('S',Name_O=CtrlVect(1),Stats_O=Ctrl%Current(1:3)))
      CALL CloseHDF()
!
      DDelta = 1.D-3
      DO AtA = 1,1!NAtoms
         DO IX = 1,1 !3
            DO II=1,2
!              Change the Geometry
               IF(II==1) THEN
                  GM%Carts%D(IX,AtA) = GMOLD%Carts%D(IX,AtA)+DDelta
               ELSEIF(II==2) THEN
                  GM%Carts%D(IX,AtA) = GMOLD%Carts%D(IX,AtA)-DDelta
               ENDIF
               CALL OpenHDF(InfFile)
               CALL Put(GM,Tag_O=IntToChar(IGeo))
               CALL CloseHDF()
!======================================================================================
!              Remake the Matrices
!               
               CALL Invoke('MakeS'  ,   CtrlVect)
               CALL Invoke('MakeT'  ,   CtrlVect)
               CALL Invoke('MakeRho',   CtrlVect)
               CALL Invoke('QCTC'   ,   CtrlVect)
               IF(HasDFT(ModelChem)) &
               CALL Invoke('HiCu'   ,   CtrlVect)       
               IF(HasHF(ModelChem))  &
               CALL Invoke('ONX'    ,   CtrlVect)       
!======================================================================================
!              Calculate the Overlap Term
! 
               CALL Get(S,TrixFile('S',Name_O=CtrlVect(1),Stats_O=Ctrl%Current(1:3)))
               CALL Multiply(P,F,T1)       
               CALL Multiply(T1,P,T2)     
               ES(II) = -Two*Trace(T2,S)
!======================================================================================
!              Recalculate the Energies
! 
               CALL Invoke('SCFstats',  CtrlVect,MPIRun_O=.TRUE.)   
!======================================================================================
!              GET ENERGIES
               CALL OpenHDF(InfFile)
               CALL Get(ET(II), 'Kin'     ,Tag_O=IntToChar(Ctrl%Current(1)))
               CALL Get(EJ(II), 'ene+eee' ,Tag_O=IntToChar(Ctrl%Current(1)))
               CALL Get(EN(II), 'enn+ene' ,Tag_O=IntToChar(Ctrl%Current(1)))
               IF(HasHF(ModelChem))  &
               CALL Get(EX(II), 'Ex'      ,Tag_O=IntToChar(Ctrl%Current(1)))
               IF(HasDFT(ModelChem)) &
               CALL Get(EXC(II),'Exc'     ,Tag_O=IntToChar(Ctrl%Current(1)))
               CALL CloseHDF() 
            ENDDO
            F_S(AtA,IX)  = (ES(1)-ES(2))/(Two*DDelta)
            F_T(AtA,IX)  = (ET(1)-ET(2))/(Two*DDelta)
            F_J(AtA,IX)  = (EJ(1)-EJ(2))/(Two*DDelta)
            F_N(AtA,IX)  = (EN(1)-EN(2))/(Two*DDelta)
            F_X(AtA,IX)  = (EX(1)-EX(2))/(Two*DDelta)
            F_XC(AtA,IX) = (EXC(1)-EXC(2))/(Two*DDelta)
            F_TOT(AtA,IX)= F_S(ATA,IX)+ F_T(AtA,IX)+F_J(AtA,IX)+F_N(AtA,IX)+F_X(AtA,IX)+F_XC(AtA,IX)
         ENDDO
      ENDDO
!     RESET Geometry
      CALL OpenHDF(InfFile)
      CALL Put(GMOLD,Tag_O=IntToChar(IGeo))
      CALL CloseHDF()

      WRITE(*,*)' --------------------------------------------------'
      WRITE(*,*)'SForce = ',F_S(1,1)
      WRITE(*,*)'TForce = ',F_T(1,1)
      WRITE(*,*)'JForce = ',F_J(1,1)+F_N(1,1)
      WRITE(*,*)'XForce = ',F_X(1,1)     
      WRITE(*,*)'XCForce= ',F_XC(1,1)     
      WRITE(*,*) ' '
      WRITE(*,*)'Total Force = ',F_TOT(1,1)
      WRITE(*,*)' --------------------------------------------------'
!
      CALL New(Frc_Num,3*NAtoms)
      DO AtA = 1,NAtoms
         A1=3*(AtA-1)+1
         Frc_Num%D(A1)   = F_TOT(AtA,1)
         Frc_Num%D(A1+1) = F_TOT(AtA,2)
         Frc_Num%D(A1+2) = F_TOT(AtA,3)
      ENDDO
!      
      CALL OpenHDF(InfFile)
      CALL Put(Frc_Num,'GradE_Num',Tag_O=IntToChar(IGeo))
      CALL CloseHDF()
!
    END SUBROUTINE Force_Num_2

#endif

 END MODULE DrvFrcs
