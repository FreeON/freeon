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
 END MODULE DrvFrcs
