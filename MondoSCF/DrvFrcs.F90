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
      TYPE(SCFControls)              :: Ctrl
      CHARACTER(LEN=DCL),OPTIONAL    :: FileName_O
      INTEGER,OPTIONAL               :: Unit_O
      INTEGER                        :: ICyc,IGeo,IBas
      TYPE(CRDS)                     :: GM,GM_MM
      TYPE(DBL_VECT)                 :: Frc
      INTEGER                        :: Modl
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
!----------------------------------------------------------------------------------------
#ifdef MMech
      IF(HasQM()) THEN
#endif
        Modl=Ctrl%Model(CBas)
!       Calculate analytic gradients piece by piece
        CALL LogSCF(Current,' Evaluating forces ')
        CtrlVect=SetCtrlVect(Ctrl,'ForceEvaluation')
!       The non-orthogonal response    
        CALL Invoke('SForce',CtrlVect)
!       Kinetic energy piece
        CALL Invoke('TForce',CtrlVect)
!       Build a density with last DM
        CALL Invoke('MakeRho',CtrlVect)
!       Coulomb part
        CALL Invoke('JForce',CtrlVect)
!       Exact Hartree-Fock exchange component
        IF(HasHF(Modl))CALL Invoke('XForce', CtrlVect)       
!       DFT exchange corrleation term
        IF(HasDFT(Modl))CALL Invoke('XCForce',CtrlVect)       
#ifdef MMech
      ENDIF
#endif
!---------------------------------------------------------------------------
!     If asking for just one force evaluation, we print forces to outfile
      IF(Ctrl%Grad==GRAD_ONE_FORCE)THEN
!        Get the current forces and atomic numbers 
#ifdef MMech
           IF(QMOnly()) THEN
#endif
             CALL Get(GM,Tag_O=CurGeom)
#ifdef MMech
           ELSE 
             CALL Get(GM,Tag_O='GM_MM'//CurGeom)
           ENDIF
#endif
           CALL New(Frc,3*GM%NAtms)
           CALL Get(Frc,'GradE',CurGeom)
!        Print the gradients
         Frc%D = -Frc%D/KJPerMolPerAngstToHPerBohr !!!in KJ/mol/Angstroem
         Mssg = '     Nuclear gradients for geometry #'//TRIM(CurGeom)//Rtrn
         CALL Print_Force(GM,Frc,Mssg,FileName_O,Unit_O,Fmat_O=2)
!        Tidy
         CALL Delete(Frc)
         CALL Delete(GM)
      ENDIF
!
    END SUBROUTINE Forces
!========================================================================================
!
!========================================================================================
    SUBROUTINE GeOp(Ctrl)
       TYPE(SCFControls)         :: Ctrl
       TYPE(CRDS)                :: GM
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
!  Compute starting energy and forces for the 
!  actual basis set and geometry
#ifdef MMech
       IF(HasQM()) Then
         IF(HasMM()) CALL MM_ENERG(Ctrl)
         CALL OneSCF(Ctrl)
         CALL Forces(Ctrl)
       ENDIF
!
       IF(MMOnly()) Then
         Ctrl%Current=(/0,0,1/)
         CALL SetGlobalCtrlIndecies(Ctrl)           
         CALL MM_ENERG(Ctrl)
       ENDIF
#else
       CALL OneSCF(Ctrl)
       CALL Forces(Ctrl)
#endif
!  Print first configuration
#ifdef MMech
       IF(QMOnly()) THEN
         CALL Get(GM,CurGeom)
         CALL Get(GM%ETotal,'Etot',Tag_O=StatsToChar(Ctrl%Current))
         GM%Confg=CGeo
         CALL PPrint(GM,GeoFile,Geo,'PDB')
         CALL Delete(GM)
       ELSE 
         CALL Get(GM_MM,'GM_MM'//CurGeom)
         CALL Get(GM_MM%ETotal,'Etot',Tag_O=StatsToChar(Ctrl%Current))
         GM_MM%Confg=CGeo
         CALL PPrint(GM_MM,GeoFile,Geo,'PDB')
         CALL Delete(GM_MM)
       ENDIF
#else
       CALL Get(GM,CurGeom)
       CALL Get(GM%ETotal,'Etot',Tag_O=StatsToChar(Ctrl%Current))
       GM%Confg=CGeo
       CALL PPrint(GM,GeoFile,Geo,'PDB')
       CALL Delete(GM)
#endif
!  First Hessian update
!      CtrlVect=SetCtrlVect(Ctrl,'BFGSHess') !!!no CtrlVect definition was given before, in front of Invoke('BFGSHess'
       CALL Invoke('BFGSHess',CtrlVect)
       PrevState=Ctrl%Current
       Ctrl%Previous=Ctrl%Current
       Ctrl%Current=(/0,CBas,CGeo+1/) ! step geometry count
       CrntState=Ctrl%Current
       CALL SetGlobalCtrlIndecies(Ctrl)           
       CALL Put(One,'StepSize')
!  Take a step 
       CtrlVect=SetCtrlVect(Ctrl,'QuNew',NoAdvance_O=.TRUE.)
       CALL Invoke('NewStep',CtrlVect)
!  Loop over geometries
       JGeo=CGeo
       DO WHILE(JGeo<=Ctrl%NGeom)
!  Evaluate energies and forces
#ifdef MMech
       IF(HasQM()) Then
         IF(HasMM()) CALL MM_ENERG(Ctrl)
         CALL OneSCF(Ctrl,Sum_O=.FALSE.)
         Ctrl%Previous=PrevState
         CALL Forces(Ctrl)
       ENDIF
!
       IF(MMOnly()) Then
         Ctrl%Current=(/0,0,JGeo/)
         CALL SetGlobalCtrlIndecies(Ctrl)           
         CALL MM_ENERG(Ctrl)
       ENDIF
#else
       CALL OneSCF(Ctrl,Sum_O=.FALSE.)
!  Set to last succesful (minimizing) state
       Ctrl%Previous=PrevState
!  Evaluate forces 
       CALL Forces(Ctrl)
#endif
          KS=KeepStep(Ctrl)
!  Decide if we will keep this step, and check convergence
          IF(KS==1)THEN
!  Update Hessian
             CALL Invoke('BFGSHess',CtrlVect)
!  Write SCF statistics for minimizing geometry
#ifdef MMech
             IF(HasQM()) CALL SCFSummry(Ctrl)
#else
             CALL SCFSummry(Ctrl)
#endif
!  Update status
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
             Ctrl%Previous=Ctrl%Current
             RETURN
          ENDIF
!  Take another step 
          CtrlVect=SetCtrlVect(Ctrl,'QuNew',NoAdvance_O=.TRUE.)
          CALL Invoke('NewStep',CtrlVect)
       ENDDO
       CLOSE(Geo)
   END SUBROUTINE GeOp
!==============================================================================
!     Simple, very naive backtracking (line search)
!==============================================================================
      FUNCTION KeepStep(Ctrl)
         INTEGER            :: KeepStep
         TYPE(SCFControls)  :: Ctrl
         TYPE(CRDS)         :: GM
         INTEGER            :: AccL
         REAL(DOUBLE)       :: StepSz,OldStep,E0,E1,ET,RelErrE,EUncert, &
                               RMSGrad,MAXGrad,RMSDisp,MaxDisp
         LOGICAL            :: ECnvrgd,GCnvrgd,XCnvrgd
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
!---------------------------------------------------------------------------
!        Open the InfFile
!        Get the current step size
         CALL Get(StepSz,'StepSize')
         OldStep=StepSz
!        Get the previous and current total energy
         CALL Get(E0,'Etot',Tag_O=StatsToChar(Ctrl%Previous))
         CALL Get(E1,'Etot',Tag_O=StatsToChar(Ctrl%Current))
         CALL Get(RMSDisp,'RMSDisp',CurGeom)
         CALL Get(MaxDisp,'MaxDisp',CurGeom)
         CALL Get(RMSGrad,'RMSGrad',CurGeom)
         CALL Get(MaxGrad,'MaxGrad',CurGeom)
!         WRITE(*,*)'E'//TRIM(StatsToChar(Ctrl%Previous))//' = '//TRIM(DblToChar(E0))
!         WRITE(*,*)'E'//TRIM(StatsToChar(Ctrl%Current))//' = '//TRIM(DblToChar(E1))
!        Check for decreasing energies and need for backtracking
         RelErrE=(E1-E0)/ABS(E0)
!        Check for convergence
         AccL=Ctrl%AccL(CBas)
         ET=ETol(AccL)
         ECnvrgd=ABS(RelErrE)<1D1*ET.AND.StepSz==One
         GCnvrgd=RMSGrad<GTol(AccL).AND.MaxGrad<GTol(AccL)
         XCnvrgd=RMSDisp<XTol(AccL).AND.MaxDisp<XTol(AccL)
         IF(ECnvrgd.OR.(GCnvrgd.AND.XCnvrgd))THEN
            KeepStep=-1
            Mssg=ProcessName('QuNew','Converged #'//TRIM(CurGeom))
          ELSEIF(RelErrE<Zero)THEN
            StepSz=One 
!           Dont back track
            KeepStep=1
            Mssg=ProcessName('QuNew','Descent   #'//TRIM(CurGeom))
         ELSE
!           Decrease step by half (could get fancy here and try interpolation...)
            Mssg=ProcessName('QuNew','Backtrack #'//TRIM(CurGeom))
            StepSz=StepSz*Half
            IF(StepSz<1.D-2) CALL MondoHalt(DRIV_ERROR,                    &
                                            ' Reached max resolution = '   &
                                           //TRIM(DblToShrtChar(RelErrE))  &
                                           //' in minimization of energy')
!           Mark for backtracking
            KeepStep=0
         ENDIF
         CALL Put(StepSz,'StepSize')
         Mssg=TRIM(Mssg)//' dE= '//TRIM(DblToShrtChar(RelErrE))    &
                        //', Grms= '//TRIM(DblToShrtChar(RMSGrad)) &
                        //', Gmax= '//TRIM(DblToShrtChar(MAXGrad)) &
                        //', Xrms= '//TRIM(DblToShrtChar(RMSDisp)) &
                        //', Xmax= '//TRIM(DblToShrtChar(MAXDisp)) &
                        //', Step= '//TRIM(DblToShrtChar(OldStep))

         WRITE(*,*)TRIM(Mssg)             
         CALL OpenASCII(OutFile,Out)
         WRITE(Out,*)TRIM(Mssg)             
         CLOSE(Out)
         IF(ABS(KeepStep)==1)THEN
#ifdef MMech
           IF(QMOnly()) THEN
             CALL Get(GM,Tag_O=CurGeom)
           ELSE IF(HasMM()) THEN
             CALL Get(GM,Tag_O='GM_MM'//CurGeom)
           ENDIF
#else
           CALL Get(GM,CurGeom)
#endif
            GM%Confg=CGeo
            GM%ETotal=E1
            CALL PPrint(GM,GeoFile,Geo,'PDB')
         ENDIF
      END FUNCTION KeepStep
!========================================================================================
!
!========================================================================================
    SUBROUTINE NForce(Ctrl)
      TYPE(SCFControls)                :: Ctrl
      TYPE(CRDS)                       :: GM,GMOLD
      TYPE(BCSR)                       :: P,F,T1,T2,T3
      INTEGER                          :: ICyc,IGeo,IBas
      INTEGER                          :: AtA,IX,II,IA,A1,A2,IS,N
      REAL(DOUBLE),DIMENSION(2)        :: ES,ET,EJ,EN,EX,EXC
      TYPE(DBL_VECT)                   :: FS,FT,FJ,FC,FX,FXC,FTot
      REAL(DOUBLE),PARAMETER           :: DDelta = 5.D-4
      LOGICAL                          :: DoMat
!------------------------------------------------------------------
      DoMat = .TRUE.
!  
      CALL SetGlobalCtrlIndecies(Ctrl)
      CtrlVect=SetCtrlVect(Ctrl,'NumForceEvaluation')
      CALL Get(GM   ,Tag_O=CurGeom)
      CALL Get(GMOLD,Tag_O=CurGeom)
!     Load globals 
      CALL Get(ModelChem,'ModelChemistry',Tag_O=CurBase)
      CALL Get(MaxAtms,  'maxatms',   Tag_O=CurBase)
      CALL Get(MaxBlks,  'maxblks',   Tag_O=CurBase)
      CALL Get(MaxNon0,  'maxnon0',   Tag_O=CurBase)
      CALL Get(NBasF,    'nbasf',     Tag_O=CurBase)
      CALL New(BSiz,NAtoms)
      CALL New(OffS,NAtoms)
      CALL Get(BSiz,'atsiz',          Tag_O=CurBase)
      CALL Get(OffS,'atoff',          Tag_O=CurBase)
      MaxBlkSize=0
      DO II=1,NAtoms; MaxBlkSize=MAX(MaxBlkSize,BSiz%I(II)); ENDDO
!
      CALL New(F)
      CALL New(P)
      CALL New(T1)
      CALL New(T2)
      CALL New(T3)
!
      CALL New(FS,3*NAtoms)
      CALL New(FT,3*NAtoms)
      CALL New(FC,3*NAtoms)
      CALL New(FX,3*NAtoms)
      CALL New(FXC,3*NAtoms)
      CALL New(FTot,3*NAtoms)
!
      FS%D=Zero
      FT%D=Zero
      FC%D=Zero
      FX%D=Zero
      FXC%D=Zero
      FTot%D=Zero
!
      IF(DoMat) THEN
!        Get Matices for computing orthogonalizaiton term
         CALL Get(F,TrixFile('F',OffSet_O=0,Name_O=CtrlVect(1),Stats_O=Current))
         Ctrl%Current(1)=Ctrl%Current(1)+1
         CALL SetGlobalCtrlIndecies(Ctrl)
         CtrlVect=SetCtrlVect(Ctrl,'NumForceEvaluation')
         CALL Get(P,TrixFile('D',OffSet_O=0,Name_O=CtrlVect(1),Stats_O=Current))
!
         DO AtA=1,NAtoms
            DO IX=1,3
               DO II=1,2
                  WRITE(*,*) 'AtA = ',AtA,' IX = ',IX,' II = ',II
!                 Change the Geometry
                  IF(II==1) THEN
                     GM%Carts%D(IX,AtA) = GMOLD%Carts%D(IX,AtA)+DDelta
                  ELSEIF(II==2) THEN
                     GM%Carts%D(IX,AtA) = GMOLD%Carts%D(IX,AtA)-DDelta
                  ENDIF
                  CALL Put(GM,Tag_O=CurGeom)
!                 Calculate the Overlap Term
                  CALL Invoke('MakeS'  ,   CtrlVect)
                  CALL Get(T1,TrixFile('S',Name_O=CtrlVect(1),Stats_O=Ctrl%Current(1:3)))
                  CALL Multiply(P,F,T2)       
                  CALL Multiply(T2,P,T3)     
                  ES(II) = -Two*Trace(T3,T1)
!                 Calcualte the Kinectic Term
                  CALL Invoke('MakeT'  ,   CtrlVect)
                  CALL Get(T1,TrixFile('T',Name_O=CtrlVect(1),Stats_O=Ctrl%Current(1:3)))
                  ET(II) =  Two*Trace(P,T1)        
!                 Calculate the Coulomb Term d(Vee + Vne)/dR
                  CALL Invoke('QCTC'   ,   CtrlVect)
                  CALL Get(T1,TrixFile('J',OffSet_O=0,Name_O=CtrlVect(1),Stats_O=Ctrl%Current(1:3)))
                  EJ(II) =  Two*Trace(P,T1)
!                 Calculate the Coulomb Term d(Vne + Vnn)/dR
                  CALL Get(EN(II), 'E_NuclearTotal'    ,Tag_O=IntToChar(Ctrl%Current(1)))
                  EN(II) =  Two*EN(II)
!                 Calculate the x term
                  IF(HasHF(ModelChem))  THEN
                     CALL Invoke('ONX'    ,   CtrlVect) 
                     CALL Get(T1,TrixFile('K',OffSet_O=0,Name_O=CtrlVect(1),Stats_O=Ctrl%Current(1:3)))
                     EX(II) = Trace(P,T1)   
                  ENDIF
!                 Calculate the xc term
                  IF(HasDFT(ModelChem)) THEN
                     CALL Invoke('HiCu'   ,   CtrlVect)
                     CALL Get(T1,TrixFile('Kxc',OffSet_O=0,Name_O=CtrlVect(1),Stats_O=Ctrl%Current(1:3)))
                     EXC(II) = Two*Trace(P,T1)
                  ENDIF
!                 Reset geometry
                  GM%Carts%D(IX,AtA) = GMOLD%Carts%D(IX,AtA)
                  CALL Put(GM,Tag_O=CurGeom)
               ENDDO
               IA=3*(AtA-1)+IX
               FS%D(IA)=(ES(1)-ES(2))/(Two*DDelta)
               FT%D(IA)=(ET(1)-ET(2))/(Two*DDelta)
               FC%D(IA)=(EJ(1)-EJ(2))/(Two*DDelta)+(EN(1)-EN(2))/(Two*DDelta)
               FX%D(IA)=(EX(1)-EX(2))/(Two*DDelta)
               FXC%D(IA)=(EXC(1)-EXC(2))/(Two*DDelta)
               FTot%D(IA)=FS%D(IA)+FT%D(IA)+FC%D(IA)+FX%D(IA)+FXC%D(IA)
            ENDDO
         ENDDO
!        Do some checksumming and IO 
         CALL Print_Force(GMOLD,FS,'   dS/dR')
         CALL Print_Force(GMOLD,FT,'   dT/dR')
         CALL Print_Force(GMOLD,FC,'   dJ/dR')
         IF(HasHF(ModelChem))   CALL Print_Force(GMOLD,FX, '   dHF/dR')
         IF(HasDFT(ModelChem))  CALL Print_Force(GMOLD,FXC,'   dXC/dR')
         CALL Print_Force(GMOLD,FTot,'   dSCF/dR')
!   
         CALL PChkSum(FS,'dS/dR',Proc_O='NForce')
         CALL PChkSum(FT,'dT/dR',Proc_O='NForce')
         CALL PChkSum(FC,'dJ/dR',Proc_O='NForce')
         IF(HasHF(ModelChem))  &
              CALL PChkSum(FX,'dHF/dR',Proc_O='NForce')
         IF(HasDFT(ModelChem))  &
              CALL PChkSum(FXC,'dXC/dR',Proc_O='NForce')
         CALL PChkSum(FTot,'dSCF/dR',Proc_O='NForce')
      ELSE     
         ICyc = Ctrl%Current(1)
         DO AtA=1,NAtoms
            DO IX=1,3
               DO II=1,2
                  WRITE(*,*) 'AtA = ',AtA,' IX = ',IX,' II = ',II
 !                Change the Geometry
                  IF(II==1) THEN
                     GM%Carts%D(IX,AtA) = GMOLD%Carts%D(IX,AtA)+DDelta
                  ELSEIF(II==2) THEN
                     GM%Carts%D(IX,AtA) = GMOLD%Carts%D(IX,AtA)-DDelta
                  ENDIF
                  CALL Put(GM,Tag_O=CurGeom)
!                 SCF
                  DO N=1,4
                     Ctrl%Current(1)=ICyc+N
                     CALL SetGlobalCtrlIndecies(Ctrl)
                     IF(N==1) CALL OneEMats(Ctrl)
                     CALL SCFCycle(Ctrl) 
                  ENDDO
                  
!                 Get Energy     
                  CALL Get(ES(II),'Etot',StatsToChar(Ctrl%Current))
!                 Reset geometry
                  GM%Carts%D(IX,AtA) = GMOLD%Carts%D(IX,AtA)
                  CALL Put(GM,Tag_O=CurGeom)
               ENDDO
               IA=3*(AtA-1)+IX
               WRITE(*,*) " Energies ",ES(1),ES(2),ES(1)-ES(2)
               FTot%D(IA)=(ES(1)-ES(2))/(Two*DDelta)
            ENDDO
         ENDDO
!        Do some checksumming and IO 
         CALL PChkSum(FTot,'dSCF/dR',Proc_O='NForce')
      ENDIF
!
      CALL Delete(GM)
      CALL Delete(GMOLD)
      CALL Delete(BSiz)
      CALL Delete(OffS)
      CALL Delete(FS)
      CALL Delete(FT)
      CALL Delete(FC)
      CALL Delete(FX)
      CALL Delete(FXC)
      CALL Delete(FTot)
!
      CALL Delete(F)
      CALL Delete(P)
      CALL Delete(T1)
      CALL Delete(T2)
      CALL Delete(T3)
!
    END SUBROUTINE NForce
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
         CALL Get(GM     ,Tag_O=IntToChar(IGeo))  
         CALL Get(GMNEW  ,Tag_O=IntToChar(IGeo+1))
         CALL Get(Frc    ,'GradE' ,Tag_O=IntToChar(IGeo))
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
         CALL Get(GM     ,Tag_O=IntToChar(IGeo))
         CALL Get(GMNEW  ,Tag_O=IntToChar(IGeo+1))
         CALL Get(Frc    ,'GradE' ,Tag_O=IntToChar(IGeo))
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
         CALL Get(GM     ,Tag_O=IntToChar(IGeo))
         CALL Get(Frc    ,'GradE' ,Tag_O=IntToChar(IGeo))
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
      CALL Put(GM     ,Tag_O=IntToChar(IGeo))
      CALL Put(GMNEW  ,Tag_O=IntToChar(IGeo+1))
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
      Volume = GM%PBC%CellVolume
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
      CALL Get(E_TOT_EL,'E_total',Tag_O=IntToChar(ICyc))
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
#endif
!------------------------------------------------------------------
  SUBROUTINE CALC_GRAD_ONE_FORCE(Ctrl)
    IMPLICIT NONE
    TYPE(SCFControls)     :: Ctrl
    INTEGER :: I,ISet,IGeo
!
#ifdef MMech
   If(HasQM()) Then
#endif
!  Loop first over basis sets 
     DO ISet=1,Ctrl%NSet
        Ctrl%Current=(/0,ISet,1/)
        CALL SetGlobalCtrlIndecies(Ctrl)           
#ifdef MMech
        IF(ISet==1.AND.HasMM()) CALL MM_ENERG(Ctrl)
#endif
        CALL OneSCF(Ctrl)
     ENDDO
!    Calculate the Force (last basis set)
     CALL Forces(Ctrl)
!
     IF(Ctrl%NGeom>1)THEN
!  Go over additional geometries at last basis set
        DO IGeo=2,Ctrl%NGeom
           Ctrl%Current=(/0,Ctrl%NSet,IGeo/)
           CALL SetGlobalCtrlIndecies(Ctrl)           
#ifdef MMech
           IF(HasMM()) CALL MM_ENERG(Ctrl) !!! temporary; overwrites energy terms calcd at prev geoms
#endif
           CALL OneSCF(Ctrl)
           CALL Forces(Ctrl)
        ENDDO
     ENDIF
#ifdef MMech
   EndIf
#endif
!
#ifdef MMech
   If(MMOnly()) Then
     Ctrl%Current=(/0,0,1/)
     CALL SetGlobalCtrlIndecies(Ctrl)           
     CALL MM_ENERG(Ctrl)
   ENDIF
#endif
  END SUBROUTINE CALC_GRAD_ONE_FORCE 
!--------------------------------------------------------
  SUBROUTINE CALC_GRAD_MD(Ctrl)
     IMPLICIT NONE
     TYPE(SCFControls)     :: Ctrl
     CALL MondoHalt(USUP_ERROR,' Look for MD in version 1.0b2. ')
  END SUBROUTINE CALC_GRAD_MD
!--------------------------------------------------------
  SUBROUTINE CALC_GRAD_QNEW_OPT(Ctrl)
     IMPLICIT NONE
     TYPE(SCFControls)     :: Ctrl
     INTEGER :: I,ISet,IGeo
!
        CALL SetGlobalCtrlIndecies(Ctrl)           
!
#ifdef MMech
     IF(MMOnly()) THEN
        Ctrl%Current=(/0,0,CGeo/)
        CALL SetGlobalCtrlIndecies(Ctrl)           
        CALL GeOp(Ctrl)
     ELSE
#endif
     DO ISet=1,Ctrl%NSet
!       Optimize geometry for each basis set
        Ctrl%Current=(/0,ISet,CGeo/)
        CALL SetGlobalCtrlIndecies(Ctrl)           
        CALL GeOp(Ctrl)
     ENDDO
#ifdef MMech
     ENDIF
#endif
  END SUBROUTINE CALC_GRAD_QNEW_OPT
!------------------------------------------------------------
  SUBROUTINE CALC_GRAD_QNEW_ONE_OPT(Ctrl)
     IMPLICIT NONE
     TYPE(SCFControls)     :: Ctrl
     INTEGER :: I,ISet,IGeo
!
!optimize only in last basis set and just one step
!
#ifdef MMech
   If(HasQM()) Then
#endif
!  Loop first over basis sets 
     DO ISet=1,Ctrl%NSet-1
        Ctrl%Current=(/0,ISet,1/)
        CALL SetGlobalCtrlIndecies(Ctrl)           
#ifdef MMech
        IF(ISet==1.AND.HasMM()) CALL MM_ENERG(Ctrl)
#endif
        CALL OneSCF(Ctrl)
     ENDDO
     Ctrl%Current=(/0,Ctrl%NSet,1/)
     CALL SetGlobalCtrlIndecies(Ctrl)           
     CALL GeOp(Ctrl)
!
!     IF(Ctrl%NGeom>1)THEN
!!  Go over additional geometries at last basis set
!        DO IGeo=2,Ctrl%NGeom
!           Ctrl%Current=(/0,Ctrl%NSet,IGeo/)
!           CALL SetGlobalCtrlIndecies(Ctrl)           
!           IF(HasMM()) CALL MM_ENERG(Ctrl) !!! temporary; overwrites energy terms calcd at prev geoms
!           CALL OneSCF(Ctrl)
!           Ctrl%Current=(/0,Ctrl%NSet,IGeo/)
!           CALL SetGlobalCtrlIndecies(Ctrl)           
!           CALL GeOp(Ctrl)
!        ENDDO
!     ENDIF
#ifdef MMech
   EndIf
#endif
!
#ifdef MMech
   If(MMOnly()) Then
     Ctrl%Current=(/0,0,1/)
     CALL SetGlobalCtrlIndecies(Ctrl)           
     CALL MM_ENERG(Ctrl)
     CALL GeOp(Ctrl)
   ENDIF
#endif
!    Optimize geometry only in last basis set
  END SUBROUTINE CALC_GRAD_QNEW_ONE_OPT
!--------------------------------------------------------
  SUBROUTINE CALC_TS_SEARCH(Ctrl)
     IMPLICIT NONE
     TYPE(SCFControls)     :: Ctrl
     INTEGER :: I,ISet,IGeo
!
     CALL MondoHalt(USUP_ERROR,' Look for transition state optimizer in version 1.0b2.')
  END SUBROUTINE CALC_TS_SEARCH
!--------------------------------------------------------
  SUBROUTINE CALC_GRAD_NO_GRAD(Ctrl)
     IMPLICIT NONE
     TYPE(SCFControls)     :: Ctrl
     INTEGER :: I,ISet,IGeo
!
#ifdef MMech
   If(HasQM()) Then
#endif
!  Loop first over basis sets 
     DO ISet=1,Ctrl%NSet
        Ctrl%Current=(/0,ISet,1/)
        CALL SetGlobalCtrlIndecies(Ctrl)           
#ifdef MMech
        IF(ISet==1.AND.HasMM()) CALL MM_ENERG(Ctrl)
#endif
        CALL OneSCF(Ctrl)
     ENDDO
     IF(Ctrl%NGeom>1)THEN
!  Go over additional geometries at last basis set
        DO IGeo=2,Ctrl%NGeom
           Ctrl%Current=(/0,Ctrl%NSet,IGeo/)
           CALL SetGlobalCtrlIndecies(Ctrl)           
#ifdef MMech
           IF(HasMM()) CALL MM_ENERG(Ctrl) !!! temporary; overwrites energy terms calcd at prev geoms, also for MM multiple geomeptries are not available yet
#endif
           CALL OneSCF(Ctrl)
        ENDDO
     ENDIF
#ifdef MMech
   EndIf
#endif
!
#ifdef MMech
   If(MMOnly()) Then
     Ctrl%Current=(/0,0,1/)
     CALL SetGlobalCtrlIndecies(Ctrl)           
     CALL MM_ENERG(Ctrl)
   ENDIF
#endif
  END SUBROUTINE CALC_GRAD_NO_GRAD
!-------------------------------------------------------------
  SUBROUTINE CALC_SinglePoint(Ctrl)
!
! Energy only calculation, at current geometry and basis set,
! as defined in Ctrl vector.
! Ctrl array is supposed to be set correctly at this point.
!
       IMPLICIT NONE
       TYPE(SCFControls)     :: Ctrl
       INTEGER :: I,ISet,IGeo
!
#ifdef MMech
       IF(HasQM()) THEN
!!! temporary; may overwrite MM energy terms calcd at prev geoms 
         IF(HasMM()) CALL MM_ENERG(Ctrl) 
         CALL OneSCF(Ctrl)
           IF(Ctrl%Grad/=GRAD_NO_GRAD) THEN
             CALL Forces(Ctrl)
           ENDIF
       ENDIF
#else
       CALL OneSCF(Ctrl)
         IF(Ctrl%Grad/=GRAD_NO_GRAD) THEN
           CALL Forces(Ctrl)
         ENDIF
#endif
!
#ifdef MMech
       IF(MMOnly()) Then
         CALL MM_ENERG(Ctrl)
       ENDIF
#endif
!
       END SUBROUTINE CALC_SinglePoint 
!-------------------------------------------------------------
  SUBROUTINE OptimizeNew(Ctrl)
     IMPLICIT NONE
     TYPE(SCFControls)     :: Ctrl
     INTEGER :: I,ISet,IGeo
!
! Geometry optimization, starting at CGeo.
! It is assumed, that CGeo is set by Ctrl at this point.
!
! Initialize file for geometry updates     
!
     IF(CBas==1.OR.MMOnly())THEN
       CALL OpenASCII(GeoFile,Geo,NewFile_O=.TRUE.)
       CLOSE(UNIT=Geo,STATUS='KEEP')
     ENDIF
!
#ifdef MMech
     IF(MMOnly()) THEN
         CALL GeOpNew(Ctrl)
     ELSE
#endif
!
!    Optimize geometry for each basis set
!
       DO ISet=1,Ctrl%NSet
         Ctrl%Current=(/0,ISet,CGeo/)
         CALL SetGlobalCtrlIndecies(Ctrl)
         CALL GeOpNew(Ctrl)
       ENDDO
#ifdef MMech
     ENDIF
#endif
  END SUBROUTINE OptimizeNew
!------------------------------------------------------------
    SUBROUTINE GeOpNew(Ctrl)
       TYPE(SCFControls)         :: Ctrl
       TYPE(CRDS)                :: GMLoc
       INTEGER                   :: JGeo,KS,BlkGeomSize,NCoordsEff,NCoordsEff0
       INTEGER                   :: Refresh,I,J,K,L,II,MaxGeOpSteps
       INTEGER                   :: NCart,NIntC
       REAL(DOUBLE)              :: RMSGrad,MaxGrad
       REAL(DOUBLE)              :: TrixThresh,AInvDistanceThresh
       REAL(DOUBLE)              :: BondConvCrit,AngleConvCrit,GradConvCrit
       REAL(DOUBLE)              :: MaxBondDispl,MaxAngleDispl,RMSIntDispl
       REAL(DOUBLE)              :: Etot
       REAL(DOUBLE),DIMENSION(3) :: PrevState,CrntState
       TYPE(INTC)                :: IntCs
       TYPE(DBL_VECT)            :: IntOld,Displ
       INTEGER                   :: GDIISMemory,FirstGeom
       INTEGER                   :: AccL
       LOGICAL                   :: GCnvrgd     
       LOGICAL                   :: LineSearchOn
       CHARACTER(LEN=DEFAULT_CHR_LEN) :: GMTag
!-------------------------------------------------------
!
       CALL SetGlobalCtrlIndecies(Ctrl)           
!
! Get some thresholds
!
         AInvDistanceThresh=1.D3
         GMTag=''
#ifdef MMech
       IF(HasQM()) THEN
         CALL Get(TrixThresh,'trixneglect',Tag_O=CurBase)
       ELSE
         GMTag='GM_MM'
         CALL Get(TrixThresh,'trixneglect',Tag_O=IntToChar(1))
       ENDIF
#else
         CALL Get(TrixThresh,'trixneglect',Tag_O=CurBase)
#endif
!
! Get GM data base
!
       CALL Get(GMLoc,TRIM(GMTag)//CurGeom)
       GMLoc%Confg=CGeo
!
! Set Convergence Criteria
!
       AccL=Ctrl%AccL(CBas)
       LineSearchOn=.FALSE.
       BondConvCrit=0.0005D0 ! in Angstroems
       AngleConvCrit=0.5D0 ! in degrees   
!
! Define maximum number of optimization steps
!
       NCart=3*GMLoc%Natms
       MaxGeOpSteps=MAX(NCart,200)
!
! Set GDIISmemory, it tells the number of 'remembered' steps in GDIIS
!
       GDIISMemory=500
!
       CALL OpenASCII(OutFile,Out)
       IF(Ctrl%Geop%CoordType==CoordType_Cartesian) THEN
         WRITE(*,*) '*****************************************'
         WRITE(*,*) 'Geometry Optimization in Cartesian Coords'
         WRITE(*,*) '*****************************************'
         WRITE(Out,*) '*****************************************'
         WRITE(Out,*) 'Geometry Optimization in Cartesian Coords'
         WRITE(Out,*) '*****************************************'
       ELSE IF(Ctrl%GeOp%Coordtype==CoordType_PrimInt) THEN
         WRITE(*,*) '****************************************'
         WRITE(*,*) 'Geometry Optimization in Internal Coords'
         WRITE(*,*) '****************************************'
         WRITE(Out,*) '****************************************'
         WRITE(Out,*) 'Geometry Optimization in Internal Coords'
         WRITE(Out,*) '****************************************'
       ENDIF
       CLOSE(Out)
!
! Define parameters of block operations for the optimizer,
! this will be the size of the vectors and matrices in the optimizer
!
!!!!       BlkGeomSize=3
!!!!       CALL Put(BlkGeomSize,'BlkGeomSize')
!
! Parametrize block matrix algebra for different cases
!
!!!!!  IF(Ctrl%Geop%CoordType==CoordType_Cartesian) THEN
!!!!!    NCoordsEff0=3*GMLoc%Natms
!!!!!  ELSE IF(Ctrl%GeOp%Coordtype==CoordType_PrimInt) THEN
!!!!!    NCoordsEff0=MAX(NIntC,3*GMLoc%Natms)
!!!!!  ENDIF
!!!!!
!!!!! General block parametrization, once NCoordsEff0 is set.
!!!!!
!!!!         CALL Put(NCoordsEff0,'NCoordsEff0') 
!!!!         NCoordsEff=(NCoordsEff0+BlkGeomSize-1)/BlkGeomSize 
!!!!         CALL Put(NCoordsEff,'NCoordsEff') 
!!!!         Natoms=NCoordsEff !!! number of blocks in a vector
!!!!         NBasF=BlkGeomSize*Natoms
!!!!         MaxAtms=Natoms+1
!!!!         CALL New(BSiz,Natoms)
!!!!         CALL New(OffS,Natoms)
!!!!         BSiz%I=BlkGeomSize
!!!!         OffS%I(1)=1
!!!!         DO I=2,Natoms
!!!!           OffS%I(I)=OffS%I(I-1)+BlkGeomSize
!!!!         ENDDO
!!!!         MaxBlkSize=MAXVAL(BSiz%I)
!!!!         PrintFlags%Mat=DEBUG_MATRICES
!
!
! Start iteration over geometries
!
       FirstGeom=Ctrl%Current(3)
       II=FirstGeom
300    CONTINUE
!
! Check for refresh of int. coord. definitions
!
       IF(II==FirstGeom) THEN !!!! refresh does not work in HDF
         CALL IntCReDef(Ctrl,Refresh)
         IF(Refresh/=0) THEN
           CALL GetIntCs(GMLoc%Carts%D,GMLoc%Natms,InfFile, &
                         IntCs,NIntC,Refresh)
         ENDIF
       ENDIF
       CALL INTCValue(IntCs,GMLoc%Carts%D)
       CALL New(IntOld,NIntC)
       IntOld%D=IntCs%Value
!      CALL PrtIntCoords(IntCs,IntCs%Value, &
!                        'Internals at step '//TRIM(IntToChar(II)))
!
! Get Energy and Cartesian gradient at current geometry
!
       CALL CALC_SinglePoint(Ctrl)
!
! Calculate single relaxation (SR) step from an inverse Hessian
!
       CALL NewDispl(Ctrl,Displ,NCart,NIntC)
       CALL SRStep(Ctrl,GMLoc,Displ,MaxGrad,RMSGrad,IntCs) 
!
! Calculate new geometry 
!
       CALL NewStructure(Ctrl,GMLoc,Displ,LineSearchOn,IntCs)
       CALL Delete(Displ)
!
! Check convergence
!
       CALL INTCValue(IntCs,GMLoc%Carts%D)
       IntOld%D=IntCs%Value-IntOld%D
       MaxBondDispl=Zero
       MaxAngleDispl=Zero
       DO I=1,NIntC 
          IF(IntCs%Def(I)(1:4)=='STRE') THEN 
            MaxBondDispl=MAX(MaxBondDispl,ABS(IntOld%D(I))) 
          ELSE
            MaxAngleDispl=MAX(MaxAngleDispl,ABS(IntOld%D(I))) 
          ENDIF
       ENDDO
       MaxBondDispl=MaxBondDispl/AngstromsToAU
       MaxAngleDispl=MaxAngleDispl*180.D0/PI
       RMSIntDispl=SQRT(DOT_PRODUCT(IntOld%D,IntOld%D)/DBLE(NIntC))
         CALL Delete(IntOld)
!
       GCnvrgd=RMSGrad<GTol(AccL).AND.MaxGrad<GTol(AccL)
!
! Review iterations
!
       CALL OpenASCII(OutFile,Out)
!
         WRITE(*,399) II-FirstGeom+1
       IF(LineSearchOn) THEN
         CALL Get(Etot,'ETot',StatsToChar(Ctrl%Current))
         WRITE(*,400) Etot
         WRITE(Out,400) Etot
       ELSE
         CALL Get(Etot,'ETot',StatsToChar(Ctrl%Previous))
         WRITE(*,401) Etot
         WRITE(Out,401) Etot
       ENDIF
!
       WRITE(*,410) MaxGrad
       WRITE(*,420) RMSGrad
       WRITE(*,430) MaxBondDispl
       WRITE(*,435) MaxAngleDispl
       WRITE(*,440) RMSIntDispl
       WRITE(Out,410) MaxGrad
       WRITE(Out,420) RMSGrad
       WRITE(Out,430) MaxBondDispl
       WRITE(Out,435) MaxAngleDispl
       WRITE(Out,440) RMSIntDispl
!
       CLOSE(Out,STATUS='KEEP')
!
399 FORMAT('GeOp step= ',I6)
400 FORMAT(' Total Energy at Current Geometry  = ',F20.8)
401 FORMAT('Total Energy at Previous Geometry  = ',F20.8)
410 FORMAT('           ','             MaxGrad= ',F12.6)
420 FORMAT('           ','             RMSGrad= ',F12.6)
430 FORMAT('           ','      Max Bond Displ= ',F12.6)
435 FORMAT('           ','     Max Angle Displ= ',F12.6)
440 FORMAT('           ','           RMS Displ= ',F12.6)
!
! Check convergence
!!!!!!!build in AccL dependent criteria for bonds and angles
!
       IF(II-FirstGeom+1<=MaxGeOpSteps) THEN
         IF(.NOT.GCnvrgd) THEN
!        IF(MaxBondDispl>BondConvCrit .OR. &
!           MaxAngleDispl>AngleConvCrit .OR. &
!           .NOT.GCnvrgd) THEN
           II=II+1
           Ctrl%Current(3)=II
           CALL SetGlobalCtrlIndecies(Ctrl)           
             CALL Put(GMLoc,TRIM(GMTag)//CurGeom)
#ifdef MMech
             IF(HasQM()) CALL Halt('put new GM coordinates separately, to be done')
#else
             CALL Halt('put new GM coordinates separately, to be done')
#endif
           GO TO 300
         ENDIF
       ELSE
       ENDIF
!
! Convergence is reached at this point, calculate final energy
! and finish optimization
!
       CALL CALC_SinglePoint(Ctrl)
       CALL Get(Etot,'ETot',StatsToChar(Ctrl%Current))
!
       CALL OpenAscii(OutFile,Out)
!
       WRITE(*,460) II
       WRITE(*,470) Etot
       WRITE(Out,460) II
       WRITE(Out,470) Etot
460  FORMAT('Geometry Optimization converged in ',I6,' steps.')
470  FORMAT('Total Energy at optimum structure= ',F20.8)
!
       CLOSE(Out,STATUS='KEEP')
!
! Tidy up
!
       CALL Delete(IntCs)
       CALL Delete(GMLoc)
!!!!!!      CALL Delete(BSiz)
!!!!!!       CALL Delete(OffS)
!
   END SUBROUTINE GeOpNew
!
!--------------------------------------------------------------------
!
      SUBROUTINE SRStep(Ctrl,GMLoc,Displ,MaxGrad,RMSGrad,IntCs)
!
! Single Relaxation step
!
      TYPE(SCFControls)              :: Ctrl 
      TYPE(CRDS)                     :: GMLoc
      TYPE(DBL_VECT)                 :: Displ
      REAL(DOUBLE)                   :: MaxGrad,RMSGrad
      TYPE(DBL_VECT)                 :: CartGrad,IntGrad,Grad
      TYPE(INTC)                     :: IntCs
      INTEGER                        :: I,J,NDim
      INTEGER                        :: NatmsLoc,NCart,NIntC
      LOGICAL                        :: DoInternals
!
      NatmsLoc=GMLoc%Natms
      NCart=3*NatmsLoc
      NDim =SIZE(Displ%D)
      DoInternals=.FALSE.
      IF(NCart/=NDim) DoInternals=.TRUE.
!
      CALL New(Grad,NDim)
!
! First, get Cartesian gradient
!
      CALL New(CartGrad,NCart)
      CALL Get(CartGrad,'GradE',Tag_O=CurGeom) 
!
! If requested, compute internal coord. gradients
!
      IF(DoInternals) THEN
        NIntC=SIZE(IntCs%Def)
        IF(NIntC/=NDim) CALL Halt('Dimension error in SRStep')
        CALL New(IntGrad,NDim)
        CALL CartToInternal(GMLoc%Carts%D,IntCs,CartGrad%D,IntGrad%D)
        Grad%D=IntGrad%D
      ELSE
        Grad%D=CartGrad%D
      ENDIF
!
! Steepest descent (may be either internal or Cartesian)
! or other inverse Hessian guess
!
      CALL SteepestDesc(Ctrl,Grad,Displ)
!
! Set constraints on the displacements
!
      CALL SetConstraint(IntCs,GMLoc%Carts%D,Displ)
!
! Check for convergence
!
      MaxGrad=Zero
      DO I=1,NDim ; MaxGrad=MAX(MaxGrad,ABS(Grad%D(I))) ; ENDDO
      RMSGrad=SQRT(DOT_PRODUCT(Grad%D,Grad%D)/DBLE(NDim))
!
! Tidy up
!
      IF(DoInternals) CALL Delete(IntGrad)
      CALL Delete(CartGrad)
      CALL Delete(Grad)
!
      END SUBROUTINE SRStep
!-------------------------------------------------------
      SUBROUTINE NewStructure(Ctrl,GMLoc,Displ,LineSearchOn,IntCs)
!
! Line search. Can be carried out either in internals
! or in Cartesian.
! Here we optimize a single parameter 'Fact',
! to get optimum for the energy at Fact*Displ displacement.
!
      TYPE(SCFControls)              :: Ctrl
      TYPE(CRDS)                     :: GMLoc
      TYPE(INTC)                     :: IntCs
      INTEGER                        :: I,J,II,MaxSteps,NDim,NInTc
      INTEGER                        :: NatmsLoc,NCart
      REAL(DOUBLE)                   :: EStart,Fact
      TYPE(DBL_VECT)                 :: Displ
      TYPE(DBL_RNK2)                 :: ActXYZ
      TYPE(DBL_VECT)                 :: ActDispl
      TYPE(DBL_VECT)                 :: Energy
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: GMTag
      LOGICAL                        :: IntSearch,LineSearchOn
      INTEGER                        :: OldGrad
!
! MaxSteps : Maximum number of Line Search Steps
! Displ    : Initial displacement vector from SR step, 
!            may be Cartesian or internal
! ActDispl : Size of actual displacement at the II-th line search step
!            maximum value is ActDispl=MaxFact*Displ
!            MaxFact must be small enough for internal coordinates
!            in order coordinate transformation can converge
!
      MaxSteps=15
!
      NatmsLoc=GMLoc%Natms
      NCart=3*NatmsLoc   
      NDim =SIZE(Displ%D)
      IntSearch=.FALSE.
      IF(NCart/=NDim) THEN
        IntSearch=.TRUE.
        NIntC =SIZE(IntCs%Def)
        IF(NIntC/=NDim) CALL Halt('Dimensionality error in LineSearch')
      ENDIF
!
! Geometry Tag
!
    GMTag=''
#ifdef MMech
    IF(HasMM()) GMTag='GM_MM'
#endif
!
! Allocate arrays
!
      CALL New(ActDispl,NDim)
      CALL New(ActXYZ,(/3,NatmsLoc/))
      CALL New(Energy,MaxSteps)
!
! Starting energy      
!
      CALL Get(EStart,'Etot',StatsToChar(Current))
!
! Increment Current and prepare for energy only calculation
!
      OldGrad=Ctrl%Grad
      Ctrl%Grad=GRAD_NO_GRAD
      Ctrl%Previous  =Ctrl%Current
      Ctrl%Current(3)=Ctrl%Current(3)+1
      CALL SetGlobalCtrlIndecies(Ctrl)
!
! Initialization
!
      Fact=One
      ActXYZ%D=GMLoc%Carts%D
!
      II=0
100   CONTINUE
      II=II+1
!
      ActDispl%D=Fact*Displ%D
!
! Construct new structure either by Cartesian or by internal
! displacements
!
      IF(IntSearch) THEN 
        CALL INTCValue(IntCs,ActXYZ%D)
        CALL InternalToCart(ActXYZ%D,IntCs,ActDispl%D)
      ELSE
        CALL CartRNK1ToCartRNK2(ActDispl%D,ActXYZ%D,.TRUE.)
      ENDIF
!
! Calculate energy at new structure; CurGeom has a new value now,
! and does not change until new call for LineSearch.
!
      GMLoc%Carts%D=ActXYZ%D
      CALL Put(GMLoc,TRIM(GMTag)//CurGeom)
      IF(.NOT.LineSearchOn) GO TO 200
! Continue here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL CALC_SinglePoint(Ctrl)
      CALL Get(Energy%D(II),'Etot',StatsToChar(Current))
!write(*,*) II,'  energy= ',Energy%D(II)
!
!
!
200   CONTINUE
!
! Set back Ctrl%Grad to its original value
!
      Ctrl%Grad=OldGrad
!
! Now call GDIIS
! and construct GDIIS step from SR step and previous displacements
! May work in either Cartesian or internal displacemets space.
! Subroutine modifies Cartesian Coordinates in GMLoc, now.
!
!      CALL GDIIS(CurGeom,GDIISMemory,GMLoc)
!
!
! Tidy up
!
      CALL Delete(ActDispl)
      CALL Delete(ActXYZ)
      CALL Delete(Energy)
!     IF(IntSearch) CALL Delete(ActIntC)
!
      END SUBROUTINE NewStructure
!
!-----------------------------------------------------------------
!
      SUBROUTINE IntCReDef(Ctrl,Refresh)
      TYPE(SCFControls) :: Ctrl
      INTEGER           :: Refresh
!
! WARNING! refresh may change the number of internal coordinates!
! This may cause problems in HDF!!!!!!
!
         IF(Ctrl%Current(3)==1) THEN
           Ctrl%GeOp%ReDefIntC=1
         ELSE
           Ctrl%GeOp%ReDefIntC=2
         ENDIF
!
         Refresh=Ctrl%GeOp%ReDefIntC
!
      END SUBROUTINE IntCReDef
!------------------------------------------------------------------
!
       SUBROUTINE NewDispl(Ctrl,Displ,NCart,NIntC)
       TYPE(SCFControls) :: Ctrl 
       TYPE(DBL_VECT) :: Displ
       INTEGER        :: NCart,NIntC
!
         IF(Ctrl%Geop%CoordType==CoordType_Cartesian) THEN
           CALL New(Displ,NCart)
           Displ%D(:)=Zero
         ELSE
           CALL New(Displ,NIntC)
           Displ%D(:)=Zero
         ENDIF
!
       END SUBROUTINE NewDispl
!
!-------------------------------------------------------
!
      SUBROUTINE SteepestDesc(Ctrl,Grad,Displ)
!
      TYPE(SCFControls) :: Ctrl
      TYPE(DBL_VECT)    :: Grad,Displ 
!
        IF(Ctrl%GeOp%CoordType==CoordType_Cartesian) THEN
          Displ%D=-1.D0*Grad%D !!!! it oscillates with 2.D0
        ELSE IF(Ctrl%GeOp%CoordType==CoordType_PrimInt) THEN
          Displ%D=-2.D0*Grad%D 
        ELSE
          Displ%D=-5.D0*Grad%D 
        ENDIF 
!
      END SUBROUTINE SteepestDesc
!-------------------------------------------------------
!
 END MODULE DrvFrcs
