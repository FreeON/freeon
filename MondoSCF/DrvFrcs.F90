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
      TYPE(CRDS)                     :: GM
      TYPE(DBL_VECT)                 :: Frc
      INTEGER                        :: Modl
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
!----------------------------------------------------------------------------------------
      Modl=Ctrl%Model(CBas)
!     Calculate analytic gradients piece by piece
      CALL LogSCF(Current,' Evaluating forces ')
      CtrlVect=SetCtrlVect(Ctrl,'ForceEvaluation')
!     The non-orthongonal response    
      CALL Invoke('SForce',CtrlVect)
!     Kinetic energy piece
      CALL Invoke('TForce',CtrlVect)
!     Build a density with last DM
      CALL Invoke('MakeRho',CtrlVect)
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
         CALL Get(Frc,'GradE',CurGeom)
         CALL Get(GM,Tag_O=CurGeom)
         CALL CloseHDF()
!        Print the gradients
         Mssg = '     Nuclear gradients for geometry #'//TRIM(CurGeom)//Rtrn
         CALL Print_Force(GM,Frc,Mssg,FileName_O,Unit_O,Fmat_O=2)
         Mssg = 'Nuclear gradients for geometry #'//TRIM(CurGeom)
         CALL Print_CheckSum_Force(Frc,Mssg,Unit_O=Unit_O)
!        Tidy
         CALL Delete(Frc)
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
!      Compute starting energy
       CALL OneSCF(Ctrl)
!      Print first configuration
       CALL OpenHDF(Ctrl%Info)
       CALL Get(GM,CurGeom)
       CALL Get(GM%ETotal,'Etot',Tag_O=StatsToChar(Ctrl%Current))
       GM%Confg=CGeo
       CALL PPrint(GM,GeoFile,Geo,'PDB')
       CALL Delete(GM)
       CALL CloseHDF()
!      Evaluate starting forces
       CALL Forces(Ctrl)
!      First Hessian update
       CALL Invoke('BFGSHess',CtrlVect)
       PrevState=Ctrl%Current
       Ctrl%Previous=Ctrl%Current
       Ctrl%Current=(/0,CBas,CGeo+1/)
       CrntState=Ctrl%Current
       CALL SetGlobalCtrlIndecies(Ctrl)           
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
             Ctrl%Previous=Ctrl%Current
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
         INTEGER            :: AccL
         REAL(DOUBLE)       :: StepSz,OldStep,E0,E1,ET,RelErrE,EUncert, &
                               RMSGrad,MAXGrad,RMSDisp,MaxDisp
         LOGICAL            :: ECnvrgd,GCnvrgd,XCnvrgd
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
!---------------------------------------------------------------------------
!        Open the InfFile
         CALL OpenHDF(Ctrl%Info)
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
            CALL Get(GM,CurGeom)
            GM%Confg=CGeo
            GM%ETotal=E1
            CALL PPrint(GM,GeoFile,Geo,'PDB')
         ENDIF
         CALL CloseHDF()
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
      CALL OpenHDF(InfFile)
      CALL Get(GM   ,Tag_O=CurGeom)
      CALL Get(GMOLD,Tag_O=CurGeom)
      CALL CloseHDF()
!     Load globals 
      CALL OpenHDF(InfFile)
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
      CALL CloseHDF()
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
                  CALL OpenHDF(InfFile)
                  CALL Put(GM,Tag_O=CurGeom)
                  CALL CloseHDF()
!                 Calculate the Overlap Term
                  CALL Invoke('MakeS'  ,   CtrlVect)
                  CALL OpenHDF(InfFile)
                  CALL Get(T1,TrixFile('S',Name_O=CtrlVect(1),Stats_O=Ctrl%Current(1:3)))
                  CALL CloseHDF()
                  CALL Multiply(P,F,T2)       
                  CALL Multiply(T2,P,T3)     
                  ES(II) = -Two*Trace(T3,T1)
!                 Calcualte the Kinectic Term
                  CALL Invoke('MakeT'  ,   CtrlVect)
                  CALL OpenHDF(InfFile)
                  CALL Get(T1,TrixFile('T',Name_O=CtrlVect(1),Stats_O=Ctrl%Current(1:3)))
                  CALL CloseHDF()
                  ET(II) =  Two*Trace(P,T1)        
!                 Calculate the Coulomb Term d(Vee + Vne)/dR
                  CALL Invoke('QCTC'   ,   CtrlVect)
                  CALL OpenHDF(InfFile)
                  CALL Get(T1,TrixFile('J',OffSet_O=0,Name_O=CtrlVect(1),Stats_O=Ctrl%Current(1:3)))
                  CALL CloseHDF()                  
                  EJ(II) =  Two*Trace(P,T1)
!                 Calculate the Coulomb Term d(Vne + Vnn)/dR
                  CALL OpenHDF(InfFile)
                  CALL Get(EN(II), 'E_NuclearTotal'    ,Tag_O=IntToChar(Ctrl%Current(1)))
                  CALL CloseHDF()
                  EN(II) =  Two*EN(II)
!                 Calculate the x term
                  IF(HasHF(ModelChem))  THEN
                     CALL Invoke('ONX'    ,   CtrlVect) 
                     CALL OpenHDF(InfFile)
                     CALL Get(T1,TrixFile('K',OffSet_O=0,Name_O=CtrlVect(1),Stats_O=Ctrl%Current(1:3)))
                     EX(II) = Trace(P,T1)   
                     CALL CloseHDF()
                  ENDIF
!                 Calculate the xc term
                  IF(HasDFT(ModelChem)) THEN
                     CALL Invoke('HiCu'   ,   CtrlVect)
                     CALL OpenHDF(InfFile)
                     CALL Get(T1,TrixFile('Kxc',OffSet_O=0,Name_O=CtrlVect(1),Stats_O=Ctrl%Current(1:3)))
                     CALL CloseHDF()
                     EXC(II) = Two*Trace(P,T1)
                  ENDIF
!                 Reset geometry
                  GM%Carts%D(IX,AtA) = GMOLD%Carts%D(IX,AtA)
                  CALL OpenHDF(InfFile)
                  CALL Put(GM,Tag_O=CurGeom)
                  CALL CloseHDF()
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
                  CALL OpenHDF(InfFile)
                  CALL Put(GM,Tag_O=CurGeom)
                  CALL CloseHDF()   
!                 SCF
                  DO N=1,4
                     Ctrl%Current(1)=ICyc+N
                     CALL SetGlobalCtrlIndecies(Ctrl)
                     IF(N==1) CALL OneEMats(Ctrl)
                     CALL SCFCycle(Ctrl) 
                  ENDDO
                  
!                 Get Energy     
                  CALL OpenHDF(InfFile)
                  CALL Get(ES(II),'Etot',StatsToChar(Ctrl%Current))
                  CALL CloseHDF()             
!                 Reset geometry
                  GM%Carts%D(IX,AtA) = GMOLD%Carts%D(IX,AtA)
                  CALL OpenHDF(InfFile)
                  CALL Put(GM,Tag_O=CurGeom)
                  CALL CloseHDF()
               ENDDO
               IA=3*(AtA-1)+IX
               WRITE(*,*) " Energys ",ES(1),ES(2),ES(1)-ES(2)
               FTot%D(IA)=(ES(1)-ES(2))/(Two*DDelta)
            ENDDO
         ENDDO
!        Do some checksumming and IO 
         CALL Print_Force(GMOLD,FTot,'   dSCF/dR')
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
#endif
!
 END MODULE DrvFrcs
