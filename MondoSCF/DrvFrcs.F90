!------------------------------------------------------------------------------
!--  This code is part of the MondoSCF suite of programs for linear scaling 
!    electronic structure theory and ab initio molecular dynamics.
!
!--  Copyright (c) 2001, the Regents of the University of California.  
!    This SOFTWARE has been authored by an employee or employees of the 
!    University of California, operator of the Los Alamos National Laboratory 
!    under Contract No. W-7405-ENG-36 with the U.S. Department of Energy.  
!    The U.S. Government has rights to use, reproduce, and distribute this 
!    SOFTWARE.  The public may copy, distribute, prepare derivative works 
!    and publicly display this SOFTWARE without charge, provided that this 
!    Notice and any statement of authorship are reproduced on all copies.  
!    Neither the Government nor the University makes any warranty, express 
!    or implied, or assumes any liability or responsibility for the use of 
!    this SOFTWARE.  If SOFTWARE is modified to produce derivative works, 
!    such modified SOFTWARE should be clearly marked, so as not to confuse 
!    it with the version available from LANL.  The return of derivative works
!    to the primary author for integration and general release is encouraged. 
!    The first publication realized with the use of MondoSCF shall be
!    considered a joint work.  Publication of the results will appear
!    under the joint authorship of the researchers nominated by their
!    respective institutions. In future publications of work performed
!    with MondoSCF, the use of the software shall be properly acknowledged,
!    e.g. in the form "These calculations have been performed using MondoSCF, 
!    a suite of programs for linear scaling electronic structure theory and
!    ab initio molecular dynamics", and given appropriate citation.  
!------------------------------------------------------------------------------
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
  USE ChkSCFs
  USE ProcessControl    
  USE InOut
  USE AtomPairs
  USE LinALg
  IMPLICIT NONE 
  CONTAINS
!========================================================================================
!
!========================================================================================
    SUBROUTINE Forces(Ctrl)
      TYPE(SCFControls)  :: Ctrl
      INTEGER            :: ICyc,IGeo,IBas,AtA,A1
      TYPE(CRDS)         :: GM
      TYPE(DBL_VECT)     :: Frc
      CHARACTER(LEN=30)  :: ForceFile
!
      ICyc=Ctrl%Current(1)
      IBas=Ctrl%Current(2)
      IGeo=Ctrl%Current(3)
      CtrlVect=SetCtrlVect(Ctrl,'Force')
!
      CALL OpenHDF(InfFile)
      CALL Get(GM,Tag_O=IntToChar(IGeo))
      CALL CloseHDF()
!
!     Calculate the Analytical forces
!
      CALL Invoke('SForce',  CtrlVect)
      CALL Invoke('TForce',  CtrlVect)
      CALL Invoke('JForce',  CtrlVect)
      IF(HasDFT(Ctrl%Model(Ctrl%Current(2)))) &
           CALL Invoke('XCForce',     CtrlVect)       
      IF(HasHF(Ctrl%Model(Ctrl%Current(2))))THEN
!         CALL Invoke('XForce',     CtrlVect)       
         CALL MondoHalt(-999,'ONX Gradients not in yet, soon... ')
      ENDIF
!
!     Output the Forces
!      
      CALL New(Frc,3*NAtoms)
      CALL OpenHDF(InfFile)
      CALL Get(Frc    ,'GradE'    ,Tag_O=IntToChar(IGeo))
      CALL CloseHDF()
!
      ForceFile   = 'Force_'// TRIM(SCF_NAME)// "_" // TRIM(PROCESS_ID) // '.out'
      OPEN(UNIT=30,FILE=ForceFile,POSITION='APPEND',STATUS='UNKNOWN')
      WRITE(30,12) IGeo
      DO AtA = 1,NAtoms
         A1 = 3*(AtA-1)+1
         WRITE(30,11) AtA,GM%Carts%D(1,AtA),GM%Carts%D(2,AtA),GM%Carts%D(3,AtA),Frc%D(A1),Frc%D(A1+1),Frc%D(A1+2)
      ENDDO
      CLOSE(30)
      CALL Delete(Frc)
!
11    FORMAT('ATOM#',I6,2X,F20.16,2X,F20.16,2X,F20.16,2X,F20.16,2X,F20.16,2X,F20.16)
12    FORMAT('Geometry #',I6)
!
    END SUBROUTINE Forces
!========================================================================================
!
!========================================================================================
    SUBROUTINE NextMDGeometry(Ctrl)
      TYPE(SCFControls)     :: Ctrl
      INTEGER               :: ICyc,IBas,IGeo,AtA,A1
      TYPE(CRDS)            :: GM,GMOLD,GMNEW
      TYPE(DBL_VECT)        :: Frc
      REAL(DOUBLE)          :: DELT,DELT2,D1MASS,D2MASS,E_MD_KIN,E_TOT_EL,VSCALE
      CHARACTER(LEN=30)     :: MDFile
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
      E_MD_KIN = Zero
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
!
            E_MD_KIN = E_MD_KIN+Half*GM%AtMss%D(AtA)*(GM%Vects%D(1,AtA)**2+GM%Vects%D(2,AtA)**2+GM%Vects%D(3,AtA)**2)
!
#ifdef PERIODIC
            CALL AtomCyclic(GMNEW,GMNEW%Carts%D(:,AtA))
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
!
            E_MD_KIN = E_MD_KIN+Half*GM%AtMss%D(AtA)*(GM%Vects%D(1,AtA)**2+GM%Vects%D(2,AtA)**2+GM%Vects%D(3,AtA)**2)
!
#ifdef PERIODIC
            CALL AtomCyclic(GMNEW,GMNEW%Carts%D(:,AtA))
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
!
            E_MD_KIN = E_MD_KIN+Half*GM%AtMss%D(AtA)*(GM%Vects%D(1,AtA)**2+GM%Vects%D(2,AtA)**2+GM%Vects%D(3,AtA)**2)
!
#ifdef PERIODIC
            CALL AtomCyclic(GMNEW,GMNEW%Carts%D(:,AtA))
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
!     Compile Statisics and Store thr Energies
!
      CALL OpenHDF(InfFile)
      CALL Get(E_TOT_EL,'E_total',Tag_O=IntToChar(ICyc))
      CALL CloseHDF()
!
      MDFile   = 'MDEnergy_'// TRIM(SCF_NAME)// "_" // TRIM(PROCESS_ID) // '.out'
      OPEN(UNIT=30,FILE=MDFile,POSITION='APPEND',STATUS='UNKNOWN')
      WRITE(30,10) IGeo,E_MD_KIN,E_TOT_EL,E_MD_KIN+E_TOT_EL
      CLOSE(30)
10    FORMAT(2X,I6,2X,F22.16,2X,F22.16,2X,F22.16)

!
!      WRITE(*,*) 'E_MD_KIN = ',E_MD_KIN
!      WRITE(*,*) 'E_TOT_EL = ',E_TOT_EL
!      WRITE(*,*) 'E_TOTAL  = ',E_MD_KIN+E_TOT_EL
!
    END SUBROUTINE NextMDGeometry
!========================================================================================
!
!========================================================================================
    SUBROUTINE NextOGGeometry(Ctrl,IChk)
      TYPE(SCFControls)     :: Ctrl
      INTEGER               :: IChk
      INTEGER               :: ICyc,IBas,IGeo
!----------------------------------------------------------------------------------------
      ICyc=Ctrl%Current(1)
      IBas=Ctrl%Current(2)
      IGeo=Ctrl%Current(3)
      CtrlVect=SetCtrlVect(Ctrl,'QuNew')
      IF(IChk==1)THEN
!        Current geometry is good, compute gradients and update Hessian
         CALL Invoke('SForce',  CtrlVect)
         CALL Invoke('TForce',  CtrlVect)
         CALL Invoke('JForce',  CtrlVect)
         IF(HasDFT(Ctrl%Model(Ctrl%Current(2)))) &
              CALL Invoke('XCForce',     CtrlVect)       
         IF(HasHF(Ctrl%Model(Ctrl%Current(2))))THEN
!              CALL Invoke('XForce',     CtrlVect)       
            CALL MondoHalt(-999,'ONX Gradients not in yet, soon... ')
         ENDIF 
         CALL Invoke('BFGSHess',   CtrlVect)
      ENDIF
      CALL Invoke('NewStep',    CtrlVect)
    END SUBROUTINE NextOGGeometry
!========================================================================================
!
!========================================================================================
    SUBROUTINE Force_Num_1(Ctrl)
      TYPE(SCFControls)                :: Ctrl
      TYPE(CRDS)                       :: GM,GMOLD      
      TYPE(BCSR)                       :: S,P,F,T1,T2
      INTEGER                          :: ICyc,IGeo,IBas
      INTEGER                          :: AtA,IX,II,A1,A2,IS
      REAL(DOUBLE)                     :: DDelta
      REAL(DOUBLE),DIMENSION(2)        :: ES,ET,EJ,EN,EXC 
      REAL(DOUBLE),DIMENSION(NAtoms,3) :: F_S,F_T,F_J,F_N,F_XC,F_TOT
      TYPE(DBL_VECT)                   :: Frc_Num
!
      ICyc=Ctrl%Current(1)
      IBas=Ctrl%Current(2)
      IGeo=Ctrl%Current(3)
      CtrlVect=SetCtrlVect(Ctrl,'Force')
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
      CALL OpenHDF(InfFile)
      CALL Get(F,TrixFile('F',OffSet_O=0,Name_O=CtrlVect(1),Stats_O=Ctrl%Current(1:3)))
      CALL CloseHDF()
!
      DDelta = 1.D-6
      DO AtA = 1,Natoms
         DO IX = 1,3
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
               DO IS = 1,2
                  Ctrl%Current(1) = Ctrl%Current(1)+1
                  IF(IS .EQ. 1) CALL OneEMats(Ctrl)
!                  CALL SCFCycle(Ctrl)
                  IF(ConvergedQ(Ctrl))GOTO 999
999               CONTINUE 
               ENDDO
!======================================================================================
!              Calculate the Overlap: Tr[F.dP]
               CALL OpenHDF(InfFile)
               CALL Get(P,TrixFile('D',OffSet_O=0,Name_O=CtrlVect(1),Stats_O=Ctrl%Current(1:3)))
               CALL CloseHDF()  
               ES(II) = Trace(F,P)
!======================================================================================
!              GET ENERGIES
               CALL OpenHDF(InfFile)
               CALL Get(ET(II), 'Kin'     ,Tag_O=IntToChar(Ctrl%Current(1)))
               CALL Get(EJ(II), 'ene+eee' ,Tag_O=IntToChar(Ctrl%Current(1)))
               CALL Get(EN(II), 'enn+ene' ,Tag_O=IntToChar(Ctrl%Current(1)))
               CALL Get(EXC(II),'Exc'     ,Tag_O=IntToChar(Ctrl%Current(1)))
               CALL CloseHDF() 
            ENDDO
            F_S(AtA,IX)  = (ES(1)-ES(2))/(Two*DDelta)
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
      CALL OpenHDF(InfFile)
      CALL Put(Frc_Num,'GradE_Num',Tag_O=IntToChar(IGeo))
      CALL CloseHDF()
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
      REAL(DOUBLE),DIMENSION(2)        :: ES,ET,EJ,EN,EXC
      REAL(DOUBLE),DIMENSION(NAtoms,3) :: F_S,F_T,F_J,F_N,F_XC,F_TOT
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
               CALL Invoke('HiCu'   ,   CtrlVect)       
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
               CALL Get(EXC(II),'Exc'     ,Tag_O=IntToChar(Ctrl%Current(1)))
               CALL CloseHDF() 
            ENDDO
            F_S(AtA,IX)  = (ES(1)-ES(2))/(Two*DDelta)
            F_T(AtA,IX)  = (ET(1)-ET(2))/(Two*DDelta)
            F_J(AtA,IX)  = (EJ(1)-EJ(2))/(Two*DDelta)
            F_N(AtA,IX)  = (EN(1)-EN(2))/(Two*DDelta)
            F_XC(AtA,IX) = (EXC(1)-EXC(2))/(Two*DDelta)
            F_TOT(AtA,IX) = F_S(ATA,IX)+ F_T(AtA,IX)+F_J(AtA,IX)+F_N(AtA,IX)+F_XC(AtA,IX)
         ENDDO
      ENDDO
      WRITE(*,*) F_S(1,1)
      WRITE(*,*) F_T(1,1)
      WRITE(*,*) F_J(1,1)+F_N(1,1)
      WRITE(*,*) F_XC(1,1)     
      WRITE(*,*) ' '
      WRITE(*,*) F_TOT(1,1)
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
!
 END MODULE DrvFrcs
