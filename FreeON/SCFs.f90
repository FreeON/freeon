!------------------------------------------------------------------------------
! This code is part of the MondoSCF suite of programs for linear scaling
! electronic structure theory and ab initio molecular dynamics.
!
! Copyright (2004). The Regents of the University of California. This
! material was produced under U.S. Government contract W-7405-ENG-36
! for Los Alamos National Laboratory, which is operated by the University
! of California for the U.S. Department of Energy. The U.S. Government has
! rights to use, reproduce, and distribute this software.  NEITHER THE
! GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
! OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the
! Free Software Foundation; either version 2 of the License, or (at your
! option) any later version. Accordingly, this program is distributed in
! the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
! the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
! PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
! While you may do as you like with this software, the GNU license requires
! that you clearly mark derivative software.  In addition, you are encouraged
! to return derivative works to the MondoSCF group for review, and possible
! disemination in future releases.
!------------------------------------------------------------------------------

#include "MondoConfig.h"

MODULE SCFs
  USE Parse
  USE InOut
  USE LinAlg
  USE GlobalObjects
  USE SCFKeys
  USE Overlay
  USE DynamicsKeys
  USE PunchHDF
  USE Numerics
  USE OptionKeys
  USE Functionals
  USE ControlStructures
  USE NEB
  USE SetXYZ
  USE PrettyPrint
  USE MondoLogger
  USE Utilities

  IMPLICIT NONE

  INTEGER HDFFileID,H5GroupID
  INTEGER,PARAMETER :: NOT_CONVERGE=0
  INTEGER,PARAMETER :: SCF_STALLED =1
  INTEGER,PARAMETER :: DIIS_NOPATH =2
  INTEGER,PARAMETER :: DID_CONVERGE=3
CONTAINS
  !===============================================================================
  !
  !===============================================================================
  SUBROUTINE SinglePoints(C)
    TYPE(Controls) :: C
    INTEGER        :: iBAS,iGEO,iBBegin

!!$    CALL MondoLog(DEBUG_MAXIMUM, "SinglePoints", "calculating energy of single geometry")
    ! Loop over geometry
    DO iGEO = 1,1 !C%Geos%NGeom
      ! Init previous state
      C%Stat%Previous%I=(/0,1,1/)
      ! Init iBBegin
      iBBegin = 1
      IF(iGEO > 1) iBBegin = C%Sets%NBSets
      ! Loop over basis sets
      DO iBAS=iBBegin,C%Sets%NBSets
        ! Archive
        CALL GeomArchive(iBAS,iGEO,C%Nams,C%Opts,C%Sets,C%Geos)
        CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
        ! Converge an SCF
        CALL SCF(iBAS,iGEO,C)
      ENDDO
    ENDDO
  END SUBROUTINE SinglePoints
  !===============================================================================
  !
  !===============================================================================
  SUBROUTINE SCF(cBAS,cGEO,C)
    TYPE(Controls)       :: C
    TYPE(DBL_RNK2),SAVE  :: ETot,DMax,DIIS
    INTEGER,PARAMETER    :: MaxSCFs = HAVE_MAX_SCF
    INTEGER              :: cBAS,cGEO,iSCF

    CALL MondoLog(DEBUG_NONE, "SCF", "doing SCF")

    ! Allocate space for action.
    CALL New(C%Stat%Action,1)

    ! Determine if there was a geomety or Basis Set Change
    CALL SameBasisSameGeom(cBAS,cGEO,C%Nams,C%Opts,C%Stat,C%Geos)

    ! Compute one-electron matrices
    CALL OneEMats(cBAS,cGEO,C%Nams,C%Sets,C%Stat,C%Opts,C%MPIs)

    ! Allocate space for convergence statistics
    CALL New(ETot,(/MaxSCFs,C%Geos%Clones/),(/0,1/))
    CALL New(DMax,(/MaxSCFs,C%Geos%Clones/),(/0,1/))
    CALL New(DIIS,(/MaxSCFs,C%Geos%Clones/),(/0,1/))

    ! For now we are forcing guess to superposition when we run the optimizer.
    ! This should get fixed in a proper way sometime.
    !
    ! DO NOT PUT INTO master YET!
    !IF(.FALSE.) THEN
    !  IF(C%Opts%Grad == GRAD_GO_DOWNHILL) THEN
    !    IF(C%Opts%Guess /= GUESS_EQ_SUPR) THEN
    !      CALL MondoLog(DEBUG_MAXIMUM, "SCF", "switching guess from "// &
    !        TRIM(IntToChar(C%Opts%Guess))// " to superposition")
    !      C%Opts%Guess = GUESS_EQ_SUPR
    !    ENDIF
    !  ENDIF
    !ENDIF
    ! End of fix.....

    DO iSCF=0,MaxSCFs
      ! Do an SCF cycle
      IF(SCFCycle(iSCF,cBAS,cGEO,C%Nams,C%Stat,C%Opts,C%Geos,C%Dyns,C%MPIs,ETot,DMax,DIIS)) THEN
        ! Free memory
        CALL Delete(ETot)
        CALL Delete(DMax)
        CALL Delete(DIIS)
        CALL Delete(C%Stat%Action)
        CALL MondoLog(DEBUG_NONE, "FreeON", "SCF converged to required accuracy")
        RETURN
      ENDIF
    ENDDO
    CALL MondoHalt(DRIV_ERROR,'Failed to converge SCF in '//TRIM(IntToChar(MaxSCFs))//' SCF iterations.')

  END SUBROUTINE SCF
  !===============================================================================
  !
  !===============================================================================
  SUBROUTINE SCFLogic(cSCF,cBAS,cGEO,SCF_STATUS,ODA_DONE,DIIS_FAIL,IConAls,N,S,O,D)
    TYPE(FileNames) :: N
    TYPE(State)     :: S
    TYPE(Options)   :: O
    TYPE(Dynamics)  :: D
    INTEGER         :: cSCF,cBAS,cGEO,IConAls,MinMDGeo,iREMOVE
    INTEGER         :: SCF_STATUS
    LOGICAL         :: DIIS_FAIL,ODA_DONE

!!$    CALL MondoLog(DEBUG_MAXIMUM, "SCFLogic", "cSCF = "//TRIM(IntToChar(cSCF))//", Action = "//TRIM(S%Action%C(1)))

    IF(cSCF == 0) THEN
      SCF_STATUS = NOT_CONVERGE
      ODA_DONE   = .FALSE.
      DIIS_FAIL  = .FALSE.
    ENDIF

    ! Logic for Algorithm Choice
    IF(O%ConAls(cBAS) == ODMIX_CONALS) THEN
      IF(ODA_DONE) THEN
        IConAls = DIIS_CONALS
      ELSE
        IF(SCF_STATUS == DIIS_NOPATH) THEN
          IConAls = DIIS_CONALS
          ODA_DONE = .TRUE.
          CALL MondoLog(DEBUG_NONE, "SCFLogic", "Turning on DIIS")
        ELSE
          IConAls = ODA_CONALS
        ENDIF
      ENDIF
    ELSE
      IConAls = O%ConAls(cBAS)
    ENDIF

    ! Parse for strict ODA or DIIS Over-Ride
    CALL OpenASCII(N%IFile,Inp)
    IF(OptKeyQ(Inp,CONALS_OVRIDE,CONALS_ODA))  IConAls = ODA_CONALS
    IF(OptKeyQ(Inp,CONALS_OVRIDE,CONALS_DIIS)) IConAls = DIIS_CONALS
    CLOSE(Inp)

    ! Defaults
    IF(cSCF < 1) THEN
      IConAls = NO_CONALS
    ENDIF

    IF(S%Action%C(1) == SCF_GUESSEQCORE .AND. cSCF < 1) THEN
      IConAls = NO_CONALS
    ENDIF

    IF(S%Action%C(1) == SCF_BASISSETSWITCH .AND. cSCF < 2) THEN
      IConAls = NO_CONALS
    ENDIF

    IF(S%Action%C(1) == SCF_RWBSS .AND. cSCF < 2) THEN
      IConAls = NO_CONALS
    ENDIF

    ! MD OverRule
    IF(D%DoingMD) THEN
      CALL CalculateMDGeo(D,iREMOVE,MinMDGeo)
      IF(cGEO > MinMDGeo) THEN
        ! IConAls = NO_CONALS
      ENDIF
    ENDIF

  END SUBROUTINE SCFLogic

  FUNCTION SCFCycle(cSCF,cBAS,cGEO,N,S,O,G,D,M,ETot,DMax,DIIS,CPSCF_O)
    TYPE(FileNames)    :: N
    TYPE(State)        :: S
    TYPE(Options)      :: O
    TYPE(Geometries)   :: G
    TYPE(Dynamics)     :: D
    TYPE(Parallel)     :: M
    TYPE(DBL_RNK2)     :: ETot,DMax,DIIS
    INTEGER            :: cSCF,cBAS,cGEO,iCLONE,Modl,IConAls
    INTEGER,SAVE       :: SCF_STATUS
    LOGICAL,OPTIONAL   :: CPSCF_O
    LOGICAL            :: SCFCycle,DoCPSCF
    LOGICAL,SAVE       :: DIIS_FAIL,ODA_DONE
    REAL(DOUBLE)       :: DIISErr

    ! Initialize
    SCFCycle=.FALSE.
    S%Current%I=(/cSCF,cBAS,cGEO/)

    ! Are we maybe solving CPSCF equations?
    IF(PRESENT(CPSCF_O))THEN
      DoCPSCF=CPSCF_O
    ELSE
      DoCPSCF=.FALSE.
    ENDIF

    ! Init and Archives the State
    CALL StateArchive(N,G,S,Init_O=.TRUE.)

    ! Decide on the Choice of convergence Algorithms
    CALL SCFLogic(cSCF,cBAS,cGEO,SCF_STATUS,ODA_DONE,DIIS_FAIL,IConAls,N,S,O,D)

    ! The options...
    IF(DoCPSCF)THEN
      CALL DensityLogic(cSCF,cBAS,cGEO,N,S,O,D,CPSCF_O=.TRUE.)
      CALL DensityBuild(N,S,M)
      IF(cSCF.EQ.0)S%Action%C(1)=CPSCF_START_RESPONSE
      IF(cSCF.GT.0)S%Action%C(1)=CPSCF_FOCK_BUILD
      IF(cSCF.GT.0)THEN
        CALL Invoke('QCTC',N,S,M)
        Modl=O%Models(cBAS)
        IF(HasHF(Modl)) CALL Invoke('ONX',N,S,M)
        IF(HasDFT(Modl)) THEN
          CALL Halt('SCFs: DFT-Response not yet supported!')
          CALL Invoke('HiCu',N,S,M)
        ENDIF
      ENDIF
      CALL Invoke('FBuild',N,S,M)
      IF(cSCF.GT.0) CALL Invoke('DDIIS',N,S,M)
      S%Action%C(1)=CPSCF_SOLVE_SCF
      CALL SolveSCF(cBAS,N,S,O,M)
      CALL Invoke('CPSCFStatus',N,S,M)
    ELSE
      CALL DensityLogic(cSCF,cBAS,cGEO,N,S,O,D)
      CALL DensityBuild(N,S,M)
      CALL FockBuild(cSCF,cBAS,N,S,O,M)

      ! Select the Case
      SELECT CASE (IConAls)

      CASE (DIIS_CONALS)
        CALL Invoke('DIIS',N,S,M)
        CALL SolveSCF(cBAS,N,S,O,M)
        CALL Invoke('SCFstats',N,S,M)

      CASE (ODA_CONALS)
        CALL SolveSCF(cBAS,N,S,O,M)
        CALL Invoke('ODA',N,S,M)
        IF(HasDFT(O%Models(cBAS)))THEN
          ! Rebuild non-linear KS matrix
          CALL DensityLogic(cSCF,cBAS,cGEO,N,S,O,D)
          CALL DensityBuild(N,S,M)
#ifdef DIPMW
          CALL Invoke('HiCu',N,S,M)
          CALL Invoke('DipMW',N,S,M)
          IF(.TRUE.) STOP
#else
          CALL Invoke('HiCu',N,S,M)
#endif
          CALL Invoke('FBuild',N,S,M)
        ENDIF
        CALL SolveSCF(cBAS,N,S,O,M)
        CALL Invoke('SCFstats',N,S,M)

      CASE (NO_CONALS)
        CALL SolveSCF(cBAS,N,S,O,M)
        CALL Invoke('SCFstats',N,S,M)

      CASE DEFAULT
        CALL MondoHalt(DRIV_ERROR,'Logic failure in SCFCycle')
      END SELECT
    ENDIF

    ! Archive and Check Status
    CALL StateArchive(N,G,S)
    SCF_STATUS=ConvergedQ(cSCF,cBAS,cGEO,N,S,O,G,D,M,ETot,DMax,DIIS,IConAls,CPSCF_O)
    S%Previous%I=S%Current%I
    IF(SCF_STATUS==DID_CONVERGE) SCFCycle=.TRUE.
  END FUNCTION SCFCycle

  FUNCTION ConvergedQ(cSCF,cBAS,cGEO,N,S,O,G,D,M,ETot,DMax,DIIS,IConAls,CPSCF_O)
    TYPE(FileNames)             :: N
    TYPE(State)                 :: S
    TYPE(Options)               :: O
    TYPE(Geometries)            :: G
    TYPE(Dynamics)              :: D
    TYPE(Parallel)              :: M
    TYPE(DBL_RNK2)              :: ETot,DMax,DIIS
    LOGICAL,OPTIONAL            :: CPSCF_O
    LOGICAL                     :: DoCPSCF,DoDIIS,DoODA,RebuildPostODA
    LOGICAL                     :: ALogic,BLogic,CLogic,DLogic,ELogic,A2Logic, &
         GLogic,QLogic,ILogic,OLogic,FLogic
    INTEGER                     :: cSCF,cBAS,cGEO,iGEO,iCLONE,MinMDGeo,iREMOVE
    REAL(DOUBLE)                :: DIISA,DIISB,DDIIS,DIISQ,       &
         DETOT,ETOTA,ETOTB,ETOTQ,ETest, &
         DDMAX,DMAXA,DMAXB,DMAXQ,DTest,ETOTO,ODAQ,DMaxMax
    INTEGER,DIMENSION(G%Clones) :: Converged
    INTEGER                     :: ConvergedQ,iSCF,IConAls,MinSCF,MaxSCF
    CHARACTER(LEN=DCL)          :: chGEO

    IF(PRESENT(CPSCF_O)) THEN
      DoCPSCF=CPSCF_O
    ELSE
      DoCPSCF=.FALSE.
    ENDIF

    ! Convergence thresholds
    CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", "entering ConvergedQ")

    IF(DoCPSCF) THEN
      ETest=RTol(O%AccuracyLevels(cBAS))
      DTest=DTol(O%AccuracyLevels(cBAS))
      IF(cSCF==0)THEN
        ConvergedQ=NOT_CONVERGE
        RETURN
      ENDIF

      ! Accumulate current statistics
      chGEO=IntToChar(iGEO)
      HDFFileID=OpenHDF(N%HFile)
      DO iCLONE=1,G%Clones
        HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
        ! Gather convergence parameters
        CALL Get(ETot%D(cSCF,iCLONE),'Prop'    )
        CALL Get(DMax%D(cSCF,iCLONE),'DPrimMax')
        CALL Get(DIIS%D(cSCF,iCLONE),'DDIISErr')
        CALL CloseHDFGroup(HDF_CurrentID)
        ! Load current energies
        G%Clone(iCLONE)%ETotal=ETot%D(cSCF,iCLONE)

        Converged(iCLONE)=NOT_CONVERGE
        IF(cSCF>1)THEN
          ETotA=ETot%D(cSCF-1,iCLONE)
          ETotB=ETot%D(cSCF  ,iCLONE)
          DMaxA=DMax%D(cSCF-1,iCLONE)
          DMaxB=DMax%D(cSCF  ,iCLONE)
          DIISA=DIIS%D(cSCF-1,iCLONE)
          DIISB=DIIS%D(cSCF  ,iCLONE)
          ! Absolute numbers
          dETot=ABS(ETotA-ETotB)
          dDMax=ABS(DMaxA-DMaxB)
          dDIIS=ABS(DIISA-DIISB)
          ! Relative numbers (Quotients)
          ETotQ=dETot/ABS(ETotB)
          !IF(CPSCF) write(*,*) 'dETot',dETot,' ETest',ETest
          !IF(CPSCF) write(*,*) 'DMaxB',DMaxB,' dTest',dTest
          DMaxQ=dDMax/ABS(DMaxB+1.D-50)
          DIISQ=dDIIS/ABS(DIISB+1.D-50)
          ! Convergence tests
          ! IF(((DMaxB<dTest.AND.ETotQ<ETest).OR.DMaxB<5D-1*dTest))THEN
          IF((DMaxB<dTest.AND.ETotQ<ETest).OR.DMaxB<dTest/10.0d0)THEN
            Converged(iCLONE)=DID_CONVERGE
            Mssg='Normal CPSCF convergence'
          ENDIF
          ! Look for stall out if we have at least one consecutive digit in the DM
          IF(DMaxB<1.D-1.AND.DMaxA<1.D-1)THEN
            ! Look for non-decreasing errors due to incomplete numerics
            IF(DIISQ<1.D-1.AND.DMaxQ<1.D-1.AND.cSCF>6)THEN
              IF(DIISB>DIISA.AND.DMaxB>DMaxA)THEN
                Mssg='CPSCF hit DDIIS & DMax increase.'
                Converged(iCLONE)=DID_CONVERGE
              ENDIF
            ELSEIF(DIISQ<1.D-2.AND.DMaxQ<1.D-2.AND.cSCF>6)THEN
              IF(DIISB>DIISA)THEN
                Mssg='CPSCF hit DDIIS increase'
                Converged(iCLONE)=DID_CONVERGE
              ELSEIF(DMaxQ<1D-1.AND.DMaxB>DMaxA)THEN
                Mssg='CPSCF hit DMax increase'
                Converged(iCLONE)=DID_CONVERGE
              ENDIF
            ELSEIF((DIISQ<1D-4.OR.DMaxQ<1D-4).AND.cSCF>6)THEN
              Mssg='CPSCF convergence due to DDIIS stagnation.'
              Converged(iCLONE)=DID_CONVERGE
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      CALL CloseHDF(HDFFileID)
      ! IF(cSCF>1)ConvergedQ=NOT_CONVERGE
      ConvergedQ = DID_CONVERGE
      DO iCLONE=1,G%Clones
        ConvergedQ=MIN(ConvergedQ,Converged(iCLONE))
      ENDDO
      ! Convergence announcement
      IF(ConvergedQ.NE.NOT_CONVERGE.AND.cSCF>2)THEN!.AND.PrintFlags%Key>DEBUG_MAXIMUM)THEN
        CALL MondoLogPlain(TRIM(Mssg))
        CALL MondoLogPlain("Normal CPSCF convergence")
      ENDIF

    ELSE ! IF(DoCPSCF) THEN

      ! NORMAL HUMANS CONVERGENCE CRITERIA

      IF(IConAls==DIIS_CONALS) THEN
        DoDIIS=.TRUE.
        DoODA =.FALSE.
      ELSEIF(IConAls==ODA_CONALS) THEN
        DoDIIS=.FALSE.
        DoODA =.TRUE.
      ELSE
        DoDIIS=.FALSE.
        DoODA =.FALSE.
      ENDIF

      ETest=ETol(O%AccuracyLevels(cBAS))
      DTest=DTol(O%AccuracyLevels(cBAS))

      ! Accumulate current statistics
      chGEO=IntToChar(iGEO)
      HDFFileID=OpenHDF(N%HFile)

      DO iCLONE=1,G%Clones

        HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
        ! Determine SCF if restarting MD
        MinSCF = O%MinSCF
        MaxSCF = O%MaxSCF

        ! Gather convergence parameters
        CALL Get(ETot%D(cSCF,iCLONE),'Etot')
        CALL Get(DMax%D(cSCF,iCLONE),'DMax')
        CALL Get(DIIS%D(cSCF,iCLONE),'DIISErr' )

        IF(DoODA.AND.cSCF>1)THEN
          CALL Get(ETotO,'ODAEnergy')
        ELSE
          ETotO=1D10
        ENDIF

        ! Load current energies
        G%Clone(iCLONE)%ETotal=ETot%D(cSCF,iCLONE)

        ! Load current energy into energy vector.
        G%Clone(iCLONE)%ETotalPerSCF%D(cSCF) = G%Clone(iCLONE)%ETotal

        Converged(iCLONE)=NOT_CONVERGE
        IF(cSCF>1)THEN
          ETotA=ETot%D(cSCF-1,iCLONE)
          ETotB=ETot%D(cSCF  ,iCLONE)
          DMaxA=DMax%D(cSCF-1,iCLONE)
          DMaxB=DMax%D(cSCF  ,iCLONE)
          DIISA=DIIS%D(cSCF-1,iCLONE)
          DIISB=DIIS%D(cSCF  ,iCLONE)
          IF(cSCF>2)THEN
            ETotA=1D10
            ! DMaxA=1D10
            DIISA=1D10
            DO iSCF=2,cSCF-1
              ETotA=MIN(ETotA,ETot%D(iSCF,iCLONE))
              ! DMaxA=MIN(DMaxA,DMax%D(iSCF,iCLONE))
              DIISA=MIN(DIISA,DIIS%D(iSCF,iCLONE))
            ENDDO
          ENDIF
          dETot=ABS(ETotA-ETotB)
          dDMax=ABS(DMaxA-DMaxB)
          dDIIS=ABS(DIISA-DIISB)

          ! Relative numbers
          ETotQ=dETot/ABS(ETotB)
          DMaxQ=dDMax/ABS(DMaxB+1.D-50)
          DIISQ=dDIIS/ABS(DIISB+1.D-50)
          IF(DoODA)THEN
            ODAQ=ABS(ETotB-ETotO)/ABS(ETotB)
          ELSE
            ODAQ=Zero
          ENDIF

    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'ODAQ  = '//TRIM(FltToChar(ODAQ)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'ETotQ = '//TRIM(FltToChar(ETotQ)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'DIISQ = '//TRIM(FltToChar(DIISQ)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'DMaxQ = '//TRIM(FltToChar(DMaxQ)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'ETOTO = '//TRIM(FltToChar(ETotO)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'ETOTA = '//TRIM(FltToChar(ETOTA)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'ETOTB = '//TRIM(FltToChar(ETOTB)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'DIISA = '//TRIM(FltToChar(DIISA)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'DIISB = '//TRIM(FltToChar(DIISB)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'DMaxA = '//TRIM(FltToChar(DMaxA)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'DMaxB = '//TRIM(FltToChar(DMaxB)))

          Converged(iCLONE)=NOT_CONVERGE
          ! Convergence from above +/- expected delta relative to historical
          ! minimum energy
          ALogic=ETotB*(One+ETest)<ETotA
          ! Convergence from above +/- expected delta simply relative to
          ! previous energy
          A2Logic=ETot%D(cSCF  ,iCLONE)*(One+ETest)<ETot%D(cSCF-1,iCLONE)
          ! Met all criteria
          CLogic=DMaxB<DTest.AND.ETotQ<ETest.AND.DMaxB.NE.Zero
          ! Exceeded density criteria
          DLogic=DMaxB<5D-2*DTest.AND.DMaxB.NE.Zero
          ! Exceeded energy criteria
          ELogic=ETotQ<3D-2*ETest.AND.DMaxB<1D-2
          ! Quasi convergence from below (bad)
          QLogic=(.NOT.ALogic).AND.DLogic.AND.ELogic
          ! Going to wrong state with DIIS
          ILogic=DoDIIS.AND.DLogic.AND.(.NOT.ELogic)
          ! DIIS is oscillating
          OLogic=DoDIIS.AND.(.NOT.ALogic).AND.( (ETotQ>1D-4.AND.(DMaxQ>1D0.AND.DIISQ>1D0)).OR. &
               (ETotQ>2D-3.AND.(DMaxQ>1D0.OR.DIISQ>1D0)) )
          ! Maybe DIIS would be a good idea
          GLogic=DoODA.AND.cSCF>4.AND.DIISB<75D-5.AND.ETotQ<5D-5.AND.DMaxB<1D-1
          ! If we are increasing with ODA and rebuild is on, we are well fucked.
          FLogic=DoODA.AND..NOT.ALogic.AND.cSCF>3

          ! Sort through logic hopefully in the conditionally correct order ...
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'ETest  = '//TRIM(FltToChar(ETest)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'DTest  = '//TRIM(FltToChar(DTest)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'ALogic = '//TRIM(LogicalToChar(ALogic)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'A2Logic= '//TRIM(LogicalToChar(A2Logic)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'ELogic = '//TRIM(LogicalToChar(ELogic)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'CLogic = '//TRIM(LogicalToChar(CLogic)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'DLogic = '//TRIM(LogicalToChar(DLogic)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'QLogic = '//TRIM(LogicalToChar(QLogic)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'ILogic = '//TRIM(LogicalToChar(ILogic)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'GLogic = '//TRIM(LogicalToChar(GLogic)))
    !!$      CALL MondoLog(DEBUG_MAXIMUM, "ConvergedQ", 'FLogic = '//TRIM(LogicalToChar(FLogic)))

          ! No message.
          Mssg=" "
          IF(ALogic.AND.CLogic)THEN
            Converged(iCLONE)=DID_CONVERGE
            Mssg='Normal SCF convergence.'
          ELSEIF(A2Logic.AND.DLogic)THEN
            Converged(iCLONE)=DID_CONVERGE
            Mssg='Convergence of density only'
          ELSEIF(ALogic.AND.ELogic)THEN
            Converged(iCLONE)=DID_CONVERGE
            Mssg='Convergence of energy only'
          ELSEIF(QLogic)THEN
            Converged(iCLONE)=DID_CONVERGE
            Mssg='Quasi convergence from wrong side.'
          ELSEIF(FLogic)THEN
            ! Converged(iCLONE)=DID_CONVERGE
            Mssg='Warning: ODA not strictly decreasing'
          ELSEIF(OLogic)THEN
            Converged(iCLONE)=SCF_STALLED
            Mssg='DIIS oscillation'
          ELSEIF(GLogic)THEN
            Converged(iCLONE)=DIIS_NOPATH
          ENDIF
        ELSE
          DIISA = 0.0D0
          DIISB = 0.0D0
        ENDIF

        IF(DoDIIS.AND.DIISB<DIISA)THEN
          ! If DIIS is making progress, then turn on archivation of the density
          CALL Put(.TRUE.,'ArchiveDensity')
        ELSE
          ! otherwise, dont archive a potential instability
          CALL Put(.FALSE.,'ArchiveDensity')
        ENDIF
        CALL CloseHDFGroup(HDF_CurrentID)
      ENDDO
      CALL CloseHDF(HDFFileID)

      ConvergedQ=DID_CONVERGE
      DO iCLONE=1,G%Clones
        ConvergedQ=MIN(ConvergedQ,Converged(iCLONE))
      ENDDO

      ! Molecular Dynamics Convergence Criteria
      !IF(D%DoingMD) THEN
      !  CALL CalculateMDGeo(D,iREMOVE,MinMDGeo)

      !  DMaxMax=Zero
      !  DO iCLONE=1,G%Clones
      !    DMaxMax = MAX(DMax%D(cSCF,iCLONE),DMaxMax)
      !  ENDDO

      !  IF(cGEO > MinMDGeo .AND. cSCF .GE. MinSCF) THEN
      !    ConvergedQ=DID_CONVERGE
      !    Mssg = "MD Verlet SCF convergence"
      !    CALL MondoLogPlain(TRIM(Mssg))
      !    RETURN
      !  ELSE
      !    IF(DMaxMax > DTest*1.D-2) THEN
      !      ConvergedQ=NOT_CONVERGE
      !      Mssg = " "
      !    ELSE
      !      ConvergedQ=DID_CONVERGE
      !      Mssg = "MD SCF convergence"
      !    ENDIF
      !  ENDIF
      !ENDIF

      ! No message.
      Mssg = " "
      IF(cSCF .LT. MinSCF) THEN
        ConvergedQ=NOT_CONVERGE
        Mssg = " "
      ENDIF
      IF(cSCF .GE. MaxSCF) THEN
        ConvergedQ=DID_CONVERGE
        Mssg = "Forced SCF convergence"
      ENDIF

      ! Convergence announcement
      IF(Mssg .NE. " " .AND. cSCF >0)THEN
        CALL MondoLogPlain(TRIM(Mssg))
      ENDIF
    ENDIF

  END FUNCTION ConvergedQ

  !===============================================================================
  ! BUILD A HGTF DENSITY BY HOOK OR BY CROOK
  !===============================================================================
  SUBROUTINE DensityBuild(N,S,M)
    TYPE(FileNames)     :: N
    TYPE(State)         :: S
    TYPE(Parallel)      :: M
    INTEGER             :: oldSCF

    CALL MondoLog(DEBUG_MAXIMUM, "DensityBuild", "Action = "//TRIM(S%Action%C(1)))
    IF(TRIM(S%Action%C(1))/=SCF_DENSITY_NORMAL   .AND. &
       TRIM(S%Action%C(1))/=SCF_BASISSETSWITCH   .AND. &
       TRIM(S%Action%C(1))/=CPSCF_START_RESPONSE .AND. &
       TRIM(S%Action%C(1))/=CPSCF_DENSITY_NORMAL) THEN

      CALL Invoke('P2Use',N,S,M)
    ENDIF

    ! Build some density ...
    IF(S%Action%C(1)/=CPSCF_START_RESPONSE) THEN
      CALL Invoke('MakeRho',N,S,M)
    ENDIF

  END SUBROUTINE DensityBuild

  SUBROUTINE DensityLogic(cSCF,cBAS,cGEO,N,S,O,D,CPSCF_O)
    TYPE(FileNames)    :: N
    TYPE(State)        :: S
    TYPE(Options)      :: O
    TYPE(Dynamics)     :: D
    INTEGER            :: cSCF,cBAS,cGEO,pBAS,I,J
    LOGICAL            :: DoCPSCF
    LOGICAL,OPTIONAL   :: CPSCF_O

    pBAS=S%Previous%I(2)

    ! S%Current%I=(/cSCF,cBAS,cGEO/)
    IF(PRESENT(CPSCF_O))THEN
      DoCPSCF=CPSCF_O
    ELSE
      DoCPSCF=.FALSE.
    ENDIF

    ! Determine the Action to be Taken
    IF(DoCPSCF)THEN
      IF(O%Guess==GUESS_EQ_DIPOLE.AND.cSCF==0)THEN
        O%Guess=0
        S%Previous%I=S%Current%I
        S%Action%C(1)=CPSCF_START_RESPONSE
      ELSE
        S%Action%C(1)=CPSCF_DENSITY_NORMAL
      ENDIF
    ELSE
      IF(O%Guess==GUESS_EQ_CORE)THEN
        O%Guess=0
        S%Previous%I  = S%Current%I
        S%Action%C(1) = SCF_GUESSEQCORE
      ELSEIF(O%Guess==GUESS_EQ_SUPR)THEN
        O%Guess=0
        S%Previous%I  = S%Current%I
        S%Action%C(1) = SCF_SUPERPOSITION
      ELSEIF(O%Guess==GUESS_EQ_NUGUESS)THEN
        O%Guess=0
        S%Previous%I  = O%RestartState%I
        S%Action%C(1) = SCF_SUPERPOSITION
      ELSEIF(O%Guess==GUESS_EQ_NEWGEOM)THEN
        O%Guess=0
        S%Previous%I  = O%RestartState%I
        S%Action%C(1) = SCF_EXTRAPOLATE
      ELSEIF(O%Guess==GUESS_EQ_RESTART)THEN
        IF(S%SameBasis .AND. .NOT. S%SameGeom) THEN
          O%Guess=0
          S%Previous%I  = O%RestartState%I
          S%Action%C(1) = SCF_EXTRAPOLATE
        ELSEIF( .NOT. S%SameBasis) THEN
          O%Guess=0
          S%Previous%I  = O%RestartState%I
          S%Action%C(1) = SCF_RWBSS
        ELSE
          O%Guess=0
          S%Previous%I  = O%RestartState%I
          S%Action%C(1) = SCF_RESTART
        ENDIF
      ELSEIF(S%SameBasis .AND. .NOT.S%SameGeom)THEN
        O%Guess=0
        S%Action%C(1) = SCF_EXTRAPOLATE
      ELSEIF(.NOT. S%SameBasis .OR. pBAS /= cBAS)THEN
        O%Guess=0
        S%Action%C(1)   = SCF_BASISSETSWITCH
        S%Previous%I(1) = S%Previous%I(1)+1
      ELSE
        O%Guess=0
        S%Action%C(1) = SCF_DENSITY_NORMAL
      ENDIF
    ENDIF

    ! If we are doing MD, Geuss to P2Use is Different
    IF(D%DoingMD .AND. cSCF==0) THEN
      S%Action%C(1)=O%GeussToP2Use
    ENDIF

    CALL MondoLog(DEBUG_MAXIMUM, "DensityLogic", "Action = "//TRIM(S%Action%C(1)))

    ! Reset
    S%SameBasis=.TRUE.
    S%SameGeom =.TRUE.
    S%SameCrds =.TRUE.
    S%SameLatt =.TRUE.

  END SUBROUTINE DensityLogic

  !===============================================================================
  ! BUILD A FOCK MATRIX
  !===============================================================================
  SUBROUTINE FockBuild(cSCF,cBAS,N,S,O,M)
    TYPE(FileNames):: N
    TYPE(State)    :: S
    TYPE(Options)  :: O
    TYPE(Parallel) :: M
    REAL(DOUBLE)   :: Lambda
    INTEGER        :: cSCF,cBAS,Modl
    LOGICAL        :: DoDIIs
    !----------------------------------------------------------------------------!
    Modl=O%Models(cBAS)

    IF(S%Action%C(1).NE.CPSCF_START_RESPONSE) CALL Invoke('QCTC',N,S,M)

    IF(S%Action%C(1)/=SCF_GUESSEQCORE.AND.S%Action%C(1).NE.CPSCF_START_RESPONSE)THEN
      IF(HasHF(Modl)) CALL Invoke('ONX',N,S,M)
#ifdef DIPMW
      IF(HasDFT(Modl)) THEN
        CALL Invoke('HiCu',N,S,M)
        CALL Invoke('DipMW',N,S,M)
        IF(.TRUE.) STOP
      ENDIF
#else
      IF(HasDFT(Modl)) THEN
        CALL Invoke('HiCu',N,S,M)
      ENDIF
#endif
    ENDIF

    CALL Invoke('FBuild',N,S,M)
  END SUBROUTINE FockBuild
  !===============================================================================
  ! EXTRAPOLATE (DIIS)
  !===============================================================================
  SUBROUTINE Xtra(cSCF,cBAS,N,S,O,M)
    TYPE(FileNames):: N
    TYPE(State)    :: S
    TYPE(Options)  :: O
    TYPE(Parallel) :: M
    REAL(DOUBLE)   :: Lambda
    INTEGER        :: cSCF,cBAS,Modl
    LOGICAL        :: DoDIIs
    !----------------------------------------------------------------------------!
    DoDIIS=cSCF>0
    IF(DoDIIS)THEN
      SELECT CASE(S%Action%C(1))
      CASE(CPSCF_SOLVE_SCF,CPSCF_START_RESPONSE,CPSCF_DENSITY_NORMAL,CPSCF_FOCK_BUILD)
        CALL Invoke('DDIIS',N,S,M)
      CASE DEFAULT
        CALL Invoke('DIIS',N,S,M)
      END SELECT
    ENDIF
  END SUBROUTINE Xtra
  !===============================================================================
  ! Solve the SCF equations
  !===============================================================================
  SUBROUTINE SolveSCF(cBAS,N,S,O,M)
    TYPE(FileNames):: N
    TYPE(State)    :: S
    TYPE(Options)  :: O
    TYPE(Parallel) :: M
    INTEGER        :: cBAS
    !----------------------------------------------------------------------------!

    IF(S%Action%C(1)==CPSCF_SOLVE_SCF.OR. &
       S%Action%C(1)==CPSCF_START_RESPONSE)THEN
      CALL Invoke('TC2Response',N,S,M)
    ELSEIF(O%Methods(cBAS)==RH_R_SCF)THEN
      CALL Invoke('RHeqs',N,S,M)
    ELSEIF(O%Methods(cBAS)==SDMM_R_SCF) THEN
      CALL Invoke('SDMM',N,S,M)
    ELSEIF(O%Methods(cBAS)==PM_R_SCF) THEN
      CALL Invoke('PM',N,S,M)
    ELSEIF(O%Methods(cBAS)==SP2_R_SCF) THEN
      CALL Invoke('SP2',N,S,M)
    ELSEIF(O%Methods(cBAS)==SP4_R_SCF) THEN
      CALL Invoke('SP4',N,S,M)
    ELSEIF(O%Methods(cBAS)==TS4_R_SCF) THEN
      CALL Invoke('TS4',N,S,M)
    ELSE
      CALL MondoHalt(99,'Unknown method key = '//TRIM(IntToChar(O%Methods(cBAS))))
    ENDIF
  END SUBROUTINE SolveSCF
  !---------------------------------------------------------------------------------
  !
  !---------------------------------------------------------------------------------
  SUBROUTINE SameBasisSameGeom(cBAS,cGEO,N,O,S,G)
    TYPE(FileNames)    :: N
    TYPE(Options)      :: O
    TYPE(State)        :: S
    TYPE(Geometries)   :: G
    TYPE(BSET),SAVE    :: BS,BS_rs
    TYPE(CRDS),SAVE    :: GM,GM_rs
    REAL(DOUBLE)       :: MaxDiff
    CHARACTER(LEN=DCL) :: chBAS,chGEO
    INTEGER            :: I,J,cBAS,cGEO,pBAS,pGEO,iCLONE
    LOGICAL,DIMENSION(G%Clones) :: SameCrds,SameLatt

    pBAS=S%Previous%I(2)
    pGEO=S%Previous%I(3)
    S%Current%I=(/0,cBAS,cGEO/)

    SameCrds=.TRUE.
    SameLatt=.TRUE.
    S%SameBasis = .TRUE.

    IF(O%Guess==GUESS_EQ_RESTART.OR.O%Guess==GUESS_EQ_NUGUESS)THEN
      DO iCLONE=1,G%Clones
        HDFFileID=OpenHDF(N%HFile)
        HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
        chBAS = IntToChar(S%Current%I(2))
        chGEO = IntToChar(S%Current%I(3))

        IF(iCLONE==1)CALL Get(BS,Tag_O=chBAS)

        CALL Get(GM,Tag_O=chGEO)
        CALL CloseHDFGroup(HDF_CurrentID)
        CALL CloseHDF(HDFFileID)

        HDFFileID=OpenHDF(N%RFile)
        HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
        chBAS = IntToChar(O%RestartState%I(2))
        chGEO = IntToChar(O%RestartState%I(3))
        IF(iCLONE==1)CALL Get(BS_rs,Tag_O=chBAS)
        CALL Get(GM_rs,Tag_O=chGEO)
        CALL CloseHDFGroup(HDF_CurrentID)
        CALL CloseHDF(HDFFileID)
        IF(iCLONE==1.AND.BS%BName/=BS_rs%BName)S%SameBasis=.FALSE.
        MaxDiff=Zero
        DO I=1,GM%Natms
          MaxDiff=MAX(MaxDiff,ABS(GM%Carts%D(1,I)-GM_rs%Carts%D(1,I)) + &
               ABS(GM%Carts%D(2,I)-GM_rs%Carts%D(2,I)) + &
               ABS(GM%Carts%D(3,I)-GM_rs%Carts%D(3,I)))
        ENDDO
        IF(MaxDiff>1D-8)SameCrds(iCLONE)=.FALSE.
        MaxDiff=Zero
        DO I=1,3
          DO J=1,3
            MaxDiff = MAX(MaxDiff,ABS(GM%PBC%BoxShape%D(I,J)-GM_rs%PBC%BoxShape%D(I,J)))
          ENDDO
        ENDDO
        IF(MaxDiff>1D-8)SameLatt(iCLONE)=.FALSE.
      ENDDO
    ELSE
      DO iCLONE=1,G%Clones

        HDFFileID=OpenHDF(N%HFile)
        HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
        chBAS = IntToChar(cBAS)
        chGEO = IntToChar(cGEO)

        CALL Get(GM,Tag_O=chGEO)

        IF(iCLONE==1) CALL Get(BS,Tag_O=chBAS)

        chBAS = IntToChar(pBAS)
        chGEO = IntToChar(pGEO)

        IF(iCLONE==1) CALL Get(BS_rs,Tag_O=chBAS)
        CALL Get(GM_rs,Tag_O=chGEO)
        CALL CloseHDFGroup(HDF_CurrentID)
        CALL CloseHDF(HDFFileID)

        IF(iCLONE==1.AND.BS%BName/=BS_rs%BName)S%SameBasis=.FALSE.
        MaxDiff=Zero
        DO I=1,GM%Natms
          MaxDiff=MAX(MaxDiff,ABS(GM%Carts%D(1,I)-GM_rs%Carts%D(1,I)) + &
               ABS(GM%Carts%D(2,I)-GM_rs%Carts%D(2,I)) + &
               ABS(GM%Carts%D(3,I)-GM_rs%Carts%D(3,I)))
        ENDDO
        IF(MaxDiff>1D-8)SameCrds(iCLONE)=.FALSE.
        MaxDiff=Zero
        DO I=1,3
          DO J=1,3
            MaxDiff = MAX(MaxDiff,ABS(GM%PBC%BoxShape%D(I,J)-GM_rs%PBC%BoxShape%D(I,J)))
          ENDDO
        ENDDO
        IF(MaxDiff>1D-8)SameLatt(iCLONE)=.FALSE.
      ENDDO
    ENDIF
    S%SameGeom=.TRUE.
    S%SameLatt=.TRUE.
    DO iCLONE=1,G%Clones
      S%SameGeom=S%SameGeom.AND.SameCrds(iCLONE).AND.SameLatt(iCLONE)
      S%SameLatt=S%SameLatt.AND.SameLatt(iCLONE)
    ENDDO

    ! A call to get() will allocate. We therefore need to delete these
    ! variables.
    CALL DELETE(GM)
    CALL DELETE(GM_rs)
    CALL DELETE(BS)
    CALL DELETE(BS_rs)

  END SUBROUTINE SameBasisSameGeom
  !===============================================================================
  !
  !===============================================================================
  SUBROUTINE OneEMats(cBAS,cGEO,N,B,S,O,M)
    TYPE(FileNames):: N
    TYPE(BasisSets):: B
    TYPE(State)    :: S
    TYPE(Options)  :: O
    TYPE(Parallel) :: M
    INTEGER        :: cBAS,cGEO,pBAS
    LOGICAL        :: DoPFFT
    !----------------------------------------------------------------------------!
    !
    S%Action%C(1)='OneElectronMatrices'
    !
    DoPFFT = .FALSE.
    IF(O%Guess==GUESS_EQ_CORE)     DoPFFT=.TRUE.
    IF(O%Guess==GUESS_EQ_SUPR)     DoPFFT=.TRUE.
    IF(O%Guess==GUESS_EQ_RESTART)  DoPFFT=.TRUE.
    IF(O%Guess==GUESS_EQ_NEWGEOM)  DoPFFT=.TRUE.
    IF(O%Guess==GUESS_EQ_NUGUESS)  DoPFFT=.TRUE.
    IF(.NOT. S%SameLatt)           DoPFFT=.TRUE.
    IF(.NOT. S%SameBasis)          DoPFFT=.TRUE.
    !
    IF(DoPFFT) THEN
      CALL Invoke('MakePFFT',N,S,M)
    ENDIF
    !
    IF((O%Guess==GUESS_EQ_RESTART .AND.(.NOT.S%SameGeom)) &
         .OR.O%Guess==GUESS_EQ_NEWGEOM)THEN
      ! Make previous geometrys S matrix
      S%Action%C(1)='RestartGeomSwitch'
      CALL Invoke('MakeS',N,S,M)
      ! now make current geometrys S matrix
      S%Action%C(1)='OneElectronMatrices'
      CALL Invoke('MakeS',N,S,M)
    ELSE
      CALL Invoke('MakeS',N,S,M)
    ENDIF
    IF(O%Methods(cBAS)==RH_R_SCF)THEN
      CALL Invoke('LowdinO',N,S,M)
      ! CALL Invoke('IRInv',N,S,M)
      ! CALL Invoke('AInv',N,S,M)
    ELSE
      CALL Invoke('AInv',N,S,M)
    ENDIF
    ! Kinetic energy matrix T
    CALL Invoke('MakeT',N,S,M)
    IF(B%BSets(1,cBAS)%HasECPs)THEN
      ! Make the ECP matrix U
      CALL Invoke('MakeU',N,S,M)
    ENDIF
  END SUBROUTINE OneEMats
  !===============================================================================
  ! COMPUTE AN ENERGY GRADIENT
  !===============================================================================
  SUBROUTINE Force(cBAS,cGEO,N,O,S,G,B,M)
    TYPE(FileNames)    :: N
    TYPE(Options)      :: O
    TYPE(State)        :: S
    TYPE(Geometries)   :: G
    TYPE(Parallel)     :: M
    TYPE(BasisSets)    :: B
    INTEGER            :: cBAS,cGEO,I,J,K,iATS,iCLONE,A1,A2
    CHARACTER(LEN=DCL) :: chGEO,chBAS
    REAL(DOUBLE)       :: GradVal,Pres,Vol,PMat
    TYPE(DBL_RNK2)     :: AuxLatF
    TYPE(DBL_VECT)     :: Ftmp

    !----------------------------------------------------------------------------!
    CALL New(S%Action,1)

    ! Initialize the force vector in HDF, clone by clone
    chGEO=IntToChar(cGEO)
    chBAS=IntToChar(cBAS)
    !
    DO iCLONE=1,G%Clones
       G%Clone(iCLONE)%Gradients%D=BIG_DBL
       G%Clone(iCLONE)%GradRMS = SQRT(G%Clone(iCLONE)%GradRMS)/DBLE(3*G%Clone(iCLONE)%NAtms)
    ENDDO
!!$    CALL MondoLog(DEBUG_MAXIMUM, "Force", "N%SCF_NAME = "//TRIM(N%SCF_NAME))
    CALL GeomArchive(cBAS,cGEO,N,O,B,G)

    ! Now evaluate the forces
!!$    CALL MondoLog(DEBUG_MAXIMUM, "Force", "State%Current = "//TRIM(IntVectToChar(S%Current)))
!!$    CALL MondoLog(DEBUG_MAXIMUM, "Force", "State%Previous = "//TRIM(IntVectToChar(S%Previous)))

    S%Action%C(1)='ForceEvaluation'
    ! The non-orthogonal response
    CALL Invoke('SForce',N,S,M)
    ! Kinetic energy piece
    CALL Invoke('TForce',N,S,M)
!!$    CALL NTHessian(cBAS,cGEO,N,G,B,S,M)
    ! Compute ECP component of the force
    IF(B%BSets(1,cBAS)%HasECPs)THEN
!!$       CALL NLATTFORCE_U(cBAS,cGEO,G,B,N,O,S,M)
       CALL Invoke('UForce',N,S,M)
    ENDIF
    ! Build density with last DM
    CALL Invoke('MakeRho',N,S,M)
    ! Coulomb part
    CALL Invoke('JForce',N,S,M)
    ! DFT exchange corrleation term
    IF(HasDFT(O%Models(cBas))) THEN
!!$       CALL NLATTFORCE_XC(cBAS,cGEO,G,B,N,S,M)
       CALL Invoke('XCForce',N,S,M)
    ENDIF
    ! Exact Hartree-Fock exchange component
    IF(HasHF(O%Models(cBas)))THEN
       CALL Invoke('GONX',N,S,M)
    ENDIF
    !
    ! Open the HDF and monkey with forces
    HDFFileID=OpenHDF(N%HFile)
    DO iCLONE=1,G%Clones
       ! Get forces and corresponding geometry straight up for each clone in HDF
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
       CALL Get(G%Clone(iCLONE),Tag_O=chGEO)
       CALL CloseHDFGroup(HDF_CurrentID)
       ! Include a Hydrostaic Presure into the lattice forces
       ! LFrc_ij = LFrc_ij + P*V*I_ii*(M^(-1))_ij
       IF(O%Pressure .NE. Zero .AND. G%Clone(iCLONE)%PBC%Dimen .GT. 0) THEN
          Pres = O%Pressure
          Vol  = G%Clone(iCLONE)%PBC%CellVolume
          DO I=1,3
             DO J=1,3
                IF(G%Clone(iCLONE)%PBC%AutoW%I(I)==1 .AND. G%Clone(iCLONE)%PBC%AutoW%I(J)==1) THEN
                   PMat = Pres*Vol*G%Clone(iCLONE)%PBC%InvBoxSh%D(I,J)
                   G%Clone(iCLONE)%PBC%LatFrc%D(I,J) = G%Clone(iCLONE)%PBC%LatFrc%D(I,J) + PMat
                ENDIF
             ENDDO
          ENDDO
       ENDIF
!      Zero forces on constrained atoms
       DO iATS=1,G%Clone(iCLONE)%NAtms
          IF(G%Clone(iCLONE)%CConstrain%I(iATS)==1 .OR. G%Clone(iCLONE)%CConstrain%I(iATS)==2)THEN
             !IF(O%Coordinates /= GRAD_INTS_OPT) THEN
             G%Clone(iCLONE)%Gradients%D(1:3,iATS)=Zero
             !ENDIF
          ENDIF
       ENDDO
!      Add additional External Forces to Atoms
       DO iATS=1,G%Clone(iCLONE)%NAtms
          IF(G%Clone(iCLONE)%CConstrain%I(iATS)==3)THEN
             G%Clone(iCLONE)%Gradients%D(1:3,iATS)=G%Clone(iCLONE)%Gradients%D(1:3,iATS)-G%Clone(iCLONE)%Fext%D(1:3,iATS)
          ENDIF
       ENDDO
!      Calculate GrandMax and GrandRMS
       G%Clone(iCLONE)%GradMax=Zero
       G%Clone(iCLONE)%GradRMS=Zero
       DO iATS=1,G%Clone(iCLONE)%NAtms
          DO J=1,3
             GradVal=G%Clone(iCLONE)%Gradients%D(J,iATS)
             G%Clone(iCLONE)%GradRMS=G%Clone(iCLONE)%GradRMS+GradVal**2
             G%Clone(iCLONE)%GradMax=MAX(G%Clone(iCLONE)%GradMax,ABS(GradVal))
          ENDDO
       ENDDO
       G%Clone(iCLONE)%GradRMS=SQRT(G%Clone(iCLONE)%GradRMS/DBLE(3*G%Clone(iCLONE)%NAtms))
    ENDDO
    ! Close up the HDF file
    CALL CloseHDF(HDFFileID)

    ! Double check zero forces on constrained atoms.
    DO iCLONE=1,G%Clones
      DO iATS=1,G%Clone(iCLONE)%NAtms
        IF(G%Clone(iCLONE)%CConstrain%I(iATS)==1 .OR. G%Clone(iCLONE)%CConstrain%I(iATS)==2)THEN
          IF((G%Clone(iCLONE)%Gradients%D(1,iATS) /= Zero) .OR. &
             (G%Clone(iCLONE)%Gradients%D(2,iATS) /= Zero) .OR. &
             (G%Clone(iCLONE)%Gradients%D(3,iATS) /= Zero)) THEN
            CALL MondoLog(DEBUG_NONE, "Force", "force on atom " &
              //TRIM(IntToChar(iATS))//" in clone " &
              //TRIM(IntToChar(iCLONE))//" is not zero!")
            G%Clone(iCLONE)%Gradients%D(1:3,iATS)=Zero
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    ! Now add in any NEB force projections
    IF(O%Grad==GRAD_TS_SEARCH_NEB) THEN
      CALL MondoLog(DEBUG_NONE, "Force", "adding in NEB forces")
      CALL NEBForce(G,O)
    ENDIF

    ! Finally, archive the whole mother
    CALL GeomArchive(cBAS,cGEO,N,O,B,G)

    ! Done with this sucka
    CALL Delete(S%Action)

    IF(O%Grad==GRAD_ONE_FORCE) THEN
       DO iCLONE=1,G%Clones
          ! Print Total Forces and Lattice Forces
          CALL New(Ftmp,3*G%Clone(iCLONE)%NAtms)
          DO iATS=1,G%Clone(iCLONE)%NAtms
             A1=3*(iATS-1)+1
             A2=3*iATS
             Ftmp%D(A1:A2) = -G%Clone(iCLONE)%Gradients%D(1:3,iATS)
          ENDDO
          PrintFlags%Key=DEBUG_MAXIMUM
          PrintFlags%MM=DEBUG_FRC
          CALL Print_Force(G%Clone(iCLONE),Ftmp,'Force')
          CALL Print_LatForce(G%Clone(iCLONE),G%Clone(iCLONE)%PBC%LatFrc%D,'Lattice Force')
          CALL Delete(Ftmp)
       ENDDO
    ENDIF

    ! Print out the positions.
    DO iCLONE=LBOUND(G%Clone, 1), UBOUND(G%Clone, 1)
      CALL MondoLog(DEBUG_NONE, "Force", "Positions (in A) and forces (in eV/A)", &
        "Clone "//TRIM(IntToChar(iCLONE)))
      DO iATS=1, G%Clone(iCLONE)%NAtms
        CALL MondoLog(DEBUG_NONE, "Force", &
          TRIM(G%Clone(iCLONE)%AtNam%C(iATS))//" "// &
          TRIM(DblToMedmChar(G%Clone(iCLONE)%Carts%D(1, iATS)*AUToAngstroms))//" "// &
          TRIM(DblToMedmChar(G%Clone(iCLONE)%Carts%D(2, iATS)*AUToAngstroms))//" "// &
          TRIM(DblToMedmChar(G%Clone(iCLONE)%Carts%D(3, iATS)*AUToAngstroms))//" "// &
          TRIM(DblToMedmChar(-G%Clone(iCLONE)%Gradients%D(1, iATS)*au2eV/AUToAngstroms))//" "// &
          TRIM(DblToMedmChar(-G%Clone(iCLONE)%Gradients%D(2, iATS)*au2eV/AUToAngstroms))//" "// &
          TRIM(DblToMedmChar(-G%Clone(iCLONE)%Gradients%D(3, iATS)*au2eV/AUToAngstroms))//" "// &
          TRIM(IntToChar(G%Clone(iCLONE)%CConstrain%I(iATS))), &
          "Clone "//TRIM(IntToChar(iCLONE)))
      ENDDO
    ENDDO

  END SUBROUTINE Force
  !===============================================================================
  ! Numerically compute Lattice Forces for J
  !===============================================================================
  SUBROUTINE NLATTFORCE_J(cBAS,cGEO,G,B,N,O,S,M)
    TYPE(FileNames)    :: N
    TYPE(Options)      :: O
    TYPE(State)        :: S
    TYPE(Geometries)   :: G
    TYPE(BasisSets)    :: B
    TYPE(Parallel)     :: M
    INTEGER            :: cBAS,cGEO,iCLONE,I,J,II
    CHARACTER(LEN=DCL) :: chGEO,chBAS,chSCF
    REAL(DOUBLE)       :: DDelta,Lat00,E1,E2
    CHARACTER(LEN=DCL) :: TrixName
    REAL(DOUBLE),DIMENSION(3,3) :: LatFrc_J
#ifdef PARRALEL
    TYPE(DBCSR)        :: P,J1,J2,J3
#else
    TYPE(BCSR)         :: P,J1,J2,J3
#endif
    !
    DDelta = 1.D-4
    DO iCLONE=1,G%Clones
      !
      LatFrc_J = Zero
      !
      chGEO=IntToChar(cGEO)
      chBAS=IntToChar(cBAS)
      chSCF=IntToChar(S%Current%I(1))
      CALL New(BSiz,G%Clone(iCLONE)%NAtms)
      CALL New(OffS,G%Clone(iCLONE)%NAtms)
      !
      HDFFileID=OpenHDF(N%HFile)
      HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
      !
      O%Thresholds(cBAS)%TwoE=O%Thresholds(cBAS)%TwoE
      CALL Put(O%Thresholds(cBAS),chBAS)
      !
      CALL Get(G%Clone(iCLONE)%PBC%LatFrc,'latfrc',Tag_O=chGEO)
      CALL Get(G%Clone(iCLONE)%Gradients,'Gradients',Tag_O=chGEO)
      CALL CloseHDFGroup(HDF_CurrentID)
      CALL CloseHDF(HDFFileID)
      !
      NAToms=G%Clone(iCLONE)%NAtms
      MaxAtms=B%MxAts(cBAS)
      MaxBlks=B%MxBlk(cBAS)
      MaxNon0=B%MxN0s(cBAS)
      NBasF=B%BSets(iCLONE,cBAS)%NBasF
      BSiz%I=B%BSiz(iCLONE,cBAS)%I
      OffS%I=B%OffS(iCLONE,cBAS)%I
      MaxBlkSize=0
      DO II=1,G%Clone(iCLONE)%NAtms
        MaxBlkSize=MAX(MaxBlkSize,BSiz%I(II))
      ENDDO
      !
      CALL New(P)
      CALL New(J1)
      CALL New(J2)
      CALL New(J3)
      !
      TrixName=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)//'_Geom#'//TRIM(chGEO)//'_Base#'//TRIM(chBAS) &
           //'_Cycl#'//TRIM(chSCF)//'_Clone#'//TRIM(IntToChar(iCLONE))//'.D'
      CALL Get(P,TrixName)
      !
      DO I=1,3
        DO J=1,3
          IF(G%Clone(iCLONE)%PBC%AutoW%I(I) == 1 .AND. G%Clone(iCLONE)%PBC%AutoW%I(J) == 1) THEN
            Lat00 = G%Clone(iCLONE)%PBC%BoxShape%D(I,J)
            !
            G%Clone(iCLONE)%PBC%BoxShape%D(I,J) =  Lat00+DDelta
            CALL MkGeomPeriodic(G%Clone(iCLONE))
            CALL GeomArchive(cBAS,cGEO,N,O,B,G)
            CALL Invoke('MakeRho' ,N,S,M)
            CALL Invoke('MakePFFT',N,S,M)
            CALL Invoke('QCTC'    ,N,S,M)
            TrixName=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)//'_Geom#'//TRIM(chGEO)//'_Base#'//TRIM(chBAS) &
                 //'_Cycl#'//TRIM(chSCF)//'_Clone#'//TRIM(IntToChar(iCLONE))//'.J'
            CALL Get(J1,TrixName)
            HDFFileID=OpenHDF(N%HFile)
            HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
            CALL Get(E1,'E_NuclearTotal')
            CALL CloseHDFGroup(HDF_CurrentID)
            CALL CloseHDF(HDFFileID)
            !
            G%Clone(iCLONE)%PBC%BoxShape%D(I,J) =  Lat00-DDelta
            CALL MkGeomPeriodic(G%Clone(iCLONE))
            CALL GeomArchive(cBAS,cGEO,N,O,B,G)
            CALL Invoke('MakeRho' ,N,S,M)
            CALL Invoke('MakePFFT',N,S,M)
            CALL Invoke('QCTC',N,S,M)
            TrixName=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)//'_Geom#'//TRIM(chGEO)//'_Base#'//TRIM(chBAS) &
                 //'_Cycl#'//TRIM(chSCF)//'_Clone#'//TRIM(IntToChar(iCLONE))//'.J'
            CALL Get(J2,TrixName)
            HDFFileID=OpenHDF(N%HFile)
            HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
            CALL Get(E2,'E_NuclearTotal')
            CALL CloseHDFGroup(HDF_CurrentID)
            CALL CloseHDF(HDFFileID)
            !
            CALL Multiply(J2,-One)
            CALL Add(J1,J2,J3)
            CALL Multiply(P,J3,J1)
            ! LatFrc_J(I,J) =  Trace(J1)/(Two*DDelta)
            ! LatFrc_J(I,J) =  (E1-E2)/(Two*DDelta)
            LatFrc_J(I,J) =  (Trace(J1) + (E1-E2))/(Two*DDelta)
            !
            G%Clone(iCLONE)%PBC%BoxShape%D(I,J) =  Lat00
            CALL MkGeomPeriodic(G%Clone(iCLONE))
            CALL GeomArchive(cBAS,cGEO,N,O,B,G)
            !
          ENDIF
        ENDDO
      ENDDO
      ! Update Lattice Forces
      G%Clone(iCLONE)%PBC%LatFrc%D = G%Clone(iCLONE)%PBC%LatFrc%D + LatFrc_J
      HDFFileID=OpenHDF(N%HFile)
      HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
      !
      O%Thresholds(cBAS)%TwoE=O%Thresholds(cBAS)%TwoE*1.D4
      CALL Put(O%Thresholds(cBAS),chBAS)
      !
      CALL Put(G%Clone(iCLONE)%PBC%LatFrc,'latfrc',Tag_O=chGEO)
      CALL Put(G%Clone(iCLONE)%Gradients,'Gradients',Tag_O=chGEO)
      CALL CloseHDFGroup(HDF_CurrentID)
      CALL CloseHDF(HDFFileID)
      ! Print J Lattice Forces
!!$       WRITE(*,*) 'LatFrc_J NUM'
!!$       DO I=1,3
!!$          WRITE(*,*) (LatFrc_J(I,J),J=1,3)
!!$       ENDDO
      !
      CALL Delete(BSiz)
      CALL Delete(OffS)
      CALL Delete(P)
      CALL Delete(J1)
      CALL Delete(J2)
      CALL Delete(J3)
    ENDDO
    !
  END SUBROUTINE NLATTFORCE_J
  !===============================================================================
  ! Numerically compute Lattice Forces for J
  !===============================================================================
  SUBROUTINE NLATTFORCE_X(cBAS,cGEO,G,B,N,O,S,M)
    TYPE(FileNames)    :: N
    TYPE(Options)      :: O
    TYPE(State)        :: S
    TYPE(Geometries)   :: G
    TYPE(BasisSets)    :: B
    TYPE(Parallel)     :: M
    INTEGER            :: cBAS,cGEO,iCLONE,I,J,II
    CHARACTER(LEN=DCL) :: chGEO,chBAS,chSCF
    REAL(DOUBLE)       :: DDelta,Lat00,KScale,Vec(6),BoxShape(3,3)
    CHARACTER(LEN=DCL) :: TrixName
    REAL(DOUBLE),DIMENSION(3,3) :: LatFrc_X
#ifdef PARRALEL
    TYPE(DBCSR)        :: P,K1,K2,K3
#else
    TYPE(BCSR)         :: P,K1,K2,K3
#endif
    !
    DDelta = 1.D-3
    DO iCLONE=1,G%Clones
      !
      LatFrc_X = Zero
      !
      chGEO=IntToChar(cGEO)
      chBAS=IntToChar(cBAS)
      chSCF=IntToChar(S%Current%I(1))
      CALL New(BSiz,G%Clone(iCLONE)%NAtms)
      CALL New(OffS,G%Clone(iCLONE)%NAtms)
      !
      HDFFileID=OpenHDF(N%HFile)
      HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
      CALL Get(G%Clone(iCLONE)%PBC%LatFrc,'latfrc',Tag_O=chGEO)
      CALL Get(G%Clone(iCLONE)%Gradients,'Gradients',Tag_O=chGEO)
      CALL CloseHDFGroup(HDF_CurrentID)
      CALL CloseHDF(HDFFileID)
      !
      NAToms=G%Clone(iCLONE)%NAtms
      MaxAtms=B%MxAts(cBAS)
      MaxBlks=B%MxBlk(cBAS)
      MaxNon0=B%MxN0s(cBAS)
      NBasF=B%BSets(iCLONE,cBAS)%NBasF
      BSiz%I=B%BSiz(iCLONE,cBAS)%I
      OffS%I=B%OffS(iCLONE,cBAS)%I
      MaxBlkSize=0
      DO II=1,G%Clone(iCLONE)%NAtms
        MaxBlkSize=MAX(MaxBlkSize,BSiz%I(II))
      ENDDO
      !
      CALL New(P)
      CALL New(K1)
      CALL New(K2)
      CALL New(K3)
      !
      TrixName=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)//'_Geom#'//TRIM(chGEO)//'_Base#'//TRIM(chBAS) &
           //'_Cycl#'//TRIM(chSCF)//'_Clone#'//TRIM(IntToChar(iCLONE))//'.D'
      CALL Get(P,TrixName)
      chSCF=IntToChar(S%Current%I(1))
      !
      BoxShape=G%Clone(iCLONE)%PBC%BoxShape%D

      DO I=1,3
        DO J=1,3
          IF(G%Clone(iCLONE)%PBC%AutoW%I(I) == 1 .AND. G%Clone(iCLONE)%PBC%AutoW%I(J) == 1) THEN
            G%Clone(iCLONE)%PBC%BoxShape%D=BoxShape
            Lat00 = G%Clone(iCLONE)%PBC%BoxShape%D(I,J)
            !
            G%Clone(iCLONE)%PBC%BoxShape%D(I,J) =  Lat00+DDelta
            CALL CalcBoxPars(Vec,G%Clone(iCLONE)%PBC%BoxShape%D)
            CALL BoxParsToCart(Vec,G%Clone(iCLONE)%PBC%BoxShape%D)
            CALL MkGeomPeriodic(G%Clone(iCLONE))
            CALL GeomArchive(cBAS,cGEO,N,O,B,G)
            CALL Invoke('ONX'    ,N,S,M)
            TrixName=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)//'_Geom#'//TRIM(chGEO)//'_Base#'//TRIM(chBAS) &
                 //'_Cycl#'//TRIM(chSCF)//'_Clone#'//TRIM(IntToChar(iCLONE))//'.K'
            CALL Get(K1,TrixName)
            !
            G%Clone(iCLONE)%PBC%BoxShape%D=BoxShape
            G%Clone(iCLONE)%PBC%BoxShape%D(I,J) =  Lat00-DDelta
            CALL CalcBoxPars(Vec,G%Clone(iCLONE)%PBC%BoxShape%D)
            CALL BoxParsToCart(Vec,G%Clone(iCLONE)%PBC%BoxShape%D)
            CALL MkGeomPeriodic(G%Clone(iCLONE))
            CALL GeomArchive(cBAS,cGEO,N,O,B,G)
            CALL Invoke('ONX',N,S,M)
            TrixName=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)//'_Geom#'//TRIM(chGEO)//'_Base#'//TRIM(chBAS) &
                 //'_Cycl#'//TRIM(chSCF)//'_Clone#'//TRIM(IntToChar(iCLONE))//'.K'
            CALL Get(K2,TrixName)
            !
            CALL Multiply(K2,-One)
            CALL Add(K1,K2,K3)
            CALL Multiply(P,K3,K1)
            LatFrc_X(I,J) = Trace(K1)/(Two*DDelta)
            !
            G%Clone(iCLONE)%PBC%BoxShape%D(I,J) =  Lat00
            CALL MkGeomPeriodic(G%Clone(iCLONE))
            CALL GeomArchive(cBAS,cGEO,N,O,B,G)
            !
          ENDIF
        ENDDO
      ENDDO
      ! Multiply By KScale
      KScale   = ExactXScale(O%Models(cBAS))
      LatFrc_X = KScale*LatFrc_X
      ! Update X Lattice Force
      G%Clone(iCLONE)%PBC%LatFrc%D = G%Clone(iCLONE)%PBC%LatFrc%D + LatFrc_X
      HDFFileID=OpenHDF(N%HFile)
      HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
!!$       CALL Put(G%Clone(iCLONE)%PBC%LatFrc,'latfrc',Tag_O=chGEO)
!!$       CALL Put(G%Clone(iCLONE)%Gradients,'Gradients',Tag_O=chGEO)
      CALL CloseHDFGroup(HDF_CurrentID)
      CALL CloseHDF(HDFFileID)
      ! Print X Lattice Forces
      WRITE(*,*) 'LatFrc_X NUM'
      DO I=1,3
        WRITE(*,*) (LatFrc_X(I,J),J=1,3)
      ENDDO
      !
      CALL Delete(BSiz)
      CALL Delete(OffS)
      CALL Delete(P)
      CALL Delete(K1)
      CALL Delete(K2)
      CALL Delete(K3)
      !
    ENDDO
    !
  END SUBROUTINE NLATTFORCE_X
  !===============================================================================
  ! Numerically compute Lattice Forces for J
  !===============================================================================
  SUBROUTINE NLATTFORCE_U(cBAS,cGEO,G,B,N,O,S,M)
    TYPE(FileNames)      :: N
    TYPE(Options)        :: O
    TYPE(State)          :: S
    TYPE(Geometries)     :: G
    TYPE(BasisSets)      :: B
    TYPE(Parallel)       :: M
    INTEGER              :: cBAS,cGEO,iCLONE,I,J,II
    CHARACTER(LEN=DCL)   :: chGEO,chBAS,chSCF
    REAL(DOUBLE)         :: DDelta,Lat00,KScale,Vec(6),BoxShape(3,3)
    CHARACTER(LEN=DCL)   :: TrixName
    REAL(DOUBLE),DIMENSION(3,3) :: LatFrc_U
#ifdef PARRALEL
    TYPE(DBCSR)        :: P,U1,U2,U3
#else
    TYPE(BCSR)         :: P,U1,U2,U3
#endif
    !
    DDelta = 1.D-3
    DO iCLONE=1,G%Clones
      !
      LatFrc_U = Zero
      !
      chGEO=IntToChar(cGEO)
      chBAS=IntToChar(cBAS)
      chSCF=IntToChar(S%Current%I(1))
      CALL New(BSiz,G%Clone(iCLONE)%NAtms)
      CALL New(OffS,G%Clone(iCLONE)%NAtms)
      !
      HDFFileID=OpenHDF(N%HFile)
      HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
      CALL Get(G%Clone(iCLONE)%PBC%LatFrc,'latfrc',Tag_O=chGEO)
      CALL Get(G%Clone(iCLONE)%Gradients,'Gradients',Tag_O=chGEO)
      CALL CloseHDFGroup(HDF_CurrentID)
      CALL CloseHDF(HDFFileID)
      !
      NAToms=G%Clone(iCLONE)%NAtms
      MaxAtms=B%MxAts(cBAS)
      MaxBlks=B%MxBlk(cBAS)
      MaxNon0=B%MxN0s(cBAS)
      NBasF=B%BSets(iCLONE,cBAS)%NBasF
      BSiz%I=B%BSiz(iCLONE,cBAS)%I
      OffS%I=B%OffS(iCLONE,cBAS)%I
      MaxBlkSize=0
      DO II=1,G%Clone(iCLONE)%NAtms
        MaxBlkSize=MAX(MaxBlkSize,BSiz%I(II))
      ENDDO
      !
      CALL New(P)
      CALL New(U1)
      CALL New(U2)
      CALL New(U3)
      !
      TrixName=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)//'_Geom#'//TRIM(chGEO)//'_Base#'//TRIM(chBAS) &
           //'_Cycl#'//TRIM(chSCF)//'_Clone#'//TRIM(IntToChar(iCLONE))//'.D'
      CALL Get(P,TrixName)
      chSCF=IntToChar(S%Current%I(1))
      !
      BoxShape=G%Clone(iCLONE)%PBC%BoxShape%D
      DO I=1,3
        DO J=1,3
          IF(G%Clone(iCLONE)%PBC%AutoW%I(I) == 1 .AND. G%Clone(iCLONE)%PBC%AutoW%I(J) == 1) THEN
            G%Clone(iCLONE)%PBC%BoxShape%D=BoxShape
            Lat00 = G%Clone(iCLONE)%PBC%BoxShape%D(I,J)
            !
            G%Clone(iCLONE)%PBC%BoxShape%D(I,J) =  Lat00+DDelta
            CALL CalcBoxPars(Vec,G%Clone(iCLONE)%PBC%BoxShape%D)
            CALL BoxParsToCart(Vec,G%Clone(iCLONE)%PBC%BoxShape%D)
            CALL MkGeomPeriodic(G%Clone(iCLONE))
            CALL GeomArchive(cBAS,cGEO,N,O,B,G)
            CALL Invoke('MakeU'    ,N,S,M)
            TrixName=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)//'_Geom#'//TRIM(chGEO)//'_Base#'//TRIM(chBAS) &
                 //'_Clone#'//TRIM(IntToChar(iCLONE))//'.U'
            CALL Get(U1,TrixName)
            !
            G%Clone(iCLONE)%PBC%BoxShape%D=BoxShape
            G%Clone(iCLONE)%PBC%BoxShape%D(I,J) =  Lat00-DDelta
            CALL CalcBoxPars(Vec,G%Clone(iCLONE)%PBC%BoxShape%D)
            CALL BoxParsToCart(Vec,G%Clone(iCLONE)%PBC%BoxShape%D)
            CALL MkGeomPeriodic(G%Clone(iCLONE))
            CALL GeomArchive(cBAS,cGEO,N,O,B,G)
            CALL Invoke('MakU',N,S,M)
            TrixName=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)//'_Geom#'//TRIM(chGEO)//'_Base#'//TRIM(chBAS) &
                 //'_Clone#'//TRIM(IntToChar(iCLONE))//'.U'
            CALL Get(U2,TrixName)
            !
            CALL Multiply(U2,-One)
            CALL Add(U1,U2,U3)
            CALL Multiply(P,U3,U1)
            LatFrc_U(I,J) = Trace(U1)/(Two*DDelta)
            !
            G%Clone(iCLONE)%PBC%BoxShape%D(I,J) =  Lat00
            CALL MkGeomPeriodic(G%Clone(iCLONE))
            CALL GeomArchive(cBAS,cGEO,N,O,B,G)
            !
          ENDIF
        ENDDO
      ENDDO
      ! Update U Lattice Force
      G%Clone(iCLONE)%PBC%LatFrc%D = G%Clone(iCLONE)%PBC%LatFrc%D + LatFrc_U
      HDFFileID=OpenHDF(N%HFile)
      HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
      CALL Put(G%Clone(iCLONE)%PBC%LatFrc,'latfrc',Tag_O=chGEO)
      CALL Put(G%Clone(iCLONE)%Gradients,'Gradients',Tag_O=chGEO)
      CALL CloseHDFGroup(HDF_CurrentID)
      CALL CloseHDF(HDFFileID)
      ! Print U Lattice Forces
      WRITE(*,*) 'LatFrc_U NUM'
      DO I=1,3
        WRITE(*,*) (LatFrc_U(I,J),J=1,3)
      ENDDO
      !
      CALL Delete(BSiz)
      CALL Delete(OffS)
      CALL Delete(P)
      CALL Delete(U1)
      CALL Delete(U2)
      CALL Delete(U3)
      !
    ENDDO
    !
  END SUBROUTINE NLATTFORCE_U
  !===============================================================================
  ! Numerically compute Lattice Forces for J
  !===============================================================================
  SUBROUTINE NLATTFORCE_XC(cBAS,cGEO,G,B,N,S,M)
    TYPE(FileNames)    :: N
    TYPE(Options)      :: O
    TYPE(State)        :: S
    TYPE(Geometries)   :: G
    TYPE(BasisSets)    :: B
    TYPE(Parallel)     :: M
    INTEGER            :: cBAS,cGEO,iCLONE,I,J,II
    CHARACTER(LEN=DCL) :: chGEO,chBAS,chSCF
    REAL(DOUBLE)       :: DDelta,Lat00,E1,E2
    CHARACTER(LEN=DCL) :: TrixName
    REAL(DOUBLE),DIMENSION(3,3) :: LatFrc_XC
#ifdef PARRALEL
    TYPE(DBCSR)        :: P,K1,K2,K3
#else
    TYPE(BCSR)         :: P,K1,K2,K3
#endif
    !
    DDelta = 1.D-3
    DO iCLONE=1,G%Clones
      !
      LatFrc_XC = Zero
      !
      chGEO=IntToChar(cGEO)
      chBAS=IntToChar(cBAS)
      chSCF=IntToChar(S%Current%I(1))
      !
      HDFFileID=OpenHDF(N%HFile)
      HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
      CALL Get(G%Clone(iCLONE)%Gradients,'Gradients',Tag_O=chGEO)
      CALL Get(G%Clone(iCLONE)%PBC%LatFrc,'latfrc',Tag_O=chGEO)
      CALL CloseHDFGroup(HDF_CurrentID)
      CALL CloseHDF(HDFFileID)
      !
      DO I=1,3
        DO J=1,3
          IF(G%Clone(iCLONE)%PBC%AutoW%I(I) == 1 .AND. G%Clone(iCLONE)%PBC%AutoW%I(J) == 1) THEN
            Lat00 = G%Clone(iCLONE)%PBC%BoxShape%D(I,J)
            !
            G%Clone(iCLONE)%PBC%BoxShape%D(I,J) =  Lat00+DDelta
            CALL MkGeomPeriodic(G%Clone(iCLONE))
            CALL GeomArchive(cBAS,cGEO,N,O,B,G)
            CALL Invoke('MakeRho' ,N,S,M)
            CALL Invoke('HiCu'    ,N,S,M)
            !
            HDFFileID=OpenHDF(N%HFile)
            HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
            CALL Get(E1,'Exc')
            CALL Get(G%Clone(iCLONE)%Gradients,'Gradients',Tag_O=chGEO)
            CALL CloseHDFGroup(HDF_CurrentID)
            CALL CloseHDF(HDFFileID)
            !
            G%Clone(iCLONE)%PBC%BoxShape%D(I,J) =  Lat00-DDelta
            CALL MkGeomPeriodic(G%Clone(iCLONE))
            CALL GeomArchive(cBAS,cGEO,N,O,B,G)
            CALL Invoke('MakeRho' ,N,S,M)
            CALL Invoke('HiCu',N,S,M)
            !
            HDFFileID=OpenHDF(N%HFile)
            HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
            CALL Get(E2,'Exc')
            CALL CloseHDFGroup(HDF_CurrentID)
            CALL CloseHDF(HDFFileID)
            !
            LatFrc_XC(I,J) =  (E1-E2)/(Two*DDelta)
            !
            G%Clone(iCLONE)%PBC%BoxShape%D(I,J) =  Lat00
            CALL MkGeomPeriodic(G%Clone(iCLONE))
            CALL GeomArchive(cBAS,cGEO,N,O,B,G)
            !
          ENDIF
        ENDDO
        G%Clone(iCLONE)%PBC%LatFrc%D=G%Clone(iCLONE)%PBC%LatFrc%D+LatFrc_XC
        CALL GeomArchive(cBAS,cGEO,N,O,B,G)
      ENDDO
      ! Print J Lattice Forces
      ! CALL Print_LatForce(G,LatFrc_XC%D,'XC Lattice Force Numerical')
      !
    ENDDO
    !
  END SUBROUTINE NLATTFORCE_XC
  !===============================================================================
  ! Numerically compute gradients of the exact HF exchange
  !===============================================================================
  SUBROUTINE NXForce(cBAS,cGEO,N,G,B,S,M)
    TYPE(FileNames)  :: N
    TYPE(Options)    :: O
    TYPE(State)      :: S
    TYPE(Geometries) :: G
    TYPE(BasisSets)  :: B
    TYPE(Parallel)   :: M
    TYPE(DBL_RNK2)   :: GradAux1,GradAux2
    INTEGER          :: cBAS,cGEO,J,iATS,iCLONE
    CHARACTER(LEN=DCL) :: chGEO,chBAS,chSCF
    TYPE(BCSR),DIMENSION(G%Clones)   :: P
    TYPE(CRDS),DIMENSION(G%Clones)   :: GTmp
    TYPE(BCSR)                       :: K
    REAL(DOUBLE),DIMENSION(G%Clones,G%Clone(1)%NAtms*3) :: FX
    REAL(DOUBLE),DIMENSION(G%Clones,2)        :: EX
    INTEGER                          :: AtA,IX,II,IA,A1,A2,IS
    REAL(DOUBLE),PARAMETER           :: DDelta = 1.D-3
    CHARACTER(LEN=DCL)               :: TrixName
    !------------------------------------------------------------------------------
    chGEO=IntToChar(cGEO)
    chBAS=IntToChar(cBAS)
    chSCF=IntToChar(S%Current%I(1)+1)
    CALL New(BSiz,G%Clone(1)%NAtms)
    CALL New(OffS,G%Clone(1)%NAtms)
    CALL New(GradAux1,(/3,G%Clone(1)%NAtms/))
    CALL New(GradAux2,(/3,G%Clone(1)%NAtms/))
    DO iCLONE=1,G%Clones
      !
      ! Load current gradients, we need that, cause we save the geo through
      ! GeoArchive and G%..%Gradients have been det to BIG_DBL previously.
      !
      HDFFileID=OpenHDF(N%HFile)
      HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
      CALL Get(G%Clone(iCLONE)%PBC%LatFrc,'latfrc',Tag_O=chGEO)
      CALL Get(G%Clone(iCLONE)%Gradients,'Gradients',Tag_O=chGEO)
      CALL CloseHDFGroup(HDF_CurrentID)
      CALL CloseHDF(HDFFileID)
      !
      ! Load globals
      NAToms=G%Clone(1)%NAtms
      MaxAtms=B%MxAts(cBAS)
      MaxBlks=B%MxBlk(cBAS)
      MaxNon0=B%MxN0s(cBAS)
      NBasF=B%BSets(iCLONE,cBAS)%NBasF
      BSiz%I=B%BSiz(iCLONE,cBAS)%I
      OffS%I=B%OffS(iCLONE,cBAS)%I
      MaxBlkSize=0
      DO II=1,G%Clone(1)%NAtms
        MaxBlkSize=MAX(MaxBlkSize,BSiz%I(II))
      ENDDO
      ! Set temporary geometries
      GTmp(iCLONE)%NAtms=G%Clone(iCLONE)%NAtms
      CALL New_CRDS(GTmp(iCLONE))
      ! Get the density matrix for this clone
      TrixName=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)//'_Geom#'//TRIM(chGEO)//'_Base#'//TRIM(chBAS)//'_Cycl#'//TRIM(chSCF) &
           //'_Clone#'//TRIM(IntToChar(iCLONE))//'.D'
      CALL Get(P(iCLONE),TrixName)
    ENDDO
    chSCF=IntToChar(S%Current%I(1))
    FX=Zero
    DO AtA=1,G%Clone(1)%NAtms
      DO IX=1,3
        DO II=1,2
          DO iCLONE=1,G%Clones
            !
            ! Move the atom.
            IF(II==1) THEN
              G%Clone(iCLONE)%Carts%D(IX,AtA)=GTmp(iCLONE)%Carts%D(IX,AtA)+DDelta
            ELSEIF(II==2) THEN
              G%Clone(iCLONE)%Carts%D(IX,AtA)=GTmp(iCLONE)%Carts%D(IX,AtA)-DDelta
            ENDIF
          ENDDO
          !
          CALL GeomArchive(cBAS,cGEO,N,O,B,G)
          ! vwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvw>>>
          ! Move back the atom.
          DO iCLONE=1,G%Clones
            IF(II==1) THEN
              G%Clone(iCLONE)%Carts%D(IX,AtA)=GTmp(iCLONE)%Carts%D(IX,AtA)
            ELSEIF(II==2) THEN
              G%Clone(iCLONE)%Carts%D(IX,AtA)=GTmp(iCLONE)%Carts%D(IX,AtA)
            ENDIF
          ENDDO
          ! vwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvw<<<
          CALL Invoke('ONX',N,S,M)
          DO iCLONE=1,G%Clones
            ! Load globals
            NAToms=G%Clone(1)%NAtms
            MaxAtms=B%MxAts(cBAS)
            MaxBlks=B%MxBlk(cBAS)
            MaxNon0=B%MxN0s(cBAS)
            NBasF=B%BSets(iCLONE,cBAS)%NBasF
            BSiz%I=B%BSiz(iCLONE,cBAS)%I
            OffS%I=B%OffS(iCLONE,cBAS)%I
            ! Get the exact HF exchange matrix from disk
            TrixName=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)//'_Geom#'//TRIM(chGEO)//'_Base#'//TRIM(chBAS)//'_Cycl#'//TRIM(chSCF) &
                 //'_Clone#'//TRIM(IntToChar(iCLONE))//'.K'
            CALL Get(K,TrixName)
            EX(iCLONE,II)=Trace(P(iCLONE),K)
          ENDDO
        ENDDO
        IA=3*(AtA-1)+IX
        DO iCLONE=1,G%Clones
          FX(iCLONE,IA)=(EX(iCLONE,1)-EX(iCLONE,2))/(Two*DDelta)
          !vwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvw
#ifdef NGONXINFO
          write(*,*) 'FX(',AtA,IX,')=',FX(iCLONE,IA)
#endif
          !vwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvw
        ENDDO
      ENDDO
    ENDDO
    !
    CALL GeomArchive(cBAS,cGEO,N,O,B,G)
    !
    !vwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvw
#ifdef NGONXINFO
    DO iCLONE=1,G%Clones
      WRITE(*,'(A,I3,A,E20.12)') ' Total XForce(',iCLONE,') =',SUM(FX(iCLONE,:))
    ENDDO
#endif
    !vwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvw
    ! Add in the forces to the global gradient and put back to HDF
    HDFFileID=OpenHDF(N%HFile)
    DO iCLONE=1,G%Clones
      HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
      CALL Get(GradAux1,'gradients',Tag_O=chGEO)
      CALL CartRNK1ToCartRNK2(FX(iCLONE,:),GradAux2%D)
      GradAux1%D=GradAux1%D+GradAux2%D
      CALL Put(GradAux1,'gradients',Tag_O=chGEO)
      ! Close the group
      CALL CloseHDFGroup(HDF_CurrentID)
    ENDDO
    CALL CloseHDF(HDFFileID)
    DO iCLONE=1,G%Clones
      G%Clone(iCLONE)%Carts%D=GTmp(iCLONE)%Carts%D
      CALL Delete(GTmp(iCLONE))
      CALL Delete(P(iCLONE))
    ENDDO
    CALL Delete(GradAux1)
    CALL Delete(GradAux2)
    CALL Delete(BSiz)
    CALL Delete(OffS)
    CALL Delete(K)
  END SUBROUTINE NXForce
  !===============================================================================
  ! Clean the Scratch
  !===============================================================================
  SUBROUTINE CleanScratch(C,iGEO,DoingMD_O)
    TYPE(Controls)                 :: C
    INTEGER                        :: iGEO, iCLONE
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: RemoveFile,chGEO
    LOGICAL,OPTIONAL               :: DoingMD_O
    LOGICAL                        :: DoingMD

    DoingMD=.FALSE.
    IF(PRESENT(DoingMD_O)) DoingMD=DoingMD_O

    !CALL MondoLog(DEBUG_MAXIMUM, "CleanScratch", "doing MD = "//TRIM(LogicalToChar(DoingMD)) &
    !  //", iGEO = "//TRIM(IntToChar(iGEO)))

    chGEO = IntToChar(iGEO)
    RemoveFile=TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)//'_*Geom#'//TRIM(chGEO)//"_*"
    CALL MondoLog(DEBUG_NONE, "CleanScratch", "removing "//TRIM(RemoveFile))
    CALL FileRemove(RemoveFile)
    IF(DoingMD) THEN
      RemoveFile=TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)//'_*G#'//TRIM(chGEO)//"_*"
      !CALL MondoLog(DEBUG_MAXIMUM, "CleanScratch", "removing "//TRIM(RemoveFile))
      CALL FileRemove(RemoveFile)
    ENDIF

    !IF(iGEO >= 2) THEN
    !  CALL MondoLog(DEBUG_MAXIMUM, "CleanScratch", "identifying objects in hdf group...")

    !  ! Open HDF file.
    !  HDFFileID = OpenHDF(C%Nams%HFile)
    !  DO iCLONE=1, C%Geos%Clones
    !    CALL HDF5DeleteObject(HDFFileID, "Clone #"//TRIM(IntToChar(iCLONE)), "^.*"//TRIM(IntToChar(iGEO))//"$")
    !  ENDDO

    !  ! Close the HDF.
    !  CALL CloseHDF(HDFFileID)
    !ENDIF

  END SUBROUTINE CleanScratch
  !===============================================================================
  !
  !===============================================================================
  SUBROUTINE NTHessian(cBAS,cGEO,N,G,B,S,M)
    TYPE(FileNames)  :: N
    TYPE(Options)    :: O
    TYPE(State)      :: S
    TYPE(Geometries) :: G
    TYPE(BasisSets)  :: B
    TYPE(Parallel)   :: M
    TYPE(DBL_RNK2)   :: HT
    INTEGER          :: cBAS,cGEO,J,iATS,iCLONE
    CHARACTER(LEN=DCL) :: chGEO,chBAS,chSCF
    TYPE(BCSR) :: P
    TYPE(CRDS) :: GTmp
    REAL(DOUBLE),DIMENSION(3,G%Clone(1)%NAtms*3,2)        :: FT
    INTEGER                          :: AtA,AtB,IX,II,IA,IB,A1,A2,IS
    REAL(DOUBLE),PARAMETER           :: DDelta = 1.D-3
    CHARACTER(LEN=DCL)               :: TrixName
    !------------------------------------------------------------------------------
    write(*,*) 'NTHessian NTHessian NTHessian NTHessian NTHessian NTHessian'
    chGEO=IntToChar(cGEO)
    chBAS=IntToChar(cBAS)
    chSCF=IntToChar(S%Current%I(1)+1)
    CALL New(BSiz,G%Clone(1)%NAtms)
    CALL New(OffS,G%Clone(1)%NAtms)
    CALL New(HT,(/G%Clone(1)%NAtms*3,G%Clone(1)%NAtms*3/))
    HT%D=0.0D0
    !
    ! Load current gradients, we need that, cause we save the geo through
    ! GeoArchive and G%..%Gradients have been set to BIG_DBL previously.
    !
    HDFFileID=OpenHDF(N%HFile)
    HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
    CALL Get(G%Clone(1)%PBC%LatFrc,'latfrc',Tag_O=chGEO)
    CALL Get(G%Clone(1)%Gradients,'Gradients',Tag_O=chGEO)
    CALL CloseHDFGroup(HDF_CurrentID)
    CALL CloseHDF(HDFFileID)
    !
    ! Load globals
    NAToms=G%Clone(1)%NAtms
    MaxAtms=B%MxAts(cBAS)
    MaxBlks=B%MxBlk(cBAS)
    MaxNon0=B%MxN0s(cBAS)
    NBasF=B%BSets(1,cBAS)%NBasF
    BSiz%I=B%BSiz(1,cBAS)%I
    OffS%I=B%OffS(1,cBAS)%I
    MaxBlkSize=0
    DO II=1,NAtoms
      MaxBlkSize=MAX(MaxBlkSize,BSiz%I(II))
    ENDDO
    ! Set temporary geometries
    GTmp%NAtms=G%Clone(1)%NAtms
    CALL New_CRDS(GTmp)
    GTmp%Carts%D=G%Clone(1)%Carts%D
    GTmp%Gradients%D=G%Clone(1)%Gradients%D
    ! Get the density matrix for this clone
    TrixName=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)//'_Geom#'//TRIM(chGEO)//'_Base#'//TRIM(chBAS)//'_Cycl#'//TRIM(chSCF) &
         //'_Clone#'//TRIM(IntToChar(1))//'.D'
    CALL Get(P,TrixName)
    chSCF=IntToChar(S%Current%I(1))
    DO AtA=1,NAtoms
      DO IX=1,3
        DO II=1,2
          !
          ! Move the atom.
          IF(II==1) THEN
            G%Clone(1)%Carts%D(IX,AtA)=GTmp%Carts%D(IX,AtA)+DDelta
          ELSEIF(II==2) THEN
            G%Clone(1)%Carts%D(IX,AtA)=GTmp%Carts%D(IX,AtA)-DDelta
          ENDIF
          !
          G%Clone(1)%Gradients%D = Zero
          CALL GeomArchive(cBAS,cGEO,N,O,B,G)
          ! vwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvw>>>
          ! Move back the atom.
          IF(II==1) THEN
            G%Clone(1)%Carts%D(IX,AtA)=GTmp%Carts%D(IX,AtA)
          ELSEIF(II==2) THEN
            G%Clone(1)%Carts%D(IX,AtA)=GTmp%Carts%D(IX,AtA)
          ENDIF
          ! vwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvwvw<<<
          CALL Invoke('TForce',N,S,M)
          ! Get the T Gradient
          HDFFileID=OpenHDF(N%HFile)
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
          CALL Get(G%Clone(1)%Gradients,'Gradients',Tag_O=chGEO)
          CALL CloseHDFGroup(HDF_CurrentID)
          CALL CloseHDF(HDFFileID)
          !copy grad
          FT(:,:,II)=0.0D0
          FT(:,:,II)=G%Clone(1)%Gradients%D(:,:)
          ! DO IA=1,NAtoms
          ! write(*,91919) 'G2(',1,',',IA,')=',G%Clone(1)%Gradients%D(1,IA),';'
          ! write(*,91919) 'G2(',2,',',IA,')=',G%Clone(1)%Gradients%D(2,IA),';'
          ! write(*,91919) 'G2(',3,',',IA,')=',G%Clone(1)%Gradients%D(3,IA),';'
          !91919        format(A,I2,A,I2,A,E26.15,A)
          ! ENDDO
        ENDDO
        IB=0
        IA=3*(AtA-1)+IX
        DO AtB=1,NAtoms
          IB=IB+1
          HT%D(IA,IB)=(FT(1,AtB,1)-FT(1,AtB,2))/(Two*DDelta)
          IB=IB+1
          HT%D(IA,IB)=(FT(2,AtB,1)-FT(2,AtB,2))/(Two*DDelta)
          IB=IB+1
          HT%D(IA,IB)=(FT(3,AtB,1)-FT(3,AtB,2))/(Two*DDelta)
        ENDDO
      ENDDO
    ENDDO
    !
    !
    G%Clone(1)%Carts%D=GTmp%Carts%D
    G%Clone(1)%Gradients%D=GTmp%Gradients%D
    !
    CALL GeomArchive(cBAS,cGEO,N,O,B,G)
    !
    CALL Print_DBL_RNK2(HT,'THessian',Unit_O=6)

    CALL Delete(HT)
    CALL Delete(GTmp)
    CALL Delete(P)
    CALL Delete(BSiz)
    CALL Delete(OffS)
  END SUBROUTINE NTHessian
  !===============================================================================
  !===============================================================================
  !===============================================================================
  !===============================================================================
  SUBROUTINE Multiply_DR2(N,A,B,C)
    INTEGER            :: I,J,K,N
    TYPE(DBL_RNK2)     :: A,B,C
    !
    DO I=1,N
      DO J=1,N
        C%D(I,J) = Zero
        DO K=1,N
          C%D(I,J) = C%D(I,J) + A%D(I,K)*B%D(K,J)
        ENDDO
      ENDDO
    ENDDO
    !
  END SUBROUTINE Multiply_DR2
  !
  FUNCTION Trace_DR2(N,A)
    INTEGER            :: I,N
    TYPE(DBL_RNK2)     :: A
    REAL(DOUBLE)       :: Trace_DR2
    !
    Trace_DR2 = Zero
    DO I=1,N
      Trace_DR2 = Trace_DR2 + A%D(I,I)
    ENDDO
    !
  END FUNCTION Trace_DR2
  !===============================================================================
  ! Numerically compute Lattice gradients
  !===============================================================================
  SUBROUTINE NLattForce0(C)
    TYPE(Controls)                   :: C
    INTEGER                          :: I1,I2,cSCF,cBAS,cGEO
    REAL(DOUBLE),PARAMETER           :: DDelta = 1.D-4
    REAL(DOUBLE)                     :: E_low,E_hig,Lat00
    TYPE(DBL_RNK2)                   :: LattF

    ! Write Out Delta
    WRITE(*,*) 'DDelta = ',DDelta
    ! Allocate
    CALL New(LattF,(/3,3/))
    LattF%D=Zero
    !
    cSCF = C%Stat%Current%I(1)
    cBAS = C%Stat%Current%I(2)
    cGEO = C%Stat%Current%I(3)
    ! Loop Over I1,I2
    DO I1=1,3
      DO I2=1,3
        Lat00 = C%Geos%Clone(1)%PBC%BoxShape%D(I1,I2)
        IF(C%Geos%Clone(1)%PBC%AutoW%I(I1)==1 .AND. C%Geos%Clone(1)%PBC%AutoW%I(I2)==1) THEN
          C%Geos%Clone(1)%PBC%BoxShape%D(I1,I2) = Lat00-DDelta
          CALL MkGeomPeriodic(C%Geos%Clone(1))
          CALL SinglePoints(C)
          !
          HDFFileID=OpenHDF(C%Nams%HFile)
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
          CALL Get(E_low,'Etot')
          CALL CloseHDF(HDFFileID)
          !
          C%Geos%Clone(1)%PBC%BoxShape%D(I1,I2) = Lat00+DDelta
          CALL MkGeomPeriodic(C%Geos%Clone(1))
          CALL SinglePoints(C)
          !
          HDFFileID=OpenHDF(C%Nams%HFile)
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
          CALL Get(E_hig,'Etot')
          CALL CloseHDF(HDFFileID)
          !
          LattF%D(I1,I2) =  (E_hig-E_low)/(Two*DDelta)
          !
          C%Geos%Clone(1)%PBC%BoxShape%D(I1,I2) = Lat00
          CALL MkGeomPeriodic(C%Geos%Clone(1))
        ENDIF
      ENDDO
    ENDDO
    !
    DO I1=1,3
      WRITE(*,99) (LattF%D(I1,I2),I2=1,3)
    END DO
    !
99  FORMAT(3(D23.14,1X))

  END SUBROUTINE NLattForce0
  !===============================================================================
  ! Numerically compute gradients
  !===============================================================================
  SUBROUTINE NLattForce1(C,I1,I2)
    TYPE(Controls)                   :: C
    INTEGER                          :: I,J,K,iSCF,I1,I2
    REAL(DOUBLE),PARAMETER           :: DDelta = 1.D-5
    REAL(DOUBLE)                     :: E_low,E_hig,Lat00
    REAL(DOUBLE)                     :: TracePFPdS,TracePdT,TracePdJ1,TracePdJ2
    REAL(DOUBLE)                     :: dNukE,TracePdK,dExc
    CHARACTER(LEN=DCL)               :: TrixName
    TYPE(DBL_RNK2)                   :: P,F,T1,T2,T3
    TYPE(BCSR)                       :: Temp
    TYPE(CellSet)                    :: CS
    ! Load global
    Lat00 = C%Geos%Clone(1)%PBC%BoxShape%D(I1,I2)
    CALL New(BSiz,C%Geos%Clone(1)%NAtms)
    CALL New(OffS,C%Geos%Clone(1)%NAtms)
    CALL New(C%Stat%Action,1)
    NAToms  = C%Geos%Clone(1)%NAtms
    MaxAtms = C%Sets%MxAts(1)
    MaxBlks = C%Sets%MxN0s(1)
    MaxNon0 = C%Sets%MxBlk(1)
    NBasF   = C%Sets%BSets(1,1)%NBasF
    BSiz%I  = C%Sets%BSiz(1,1)%I
    OffS%I  = C%Sets%OffS(1,1)%I
    MaxBlkSize=0
    DO I=1,C%Geos%Clone(1)%NAtms
      MaxBlkSize=MAX(MaxBlkSize,BSiz%I(I))
    ENDDO
    ! Write Out Delta
    WRITE(*,*) 'DDelta = ',DDelta
    WRITE(*,*) 'Lattice'
    WRITE(*,*) C%Geos%Clone(1)%PBC%BoxShape%D(1,1:3)
    WRITE(*,*) C%Geos%Clone(1)%PBC%BoxShape%D(2,1:3)
    WRITE(*,*) C%Geos%Clone(1)%PBC%BoxShape%D(3,1:3)
    WRITE(*,*)
    ! Intitialize
    CALL New(F ,(/NBasF,NBasF/))
    CALL New(P ,(/NBasF,NBasF/))
    CALL New(T1,(/NBasF,NBasF/))
    CALL New(T2,(/NBasF,NBasF/))
    CALL New(T3,(/NBasF,NBasF/))
    iSCF = C%Stat%Current%I(1)
    ! Get F
    HDFFileID=OpenHDF(C%Nams%HFile)
    HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
    TrixName=TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)//'_Geom#1_Base#1_Cycl#'//TRIM(IntToChar(iSCF))//'_Clone#1.F'
    CALL Get(Temp,TrixName)
    CALL SetEq(F,Temp)
    ! Get P
    TrixName=TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)//'_Geom#1_Base#1_Cycl#'//TRIM(IntToChar(iSCF))//'_Clone#1.D'
    CALL Get(Temp,TrixName)
    CALL SetEq(P,Temp)
    CALL CloseHDF(HDFFileID)
    ! Compute PFP
    CALL Multiply_DR2(NBasF,F,P,T1)
    CALL Multiply_DR2(NBasF,P,T1,F)
    !
    ! Calculate -2*Trace[P*F*P*dS]
    !
    C%Stat%Action%C(1)='OneElectronMatrices'
    C%Geos%Clone(1)%PBC%BoxShape%D(I1,I2) = Lat00-DDelta
    CALL MkGeomPeriodic(C%Geos%Clone(1))
    CALL GeomArchive(1,1,C%Nams,C%Opts,C%Sets,C%Geos)
    CALL Invoke('MakeS',C%Nams,C%Stat,C%MPIs)
    !
    HDFFileID=OpenHDF(C%Nams%HFile)
    HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
    TrixName=TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)//'_Geom#1_Base#1_Clone#1.S'
    CALL Get(Temp,TrixName)
    CALL SetEq(T1,Temp)
    CALL CloseHDF(HDFFileID)
    !
    C%Geos%Clone(1)%PBC%BoxShape%D(I1,I2) = Lat00+DDelta
    CALL MkGeomPeriodic(C%Geos%Clone(1))
    CALL GeomArchive(1,1,C%Nams,C%Opts,C%Sets,C%Geos)
    CALL Invoke('MakeS'   ,C%Nams,C%Stat,C%MPIs)
    !
    HDFFileID=OpenHDF(C%Nams%HFile)
    HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
    TrixName=TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)//'_Geom#1_Base#1_Clone#1.S'
    CALL Get(Temp,TrixName)
    CALL SetEq(T2,Temp)
    !
    T3%D = (T2%D-T1%D)/(Two*DDelta)
    !
    CALL Multiply_DR2(NBasF,F,T3,T1)
    TracePFPdS = -Two*Trace_DR2(NBasF,T1)
    ! WRITE(*,*) 'LatFrc_S'
    ! WRITE(*,*) '(',I1,',',I2,')=',TracePFPdS
    !
    ! Calculate 2*Trace[P*dT]
    !
    C%Geos%Clone(1)%PBC%BoxShape%D(I1,I2) = Lat00-DDelta
    CALL MkGeomPeriodic(C%Geos%Clone(1))
    CALL GeomArchive(1,1,C%Nams,C%Opts,C%Sets,C%Geos)
    CALL Invoke('MakeT',C%Nams,C%Stat,C%MPIs)
    !
    HDFFileID=OpenHDF(C%Nams%HFile)
    HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
    TrixName=TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)//'_Geom#1_Base#1_Clone#1.T'
    CALL Get(Temp,TrixName)
    CALL SetEq(T1,Temp)
    CALL CloseHDF(HDFFileID)
    !
    C%Geos%Clone(1)%PBC%BoxShape%D(I1,I2) = Lat00+DDelta
    CALL MkGeomPeriodic(C%Geos%Clone(1))
    CALL GeomArchive(1,1,C%Nams,C%Opts,C%Sets,C%Geos)
    CALL Invoke('MakeT',C%Nams,C%Stat,C%MPIs)
    !
    HDFFileID=OpenHDF(C%Nams%HFile)
    HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
    TrixName=TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)//'_Geom#1_Base#1_Clone#1.T'
    CALL Get(Temp,TrixName)
    CALL SetEq(T2,Temp)
    CALL CloseHDF(HDFFileID)
    !
    T3%D = (T2%D-T1%D)/(Two*DDelta)
    CALL Multiply_DR2(NBasF,P,T3,T1)
    TracePdT = Two*Trace_DR2(NBasF,T1)
    ! WRITE(*,*) 'LatFrc_T '
    ! WRITE(*,*) '(',I1,',',I2,')=',TracePdT
    !
    ! Calculate d(Exc)
    !
    C%Geos%Clone(1)%PBC%BoxShape%D(I1,I2) = Lat00-DDelta
    CALL MkGeomPeriodic(C%Geos%Clone(1))
    CALL GeomArchive(1,1,C%Nams,C%Opts,C%Sets,C%Geos)
    !
    CALL Invoke('MakeRho'     ,C%Nams,C%Stat,C%MPIs)
    CALL Invoke('HiCu'     ,C%Nams,C%Stat,C%MPIs)
    HDFFileID=OpenHDF(C%Nams%HFile)
    HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
    CALL Get(E_low,'Exc')
    CALL CloseHDF(HDFFileID)
    !
    C%Geos%Clone(1)%PBC%BoxShape%D(I1,I2) = Lat00+DDelta
    CALL MkGeomPeriodic(C%Geos%Clone(1))
    CALL GeomArchive(1,1,C%Nams,C%Opts,C%Sets,C%Geos)
    !
    CALL Invoke('MakeRho'     ,C%Nams,C%Stat,C%MPIs)
    CALL Invoke('HiCu'     ,C%Nams,C%Stat,C%MPIs)
    HDFFileID=OpenHDF(C%Nams%HFile)
    HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
    CALL Get(E_hig,'Exc')
    CALL CloseHDF(HDFFileID)
    !
    dExc = (E_hig-E_low)/(Two*DDelta)

    WRITE(*,*) 'LatFrc_XC'
    WRITE(*,*) '(',I1,',',I2,')=',dExc
    WRITE(*,*)

    C%Geos%Clone(1)%PBC%BoxShape%D(I1,I2) = Lat00
    CALL MkGeomPeriodic(C%Geos%Clone(1))
    CALL GeomArchive(1,1,C%Nams,C%Opts,C%Sets,C%Geos)
    CALL Invoke('MakePFFT' ,C%Nams,C%Stat,C%MPIs)
    CALL Invoke('MakeRho'     ,C%Nams,C%Stat,C%MPIs)
    !
    ! Calculate 2*(Trace[P*dJ]+NukeE)
    !
    C%Geos%Clone(1)%PBC%BoxShape%D(I1,I2) = Lat00-DDelta
    CALL MkGeomPeriodic(C%Geos%Clone(1))
    CALL GeomArchive(1,1,C%Nams,C%Opts,C%Sets,C%Geos)
    CALL Invoke('MakePFFT' ,C%Nams,C%Stat,C%MPIs)! CALL Invoke('MakeRho'     ,C%Nams,C%Stat,C%MPIs)
    CALL Invoke('MakeRho'     ,C%Nams,C%Stat,C%MPIs)
    CALL Invoke('QCTC'     ,C%Nams,C%Stat,C%MPIs)
    HDFFileID=OpenHDF(C%Nams%HFile)
    HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
    CALL Get(E_low,'E_NuclearTotal')
    TrixName=TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)//'_Geom#1_Base#1_Cycl#'//TRIM(IntToChar(iSCF))//'_Clone#1.J'
    CALL Get(Temp,TrixName)
    CALL SetEq(T1,Temp)
    CALL CloseHDF(HDFFileID)
    !
    CALL Multiply_DR2(NBasF,P,T1,T3)
    E_low = Trace_DR2(NBasF,T3)+E_low
    !
    !
    C%Geos%Clone(1)%PBC%BoxShape%D(I1,I2) = Lat00+DDelta
    CALL MkGeomPeriodic(C%Geos%Clone(1))
    CALL GeomArchive(1,1,C%Nams,C%Opts,C%Sets,C%Geos)
    CALL Invoke('MakePFFT' ,C%Nams,C%Stat,C%MPIs)
    CALL Invoke('MakeRho'     ,C%Nams,C%Stat,C%MPIs)
    CALL Invoke('QCTC'     ,C%Nams,C%Stat,C%MPIs)
    HDFFileID=OpenHDF(C%Nams%HFile)
    HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
    CALL Get(E_hig,'E_NuclearTotal')
    TrixName=TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)//'_Geom#1_Base#1_Cycl#'//TRIM(IntToChar(iSCF))//'_Clone#1.J'
    CALL Get(Temp,TrixName)
    CALL SetEq(T2,Temp)
    CALL CloseHDF(HDFFileID)
    !
    CALL Multiply_DR2(NBasF,P,T2,T3)
    E_hig = Trace_DR2(NBasF,T3)+E_hig
    !
    !
    !
    dNukE = (E_hig-E_low)/(Two*DDelta)
    WRITE(*,*) 'LatFrc_QCTC'
    WRITE(*,*) '(',I1,',',I2,')=',dNukE
    WRITE(*,*)
    !
    !
    C%Geos%Clone(1)%PBC%BoxShape%D(I1,I2) = Lat00
    CALL MkGeomPeriodic(C%Geos%Clone(1))
    CALL GeomArchive(1,1,C%Nams,C%Opts,C%Sets,C%Geos)
    CALL Invoke('MakePFFT' ,C%Nams,C%Stat,C%MPIs)
    CALL Invoke('MakeRho'     ,C%Nams,C%Stat,C%MPIs)
    !
    ! Delete
    CALL Delete(BSiz)
    CALL Delete(OffS)
    CALL Delete(C%Stat%Action)
    CALL Delete(F)
    CALL Delete(P)
    CALL Delete(T1)
    CALL Delete(T2)
    CALL Delete(T3)
    !
  END SUBROUTINE NLattForce1
  !--------------------------------------------------------------
  ! Determin iREMOVE and MinMDGeo
  !--------------------------------------------------------------
  SUBROUTINE CalculateMDGeo(D,iREMOVE,MinMDGeo)
    TYPE(Dynamics)  :: D
    INTEGER         :: iREMOVE,MinMDGeo

    iREMOVE = D%MDMaxSteps+1
    SELECT CASE(D%MDGeuss)
    CASE("DMLinear","FMVerlet0")
      iREMOVE  = 2
      MinMDGeo = 2
    CASE("DMTRBO")
      iREMOVE  = 3
      MinMDGeo = 3
    CASE("DMTRBO_Damp_dt3")
      iREMOVE  = 4
      MinMDGeo = 4
    CASE("DMTRBO_Damp_dt5")
      iREMOVE  = 5
      MinMDGeo = 5
    CASE("DMTRBO_Damp_dt7")
      iREMOVE  = 6
      MinMDGeo = 6
    CASE("DMTRBO_Damp_dt9")
      iREMOVE  = 7
      MinMDGeo = 7
    CASE("DMTRBO_Damp_dt11")
      iREMOVE  = 8
      MinMDGeo = 8
    CASE("DMSymplectic")
      iREMOVE  = 6
      MinMDGeo = 6
    CASE("FMVerlet1")
      iREMOVE  = 4
      MinMDGeo = 4
    CASE('DMProj0')
      iREMOVE  = 1
      MinMDGeo = 1
    CASE('DMProj1')
      iREMOVE  = 2
      MinMDGeo = 2
    CASE('DMProj2')
      iREMOVE  = 3
      MinMDGeo = 3
    CASE('DMProj3')
      iREMOVE  = 4
      MinMDGeo = 4
    CASE('DMProj4')
      iREMOVE  = 5
      MinMDGeo = 5
    CASE('DMDGeuss')
      iREMOVE  = 1
      MinMDGeo = 1
    CASE DEFAULT
      CALL Halt("[SCFs.CalculateMDGeo] unknown MDGeuss "//TRIM(D%MDGeuss))
    END SELECT
  END SUBROUTINE CalculateMDGeo

END MODULE SCFs
