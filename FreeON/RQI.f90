#define KWERQY
!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has

!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
MODULE RayleighQuotientIteration
  USE InOut
  USE Macros
  USE MemMan
  USE Overlay
  USE Indexing
  USE AtomPairs
  USE PunchHDF
  USE PrettyPrint
  USE GlobalCharacters
  USE GlobalScalars
  USE ControlStructures
  USE OptionKeys
  USE Response
  !  USE McMurchie
  USE MondoLogger
  !!  USE CWrappers
  IMPLICIT NONE
  TYPE(TIME)   :: TimeTotal,TimeONX,TimeQCTC,TimeBCSR, &
       TimeGRAD,TimeSRCH,TimeLON2,TimeRECO,TimeSPLT
  LOGICAL      :: DoKoopmans=.FALSE.
  REAL(DOUBLE) :: OneERPAScaling

  CHARACTER(LEN=DEFAULT_CHR_LEN) :: KoopmansName
  ! STO-2G/RPA/p20
  REAL(DOUBLE) :: E_RPA=0.3476049424338692D+01 !Tight
  !   REAL(DOUBLE) :: E_RPA=0.3476176935406602D+01 !VeryTight
  !                Ev = 0.3476176141052899D+01 !TwoE*1D2


  ! 6-31G**/RPA/p5
  ! QUIRQI                 ::  Ev = 0.4958423523032192D+01, ErrRel = 0.55D+03
  !  REAL(DOUBLE) :: E_RPA=0.3235800692514660D+01

  ! 3-21G/RPA/p5
  ! QUIRQI                 ::  Ev = 0.4958423523032192D+01, ErrRel = 0.55D+03
  !    REAL(DOUBLE) :: E_RPA=0.3399129180718000D+01


  !   REAL(DOUBLE) ::  E_RPA=0.406088110418702D0*27.21139613182D0

  TYPE PQName
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: AO
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: Ortho
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: PosP
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: MomQ
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: PosP_Old
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: MomQ_Old
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: Act2
  END TYPE PQName
  !
  TYPE QUIRQIKontrol
     CHARACTER(LEN=2) :: chBas
     TYPE(TOLS)   :: Current
     TYPE(TOLS)   :: Ultimate
     INTEGER,DIMENSION(3) :: Status
     LOGICAL      :: CholFact
     REAL(DOUBLE) :: Threshold
     REAL(DOUBLE) :: Konvergence
     REAL(DOUBLE) :: UltimateThreshold
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: P
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: J
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: K
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: F
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: Z
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: ZT
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: ASCII_Xk
     TYPE(PQName)                   :: Xk
     TYPE(PQName)                   :: Gk
     TYPE(PQName)                   :: Pk
     TYPE(PQName)                   :: LXk
     TYPE(PQName)                   :: LPk
  END TYPE QUIRQIKontrol
  !
CONTAINS

  SUBROUTINE TDSCF(C)
    IMPLICIT NONE
    TYPE(Controls)     :: C
    TYPE(QUIRQIKontrol):: QN
    IF(.NOT. C%POpt%Resp%TD_SCF) RETURN
    CALL SetFrontEndMacros(C%Nams,C%Geos,C%Sets,C%Opts,QN)
    CALL xQUIRQI(C,C%Nams,C%Opts,C%Stat,C%MPIs,C%Geos,C%Sets%BSets(1,C%Sets%NBSets),QN)
    !    WRITE(*,*)'=================== QQQQQQQQQQ ==============='
    !    CALL QUIRQI(NBasF,4,C%Nams,C%Opts,C%Stat,C%MPIs,C%Geos,C%Sets%BSets(1,C%Sets%NBSets))
    !    WRITE(*,*)'=================== QQQQQQQQQQ ==============='

    !    CALL sQUIRQI(C%Nams,C%Opts,C%Stat,C%MPIs,C%Geos,C%Sets%BSets(1,C%Sets%NBSets),QN)
    !    CALL RQI(NBasF,4,C%Nams,C%Opts,C%Stat,C%MPIs,C%Geos,C%Sets%BSets(1,C%Sets%NBSets))
  END SUBROUTINE TDSCF
  !===============================================================================
  !   Subroutine to load the global macro parameters used mostly in the backend
  !   into the front end, so that we can do kluge work. Assumes that we are
  !   at the end of our basis set list, first clone etc. NOT APPROPRIATE for
  !   general work.
  !===============================================================================

  SUBROUTINE SetFrontEndMacros(N,G,B,O,QN)
    TYPE(FileNames)    :: N
    TYPE(BasisSets)  :: B
    TYPE(Geometries) :: G
    TYPE(Options)    :: O
    TYPE(QUIRQIKontrol):: QN
    INTEGER          :: II
    !

    PrintFlags%Mat=DEBUG_MATRICES
    !   PrintFlags%Key=DEBUG_MAXIMUM
    !   PrintFlags%Chk=DEBUG_CHKSUMS

    MyClone=1
    NBasF=B%BSets(1,B%NBSets)%NBasF
    NAToms=G%Clone(1)%NAtms
    MaxAtms=B%MxAts(B%NBSets)
    MaxBlks=B%MxBlk(B%NBSets)
    MaxNon0=B%MxN0s(B%NBSets)
    CALL New(BSiz,NAtoms)
    CALL New(OffS,NAtoms)
    BSiz%I=B%BSiz(1,B%NBSets)%I
    OffS%I=B%OffS(1,B%NBSets)%I
    MaxBlkSize=0
    DO II=1,G%Clone(1)%NAtms
       MaxBlkSize=MAX(MaxBlkSize,BSiz%I(II))
    ENDDO
!!$
    QN%chBas=IntToChar(B%NBSets)
    HDFFileID=OpenHDF(N%HFile)
    HDF_CurrentID=HDFFileID
    HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
    CALL Get(QN%Current,QN%chBas)
    CALL Get(QN%Ultimate,QN%chBas)
    CALL CloseHDFGroup(HDF_CurrentID)
    CALL CloseHDF(HDFFileID)

    QN%Threshold=O%Thresholds(B%NBSets)%Trix*1D2
    QN%Konvergence=O%Thresholds(B%NBSets)%ETol*1D2

    !
  END SUBROUTINE SetFrontEndMacros
  ! ===================================================================================
  ! QUASI-INDEPENDENT RAYLEIGHT QUOTIENT ITERATION (QUIRQI)
  ! EOM FORMULATION OF TD-SCF IS NEARLY HERMETIAN, COUPLED ONLY THROUGH
  ! THE ORTHOGONALITY RELATIONSHIP BETWEEN EOM VARIABLES P AND Q (NOT PROJECTORS,
  ! BUT COMPONENTS OF THE TRANSITION DENSITY MATRIX).  NEW IDEA IS TO SOLVE THE
  ! TD-SCF IN EOM REPRESENTATION ASSUMING P AND Q ARE NEARLY INDEPENDENT.  WEAK
  ! COUPLING GIVES RISE TO NONLINEARITY, SOLVED BY USE OF DUAL TRACK POLAK-RIBIERE
  ! NLCG AND EXACT, ANALYTIC LINE SEARCH.
  ! ===================================================================================

  SUBROUTINE PolarizationGuess(C,qn)
    IMPLICIT NONE
    TYPE(Controls)                 :: C
    TYPE(QUIRQIKontrol)            :: QN
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: Guess
    REAL(DOUBLE)                   :: Norm
    TYPE(BCSR)                     :: sP,sT1,sT2,sT3
    INTEGER,DIMENSION(3)           :: Status
    !
    C%POpt%Resp%StcAlpha=.TRUE.
    C%POpt%Resp%StcBeta=.FALSE.
    C%POpt%Resp%StcGamma=.FALSE.
    !
    C%POpt%Resp%AlphaAxis=.FALSE.
    IF(TRIM(C%Opts%RQIGuess)=='PolarX')THEN
       C%POpt%Resp%AlphaAxis(1)=.TRUE.
    ELSEIF(TRIM(C%Opts%RQIGuess)=='PolarY')THEN
       C%POpt%Resp%AlphaAxis(2)=.TRUE.
    ELSEIF(TRIM(C%Opts%RQIGuess)=='PolarZ')THEN
       C%POpt%Resp%AlphaAxis(3)=.TRUE.
    ELSE
       CALL MondoHalt(PRSE_ERROR,' Bad logic in CPSCF guess ')
    ENDIF
    !
    HDFFileID=OpenHDF(C%Nams%HFile)
    HDF_CurrentID=HDFFileID
    HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
    QN%Current%TwoE=1D-6
    QN%Current%Trix=1D-5
    QN%Current%Dist=1D-9
    CALL Put(QN%Current,QN%chBas)
    CALL CloseHDFGroup(HDF_CurrentID)
    CALL CloseHDF(HDFFileID)
    !
    ! Archive SCF state
    Status(:)=C%Stat%Current%I
    ! Solve the CPSCF
    CALL CPSCF(C)
    !
    HDFFileID=OpenHDF(C%Nams%HFile)
    HDF_CurrentID=HDFFileID
    HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
    CALL Put(QN%Ultimate,QN%chBas)
    CALL CloseHDFGroup(HDF_CurrentID)
    CALL CloseHDF(HDFFileID)
    ! Rename the response
    IF(TRIM(C%Opts%RQIGuess)=='PolarX')THEN
       Guess=TrixFile("OrthoDPrimeX",PWD_O=C%Nams%M_SCRATCH,Name_O=C%Nams%SCF_NAME,Stats_O=C%Stat%Current%I,OffSet_O=0)
    ELSEIF(TRIM(C%Opts%RQIGuess)=='PolarY')THEN
       Guess=TrixFile("OrthoDPrimeY",PWD_O=C%Nams%M_SCRATCH,Name_O=C%Nams%SCF_NAME,Stats_O=C%Stat%Current%I,OffSet_O=0)
    ELSEIF(TRIM(C%Opts%RQIGuess)=='PolarZ')THEN
       Guess=TrixFile("OrthoDPrimeZ",PWD_O=C%Nams%M_SCRATCH,Name_O=C%Nams%SCF_NAME,Stats_O=C%Stat%Current%I,OffSet_O=0)
    ENDIF
    CALL FileCopy(Guess,QN%Xk%Ortho)
    ! Reset to reflect SCF, not CPSCF reality.
    C%Stat%Current%I=Status
    ! Normalize
    CALL Get(sT1,QN%Xk%Ortho)
    Norm=SQRT(Two*ABS(Trace(sT1,sT1)))
    CALL Multiply(sT1,One/Norm)
    CALL Put(sT1,QN%Xk%Ortho)
    CALL Delete(sT1)
    ! Split, replicate and recombine
    CALL SplitPQ(QN,QN%Xk,.TRUE.)
    CALL Get(sT1,QN%Xk%PosP)
    CALL Multiply(sT1,-One)
    CALL Put(sT1,QN%Xk%MomQ)
    CALL Recomb(QN,QN%Xk)
    CALL Delete(sT1)
    !
  END SUBROUTINE PolarizationGuess

  SUBROUTINE xQUIRQI(C,Nam,O,S,MPI,G,B,QN)
    TYPE(Controls)     :: C
    TYPE(QUIRQIKontrol):: QN
    INTEGER            :: I,J,K
    TYPE(FileNames)    :: Nam,RQINams
    TYPE(State)        :: S,RQIStat
    TYPE(Parallel)     :: MPI
    TYPE(Options)      :: O
    TYPE(Geometries)   :: G
    TYPE(BSET)         :: B
    LOGICAL            :: DoTDA
    REAL(DOUBLE)       :: Ek_P,Ek_Q,Ek,ErrAbs,EkOld0,EkOld1,ErrRel,Delta,BetaP,BetaQ, &
         XkPcntP,XkPcntQ,XkPcnt,PkPcnt,LXkPcnt,LPkPcnt,MaxGP,MaxGQ,KTau,XTau,GTau,SqrRel,MinRel
    TYPE(BCSR)         :: sX
    TYPE(DBL_RNK2)     :: dX
    CHARACTER(LEN=DCL) :: Iteration,Statistics1,Statistics2,Statistics3
    REAL(DOUBLE),DIMENSION(0:O%MaxRQI) :: EConverge

!    TYPE(DBL_VECT)    :: EigenEs
!    TYPE(DBL_RNK2)    :: X,P,Q,T1,T2,T3
    !
    DoTDA=.TRUE.
    DoTDA=.FALSE.
    !
    CALL SetQName(Nam,S,QN,RQIStat)
    !
    IF(TRIM(O%RQIGuess)=='Koopmans')THEN
       QN%Threshold=1D-8
       CALL Koopmans(QN)
    ELSE
       CALL PolarizationGuess(C,qn)
    ENDIF
    !
    CALL Elapsed_TIME(TimeTotal,Init_O='Init')
    CALL Elapsed_TIME(TimeONX,Init_O='Init')
    CALL Elapsed_TIME(TimeQCTC,Init_O='Init')
    CALL Elapsed_TIME(TimeGRAD,Init_O='Init')
    CALL Elapsed_TIME(TimeSRCH,Init_O='Init')
    CALL Elapsed_TIME(TimeLON2,Init_O='Init')
    CALL Elapsed_TIME(TimeRECO,Init_O='Init')
    CALL Elapsed_TIME(TimeSPLT,Init_O='Init')
    CALL Elapsed_TIME(TimeTotal,Init_O='Start')
    !
    EkOld0=1D2
    EkOld1=1D2
    XkPcnt=Zero
    PkPcnt=Zero
    XTau=1D-4
    KTau=1D-5
    GTau=1D-4
    !
    CALL LOn2AO(Nam,MPI,S,QN,RQIStat,QN%Xk,QN%LXk)
    CALL SNihilate(QN,QN%LXk,DoTDA)
    CALL OrthoPrune(QN,QN%LXk%Ortho,.FALSE.,XTau,XkPcnt)
    CALL SplitPQ(QN,QN%LXk,PQThreshold_O=XTau)
    CALL OrthoPrune(QN,QN%LXk%PosP ,.FALSE.,XTau,XkPcnt)
    CALL OrthoPrune(QN,QN%LXk%MomQ ,.FALSE.,XTau,XkPcnt)
    !
    DO K=0,O%MaxRQI-1
       ! - - - - - - - - - - - - - - - - - - - - -
       CALL GradGrad(QN,Ek_P,Ek_Q,MaxGP,MaxGQ)
       CALL Recomb(QN,QN%Gk)
       CALL SNihilate(QN,QN%Gk,DoTDA)
       CALL SplitPQ(QN,QN%Gk,Switch_O=.TRUE.,PQThreshold_O=GTau)
       !
       Ek=Half*(Ek_P+Ek_Q)
       EConverge(K)=Ek
       IF(MOD(K,2)==0)THEN
          EkOld0=MIN(EkOld0,Ek)
          ErrAbs=EkOld1-Ek
       ELSE
          EkOld1=MIN(EkOld1,Ek)
          ErrAbs=EkOld0-Ek
       ENDIF
       !
       ErrRel=ErrAbs/Ek
       IF(ErrRel<5D-4.AND.K>15)EXIT
       !
       MinRel=MIN(MinRel,SQRT(ABS(ErrRel)))

       !       XTau=MAX(MIN(XTau,MinRel),QN%Threshold)
       GTau=MIN(GTau,1D-1*MAX(MaxGP,MaxGQ))
       !
       CALL Elapsed_TIME(TimeTOTAL,Init_O='Accum')
       CALL QUIRQIPrint(QN,K,Ek,Ek_P,Ek_Q,ErrRel,QN%Current%TwoE,XTau,GTau,XkPcnt,PkPcnt,EConverge,Nam,.FALSE.)
       XkPcnt=Zero
       PkPcnt=Zero
       !
       IF(K==0)THEN
          CALL Recomb(QN,QN%Gk)
          CALL FileCopy(QN%Gk%Ortho,QN%Pk%Ortho)
       ENDIF
       ! - - - - - - - - - - - - - - - - - - - - -
       ! Polak Ribierre step to generate P
       ! - - - - - - - - - - - - - - - - - - - - -
       CALL PRStep(QN,K)
       !
       IF(K>0)THEN
          CALL Recomb(QN,QN%Pk)
       ENDIF
       CALL OrthoPrune(QN,QN%Pk%PosP ,.FALSE.,GTau,PkPcnt)
       CALL OrthoPrune(QN,QN%Pk%MomQ ,.FALSE.,GTau,PkPcnt)
       CALL OrthoPrune(QN,QN%Pk%Ortho,.FALSE.,GTau,PkPcnt)
       !
       CALL SNihilate(QN,QN%Pk,DoTDA)

       CALL LOn2AO(Nam,MPI,S,QN,RQIStat,QN%Pk,QN%LPk,KTau_O=KTau)
       CALL SNihilate(QN,QN%LPk,DoTDA)
       CALL OrthoPrune(QN,QN%LPk%Ortho,.FALSE.,GTau,PkPcnt)
       CALL SplitPQ(QN,QN%LPk,PQThreshold_O=GTau)
       CALL OrthoPrune(QN,QN%LPk%PosP,.FALSE.,GTau,PkPcnt)
       CALL OrthoPrune(QN,QN%LPk%MomQ,.FALSE.,GTau,PkPcnt)
       ! - - - - - - - - - - - - - - - - - - - - -
       ! Line search
       ! - - - - - - - - - - - - - - - - - - - - -
       IF(.NOT. &
            LinSrch(QN)    )EXIT
       CALL Recomb(QN,QN%Xk)
       CALL OrthoPrune(QN,QN%Xk%PosP ,.FALSE.,XTau,XkPcnt)
       CALL OrthoPrune(QN,QN%Xk%MomQ ,.FALSE.,XTau,XkPcnt)
       CALL OrthoPrune(QN,QN%Xk%Ortho,.FALSE.,XTau,XkPcnt)
       !
       CALL SNihilate(QN,QN%Xk,DoTDA)

       CALL LOn2AO(Nam,MPI,S,QN,RQIStat,QN%Xk,QN%LXk,KTau_O=KTau)
       CALL SNihilate(QN,QN%LXk,DoTDA)
       CALL OrthoPrune(QN,QN%LXk%Ortho,.FALSE.,XTau,XkPcnt)
       CALL SplitPQ(QN,QN%LXk,PQThreshold_O=XTau)
       CALL OrthoPrune(QN,QN%LXk%PosP ,.FALSE.,XTau,XkPcnt)
       CALL OrthoPrune(QN,QN%LXk%MomQ ,.FALSE.,XTau,XkPcnt)
       !
    ENDDO
    !
    EConverge(K)=E_RPA
    CALL QUIRQIPrint(QN,K,Ek,Ek_P,Ek_Q,ErrRel,QN%Current%TwoE,XTau,GTau,XkPcnt,PkPcnt,EConverge,Nam,.TRUE.)
    CALL Get(sX,QN%Xk%Ortho)
    CALL SetEq(dX,sX)

!    CALL OpenASCII(QN%ASCII_Xk,77)
!    WRITE(77,*)' Atoms = ',NAtoms
!    WRITE(77,*)' Basis Functions = ',NBasF
!    Close(77)
!    CALL PPrint(BSiz,'Block Size ',FileName_O=QN%ASCII_Xk,Unit_O=77)
!    CALL PPrint(OffS,'Block Offset',FileName_O=QN%ASCII_Xk,Unit_O=77)
    CALL OpenASCII(QN%ASCII_Xk,77)
    DO I=1,NBasF
       DO J=1,NBasF
          WRITE(77,*)I,J,dX%D(I,J)
       ENDDO
    ENDDO
    CLOSE(Unit=77)


    !
  END SUBROUTINE xQUIRQI


  SUBROUTINE QUIRQIPrint(QN,K,Ek,Ek_P,Ek_Q,ErrRel,ETau,XTau,GTau,XkPcnt,PkPcnt,EConverge,Nam,FinalTally_O)
    TYPE(QUIRQIKontrol):: QN
    INTEGER :: I,K
    REAL(DOUBLE)       :: Ek_P,Ek_Q,Ek,ErrAbs,EkOld,ErrRel,XTau,GTau,ETau, &
         XkPcntP,XkPcntQ,PkPcntP,PkPcntQ,XkPcnt,PkPcnt,MaxGP,MaxGQ,thresh
    CHARACTER(LEN=DCL) :: Iteration,Statistics1,Statistics2,Statistics3
    REAL(DOUBLE),DIMENSION(:) :: EConverge
    TYPE(FileNames)    :: Nam
    LOGICAL, OPTIONAL  :: FinalTally_O
    CHARACTER(LEN=DCL) :: SCF_NAME
    !
    Iteration="RPA(1,"//TRIM(IntToChar(K))//')'
    Statistics1="Ev = "//TRIM(DblToMedmChar(Ek*27.21139613182D0))   &
         //", Err = "//TRIM(DblToShrtChar(ErrRel))                  &
         //", Torq = "//TRIM(DblToShrtChar(ABS(Ek_Q-Ek_P)/ABS(Ek)))
    !
    !       CALL MondoLog(DEBUG_NONE, "QUIRQI",Statistics1,Iteration)
    !
    Statistics2="tK = "//TRIM(DblToMedmChar(TimeONX%Wall))  &
         //", tQ = "//TRIM(DblToMedmChar(TimeQCTC%Wall)) &
         //", tM = "//TRIM(DblToMedmChar(TimeBCSR%Wall))
    !       CALL MondoLog(DEBUG_NONE, "QUIRQI",Statistics2)
    Statistics3=" Ev = "//TRIM(DblToChar(Ek*27.21139613182D0))   &
         //", ErrRel = "//TRIM(DblToShrtChar(ErrRel)) &
         //", TrueEr = "//TRIM(DblToShrtChar(         &
         ABS(E_RPA -Ek*27.21139613182D0)/E_RPA))  &
         //", ETau = "//TRIM(DblToShrtChar(ETau)) &
         //", XTau = "//TRIM(DblToShrtChar(XTau)) &
         //", GTau = "//TRIM(DblToShrtChar(GTau)) &
         //", %X = "//TRIM(DblToShrtChar(XkPcnt)) &
         //", %P = "//TRIM(DblToShrtChar(PkPcnt))

    CALL MondoLog(DEBUG_NONE, "QUIRQI",Statistics3,Iteration)

    SCF_NAME=Nam%SCF_NAME
    !       SCF_NAME='p20_X_FULL_RPA_PRMOD'


    CALL OpenASCII( TRIM(Nam%M_PWD)//'/'//TRIM(SCF_NAME)//'_XkPcnt_'//TRIM(DblToShrtChar(QN%Threshold)),77)
    WRITE(77,*)K,XkPcnt
    CLOSE(Unit=77)



    IF(PRESENT(FinalTally_O))THEN
       IF(FinalTally_O)THEN
          CALL OpenASCII( TRIM(Nam%M_PWD)//'/'//TRIM(SCF_NAME)//'_TD-SCF_Convergence_' &
               //TRIM(DblToShrtChar(QN%Threshold)),77)
          DO I=1,K-1
             WRITE(77,*)I,LOG10((MAX(1D-10,ABS(EConverge(K)-EConverge(I))/EConverge(K)))),EConverge(I)
          ENDDO
          CLOSE(Unit=77)
          !
          CALL OpenASCII( TRIM(Nam%M_PWD)//'TD-SCF_ONX',77)
          WRITE(77,*)NAtoms,TimeONX%Wall
          CLOSE(Unit=77)
          CALL OpenASCII( TRIM(Nam%M_PWD)//'TD-SCF_QCTC',77)
          WRITE(77,*)NAToms,TimeQCTC%Wall
          !
          CLOSE(Unit=77)
          CALL OpenASCII( TRIM(Nam%M_PWD)//'TD-SCF_GRAD',77)
          WRITE(77,*)NAtoms,TimeGRAD%Wall
          CLOSE(Unit=77)
          CALL OpenASCII( TRIM(Nam%M_PWD)//'TD-SCF_SRCH',77)
          WRITE(77,*)NAtoms,TimeSRCH%Wall
          CLOSE(Unit=77)
          CALL OpenASCII( TRIM(Nam%M_PWD)//'TD-SCF_LON2',77)
          WRITE(77,*)NAtoms,TimeLON2%Wall
          CLOSE(Unit=77)
          CALL OpenASCII( TRIM(Nam%M_PWD)//'TD-SCF_RECO',77)
          WRITE(77,*)NAtoms,TimeRECO%Wall
          CLOSE(Unit=77)
          CALL OpenASCII( TRIM(Nam%M_PWD)//'TD-SCF_SPLT',77)
          WRITE(77,*)NAtoms,TimeSPLT%Wall
          CLOSE(Unit=77)
          !
          CALL OpenASCII( TRIM(Nam%M_PWD)//'TD-SCF_BCSR',77)
          WRITE(77,*)NAtoms,TimeGRAD%Wall+TimeSRCH%Wall+ &
               TimeLON2%Wall+TimeRECO%Wall+TimeSPLT%Wall
          CLOSE(Unit=77)
          !
          CALL OpenASCII( TRIM(Nam%M_PWD)//'TD-SCF_ENERGIES',77)
          WRITE(77,*)NAtoms,Ek
          CLOSE(Unit=77)
          !
       ENDIF
    ENDIF

  END SUBROUTINE QUIRQIPrint


  SUBROUTINE PrunePH(QN,Xk,NormQ,Threshold,XkPcnt)
    REAL(DOUBLE)          :: Norm,Threshold
    CHARACTER(LEN=*)      :: Xk
    TYPE(QUIRQIKontrol)   :: QN
    LOGICAL               :: NormQ
    REAL(DOUBLE)          :: XkPcnt
    TYPE(BCSR)            :: sXk,sP,sQ,sT1,sT2
    INTEGER               :: Kount
    !-----------------------------------------------------------------------------

    CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
    !
    CALL New(sXk) ! Kluge, shouldn't have to do this!!
    Kount=1
111 CALL Get(sXk,Xk)

    CALL Get(sP,QN%P)
    CALL Multiply(sP,Two)
    CALL SetEq(sQ,sP)
    CALL Multiply(sQ,-One)
    CALL Add(sQ,Two)
    !
    CALL Multiply(sP,sXk,sT1)
    CALL Multiply(sT1,sQ,sT2)
    !
    CALL Filter(sXk,sT2,Tol_O=Threshold)
    !
    CALL Delete(sP)
    CALL Delete(sQ)
    CALL Delete(sT1)
    CALL Delete(sT2)
    !
!!$    IF(sXk%NNon0==0)THEN
!!$       WRITE(*,*)' NNon0 = ',sXk%NNon0
!!$       Kount=Kount+1
!!$       Threshold=5D-1*Threshold
!!$       WRITE(*,*)TRIM(Xk)
!!$       WRITE(*,*)'NNon0 = ',sXk%NNon0
!!$       IF(Kount>8) &
!!$            CALL Halt(' Logical error in RQI: OrthoPrune ')
!!$       GOTO 111
!!$    ENDIF
    !
    XkPcnt=1D2*DBLE(sXk%NNon0)/DBLE(NBasF**2)
    !    WRITE(*,*)Threshold,XkPcnt
    !
    CALL Put(sXk,Xk)
    ! Tidy up
    CALL Delete(sXk)
    ! All done
    CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
    !
  END SUBROUTINE PrunePH

  !
  SUBROUTINE OrthoPrune(QN,Xk,NormQ,Threshold,XkPcnt)
    REAL(DOUBLE)          :: Norm,Threshold
    CHARACTER(LEN=*)      :: Xk
    TYPE(QUIRQIKontrol)   :: QN
    LOGICAL               :: NormQ
    REAL(DOUBLE)          :: XkPcnt
    TYPE(BCSR)            :: sXk,sT
    INTEGER               :: Kount
    !-----------------------------------------------------------------------------
    CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
    !
    CALL New(sXk) ! Kluge, shouldn't have to do this!!
    Kount=1
111 CALL Get(sT,Xk)
    !
    CALL Filter(sXk,sT,Tol_O=Threshold)
    CALL Delete(sT)
    !
    IF(sXk%NNon0==0)THEN
       Kount=Kount+1
       Threshold=75D-2*Threshold
       WRITE(*,*)TRIM(Xk)
       WRITE(*,*)'NNon0 = ',sXk%NNon0
       IF(Kount>8) &
            CALL Halt(' Logical error in RQI: OrthoPrune ')
       GOTO 111
    ENDIF
    !
    XkPcnt=MAX(XkPcnt,1D2*DBLE(sXk%NNon0)/DBLE(NBasF**2))
    !
    CALL Put(sXk,Xk)
    ! Tidy up
    CALL Delete(sXk)
    ! All done
    CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
  END SUBROUTINE OrthoPrune

  SUBROUTINE EigenResolution(QN,Xk,LXk,Delta)
    TYPE(FileNames)               :: Nam
    TYPE(QUIRQIKontrol)           :: QN
    TYPE(PQName)                  :: Xk,LXk
    REAL(DOUBLE)                  :: EkRaw,EkFiltered,Delta,dEk,Ek,XkPcnt,Norm
    TYPE(BCSR)                    :: sP,sXk,sLXk,sT1
    INTEGER :: L
    !
    CALL Get(sP,QN%P)
    CALL Get(sXk,QN%Xk%Ortho)
    CALL Get(sLXk,QN%LXk%Ortho)
    !
    Ek=OneDotT(sP,sXk,sLXk)
    CALL Multiply(sXk,-Two*Ek)
    CALL Add(sXk,sLXk,sT1)
    !
    dEk=-1D10
    DO L=1,sT1%NNon0
       dEk=MAX(dEk,ABS(sT1%MTrix%D(L)))
    ENDDO
    Delta=MIN(Delta,ABS(dEk/Ek))
    !
    CALL Delete(sP)
    CALL Delete(sXk)
    CALL Delete(sLXk)
    CALL Delete(sT1)
  END SUBROUTINE EigenResolution

  SUBROUTINE PRStep(QN,K)!,BetaP,BetaQ)
    INTEGER             :: K
    TYPE(QUIRQIKontrol)   :: QN
    REAL(DOUBLE)        :: NumDot,DenDot,Beta,BetaP,BetaQ
    TYPE(BCSR)          :: sP,sPk,sGk,sGk_Old,sPk_Old,sT1,sT
    !
    IF(K==0)THEN
       !       CALL FileCopy(QN%Gk%Ortho,QN%Pk%Ortho)
       CALL FileCopy(QN%Gk%PosP, QN%Pk%PosP)
       CALL FileCopy(QN%Gk%MomQ, QN%Pk%MomQ)
    ELSE
       CALL New(sT1)
       CALL Get(sP,QN%P)
       CALL Multiply(sP,Two)
       !--------------------------------------------------
       CALL Get(sGk,QN%Gk%PosP)
       CALL Get(sGk_Old,QN%Gk%PosP_Old)
       !
       DenDot=OneDotT(sP,sGk_Old,sGk_Old)
       CALL Multiply(sGk_Old,-One)
       CALL Add(sGk,sGk_Old,sT1)
       NumDot=OneDotT(sP,sGk,sT1)
       Beta=MAX(Zero,NumDot/DenDot)

!!$       IF(MOD(K,10)==0)THEN
!!$          BetaP=Zero
!!$       ELSE
       BetaP=Beta
!!$       ENDIF

       !
       CALL Delete(sGk_Old)
       CALL Get(sPk_Old,QN%Pk%PosP_Old)
       CALL Multiply(sPk_Old,BetaP)
       CALL Add(sGk,sPk_Old,sPk)
       CALL Put(sPk,QN%Pk%PosP)
       !
       CALL Delete(sGk)
       CALL Delete(sT1)
       CALL Delete(sPk)
       CALL Delete(sPk_Old)
       !--------------------------------------------------
       CALL Get(sGk,QN%Gk%MomQ)
       CALL Get(sGk_Old,QN%Gk%MomQ_Old)
       !
       DenDot=OneDotT(sP,sGk_Old,sGk_Old)

       CALL Multiply(sGk_Old,-One)
       CALL Add(sGk,sGk_Old,sT1)
       NumDot=OneDotT(sP,sGk,sT1)
       Beta=MAX(Zero,NumDot/DenDot)
!!$
!!$       IF(MOD(K,10)==0)THEN
!!$          BetaQ=Zero
!!$       ELSE
       BetaQ=Beta
!!$       ENDIF
!!$
       CALL Delete(sGk_Old)
       CALL Get(sPk_Old,QN%Pk%MomQ_Old)
       CALL Multiply(sPk_Old,BetaQ)
       CALL Add(sGk,sPk_Old,sPk)
       CALL Put(sPk,QN%Pk%MomQ)
       !
       CALL Delete(sP)
       CALL Delete(sGk)
       CALL Delete(sT1)
       CALL Delete(sPk)
       CALL Delete(sPk_Old)
    ENDIF
    !
    CALL FileCopy(QN%Gk%PosP, QN%Gk%PosP_Old)
    CALL FileCopy(QN%Gk%MomQ, QN%Gk%MomQ_Old)
    IF(K==0)THEN
       CALL FileCopy(QN%Gk%PosP, QN%Pk%PosP_Old)
       CALL FileCopy(QN%Gk%MomQ, QN%Pk%MomQ_Old)
    ELSE
       CALL FileCopy(QN%Pk%PosP, QN%Pk%PosP_Old)
       CALL FileCopy(QN%Pk%MomQ, QN%Pk%MomQ_Old)
    ENDIF
    !
  END SUBROUTINE PRStep

  SUBROUTINE GradGrad(QN,Ek_P,Ek_Q,MaxGP,MaxGQ,NoQGrad_O)
    TYPE(FileNames)               :: Nam
    TYPE(State)                   :: S,RQIStat
    TYPE(Parallel)                :: MPI
    TYPE(QUIRQIKontrol)           :: QN
    TYPE(BCSR)                    :: sP,sXk_PosP,sXk_MomQ,sLXk_PosP,sLXk_MomQ, &
         sGk_PosP,sGk_MomQ,sXk,sLXk,sGk
    REAL(DOUBLE)                  :: DenDot,Ek,Ek_P,Ek_Q,MaxGP,MaxGQ
    INTEGER :: I
    LOGICAL,OPTIONAL :: NoQGrad_O
    !----------------------------------------------------------------------------
    CALL Elapsed_TIME(TimeGRAD,Init_O='Start')
    IF(PRESENT(NoQGrad_O))THEN
       ! Compute composit RQI gradients
       CALL Get(sP,QN%P)
       CALL Get(sXk,QN%Xk%Ortho)
       CALL Get(sLXk,QN%LXk%Ortho)
       !       Ek=Two*OneDotT(sP,sLXk,sXk)

       Ek=Two*OneDotT(sP,sXk,sLXk)
       Ek_P=Ek
       Ek_Q=Ek
       ! Gk=Two*(LXk-Ek*Xk)
       CALL Multiply(sXk,-Ek)
       CALL Add(sLXk,sXk,sGk)
       CALL Multiply(sGk,Two)
       !
       MaxGP=Zero
       DO I=1,sGk%NNon0
          MaxGP=MAX(MaxGP,ABS(sGk%MTrix%D(I)))
       ENDDO
       MaxGQ=MaxGP
       !
       CALL Put(sGk,QN%Gk%Ortho)
       !
       CALL Delete(sP)
       CALL Delete(sXk)
       CALL Delete(sLXk)
       CALL Delete(sGk)
    ELSE
       ! Compute split QUIRQI gradients
       CALL Get(sP,QN%P)
       CALL Multiply(sP,Two)
       CALL Get(sXk_PosP,QN%Xk%PosP)
       CALL Get(sXk_MomQ,QN%Xk%MomQ)
       CALL Get(sLXk_PosP,QN%LXk%PosP)
       CALL Get(sLXk_MomQ,QN%LXk%MomQ)
       DenDot=OneDotT(sP,sXk_PosP,sXk_MomQ)
       Ek_P=OneDotT(sP,sXk_PosP,sLXk_PosP)/DenDot
       Ek_Q=OneDotT(sP,sXk_MomQ,sLXk_MomQ)/DenDot
       CALL Delete(sP)
       CALL Multiply(sLXk_PosP,-One)
       CALL Multiply(sLXk_MomQ,-One)
       CALL Multiply(sXk_MomQ,Ek_Q)
       CALL Multiply(sXk_PosP,Ek_P)
       CALL Add(sLXk_PosP,sXk_MomQ,sGk_PosP)
       CALL Add(sLXk_MomQ,sXk_PosP,sGk_MomQ)
       !
       MaxGP=Zero
       MaxGQ=Zero
       DO I=1,sGk_PosP%NNon0
          MaxGP=MAX(MaxGP,ABS(sGk_PosP%MTrix%D(I)))
       ENDDO
       !
       DO I=1,sGk_MomQ%NNon0
          MaxGQ=MAX(MaxGQ,ABS(sGk_MomQ%MTrix%D(I)))
       ENDDO
       !
       CALL Put(sGk_PosP,QN%Gk%PosP)
       CALL Put(sGk_MomQ,QN%Gk%MomQ)
       !
       CALL Delete(sXk_PosP)
       CALL Delete(sXk_MomQ)
       CALL Delete(sGk_PosP)
       CALL Delete(sGk_MomQ)
       CALL Delete(sLXk_PosP)
       CALL Delete(sLXk_MomQ)
    ENDIF

    CALL Elapsed_TIME(TimeGRAD,Init_O='Accum')
  END SUBROUTINE GradGrad

  FUNCTION LinSrch(QN) !,DoTDA,Ek_P,Ek_Q)
    LOGICAL             :: LinSrch
    TYPE(QUIRQIKontrol) :: QN
    INTEGER             :: I,JJ,KK
    LOGICAL             :: DoTDA
    REAL(DOUBLE)        :: Ap,Aq,Bp,Bq,Cp,Cq,Rpq,Spq,Tpq,Upq, &
         LamP,LamQ,PosP_MIN,MomQ_MIN,dLamP,dLamQ,LamP_Old,LamQ_Old, &
         E_EOM,dEOM,EOM_Old,Ek_P,Ek_Q,Den,Tmp_EOM,E_EOM_AT_A
    TYPE(BCSR)          :: sP,sT1,sPosP_Xk,sMomQ_Xk,sL_PosP_Xk, &
         sL_MomQ_Xk,sPosP_Pk,sMomQ_Pk,sL_PosP_Pk,sL_MomQ_Pk
    !----------------------------------------------------------------------------------------
    CALL Elapsed_TIME(TimeSRCH,Init_O='Start')
    CALL Get(sP,QN%P)
    CALL Multiply(sP,Two)
    CALL Get(sPosP_Xk,QN%Xk%PosP)
    CALL Get(sPosP_Pk,QN%Pk%PosP)
    CALL Get(sL_PosP_Xk,QN%LXk%PosP)
    CALL Get(sL_PosP_Pk,QN%LPk%PosP)
    !
    Ap=OneDotT(sP,sPosP_Xk,sL_PosP_Xk)
    Bp=OneDotT(sP,sPosP_Pk,sL_PosP_Xk)+OneDotT(sP,sPosP_Xk,sL_PosP_Pk)
    Cp=OneDotT(sP,sPosP_Pk,sL_PosP_Pk)
    !
    CALL Delete(sPosP_Xk)
    CALL Delete(sPosP_Pk)
    CALL Delete(sL_PosP_Xk)
    CALL Delete(sL_PosP_Pk)
    !
    CALL Get(sMomQ_Xk,QN%Xk%MomQ)
    CALL Get(sMomQ_Pk,QN%Pk%MomQ)
    CALL Get(sL_MomQ_Xk,QN%LXk%MomQ)
    CALL Get(sL_MomQ_Pk,QN%LPk%MomQ)
    !
    Aq=OneDotT(sP,sMomQ_Xk,sL_MomQ_Xk)
    Bq=OneDotT(sP,sMomQ_Pk,sL_MomQ_Xk)+OneDotT(sP,sMomQ_Xk,sL_MomQ_Pk)
    Cq=OneDotT(sP,sMomQ_Pk,sL_MomQ_Pk)
    !
    CALL Delete(sL_MomQ_Xk)
    CALL Delete(sL_MomQ_Pk)
    CALL Get(sPosP_Xk,QN%Xk%PosP)
    CALL Get(sPosP_Pk,QN%Pk%PosP)
    !
    Rpq=OneDotT(sP,sPosP_Xk,sMomQ_Xk)
    Spq=OneDotT(sP,sPosP_Pk,sMomQ_Xk)
    Tpq=OneDotT(sP,sPosP_Xk,sMomQ_Pk)
    Upq=OneDotT(sP,sPosP_Pk,sMomQ_Pk)

!!$    WRITE(*,*)' Ap = ',Ap
!!$    WRITE(*,*)' Aq = ',Aq
!!$    WRITE(*,*)' Bp = ',Bp
!!$    WRITE(*,*)' Bq = ',Bq
!!$    WRITE(*,*)' Cp = ',Cp
!!$    WRITE(*,*)' Cq = ',Cq
!!$    WRITE(*,*)' Rpq = ',Rpq
!!$    WRITE(*,*)' Spq = ',Spq
!!$    WRITE(*,*)' Tpq = ',Tpq
!!$    WRITE(*,*)' Upq = ',Upq
    !
    CALL Delete(sP)
    !
    E_EOM=Half*(Ap+Aq)/Rpq
    !!    WRITE(*,44)E_EOM,Spq,Tpq,Upq
44  FORMAT(' E_EOM TO START = ',F20.10,' Spq = ',D20.10,' Tpq = ',D20.10,' Upq = ',D20.10)

    PosP_MIN=Zero
    MomQ_MIN=Zero
    !
    DO JJ=-100,100
       LamP=JJ*0.1D0
       DO KK=-100,100
          LamQ=KK*0.1D0
          Tmp_EOM=Half*(Ap+Bp*LamP+Cp*LamP**2+Aq+Bq*LamQ+Cq*LamQ**2)/(Rpq+LamP*Spq+LamQ*Tpq+LamP*LamQ*Upq)
          IF(Tmp_EOM<E_EOM.AND.Tmp_EOM>0D0)THEN
             E_EOM=Tmp_EOM
             PosP_MIN=LamP
             MomQ_MIN=LamQ
          ENDIF
       ENDDO
    ENDDO
    LamP=PosP_MIN
    LamQ=MomQ_MIN
    !
    E_EOM_AT_A=E_EOM

    !!    WRITE(*,46)E_EOM,LamP,LamQ
46  FORMAT(' E_EOM AT PT  A = ',F20.10,' LapP = ',D20.10,' LamQ = ',D20.10)
    DO I=1,1000
       LamQ_Old=LamQ
       LamQ=(-2D0*Cq*Rpq - 2D0*Cq*LamP*Spq +  Sqrt((2D0*Cq*Rpq + 2D0*Cq*LamP*Spq)**2 - &
            4D0*(Cq*Tpq + Cq*LamP*Upq)*(Bq*Rpq + Bq*LamP*Spq - Ap*Tpq - Aq*Tpq - Bp*LamP*Tpq - Cp*LamP**2*Tpq - Ap*LamP*Upq - &
            Aq*LamP*Upq - Bp*LamP**2*Upq - Cp*LamP**3*Upq)))/(2D0*(Cq*Tpq + Cq*LamP*Upq))
       LamP_Old=LamP
       LamP= (-2D0*Cp*Rpq - 2D0*Cp*LamQ*Tpq + Sqrt((2D0*Cp*Rpq + 2D0*Cp*LamQ*Tpq)**2 -          &
            4d0*(Cp*Spq + Cp*LamQ*Upq)*(Bp*Rpq - Ap*Spq - Aq*Spq - Bq*LamQ*Spq - Cq*LamQ**2*Spq + Bp*LamQ*Tpq - Ap*LamQ*Upq - &
            Aq*LamQ*Upq - Bq*LamQ**2*Upq - Cq*LamQ**3*Upq)))/(2D0*(Cp*Spq + Cp*LamQ*Upq))
       EOM_Old=E_EOM
       E_EOM=Half*(Ap+Bp*LamP+Cp*LamP**2+Aq+Bq*LamQ+Cq*LamQ**2)/(Rpq+LamP*Spq+LamQ*Tpq+LamP*LamQ*Upq)
       IF(E_EOM_AT_A<E_EOM)EXIT
       dLamQ=ABS((LamQ_Old-LamQ)/LamQ)
       dLamP=ABS((LamP_Old-LamP)/LamP)
       dEOM=ABS((EOM_Old-E_EOM)/E_EOM)
       IF(dEOM<1D-8.AND.I>4)EXIT
    ENDDO
    !
    IF( (I==1001).OR.(ABS(LamP)<1D-8).OR.(ABS(LamQ)<1D-8) )THEN
       LamP=-Bp/(2D0*Cp)
       LamQ=-Bq/(2D0*Cq)
    ENDIF
    !
    Den=(Rpq+LamP*Spq+LamQ*Tpq+LamP*LamQ*Upq)
    Ek_P=Half*(Ap+Bp*LamP+Cp*LamP**2)/Den
    Ek_Q=Half*(Aq+Bq*LamQ+Cq*LamQ**2)/Den
    E_EOM=Ek_P+Ek_Q
    !
    IF(E_EOM<Zero.OR.E_EOM>E_EOM_AT_A)THEN
       LamP=PosP_MIN
       LamQ=MomQ_MIN
       Den=(Rpq+LamP*Spq+LamQ*Tpq+LamP*LamQ*Upq)
       Ek_P=Half*(Ap+Bp*LamP+Cp*LamP**2)/Den
       Ek_Q=Half*(Aq+Bq*LamQ+Cq*LamQ**2)/Den
       E_EOM=Ek_P+Ek_Q
       WRITE(*,46)Half*(Ap+Aq)/Rpq,0D0,0D0
       WRITE(*,46)E_EOM,LamP,LamQ
    ENDIF

    !!    WRITE(*,47)E_EOM,LamP,LamQ,I
47  FORMAT(' E_EOM AT PT  B = ',F20.10,' LapP = ',D20.10,' LamQ = ',D20.10,' ICONVERGE = ',I8)

    !
    CALL FileCopy(QN%Xk%PosP,QN%Xk%PosP_Old)
    CALL FileCopy(QN%Xk%MomQ,QN%Xk%MomQ_Old)
    !
    Den=SQRT(Den)
    !
    CALL Multiply(sPosP_Pk,LamP)
    CALL Add(sPosP_Xk,sPosP_Pk,sT1)
    CALL Multiply(sT1,One/Den)
    CALL Put(sT1,QN%Xk%PosP)
    !
    CALL Delete(sT1)
    CALL Delete(sPosP_Xk)
    CALL Delete(sPosP_Pk)
    !
    CALL Multiply(sMomQ_Pk,LamQ)
    CALL Add(sMomQ_Xk,sMomQ_Pk,sT1)
    CALL Multiply(sT1,One/Den)

    CALL Put(sT1,QN%Xk%MomQ)
    !
    CALL Delete(sT1)
    CALL Delete(sMomQ_Xk)
    CALL Delete(sMomQ_Pk)
    !
    CALL Elapsed_TIME(TimeSRCH,Init_O='Accum')
    LinSrch=.TRUE.
    RETURN
111 CONTINUE
    LinSrch=.FALSE.
    CALL Warn('Line search failed. QUIRQI is likely converged.')
    CALL Elapsed_TIME(TimeSRCH,Init_O='Accum')
    !-----------------------------------------------------------------
    WRITE(*,*)' Ap = ',Ap
    WRITE(*,*)' Aq = ',Ap
    WRITE(*,*)' Bp = ',Bp
    WRITE(*,*)' Bq = ',Bq
    WRITE(*,*)' Cp = ',Cp
    WRITE(*,*)' Cq = ',Cq
    WRITE(*,*)' Rpq = ',Rpq
    WRITE(*,*)' Spq = ',Spq
    WRITE(*,*)' Tpq = ',Tpq
    WRITE(*,*)' Upq = ',Upq
    LamP=PosP_MIN
    LamQ=MomQ_MIN
    WRITE(*,446)E_EOM_AT_A,LamP,LamQ
446 FORMAT(' E_EOM AT PT  A = ',F20.10,' LapP = ',D20.10,' LamQ = ',D20.10)
    DO I=1,10
       LamQ_Old=LamQ
       LamQ=(-2D0*Cq*Rpq - 2D0*Cq*LamP*Spq +  Sqrt((2D0*Cq*Rpq + 2D0*Cq*LamP*Spq)**2 - &
            4D0*(Cq*Tpq + Cq*LamP*Upq)*(Bq*Rpq + Bq*LamP*Spq - Ap*Tpq - Aq*Tpq - Bp*LamP*Tpq - Cp*LamP**2*Tpq - Ap*LamP*Upq - &
            Aq*LamP*Upq - Bp*LamP**2*Upq - Cp*LamP**3*Upq)))/(2D0*(Cq*Tpq + Cq*LamP*Upq))
       LamP_Old=LamP
       LamP= (-2D0*Cp*Rpq - 2D0*Cp*LamQ*Tpq + Sqrt((2D0*Cp*Rpq + 2D0*Cp*LamQ*Tpq)**2 -          &
            4d0*(Cp*Spq + Cp*LamQ*Upq)*(Bp*Rpq - Ap*Spq - Aq*Spq - Bq*LamQ*Spq - Cq*LamQ**2*Spq + Bp*LamQ*Tpq - Ap*LamQ*Upq - &
            Aq*LamQ*Upq - Bq*LamQ**2*Upq - Cq*LamQ**3*Upq)))/(2D0*(Cp*Spq + Cp*LamQ*Upq))
       EOM_Old=E_EOM
       E_EOM=Half*(Ap+Bp*LamP+Cp*LamP**2+Aq+Bq*LamQ+Cq*LamQ**2)/(Rpq+LamP*Spq+LamQ*Tpq+LamP*LamQ*Upq)
       dLamQ=ABS((LamQ_Old-LamQ)/LamQ)
       dLamP=ABS((LamP_Old-LamP)/LamP)
       dEOM=ABS((EOM_Old-E_EOM)/E_EOM)
       WRITE(*,223)LamP,LamQ,E_EOM
223    FORMAT(' Lp = ',F12.4,', Lq = ',F12.4,' E = ',D16.6)
    ENDDO
    STOP ' Bad result in LinSrch '
  END FUNCTION LinSrch

  SUBROUTINE LOn2AO(Nam,MPI,S,QN,RQIStat,Xk,LXk,KTau_O,JTau_O)
    TYPE(FileNames)               :: Nam
    TYPE(State)                   :: S,RQIStat
    TYPE(Parallel)                :: MPI
    TYPE(QUIRQIKontrol)             :: QN
    TYPE(PQName)                  :: Xk,LXk
    TYPE(BCSR)                    :: sF,sX,sJ,sK,sZ,sP,sJK,sT1,sT2,sT3
    REAL(DOUBLE),OPTIONAL         :: JTau_O,KTau_O
    !----------------------------------------------------------------------------
    CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
    CALL New(sX) ! Kluge for now.  Shouldn't have to do this
    CALL Get(sX,Xk%Ortho)
    CALL Get(sZ,QN%Z)
    ! Xao = Z^t.Xor.Z
    CALL Multiply(sZ,sX,sT1)
    CALL Get(sZ,QN%ZT)
    CALL Multiply(sT1,sZ,sX)
    CALL Put(sX,Xk%AO) ! AO transition density matrix (or CG gradient)
    ! Done with sX and sZ for now
    CALL Delete(sZ)
    CALL Delete(sX)
    CALL Delete(sT1)
    ! To here, all big memory should have been purged
    CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
    ! Build JK[X] in the AO basis
    RQIStat%Action%C(2)=TRIM(Xk%Act2)
    !---------------------------------------------------------------------------------------
    ! ONX CALL:
    IF(PRESENT(KTau_O))THEN
       HDFFileID=OpenHDF(Nam%HFile)
       HDF_CurrentID=HDFFileID
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #1")
       QN%Current%Dist=KTau_O*1D-2
       QN%Current%TwoE=KTau_O
       CALL Put(QN%Current,QN%chBas)
       CALL CloseHDFGroup(HDF_CurrentID)
       CALL CloseHDF(HDFFileID)
    ENDIF
    CALL Elapsed_TIME(TimeONX,Init_O='Start');  CALL Invoke('ONX',Nam,RQIStat,MPI);    CALL Elapsed_TIME(TimeONX,Init_O='Accum')
    IF(PRESENT(KTau_O))THEN
       HDFFileID=OpenHDF(Nam%HFile)
       HDF_CurrentID=HDFFileID
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #1")
       CALL Put(QN%Ultimate,QN%chBas)
       CALL CloseHDFGroup(HDF_CurrentID)
       CALL CloseHDF(HDFFileID)
    ENDIF
    !---------------------------------------------------------------------------------------
    CALL Elapsed_TIME(TimeQCTC,Init_O='Start'); CALL Invoke('QCTC',Nam,RQIStat,MPI);   CALL Elapsed_TIME(TimeQCTC,Init_O='Accum')
    CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
    ! Pick up J and K
    CALL Get(sJ,QN%J)
    CALL Get(sK,QN%K)
    ! JK=Jao[X]+Kao[X]
    CALL Add(sJ,sK,sJK)
    CALL Delete(sJ)
    CALL Delete(sK)
    ! JK[X]=Zt.JKao[X].Z==JKor
    CALL Get(sZ,QN%ZT)
    CALL Multiply(sZ,sJK,sT1)
    CALL Get(sZ,QN%Z)
    CALL Multiply(sT1,sZ,sJK)
    ! Done with Z
    CALL Delete(sZ)
    CALL Get(sP,QN%P)
    ! T1=[JKor,Por]
    CALL Multiply(sP,sJK,sT1)
    CALL Multiply(sJK,sP,sT1,-One)
    !
    CALL Delete(sJK)
    !
    CALL Get(sF,QN%F)
    CALL Get(sX,Xk%Ortho)
    ! T2=[Xor,For]
    CALL Multiply(sX,sF,sT2)
    CALL Multiply(sF,sX,sT2,-One)
    !    CALL Multiply(sT2,OneERPAScaling)
    ! Done with F
    CALL Delete(sF)
    CALL Delete(sP)
    CALL Delete(sX)
    ! L[Xk]=[F,Xk]+[P,JK[X]] (orthogonal)
    CALL Add(sT1,sT2,sT3)
    ! Put orthogonal L[Xk] or L[Pk] to disk
    CALL Put(sT3,LXk%Ortho)
    ! Tidy up
    CALL Delete(sT1)
    CALL Delete(sT2)
    CALL Delete(sT3)
    !
    CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
  END SUBROUTINE LOn2AO

  SUBROUTINE LOn2AOpq(Nam,MPI,S,QN,RQIStat,Xk,LXk,JTau_O,KTau_O)
    TYPE(FileNames)               :: Nam
    TYPE(State)                   :: S,RQIStat
    TYPE(Parallel)                :: MPI
    TYPE(QUIRQIKontrol)           :: QN
    TYPE(PQName)                  :: Xk,LXk
    TYPE(BCSR)                    :: sF,sX,sJ,sK,sZ,sP,sQ,sJK,sT1,sT2,sT3
    REAL(DOUBLE),OPTIONAL         :: JTau_O,KTau_O
    REAL(DOUBLE)                  :: JTmp,KTmp
    !----------------------------------------------------------------------------
    CALL Elapsed_TIME(TimeLON2,Init_O='Start')
    !===================================================================
    ! DO THE PosP part
    !
    CALL New(sX) ! Kluge for now.  Shouldn't have to do this
    CALL Get(sX,Xk%PosP)

    !    IF(QN%CholFact)THEN
    !       CALL Get(sZ,QN%ZT)
    !    ELSE
    CALL Get(sZ,QN%Z)

    !    ENDIF
    ! Xao = Z^t.Xor.Z
    CALL Multiply(sZ,sX,sT1)
    IF(QN%CholFact) &
         CALL Get(sZ,QN%ZT)

    CALL Multiply(sT1,sZ,sX)
    CALL Put(sX,Xk%AO) ! AO transition density matrix (or CG gradient)
    ! Done with sX and sZ for now
    CALL Delete(sZ)
    CALL Delete(sX)
    CALL Delete(sT1)
    ! To here, all big memory should have been purged
    CALL Elapsed_TIME(TimeLON2,Init_O='Accum')
    ! Build JK[X] in the AO basis
    RQIStat%Action%C(2)=TRIM(Xk%Act2)
    !---------------------------------------------------------------------------------------
    CALL Elapsed_TIME(TimeONX,Init_O='Start');  CALL Invoke('ONX',Nam,RQIStat,MPI);    CALL Elapsed_TIME(TimeONX,Init_O='Accum')
    !---------------------------------------------------------------------------------------
    CALL Elapsed_TIME(TimeQCTC,Init_O='Start'); CALL Invoke('QCTC',Nam,RQIStat,MPI);   CALL Elapsed_TIME(TimeQCTC,Init_O='Accum')
    CALL Elapsed_TIME(TimeLON2,Init_O='Start')
    ! Pick up J and K
    CALL Get(sJ,QN%J)
    CALL Get(sK,QN%K)
    ! JK=Jao[X]+Kao[X]
    CALL Add(sJ,sK,sJK)
    CALL Delete(sJ)
    CALL Delete(sK)
    ! JK[X]=Zt.JKao[X].Z==JKor
    CALL Get(sZ,QN%Z)
    CALL Multiply(sZ,sJK,sT1)
    IF(QN%CholFact) &
         CALL Get(sZ,QN%ZT)
    ! Xao = Z^t.Xor.Z
    CALL Multiply(sZ,sX,sT1)
    IF(QN%CholFact) &
         CALL Get(sZ,QN%Z)
    CALL Multiply(sT1,sZ,sX)

    CALL Multiply(sT1,sZ,sJK)
    ! Done with Z
    CALL Delete(sZ)
    CALL Get(sP,QN%P)
    ! T1=[JKor,Por]
    CALL Multiply(sP,sJK,sT1)
    CALL Multiply(sJK,sP,sT1,-One)
    CALL Delete(sJK)
    CALL Get(sF,QN%F)
    CALL Get(sX,Xk%PosP)
    ! T2=[Xor,For]
    CALL Multiply(sX,sF,sT2)
    CALL Multiply(sF,sX,sT2,-One)
    ! Done with F
    CALL Delete(sF)
    CALL Delete(sP)
    CALL Delete(sX)
    ! L[Xk]=[F,Xk]+[P,JK[X]] (orthogonal)
    CALL Add(sT1,sT2,sT3)
    ! Put orthogonal L[Xk] or L[Pk] to disk
    CALL Put(sT3,LXk%PosP)
    ! Tidy up
    CALL Delete(sT1)
    CALL Delete(sT2)
    CALL Delete(sT3)
    !===================================================================
    ! DO THE MomQ part
    !
    CALL New(sX) ! Kluge for now.  Shouldn't have to do this
    CALL Get(sX,Xk%MomQ)
    CALL Get(sZ,QN%Z)
    ! Xao = Z^t.Xor.Z
    CALL Multiply(sZ,sX,sT1)
    CALL Multiply(sT1,sZ,sX)
    CALL Put(sX,Xk%AO) ! AO transition density matrix (or CG gradient)
    ! Done with sX and sZ for now
    CALL Delete(sZ)
    CALL Delete(sX)
    CALL Delete(sT1)
    ! To here, all big memory should have been purged
    CALL Elapsed_TIME(TimeLON2,Init_O='Accum')
    ! Build JK[X] in the AO basis
    RQIStat%Action%C(2)=TRIM(Xk%Act2)
    CALL Elapsed_TIME(TimeONX,Init_O='Start');  CALL Invoke('ONX',Nam,RQIStat,MPI);    CALL Elapsed_TIME(TimeONX,Init_O='Accum')
    CALL Elapsed_TIME(TimeQCTC,Init_O='Start'); CALL Invoke('QCTC',Nam,RQIStat,MPI);   CALL Elapsed_TIME(TimeQCTC,Init_O='Accum')
    CALL Elapsed_TIME(TimeLON2,Init_O='Start')
    ! Pick up J and K
    CALL Get(sJ,QN%J)
    CALL Get(sK,QN%K)
    ! JK=Jao[X]+Kao[X]
    CALL Add(sJ,sK,sJK)
    CALL Delete(sJ)
    CALL Delete(sK)
    ! JK[X]=Zt.JKao[X].Z==JKor
    CALL Get(sZ,QN%Z)
    CALL Multiply(sZ,sJK,sT1)
    CALL Multiply(sT1,sZ,sJK)
    ! Done with Z
    CALL Delete(sZ)
    CALL Get(sP,QN%P)
    ! T1=[JKor,Por]
    CALL Multiply(sP,sJK,sT1)
    CALL Multiply(sJK,sP,sT1,-One)
    CALL Delete(sJK)
    CALL Get(sF,QN%F)
    CALL Get(sX,Xk%MomQ)
    ! T2=[Xor,For]
    CALL Multiply(sX,sF,sT2)
    CALL Multiply(sF,sX,sT2,-One)
    ! Done with F
    CALL Delete(sF)
    CALL Delete(sP)
    CALL Delete(sX)
    ! L[Xk]=[F,Xk]+[P,JK[X]] (orthogonal)
    CALL Add(sT1,sT2,sT3)
    ! Put orthogonal L[Xk] or L[Pk] to disk
    !    CALL PPrint(sT3,'pq L[MomQ] RAW ',Unit_O=6)
    CALL Put(sT3,LXk%MomQ)
    ! Tidy up
    CALL Delete(sT1)
    CALL Delete(sT2)
    CALL Delete(sT3)
    !===================================================================
    ! Now, do the folding
    !
    CALL Get(sP,QN%P)
    CALL Multiply(sP,Two)
    CALL SetEq(sQ,sP)
    CALL Multiply(sQ,-One)
    CALL Add(sQ,Two)
    ! - - - - - - - - - - - - - - - - - - - - - - - -
    ! First the PosP part
    CALL Get(sX,LXk%PosP)

    ! TRANSPOSE(MATMUL(Q,MATMUL(X,P)))
    CALL Multiply(sQ,sX,sT1)
    CALL Multiply(sT1,sP,sT2)
    CALL XPose(sT2,sT1)
    ! MATMUL(P,MATMUL(X,Q))
    CALL Multiply(sP,sX,sT2)
    CALL Multiply(sT2,sQ,sT3)
    ! PosP=P.X.Q-TRANSPOSE(Q.X.P)
    CALL Multiply(sT1,-One)
    CALL Add(sT1,sT3,sT2)
    CALL Multiply(sT2,25D-2)
    !
    !    CALL PPrint(sT2,'pq L[MomQ]',Unit_O=6)
    !
    CALL Put(sT2,LXk%PosP)
    CALL Delete(sX)
    CALL Delete(sT1)
    CALL Delete(sT2)
    CALL Delete(sT3)
    ! - - - - - - - - - - - - - - - - - - - - - - - -
    CALL Get(sX,LXk%MomQ)
    ! TRANSPOSE(MATMUL(Q,MATMUL(X,P)))
    CALL Multiply(sQ,sX,sT1)
    CALL Multiply(sT1,sP,sT2)
    CALL XPose(sT2,sT1)
    ! MATMUL(P,MATMUL(X,Q))
    CALL Multiply(sP,sX,sT2)
    CALL Multiply(sT2,sQ,sT3)
    ! PosP=P.X.Q+TRANSPOSE(Q.X.P)
    CALL Add(sT1,sT3,sT2)
    CALL Multiply(sT2,25D-2)
    !    CALL PPrint(sT2,'pq L[PosP]',Unit_O=6)
    CALL Put(sT2,LXk%MomQ)
    ! - - - - - - - - - - - - - - - - - - - - - - - -
    CALL Delete(sX)
    CALL Delete(sT1)
    CALL Delete(sT2)
    CALL Delete(sT3)
    CALL Delete(sP)
    CALL Delete(sQ)
    CALL Elapsed_TIME(TimeLON2,Init_O='Accum')
  END SUBROUTINE LOn2AOpq

  SUBROUTINE Recomb(QN,Xk)
    TYPE(QUIRQIKontrol)     :: QN
    TYPE(PQName)          :: Xk
    TYPE(BCSR)            :: sXk,sPosP,sMomQ,sT1,sT2,sT3
    CALL Elapsed_TIME(TimeRECO,Init_O='Start')
    ! Xk=Half*(PosP_Xk+MomQ_Xk+TRANSPOSE(PosP_Xk-MomQ_Xk))
    CALL Get(sPosP,Xk%PosP)
    CALL Get(sMomQ,Xk%MomQ)
    !    WRITE(*,*)'recomb recomb recomb recomb recomb '
    !    CALL PChkSum(sPosP,Xk%PosP,Unit_O=6)
    !    CALL PChkSum(sMomQ,Xk%MomQ,Unit_O=6)
    !    WRITE(*,*)'- - - - - - - - - - - - - - - - - - '
    CALL Multiply(sMomQ,-One)
    CALL Add(sPosP,sMomQ,sT1)
    CALL XPose(sT1,sT2)
    CALL Multiply(sMomQ,-One)
    CALL Add(sPosP,sMomQ,sT1)
    CALL Add(sT1,sT2,sT3)
    CALL Multiply(sT3,Half)
    CALL Put(sT3,Xk%Ortho)
    ! Tidy
    CALL Delete(sT1)
    CALL Delete(sT2)
    CALL Delete(sT3)
    CALL Delete(sPosP)
    CALL Delete(sMomQ)
    ! All done
    CALL Elapsed_TIME(TimeRECO,Init_O='Accum')
  END SUBROUTINE Recomb

  SUBROUTINE SplitPQ(QN,Xk,Switch_O,PQThreshold_O)
    TYPE(QUIRQIKontrol)   :: QN
    TYPE(PQName)          :: Xk
    TYPE(BCSR)            :: sP,sXk,sQ,sT1,sT2,sT3
    LOGICAL,OPTIONAL      :: Switch_O
    REAL(DOUBLE),OPTIONAL :: PQThreshold_O
    !
    CALL Elapsed_TIME(TimeSPLT,Init_O='Start')
    CALL Get(sP,QN%P)
    CALL Multiply(sP,Two)
    CALL SetEq(sQ,sP)
    CALL Multiply(sQ,-One)
    CALL Add(sQ,Two)
    CALL Get(sXk,Xk%Ortho)
    !
    IF(PRESENT(PQThreshold_O))THEN
       CALL Filter(sP,PQThreshold_O)
       CALL Filter(sQ,PQThreshold_O)
    ENDIF
    !
    ! TRANSPOSE(MATMUL(Q,MATMUL(X,P)))
    CALL Multiply(sQ,sXk,sT1)
    CALL Multiply(sT1,sP,sT2)
    CALL XPose(sT2,sT1)
    ! MATMUL(P,MATMUL(X,Q))
    CALL Multiply(sP,sXk,sT2)
    CALL Multiply(sT2,sQ,sT3)
    ! PosP=P.X.Q+TRANSPOSE(Q.X.P)
    CALL Add(sT1,sT3,sT2)
    CALL Multiply(sT2,25D-2)
    !
    IF(PRESENT(Switch_O))THEN
       IF(Switch_O)THEN
          CALL Put(sT2,Xk%PosP)
       ELSE
          CALL Put(sT2,Xk%MomQ)
       ENDIF
    ELSE
       CALL Put(sT2,Xk%MomQ)
    ENDIF
    ! MomQ=P.X.Q-TRANSPOSE(Q.X.P)
    CALL Multiply(sT1,-One)
    CALL Add(sT1,sT3,sT2)
    CALL Multiply(sT2,25D-2)
    IF(PRESENT(Switch_O))THEN
       IF(Switch_O)THEN
          CALL Put(sT2,Xk%MomQ)
       ELSE
          CALL Put(sT2,Xk%PosP)
       ENDIF
    ELSE
       CALL Put(sT2,Xk%PosP)
    ENDIF
    ! Tidy up
    CALL Delete(sP)
    CALL Delete(sQ)
    CALL Delete(sXk)
    CALL Delete(sT1)
    CALL Delete(sT2)
    CALL Delete(sT3)
    ! All done
    CALL Elapsed_TIME(TimeSPLT,Init_O='Accum')
  END SUBROUTINE SplitPQ

  SUBROUTINE SNihilate(QN,Xk,DoTDA)
    REAL(DOUBLE)          :: Norm
    TYPE(PQName)          :: Xk
    TYPE(QUIRQIKontrol)     :: QN
    CHARACTER(LEN=DCL)    :: XName
    LOGICAL               :: DoNorm,DoTDA
    TYPE(BCSR)            :: sP,sQ,sXk,sI,sT1,sT2,sT3
    !!       TYPE(DBL_RNK2) :: PP,QQ,XX,BB
    !-----------------------------------------------------------------------------
    CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
    !
    CALL New(sXk) ! Kluge, shouldn't have to do this!!
    CALL Get(sXk,Xk%Ortho)
    CALL Get(sP,QN%P)
    CALL Get(sQ,QN%P)
    CALL Multiply(sP,Two)
    CALL Multiply(sQ,-Two)
    CALL Add(sQ,Two)
    IF(DoTDA)THEN
       CALL Multiply(sQ,sXk,sT1)
       CALL Multiply(sT1,sP,sXk)
       CALL Delete(sP)
       CALL Delete(sQ)
       CALL Delete(sT1)
    ELSE
       ! anihilate particle-hole and hole-particle
       ! Xk=(P.Xk.Q+Q.Xk.P) (both h-p and p-h zapped)
       CALL Multiply(sP,sXk,sT1)
       CALL Multiply(sT1,sQ,sT2)
       CALL Multiply(sQ,sXk,sT1)
       CALL Multiply(sT1,sP,sT3)
       CALL Add(sT2,sT3,sXk)
       CALL Delete(sP)
       CALL Delete(sQ)
       CALL Delete(sT1)
       CALL Delete(sT2)
       CALL Delete(sT3)
    ENDIF
    CALL Multiply(sXk,25D-2)
    CALL Put(sXk,Xk%Ortho)
    ! Tidy up
    CALL Delete(sXk)
    ! All done
    CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
  END SUBROUTINE SNihilate
  !==========================================================================
  !  END MAIN QUIRQI ROUTINES
  !==========================================================================

  SUBROUTINE sQUIRQI(Nam,O,S,MPI,G,B,QN)
    INTEGER            :: I,J,K
    TYPE(QUIRQIKontrol):: QN
    TYPE(FileNames)    :: Nam,RQINams
    TYPE(State)        :: S,RQIStat
    TYPE(Parallel)     :: MPI
    TYPE(Options)      :: O
    TYPE(Geometries)   :: G
    TYPE(BSET)         :: B
    LOGICAL            :: DoTDA
    REAL(DOUBLE)       :: Ek_P,Ek_Q,Ek,ErrAbs,EkOld,ErrRel,MaxGP,MaxGQ, &
         XkPcntP,XkPcntQ,PkPcntP,PkPcntQ,XkPcnt,PkPcnt,GTau,XTau,ETau,MinRel,SqrRel
    TYPE(BCSR)         :: sP,sXk,sLXk,sT
    CHARACTER(LEN=DCL) :: Iteration,Statistics1,Statistics2,Statistics3
    INTEGER,PARAMETER  :: MaxQIts=100
    REAL(DOUBLE),DIMENSION(0:MaxQIts) :: EConverge

    WRITE(*,*)'SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSs'
    !
    CALL SetQName(Nam,S,QN,RQIStat)
    !
    QN%Threshold=1D-10
    CALL Koopmans(QN)
    !
    CALL Elapsed_TIME(TimeTotal,Init_O='Init')
    CALL Elapsed_TIME(TimeONX,Init_O='Init')
    CALL Elapsed_TIME(TimeQCTC,Init_O='Init')
    CALL Elapsed_TIME(TimeBCSR,Init_O='Init')
    CALL Elapsed_TIME(TimeTotal,Init_O='Start')
    !
    CALL SNihilate(QN,QN%Xk,.TRUE.)
    CALL SplitPQ(QN,QN%Xk,.TRUE.)
    CALL LOn2AOpq(Nam,MPI,S,QN,RQIStat,QN%Xk,QN%LXk)
    !
    XTau=1D-10
    GTau=1D-10
    !
    EkOld=1D2
    MinRel=1D2
    SqrRel=1D2
    !
    DO K=0,MaxQIts-1
       !
!!$       IF(K<6)THEN
!!$          XTau=1D-1
!!$!          DoTDA=.TRUE.
!!$       ELSE
!!$          XTau=1D-3
!!$!          DoTDA=.FALSE.
!!$       ENDIF
!!$
       CALL GradGrad(QN,Ek_P,Ek_Q,MaxGP,MaxGQ)
       !
       Ek=Half*(Ek_P+Ek_Q)
       !
       EConverge(K)=Ek
       ErrAbs=EkOld-Ek
       ErrRel=ErrAbs/Ek
       IF(MOD(K,2)==0) &
            EkOld=MIN(EkOld,Ek)
       !
       MinRel=MIN(MinRel,ABS(ErrRel))
       SqrRel=MIN(SqrRel,SQRT(ABS(ErrRel)))
       !
       !       XTau=MAX(MIN(XTau,1D-1*SqrRel),QN%Threshold)
       !       GTau=MIN(GTau,5D-1*MAX(MaxGP,MaxGQ)) !,QN%Threshold)

       CALL Elapsed_TIME(TimeTOTAL,Init_O='Accum')
       !
       !       IF(ErrRel<QN%Konvergence.OR.(ErrAbs<Zero.AND.XTau<=QN%Threshold))EXIT
       !
       IF(ErrAbs<Zero)THEN
          !          XTau=MAX(XTau*1D-1,QN%Threshold)
          !         GTau=MAX(GTau*1D-1,QN%Threshold)
       ENDIF
       !
       CALL QUIRQIPrint(QN,K,Ek,Ek_P,Ek_Q,ErrRel,QN%Current%TwoE,XTau,GTau,XkPcnt,PkPcnt,EConverge,Nam,.FALSE.)
       CALL PRStep(QN,K)
!!$       !
       CALL Recomb(QN,QN%Pk)
       CALL SNihilate(QN,QN%Pk,DoTDA)
       CALL SplitPQ(QN,QN%Pk,.TRUE.)
!!$       CALL PrunePH(QN,QN%Pk%PosP,.FALSE.,GTau,PkPcntP)
!!$       CALL PrunePH(QN,QN%Pk%MomQ,.FALSE.,GTau,PkPcntQ)
       !
       CALL LOn2AOpq(Nam,MPI,S,QN,RQIStat,QN%Pk,QN%LPk)
       !
       CALL Recomb(QN,QN%LPk)
       CALL SNihilate(QN,QN%LPk,DoTDA)
       CALL SplitPQ(QN,QN%LPk,.TRUE.)
!!$       CALL PrunePH(QN,QN%LPk%PosP,.FALSE.,1D-10,PkPcntP)
!!$       CALL PrunePH(QN,QN%LPk%MomQ,.FALSE.,1D-10,PkPcntQ)
       !
       IF(.NOT.LinSrch(QN))EXIT
       !
       CALL Recomb(QN,QN%Xk)
       CALL SNihilate(QN,QN%Xk,DoTDA)
       CALL SplitPQ(QN,QN%Xk,.TRUE.)

       CALL OrthoPrune(QN,QN%Xk%PosP,.FALSE.,XTau,XkPcntP)
       CALL OrthoPrune(QN,QN%Xk%MomQ,.FALSE.,XTau,XkPcntQ)
       XkPcnt=MAX(XkPcntP,XkPcntQ)

       !!       CALL PrunePH(QN,QN%Xk%PosP,.FALSE.,XTau,XkPcntP)
!!$       CALL PrunePH(QN,QN%Xk%MomQ,.FALSE.,XTau,XkPcntQ)

!!$!       WRITE(*,*)' XKPCNT = ',XkPcntP,XkPcntQ
       !
       CALL LOn2AOpq(Nam,MPI,S,QN,RQIStat,QN%Xk,QN%LXk)
       !
       CALL Recomb(QN,QN%LXk)
       CALL SNihilate(QN,QN%LXk,DoTDA)
       CALL SplitPQ(QN,QN%LXk,.TRUE.)
!!$       CALL PrunePH(QN,QN%LXk%PosP,.FALSE.,1D-10,PkPcntP)
!!$       CALL PrunePH(QN,QN%LXk%MomQ,.FALSE.,1D-10,PkPcntQ)
       !
    ENDDO
    !
    EConverge(MaxQIts)=E_RPA
    CALL QUIRQIPrint(QN,K,Ek,Ek_P,Ek_Q,ErrRel,QN%Current%TwoE,XTau,GTau,XkPcnt,PkPcnt,EConverge,Nam,.TRUE.)
    !
    !    WRITE(*,*)'Xk%Ortho = ',QN%Xk%Ortho
    !    CALL Recomb(QN,QN%Xk)
!!!    OrthoPrune(QN,Xk,NormQ,Threshold,XkPcnt)




  END SUBROUTINE sQUIRQI


  SUBROUTINE Split2PQ(QN,Xk)
    TYPE(QUIRQIKontrol)     :: QN
    TYPE(PQName)          :: Xk
    TYPE(BCSR)            :: sP,sXk,sQ,sT1,sT2,sT3
    TYPE(DBL_RNK2)        :: P,Q,X,BB

    CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
    CALL Get(sXk,Xk%Ortho)
    !
    CALL Get(sP,QN%P)
    CALL Get(sQ,QN%P)
    CALL Multiply(sP,Two)
    CALL Multiply(sQ,-Two)
    CALL Add(sQ,Two)
    ! T1=TRANSPOSE(MATMUL(Q,MATMUL(X,P)))
    CALL Multiply(sQ,sXk,sT1)
    CALL Multiply(sT1,sP,sT2)
    CALL XPose(sT2,sT1)
    ! T3=MATMUL(P,MATMUL(X,Q))
    CALL Multiply(sP,sXk,sT2)
    CALL Multiply(sT2,sQ,sT3)
    ! T2=PosP=P.X.Q+TRANSPOSE(Q.X.P)
    CALL Add(sT1,sT3,sT2)
    CALL Multiply(sT2,25D-2)
    CALL Put(sT2,Xk%PosP)
    ! T2=MomQ=P.X.Q-TRANSPOSE(Q.X.P)
    CALL Multiply(sT1,-One)
    CALL Add(sT3,sT1,sT2)
    CALL Multiply(sT2,25D-2)
    !-----------------------------------
    ! MomQ->PosP
    CALL SetEq(sXk,sT2)
    ! TRANSPOSE(MATMUL(Q,MATMUL(X,P)))
    CALL Multiply(sQ,sXk,sT1)
    CALL Multiply(sT1,sP,sT2)
    CALL XPose(sT2,sT1)
    ! MATMUL(P,MATMUL(X,Q))
    CALL Multiply(sP,sXk,sT2)
    CALL Multiply(sT2,sQ,sT3)
    ! PosP=P.X.Q+TRANSPOSE(Q.X.P)
    CALL Add(sT1,sT3,sT2)
    CALL Multiply(sT2,25D-2)
    !
    CALL Get(sXk,Xk%PosP)
    CALL Put(sT2,Xk%PosP)
    !-----------------------------------
    ! PosP->MomQ
    ! TRANSPOSE(MATMUL(Q,MATMUL(X,P)))
    CALL Multiply(sQ,sXk,sT1)
    CALL Multiply(sT1,sP,sT2)
    CALL XPose(sT2,sT1)
    CALL Multiply(sT1,-One)
    ! MATMUL(P,MATMUL(X,Q))
    CALL Multiply(sP,sXk,sT2)
    CALL Multiply(sT2,sQ,sT3)
    ! MomQ=P.X.Q-TRANSPOSE(Q.X.P)
    CALL Add(sT3,sT1,sT2)
    CALL Multiply(sT2,25D-2)
    CALL Put(sT2,Xk%MomQ)
    !       CALL PPrint(sT3,'Xk_Q = ',Unit_O=6)
    ! Tidy up
    CALL Delete(sP)
    CALL Delete(sQ)
    CALL Delete(sXk)
    CALL Delete(sT1)
    CALL Delete(sT2)
    CALL Delete(sT3)
    ! All done
    CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
  END SUBROUTINE Split2PQ






  SUBROUTINE Koopmans(QN,RelErrTol_O)
    !
    TYPE(QUIRQIKontrol)  :: QN
    INTEGER            :: I,J,K,L
    REAL(DOUBLE),OPTIONAL :: RelErrTol_O
    REAL(DOUBLE)       :: MaxP,Ek,EkOld,dEk,Beta,Lambda,ErrRel,ErrAbs,Err,Norm,Num,Den, &
         XkPcnt,PkPcnt,GkPcnt,LXkPcnt,LPkPcnt
    TYPE(BCSR)         :: sP,sQ,sF,sXk,sGk,sPk,sXkOld,sGkOld,sPkOld,sLXk,sLPk,sT1,sT2,sT3
    CHARACTER(LEN=DCL) :: Iteration,Statistics1,Statistics2,Statistics3

    TYPE(DBL_VECT)    :: EigenEs
    TYPE(DBL_RNK2)    :: X


    !

    CALL Get(sP,QN%P)
    CALL Get(sF,QN%F)
    CALL SetEq(sQ,sP)
    !
    CALL Multiply(sP,Two)
    CALL Multiply(sQ,-Two)
    CALL Add(sQ,Two)
    !
    CALL New(sT1)
    CALL New(sT2)
    CALL New(sXk)
    CALL New(sGk)
    CALL New(sPk)
    CALL New(sLXk)
    CALL New(sLPk)
    CALL New(sXkOld)
    CALL New(sGkOld)
    CALL New(sPkOld)
    !
    ! If we use random start we get non-deterministic behavior due to sigularities.
    !
    CALL SetEq(sXk,sP)
    !
    ! LXk=MATMUL(F%D,Xk)-MATMUL(Xk,F%D)
    DO I=1,sXk%NNon0
       MaxP=ABS(sXk%MTrix%D(I))
       sXk%MTrix%D(I)=RANDOM_DBL((/-MaxP,MaxP/))
    ENDDO

    !    DO I=1,sXk%NNon0
    !       MaxP=ABS(sXk%MTrix%D(I))/1D2
    !       sXk%MTrix%D(I)=sXk%MTrix%D(I)+RANDOM_DBL((/-MaxP,MaxP/))
    !    ENDDO
!!$
!!$    ! CONSIDER A LEVEL SHIFT TO ENHANCE DECAY OF THE GUESS!
!!$    ! The shifted Fockian is F[Lambda] = P.F + (1+Lambda)*Q.F
!!$    CALL Multiply(sQ,sF,sT1)
!!$    CALL Multiply(sT1,sQ,sT2)
!!$    CALL Multiply(sP,sF,sT1)
!!$    CALL Multiply(sT1,sP,sT3)
!!$    CALL Multiply(sT3,1D2)
!!$    CALL Delete(sF)
!!$    CALL Add(sT2,sT3,sF)


    !QUIRQI : Koopmans' 42  :: Ev = 0.47801931D+01, Err = -.59D-05, %X = 89.83703, %LX = 0.00000, %P = 90.75373, %LP = 0.00000
    !QUIRQI : RPA(1,0)      ::  Ev = 0.5773009627919349D+01, ErrRel = 0.47D+03, TrueEr = 0.66D+00, ETau = 0.10D-09,

    !QUIRQI : Koopmans' 178 :: Ev = 0.35769787D+04, Err = 0.12D-06, %X = 91.67044, %LX = 0.00000, %P = 91.67044, %LP = 0.00000
    !QUIRQI : RPA(1,0)      ::  Ev = 0.5216907115502461D+01, ErrRel = 0.52D+03, TrueEr = 0.50D+00, ETau = 0.10D-09, XTau = 0.10D-09, GTau = 0.10D-09, %X = 0.44-322, %P = 0.45-316



    ! LXk=MATMUL(F%D,Xk)-MATMUL(Xk,F%D)
    CALL Multiply(sF,sXk,sT1)
    CALL Multiply(sXk,sF,sT2)
    CALL Multiply(sT2,-One)
    CALL Add(sT1,sT2,sLXk)
    ! LXk=0.25D0*(MATMUL(MATMUL(Q%D,LXk),P%D))
    CALL Multiply(sQ,sLXk,sT1)
    CALL Multiply(sT1,sP,sLXk)
    CALL Multiply(sLXk,25D-2)
    ! Ek=ThoulessQ(N,P%D,Xk,LXk)
    EkOld=1D100
    Ek=OneDotT(sP,sXk,sLXk)
    DO K=0,300
       ! Gk=Two*(LXk-Ek*Xk)
       CALL SetEq(sT1,sXk)
       CALL Multiply(sT1,-Ek)
       CALL Add(sLXk,sT1,sGk)
       CALL Multiply(sGk,Two)
       !
       IF(K==0)THEN
          Beta=Zero
       ELSE
          ! Beta=Pdot1(N,P%D,Gk,Gk-Gkold)/Pdot1(N,P%D,GkOld,GkOld)
          CALL Multiply(sGkOld,-One)
          CALL Add(sGk,sGkOld,sT1)
          CALL Multiply(sGkOld,-One)
          Num=OneDotT(sP,sT1,sGk)
          Den=OneDotT(sP,sGkOld,sGkOld)
          Beta=Num/Den
       ENDIF
       IF(K==0)THEN
          ! Pk=Gk
          CALL SetEq(sPk,sGk)
       ELSE
          ! Pk=Gk+Beta*PkOld
          CALL Multiply(sPkOld,Beta)
          CALL Add(sPkOld,sGk,sPk)
       ENDIF

       ! CALL ReNorm(N,P%D,Xk)
       Norm=SQRT(ABS(OneDotT(sP,sPk,sPk)))
       CALL Multiply(sPk,One/Norm)
       CALL Filter(sPk,Tol_O=QN%Threshold)
       PkPcnt=1D2*DBLE(sPk%NNon0)/DBLE(NBasF**2)
       ! LPk=MATMUL(F%D,Pk)-MATMUL(Pk,F%D)
       CALL Multiply(sF,sPk,sT1)
       CALL Multiply(sPk,sF,sT2)
       CALL Multiply(sT2,-One)
       CALL Add(sT1,sT2,sLPk)
       ! LPk=0.25D0*(MATMUL(MATMUL(Q%D,LPk),P%D))
       CALL Multiply(sQ,sLPk,sT1)
       CALL Multiply(sT1,sP,sLPk)
       CALL Multiply(sLPk,25D-2)
!!$
!!$
!!$
!!$       ! Filter the LPk cooresponding to the NORMALIZED gradient
!!$       CALL Filter(sLPk,Tol_O=QN%Threshold)
!!$       LPkPcnt=1D2*DBLE(sLPk%NNon0)/DBLE(NBasF**2)
!!$       ! Now, de-normalized the gradient and cooresponding LPk
!!$
       CALL Multiply(sPk,Norm)
       CALL Multiply(sLPk,Norm)
       !
       CALL SparseRQILineSearch(sP,sPk,sXk,sLXk,sLPk,Lambda)
       !
       EkOld=Ek
       CALL SetEq(sXkOld,sXk)
       CALL SetEq(sGkOld,sGk)
       CALL SetEq(sPkOld,sPk)
       CALL Multiply(sPk,Lambda)
       ! Xk=XkOld+Lambda*Pk
       CALL Add(sXkOld,sPk,sXk)
       CALL Multiply(sPk,One/Lambda)  ! Multiply back
       ! Xk=0.25D0*(MATMUL(MATMUL(Q%D,Xk),P%D))
       CALL Multiply(sQ,sXk,sT1)
       CALL Multiply(sT1,sP,sXk)
       CALL Multiply(sXk,25D-2)
       ! CALL ReNorm(N,P%D,Xk)
       Norm=SQRT(ABS(OneDotT(sP,sXk,sXk)))
       CALL Multiply(sXk,One/Norm)
       CALL Filter(sXk,Tol_O=QN%Threshold)
       XkPcnt=1D2*DBLE(sXk%NNon0)/DBLE(NBasF**2)
       !
       ! LXk=MATMUL(F%D,Xk)-MATMUL(Xk,F%D)
       CALL Multiply(sF,sXk,sT1)
       CALL Multiply(sXk,sF,sT2)
       CALL Multiply(sT2,-One)
       CALL Add(sT1,sT2,sLXk)

       ! LXk=0.25D0*(MATMUL(MATMUL(Q%D,LXk),P%D))
       CALL Multiply(sQ,sLXk,sT1)
       CALL Multiply(sT1,sP,sLXk)
       CALL Multiply(sLXk,25D-2)


       !       CALL Filter(sLXk,Tol_O=QN%Threshold)
       !       LXkPcnt=1D2*DBLE(sLXk%NNon0)/DBLE(NBasF**2)
       !
       ! Ek=ThoulessQ(N,P%D,Xk,LXk)
       Ek=OneDotT(sP,sXk,sLXk)
       ErrAbs=EkOld-Ek
       ErrRel=ErrAbs/Ek
       IF(Mod(K,2)==0)THEN
          EkOld=MIN(EkOld,Ek)
       ENDIF



!!$
!!$       CALL Multiply(sXk,-Ek)
!!$       CALL Add(sXk,sLXk,sT1)
!!$       CALL Multiply(sXk,-One/Ek) ! Multiply back
!!$       dEk=-1D10
!!$       DO L=1,sT1%NNon0
!!$          dEk=MAX(dEk,ABS(sT1%MTrix%D(L)))
!!$       ENDDO


       Iteration="Koopmans' "//TRIM(IntToChar(K))

       Statistics1="Ev = "//TRIM(DblToMedmChar(Ek*27.21139613182D0))   &
            //", Err = "//TRIM(DblToShrtChar(ErrRel)) &
            //", %X = "//TRIM(FltToShrtChar(XkPcnt))   &
            //", %LX = "//TRIM(FltToShrtChar(LXkPcnt)) &
            //", %P = "//TRIM(FltToShrtChar(PkPcnt))   &
            //", %LP = "//TRIM(FltToShrtChar(LPkPcnt))

       CALL MondoLog(DEBUG_NONE, "QUIRQI",Statistics1,Iteration)

       IF(PRESENT(RelErrTol_O))THEN
          IF(ErrRel<RelErrTol_O)THEN
             CALL Elapsed_TIME(TimeTOTAL,Init_O='Accum')
             CALL OpenASCII('k_scaling',77)
             WRITE(77,*)NAtoms,Ek*27.21139613182D0,TimeTOTAL%Wall,K,XkPcnt,PkPcnt
             CLOSE(UNIT=77)
             EXIT
          ENDIF
       ENDIF

       IF(K>4.AND.(ErrRel<QN%Konvergence.OR.ErrAbs<Zero))THEN
          CALL Elapsed_TIME(TimeTOTAL,Init_O='Accum')
          CALL OpenASCII('k_scaling',77)
          WRITE(77,*)NAtoms,Ek*27.21139613182D0,TimeTOTAL%Wall,K,XkPcnt,PkPcnt
          CLOSE(UNIT=77)
          !STOP
          EXIT
       ENDIF

    ENDDO


    !    DO I=1,sXk%NNon0
    !       MaxP=ABS(sXk%MTrix%D(I))*1D-1
    !       sXk%MTrix%D(I)=sXk%MTrix%D(I)+RANDOM_DBL((/-MaxP,MaxP/))
    !    ENDDO
!!$
111 CONTINUE
!!$
!!$    CALL New(X,(/NBasF,NBasF/))
!!$    CALL New(EigenEs,NBasF)
!!$    CALL SetEq(X,sXk)
!!$    CALL MDiag_DSYEVD(X,NBasF,EigenEs,0)
!!$
!!$    DO I=NBasF/2+1,NBasF
!!$       WRITE(*,*)I,EigenEs%D(I)
!!$    ENDDO
!!$
!!$



    CALL Put(sXk,QN%Xk%Ortho)


    CALL Delete(sQ)
    CALL Delete(sP)
    CALL Delete(sF)
    CALL Delete(sT1)
    CALL Delete(sT2)
    CALL Delete(sXk)
    CALL Delete(sGk)
    CALL Delete(sPk)
    CALL Delete(sLXk)
    CALL Delete(sLPk)
    CALL Delete(sXkOld)
    CALL Delete(sGkOld)
    CALL Delete(sPkOld)

    ! Split the guess
    CALL SplitPQ(QN,QN%Xk,.TRUE.)
    !

    !
  END SUBROUTINE Koopmans

  !QUIRQI                 ::  Ev = 0.60184723D+01, ErrRel = 0.45D+03, TrueEr = 0.63D+00, Delta = 0.10D+03, BetaP = 0.00D+00, BetaQ = 0.19-311, Tau = 0.10D-19,%X = 0.00D+00,%P = 0.28-316

  SUBROUTINE SetQName(Nam,S,QN,RQIStat)
    !
    TYPE(State)          :: S,RQIStat
    TYPE(QUIRQIKontrol)    :: QN
    TYPE(FileNames)      :: Nam
    LOGICAL              :: Present
    !
    CALL New(RQIStat%Action,2)
    CALL New(RQIStat%Current,3)
    CALL New(RQIStat%Previous,3)
    RQIStat%Current%I=S%Current%I
    RQIStat%Previous%I=S%Previous%I
    ! Action is TD-SCF with secondary parameter the product LX or LP (L[Xk], L[Pk])
    RQIStat%Action%C(1)="TD-SCF"
    RQIStat%Action%C(2)="" ! To be determined throughout RQI cycle
    !
    QN%P=          TrixFile("OrthoD",     PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%F=          TrixFile("F_DIIS",     PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    INQUIRE(FILE=QN%F,EXIST=Present)
    IF(.NOT.Present) &
         QN%F=          TrixFile("OrthoF",     PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%J=          TrixFile("J",          PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%K=          TrixFile("K",          PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    ! INQUIRE, Z Z^t or X:
    QN%Z=          TrixFile("X",          PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I)
    INQUIRE(FILE=QN%Z,EXIST=Present)
    IF(.NOT.Present)THEN
       QN%CholFact=.TRUE.
       QN%Z=          TrixFile("Z",           PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I)
       QN%ZT=         TrixFile("ZT",          PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I)
    ELSE
       QN%ZT=QN%Z
    ENDIF
    !
    QN%ASCII_Xk=TRIM(Nam%M_PWD)//TRIM(Nam%SCF_NAME)//'.'//'ASCII_Xk'
    !
    QN%Xk%Ortho=   TrixFile('OrthoXk',    PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Xk%AO=      TrixFile('Xk',         PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Xk%MomQ=    TrixFile('Xk_MomQ',    PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Xk%PosP=    TrixFile('Xk_PosP',    PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Xk%MomQ_Old=TrixFile('Xk_MomQ_Old',PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Xk%PosP_Old=TrixFile('Xk_PosP_Old',PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Xk%Act2='Xk'
    !
    QN%LXk%Ortho=  TrixFile('LXk',        PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%LXk%MomQ=   TrixFile('LXk_MomQ',   PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%LXk%PosP=   TrixFile('LXk_PosP',   PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    !
    QN%Gk%Ortho=   TrixFile('OrthoGk',    PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Gk%AO=      TrixFile('Gk',         PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Gk%MomQ=    TrixFile('Gk_MomQ',    PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Gk%PosP=    TrixFile('Gk_PosP',    PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Gk%MomQ_Old=TrixFile('Gk_MomQ_Old',PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Gk%PosP_Old=TrixFile('Gk_PosP_Old',PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Gk%Act2='Gk'
    !
    QN%Pk%Ortho=   TrixFile('OrthoPk',    PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Pk%AO=      TrixFile('Pk',         PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Pk%MomQ=    TrixFile('Pk_MomQ',    PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Pk%PosP=    TrixFile('Pk_PosP',    PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Pk%MomQ_Old=TrixFile('Pk_MomQ_Old',PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Pk%PosP_Old=TrixFile('Pk_PosP_Old',PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%Pk%Act2='Pk'
    !
    QN%LPk%Ortho=  TrixFile('LPk',        PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%LPk%MomQ=   TrixFile('LPk_MomQ',   PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QN%LPk%PosP=   TrixFile('LPk_PosP',   PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    !
  END SUBROUTINE SetQName



  SUBROUTINE QUIRQI(N,M,Nam,O,S,MPI,G,B)
    INTEGER            :: N,M,I,II,J,K,L,U,V,JTDA,cBAS
    TYPE(FileNames)    :: Nam,RQINams
    TYPE(State)        :: S,RQIStat
    TYPE(Parallel)     :: MPI
    TYPE(Options)      :: O
    TYPE(Geometries)   :: G
    TYPE(BSET)         :: B
    LOGICAL            :: DoTDA
    !
    TYPE(BCSR)                      :: sP,sQ,sF,sZ,sXk,sCPSCF
    TYPE(DBL_VECT)                  :: EigenEs
    TYPE(DBL_RNK2)                  :: P,Q,F,Z,Orbitals,Xk_Guess,CPSCF
    REAL(DOUBLE)                    :: Ek,EkOld,dEk,Beta,Lambda,ErrRel,ErrAbs,Shift,dNorm, &
         Tmp_EOM_P,Tmp_EOM_Q,Tmp_EOM,LamP,LamQ,EkT
    INTEGER                         :: XkNon0s,PkNon0s,KK,JJ,KStop
    REAL(DOUBLE)                    :: XkThreshold,PkThreshold,PMax,NTr,DTr,E_EOM,E_EOM_P,E_EOM_Q,Lambda_PosP,Lambda_MomQ,Lam,Err,E2

    REAL(DOUBLE),DIMENSION(N,N)     :: Xk,Gk,Pk,LXk,LPk,LQk,PositionP,MomentumQ,XXk,GkTmp,GkTmpOld
    REAL(DOUBLE),DIMENSION(N,N)     :: XkOld,PkOld,GkOld,GPk,GQk,LX_1E,LX_2E,LP_1E,LP_2E
    REAL(DOUBLE),DIMENSION(N)       :: Values
    REAL(DOUBLE),DIMENSION(N,N,M)   :: Vectors
    REAL(DOUBLE),DIMENSION(N,N,N,N) :: TwoE,DoubleSlash



    REAL(DOUBLE), DIMENSION(N,N)    :: PosP,MomQ,PosP_Xk,MomQ_Xk, PosP_Pk,MomQ_Pk,L_PosP_Xk, &
         L_MomQ_Xk,L_PosP_Pk,L_MomQ_Pk,PosP_Gk,MomQ_Gk, &
         PosP_GkOld,MomQ_GkOld,PosP_PkOld,MomQ_PkOld,LGk, &
         PosP_XkOld,MomQ_XkOld

    REAL(DOUBLE) :: Ap,Aq,Bp,Bq,Cp,Cq,Rpq,Spq,Tpq,Upq,TTmp_EOM,PosP_Beta,MomQ_Beta,PosP_Ek,MomQ_Ek,Norm,Tmp,EDiff1,EDiff2,ASwitch
    !
    RQIStat=S
    RQINams=Nam
    cBAS=S%Current%I(2)
    !
    CALL Elapsed_TIME(TimeTotal,Init_O='Init')
    CALL Elapsed_TIME(TimeONX,Init_O='Init')
    CALL Elapsed_TIME(TimeQCTC,Init_O='Init')
    CALL Elapsed_TIME(TimeBCSR,Init_O='Init')
    CALL Elapsed_TIME(TimeTotal,Init_O='Start')
    !
    CALL Get(sP,TrixFile("OrthoD",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0))
    CALL Get(sF,TrixFile("OrthoF",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0))
    CALL Get(sZ,TrixFile("X",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I))

    CALL Get(sCPSCF,'/scratch/mchalla/MONDO_SCRATCH/cpscf.guess')

    CALL SetEq(sQ,sP)
    !
    CALL SetEq(P,sP)
    CALL SetEq(Q,sP)
    CALL SetEq(F,sF)
    CALL SetEq(Z,sZ)
    CALL SetEq(CPSCF,sCPSCF)
    !
    CALL Multiply(sP,Two)
    CALL Multiply(sQ,-Two)
    CALL Add(sQ,Two)
    !

    P%D=Two*P%D
    Q%D=-Two*Q%D
    DO I=1,N
       Q%D(I,I)=Q%D(I,I)+Two
    ENDDO
    !
    !
    !
    ! BEGIN LOGIC TO DO CALCS ALL IN MO REPRESENTATION
    !
    CALL New(Orbitals,(/NBasF,NBasF/))
    Orbitals%D=F%D
    !    CALL SetEq(Orbitals,F)
    CALL New(EigenEs,NBasF)
    CALL MDiag_DSYEVD(Orbitals,NBasF,EigenEs,0)

    CALL F2Rot(NBasF,Orbitals%D,P%D)
    CALL F2Rot(NBasF,Orbitals%D,Q%D)
    CALL F2Rot(NBasF,Orbitals%D,F%D)
    CALL F2Rot(NBasF,Orbitals%D,CPSCF%D)

    CALL Integrals2E(B,G%Clone(1),TwoE)
    DO I=1,N
       DO J=1,N
          DO K=1,N
             DO L=1,N
                DoubleSlash(I,J,K,L)=TwoE(I,J,K,L)-TwoE(I,K,J,L)/2D0
             ENDDO
          ENDDO
       ENDDO
    ENDDO
!!$
    ! First trapose AO->OR
    CALL F4Rot(NBasF,Z%D,DoubleSlash)
    ! Then OR->MO
    CALL F4Rot(NBasF,Orbitals%D,DoubleSlash)
    ! Now set Z->I, avoiding any further AO->OR action
    Z%D=0D0
    DO I=1,NBasF; Z%D(I,I)=1D0; ENDDO
       !
!!$       CALL PPrint(F,'F_MO',Unit_O=6)
!!$       CALL PPrint(P,'P_MO',Unit_O=6)
!!$       CALL PPrint(Q,'Q_MO',Unit_O=6)
       CALL PPrint(CPSCF,'CPSCF_MO',Unit_O=6)
!!$       CALL KoopmansGuess(N,P,Q,F,Xk)

!!$       CALL PPrint(Xk,'KOOP_MO',Unit_O=6)

       !
       ! END LOGIC TO DO CALCS ALL IN MO REPRESENTATION
       !

       !

33     FORMAT('St=',I2,', It=',I3,', Ev=',F10.6,', dE=',D8.2,', dN=',D7.2, &
            ', Tk=',D10.4,', Tj=',D10.4,', Tm=',D10.4,', |Gk|=',D8.2,', %Xk=',F6.2,', %Pk=',F6.2)


44     FORMAT(A4,' State=',I2,', Nk=',I3,', Ev=',F9.6,', dE=',D7.2,', WallSec=',D12.4)


       WRITE(*,*)'QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ'

       DO I=0,0
          !
          DO JTDA=1,1
             IF(JTDA==0)THEN
                DoTDA=.TRUE.
                CALL KoopmansGuess(N,P,Q,F,Xk)
             ELSEIF(JTDA==1)THEN
                DoTDA=.FALSE.
                CALL KoopmansGuess(N,P,Q,F,Xk)
             ENDIF


             PosP_Xk=pOp(NBasF,P%D,Q%D,Xk)
             MomQ_Xk=qOp(NBasF,P%D,Q%D,Xk)


             !    CALL PPrint(PosP_Xk,'PosPPP',Unit_O=6)
             !    CALL PPrint(MomQ_Xk,'MomQQQ',Unit_O=6)


             !

             CALL LOn2(N,I,Shift,F%D,P%D,Z%D,DoubleSlash,Values,Vectors,Xk,LXk)
             CALL Anihilate(N,P%D,Q%D,LXk,TDA_O=.FALSE. ) !=DoTDA)
             L_PosP_Xk=qOp(NBasF,P%D,Q%D,LXk)
             L_MomQ_Xk=pOp(NBasF,P%D,Q%D,LXk)

             !    CALL PChkSum(L_PosP_Xk,'L_PosP',Unit_O=6)
             !    CALL PChkSum(L_MomQ_Xk,'L_MomQ',Unit_O=6)
!!$

             PosP_XkOld=Zero
             MomQ_XkOld=Zero
             PosP_GkOld=Zero
             MomQ_GkOld=Zero
             PosP_PkOld=Zero
             MomQ_PkOld=Zero
             !
             ErrAbs=1D10
             ErrRel=1D10
             Ek=ThoulessQ(N,P%D,Xk,LXk)
             PosP_Ek=Trace_EOM(N,PosP_Xk,L_PosP_Xk)/Trace_EOM(N,PosP_Xk,MomQ_Xk)
             MomQ_Ek=Trace_EOM(N,MomQ_Xk,L_MomQ_Xk)/Trace_EOM(N,PosP_Xk,MomQ_Xk)

             WRITE(*,*)' Ek = ',Ek !*27.21139613182D0

             !          RETURN

             DO K=0,10
                !------------------------------------------------------------------------------
                ! QUASI-INDEPENDENT GRADIENTS OF THE EOM FUNCTIONAL
                !------------------------------------------------------------------------------
                PosP_Gk=-L_PosP_Xk+MomQ_Ek*MomQ_Xk
                MomQ_Gk=-L_MomQ_Xk+PosP_Ek*PosP_Xk
                !------------------------------------------------------------------------------
                ! QUASI-INDEPENDENT POLAK-RIBIERE NLCG UPDATE CONSTANTS
                !------------------------------------------------------------------------------
                IF(K==0)THEN
                   PosP_Beta=Zero
                   MomQ_Beta=Zero
                ELSE
                   PosP_Beta=Trace_EOM(N,PosP_Gk,PosP_Gk-PosP_Gkold)/Trace_EOM(N,PosP_GkOld,PosP_GkOld)
                   MomQ_Beta=Trace_EOM(N,MomQ_Gk,MomQ_Gk-MomQ_Gkold)/Trace_EOM(N,MomQ_GkOld,MomQ_GkOld)
                   PosP_Beta=MAX(Zero,PosP_Beta)
                   MomQ_Beta=MAX(Zero,MomQ_Beta)
                ENDIF
                !------------------------------------------------------------------------------
                ! QUASI-INDEPENDENT UPDATE OF THE CONJUGATE GRADIENTS
                !------------------------------------------------------------------------------
                IF(K==0)THEN
                   PosP_Pk=PosP_Gk
                   MomQ_Pk=MomQ_Gk
                ELSE
                   PosP_Pk=PosP_Gk+PosP_Beta*PosP_PkOld
                   MomQ_Pk=MomQ_Gk+MomQ_Beta*MomQ_PkOld
                ENDIF
                !------------------------------------------------------------------------------
                ! L ONTO Pk, THEN SPLIT
                !------------------------------------------------------------------------------
                Pk=Half*(PosP_Pk+MomQ_Pk+TRANSPOSE(PosP_Pk-MomQ_Pk))
                CALL LOn2(N,I,Shift,F%D,P%D,Z%D,DoubleSlash,Values,Vectors,Pk,LPk)
                CALL Anihilate(N,P%D,Q%D,LPk,TDA_O=DoTDA)
                L_PosP_Pk=qOp(NBasF,P%D,Q%D,LPk)
                L_MomQ_Pk=pOp(NBasF,P%D,Q%D,LPk)
!!$
!!$             CALL PChkSum(L_PosP_Pk,'PosP_PK',Unit_O=6)
!!$             CALL PChkSum(L_MomQ_Pk,'MomQ_PK',Unit_O=6)

                !------------------------------------------------------------------------------
                ! QUASI-INDEPENDENT LINE SEARCH OF EOM FUNCTIONAL
                !------------------------------------------------------------------------------
                CALL EOMLineSearch(N,P%D,Q%D,F%D,Z%D,DoubleSlash,DoTDA, &
                     PosP_Xk,MomQ_Xk,L_PosP_Xk,L_MomQ_Xk, &
                     PosP_Pk,MomQ_Pk,L_PosP_Pk,L_MomQ_Pk, &
                     LamP,LamQ,E_EOM)
                !------------------------------------------------------------------------------
                ! SET OLD VARIABLES FOR NLCG ETC
                !------------------------------------------------------------------------------
                EkOld=Ek
                PosP_XkOld=PosP_Xk
                MomQ_XkOld=MomQ_Xk
                PosP_GkOld=PosP_Gk
                MomQ_GkOld=MomQ_Gk
                PosP_PkOld=PosP_Pk
                MomQ_PkOld=MomQ_Pk
                !------------------------------------------------------------------------------
                ! APPLY CONJUGATE GRADIENTS
                !------------------------------------------------------------------------------
                PosP_Xk=PosP_Xk+LamP*PosP_Pk
                MomQ_Xk=MomQ_Xk+LamQ*MomQ_Pk
                !------------------------------------------------------------------------------
                ! RENORM
                !------------------------------------------------------------------------------
                Norm=SQRT(ABS(Trace_EOM(N,PosP_Xk,MomQ_Xk)))
                PosP_Xk=PosP_Xk/Norm
                MomQ_Xk=MomQ_Xk/Norm
                Xk=Half*(PosP_Xk+MomQ_Xk+TRANSPOSE(PosP_Xk-MomQ_Xk))
                !------------------------------------------------------------------------------
                ! L ONTO Xk
                !------------------------------------------------------------------------------
                CALL LOn2(N,I,Shift,F%D,P%D,Z%D,DoubleSlash,Values,Vectors,Xk,LXk)
                CALL Anihilate(N,P%D,Q%D,LXk,TDA_O=DoTDA)
                EkT=ThoulessQ(N,P%D,Xk,LXk)
                L_PosP_Xk=qOp(NBasF,P%D,Q%D,LXk)
                L_MomQ_Xk=pOp(NBasF,P%D,Q%D,LXk)
!!$
!!$             CALL PChkSum(L_PosP_Xk,'PosP_XK',Unit_O=6)
!!$             CALL PChkSum(L_MomQ_Xk,'MomQ_XK',Unit_O=6)
!!$             RETURN

                !------------------------------------------------------------------------------
                ! ENERGIES AND ERRORS
                !------------------------------------------------------------------------------
                PosP_Ek=Trace_EOM(N,PosP_Xk,L_PosP_Xk)/Trace_EOM(N,PosP_Xk,MomQ_Xk)
                MomQ_Ek=Trace_EOM(N,MomQ_Xk,L_MomQ_Xk)/Trace_EOM(N,PosP_Xk,MomQ_Xk)
                Ek=Half*(PosP_Ek+MomQ_Ek)
                PosP_Xk=pOp(NBasF,P%D,Q%D,Xk)
                MomQ_Xk=qOp(NBasF,P%D,Q%D,Xk)

                ErrAbs=Ek-EkOld
                ErrRel=-1D10
                DO U=1,N
                   DO V=1,N
                      ErrRel=MAX(ErrRel,ABS(Ek*Xk(U,V)-LXk(U,V))/Ek)
                   ENDDO
                ENDDO

                IF(JTDA==-1)THEN
                   Err=LOG10(ABS(Ek-0.663332701870816D0)+1D-20)
                   WRITE(*,*)K,Err,ErrRel,EkT
#ifdef KWERQY
                   WRITE(46,*)K,Err
#else
                   WRITE(36,*)K,Err
#endif
                ELSEIF(JTDA==0)THEN
                   Err=LOG10(ABS(Ek-0.432987734557198D0)+1D-20)
                   WRITE(*,*)K,Err,ErrRel,EkT
#ifdef KWERQY
                   WRITE(45,*)K,Err
#else
                   WRITE(35,*)K,Err
#endif
                ELSEIF(JTDA==1)THEN
                   Err=LOG10(ABS(Ek-0.406088110418702D0)+1D-20)
                   WRITE(*,*)K,Err,ErrRel,Ek*27.21139613182D0
                   WRITE(34,*)K,Err
                ELSE
                   WRITE(*,*)K,ErrRel,Ek*27.21139613182D0
                ENDIF

             ENDDO


             CALL PPrint(Xk,'Xk',Unit_O=6)
             !          CALL PPrint(MomQ_Xk,'MomQ',Unit_O=6)



             RETURN

             IF(JTDA==-1)THEN
                Xk_Guess%D=Xk
             ENDIF

             Values(I)=Ek
             Vectors(:,:,I)=Xk
          ENDDO
       ENDDO

     END SUBROUTINE QUIRQI



     SUBROUTINE KoopmansGuess(N,P,Q,F,Xk,DoTDA,E2)
       INTEGER            :: N,M,I,II,J,K,L,U,V,JTDA,cBAS
       TYPE(DBL_RNK2)                  :: P,Q,F
       REAL(DOUBLE)                    :: Ek,EkOld,dEk,Beta,Lambda,ErrRel,ErrAbs
       INTEGER                         :: XkNon0s,PkNon0s,KK,JJ,KStop
       REAL(DOUBLE)                    :: XkThreshold,PkThreshold,PMax,NTr,DTr,E_EOM,E_EOM_P,E_EOM_Q,Lambda_PosP,Lambda_MomQ,Lam,Err
       REAL(DOUBLE), OPTIONAL :: E2

       REAL(DOUBLE),DIMENSION(N,N)     :: Xk,Gk,Pk,LXk,LPk
       REAL(DOUBLE),DIMENSION(N,N)     :: XkOld,PkOld,GkOld
       LOGICAL, OPTIONAL               :: DoTDA
       INTEGER :: KCycle
       !

       !    WRITE(*,*)' BEG KOOP BEG KOOP BEG KOOP BEG KOOP BEG KOOP '
       IF(.NOT.PRESENT(E2))THEN
          DO I=1,N
             DO J=1,N
                Xk(I,J)=I+J
             ENDDO
          ENDDO
       ENDIF
       !
       LXk=MATMUL(F%D,Xk)-MATMUL(Xk,F%D)
       !
       IF(.NOT.PRESENT(DoTDA))THEN
          LXk=0.25D0*(MATMUL(MATMUL(Q%D,LXk),P%D))
       ELSE
          LXk=0.25D0*( MATMUL(MATMUL(Q%D,LXk),P%D)+MATMUL(MATMUL(P%D,LXk),Q%D))
       ENDIF
       !
       XkOld=Zero
       PkOld=Zero
       GkOld=Zero
       IF(PRESENT(E2))THEN
          Ek=ThoulessQ(N,P%D,Xk,LXk)
       ELSE
          Ek=ThoulessQ(N,P%D,Xk,LXk)
       ENDIF

       !    WRITE(*,*) ' EK KOOP = ',Ek

       IF(PRESENT(E2))THEN
          KCycle=5
       ELSE
          KCycle=100
       END IF

       DO K=0,KCycle

          Gk=Two*(LXk-Ek*Xk)

          IF(PRESENT(E2))THEN
             Gk=0.25D0*(MATMUL(MATMUL(Q%D,Gk),P%D)+MATMUL(MATMUL(P%D,Gk),Q%D))
             CALL PPrint(Gk,'Gk',Unit_O=6)
             STOP
          ENDIF

          IF(K==0)THEN
             Beta=Zero
          ELSE
             Beta=Pdot1(N,P%D,Gk,Gk-Gkold)/Pdot1(N,P%D,GkOld,GkOld)
          ENDIF

          IF(K==0)THEN
             Pk=Gk
          ELSE
             Pk=Gk+Beta*PkOld
          ENDIF
          !
          LPk=MATMUL(F%D,Pk)-MATMUL(Pk,F%D)
          !
          IF(.NOT.PRESENT(DoTDA))THEN
             LPk=0.25D0*(MATMUL(MATMUL(Q%D,LPk),P%D))
          ELSE
             LPk=0.25D0*(MATMUL(MATMUL(Q%D,LPk),P%D)+MATMUL(MATMUL(P%D,LPk),Q%D))
          ENDIF
          !
          IF(PRESENT(E2))THEN
             CALL RQILineSearch2E(N,P%D,F%D,Pk,Xk,E2,Lambda)

             !          CALL RQILineSearch(N,P%D,Pk,Xk,LXk,LPk,Lambda,E2)
          ELSE
             CALL RQILineSearch(N,P%D,Pk,Xk,LXk,LPk,Lambda)
          ENDIF

          EkOld=Ek
          XkOld=Xk
          EkOld=Ek
          GkOld=Gk
          PkOld=Pk
          !
          Xk=XkOld+Lambda*Pk

          !
          IF(.NOT.PRESENT(DoTDA))THEN
             Xk=0.25D0*(MATMUL(MATMUL(Q%D,Xk),P%D))
          ELSE
             Xk=0.25D0*( MATMUL(MATMUL(Q%D,Xk),P%D)+MATMUL(MATMUL(P%D,Xk),Q%D))
          ENDIF
          !


          CALL ReNorm(N,P%D,Xk)


          !
          LXk=MATMUL(F%D,Xk)-MATMUL(Xk,F%D)

          !      IF(PRESENT(E2)) &
          !      WRITE(*,*)' AFTER UPDATE 1 IN KOOP, EK = ',ThoulessQ(N,P%D,Xk,LXk)+E2


          IF(.NOT.PRESENT(DoTDA))THEN
             LXk=0.25D0*(MATMUL(MATMUL(Q%D,LXk),P%D))
          ELSE
             LXk=0.25D0*( MATMUL(MATMUL(Q%D,LXk),P%D)+MATMUL(MATMUL(P%D,LXk),Q%D))
          ENDIF

          IF(PRESENT(E2))THEN
             Ek=ThoulessQ(N,P%D,Xk,LXk)+E2
          ELSE
             Ek=ThoulessQ(N,P%D,Xk,LXk)
          ENDIF

          ErrAbs=Ek-EkOld
          ErrRel=-1D10
          DO U=1,N
             DO V=1,N
                ErrRel=MAX(ErrRel,ABS(Ek*Xk(U,V)-LXk(U,V))/Ek)
             ENDDO
          ENDDO


          Err=LOG10(ABS(Ek-0.663332701870816D0)+1D-20)

          !       WRITE(*,*)K,Ek,Err !,ErrRel
          !       WRITE(46,*)K,Err
       ENDDO

       !    WRITE(*,*)' END KOOP END KOOP END KOOP END KOOP END KOOP '
     END SUBROUTINE KOOPMANSGUESS


     Subroutine RQILineSearch2E(N,P,F,Pk,Xk,E2,Lambda)
       INTEGER :: N,I
       REAL(DOUBLE) :: L,Lambda,Lambda_p,Lambda_m,EkMin,Ek,E2
       REAL (DOUBLE),DIMENSION(N,N)::P,F,Pk,Xk,Tmp1
       REAL(DOUBLE) :: XFX,PFX,XFP,PFP,XXF,XPF,PXF,PPF, &
            XXL2,XPL2,PXL2,PPL2,XX,PP,XP,PX,AA,BB,CC

       XX =Pdot1(N,P,Xk,Xk) ! 1 - normalized by definition
       PP =Pdot1(N,P,Pk,Pk)
       XP =Pdot1(N,P,Xk,Pk)
       PX =Pdot1(N,P,Pk,Xk)

       XFX=Pdot1(N,P,Xk,MATMUL(F,Xk))
       PFX=Pdot1(N,P,Pk,MATMUL(F,Xk))
       XFP=Pdot1(N,P,Xk,MATMUL(F,Pk))
       PFP=Pdot1(N,P,Pk,MATMUL(F,Pk))

       XXF=Pdot1(N,P,Xk,MATMUL(Xk,F))
       XPF=Pdot1(N,P,Xk,MATMUL(Pk,F))
       PXF=Pdot1(N,P,Pk,MATMUL(Xk,F))
       PPF=Pdot1(N,P,Pk,MATMUL(Pk,F))

       !    AA=XFX
       !    BB=(PX+XP)*(PF-FP+PLX)-Two*PP !2.0*PLP*XX-2.0*XLX*PP
       !    CC=PP*(PF-FP+PLX)             !XX*(PLX+XLP)-XLX*(PX+XP)
!!$
!!$    lambda_p=(-BB+SQRT(BB*BB-4.0*AA*CC))/(2.0*AA)
!!$    lambda_m=(-BB-SQRT(BB*BB-4.0*AA*CC))/(2.0*AA)
!!$    Lambda=Lambda_P

       Ek=XFX-XXF+E2
       WRITE(*,*)'Before, E1 = ',XFX-XXF,' E2 = ',E2


       EkMin=1D5
       DO I=-10,10
          L=DBLE(I)*1D-2
          Ek=( (XFX+L*(PFX+XFP)+L**2*PFP) &
               -(XXF+L*(PXF+XPF)+L**2*PPF) &
               +E2) &
               /(XX+L*(PX+XP)+L**2*PP)
          !       WRITE(*,*)L,Ek,(XX+L*(PX+XP)+L**2*PP)
          IF(Ek<EkMin)THEN
             Lambda=L
             EkMin=Ek

          ENDIF
       ENDDO

       WRITE(*,*)'After ',Lambda,EkMin
       !    STOP

       RETURN
!!$!    WRITE(*,*)lambda_p,' AFTER EK P= ',(XLX+lambda_p*XLP+lambda_p*PLX+lambda_p**2 * PLP) /(XX+lambda_p*XP+lambda_p*PX+lambda_p**2 * PP)
!!$!    WRITE(*,*)lambda_q,' AFTER EK M= ',(XLX+lambda_m*XLP+lambda_m*PLX+lambda_m**2 * PLP) /(XX+lambda_m*XP+lambda_m*PX+lambda_m**2 * PP)
!!$
!!$    WRITE(*,*)' Lambda p/m = ',Lambda_P,Lambda_
!!$    WRITE(*,*)' Den = ',XX+Lambda*(PX+XP)+Lambda**2 *PP
!!$    WRITE(*,*)' Num = ',XF-FX+XLX+Lambda_P*(PF-FP+PLX)
!!$    WRITE(*,*)' Num = ',XF-FX+XLX+Lambda_M*(PF-FP+PLX)

     END Subroutine RQILineSearch2E


     Subroutine RQILineSearch(N,P,Pk,Xk,LXk,LPk,Lambda,E2)
       INTEGER :: N
       REAL(DOUBLE) :: Lambda,Lambda_p,Lambda_m
       REAL (DOUBLE),DIMENSION(N,N)::P,Pk,Xk,LXk,LPk,Tmp1
       REAL(DOUBLE) :: XX,PP,XP,PX,PLP,XLP,XLX,PLX,AA,BB,CC
       REAL(DOUBLE), OPTIONAL :: E2


       XX =Pdot1(N,P,Xk,Xk) ! 1 - normalized by definition
       PP =Pdot1(N,P,Pk,Pk)
       XP =Pdot1(N,P,Xk,Pk)
       PX =Pdot1(N,P,Pk,Xk)
       PLP=Pdot1(N,P,Pk,LPk)
       XLX=Pdot1(N,P,Xk,LXk)
       XLP=Pdot1(N,P,Xk,LPk)
       PLX=Pdot1(N,P,Pk,LXk)

       IF(PRESENT(E2))THEN
          AA=PLP*(PX+XP)-PP*(PLX+XLP)+E2
       ELSE
          AA=PLP*(PX+XP)-PP*(PLX+XLP)
       ENDIF
       BB=2.0*PLP*XX-2.0*XLX*PP
       CC=XX*(PLX+XLP)-XLX*(PX+XP)
       lambda_p=(-BB+SQRT(BB*BB-4.0*AA*CC))/(2.0*AA)
       lambda_m=(-BB-SQRT(BB*BB-4.0*AA*CC))/(2.0*AA)
       Lambda=Lambda_P

       !    WRITE(*,*)'L = ',Lambda,' AFTER EK P= ',(XLX+lambda_p*XLP+lambda_p*PLX+lambda_p**2 * PLP) /(XX+lambda_p*XP+lambda_p*PX+lambda_p**2 * PP)
       !    WRITE(*,*)' AFTER EK M= ',(XLX+lambda_m*XLP+lambda_m*PLX+lambda_m**2 * PLP) /(XX+lambda_m*XP+lambda_m*PX+lambda_m**2 * PP)


     END Subroutine RQILineSearch



     SUBROUTINE EOMLineSearch(N,P,Q,F,Z,DoubleSlash,DoTDA, &
          PosP_Xk,MomQ_Xk,L_PosP_Xk,L_MomQ_Xk, &
          PosP_Pk,MomQ_Pk,L_PosP_Pk,L_MomQ_Pk, &
          Lambda_PosP,Lambda_MomQ,E_EOM)
       INTEGER                      :: N,II,JJ,KK,I
       LOGICAL                      :: DoTDA
       REAL(DOUBLE), DIMENSION(N,N) :: P,Q,F,Z,PosP_Xk,MomQ_Xk,L_PosP_Xk,L_MomQ_Xk, &
            PosP_Pk,MomQ_Pk,L_PosP_Pk,L_MomQ_Pk, &
            Xk,PosP,MomQ,LXk
       REAL(DOUBLE), DIMENSION(N,N,N,N) :: DoubleSlash
       REAL(DOUBLE)                 :: Lambda_PosP,Lambda_MomQ,E_EOM,Tmp_EOM,Ap,Aq, &
            Bp,Bq,Cp,Cq,Rpq,Spq,Tpq,Upq,LamP,LamQ,PosP_MIN,MomQ_MIN,dLamP,dLamQ,LamP_Old,LamQ_Old, &
            dEOM,EOM_Old


       Ap=Trace_EOM(N,PosP_Xk,L_PosP_Xk)
       Aq=Trace_EOM(N,MomQ_Xk,L_MomQ_Xk)
       Bp=Trace_EOM(N,PosP_Pk,L_PosP_Xk)+Trace_EOM(N,PosP_Xk,L_PosP_Pk)
       Bq=Trace_EOM(N,MomQ_Pk,L_MomQ_Xk)+Trace_EOM(N,MomQ_Xk,L_MomQ_Pk)
       Cp=Trace_EOM(N,PosP_Pk,L_PosP_Pk)
       Cq=Trace_EOM(N,MomQ_Pk,L_MomQ_Pk)
       Rpq=Trace_EOM(N,PosP_Xk,MomQ_Xk)
       Spq=Trace_EOM(N,PosP_Pk,MomQ_Xk)
       Tpq=Trace_EOM(N,PosP_Xk,MomQ_Pk)
       Upq=Trace_EOM(N,PosP_Pk,MomQ_Pk)


       !    WRITE(*,*)' Dot1 = ',Trace_EOM(N,PosP_Pk,L_PosP_Xk)
       !    WRITE(*,*)' Dot2 = ',Trace_EOM(N,PosP_Xk,L_PosP_Pk)
       !    WRITE(*,*)' Bp = ',Bp

       !    CALL PPrint(PosP_Pk,' EOM PosPk ',Unit_O=6)
       !    CALL PPrint(L_PosP_Pk,' EOM L_PosPk ',Unit_O=6)
       !    STOP


!!$
!!$    WRITE(*,*)' Ap = ',Ap
!!$    WRITE(*,*)' Aq = ',Aq
!!$    WRITE(*,*)' Bp = ',Bp
!!$    WRITE(*,*)' Bq = ',Bq
!!$    WRITE(*,*)' Cp = ',Cp
!!$    WRITE(*,*)' Cq = ',Cq
!!$    WRITE(*,*)' Rpq = ',Rpq
!!$    WRITE(*,*)' Spq = ',Spq
!!$    WRITE(*,*)' Tpq = ',Tpq
!!$    WRITE(*,*)' Upq = ',Upq

       E_EOM=Half*(Ap+Aq)/Rpq
       !     WRITE(*,44)E_EOM,Spq,Tpq,Upq
44     FORMAT(' E_EOM TO START = ',F20.10,' Spq = ',D20.10,' Tpq = ',D20.10,' Upq = ',D20.10)

       PosP_MIN=Zero
       MomQ_MIN=Zero
       !
       DO JJ=-100,100
          LamP=JJ*0.1D0
          DO KK=-100,100
             LamQ=KK*0.1D0
             Tmp_EOM=Half*(Ap+Bp*LamP+Cp*LamP**2+Aq+Bq*LamQ+Cq*LamQ**2)/(Rpq+LamP*Spq+LamQ*Tpq+LamP*LamQ*Upq)
             IF(Tmp_EOM<E_EOM.AND.Tmp_EOM>0D0)THEN
                E_EOM=Tmp_EOM
                PosP_MIN=LamP
                MomQ_MIN=LamQ
             ENDIF
          ENDDO
       ENDDO

       LamP=PosP_MIN
       LamQ=MomQ_MIN

       !    WRITE(*,46)E_EOM,LamP,LamQ
46     FORMAT(' E_EOM AT PT  A = ',F20.10,' LapP = ',D20.10,' LamQ = ',D20.10)

       DO I=1,5000
          LamQ_Old=LamQ
          LamQ=(-2D0*Cq*Rpq - 2D0*Cq*LamP*Spq +  Sqrt((2D0*Cq*Rpq + 2D0*Cq*LamP*Spq)**2 - &
               4D0*(Cq*Tpq + Cq*LamP*Upq)*(Bq*Rpq + Bq*LamP*Spq - Ap*Tpq - Aq*Tpq - Bp*LamP*Tpq - Cp*LamP**2*Tpq - Ap*LamP*Upq - &
               Aq*LamP*Upq - Bp*LamP**2*Upq - Cp*LamP**3*Upq)))/(2D0*(Cq*Tpq + Cq*LamP*Upq))
          LamP_Old=LamP
          LamP= (-2D0*Cp*Rpq - 2D0*Cp*LamQ*Tpq + Sqrt((2D0*Cp*Rpq + 2D0*Cp*LamQ*Tpq)**2 -          &
               4d0*(Cp*Spq + Cp*LamQ*Upq)*(Bp*Rpq - Ap*Spq - Aq*Spq - Bq*LamQ*Spq - Cq*LamQ**2*Spq + Bp*LamQ*Tpq - Ap*LamQ*Upq - &
               Aq*LamQ*Upq - Bq*LamQ**2*Upq - Cq*LamQ**3*Upq)))/(2D0*(Cp*Spq + Cp*LamQ*Upq))

          EOM_Old=E_EOM
          E_EOM=Half*(Ap+Bp*LamP+Cp*LamP**2+Aq+Bq*LamQ+Cq*LamQ**2)/(Rpq+LamP*Spq+LamQ*Tpq+LamP*LamQ*Upq)
          dLamQ=ABS((LamQ_Old-LamQ)/LamQ)
          dLamP=ABS((LamP_Old-LamP)/LamP)
          dEOM=ABS((EOM_Old-E_EOM)/E_EOM)
          IF(dEOM<1D-10.AND.I>5)EXIT
       ENDDO
       Lambda_PosP=LamP
       Lambda_MomQ=LamQ

       IF( (I==5001).OR.(ABS(LamP)<1D-8).OR.(ABS(LamQ)<1D-8) )THEN
          WRITE(*,48)E_EOM,LamP,LamQ,I

          Lambda_PosP=-Bp/(2D0*Cp)
          Lambda_MomQ=-Bq/(2D0*Cq)

          E_EOM=Half*(Ap+Bp*LamP+Cp*LamP**2+Aq+Bq*LamQ+Cq*LamQ**2)/(Rpq+LamP*Spq+LamQ*Tpq+LamP*LamQ*Upq)
          !       WRITE(*,*)' Bp = ',Bp,' Cp = ',Cp
          !       WRITE(*,*)' Bq = ',Bq,' Cq = ',Cq
          !       WRITE(*,48)E_EOM,LamP,LamQ,I
48        FORMAT(' E_EOM AT PT  F = ',F20.10,' LapP = ',D20.10,' LamQ = ',D20.10,' ICONVERGE = ',I8)


       ENDIF




       !!             IF(I==5001)THEN
!!$                PrintFlags%Fmt=DEBUG_MMASTYLE
!!$                CALL Print_DBL_SCLR(Ap,'Ap',Unit_O=6)
!!$                CALL Print_DBL_SCLR(Aq,'Aq',Unit_O=6)
!!$                CALL Print_DBL_SCLR(Bp,'Bp',Unit_O=6)
!!$                CALL Print_DBL_SCLR(Bq,'Bq',Unit_O=6)
!!$                CALL Print_DBL_SCLR(Cp,'Cp',Unit_O=6)
!!$                CALL Print_DBL_SCLR(Cq,'Cq',Unit_O=6)
!!$                CALL Print_DBL_SCLR(Rpq,'Rpq',Unit_O=6)
!!$                CALL Print_DBL_SCLR(Spq,'Spq',Unit_O=6)
!!$                CALL Print_DBL_SCLR(Tpq,'Tpq',Unit_O=6)
!!$                CALL Print_DBL_SCLR(Upq,'Upq',Unit_O=6)
       !!                E_EOM=Half*(Ap+Bp*LamP+Cp*LamP**2+Aq+Bq*LamQ+Cq*LamQ**2)/(Rpq+LamP*Spq+LamQ*Tpq+LamP*LamQ*Upq)
       !!                WRITE(*,48)E_EOM,LamP,LamQ,I
       !!48              FORMAT(' E_EOM AT PT  F = ',F20.10,' LapP = ',D20.10,' LamQ = ',D20.10,' ICONVERGE = ',I8)
!!$                STOP
       !!             ENDIF



       !    WRITE(*,47)E_EOM,LamP,LamQ,I

47     FORMAT(' E_EOM AT PT  B = ',F20.10,' LapP = ',D20.10,' LamQ = ',D20.10,' ICONVERGE = ',I8)


       RETURN

     END SUBROUTINE EOMLineSearch

     FUNCTION EOMThouless(N,P,Q,F,Z,DoubleSlash,DoTDA,PosP_Xk,MomQ_Xk)
       INTEGER                      :: N,II,JJ,KK
       LOGICAL                      :: DoTDA
       REAL(DOUBLE), DIMENSION(N,N) :: P,Q,F,Z,PosP_Xk,MomQ_Xk,L_PosP_Xk,L_MomQ_Xk, &
            PosP_Pk,MomQ_Pk,L_PosP_Pk,L_MomQ_Pk, &
            Xk,PosP,MomQ,LXk
       REAL(DOUBLE), DIMENSION(N,N,N,N) :: DoubleSlash
       REAL(DOUBLE)                 :: Lambda_PosP,Lambda_MomQ,E_EOM,Tmp_EOM,Ap,Aq, &
            Bp,Bq,Cp,Cq,Rpq,Spq,Tpq,Upq,LamP,LamQ,EOMThouless

       L_PosP_Xk=LiouvAO(N,F,P,Z,DoubleSlash,PosP_Xk)
       L_MomQ_Xk=LiouvAO(N,F,P,Z,DoubleSlash,MomQ_Xk)
       Ap=Trace_EOM(N,PosP_Xk,L_PosP_Xk)
       Aq=Trace_EOM(N,MomQ_Xk,L_MomQ_Xk)
       Rpq=Trace_EOM(N,PosP_Xk,MomQ_Xk)
       EOMThouless=Half*(Ap+Aq)/Rpq

     END FUNCTION EOMThouless




     SUBROUTINE RQI(N,M,Nam,O,S,MPI,G,B)
       INTEGER            :: N,M,I,J,K,L,U,V,JTDA,cBAS
       TYPE(FileNames)    :: Nam,RQINams
       TYPE(State)        :: S,RQIStat
       TYPE(Parallel)     :: MPI
       TYPE(Options)      :: O
       TYPE(Geometries)   :: G
       TYPE(BSET)         :: B
       LOGICAL            :: DoTDA
       INTEGER,DIMENSION(3) :: Cur
       !
       TYPE(BCSR)                      :: sP,sQ,sT,sZ
       TYPE(DBL_RNK2)                  :: P,Q,F,Z
       REAL(DOUBLE)                    :: Ek,EkOld,dEk,Beta,Lambda,ErrRel,ErrAbs,Shift,dNorm
       INTEGER                         :: XkNon0s,PkNon0s
       REAL(DOUBLE)                    :: XkThreshold,PkThreshold,PMax,Err

!!$
!!$    REAL(DOUBLE),DIMENSION(N,N)     :: Xk,Gk,Pk,LXk,LPk
!!$    REAL(DOUBLE),DIMENSION(N,N)     :: XkOld,PkOld,GkOld
!!$    REAL(DOUBLE),DIMENSION(N)       :: Values
!!$    REAL(DOUBLE),DIMENSION(N,N,M)   :: Vectors
!!$    REAL(DOUBLE),DIMENSION(N,N,N,N) :: TwoE,DoubleSlash
       !
       RQIStat=S
       RQINams=Nam
       cBAS=RQIStat%Current%I(2)
       !
       CALL Elapsed_TIME(TimeTotal,Init_O='Init')
       CALL Elapsed_TIME(TimeONX,Init_O='Init')
       CALL Elapsed_TIME(TimeQCTC,Init_O='Init')
       CALL Elapsed_TIME(TimeBCSR,Init_O='Init')
       CALL Elapsed_TIME(TimeTotal,Init_O='Start')

!!$    !
!!$    CALL Get(sP,TrixFile("OrthoD",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0))
!!$    CALL Get(sF,TrixFile("OrthoF",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0))
!!$    CALL Get(sZ,TrixFile("X",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I))
!!$    !
!!$    CALL SetEq(P,sP)
!!$    CALL SetEq(Q,sP)
!!$    P%D=Two*P%D
!!$    Q%D=-Two*Q%D
!!$    DO I=1,N
!!$       Q%D(I,I)=Q%D(I,I)+Two
!!$    ENDDO
!!$    !
!!$    CALL SetEq(F,sF)
!!$    CALL SetEq(Z,sZ)
!!$    !
!!$    CALL Integrals2E(B,G%Clone(1),TwoE)
!!$    DO I=1,N
!!$       DO J=1,N
!!$          DO K=1,N
!!$             DO L=1,N
!!$                DoubleSlash(I,J,K,L)=TwoE(I,J,K,L)-TwoE(I,K,J,L)/2D0
!!$             ENDDO
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
!!$
!!$    goto 111

       DO I=1,1
          !       write(*,*)' ek = ',EK
          !       GOTO 111 ! STOP
          DO JTDA=1,1
             !
             IF(JTDA==0)THEN
                DoTDA=.TRUE.
             ELSE
                DoTDA=.FALSE.
             ENDIF

             XkThreshold=1D-10
             PkThreshold=1D-10
             !
             CALL Nihilate0(N,I,Nam,S,MPI)
             CALL LOn2BakEnd(N,I,0,Shift,'Xk',Nam,S,MPI,XkThreshold,XkNon0s)
             CALL Nihilate1(I,Nam,S,MPI,Ek,TDA_O=DoTDA)
             !
             EkOld=BIG_DBL
             DO K=0,500
                ! The non-linear congjugate gradient
                CALL NLCGBakEnd(I,K,Ek,Nam,S,MPI,Beta)
                ! Compute L[Pk]
                CALL LOn2BakEnd(N,I,K,Shift,'Pk',Nam,S,MPI,PkThreshold,PkNon0s,PMax)
                ! Line Search: Min_Lambda{ E[Xk+Lambda*Pk] }
                CALL RQLSBakEnd(I,Nam,S,Lambda)
                ! Anhiliate and renorm Xk
                CALL NihilateXk(I,Nam,S,MPI,Lambda,dNorm,TDA_O=DoTDA)
                ! Compute L[Xk]
                CALL LOn2BakEnd(N,I,K,Shift,'Xk',Nam,S,MPI,XkThreshold,XkNon0s)
                ! Anihilate L[Xk], compute Ek and its relative error
                CALL NihilateLXk(I,Nam,S,MPI,Ek,dEk,TDA_O=DoTDA)
                !
                CALL OpenASCII(Nam%OFile,Out)
                WRITE(*  ,33)I,K,Ek*27.21139613182D0,dEk,ABS(dNorm),TimeONX%Wall,TimeQCTC%Wall,TimeBCSR%Wall, &
                     PMax,100D0*DBLE(XkNon0s)/DBLE(N*N),100D0*DBLE(PkNon0s)/DBLE(N*N)
                WRITE(Out,33)I,K,Ek*27.21139613182D0,dEk,ABS(dNorm),TimeONX%Wall,TimeQCTC%Wall,TimeBCSR%Wall, &
                     PMax,100D0*DBLE(XkNon0s)/DBLE(N*N),100D0*DBLE(PkNon0s)/DBLE(N*N)
                CLOSE(Out,STATUS='KEEP')
                !


                CALL OpenASCII('RQI_TauScaling2',77)
                WRITE(77,*)K,LOG10((MAX(1D-10,ABS(E_RPA-Ek*27.21139613182D0)/E_RPA))),XkNon0s,PkNon0s
                CLOSE(Unit=77)

                !             IF(JTDA==0.AND.K>5)EXIT

                IF(K>3.AND.dNorm<1D-2)THEN
                   ! Look for bad behavior
!!$                IF( Ek > EkOld .AND. ABS((Ek-EkOld)/Ek) > O%Thresholds(cBAS)%ETol )THEN
!!$
!!$                   ! Sign of variational principle broken, ostensibly due to N-scaling
!!$                   ! approximaitons.  If this happens, we are DONE!
!!$                   WRITE(*,*)' Converged due to variational violation ',Ek,EkOld,  &
!!$                        (Ek-EkOld)/Ek, O%Thresholds(cBAS)%ETol*1D2
!!$                   Ek=EkOld
!!$                   EXIT
!!$                ENDIF
                   ! Look for convergence (may be too tight)
                   IF(ABS(dEk)  <O%Thresholds(cBAS)%ETol*1D2.AND. &
                        ABS(dNorm)<O%Thresholds(cBAS)%DTol )THEN
                      WRITE(*,*)' Met convergence criteria ',dEk,dNorm,O%Thresholds(cBAS)%ETol, &
                           O%Thresholds(cBAS)%DTol
                      EXIT
                   ENDIF
                ENDIF
                EkOld=MIN(EkOld,Ek)

                !
             ENDDO
             CALL Elapsed_TIME(TimeTOTAL,Init_O='Accum')
             IF(DoTDA)THEN
                WRITE(*,44)'TDA:',I,K,Ek*27.21139613182D0,dEk,TimeTotal%Wall
             ELSE
                WRITE(*,44)'RPA:',I,K,Ek*27.21139613182D0,dEk,TimeTotal%Wall
             ENDIF
          ENDDO
       ENDDO
       !




33     FORMAT('St=',I2,', It=',I3,', Ev=',F10.6,', dE=',D8.2,', dN=',D7.2, &
            ', Tk=',D10.4,', Tj=',D10.4,', Tm=',D10.4,', |Gk|=',D8.2,', %Xk=',F6.2,', %Pk=',F6.2)


44     FORMAT(A4,' State=',I2,', Nk=',I3,', Ev=',F9.6,', dE=',D7.2,', WallSec=',D12.4)

!!$
!!$
!!$111 CONTINUE
!!$    DO I=1,1
!!$       !
!!$       DO JTDA=0,1
!!$
!!$          IF(JTDA==0)THEN
!!$             DoTDA=.TRUE.
!!$             CALL RPAGuess(N,Xk)
!!$          ELSE
!!$             DoTDA=.FALSE.
!!$!             CALL RPAGuess(N,Xk)
!!$
!!$          ENDIF
!!$
!!$          CALL Anihilate(N,P%D,Q%D,Xk,TDA_O=DoTDA)
!!$          CALL Renorm(N,P%D,Xk)
!!$          CALL LOn2(N,I,Shift,F%D,P%D,Z%D,DoubleSlash,Values,Vectors,Xk,LXk)
!!$          CALL Anihilate(N,P%D,Q%D,LXk,TDA_O=DoTDA)
!!$
!!$          Beta=Zero
!!$          XkOld=Zero
!!$          PkOld=Zero
!!$          Ek=ThoulessQ(N,P%D,Xk,LXk)
!!$          DO K=0,200
!!$             !
!!$             Gk=Two*(LXk-Ek*Xk)
!!$             IF(K>0)Beta=Pdot1(N,P%D,Gk,Gk-Gkold)/Pdot1(N,P%D,GkOld,GkOld)
!!$
!!$             Pk=Gk+Beta*PkOld
!!$             !
!!$             CALL LOn2(N,I,Shift,F%D,P%D,Z%D,DoubleSlash,Values,Vectors,Pk,LPk)
!!$             CALL RQILineSearch(N,P%D,Pk,Xk,LXk,LPk,Lambda)
!!$             !
!!$             EkOld=Ek
!!$             XkOld=Xk
!!$             EkOld=Ek
!!$             GkOld=Gk
!!$             PkOld=Pk
!!$             !
!!$             Xk=XkOld+Lambda*Pk
!!$             !
!!$             CALL Anihilate(N,P%D,Q%D,Xk,TDA_O=DoTDA)
!!$             dNorm=One-sqrt(abs(Pdot1(N,P%D,Xk,Xk)))
!!$             CALL ReNorm(N,P%D,Xk)
!!$             CALL LOn2(N,I,Shift,F%D,P%D,Z%D,DoubleSlash,Values,Vectors,Xk,LXk)
!!$             CALL Anihilate(N,P%D,Q%D,LXk,TDA_O=DoTDA)
!!$             Ek=ThoulessQ(N,P%D,Xk,LXk)
!!$             ErrAbs=Ek-EkOld
!!$             ErrRel=-1D10
!!$             DO U=1,N
!!$                DO V=1,N
!!$                   ErrRel=MAX(ErrRel,ABS(Ek*Xk(U,V)-LXk(U,V))/Ek)
!!$                ENDDO
!!$             ENDDO
!!$
!!$             IF(JTDA==0)THEN
!!$                Err=LOG10(Ek-0.432987734557198D0)
!!$!0.432987752380454D0)
!!$                WRITE(*,*)K,Err,Ek
!!$                WRITE(35,*)K,Err
!!$             ELSE
!!$                Err=LOG10(Ek- 0.406088082593828D0)!
!!$!0.406088110418713D0)
!!$!                WRITE(*,*)K,Err,Ek,Ek*27.21139613182D0
!!$                WRITE(*,*)K,Err,Ek
!!$                WRITE(34,*)K,Err
!!$             ENDIF
!!$
!!$
!!$
!!$!             WRITE(*,33)I,K,Ek*27.21139613182D0,Beta,Lambda,ErrRel,dNorm
!!$
!!$!             IF(ErrRel<1D-4)EXIT
!!$             IF(JTDA==0.AND.K==10)EXIT
!!$             !
!!$          ENDDO
!!$          WRITE(*,*)I,K,Ek*27.21139613182D0,ErrRel,ErrAbs
!!$          Values(I)=Ek
!!$          Vectors(:,:,I)=Xk
!!$
!!$       ENDDO
!!$    ENDDO

     END SUBROUTINE RQI

     FUNCTION TRACE_EOM(N,P,Q)
       INTEGER :: N,II
       REAL(DOUBLE), DIMENSION(N,N) :: P,Q,Tmp
       REAL(DOUBLE) :: TRACE_EOM
       Tmp=MATMUL(TRANSPOSE(P),Q)
       TRACE_EOM=0D0

       !    CALL PPrint(TMP,'TMP',Unit_O=6)

       DO II=1,N
          TRACE_EOM=TRACE_EOM+Tmp(II,II)
       ENDDO
       !    WRITE(*,*)' TRACE_EOM = ',TRACE_EOM
     END FUNCTION TRACE_EOM

     SUBROUTINE Nihilate0(N,I,Nam,S,MPI)
       !
       INTEGER               :: I,N,  ii,jj
       REAL(DOUBLE)          :: Norm
       TYPE(FileNames)       :: Nam
       TYPE(State)           :: S
       TYPE(Parallel)        :: MPI
       TYPE(BCSR)            :: sP,sQ,sXk,sI,sT,sT1 ! Nihilate0 delete list
       INTEGER, DIMENSION(3) :: Cur
       CHARACTER(LEN=DCL)    :: XkName,PName,QName
       !
       TYPE(DBL_RNK2)        :: X

       !-----------------------------------------------------------------------------
       CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
       ! Status same as ground state, but now SCF cycle# is RQI state#
       Cur=S%Current%I
       Cur(1)=I
       ! Naming of things
       XkName=TrixFile('OrthoXk',PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=Cur,OffSet_O=0)
       PName= TrixFile("OrthoD", PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
       QName= TrixFile("OrthoQ", PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
       ! Get occupied projector
       CALL Get(sP,PName)
       ! Generate virtual space projector if first time through
       IF(I==1)THEN
          CALL SetEq(sQ,sP)
          CALL Multiply(sQ,-One)
          CALL Add(sQ,One)
          CALL Put(sQ,QName)
       ELSE
          CALL Get(sQ,QName)
       ENDIF
       ! Start guess with particle-hole space ONLY. Starting in this space
       ! GUARANTEES that we will only ever have positive line search solutions.
       ! This is equivalent to the TDA.
       !
       CALL New(X,(/N,N/))
       X%D=Zero
       do ii=1,N
          do jj=MAX(1,ii-5),MIN(N,II+5)
             X%D(ii,jj)=RANDOM_DBL((/-One,One/))
          enddo
       enddo
       CALL SetEq(sXk,X)
       CALL Delete(X)
!!$    !!
       CALL SetEq(sXk,sP)
       DO II=1,sXk%NNon0
          sXk%MTrix%D(II)=sXk%MTrix%D(II)*Zero + sXk%MTrix%D(II) !*RANDOM_DBL((/-One,One/))
       ENDDO

       CALL Multiply(sQ,sXk,sT1)
       CALL Multiply(sT1,sP,sXk)
       CALL Delete(sQ)
       CALL Delete(sT1)

       ! Normalize the guess transition density.  Note factor
       ! of two has to do with normalization of P and Q.  Here,
       ! these projectors have eigenvalues == 1, not 2.
       CALL XPose(sXk,sT)
       Norm=SQRT(Two*ABS(OneDot(sP,sXk,sT)))
       sXk%MTrix%D=sXk%MTrix%D*(One/Norm)

       !    CALL PPrint(sXk,"NORMED",Unit_O=6)

       ! Put guess Xk to disk
!!$
!!$    CALL Get(sXk,KoopmansName)

       CALL Put(sXk,XkName)



       ! Tidy up
       CALL Delete(sT)
       CALL Delete(sP)
       CALL Delete(sXk)
       ! All done
       CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
     END SUBROUTINE Nihilate0


     FUNCTION TRACE_EOM2(N,P,Q)
       INTEGER :: N,II
       REAL(DOUBLE), DIMENSION(N,N) :: P,Q,Tmp
       REAL(DOUBLE) :: TRACE_EOM2
       CALL PPrint(P,'P',Unit_O=6)
       CALL PPrint(Q,'Q',Unit_O=6)


       Tmp=MATMUL(TRANSPOSE(P),Q)
       TRACE_EOM2=0D0

       CALL PPrint(TMP,'TMP',Unit_O=6)

       DO II=1,N
          TRACE_EOM2=TRACE_EOM2+Tmp(II,II)
       ENDDO
       WRITE(*,*)' TRACE_EOM = ',TRACE_EOM2
     END FUNCTION TRACE_EOM2

     FUNCTION pOp(N,P,Q,X)
       INTEGER :: N
       REAL(DOUBLE), DIMENSION(N,N) :: P,Q,X,pOp
       pOp=MATMUL(P,MATMUL(X,Q))+TRANSPOSE(MATMUL(Q,MATMUL(X,P)))
       pOp=25D-2*pOp
     END FUNCTION POp

     FUNCTION qOp(N,P,Q,X)
       INTEGER :: N
       REAL(DOUBLE), DIMENSION(N,N) :: P,Q,X,qOp
       qOp=MATMUL(P,MATMUL(X,Q))-TRANSPOSE(MATMUL(Q,MATMUL(X,P)))
       qOp=25D-2*qOp
     END FUNCTION QOp

     SUBROUTINE F2Rot(N,C,Two2)
       IMPLICIT NONE
       INTEGER :: N,I,J,K,L,II,JJ,KK,LL,alpha
       REAL(DOUBLE), DIMENSION(N,N) :: C
       REAL(DOUBLE), DIMENSION(N,N) :: Two2,Tmp2
       alpha=0
       DO I=1,N
          DO J=1,N
             Tmp2(I,J)=0D0
             DO II=1,N
                DO JJ=1,N
                   Tmp2(I,J)=Tmp2(I,J)+Two2(II,JJ)*C(II,I)*C(JJ,J)
                END DO
             END DO
          END DO
       END DO
       Two2=Tmp2
     END SUBROUTINE F2ROT

     SUBROUTINE F4Rot(N,C,Four4)
       IMPLICIT NONE
       INTEGER :: N,I,J,K,L,II,JJ,KK,LL,alpha
       REAL(DOUBLE), DIMENSION(N,N) :: C
       REAL(DOUBLE), DIMENSION(N,N,N,N) :: Four4,Tmp4
       alpha=0
       DO I=1,N
          DO J=1,N
             DO K=1,N
                DO L=1,N
                   Tmp4(I,J,K,L)=0D0
                   DO II=1,N
                      DO JJ=1,N
                         DO KK=1,N
                            DO LL=1,N
                               Tmp4(I,J,K,L)=Tmp4(I,J,K,L)+Four4(II,JJ,KK,LL)*C(II,I)*C(JJ,J)*C(KK,K)*C(LL,L)
                            END DO
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO
       Four4=Tmp4
     END SUBROUTINE F4ROT


     SUBROUTINE MDiag_DSYEVD(A,N,EigVal,OffSet)
       TYPE(DBL_RNK2)   :: A
       INTEGER          :: N
       TYPE(DBL_VECT)   :: EigVal
       TYPE(DBL_VECT)   :: Work
       TYPE(INT_VECT)   :: IWork
       INTEGER          :: K,LWORK,Info,LIWORK,LgN,OffSet
       DO K=4,10000
          IF(2**K>=N)THEN
             LgN=K
             EXIT
          ENDIF
       ENDDO
       LWORK=2*(1+5*N+2*N*LgN+3*N**2)
       LIWORK=2*(2+5*N)
       CALL New(Work,LWork)
       CALL New(IWork,LIWork)
       CALL DSYEVD('V','U',N,A%D(1,1+OffSet),N,EigVal%D(1+OffSet),Work%D(1),LWORK, &
            IWork%I(1),LIWORK,Info)
       IF(Info/=SUCCEED)CALL Halt('DSYEVD flaked in RHEqs. INFO='//TRIM(IntToChar(Info)))
       CALL Delete(IWork)
       CALL Delete(Work)
     END SUBROUTINE MDiag_DSYEVD








     SUBROUTINE Nihilate1(I,N,S,MPI,Ek,TDA_O)
       LOGICAL, OPTIONAL     :: TDA_O
       LOGICAL               :: TDA
       INTEGER               :: I
       REAL(DOUBLE)          :: Ek
       TYPE(FileNames)       :: N
       TYPE(State)           :: S
       TYPE(Parallel)        :: MPI
       TYPE(BCSR)            :: sP,sQ,sXk,sLXk,sT,sT1,sT2,sT3 ! Nihilate1 delete list
       INTEGER, DIMENSION(3) :: Cur
       CHARACTER(LEN=DCL)    :: XkName,LXkName,PName,QName
       !-----------------------------------------------------------------------------
       CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
       IF(PRESENT(TDA_O))THEN
          TDA=TDA_O
       ELSE
          TDA=.FALSE.
       ENDIF
       ! Status same as ground state, but now SCF cycle# is RQI state#
       Cur=S%Current%I
       Cur(1)=I
       ! Naming of things
       XkName=TrixFile('OrthoXk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
       LXkName=TrixFile('LXk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
       PName= TrixFile("OrthoD", PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
       QName= TrixFile("OrthoQ", PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
       ! Get occupied projector
       CALL Get(sP,PName)
       CALL Get(sQ,QName)
       ! Kluge due to poor memory management:
       CALL New(ST1)
       CALL New(SLXK)
       ! Get that should have good management (but doesnt)
       CALL Get(sLXk,LXkName)
       ! Anihilate via TDA or TD-SCF symmetry.
       ! Note that factor of (1/4) is not present, due to the
       ! normalization of P to eigenvalues with 1s or 0s.
       IF(TDA)THEN
          ! LXk=Q.LXk.P (only p-h zapped)
          CALL Multiply(sQ,sLXk,sT1)
          CALL Multiply(sT1,sP,sLXk)
          CALL Delete(sT1)
       ELSE
          ! LXk=(P.LXk.Q+Q.LXk.P) (both h-p and p-h zapped)
          CALL Multiply(sP,sLXk,sT1)
          CALL Multiply(sT1,sQ,sT2)
          CALL Multiply(sQ,sLXk,sT1)
          CALL Multiply(sT1,sP,sT3)
          CALL Add(sT2,sT3,sLXk)
          CALL Delete(sT1)
          CALL Delete(sT2)
          CALL Delete(sT3)
       ENDIF
       !
       CALL Delete(sQ)
       !
       CALL Put(sLXk,LXkName)
       !
       CALL Get(sXk,XkName)
       CALL XPose(sLXk,sT)
       !
       ! Again a factor of two for normalization of P
       Ek=Two*OneDot(sP,sXk,sT)
       !
       CALL Delete(sP)
       CALL Delete(sT)
       CALL Delete(sXk)
       CALL Delete(sLXk)
       ! All done
       CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
     END SUBROUTINE Nihilate1

     SUBROUTINE NihilateXk(I,N,S,MPI,Lambda,dNorm,TDA_O)
       LOGICAL, OPTIONAL     :: TDA_O
       LOGICAL               :: TDA
       INTEGER               :: I
       REAL(DOUBLE)          :: Lambda,Norm,dNorm
       TYPE(FileNames)       :: N
       TYPE(State)           :: S
       TYPE(Parallel)        :: MPI
       TYPE(BCSR)            :: sP,sQ,sXk,sPk,sLXk,sT,sT1,sT2,sT3 ! NihilateXk delete list
       INTEGER, DIMENSION(3) :: Cur
       CHARACTER(LEN=DCL)    :: XkName,PkName,LXkName,PName,QName
       !-----------------------------------------------------------------------------
       CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
       IF(PRESENT(TDA_O))THEN
          TDA=TDA_O
       ELSE
          TDA=.FALSE.
       ENDIF
       ! Status same as ground state, but now SCF cycle# is RQI state#
       Cur=S%Current%I
       Cur(1)=I

       ! Naming of things

       XkName=TrixFile('OrthoXk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
       PkName=TrixFile('OrthoPk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
       PName= TrixFile("OrthoD", PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
       QName= TrixFile("OrthoQ", PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)

       CALL New(sT1)
       CALL New(sPk)
       CALL New(sXk)

       CALL Get(sT1,XkName)
       CALL Get(sPk,PkName)
       sPk%MTrix%D=sPk%MTrix%D*Lambda
       CALL Add(sT1,sPk,sXk)
       !
       CALL Delete(sPk)

       CALL Get(sP,PName)
       CALL Get(sQ,QName)

       ! Anihilate via TDA or TD-SCF symmetry.
       ! Note that factor of (1/4) is not present, due to the
       ! normalization of P to eigenvalues with 1s or 0s.
       IF(TDA)THEN
          ! Xk=Q.Xk.P (only p-h zapped)
          CALL Multiply(sQ,sXk,sT1)
          CALL Multiply(sT1,sP,sXk)
          CALL Delete(sT1)
       ELSE
          ! Xk=(P.Xk.Q+Q.Xk.P) (both h-p and p-h zapped)
          CALL Multiply(sP,sXk,sT1)
          CALL Multiply(sT1,sQ,sT2)
          CALL Multiply(sQ,sXk,sT1)
          CALL Multiply(sT1,sP,sT3)
          CALL Add(sT2,sT3,sXk)
          CALL Delete(sT1)
          CALL Delete(sT2)
          CALL Delete(sT3)
       ENDIF
       !
       CALL Delete(sQ)
       ! Normalize the guess transition density.  Note factor
       ! of two has to do with normalization of P and Q.  Here,
       ! these projectors have eigenvalues == 1, not 2.
       CALL XPose(sXk,sT)
       Norm=SQRT(Two*ABS(OneDot(sP,sXk,sT)))
       dNorm=One-Norm
       sXk%MTrix%D=sXk%MTrix%D*(One/Norm)
       ! Put guess Xk to disk
       CALL Put(sXk,XkName)
       !
       CALL Delete(sT)
       CALL Delete(sP)
       CALL Delete(sXk)
       ! All done
       CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
     END SUBROUTINE NihilateXk
     !
     SUBROUTINE NihilateLXk(I,N,S,MPI,Ek,dEk,TDA_O)
       LOGICAL, OPTIONAL     :: TDA_O
       LOGICAL               :: TDA
       INTEGER               :: I,L
       REAL(DOUBLE)          :: dEk,Ek
       TYPE(FileNames)       :: N
       TYPE(State)           :: S
       TYPE(Parallel)        :: MPI
       TYPE(BCSR)            :: sP,sQ,sXk,sPk,sLXk,sT,sT1,sT2,sT3 ! NihilateLXk delete list
       INTEGER, DIMENSION(3) :: Cur
       CHARACTER(LEN=DCL)    :: XkName,PkName,LXkName,PName,QName
       !-----------------------------------------------------------------------------
       CALL Elapsed_TIME(TimeBCSR,Init_O='Start')

       IF(PRESENT(TDA_O))THEN
          TDA=TDA_O
       ELSE
          TDA=.FALSE.
       ENDIF
       ! Status same as ground state, but now SCF cycle# is RQI state#
       Cur=S%Current%I
       Cur(1)=I

       ! Naming of things

       LXkName=TrixFile('LXk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
       XkName=TrixFile('OrthoXk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
       PkName=TrixFile('OrthoPk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
       PName= TrixFile("OrthoD", PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
       QName= TrixFile("OrthoQ", PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
       !

       CALL New(sLXk)
       !
       CALL Get(sP,PName)
       CALL Get(sQ,QName)
       CALL Get(sLXk,LXkName)
       ! Anihilate via TDA or TD-SCF symmetry.
       ! Note that factor of (1/4) is not present, due to the
       ! normalization of P to eigenvalues with 1s or 0s.
       IF(TDA)THEN
          ! Xk=Q.Xk.P (only p-h zapped)
          CALL Multiply(sQ,sLXk,sT1)
          CALL Multiply(sT1,sP,sLXk)
          CALL Delete(sT1)
       ELSE
          ! Xk=(P.Xk.Q+Q.Xk.P) (both h-p and p-h zapped)
          CALL Multiply(sP,sLXk,sT1)
          CALL Multiply(sT1,sQ,sT2)
          CALL Multiply(sQ,sLXk,sT1)
          CALL Multiply(sT1,sP,sT3)
          CALL Add(sT2,sT3,sLXk)
          CALL Delete(sT1)
          CALL Delete(sT2)
          CALL Delete(sT3)
       ENDIF
       !
       CALL Delete(sQ)
       ! Put the anihilated LXk to disk
       CALL Put(sLXk,LXkName)
       !
       CALL Get(sXk,XkName)
       CALL XPose(sLXk,sT)
       ! Again a factor of two for normalization of P
       Ek=Two*OneDot(sP,sXk,sT)
       CALL Delete(sP)
       CALL Delete(sT)
       !
       sXk%MTrix%D=-Ek*sXk%MTrix%D
       CALL Add(sXk,sLXk,sT1)
       !
       dEk=-1D10
       DO L=1,sT1%NNon0
          dEk=MAX(dEk,ABS(sT1%MTrix%D(L)))
       ENDDO
       dEk=dEk/Ek
       !
       CALL Delete(sT1)
       CALL Delete(sXk)
       CALL Delete(sLXk)
       ! All done

       CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')

     END SUBROUTINE NihilateLXk

     SUBROUTINE NLCGBakEnd(I,K,Ek,N,S,MPI,Beta)
       !
       INTEGER               :: I,K
       REAL(DOUBLE)          :: Ek
       TYPE(FileNames)       :: N
       TYPE(State)           :: S
       TYPE(Parallel)        :: MPI
       TYPE(BCSR),SAVE       :: sXk,sLXk,sGk,sGkOld,sPkOld,sP,sT,sT1 ! NLCGradient delete list
       INTEGER, DIMENSION(3) :: Cur
       REAL(DOUBLE)          :: Num,Den,Beta
       CHARACTER(LEN=DCL)    :: XkName,LXkName,GkName,PName,GkOldName,PkOldName,PkName

       CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
       !
       Cur=S%Current%I
       Cur(1)=I
       !
       XkName=   TrixFile('OrthoXk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
       LXkName=  TrixFile('LXk',    PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
       GkName=   TrixFile('Gk',     PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
       GkOldName=TrixFile('GkOld',  PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
       PkName=   TrixFile('OrthoPk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)

       PkOldName=PkName

       PName=    TrixFile("OrthoD", PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
       !
       ! Get Xk and LXk
       CALL Get(sXk,XkName)
       CALL Get(sLXk,LXkName)

       ! Gk=Two*(LXk-Ek*Xk)

       sXk%MTrix%D=-Ek*sXk%MTrix%D

       CALL Add(sLXk,sXk,sGk)

       sGk%MTrix%D=Two*sGk%MTrix%D

       !    CALL PPrint(sGk,   ' Gk    ',Unit_O=6)

       ! Gk to disk
       CALL Put(sGk,GkName)
       ! Done with Xk and LXk
       CALL Delete(sXk)
       CALL Delete(sLXk)
       ! Beta=0 for K=0
       IF(K==0)THEN
          CALL Put(sGk,PkName)
          CALL Put(sGk,GkOldName)
          CALL Delete(sGk)
          RETURN
       ENDIF
       ! Get Pk and GkOld
       CALL Get(sP,PName)
       CALL Get(sGkOld,GkOldName)

       !    CALL PPrint(sGkOld,' GkOld ',Unit_O=6)

       !
       ! Beta=(Gk-Gkold,Gk)_p/(GkOld,GkOld)_p
       sGkOld%MTrix%D=-sGkOld%MTrix%D
       CALL Add(sGk,sGkOld,sT1)
       sGkOld%MTrix%D=-sGkOld%MTrix%D
       CALL XPose(sGk,sT)
       Num=OneDot(sP,sT1,sT)
       CALL XPose(sGkOld,sT)
       Den=OneDot(sP,sGkOld,sT)
       Beta=Num/Den
       !    WRITE(*,*)' Beta = ',Beta


       ! Done with P and T
       CALL Delete(sT)
       CALL Delete(sP)
       CALL Delete(sGkOld)
       ! Pk=Gk+Beta*PkOld
       CALL Get(sPkOld,PkOldName)

       sPkOld%MTrix%D=Beta*sPkOld%MTrix%D

       CALL Add(sGk,sPkOld,sT1)
       ! Put GkOld to disk
       CALL Put(sGk,GkOldName)
       ! Put Pk to disk
       CALL Put(sT1,PkName)
       ! Clean up and done
       CALL Delete(sT1)
       CALL Delete(sGk)
       CALL Delete(sPkOld)

       CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
       !
     END SUBROUTINE NLCGBAKEND

     SUBROUTINE LOn2BakEnd(N,I,K,Shift,Trgt,Nam,S,MPI,MatrixThreshold,MatrixNon0s,PMax_O)
       INTEGER                       :: N,I,J,K,MatrixNon0s
       TYPE(FileNames)               :: Nam
       TYPE(State)                   :: S,RQIStat
       TYPE(Parallel)                :: MPI
       REAL(DOUBLE),DIMENSION(N,N)   :: LX,X
       TYPE(BCSR)                    :: sF,sX,sJ,sK,sZ,sP,sJK,sT1,sT2,sT3
       REAL(DOUBLE)                  :: Shift,MatrixThreshold,LocalThreshold,PMax
       CHARACTER(LEN=*)              :: Trgt
       REAL(DOUBLE),OPTIONAL         :: PMax_O
       !----------------------------------------------------------------------------
       CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
       ! Set up local Invokation parameters based on ground state values
       CALL New(RQIStat%Action,2)
       CALL New(RQIStat%Current,3)
       CALL New(RQIStat%Previous,3)
       RQIStat%Current%I=S%Current%I
       RQIStat%Previous%I=S%Previous%I
       ! Action is TD-SCF with secondary parameter the product LX or LP (L[Xk], L[Pk])
       RQIStat%Action%C(1)="TD-SCF"
       RQIStat%Action%C(2)=TRIM(Trgt) !//TRIM(IntToChar(I))
       ! "SCF cycle" is the RQI State number
       RQIStat%Current%I(1)=I
       RQIStat%Previous%I(1)=I
       ! Get the orthogonal transition density matrix Xk (or CG gradient Pk) corresponding to the
       ! resultant of this subroutine, namely L[Xk] or L[Pk] in an orthongal representation
       CALL New(sX) ! Kluge.  Somehow, dimensioning not quite right here:
       CALL Get(sX,TrixFile('Ortho'//TRIM(Trgt),PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=RQIStat%Current%I,OffSet_O=0))
       ! Ground state fockian in an orthogonal representation
       CALL Get(sF,TrixFile("OrthoF",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0))
       ! T2=[Xor,For]
       CALL Multiply(sX,sF,sT2)
       CALL Multiply(sF,sX,sT2,-One)
       ! Done with F
       CALL Delete(sF)
       ! Z is the sparse inverse factor of S
       CALL Get(sZ,TrixFile("X",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I))
       ! Xao = Z^t.Xor.Z
       CALL Multiply(sZ,sX,sT1)
       CALL Multiply(sT1,sZ,sX)

       IF(PRESENT(PMax_O))THEN
          !       PMax_O=MAX(sX)
          PMax_O=FNorm(sX)
          LocalThreshold=1D1*MatrixThreshold*PMax_O
          LocalThreshold=MatrixThreshold
       ELSE
          LocalThreshold=MatrixThreshold
       ENDIF
       ! Filter small blocks and return the # of non zero elements
       CALL Filter(sX,Tol_O=LocalThreshold)
       MatrixNon0s=sX%NNon0
       ! This is the AO transition density matrix (or CG gradient)
       CALL Put(sX,TrixFile(Trgt,PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=RQIStat%Current%I,OffSet_O=0))
       ! Done with sX
       CALL Delete(sX)
       CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
       ! Build JK[X] in the AO basis
       CALL Elapsed_TIME(TimeONX,Init_O='Start')
       CALL Invoke('ONX',Nam,RQIStat,MPI)
       CALL Elapsed_TIME(TimeONX,Init_O='Accum')
       CALL Elapsed_TIME(TimeQCTC,Init_O='Start')
       CALL Invoke('QCTC',Nam,RQIStat,MPI)
       CALL Elapsed_TIME(TimeQCTC,Init_O='Accum')
       CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
       ! Pick up J and K
       CALL Get(sJ,TrixFile("J",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=RQIStat%Current%I,OffSet_O=0))
       CALL Get(sK,TrixFile("K",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=RQIStat%Current%I,OffSet_O=0))
       ! JK=Jao[X]+Kao[X]
       CALL Add(sJ,sK,sJK)
       ! Done with J and K
       CALL Delete(sJ)
       CALL Delete(sK)
       ! JK[X]=Zt.JKao[X].Z==JKor
       CALL Multiply(sZ,sJK,sT1)
       CALL Multiply(sT1,sZ,sJK)
       ! Done with Z
       CALL Delete(sZ)
       ! Get some P
       CALL Get(sP,TrixFile("OrthoD",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0))
       ! T1=[JKor,Por]
       CALL Multiply(sP,sJK,sT1)
       CALL Multiply(sJK,sP,sT1,-One)
       ! Done with P and JK
       CALL Delete(sP)
       CALL Delete(sJK)
       ! L[Xk]=[F,Xk]+[P,JK[X]] (orthogonal)
       CALL Add(sT1,sT2,sT3)
       ! Done with temporaries 1 and 2
       CALL Delete(sT1)
       CALL Delete(sT2)
       ! Put orthogonal L[Xk] or L[Pk] to disk (*.LX or *.LP)
       CALL Put(sT3,TrixFile("L"//TRIM(Trgt),PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=RQIStat%Current%I,OffSet_O=0))
       ! Done with temp #3
       CALL Delete(sT3)
       ! Done with invokation parameters
       CALL Delete(RQIStat%Action)
       CALL Delete(RQIStat%Current)
       CALL Delete(RQIStat%Previous)
       CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
       !----------------------------------------------------------------------------
       IF(I>1)STOP
     END SUBROUTINE LOn2BakEnd

     SUBROUTINE RQLSBakEnd(I,N,S,Lambda)
       INTEGER               :: I
       TYPE(FileNames)       :: N
       TYPE(State)           :: S
       TYPE(Parallel)        :: MPI
       TYPE(BCSR)            :: sP,sXk,sPk,sLXk,sLPk,sT
       INTEGER, DIMENSION(3) :: Cur
       REAL(DOUBLE)          :: Lambda,Lambda_p,Lambda_m
       REAL(DOUBLE)          :: XX,PP,XP,PX,PLP,XLP,XLX,PLX,AA,BB,CC
       !
       CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
       !

       Cur=S%Current%I
       Cur(1)=I
       !
       CALL Get(sP ,TrixFile("OrthoD",PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0))
       CALL Get(sXk,TrixFile('OrthoXk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0))
       CALL Get(sPk,TrixFile('OrthoPk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0))
       !
       CALL XPose(sXk,sT)
       XX =OneDot(sP,sXk,sT)
       PX =OneDot(sP,sPk,sT)
       !
       CALL XPose(sPk,sT)
       PP =OneDot(sP,sPk,sT)
       XP =OneDot(sP,sXk,sT)
       !
       CALL Get(sLXk,TrixFile('LXk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0))
       CALL XPose(sLXk,sT)
       CALL Delete(sLXk)
       XLX=OneDot(sP,sXk,sT)
       PLX=OneDot(sP,sPk,sT)
       !
       CALL Get(sLPk,TrixFile('LPk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0))
       CALL XPose(sLPk,sT)
       CALL Delete(sLPk)
       PLP=OneDot(sP,sPk,sT)
       XLP=OneDot(sP,sXk,sT)
       !
       CALL Delete(sP)
       CALL Delete(sT)
       CALL Delete(sXk)
       CALL Delete(sPk)
       !
       AA=PLP*(PX+XP)-PP*(PLX+XLP)
       BB=2.0*PLP*XX-2.0*XLX*PP
       CC=XX*(PLX+XLP)-XLX*(PX+XP)
       !
       lambda_p=(-BB+SQRT(BB*BB-4.0*AA*CC))/(2.0*AA)
       lambda_m=(-BB-SQRT(BB*BB-4.0*AA*CC))/(2.0*AA)
       !
       Lambda=Lambda_P
       !
       CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
       !
     END SUBROUTINE RQLSBakEnd
     !===============================================================================
     ! HERE ARE THE CONVENTIONAL DENSE MATRIX RQI ROUTINES
     !===============================================================================
     FUNCTION OneDotT(sP,sA,sB,TDA) RESULT(Tr)
       TYPE(BCSR) :: sP,sA,sB,sT1,sT2
       REAL(DOUBLE) :: Tr
       LOGICAL,OPTIONAL  :: TDA
       CALL XPose(sB,sT1)
       CALL Multiply(sP,sA,sT2)
       CALL Multiply(sA,sP,sT2,-One)
       Tr=Half*Trace(sT1,sT2)
       CALL Delete(sT1)
       CALL Delete(sT2)
     END FUNCTION OneDotT
     !---------------------------------------------------------------------
     ! Calculates TD-SCF scalar product Tr=Tr([A,P],B^+)
     ! Assumes B^t on input
     !---------------------------------------------------------------------
     FUNCTION OneDot(sP,sA,sBt) RESULT(Tr)
       TYPE(BCSR) :: sP,sA,sBt,sT1
       REAL(DOUBLE) :: Tr
       CALL Multiply(sP,sA,sT1)
       CALL Multiply(sA,sP,sT1,-One)
       Tr=Half*Trace(sT1,sBt)
       CALL Delete(sT1)
     END FUNCTION OneDot


     !
     SUBROUTINE Anihilate(N,P,Q,X,TDA_O)
       LOGICAL, OPTIONAL :: TDA_O
       LOGICAL :: TDA
       INTEGER :: N
       REAL(DOUBLE),DIMENSION(N,N) :: P,Q,X
       IF(PRESENT(TDA_O))THEN
          TDA=TDA_O
       ELSE
          TDA=.FALSE.
       ENDIF
       IF(TDA)THEN
          X=ProjectPH(N,P,Q,X)
       ELSE
          X=Project(N,P,Q,X)
       ENDIF
     END SUBROUTINE Anihilate

     SUBROUTINE LOn2(N,M,Shift,F,P,Z,TwoE,Values,Vectors,X,LX)
       INTEGER :: N,M,J
       REAL(DOUBLE),DIMENSION(N,N)     :: F,P,Z,X,LX,AA,BB,Temp,Com
       REAL(DOUBLE),DIMENSION(N,N,N,N) :: TwoE
       REAL(DOUBLE)                  :: Shift,OmegaPls,OmegaMns,WS
       REAL(DOUBLE),DIMENSION(:)     :: Values
       REAL(DOUBLE),DIMENSION(:,:,:) :: Vectors
       LX=LiouvAO(N,F  ,P  ,Z,TwoE,X )
       IF(M==1)RETURN
       WS=Values(M-1)-Values(1)+Shift
       Com=MATMUL(X,P)-MATMUL(P,X)
       DO J=1,M-1
          OmegaPls=Trace2(MATMUL(TRANSPOSE(Vectors(:,:,J)),Com),N)
          OmegaMns=Trace2(MATMUL(Vectors(:,:,J),Com),N)
          LX=LX+WS*(OmegaMns*TRANSPOSE(Vectors(:,:,J))+OmegaPls*Vectors(:,:,J))
       ENDDO
     END SUBROUTINE LOn2
     !

     Subroutine SparseRQILineSearch(P,Pk,Xk,LXk,LPk,Lambda,TDA)
       REAL(DOUBLE) :: Lambda,Lambda_p,Lambda_m
       TYPE(BCSR)   :: P,Pk,Xk,LXk,LPk,Tmp1
       REAL(DOUBLE) :: XX,PP,XP,PX,PLP,XLP,XLX,PLX,AA,BB,CC
       LOGICAL,OPTIONAL  :: TDA
       XX =OneDotT(P,Xk,Xk,TDA)
       PP =OneDotT(P,Pk,Pk,TDA)
       XP =OneDotT(P,Xk,Pk,TDA)
       PX =OneDotT(P,Pk,Xk,TDA)
       PLP=OneDotT(P,Pk,LPk,TDA)
       XLX=OneDotT(P,Xk,LXk,TDA)
       XLP=OneDotT(P,Xk,LPk,TDA)
       PLX=OneDotT(P,Pk,LXk,TDA)

!!$       WRITE(*,*)' XX = ',XX
!!$       WRITE(*,*)' PP = ',PP
!!$       WRITE(*,*)' XP = ',XP
!!$       WRITE(*,*)' PX = ',PX
!!$       WRITE(*,*)' PLP= ',PLP
!!$       WRITE(*,*)' XLX= ',XLX
!!$       WRITE(*,*)' XLP= ',XLP
!!$       WRITE(*,*)' PLX= ',PLX


       ! WRITE(*,*)' BEFORE EK = ',XLX/XX
       AA=PLP*(PX+XP)-PP*(PLX+XLP)
       BB=2.0*PLP*XX-2.0*XLX*PP
       CC=XX*(PLX+XLP)-XLX*(PX+XP)
       lambda_p=(-BB+SQRT(BB*BB-4.0*AA*CC))/(2.0*AA)
       lambda_m=(-BB-SQRT(BB*BB-4.0*AA*CC))/(2.0*AA)
       Lambda=Lambda_P
       !    WRITE(*,*)' AFTER EK P= ',(XLX+lambda_p*XLP+lambda_p*PLX+lambda_p**2 * PLP) /(XX+lambda_p*XP+lambda_p*PX+lambda_p**2 * PP)
       !    WRITE(*,*)' AFTER EK M= ',(XLX+lambda_m*XLP+lambda_m*PLX+lambda_m**2 * PLP) /(XX+lambda_m*XP+lambda_m*PX+lambda_m**2 * PP)
     END Subroutine SparseRQILineSearch


     FUNCTION ThoulessQ(N,P,X,LX) RESULT(Ek)
       !! So, where is  the denominator??
       INTEGER :: N
       REAL(DOUBLE) :: Ek
       REAL(DOUBLE),DIMENSION(N,N)     :: F,P,X,LX,Tmp1
       Ek=Pdot1(N,P,X,LX)
     END FUNCTION ThoulessQ

     FUNCTION LiouvDot(N,BB,DSao,temp2)  RESULT(temp1)
       ! Calculates action of the Coulomb operator in AO space temp1=BB * (ij||kl)
       IMPLICIT NONE
       INTEGER :: I,J,K,N,one
       REAL (DOUBLE),DIMENSION(N*N):: BB,temp2
       REAL (DOUBLE),DIMENSION(N,N):: temp1
       REAL(DOUBLE),DIMENSION(N*N,N*N):: DSao
       REAL(DOUBLE) :: ddot

       one=1
       K=0
       DO I=1,N
          DO J=1,N
             K=K+1
             temp2=DSao(:,K)
             ! temp1(J,I)= ddot(N*N,BB,one,temp2,one)     ! This line is
             temp1(J,I)=DOT_PRODUCT(BB,Temp2)           ! the most CPU consuming step
          ENDDO
       END DO

     END FUNCTION LiouvDot

     FUNCTION LiouvAO(N,For,Por,X,DSao,AA)  RESULT(BB)
       ! Calculates action of the Liouville operator in AO space BB=L AA, (ij||kl)
       IMPLICIT NONE
       INTEGER :: I,J,M,K,L,N,one
       REAL (DOUBLE),DIMENSION(N,N)::For,Por,AA,BB,temp1,temp2,X
       REAL(DOUBLE),DIMENSION(N,N,N,N):: DSao
       REAL(DOUBLE) :: E,ddot

       ! AA to AO
       one=1
       BB=MATMUL(TRANSPOSE(X),(MATMUL(AA,X)))
       !    CALL PPrint(BB,'INPUT 2',Unit_O=6)

       temp1=Zero
       IF(.NOT.DoKoopmans)THEN
          DO I=1,N
             DO J=1,N
                DO K=1,N
                   DO L=1,N
                      temp1(I,J)=temp1(I,J)+BB(K,L)*DSao(K,L,I,J)
                   END DO
                END DO
             END DO
          END DO
       ENDIF
       !   CALL PPrint(temp1,'AO_JK[X]',Unit_O=6)
       BB=MATMUL(For,AA)-MATMUL(AA,For)
       ! temp back to orthog
       temp2=MATMUL(TRANSPOSE(X),(MATMUL(temp1,X)))
       BB=BB+MATMUL(temp2,Por)-MATMUL(Por,temp2)
     END FUNCTION LiouvAO


     FUNCTION LiouvAO2E(N,For,Por,X,DSao,AA)  RESULT(BB)
       ! Calculates action of the Liouville operator in AO space BB=L AA, (ij||kl)
       IMPLICIT NONE
       INTEGER :: I,J,M,K,L,N,one
       REAL (DOUBLE),DIMENSION(N,N)::For,Por,AA,BB,temp1,temp2,X
       REAL(DOUBLE),DIMENSION(N,N,N,N):: DSao
       REAL(DOUBLE) :: E,ddot

       ! AA to AO
       one=1
       BB=MATMUL(TRANSPOSE(X),(MATMUL(AA,X)))
       !    CALL PPrint(BB,'INPUT 2',Unit_O=6)

       temp1=Zero
       IF(.NOT.DoKoopmans)THEN
          DO I=1,N
             DO J=1,N
                DO K=1,N
                   DO L=1,N
                      temp1(I,J)=temp1(I,J)+BB(K,L)*DSao(K,L,I,J)
                   END DO
                END DO
             END DO
          END DO
       ENDIF
       !   CALL PPrint(temp1,'AO_JK[X]',Unit_O=6)
       !    BB=MATMUL(For,AA)-MATMUL(AA,For)
       ! temp back to orthog
       temp2=MATMUL(TRANSPOSE(X),(MATMUL(temp1,X)))
       BB=MATMUL(temp2,Por)-MATMUL(Por,temp2)
     END FUNCTION LiouvAO2E

     SUBROUTINE RPAGuess(N,X)
       INTEGER :: N,I,J
       REAL(DOUBLE), DIMENSION(N,N) :: X
       X=One
!!$
       do i=1,N
          do j=1,N
             X(i,j)=RANDOM_DBL((/-One,One/))
          enddo
       enddo

       !    X=HALF*(X-TRANSPOSE(X))

       !    temp1 = ProjectPH(N,Qor,Por,Xk)
       !    Xk = temp1/sqrt(abs(Pdot1(N,Por,temp1,temp1,tmp1)))
     END SUBROUTINE RPAGuess

     SUBROUTINE RANDOMIZE(N,X,Tau)
       INTEGER :: N,I,J
       REAL(DOUBLE),DIMENSION(N,N) :: X
       REAL(DOUBLE) :: Tau,XABS
       do i=1,N
          do j=1,N
             XABS=X(I,J)
             X(i,j)=X(i,j)+RANDOM_DBL((/-Tau,Tau/))
          enddo
       enddo
     END SUBROUTINE RANDOMIZE


     SUBROUTINE ReNorm(N,P,X)
       INTEGER :: N
       REAL(DOUBLE) :: Norm
       REAL(DOUBLE),DIMENSION(N,N) :: P,X
       Norm=sqrt(abs(Pdot1(N,P,X,X)))
       X=X/Norm
     END SUBROUTINE ReNorm

     !************************************************************************
     FUNCTION Pdot(N,P,AA,BB,CC) RESULT(Tr)
       ! Calculates RPA scalar product Tr=Tr([AA^+,P],BB)

       IMPLICIT NONE
       INTEGER :: N
       REAL(DOUBLE) :: Tr
       REAL (DOUBLE),DIMENSION(N,N)::P,AA,BB,CC

       CC=MATMUL((MATMUL(TRANSPOSE(AA),P)-MATMUL(P,TRANSPOSE(AA))),BB)
       Tr=0.5*Trace2(CC,N)

     END FUNCTION Pdot

     !************************************************************************
     FUNCTION Pdot1(N,P,AA,BB) RESULT(Tr)
       ! Calculates RPA scalar product Tr=Tr([AA,P],BB^+)

       IMPLICIT NONE
       INTEGER :: N
       REAL(DOUBLE) :: Tr
       REAL (DOUBLE),DIMENSION(N,N)::P,AA,BB,CC

       CC=MATMUL((MATMUL(AA,P)-MATMUL(P,AA)),TRANSPOSE(BB))
       Tr=0.5*Trace2(CC,N)

     END FUNCTION Pdot1
     !-------------------------------------------------------------------------------
     FUNCTION Project(N,P,Q,AA)  RESULT(BB)
       ! BB=P AA Q + Q AA P
       ! calculates  projection to p-h an h-p space using Q and P

       IMPLICIT NONE
       INTEGER :: N
       REAL (DOUBLE),DIMENSION(N,N)::P,Q,AA,BB
       BB=0.25*(MATMUL(MATMUL(P,AA),Q)+MATMUL(MATMUL(Q,AA),P))
     END FUNCTION Project

     !-------------------------------------------------------------------------------
     FUNCTION ProjectPH(N,P,Q,AA)  RESULT(BB)
       ! BB=Q AA P (X-component, large)   (0  Y)
       ! BB=P AA Q (Y-component, small)   (X  0)
       ! calculates  projection to p-h OR h-p space using Q and P

       IMPLICIT NONE
       INTEGER :: N
       REAL (DOUBLE),DIMENSION(N,N)::P,Q,AA,BB
       BB=0.25D0*(MATMUL(MATMUL(Q,AA),P))
     END FUNCTION ProjectPH

     FUNCTION Trace2(Matrix,N) RESULT(Tr)
       IMPLICIT NONE
       INTEGER :: I,N
       REAL(DOUBLE) :: Tr
       REAL(DOUBLE),DIMENSION(N,N)::Matrix
       Tr=0D0
       DO I=1,N
          Tr=Tr+Matrix(I,I)
       END DO
     END FUNCTION Trace2

   END MODULE RayleighQuotientIteration
