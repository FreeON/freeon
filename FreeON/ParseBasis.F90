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
MODULE ParseBasis
  USE Order
  USE InOut
  USE PrettyPrint
  USE ControlStructures
  USE BasisSetParameters
  USE MondoLogger

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: BASIS_SETS     = 'BasisSets'

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: BCSR_MAXNON0S   = 'maxnon0s'
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: BCSR_MAXNBLKS   = 'maxnblks'
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: BCSR_MAXDEFAULT = 'default'

CONTAINS
  !============================================================================
  !
  !
  !============================================================================
  SUBROUTINE LoadBasisSets(N,O,G,B)
    TYPE(FileNames)    :: N
    TYPE(Options)      :: O
    TYPE(Geometries)   :: G
    TYPE(BasisSets)    :: B
    TYPE(INT_VECT)     :: ArrMaxNon0s,ArrMaxNBlks
    INTEGER            :: I,J,Ntot,BasU
    CHARACTER(LEN=DCL) :: BaseFile,GhostFile
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: Line
    LOGICAL            :: Exists
    !-------------------------------------------------------------------------!
    CALL OpenASCII(N%IFile,Inp)
    ! Find out how many basis sets we are going over, and their names
    CALL ParseBasisNames(B)
    !--------------------------------------
    ! Look for MaxNon0s and MaxNBlks.
    ! Allocate new array.
    CALL New(ArrMaxNon0s,B%NBSets);CALL New(ArrMaxNBlks,B%NBSets)
    ArrMaxNon0s%I=-BIG_INT;ArrMaxNBlks%I=-BIG_INT
    CALL ParseMaxElemBCSR(B%NBSets,Inp,ArrMaxNon0s,ArrMaxNBlks)
    !--------------------------------------
    ! Allocate aux stuff...
    ALLOCATE(B%BSiz(1:G%Clones,1:B%NBSets))
    ALLOCATE(B%OffS(1:G%Clones,1:B%NBSets))
    ALLOCATE(B%BSets(1:G%Clones,1:B%NBSets))
    ALLOCATE(B%LnDex(1:G%Clones,1:B%NBSets))
    ALLOCATE(B%DExpt(1:G%Clones,1:B%NBSets))
    ALLOCATE(B%AtomPairThresh(1:G%Clones,1:B%NBSets))
    ALLOCATE(B%PrimPairThresh(1:G%Clones,1:B%NBSets))

    DO J=1,B%NBSets
      BaseFile=TRIM(N%M_HOME)//'BasisSets/'//TRIM(B%BName(J))//BasF
      CALL MondoLog(DEBUG_MAXIMUM, "FreeON", "trying to open basis set file: "//TRIM(BaseFile), "LoadBasisSets")
      INQUIRE(FILE=BaseFile,EXIST=Exists)
      IF(Exists)THEN
        CALL MondoLog(DEBUG_MAXIMUM, "FreeON","loading external basis set  "//TRIM(BaseFile), "LoadBasisSets")
        CALL OpenASCII(BaseFile,Bas,OldFileQ_O=.TRUE.)
        BasU=Bas
      ELSE
        CALL MondoLog(DEBUG_MAXIMUM, "FreeON","looking for <BeginBasisSet"//TRIM(B%BName(J))//"> delimiter in inputfile ","LoadBasisSets")
        CALL LowCase(B%BName(J))
        BasU=Inp
      ENDIF
      DO I=1,G%Clones
        IF(Exists)THEN
          REWIND(Bas)
        ELSE
          REWIND(Inp)
          CALL AlignLowCase('<beginbasisset'//TRIM(B%BName(J))//'>',Inp)
        ENDIF
        ! Reset Allocation flag on BSets.
        CALL Initialize(B%BSets(I,J))
        IF(.NOT.ParseBasisSets(G%Clone(I),B%BSets(I,J),BasU,B%BSiz(I,J),B%OffS(I,J)))THEN
          CALL MondoHalt(PRSE_ERROR,'ParseBasisSets failed for basis set '//TRIM(B%BName(J)))
        ENDIF
        B%BSets(I,J)%BName=B%BName(J)
        CALL BCSRDimensions(G%Clone(I),B%BSets(I,J),O%AccuracyLevels(J), &
             B%MxAts(J),B%MxBlk(J),B%MxN0s(J),ArrMaxNon0s%I(J),ArrMaxNBlks%I(J))
        ! Compute primitive distribution statistics
        CALL PrimitiveReDistribution(B%BSets(I,J),B%NExpt(J),B%DExpt(I,J),B%Lndex(I,J))
        ! Compute basis set dependent distance thresholds
        CALL DistanceThresholdSetUp(O%Thresholds(J)%Dist,B%DExpt(I,J)%D(1), &
             B%AtomPairThresh(I,J),B%PrimPairThresh(I,J))
#ifdef FULL_ON_FRONT_END_DEBUG
        PrintFlags%Set=DEBUG_BASISSET
        CALL Print_BSET(B%BSets(I,J))
#endif
      ENDDO
      CLOSE(Bas,STATUS='KEEP')
    ENDDO
    CLOSE(Inp,STATUS='KEEP')
    !Delete Arrays.
    CALL Delete(ArrMaxNon0s)
    CALL Delete(ArrMaxNBlks)
  END SUBROUTINE LoadBasisSets
  !============================================================================
  ! THIS ROUTINE IS A CLONE OF THOSE IN Modules/Thresholding.f90, BUT DOES
  ! NOT USE GLOBAL VARIABLES, WHICH ARE VERBOTEN IN THE FRONT END
  !============================================================================
  SUBROUTINE DistanceThresholdSetUp(Tau,MinZab,AtomPairThresh,PrimPairThresh)
    REAL(DOUBLE) :: Tau,MinZab,MinXab,AtomPairThresh,PrimPairThresh
    !-------------------------------------------------------------------------!
    ! MinXab=MinZa*MinZb/(MinZa+MinZb)=MinZab/4
    MinXab=MinZab/Four
    ! Set Atom-Atom thresholds
    AtomPairThresh=-LOG(Tau)/MinXab
    ! Set Prim-Prim thresholds
    PrimPairThresh=-LOG(Tau)
  END SUBROUTINE DistanceThresholdSetUp
  !============================================================================
  ! WEAK ATTEMPT TO ESTIMATE NUMBER OF NONZEROS FOR SPARSE MATRICES, NEADS TWIDDLING
  ! AND SHOULD BECOME OBSOLETE IN NEAR FURTURE AS FASTMAT TAKES OVER...
  !============================================================================
  SUBROUTINE BCSRDimensions(G,B,AccL,MaxAtms,MaxBlks,MaxNon0,MaxNon0s,MaxNBlks)
    TYPE(CRDS) :: G
    TYPE(BSET) :: B
    INTEGER    :: AccL,MaxAtms,MaxBlks,MaxNon0
    INTEGER, INTENT(IN) :: MaxNon0s,MaxNBlks ! From input.
    REAL(DOUBLE) :: BWEstim,ShrinkFactor
    REAL(DOUBLE),PARAMETER,DIMENSION(4) :: BandWidth=(/ 1.D3, 1.3D3, 1.3D3,1.6D3/)
    REAL(DOUBLE),PARAMETER,DIMENSION(4) :: BWDecay  =(/ 1.D-4,1.D-4,1.D-3,1.D-2/)

    ! Use assymptotics to set the max matrix dimensions
    MaxAtms=1+G%NAtms
    BWEstim=MIN(G%NAtms,CEILING((DBLE(G%NAtms) &
         +BandWidth(AccL)*BWDecay(AccL)*DBLE(G%NAtms)**2) &
         /(One+BWDecay(AccL)*DBLE(G%NAtms)**2) ) )
    MaxBlks=1+G%NAtms*BWEstim
    MaxNon0=1+B%NBasF*(DBLE(B%NBasF)*DBLE(BWEstim)/DBLE(G%NAtms))

    ! Print out what we found.
    CALL MondoLog(DEBUG_MAXIMUM, "BCSRDimensions", "estimated MaxBlks = "//TRIM(IntToChar(MaxBlks)))
    CALL MondoLog(DEBUG_MAXIMUM, "BCSRDimensions", "estimated MaxNon0 = "//TRIM(IntToChar(MaxNon0)))

    ! Set the variable if def in the input.
    IF(MaxNon0s.GT.0) MaxNon0=MaxNon0s !if def in the input.
    IF(MaxNBlks.GT.0) MaxBlks=MaxNBlks !if def in the input.

    ! Print out result.
    CALL MondoLog(DEBUG_MAXIMUM, "BCSRDimensions", "MaxBlks = "//TRIM(IntToChar(MaxBlks)))
    CALL MondoLog(DEBUG_MAXIMUM, "BCSRDimensions", "MaxNon0 = "//TRIM(IntToChar(MaxNon0)))

  END SUBROUTINE BCSRDimensions
  !============================================================================
  ! PARSE A BASIS SET, SET UP ITS INDECIES, NORMALIZE THE PRIMITIVES AND
  ! COMPUTE BLOCKING FOR SPARSE BLOCKED LINEAR ALGEBRA
  !============================================================================
  FUNCTION ParseBasisSets(G,B,BasU,BlkSiz,OffSet)
    TYPE(BSET)                            :: BS,B
    TYPE(CRDS)                            :: G
    TYPE(INT_VECT)                        :: BlkSiz,OffSet
    TYPE(CHR_VECT)                        :: C
    REAL(DOUBLE), DIMENSION(1:MaxASymt+2) :: Dbls
    INTEGER, DIMENSION(2)                 :: Ints
    LOGICAL                               :: ParseBasisSets
    CHARACTER(LEN=DCL)                    :: Line
    INTEGER                               :: I,J,K,L,N,NC,NK,NP,NS,MinL,MaxL,KFound,Prim,Ell,BasU

    ParseBasisSets=.FALSE.

    ! Allocate temporary set
    BS%LMNLen=LHGTF(MaxAsymt)
    BS%NCtrt=MaxCntrx
    BS%NPrim=MaxPrmtv
    BS%NAtms=G%NAtms
    BS%NKind=G%NAtms
    BS%HasECPs=.TRUE.
    BS%MxProjL=5
    BS%Typ1Fnk=MaxPrmtv
    BS%Typ2Fnk=MaxPrmtv
    CALL New(BS)
    ! Count kinds and load basis set kind index
    BS%NKind=1
    BS%Kinds%I(1)=G%AtNum%D(1)
    BS%AtNam%C(1)=G%AtNam%C(1)
    DO I=2,BS%NAtms
      DO J=1,BS%NKind
        IF(BS%AtNam%C(J)==G%AtNam%C(I))GOTO 10
        !vwIF(BS%Kinds%I(J)==G%AtNum%D(I))GOTO 10
      ENDDO
      BS%NKind=BS%NKind+1
      BS%Kinds%I(BS%NKind)=G%AtNum%D(I)
      BS%AtNam%C(BS%NKind)=G%AtNam%C(I)
10    CONTINUE
    ENDDO
    ! Now load geometry atom type (kinds pointer) array.
    G%NKind=BS%NKind
    DO K=1,BS%NKind
      DO I=1,G%NAtms
        IF(BS%AtNam%C(K)==G%AtNam%C(I))G%AtTyp%I(I)=K
        !vwIF(BS%Kinds%I(K)==G%AtNum%D(I))G%AtTyp%I(I)=K
      ENDDO
    ENDDO
    ! Zero counters
    BS%NASym=0
    BS%NCtrt=0
    BS%NPrim=0
    DO I=1,BS%NKind
      BS%NCFnc%I(I)=0
      BS%NTyp1PF%I(I)=0
      DO J=1,MaxCntrx
        BS%NPFnc%I(J,I)=0
      ENDDO
    ENDDO
    BS%HasECPs=.False.
    BS%MxProjL=0
    BS%Typ1Fnk=0
    BS%Typ2Fnk=0
    BS%NCoreEl%D=Zero
    ! Parse basis set
    KFound=0
    DO
      READ(BasU,DEFAULT_CHR_FMT,END=99)Line
      DO NK=1,BS%NKind
        ! Look for the basis set
        IF(KeyQ(Line,BS%AtNam%C(NK)).AND.KeyQ(Line,'0'))THEN
          !vwIF(KeyQ(Line,Ats(BS%Kinds%I(NK))).AND.KeyQ(Line,'0'))THEN
          NC=0
          KFound=KFound+1
          DO
            READ(BasU,DEFAULT_CHR_FMT,END=99)Line
            IF(INDEX(Line,'<End')/=0)RETURN
            IF(KeyQ(TRIM(Line),Stars))GOTO 100
            NC=NC+1
            DO K=1,MaxLTyps
              IF(KeyQ(Line,CLTyps(K)))THEN
                BS%ASymm%I(1,NC,NK)=LTyps(1,K)
                BS%ASymm%I(2,NC,NK)=LTyps(2,K)
                DO L=1,MaxPrmtv
                  IF(KeyQ(Line,TRIM(IntToChar(L))))THEN
                    NP=L
                    GOTO 101
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
            RETURN
101         CONTINUE
            BS%NCFnc%I(NK)=NC
            BS%NPFnc%I(NC,NK)=NP
            BS%NCtrt=Max(BS%NCtrt,NC)
            BS%NPrim=Max(BS%NPrim,NP)
            MinL=BS%ASymm%I(1,NC,NK)
            MaxL=BS%ASymm%I(2,NC,NK)
            BS%NAsym=Max(BS%NAsym,MaxL)
            DO NP=1,BS%NPFnc%I(NC,NK)
              READ(BasU,DEFAULT_CHR_FMT,END=99)Line
              N=MaxL-MinL+2
              CALL LineToDbls(Line,N,Dbls)
              BS%Expnt%D(NP,NC,NK)=Dbls(1)
              K=1
              BS%CCoef%D(1:,NP,NC,NK)=Zero
              DO NS=MinL,MaxL
                K=K+1
                BS%CCoef%D(NS+1,NP,NC,NK)=Dbls(K)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        ! Look for an ECP in the basis set
        IF(INDEX(Line,'ECP')/=0)THEN
          ! Add a delimeter (space) to allow correct parsing
          Line=" "//Line
          IF(KeyQ(Line,TRIM(Ats(BS%Kinds%I(NK)))//'-ecp'))THEN
            BS%HasECPs=.TRUE.
            ! Parse in the ECP info
            CALL LineToChars(Line,C)
            BS%ProjEll%I(NK)=CharToInt(C%C(2))-1
            BS%NCoreEl%D(NK)=CharToInt(C%C(3))
            BS%MxProjL=MAX(BS%MxProjL,BS%ProjEll%I(NK))
            CALL Delete(C)
            ! Read a dummy line, something about potentials
            READ(BasU,DEFAULT_CHR_FMT,END=99)Line
            ! This line should contain number of primitives in this ECP
            READ(BasU,DEFAULT_CHR_FMT,END=99)Line
            CALL LineToInts(Line,1,Ints)
            BS%NTyp1PF%I(NK)=Ints(1)
            BS%Typ1Fnk=MAX(BS%Typ1Fnk,BS%NTyp1PF%I(NK))
            ! Parse in the type one potential
            DO Prim=1,BS%NTyp1PF%I(NK)
              READ(BasU,DEFAULT_CHR_FMT,END=99)Line
              CALL LineToInts(Line,1,Ints)
              CALL LineToDbls(Line,3,Dbls)
              ! Load this ECPs primitive angular symmetries, exponents and coefficients
              BS%Typ1Ell%I(Prim,NK)=Ints(1)
              BS%Typ1Exp%D(Prim,NK)=Dbls(2)
              BS%Typ1CCo%D(Prim,NK)=Dbls(3)
            ENDDO
            ! Parse in the type two potential
            DO Ell=0,BS%ProjEll%I(NK)
              ! Read a dummy line, something about potentials
              READ(BasU,DEFAULT_CHR_FMT,END=99)Line
              ! This line should contain number of primitives in this ECP
              READ(BasU,DEFAULT_CHR_FMT,END=99)Line
              ! Parse in the type two potential for this ell value
              CALL LineToInts(Line,1,Ints)
              BS%NTyp2PF%I(Ell,NK)=Ints(1)
              BS%Typ2Fnk=MAX(BS%Typ2Fnk,BS%NTyp2PF%I(Ell,NK))
              DO Prim=1,BS%NTyp2PF%I(Ell,NK)
                READ(BasU,DEFAULT_CHR_FMT,END=99)Line
                CALL LineToDbls(Line,3,Dbls)
                CALL LineToInts(Line,1,Ints)
                ! Load this ECPs primitive angular symmetries, exponents and coefficients
                BS%Typ2Ell%I(Prim,Ell,NK)=Ints(1)
                BS%Typ2Exp%D(Prim,Ell,NK)=Dbls(2)
                BS%Typ2CCo%D(Prim,Ell,NK)=Dbls(3)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDDO
100   CONTINUE
    ENDDO
99  CONTINUE
    !   Continue On
    IF(KFound/=BS%NKind) THEN
      CALL MondoLog(DEBUG_NONE, "ParseBasis", "Could not find a basis for all atoms requested.")
      CALL MondoLog(DEBUG_NONE, "ParseBasis", "The basis set file may not contain the desired atoms.")
      RETURN
    ENDIF
    ! Computing basis set indexing
    CALL BSetIndx(BS)
    ! Normalize the primitives
    CALL ReNormalizePrimitives(BS)
    ! Set the correct dimensions
    B%NAsym=BS%NAsym
    B%LMNLen=BS%LMNLen
    B%NCtrt=BS%NCtrt
    B%NPrim=BS%NPrim
    B%NAtms=BS%NAtms
    B%NKind=BS%NKind
    B%HasECPs=BS%HasECPs
    B%MxProjL=BS%MxProjL
    B%Typ1Fnk=BS%Typ1Fnk
    B%Typ2Fnk=BS%Typ2Fnk
    ! New basis set
    CALL New(B)
    ! Set old eq to new (should go to seteq sometime soon...)
    B%Kinds%I(1:B%NKind)=BS%Kinds%I(1:B%NKind)
    B%NCFnc%I(1:B%NKind)=BS%NCFnc%I(1:B%NKind)
    B%BFKnd%I(1:B%NAtms)=BS%BFKnd%I(1:B%NAtms)
    B%LxDex%I(1:B%LMNLen)=BS%LxDex%I(1:B%LMNLen)
    B%LyDex%I(1:B%LMNLen)=BS%LyDex%I(1:B%LMNLen)
    B%LzDex%I(1:B%LMNLen)=BS%LzDex%I(1:B%LMNLen)
    B%LStrt%I(1:B%NCtrt,1:B%NKind)=BS%LStrt%I(1:B%NCtrt,1:B%NKind)
    B%LStop%I(1:B%NCtrt,1:B%NKind)=BS%LStop%I(1:B%NCtrt,1:B%NKind)
    B%NPFnc%I(1:B%NCtrt,1:B%NKind)=BS%NPFnc%I(1:B%NCtrt,1:B%NKind)
    B%ASymm%I(:,1:B%NCtrt,1:B%NKind)=BS%ASymm%I(:,1:B%NCtrt,1:B%NKind)
    B%Expnt%D(1:B%NPrim,1:B%NCtrt,1:B%NKind)=BS%Expnt%D(1:B%NPrim,1:B%NCtrt,1:B%NKind)
    B%CCoef%D(1:B%LMNLen,1:B%NPrim,1:B%NCtrt,1:B%NKind)=BS%CCoef%D(1:B%LMNLen,1:B%NPrim,1:B%NCtrt,1:B%NKind)
    IF(B%HasECPs)THEN
      ! If we have pseudopotentials, copy them over too
      B%NCoreEl%D(1:B%NKind)=BS%NCoreEl%D(1:B%NKind)
      B%NTyp1PF%I(1:B%NKind)=BS%NTyp1PF%I(1:B%NKind)
      B%ProjEll%I(1:B%NKind)=BS%ProjEll%I(1:B%NKind)
      B%NTyp2PF%I(0:B%MxProjL,1:B%NKind)=BS%NTyp2PF%I(0:B%MxProjL,1:B%NKind)
      B%Typ1Ell%I(1:B%Typ1Fnk,1:B%NKind)=BS%Typ1Ell%I(1:B%Typ1Fnk,1:B%NKind)
      B%Typ1Exp%D(1:B%Typ1Fnk,1:B%NKind)=BS%Typ1Exp%D(1:B%Typ1Fnk,1:B%NKind)
      B%Typ1CCo%D(1:B%Typ1Fnk,1:B%NKind)=BS%Typ1CCo%D(1:B%Typ1Fnk,1:B%NKind)
      B%Typ2Ell%I(1:B%Typ2Fnk,0:B%MxProjL,1:B%NKind)=BS%Typ2Ell%I(1:B%Typ2Fnk,0:B%MxProjL,1:B%NKind)
      B%Typ2Exp%D(1:B%Typ2Fnk,0:B%MxProjL,1:B%NKind)=BS%Typ2Exp%D(1:B%Typ2Fnk,0:B%MxProjL,1:B%NKind)
      B%Typ2CCo%D(1:B%Typ2Fnk,0:B%MxProjL,1:B%NKind)=BS%Typ2CCo%D(1:B%Typ2Fnk,0:B%MxProjL,1:B%NKind)
    ENDIF
    !
    ! We need to reorder the Expnt for any case. It is important to have
    ! the Expnt in a decreassing order (skip out of the primitive loops
    ! for integral evaluations).
    CALL OrderExpnt(B)
    !
    !    IF(B%HasECPs)THEN
    !       PrintFlags%Set=DEBUG_BASISSET
    !       CALL PPrint(B,Unit_O=6)
    !    ENDIF
    ! Done with the temp BS
    BS%HasECPs=.TRUE.
    CALL Delete(BS)
    ! Compute blocking for sparse matrix methods
    CALL New(BlkSiz,B%NAtms)
    CALL New(OffSet,B%NAtms)
    CALL BlockBuild2(G,B,BlkSiz,OffSet)
    ParseBasisSets=.TRUE.
  END FUNCTION ParseBasisSets
  !
  SUBROUTINE OrderExpnt(B)
    TYPE(BSET)     :: B
    TYPE(INT_VECT) :: IA
    TYPE(DBL_RNK2) :: RA
    INTEGER        :: NK,CF,PF,LMN
    REAL(DOUBLE)   :: Tmp
    LOGICAL        :: IsOrdred,GlbOrder
    !-------------------------------------------------------------------------!
    CALL New(RA,(/B%LMNLen,B%NPrim/))
    CALL New(IA,B%NPrim)
    GlbOrder=.FALSE.
    DO NK=1,B%NKind
      DO CF=1,B%NCFnc%I(NK)
        Tmp=BIG_DBL
        IsOrdred=.TRUE.
        DO PF=1,B%NPFnc%I(CF,NK)
          IA%I(PF)=PF
          IF(B%Expnt%D(PF,CF,NK).GT.Tmp.AND.IsOrdred) IsOrdred=.FALSE.
          Tmp=B%Expnt%D(PF,CF,NK)
          !write(*,*) 'Before',B%Expnt%D(PF,CF,NK), &
          !       B%CCoef%D(B%LStrt%I(CF,NK):B%LStop%I(CF,NK),PF,CF,NK)
        ENDDO
        IF(IsOrdred) CYCLE
        GlbOrder=.TRUE.
        CALL DblIntSort77(B%NPFnc%I(CF,NK),B%Expnt%D(1,CF,NK),IA%I(1),-2)
        CALL DCOPY(B%LMNLen*B%NPrim,B%CCoef%D(1,1,CF,NK),1,RA%D(1,1),1)
        DO PF=1,B%NPFnc%I(CF,NK)
          DO LMN=B%LStrt%I(CF,NK),B%LStop%I(CF,NK)
            B%CCoef%D(LMN,PF,CF,NK)=RA%D(LMN,IA%I(PF))
          ENDDO
        ENDDO
        !do PF=1,B%NPFnc%I(CF,NK)
        !   write(*,*) 'After',B%Expnt%D(PF,CF,NK), &
        !          B%CCoef%D(B%LStrt%I(CF,NK):B%LStop%I(CF,NK),PF,CF,NK)
        !enddo
      ENDDO
    ENDDO
    IF(GlbOrder) THEN
      CALL MondoLog(DEBUG_NONE, "OrderExpnt", " ... ning Warning  Warning  Warning Wa ... ")
      CALL MondoLog(DEBUG_NONE, "OrderExpnt", "The primitive functions have been automatically reordered in a decreasing way!")
    ENDIF
    CALL Delete(IA)
    CALL Delete(RA)
  END SUBROUTINE OrderExpnt

  SUBROUTINE BlockBuild2(G,B,BS,OS)
    TYPE(CRDS)      :: G
    TYPE(BSET)      :: B
    TYPE(INT_VECT)  :: BS,OS
    INTEGER         :: NA,NK,NC,Stride
    !-------------------------------------------------------------------------------------!
    B%NBasF=0
    ! Off set starts at 1
    OS%I(1)=1
    DO NA=1,G%NAtms
      ! Block size is total number of basis functions per atom
      BS%I(NA)=0
      NK=G%AtTyp%I(NA)
      ! Go over contracted functions
      DO NC=1,B%NCFnc%I(NK)
        ! Add in size of each contraction
        Stride=B%LStop%I(NC,NK)-B%LStrt%I(NC,NK)+1
        BS%I(NA)=BS%I(NA)+Stride
        ! Oh yeah, accumulate basis function counter too...
        B%NBasF=B%NBasF+Stride
      ENDDO
      ! Off set counter from block sizes
      IF(NA.GE.2)OS%I(NA)=OS%I(NA-1)+BS%I(NA-1)
    ENDDO
  END SUBROUTINE BlockBuild2



  SUBROUTINE ReNormalizePrimitives(A)
    TYPE(BSET)       :: A
    TYPE(DBL_RNK2)   :: AuxCoef
    REAL(DOUBLE)     :: Z,C,Expnt,RNorm,ZA,ZB,CA,CB,SQNrm
    INTEGER          :: K,L,M,N,NC,NK,NP,LMN,MinL,MaxL,PFA,PFB,MaxSym,MaxPrm
    REAL(DOUBLE),PARAMETER, &      ! Fact=Sqrt[Pi](2*L-1)!! 2^(-L)
         DIMENSION(0:10)  :: Fact=(/0.17724538509055160273D1, &
         0.8862269254527580136D0,  &
         0.13293403881791370205D1, &
         0.3323350970447842551D1,  &
         0.11631728396567448929D2, &
         0.5234277778455352018D2,  &
         0.28788527781504436100D3, &
         0.1871254305797788347D4,  &
         0.14034407293483412599D5, &
         0.1192924619946090071D6,  &
         0.11332783889487855673D7/)
    !-------------------------------------------------------------------------------------!
    ! Allocate temp space
    MaxSym=0
    MaxPrm=0
    DO K=1,A%NKind
      DO NC=1,A%NCFnc%I(K)
        MaxSym=MAX(MaxSym,A%ASymm%I(2,NC,K))
        MaxPrm=MAX(MaxPrm,A%NPFnc%I(NC,K))
      ENDDO
    ENDDO
    CALL New(AuxCoef,(/MaxSym+1,MaxPrm/))
    DO K=1,A%NKind
      DO NC=1,A%NCFnc%I(K)
        MinL=A%ASymm%I(1,NC,K)
        MaxL=A%ASymm%I(2,NC,K)
        DO NP=1,A%NPFnc%I(NC,K)
          Z=A%Expnt%D(NP,NC,K)
          DO L=MinL,MaxL
            AuxCoef%D(L+1,NP)=A%CCoef%D(L+1,NP,NC,K)*Z**(Half*DBLE(L)+0.75D0)
          ENDDO
        ENDDO
        DO LMN=A%LStrt%I(NC,K),A%LStop%I(NC,K)
          RNorm=0.0D0
          L=A%LxDex%I(LMN)
          M=A%LyDex%I(LMN)
          N=A%LzDex%I(LMN)
          DO PFA=1,A%NPFnc%I(NC,K)
            ZA=A%Expnt%D(PFA,NC,K)
            CA=AuxCoef%D(L+M+N+1,PFA)
            DO PFB=1,A%NPFnc%I(NC,K)
              ZB=A%Expnt%D(PFB,NC,K)
              CB=AuxCoef%D(L+M+N+1,PFB)
              RNorm=RNorm+CA*CB*(ZA+ZB)**(-DBLE(L+M+N)-1.5D0)
            ENDDO
          ENDDO
          SqNrm=1.0D0/SQRT(RNorm*Fact(L)*Fact(M)*Fact(N))
          DO NP=1,A%NPFnc%I(NC,K)
            A%CCoef%D(LMN,NP,NC,K)=AuxCoef%D(L+M+N+1,NP)*SqNrm
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    CALL Delete(AuxCoef)
  END SUBROUTINE ReNormalizePrimitives

  SUBROUTINE ParseBasisNames(B)
    TYPE(BasisSets)                 :: B
    INTEGER                         :: I,J
    TYPE(CHR_VECT)                  :: Chars
    CHARACTER(LEN=DCL)              :: Line
    CHARACTER(LEN=BASESET_CHR_LEN)  :: LowC1,LowC2
    LOGICAL                         :: Found

    CALL Align(BASIS_SETS,Inp)
    BACKSPACE(Inp)
    READ(Inp,DEFAULT_CHR_FMT)Line
    CALL RemoveComments(Line)
    CALL MondoLog(DEBUG_MAXIMUM, "ParseBasisNames", "parsing "//TRIM(Line))
    CALL LineToChars(Line,Chars)
    B%NBSets=SIZE(Chars%C)-1
    DO I=2,SIZE(Chars%C)
      Found=.FALSE.
      DO J=1,NSupSets
        LowC1=Chars%C(I)
        LowC2=CSets(1,J)
        CALL LowCase(LowC1)
        CALL LowCase(LowC2)
        IF(TRIM(LowC1)==TRIM(LowC2))THEN
          B%BName(I-1)=ADJUSTL(CSets(2,J))
          Found=.TRUE.
          EXIT
        ENDIF
      ENDDO
      IF(.NOT.Found)THEN
        B%BName(I-1)=TRIM(Chars%C(I))
      ENDIF
    ENDDO
    CALL Delete(Chars)
  END SUBROUTINE ParseBasisNames

  !=============================================================================
  ! PRECOMPUTE PRIMITIVE DISTRIBUTION INFORMATION
  !=============================================================================
  SUBROUTINE PrimitiveReDistribution(BS,NExpt,DExpt,Lndex)
    TYPE(BSET),     INTENT(INOUT):: BS
    TYPE(DBL_VECT), INTENT(OUT)  :: DExpt
    TYPE(INT_VECT), INTENT(OUT)  :: Lndex
    INTEGER,        INTENT(OUT)  :: NExpt
    TYPE(INT_VECT)               :: ITmp,IPnt
    INTEGER                      :: I,K,KA,KB,CFA,CFB,PFA,PFB,NPB
    !---------------------------------------------------------------------------!
    ! First count the types of primitive distributions(K=1), then compute
    ! their primitive distribution exponents and associated values of
    ! angular symmetry (K=2)
    DO K=1,2
      IF(K==1)THEN
        NExpt=1
      ELSE
        CALL New(DExpt,NExpt)
        CALL New(Lndex,NExpt)
        NExpt=1
      ENDIF
      DO KA=1,BS%NKind
        DO KB=1,KA-1
          DO CFA=1,BS%NCFnc%I(KA)
            DO CFB=1,BS%NCFnc%I(KB)
              DO PFA=1,BS%NPFnc%I(CFA,KA)
                DO PFB=1,BS%NPFnc%I(CFB,KB)
                  IF(K==2)THEN
                    DExpt%D(NExpt)=BS%Expnt%D(PFA,CFA,KA)  &
                         +BS%Expnt%D(PFB,CFB,KB)
                    Lndex%I(NExpt)=BS%ASymm%I(2,CFA,KA)    &
                         +BS%ASymm%I(2,CFB,KB)
                  ENDIF
                  NExpt=NExpt+1
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO KA=1,BS%NKind
        DO CFA=1,BS%NCFnc%I(KA)
          DO CFB=1,CFA
            DO PFA=1,BS%NPFnc%I(CFA,KA)
              NPB=BS%NPFnc%I(CFB,KA)
              IF(CFA.EQ.CFB)NPB=PFA
              DO PFB=1,NPB
                IF(K==2)THEN
                  DExpt%D(NExpt)=BS%Expnt%D(PFA,CFA,KA)  &
                       +BS%Expnt%D(PFB,CFB,KA)
                  Lndex%I(NExpt)=BS%ASymm%I(2,CFA,KA)    &
                       +BS%ASymm%I(2,CFB,KA)
                ENDIF
                NExpt=NExpt+1
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    ! Add exponent for nuclear delta function
    DExpt%D(NExpt)=NuclearExpnt
    Lndex%I(NExpt)=0
    ! Sort the exponents in assending order and carry
    ! allong angular symetry indecies
    CALL New(IPnt,NExpt)
    CALL New(ITmp,NExpt)
    DO I=1,NExpt
      IPnt%I(I)=I
    ENDDO
    CALL Sort(DExpt,IPnt,NExpt,2)
    DO I=1,NExpt
      ITmp%I(I)=Lndex%I(IPnt%I(I))
    ENDDO
    DO I=1,NExpt
      Lndex%I(I)=ITmp%I(I)
    ENDDO
    ! Tidy up
    CALL Delete(ITmp)
    CALL Delete(IPnt)
  END SUBROUTINE PrimitiveReDistribution
  !
  !
  SUBROUTINE ParseMaxElemBCSR(NBSets,Unit,ArrMaxNon0s,ArrMaxNBlks)
    !H---------------------------------------------------------------------------------
    !H SUBROUTINE ParseMaxElemBCSR(NBSets,Unit,ArrMaxNon0s,ArrMaxNBlks)
    !H  This routine look in the unit <Unit> for the keywords <MaxNon0s> and
    !H  <MaxNBlks>. The new value will replace the value of the estimate
    !H  dimension in <BCSRDimensions> if the key <Default> is not used.
    !H   e.g. MaxNon0s=(Default,Default,100000).
    !H        MaxNBlks=(Default,Default,   900).
    !H        MaxNon0s=(100000) or MaxNon0s=100000.
    !H        MaxNBlks=(   900) or MaxNBlks=   900.
    !H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    INTEGER       , INTENT(IN   ) :: NBSets,Unit
    TYPE(INT_VECT), INTENT(INOUT) :: ArrMaxNon0s,ArrMaxNBlks
    !-------------------------------------------------------------------
    INTEGER                       :: iBS
    TYPE(CHR_VECT)                :: Arg
    !-------------------------------------------------------------------
    !
    ! Look for MaxNon0s
    IF(OptGetKeyArg(Unit,BCSR_MAXNON0S,Arg)) THEN
      IF(SIZE(Arg%C).NE.NBSets) &
           CALL Halt('The number of MaxNon0s arguments <'  //TRIM(IntToChar(SIZE(Arg%C)))// &
           '> must be the same as the number of BasisSets <'//TRIM(IntToChar(NBSets))//'>.')

      ! Copy info.
      DO iBS=1,NBSets
        IF(TRIM(Arg%C(iBS)).EQ.BCSR_MAXDEFAULT) THEN
          ArrMaxNon0s%I(iBS)=-BIG_INT
        ELSE

          ! Check for strange char.
          IF(.NOT.ChrCkkIfInt(Arg%C(iBS))) &
               CALL Halt('The argument <'//TRIM(Arg%C(iBS))// &
               '> is not a positive integer nor <Default> in the Key <MaxNon0s>.')
          ArrMaxNon0s%I(iBS)=CharToInt(TRIM(Arg%C(iBS)))
        ENDIF
        !write(*,*) 'ArrMaxNon0s%I(iBS)',ArrMaxNon0s%I(iBS)
      ENDDO

      CALL Delete(Arg)
    ENDIF

    ! Look for MaxNBlks
    IF(OptGetKeyArg(Unit,BCSR_MAXNBLKS,Arg)) THEN
      IF(SIZE(Arg%C).NE.NBSets) &
           CALL Halt('The number of MaxNBlks arguments <'  //TRIM(IntToChar(SIZE(Arg%C)))// &
           '> must be the same as the number of BasisSets <'//TRIM(IntToChar(NBSets))//'>.')

      ! Copy info.
      DO iBS=1,NBSets
        IF(TRIM(Arg%C(iBS)).EQ.BCSR_MAXDEFAULT) THEN
          ArrMaxNon0s%I(iBS)=-BIG_INT
        ELSE

          ! Check for strange char.
          IF(.NOT.ChrCkkIfInt(Arg%C(iBS))) &
               CALL Halt('The argument <'//TRIM(Arg%C(iBS))// &
               '> is not a positive integer nor <Default> in the Key <MaxNBlks>.')
          ArrMaxNBlks%I(iBS)=CharToInt(TRIM(Arg%C(iBS)))
        ENDIF
        !write(*,*) 'ArrMaxNBlks%I(iBS)',ArrMaxNBlks%I(iBS)
      ENDDO

      CALL Delete(Arg)
    ENDIF

  END SUBROUTINE ParseMaxElemBCSR

END MODULE ParseBasis
