MODULE ParsePeriodic
  USE Parse
  USE InOut
  USE AtomPairs 
  USE PrettyPrint
  USE OptionKeys
  USE PeriodicKeys
  USE ControlStructures
CONTAINS
!=========================================================================
!
!=========================================================================
  SUBROUTINE LoadPeriodic(N,O,G,P)
    TYPE(FileNames)  :: N
    TYPE(Options)    :: O
    TYPE(Geometries) :: G
    TYPE(Periodics)  :: P
    TYPE(PBCInfo)    :: PBC
    INTEGER          :: I,GBeg,GEnd
!-----------------------------------------------------------------------!
!   If we are restarting, just use values read from HDF ...
    CALL OpenASCII(N%IFile,Inp)
    IF(O%Guess==GUESS_EQ_RESTART .AND. .NOT.OptKeyQ(Inp,RESTART_OPTION,RESTART_NEWGEOM))THEN
       CLOSE(Inp)
       RETURN
    ENDIF
    IF(O%Guess==GUESS_EQ_RESTART .AND. OptKeyQ(Inp,RESTART_OPTION,RESTART_NEWGEOM)) THEN
!       WRITE(*,*)' REPARSING PERIODIC ON RESTART'
    ENDIF
    CALL New(PBC)
    CALL LoadPeriodicOptions(PBC)
    CALL LoadLattice(PBC)
    CALL UnitCellSetUp(PBC)
    CLOSE(UNIT=Inp,STATUS='KEEP')
    IF(O%Grad==GRAD_TS_SEARCH_NEB)THEN
       GBeg=0
       GEnd=G%Clones+1
    ELSE	
       GBeg=1
       GEnd=G%Clones
    ENDIF
    DO I=GBeg,GEnd
       G%Clone(I)%PBC=PBC
       CALL CalculateCoordArrays(G%Clone(I))
    ENDDO
  END SUBROUTINE LoadPeriodic
!=========================================================================
!
!=========================================================================
  SUBROUTINE LoadPeriodicOptions(PBC)
    TYPE(PBCInfo)    :: PBC
    INTEGER          :: I,J,NTot,MaxEll
!-----------------------------------------------------------------------!
!   Parse the coordinate type 
!
    IF(OptKeyQ(Inp,PBOUNDRY,CRT_FRAC))THEN
       PBC%InAtomCrd=.FALSE.
    ELSE
       PBC%InAtomCrd=.TRUE.
    ENDIF
!   Parse Periodic Directions
    Ntot = 0
    PBC%AutoW%I=BIG_INT
    IF(FindKey(PBCWRAP,Inp))THEN
       IF(OptKeyLocQ(Inp,PBCWRAP,PBC_TRUE,MaxSets,NLoc,Location)) THEN
          Ntot = NLoc
          DO I=1,NLoc
             PBC%AutoW%I(Location(I)) = 1
          ENDDO
       ENDIF
       PBC%Dimen=NLoc
       IF(OptKeyLocQ(Inp,PBCWRAP,PBC_FALSE,MaxSets,NLoc,Location)) THEN
          Ntot = NTot+NLoc
          DO I=1,NLoc
             PBC%AutoW%I(Location(I)) = 0
          ENDDO
       ENDIF
    ELSE
       PBC%AutoW%I(1:3)=0
       PBC%Dimen=0
    ENDIF
!   Parse Translate
    IF(OptKeyQ(Inp,PBOUNDRY,CENTERATOMS))THEN
       PBC%Translate=.TRUE.
    ELSE
       PBC%Translate=.FALSE.
    ENDIF
!   Intput over-ride
    IF(OptKeyQ(Inp,PBOUNDRY,PFFOVRDE))THEN
       PBC%PFFOvRide=.TRUE.
    ELSE
       PBC%PFFOvRide=.FALSE.
    ENDIF
!   Intput maxium Ell and number of layers
    IF(PBC%PFFOvRide)THEN
       IF(.NOT.OptIntQ(Inp,PFFMXLAY,PBC%PFFMaxLay))THEN
          CALL MondoHalt(PRSE_ERROR,'PFFOverRide is on, please provide PFFMaxLay')
       ENDIF
       IF(.NOT.OptIntQ(Inp,PFFMXELL,PBC%PFFMaxEll))THEN
          CALL MondoHalt(PRSE_ERROR,'PFFOverRide is on, please provide PFFMaxEll')
       ENDIF
    ELSE      
       PBC%PFFMaxLay=1
!      Look for MaxEll in periodic boundary options
       IF(.NOT.OptIntQ(Inp,PFFMXELL,PBC%PFFMaxEll)) THEN 
          PBC%PFFMaxEll=10
       ENDIF  
    ENDIF
!   Parse permeability 
    IF(.NOT.OptDblQ(Inp,EPSILON,PBC%Epsilon))THEN
       PBC%Epsilon=BIG_DBL !1.D32
    ENDIF
!   Parse Atom Wrap, Default is on
    IF(OptKeyQ(Inp,PBOUNDRY,ATOMW_OFF))THEN
       PBC%AtomW=.FALSE.
    ELSE
       PBC%AtomW=.TRUE.
    ENDIF
!
!    WRITE(*,*)   'PBC%PFFMaxEll = ',PBC%PFFMaxEll
!    WRITE(*,*)   'PBC%PFFMaxLay = ',PBC%PFFMaxLay
!
  END SUBROUTINE LoadPeriodicOptions
!============================================================================
!
!============================================================================
  SUBROUTINE LoadLattice(PBC)
    TYPE(PBCInfo)                   :: PBC
    INTEGER                         :: NLvec,NTvec,Dimen,I,J
    CHARACTER(LEN=2)                :: At
    CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Line,LowLine
    REAL(DOUBLE),PARAMETER          :: DegToRad = 1.745329251994329576923D-2
    REAL(DOUBLE),DIMENSION(6)       :: Vec
    REAL(DOUBLE)                    :: AngAB,AngAC,AngBC,Error
!

    PBC%TransVec%D=Zero
    PBC%BoxShape%D=Zero  
    DO I=1,3
       PBC%BoxShape%D(I,I)=One
    ENDDO
    IF(PBC%Dimen==0) RETURN
!
    IF(.NOT. FindKey(BEGIN_PERIODIC,Inp)) THEN
       CALL MondoHalt(PRSE_ERROR,'Lattice Vectors Must Be Suppied for Periodic Systems')
    ENDIF
!
    NLvec=0 
    CALL Align(BEGIN_PERIODIC,Inp)
    DO 
       READ(Inp,DEFAULT_CHR_FMT,END=1) Line
       IF(INDEX(Line,END_PERIODIC)==0) THEN
          LowLine=Line
          CALL LowCase(LowLine) 
          J=SCAN(LowLine,Lower)
          IF(J/=0)THEN
             CALL LineToGeom(Line,At,Vec)
             IF(At==ALAT_VEC) THEN
                NLvec = NLvec+1
                PBC%BoxShape%D(1:3,1)=Vec(1:3)
             ELSEIF(At==BLAT_VEC) THEN
                NLvec = NLvec+1
                PBC%BoxShape%D(1:3,2)=Vec(1:3)
             ELSEIF(At==CLAT_VEC) THEN
                NLvec = NLvec+1
                PBC%BoxShape%D(1:3,3)=Vec(1:3)
             ENDIF
          ELSE
             NLvec=3
             CALL LineToDbls(Line,6,Vec)
             J=0
             DO I=1,6
                IF(ABS(Vec(I)) .GT. 1.D-14) THEN
                   J=J+1
                ENDIF
             ENDDO
             IF(PBC%InAtomCrd) THEN
                IF(J==1) THEN
                   IF(PBC%Dimen > 1) THEN
                      CALL MondoHalt(PRSE_ERROR,'Number of Magnitudes and Angles supplied is incorrect')
                   ENDIF
                   IF(PBC%AutoW%I(1)==1) PBC%BoxShape%D(1,1) = Vec(1)
                   IF(PBC%AutoW%I(2)==1) PBC%BoxShape%D(2,2) = Vec(1)
                   IF(PBC%AutoW%I(3)==1) PBC%BoxShape%D(3,3) = Vec(1)
                ELSEIF(J==3) THEN 
                   IF(PBC%Dimen > 2) THEN
                      CALL MondoHalt(PRSE_ERROR,'Number of Magnitudes and Angles supplied is incorrect')
                   ENDIF
                   IF(PBC%AutoW%I(3)==0) THEN
                      PBC%BoxShape%D(1,1)=Vec(1)
                      PBC%BoxShape%D(1,2)=Vec(2)*COS(DegToRad*Vec(3))
                      PBC%BoxShape%D(2,2)=Vec(2)*SIN(DegToRad*Vec(3)) 
                   ELSEIF(PBC%AutoW%I(2)==0) THEN
                      PBC%BoxShape%D(1,1)=Vec(1)
                      PBC%BoxShape%D(1,3)=Vec(2)*COS(DegToRad*Vec(3))
                      PBC%BoxShape%D(3,3)=Vec(2)*SIN(DegToRad*Vec(3))
                   ELSEIF(PBC%AutoW%I(1)==0) THEN
                      PBC%BoxShape%D(2,2)=Vec(1)
                      PBC%BoxShape%D(2,3)=Vec(2)*COS(DegToRad*Vec(3))
                      PBC%BoxShape%D(3,3)=Vec(2)*SIN(DegToRad*Vec(3))
                   ENDIF
                ELSEIF(J==6) THEN
                   PBC%BoxShape%D(1,1)=Vec(1)
                   PBC%BoxShape%D(1,2)=Vec(2)*COS(DegToRad*Vec(6))
                   PBC%BoxShape%D(2,2)=Vec(2)*SIN(DegToRad*Vec(6))
                   PBC%BoxShape%D(1,3)=Vec(3)*COS(DegToRad*Vec(5))
                   PBC%BoxShape%D(2,3)=(Vec(2)*Vec(3)*COS(DegToRad*Vec(4)) &
                        -PBC%BoxShape%D(1,2)*PBC%BoxShape%D(1,3))/PBC%BoxShape%D(2,2)
                   PBC%BoxShape%D(3,3)=SQRT(Vec(3)**2-PBC%BoxShape%D(1,3)**2-PBC%BoxShape%D(2,3)**2)
                ELSE
                   CALL MondoHalt(PRSE_ERROR,'Number of Magnitudes and Angles supplied is incorrect')
                ENDIF
             ELSE
                IF(J==6) THEN
                   PBC%BoxShape%D(1,1)=Vec(1)
                   PBC%BoxShape%D(1,2)=Vec(2)*COS(DegToRad*Vec(6))
                   PBC%BoxShape%D(2,2)=Vec(2)*SIN(DegToRad*Vec(6))
                   PBC%BoxShape%D(1,3)=Vec(3)*COS(DegToRad*Vec(5))
                   PBC%BoxShape%D(2,3)=(Vec(2)*Vec(3)*COS(DegToRad*Vec(4)) &
                        -PBC%BoxShape%D(1,2)*PBC%BoxShape%D(1,3))/PBC%BoxShape%D(2,2)
                   PBC%BoxShape%D(3,3)=SQRT(Vec(3)**2-PBC%BoxShape%D(1,3)**2-PBC%BoxShape%D(2,3)**2)
                ELSE
                   CALL MondoHalt(PRSE_ERROR,'In fractional coordinates all lattice vectors must be supplied')
                ENDIF
             ENDIF
          ENDIF
       ELSE
          EXIT
       ENDIF
    ENDDO
    IF(NLVec .NE. 3) THEN
       CALL MondoHalt(PRSE_ERROR,'Lattice Vectors are incorrect')
    ENDIF
    DO I=1,3
       DO J=1,3
          IF(ABS(PBC%BoxShape%D(I,J)).LT. 1.D-12) PBC%BoxShape%D(I,J)=Zero
       ENDDO
    ENDDO
    RETURN
1   CALL Halt('While parsing '//TRIM(InpFile)//', failed to find '     &
         //TRIM(END_PERIODIC)//'. You may be missing blank '  &
         //'line at the end of the inPut file.')
  END SUBROUTINE LoadLattice
!============================================================================
!
!============================================================================
  SUBROUTINE UnitCellSetUp(PBC)
    TYPE(PBCInfo)       :: PBC
    INTEGER             :: I,J,K,NLvec,NTvec
!   CalculatetheBoxVolume
    PBC%CellVolume=ABS(CellVolume(PBC%BoxShape%D,PBC%AutoW%I))
!   Calculate dipole and quadripole factors
    IF(PBC%Dimen<2)THEN
       PBC%DipoleFAC=Zero
       PBC%QupoleFAC=Zero
    ELSEIF(PBC%Dimen==2)THEN
       PBC%DipoleFAC=(Four*Pi/PBC%CellVolume)*(One/(PBC%Epsilon+One))
       PBC%QupoleFAC=Zero
       IF(ABS(PBC%DipoleFAC).LT.1.D-14) PBC%DipoleFAC=Zero
    ELSEIF(PBC%Dimen==3)THEN
       PBC%DipoleFAC=-(Four*Pi/PBC%CellVolume)*(One/Three-One/(Two*PBC%Epsilon+One))
       PBC%QupoleFAC=(Two*Pi/PBC%CellVolume)*(One/Three-One/(Two*PBC%Epsilon+One))
       IF(ABS(PBC%DipoleFAC).LT.1.D-14) PBC%DipoleFAC=Zero
       IF(ABS(PBC%QupoleFAC).LT.1.D-14) PBC%QupoleFAC=Zero
    ENDIF
!   Find the center of the cell
    DO I=1,3
       PBC%CellCenter%D(I)=Zero
       IF(PBC%AutoW%I(I)==1)THEN
          DO J=1,3
             IF(PBC%AutoW%I(J)==1)THEN
                PBC%CellCenter%D(I)=PBC%CellCenter%D(I)+Half*PBC%BoxShape%D(I,J)
             ENDIF
          ENDDO
       ENDIF
    ENDDO
!   Compute the inverse box shape
    PBC%InvBoxSh%D = InverseMatrix(PBC%BoxShape%D)
  END SUBROUTINE UnitCellSetUp
!=========================================================================
! ACCOUNT FOR ALL FOUR (YEP COUNT EM FOUR) COORDINATE ARRAYS... YEESH!
!=========================================================================
  SUBROUTINE CalculateCoordArrays(G)
    TYPE(CRDS) :: G
    INTEGER    :: I
!-----------------------------------------------------------------------!
    IF(G%PBC%InAtomCrd)THEN
       G%Carts%D=G%AbCarts%D
       CALL CalFracCarts(G)
       IF(G%PBC%AtomW) THEN 
          CALL WrapAtoms(G)
       ENDIF
    ELSE
       G%BoxCarts%D=G%AbCarts%D
       CALL CalAtomCarts(G)
       G%AbCarts%D=G%Carts%D
       IF(G%PBC%AtomW) THEN 
          CALL WrapAtoms(G)
       ENDIF
!      Convert the Velocities from Fractional to Atomic
       DO I=1,G%NAtms
          G%Velocity%D(:,I)   = FracToAtom(G,G%Velocity%D(:,I))
       ENDDO
    ENDIF
  END SUBROUTINE CalculateCoordArrays
!-----------------------------------------------------------------------!
END MODULE ParsePeriodic
