MODULE ParsePeriodic
  USE Parse
  USE InOut
  USE AtomPairs 
  USE PrettyPrint
  USE PeriodicKeys
  USE ControlStructures
CONTAINS
!=========================================================================
!
!=========================================================================
  SUBROUTINE LoadPeriodic(N,G,P)
    TYPE(FileNames)  :: N
    TYPE(Geometries) :: G
    TYPE(Periodics)  :: P
    TYPE(PBCInfo)    :: PBC
    INTEGER          :: I
!-----------------------------------------------------------------------!
    CALL OpenASCII(N%IFile,Inp)
    CALL LoadPeriodicOptions(PBC)
    CALL LoadLattice(PBC)
    CALL UnitCellSetUp(PBC)
    CLOSE(UNIT=Inp,STATUS='KEEP')
    DO I=1,G%Clones
       G%Clone(I)%PBC=PBC
       CALL YeeshFourCoordinateArrays(G%Clone(I))
    ENDDO
  END SUBROUTINE LoadPeriodic
!=========================================================================
!
!=========================================================================
  SUBROUTINE LoadPeriodicOptions(PBC)
    TYPE(PBCInfo)    :: PBC
    INTEGER          :: I,J,MaxEll
!-----------------------------------------------------------------------!
!   Parse the coordinate type 
    IF(OptKeyQ(Inp,PBOUNDRY,CRT_FRAC))THEN
       PBC%InAtomCrd=.FALSE.
    ELSE
       PBC%InAtomCrd=.TRUE.
    ENDIF
!   Parse Periodic Directions
    Ntot = 0
    IF(FindKey(PBCWRAP,Inp))THEN
       IF(OptKeyLocQ(Inp,PBCWRAP,PBC_TRUE,MaxSets,NLoc,Location)) THEN
          Ntot = NLoc
          DO I=1,NLoc
             PBC%AutoW(Location(I)) = .TRUE.
          ENDDO
       ENDIF
       PBC%Dimen=NLoc
       IF(OptKeyLocQ(Inp,PBCWRAP,PBC_FALSE,MaxSets,NLoc,Location)) THEN
          Ntot = NTot+NLoc
          DO I=1,NLoc
             PBC%AutoW(Location(I)) = .FALSE.
          ENDDO
       ENDIF
       IF(NTot .NE. 3) THEN
          CALL MondoHalt(PRSE_ERROR,'PBC = (?,?,?): Three Logicals must be Specified')
       ENDIF
    ELSE
       PBC%AutoW(:) = .FALSE.
       PBC%Dimen=0
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
       IF(.NOT.OptIntQ(Inp,PBOUNDRY,PBC%PFFMaxEll))  PBC%PFFMaxEll=16
    ENDIF
!   Parse permeability 
    IF(.NOT.OptDblQ(Inp,EpsILON,PBC%Epsilon))THEN
       PBC%Epsilon=BIG_DBL !1.D32
    ENDIF
  END SUBROUTINE LoadPeriodicOptions
!============================================================================
!
!============================================================================
  SUBROUTINE LoadLattice(PBC)!,NLvec,NTvec)
    TYPE(PBCInfo)                   :: PBC
    INTEGER                         :: NLvec,NTvec,Dimen,I,J
    CHARACTER(LEN=2)                :: At
    CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Line,LowLine
    REAL(DOUBLE),PARAMETER          :: DegToRad = 1.745329251994329576923D-2
    REAL(DOUBLE),DIMENSION(6)       :: Vec
    REAL(DOUBLE)                    :: AngAB,AngAC,AngBC,Error
!

    PBC%TransVec=Zero
    PBC%BoxShape=Zero  
    DO I=1,3
       PBC%BoxShape(I,I)=One
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
                PBC%BoxShape(1:3,1)=Vec(1:3)
             ELSEIF(At==BLAT_VEC) THEN
                NLvec = NLvec+1
                PBC%BoxShape(1:3,2)=Vec(1:3)
             ELSEIF(At==CLAT_VEC) THEN
                NLvec = NLvec+1
                PBC%BoxShape(1:3,3)=Vec(1:3)
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
                   IF(PBC%AutoW(1)) PBC%BoxShape(1,1) = Vec(1)
                   IF(PBC%AutoW(2)) PBC%BoxShape(2,2) = Vec(1)
                   IF(PBC%AutoW(3)) PBC%BoxShape(3,3) = Vec(1)
                ELSEIF(J==3) THEN 
                   IF(PBC%Dimen > 2) THEN
                      CALL MondoHalt(PRSE_ERROR,'Number of Magnitudes and Angles supplied is incorrect')
                   ENDIF
                   IF(.NOT. PBC%AutoW(3)) THEN
                      PBC%BoxShape(1,1)=Vec(1)
                      PBC%BoxShape(1,2)=Vec(2)*COS(DegToRad*Vec(3))
                      PBC%BoxShape(2,2)=Vec(2)*SIN(DegToRad*Vec(3)) 
                   ELSEIF(.NOT. PBC%AutoW(2)) THEN
                      PBC%BoxShape(1,1)=Vec(1)
                      PBC%BoxShape(1,3)=Vec(2)*COS(DegToRad*Vec(3))
                      PBC%BoxShape(3,3)=Vec(2)*SIN(DegToRad*Vec(3))
                   ELSEIF(.NOT. PBC%AutoW(1)) THEN
                      PBC%BoxShape(2,2)=Vec(1)
                      PBC%BoxShape(2,3)=Vec(2)*COS(DegToRad*Vec(3))
                      PBC%BoxShape(3,3)=Vec(2)*SIN(DegToRad*Vec(3))
                   ENDIF
                ELSEIF(J==6) THEN
                   PBC%BoxShape(1,1)=Vec(1)
                   PBC%BoxShape(1,2)=Vec(2)*COS(DegToRad*Vec(6))
                   PBC%BoxShape(2,2)=Vec(2)*SIN(DegToRad*Vec(6))
                   PBC%BoxShape(1,3)=Vec(3)*COS(DegToRad*Vec(5))
                   PBC%BoxShape(2,3)=(Vec(2)*Vec(3)*COS(DegToRad*Vec(4)) &
                        -PBC%BoxShape(1,2)*PBC%BoxShape(1,3))/PBC%BoxShape(2,2)
                   PBC%BoxShape(3,3)=SQRT(Vec(3)**2-PBC%BoxShape(1,3)**2-PBC%BoxShape(2,3)**2)
                ELSE
                   CALL MondoHalt(PRSE_ERROR,'Number of Magnitudes and Angles supplied is incorrect')
                ENDIF
             ELSE
                IF(J==6) THEN
                   PBC%BoxShape(1,1)=Vec(1)
                   PBC%BoxShape(1,2)=Vec(2)*COS(DegToRad*Vec(6))
                   PBC%BoxShape(2,2)=Vec(2)*SIN(DegToRad*Vec(6))
                   PBC%BoxShape(1,3)=Vec(3)*COS(DegToRad*Vec(5))
                   PBC%BoxShape(2,3)=(Vec(2)*Vec(3)*COS(DegToRad*Vec(4)) &
                        -PBC%BoxShape(1,2)*PBC%BoxShape(1,3))/PBC%BoxShape(2,2)
                   PBC%BoxShape(3,3)=SQRT(Vec(3)**2-PBC%BoxShape(1,3)**2-PBC%BoxShape(2,3)**2)
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
          IF(ABS(PBC%BoxShape(I,J)).LT. 1.D-12) PBC%BoxShape(I,J)=Zero
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
    PBC%CellVolume=One
    DO I=1,3
       IF(PBC%AutoW(I))THEN
          PBC%CellVolume=PBC%CellVolume*PBC%BoxShape(I,I)
       ENDIF
    ENDDO
!   Calculate dipole and quadripole factors
    IF(PBC%Dimen<2)THEN
       PBC%DipoleFAC=Zero
       PBC%QupoleFAC=Zero
    ELSEIF(PBC%Dimen==2)THEN
       PBC%DipoleFAC=(Four*Pi/PBC%CellVolume)*(One/(PBC%Epsilon+One))
       PBC%QupoleFAC=Zero
       IF(ABS(PBC%DipoleFAC).LT.1.D-14)PBC%DipoleFAC=Zero
    ELSEIF(PBC%Dimen==3)THEN
       PBC%DipoleFAC=-(Four*Pi/PBC%CellVolume)*(One/Three-One/(Two*PBC%Epsilon+One))
       PBC%QupoleFAC=(Two*Pi/PBC%CellVolume)*(One/Three-One/(Two*PBC%Epsilon+One))
       IF(ABS(PBC%DipoleFAC).LT.1.D-14)PBC%DipoleFAC=Zero
       IF(ABS(PBC%QupoleFAC).LT.1.D-14)PBC%QupoleFAC=Zero
    ENDIF
!   Find the center of the cell
    DO I=1,3
       PBC%CellCenter(I)=Zero
       IF(PBC%AutoW(I))THEN
          DO J=1,3
             IF(PBC%AutoW(J))THEN
                PBC%CellCenter(I)=PBC%CellCenter(I)+Half*PBC%BoxShape(I,J)
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    PBC%InvBoxSh=Zero
!   Compute the inverse box shape 
    PBC%InvBoxSh(1,1)=One/PBC%BoxShape(1,1)
    PBC%InvBoxSh(2,1)=Zero
    PBC%InvBoxSh(3,1)=Zero
    PBC%InvBoxSh(1,2)=-PBC%BoxShape(1,2)/(PBC%BoxShape(1,1)*PBC%BoxShape(2,2))
    PBC%InvBoxSh(2,2)=One/PBC%BoxShape(2,2)
    PBC%InvBoxSh(3,2)=Zero
    PBC%InvBoxSh(1,3)=(PBC%BoxShape(1,2)*PBC%BoxShape(2,3)&
         -PBC%BoxShape(2,2)*PBC%BoxShape(1,3))&
         /(PBC%BoxShape(1,1)*PBC%BoxShape(2,2)*PBC%BoxShape(3,3))
    PBC%InvBoxSh(2,3)=-PBC%BoxShape(2,3)/(PBC%BoxShape(2,2)*PBC%BoxShape(3,3))
    PBC%InvBoxSh(3,3)=One/PBC%BoxShape(3,3)
!
  END SUBROUTINE UnitCellSetUp
!=========================================================================
! ACCOUNT FOR ALL FOUR (YEP COUNT EM FOUR) COORDINATE ARRAYS... YEESH!
!=========================================================================
  SUBROUTINE YeeshFourCoordinateArrays(G)
    TYPE(CRDS) :: G
    INTEGER    :: I
!-----------------------------------------------------------------------!
    IF(G%PBC%InAtomCrd)THEN
       CALL WrapAtoms(G)
    ELSE
! These are the two fractional coordinate arrays ... 
       G%BoxCarts%D=G%AbCarts%D
       G%AbBoxCarts%D=G%AbCarts%D
! ... and here are the two Cartesian coordinate arrays
       DO I=1,G%NAtms
          G%Carts%D(:,I)=FracToAtom(G,G%BoxCarts%D(:,I))
          G%AbCarts%D(:,I)=FracToAtom(G,G%AbBoxCarts%D(:,I))          
       ENDDO
    ENDIF
  END SUBROUTINE YeeshFourCoordinateArrays
!-----------------------------------------------------------------------!
END MODULE ParsePeriodic
