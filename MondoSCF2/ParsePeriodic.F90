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
    INTEGER          :: I,MaxEll
    !-----------------------------------------------------------------------!
    ! Parse the coordinate type 
    IF(OptKeyQ(Inp,PBOUNDRY,CRT_FRAC))THEN
       PBC%InAtomCrd=.FALSE.
    ELSE
       PBC%InAtomCrd=.TRUE.
    ENDIF
    ! Intput over-ride
    IF(OptKeyQ(Inp,PBOUNDRY,PFFOVRDE))THEN
       PBC%PFFOvRide=.TRUE.
    ELSE
       PBC%PFFOvRide=.FALSE.
    ENDIF
    !Intput maxium Ell and number of layers
    IF(PBC%PFFOvRide)THEN
       IF(.NOT.OptIntQ(Inp,PFFMXLAY,PBC%PFFMaxLay))THEN
          CALL MondoHalt(PRSE_ERROR,'PFFOverRide is on, please provide PFFMaxLay')
       ENDIF
       IF(.NOT.OptIntQ(Inp,PFFMXELL,PBC%PFFMaxEll))THEN
          CALL MondoHalt(PRSE_ERROR,'PFFOverRide is on, please provide PFFMaxEll')
       ENDIF
    ELSE      
       PBC%PFFMaxLay=1
       ! Look for MaxEll in periodic boundary options
       IF(.NOT.OptIntQ(Inp,PBOUNDRY,PBC%PFFMaxEll))PBC%PFFMaxEll=16
    ENDIF
    ! Parse permeability 
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
    REAL(DOUBLE),PARAMETER          :: DegToRad =  1.745329251994329576923D-2
    REAL(DOUBLE),DIMENSION(6)       :: Vec
    REAL(DOUBLE)                    :: AngAB,AngAC,AngBC,Error
    !
    PBC%Dimen=0
    PBC%TransVec=Zero
    PBC%BoxShape=Zero
    NLvec=0 
    NTvec=0   
    PBC%AutoW=.FALSE.  ! This is the inapropriatly named periodic is on flag
    PBC%TransVec=Zero
    PBC%BoxShape=Zero
    DO I=1,3
       PBC%BoxShape(I,I)=One
    ENDDO
    IF(.NOT.FindKey(BEGIN_PERIODIC,Inp))RETURN
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
             CALL LineToDbls(Line,6,Vec)
             J=0
             DO I=1,6
                IF(Vec(I)==Zero)EXIT
                J=J+1
             ENDDO
             IF(J==1)THEN
                PBC%Dimen=1
                DO I=1,3
                   IF(Vec(I)/=Zero)THEN
                      PBC%AutoW(I)=.TRUE.
                      PBC%BoxShape(I,1)=Vec(I)                               
                   ENDIF
                ENDDO
             ELSEIF(J==3) THEN                 
                PBC%Dimen=2
                IF(Vec(1)/=0.AND.Vec(2)/=0)THEN
                   PBC%AutoW(1)=.TRUE.
                   PBC%AutoW(2)=.TRUE.
                   PBC%BoxShape(1,1)=Vec(1)
                   PBC%BoxShape(1,2)=Vec(2)*COS(DegToRad*Vec(3))
                   PBC%BoxShape(2,2)=Vec(2)*SIN(DegToRad*Vec(3)) 
                ELSEIF(Vec(1)/=0 .AND.Vec(3)/=0)THEN
                   PBC%AutoW(1)=.TRUE.
                   PBC%AutoW(3)=.TRUE.
                   PBC%BoxShape(1,1)=Vec(1)
                   PBC%BoxShape(1,3)=Vec(2)*COS(DegToRad*Vec(3))
                   PBC%BoxShape(3,3)=Vec(2)*SIN(DegToRad*Vec(3))
                ELSEIF(Vec(2)/=0.AND.Vec(3)/=0)THEN
                   PBC%AutoW(2)=.TRUE.
                   PBC%AutoW(3)=.TRUE.
                   PBC%BoxShape(2,2)=Vec(1)
                   PBC%BoxShape(2,3)=Vec(2)*COS(DegToRad*Vec(3))
                   PBC%BoxShape(3,3)=Vec(2)*SIN(DegToRad*Vec(3))
                ENDIF
             ELSEIF(J==6) THEN
                PBC%Dimen=3
                PBC%AutoW=.TRUE.
                PBC%BoxShape(1,1)=Vec(1)
                PBC%BoxShape(1,2)=Vec(2)*COS(DegToRad*Vec(6))
                PBC%BoxShape(2,2)=Vec(2)*SIN(DegToRad*Vec(6))
                PBC%BoxShape(1,3)=Vec(3)*COS(DegToRad*Vec(5))
                PBC%BoxShape(2,3)=(Vec(2)*Vec(3)*COS(DegToRad*Vec(4)) &
                     -PBC%BoxShape(1,2)*PBC%BoxShape(1,3))/ PBC%BoxShape(2,2)
                PBC%BoxShape(3,3)=SQRT(Vec(3)**2-PBC%BoxShape(1,3)**2-PBC%BoxShape(2,3)**2)
                AngAB = ACOS(PBC%BoxShape(1,1)*PBC%BoxShape(1,2)/(Vec(1)*Vec(2)))/DegToRad
                AngAC = ACOS(PBC%BoxShape(1,1)*PBC%BoxShape(1,3)/(Vec(1)*Vec(3)))/DegToRad
                AngBC = PBC%BoxShape(1,2)*PBC%BoxShape(1,3)+PBC%BoxShape(2,2)*PBC%BoxShape(2,3)
                AngBC = ACOS(AngBC/(Vec(2)*Vec(3)))/DegToRad
                Error = ABS(AngAB-Vec(6))+ABS(AngAC-Vec(5))+ ABS(AngBC-Vec(4))
                IF(Error .GT. 1.D-6 ) THEN
                   CALL MondoHalt(PRSE_ERROR,'Angles Are Inccorect')
                ENDIF
             ELSE
                CALL MondoHalt(PRSE_ERROR,'Number of Magnitudes and Angles supplied is incorrect')
             ENDIF
          ENDIF
       ELSE
          EXIT
       ENDIF
    ENDDO
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

  SUBROUTINE UnitCellSetUp(PBC)
    TYPE(PBCInfo)       :: PBC
    INTEGER             :: I,J,K,NLvec,NTvec
    !CalculatetheBoxVolume
    PBC%CellVolume=One
    DO I=1,3
       IF(PBC%AutoW(I))THEN
          PBC%CellVolume=PBC%CellVolume*PBC%BoxShape(I,I)
       ENDIF
    ENDDO
    ! Calculate dipole and quadripole factors
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
    ! Find the center of the cell
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
    ! Compute the inverse box shape 
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
  END SUBROUTINE UnitCellSetUp

END MODULE ParsePeriodic
