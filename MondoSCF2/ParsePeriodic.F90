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
    ! If we are restarting or reguessing, just use values read from HDF ...
    IF(O%Guess==GUESS_EQ_RESTART.OR.O%Guess==GUESS_EQ_NUGUESS)RETURN
    ! ... otherwise, we are reading in new PBC info
    CALL OpenASCII(N%IFile,Inp)
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
       ! Do we need to supercell?
       CALL SuperCellGeom(G%Clone(I))
       !
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
          PBC%PFFMaxEll=14
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
    ! Do we need to supercell?
    CALL ParseSuperCell(PBC)
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
!        
!!$    DO I=1,3
!!$       WRITE(Out,*) (PBC%BoxShape%D(I,J),J=1,3) 
!!$       WRITE(*,*)   (PBC%BoxShape%D(I,J),J=1,3) 
!!$    ENDDO
!!$    IF(.TRUE.) STOP
!
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
       CALL CalFracCarts(G)
       IF(G%PBC%AtomW) THEN 
          CALL WrapAtoms(G)
       ENDIF
    ELSE
       G%BoxCarts%D=G%Carts%D
       CALL CalAtomCarts(G)
       IF(G%PBC%AtomW) THEN 
          CALL WrapAtoms(G)
       ENDIF
!      Convert the Velocities from Fractional to Atomic
       DO I=1,G%NAtms
          G%Velocity%D(:,I)   = FracToAtom(G,G%Velocity%D(:,I))
       ENDDO
    ENDIF
  END SUBROUTINE CalculateCoordArrays
!------------------------------------------------------------------------!
!
!------------------------------------------------------------------------!
  SUBROUTINE ParseSuperCell(PBC)
    INTEGER,DIMENSION(3) :: SC 
    TYPE(PBCInfo)        :: PBC
    INTEGER              :: I,J,K,NTot
    !
    SC(:)=0
    IF(FindKey(SUPERC,Inp))THEN
       DO K=1,10
          IF(OptKeyLocQ(Inp,SUPERC,IntToChar(K),MaxSets,NLoc,Location)) THEN
             Ntot = NLoc
             DO I=1,NLoc
                SC(Location(I)) = K
             ENDDO
          ENDIF
       ENDDO
    ELSE
       SC(:)=1
    ENDIF   
    IF(SC(1).LE.0 .OR. SC(2).LE.0 .OR. SC(3).LE.0) THEN
       CALL MondoHalt(PRSE_ERROR,'ParseSuperCell: SupreCell Number is incorrect!')
    ENDIF
    PBC%SuperCell%I(:)=SC(:)
    IF(PBC%AutoW%I(1)==0.AND.PBC%SuperCell%I(1).NE.1) &
         & CALL MondoHalt(PRSE_ERROR,'ParseSuperCell: SupreCell Number should be 1 along the a-direction!')
    IF(PBC%AutoW%I(2)==0.AND.PBC%SuperCell%I(2).NE.1) &
         & CALL MondoHalt(PRSE_ERROR,'ParseSuperCell: SupreCell Number should be 1 along the b-direction!')
    IF(PBC%AutoW%I(3)==0.AND.PBC%SuperCell%I(3).NE.1) &
         & CALL MondoHalt(PRSE_ERROR,'ParseSuperCell: SupreCell Number should be 1 along the c-direction!')
    !
    IF(SUM(PBC%SuperCell%I).NE.3) THEN
       WRITE(*,*)'ParseSuperCell: Matt be ready to blame CJ... SuperCell=(' &
            & //TRIM(IntToChar(PBC%SuperCell%I(1)))//',' &
            & //TRIM(IntToChar(PBC%SuperCell%I(2)))//',' &
            & //TRIM(IntToChar(PBC%SuperCell%I(2)))//')'
       CALL Warn('ParseSuperCell: Matt be ready to blame CJ... SuperCell=(' &
            & //TRIM(IntToChar(PBC%SuperCell%I(1)))//',' &
            & //TRIM(IntToChar(PBC%SuperCell%I(2)))//',' &
            & //TRIM(IntToChar(PBC%SuperCell%I(2)))//')' )
    ENDIF
    !
  END SUBROUTINE ParseSuperCell
!=========================================================================
!
!=========================================================================
  SUBROUTINE SuperCellGeom(G)
    TYPE(CRDS)                :: G,Gtmp
    REAL(DOUBLE)              :: SCFacX,SCFacY,SCFacZ
    INTEGER                   :: I,J,K,ISCX,ISCY,ISCZ,AT,NC
    REAL(DOUBLE),DIMENSION(3) :: RVec,IVec
!    
    IF(G%PBC%SuperCell%I(1)+G%PBC%SuperCell%I(2)+G%PBC%SuperCell%I(3)==3) RETURN
!
    ISCX = G%PBC%SuperCell%I(1)
    ISCY = G%PBC%SuperCell%I(2)
    ISCZ = G%PBC%SuperCell%I(3)
    SCFacX = DBLE(ISCX)
    SCFacY = DBLE(ISCY)    
    SCFacZ = DBLE(ISCZ)
!
    Gtmp%NAtms=G%NAtms*ISCX*ISCY*ISCZ
    CALL New(Gtmp)
    Gtmp%Nkind=G%Nkind
    Gtmp%NElec=G%NElec*ISCX*ISCY*ISCZ
    GTmp%TotCh=G%TotCh*ISCX*ISCY*ISCZ
    Gtmp%NAlph=G%NAlph*ISCX*ISCY*ISCZ
    Gtmp%NBeta=G%NBeta*ISCX*ISCY*ISCZ
    Gtmp%Multp=G%Multp
    Gtmp%InAU =G%InAU
    Gtmp%Ordrd=G%Ordrd
    Gtmp%Confg=G%Confg
    Gtmp%ETotal  =G%ETotal
    Gtmp%GradRMS =G%GradRMS
    Gtmp%GradMax =G%GradMax
    Gtmp%Unstable=G%Unstable
!    
    Gtmp%BndBox%D(1,1)= G%BndBox%D(1,1)
    Gtmp%BndBox%D(1,2)= SCFacX*G%BndBox%D(1,2)
    Gtmp%BndBox%D(2,1)= G%BndBox%D(2,1)
    Gtmp%BndBox%D(2,2)= SCFacY*G%BndBox%D(2,2)
    Gtmp%BndBox%D(3,1)= G%BndBox%D(3,1)
    Gtmp%BndBox%D(3,2)= SCFacZ*G%BndBox%D(3,2)
!
    Gtmp%PBC%Dimen     = G%PBC%Dimen
    Gtmp%PBC%PFFMaxEll = G%PBC%PFFMaxEll
    Gtmp%PBC%PFFMaxLay = G%PBC%PFFMaxLay
    Gtmp%PBC%PFFOvRide = G%PBC%PFFOvRide
    Gtmp%PBC%AtomW     = G%PBC%AtomW
    Gtmp%PBC%InVecForm = G%PBC%InVecForm 
    Gtmp%PBC%InAtomCrd = G%PBC%InAtomCrd
    Gtmp%PBC%Translate = G%PBC%Translate
!
    Gtmp%PBC%CellVolume= SCFacX*SCFacY*SCFacZ*G%PBC%CellVolume
    Gtmp%PBC%Epsilon   = G%PBC%Epsilon
    Gtmp%PBC%DipoleFAC = G%PBC%DipoleFAC/SCFacX*SCFacY*SCFacZ
    Gtmp%PBC%QupoleFAC = G%PBC%QupoleFAC/SCFacX*SCFacY*SCFacZ
!
    Gtmp%PBC%AutoW%I     = G%PBC%AutoW%I
    Gtmp%PBC%SuperCell%I = G%PBC%SuperCell%I
!
    Gtmp%PBC%CellCenter%D(1)=SCFacX*G%PBC%CellCenter%D(1)
    Gtmp%PBC%CellCenter%D(2)=SCFacY*G%PBC%CellCenter%D(2)
    Gtmp%PBC%CellCenter%D(3)=SCFacZ*G%PBC%CellCenter%D(3)
    Gtmp%PBC%TransVec%D(1)  =SCFacX*G%PBC%TransVec%D(1)
    Gtmp%PBC%TransVec%D(2)  =SCFacY*G%PBC%TransVec%D(2)
    Gtmp%PBC%TransVec%D(3)  =SCFacZ*G%PBC%TransVec%D(3)
!
    Gtmp%PBC%BoxShape%D(1:3,1)=SCFacX*G%PBC%BoxShape%D(1:3,1)
    Gtmp%PBC%BoxShape%D(1:3,2)=SCFacY*G%PBC%BoxShape%D(1:3,2)
    Gtmp%PBC%BoxShape%D(1:3,3)=SCFacZ*G%PBC%BoxShape%D(1:3,3)
    Gtmp%PBC%LatFrc%D(1:3,1)  =SCFacX*G%PBC%LatFrc%D(1:3,1)
    Gtmp%PBC%LatFrc%D(1:3,2)  =SCFacY*G%PBC%LatFrc%D(1:3,2)
    Gtmp%PBC%LatFrc%D(1:3,3)  =SCFacZ*G%PBC%LatFrc%D(1:3,3)
!
    Gtmp%PBC%InvBoxSh%D = InverseMatrix(Gtmp%PBC%BoxShape%D)
!
    NC = 0
    DO I=1,ISCX
       DO J=1,ISCY
          DO K=1,ISCZ
             RVec(1) = G%PBC%BoxShape%D(1,1)*DBLE(I-1)+G%PBC%BoxShape%D(1,2)*DBLE(J-1)+G%PBC%BoxShape%D(1,3)*DBLE(K-1)
             RVec(2) = G%PBC%BoxShape%D(2,1)*DBLE(I-1)+G%PBC%BoxShape%D(2,2)*DBLE(J-1)+G%PBC%BoxShape%D(2,3)*DBLE(K-1)
             RVec(3) = G%PBC%BoxShape%D(3,1)*DBLE(I-1)+G%PBC%BoxShape%D(3,2)*DBLE(J-1)+G%PBC%BoxShape%D(3,3)*DBLE(K-1)
             IVec(1) = DBLE(I-1)
             IVec(2) = DBLE(J-1)
             IVec(3) = DBLE(K-1)
             DO AT=1,G%NAtms
                NC = NC+1
                Gtmp%AtNum%D(NC)      = G%AtNum%D(AT)
                Gtmp%AtTyp%I(NC)      = G%AtTyp%I(AT)
                Gtmp%AtNam%C(NC)      = G%AtNam%C(AT)
                Gtmp%AtMMTyp%C(NC)    = G%AtMMTyp%C(AT)
                Gtmp%AtMss%D(NC)      = G%AtMss%D(AT)
                Gtmp%CConstrain%I(NC) = G%CConstrain%I(AT)
!
                Gtmp%Velocity%D(1:3,NC) = G%Velocity%D(1:3,AT)  
                Gtmp%Gradients%D(1:3,NC)= G%Gradients%D(1:3,AT)
!
                Gtmp%Carts%D(1:3,NC)    = G%Carts%D(1:3,AT)+RVec(1:3) 
                Gtmp%BoxCarts%D(1:3,NC) = G%BoxCarts%D(1:3,AT)+IVec(1:3) 
             ENDDO
          ENDDO
       ENDDO
    ENDDO
!

    WRITE(*,*) 'SuperCell Option is On:'
    WRITE(*,*) 'No. Atoms in PrimCell  = ',G%NAtms
    WRITE(*,*) 'No. Atoms in SuperCell = ',Gtmp%NAtms
!
    CALL Delete(G)
    G%NAtms=Gtmp%NAtms
    CALL New(G)

    G%Nkind=Gtmp%Nkind
    G%NElec=Gtmp%NElec
    G%TotCh=Gtmp%TotCh
    G%NAlph=Gtmp%NAlph
    G%NBeta=Gtmp%NBeta
    G%Multp=Gtmp%Multp
    G%InAU =Gtmp%InAU
    G%Ordrd=Gtmp%Ordrd
    G%Confg=Gtmp%Confg
    G%ETotal  =Gtmp%ETotal
    G%GradRMS =Gtmp%GradRMS
    G%GradMax =Gtmp%GradMax
    G%Unstable=Gtmp%Unstable
    G%BndBox%D=Gtmp%BndBox%D
!
    G%PBC%Dimen     = Gtmp%PBC%Dimen
    G%PBC%PFFMaxEll = Gtmp%PBC%PFFMaxEll
    G%PBC%PFFMaxLay = Gtmp%PBC%PFFMaxLay
    G%PBC%PFFOvRide = Gtmp%PBC%PFFOvRide
    G%PBC%AtomW     = Gtmp%PBC%AtomW
    G%PBC%InVecForm = Gtmp%PBC%InVecForm 
    G%PBC%InAtomCrd = Gtmp%PBC%InAtomCrd
    G%PBC%Translate = Gtmp%PBC%Translate
!
    G%PBC%CellVolume= Gtmp%PBC%CellVolume
    G%PBC%Epsilon   = Gtmp%PBC%Epsilon
    G%PBC%DipoleFAC = Gtmp%PBC%DipoleFAC
    G%PBC%QupoleFAC = Gtmp%PBC%QupoleFAC
!
    G%PBC%AutoW%I     = Gtmp%PBC%AutoW%I
    G%PBC%SuperCell%I = Gtmp%PBC%SuperCell%I
!
    G%PBC%CellCenter%D=Gtmp%PBC%CellCenter%D
    G%PBC%TransVec%D  =Gtmp%PBC%TransVec%D
!
    G%PBC%BoxShape%D = Gtmp%PBC%BoxShape%D
    G%PBC%InvBoxSh%D = Gtmp%PBC%InvBoxSh%D
    G%PBC%LatFrc%D   = Gtmp%PBC%LatFrc%D
!
    G%AtNum%D      = Gtmp%AtNum%D
    G%AtTyp%I      = Gtmp%AtTyp%I
    G%AtNam%C      = Gtmp%AtNam%C
    G%AtMMTyp%C    = Gtmp%AtMMTyp%C
    G%AtMss%D      = Gtmp%AtMss%D
    G%CConstrain%I = Gtmp%CConstrain%I
!
    G%Velocity%D   = Gtmp%Velocity%D 
    G%Gradients%D  = Gtmp%Gradients%D
!
    G%Carts%D      = Gtmp%Carts%D
    G%BoxCarts%D   = Gtmp%BoxCarts%D
!
  END SUBROUTINE SuperCellGeom
!-----------------------------------------------------------------------!
END MODULE ParsePeriodic
