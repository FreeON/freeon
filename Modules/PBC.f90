!
!    ROUTINES FOR THE IMPLIMENTATION OF PERIODIC BOUNDARY CONDITIONS
!
!    Matt Challacombe and C.J. Tymczak
!
MODULE PBC
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE GlobalObjects
  USE PrettyPrint
  USE Thresholding
  USE LinAlg
  USE Order
  USE MemMan
  IMPLICIT NONE

  REAL(DOUBLE) :: PwMAX

CONTAINS
  !--------------------------------------------------------------------
  SUBROUTINE MkGeomPeriodic(GM,BS_O,Dist_O,TwoE_O,Rescale_O,SoftReset_O)
    TYPE(CRDS)                     :: GM
    TYPE(BSET),OPTIONAL            :: BS_O
    REAL(DOUBLE),OPTIONAL          :: Dist_O,TwoE_O
    LOGICAL, OPTIONAL              :: Rescale_O,SoftReset_O
    REAL(DOUBLE)                   :: X,Y,Z,MinExpt
    REAL(DOUBLE),DIMENSION(3,3)    :: OldInvBox
    REAL(DOUBLE),DIMENSION(3)      :: Cnt
    INTEGER                        :: I,J,K,N,CF,PF,NK
    LOGICAL                        :: Rescale,SoftReset
    !====================================================================================
    !  VOLUME AND ASSOCIATED CONSTANTS RELATED TO TIN-FOIL BOUNDARY CONDITIONS
    !====================================================================================
    ! Recompute the PBC cell volume and cell center (may be different than the bounding box center)
    GM%PBC%CellVolume=CellVolume(GM%PBC%BoxShape%D,GM%PBC%AutoW%I)
    ! Calculate the Dipole and Quadripole Factors
    IF(GM%PBC%Dimen < 2) THEN
       GM%PBC%DipoleFAC = Zero
       GM%PBC%QupoleFAC = Zero
    ELSEIF(GM%PBC%Dimen ==2) THEN
       GM%PBC%DipoleFAC = (Four*Pi/GM%PBC%CellVolume)*(One/(GM%PBC%Epsilon+One))
       GM%PBC%QupoleFAC =  Zero
       IF(ABS(GM%PBC%DipoleFAC) .LT. 1.D-14) GM%PBC%DipoleFAC = Zero
    ELSEIF(GM%PBC%Dimen ==3) THEN
       GM%PBC%DipoleFAC = -(Four*Pi/GM%PBC%CellVolume)*(One/Three - One/(Two*GM%PBC%Epsilon+One))
       GM%PBC%QupoleFAC =  ( Two*Pi/GM%PBC%CellVolume)*(One/Three - One/(Two*GM%PBC%Epsilon+One))
       IF(ABS(GM%PBC%DipoleFAC) .LT. 1.D-14) GM%PBC%DipoleFAC = Zero
       IF(ABS(GM%PBC%QupoleFAC) .LT. 1.D-14) GM%PBC%QupoleFAC = Zero
    ENDIF
    !====================================================================================
    !  ATOMS
    !====================================================================================
    IF(GM%PBC%Dimen>0)THEN
       !====================================================================================
       !  IF RESCALE==TRUE, THEN IT IS GENERALLY ASSUMED THAT (1) THE COORDINATES
       !  DEPEND ON THE LATTICE AND (2) THE LATTICE MAY HAVE CHANGED WHEN THIS ROUTINE IS CALLED,
       !  NESSITATING THAT THE ATOMIC COORDINATES BE RECOMPUTED WITH THE NEW, PUTATIVE LATTICE.
       !  AND EXAMPLE OF WHERE THIS IS IMPORTANT IS IN THE COMPUTATION OF LATTICE FINITE DIFFERENCES.
       !====================================================================================
       IF(PRESENT(Rescale_O))THEN
          Rescale=Rescale_O
       ELSE
          Rescale=.FALSE.
       ENDIF
       IF(Rescale)THEN
          !
          ! Rescale: It is assumed that we shold first compute fractional atomic coordinates
          ! cooresponding to (1) an OLD INVERSE and (2) NEW ATOMIC POSITIONS.  Then compute a
          ! NEW INVERSE cooresponding to a potentially NEW LATTICE, and then recompute atomic
          ! positions.  This is useful for, eg finite differences of lattice coordinates.
          ! However, this option will severely hose the geometry optimizer.
          !
          ! For finite difference schemes and analytic derivatives, this cooresponds to variation
          ! of the coordinates as follows:  A_new = Box_new.[InverseBox_old.A].  Thus, the
          ! InverseBox should be assumed constant, while Box_new and A hold all the active variables,
          ! FOR THE ATOMIC POSITIONS (ONLY).  For other items, such as the wrapped distribution center,
          ! we must consider that the wrapping function is entirely variational with respect to the
          ! lattice.  That is,  WRAPED P = Box_new.MOD[InverseBox_new.P].
          !
          DO I=1,GM%NAtms
             ! Overwrite the atomic positions with fractionals based on OLD INVERSE AND NEW LATTICE/COORDINATES
             GM%Carts%D(:,I)=AtomToFrac(GM,GM%Carts%D(:,I))
             ! wrap the fractionals
             !          WRITE(*,*)I,'FRACTED CARTS = ',GM%Carts%D(:,I)
             CALL FracCyclic(GM,GM%Carts%D(:,I))
             ! Boxcarts are the fractional coordinates
             GM%BoxCarts%D(:,I)=GM%Carts%D(:,I)
          ENDDO
          ! First save the inverse lattice vectors for later (maybe)...
          OldInvBox=GM%PBC%InvBoxSh%D
          ! ... then generate the possibly new inverse lattice coordinates
          GM%PBC%InvBoxSh%D=InverseBoxShape(GM%PBC%BoxShape%D,GM%PBC%Dimen)
          ! Now, regenerate the possibly dialated or contracted atomic positions from the fractionals
          DO I=1,GM%NAtms
             GM%Carts%D(:,I)=FracToAtom(GM,GM%Carts%D(:,I))
          ENDDO
       ELSE
          ! No rescale: Here, we just recompute the inverse lattice vectors
          ! and wrap atoms back into the unit cell.
          OldInvBox=GM%PBC%InvBoxSh%D
          GM%PBC%InvBoxSh%D=InverseBoxShape(GM%PBC%BoxShape%D,GM%PBC%Dimen)
          DO I=1,GM%NAtms
             GM%Carts%D(:,I)=AtomToFrac(GM,GM%Carts%D(:,I))
             CALL FracCyclic(GM,GM%Carts%D(:,I))
             GM%BoxCarts%D(:,I)=GM%Carts%D(:,I)
             GM%Carts%D(:,I)=FracToAtom(GM,GM%Carts%D(:,I))
          ENDDO
       ENDIF
    ENDIF
    !====================================================================================
    ! CELL CENTER
    !====================================================================================
    ! Compute the cell center with the modified positions (Zero for no PBCs)
    GM%PBC%CellCenter%D=CellCenter(GM%PBC)
    ! Recompute the bounding box and ITS cell center (May be different from the PBC cell center)
    GM%BndBox%D(1:3,1)=GM%Carts%D(1:3,1)
    GM%BndBox%D(1:3,2)=GM%Carts%D(1:3,1)
    DO J=2,GM%Natms
       GM%BndBox%D(1:3,1)=MIN(GM%BndBox%D(1:3,1),GM%Carts%D(1:3,J))
       GM%BndBox%D(1:3,2)=MAX(GM%BndBox%D(1:3,2),GM%Carts%D(1:3,J))
    ENDDO
    !====================================================================================
    ! CELL SETS
    !====================================================================================
    ! To here we have set (or reset) the unit cell.  Now we set (or reset) the different
    ! sets of periodic images (cellsets) for the sum over basis functions (outer==overlap) and
    ! the inner sum over Coulomb boxes (inner==penetration).
    !====================================================================================
    !  IF SoftReset==TRUE, THEN WE WILL NOT RECOMPUTE THE NUMBER OF CELL SETS.
    !  WITH SoftReset==TRUE, ONLY THE NUMERICAL VALUES OF THE CELL SETS WILL CHANGE.
    !  This option is not fully implemented, as it currently needs the optional basis
    !  set info, so the hard/soft is controlled by that option.  This should be ok,
    !  except in situations where the lattice changes a lot.
    !====================================================================================
    IF(PRESENT(SoftReset_O))THEN
       SoftReset=SoftReset_O
    ELSE
       SoftReset=.FALSE.
    ENDIF
!
    IF(PRESENT(BS_O))THEN
!       CALL Warn(' Hard reset of the PBC cell-sets')
       ! Logic check for optionals
       IF(.NOT.PRESENT(Dist_O))CALL MondoHalt(DRIV_ERROR,' Missing Dist_O threshold option in MkGeomPeriodic ')
       IF(.NOT.PRESENT(TwoE_O))CALL MondoHalt(DRIV_ERROR,' Missing TwoE_O threshold option in MkGeomPeriodic ')
       IF(AllocQ(GM%OvCells%Alloc))CALL Delete(GM%OvCells)
       IF(AllocQ(GM%InCells%Alloc))CALL Delete(GM%InCells)
       ! Find the min exponents for this basis set
       MinExpt=1D25
       DO NK=1,BS_O%NKind
          DO CF=1,BS_O%NCFnc%I(NK)
             DO PF=1,BS_O%NPFnc%I(CF,NK)
                MinExpt=MIN(MinExpt,BS_O%Expnt%D(PF,CF,NK))
             ENDDO
          ENDDO
       ENDDO
       ! Generally, the minimum Gaussian exponent (MinExpt)in this basis set determines the periodic
       ! images to be summed over in the Coulomb sums and in the periodic sum over basis functions:
       CALL SetCellSets(GM%OvCells,GM%PBC%AutoW%I,GM%PBC%BoxShape%D,'FreeON','Overlap',Dist_O,MinExpt,GM%NElec)
       !       CALL SetCellSets(GM%OvCells,GM%PBC%AutoW%I,GM%PBC%BoxShape%D,'HardCell',HardCell_O=1)
       ! If WellSeperated (WS) criteria is not zero, then we use a strict FMM
       ! definition to determine the penetration cell sets.  NOTE: use of WS
       ! criteria can lead to errors for small cells, and also for non-cubic
       ! cells, but is usefull in debugging Coulomb sums.
       IF(GM%PBC%PFFWelSep.NE.0)THEN
          CALL SetCellSets(GM%InCells,GM%PBC%AutoW%I,GM%PBC%BoxShape%D,'FreeON','HardCell',HardCell_O=GM%PBC%PFFWelSep)
          CALL Warn('Fixed FMM style WS = '//TRIM(IntToChar(GM%PBC%PFFWelSep))// &
            ', may lead to errors for small or non-cubic cells')
       ELSE
          !          CALL SetCellSets(GM%InCells,GM%PBC%AutoW%I,GM%PBC%BoxShape%D,'HardCell',HardCell_O=2)
          CALL SetCellSets(GM%InCells,GM%PBC%AutoW%I,GM%PBC%BoxShape%D,'FreeON','Penetration', &
               CellSetsThreshold_O=TwoE_O,MinExpt_O=MinExpt,NElect_O=GM%NElec,MaxEll_O=GM%PBC%PFFMaxEll)
       ENDIF
    ELSE
       CALL Warn('Soft reset of the PBC cell-sets')
       ! Check to see if the cell set is already allocated
       IF(.NOT.AllocQ(GM%OvCells%Alloc))CALL MondoHalt(DRIV_ERROR,'Resetting of non allocated cell sets')
       IF(.NOT.AllocQ(GM%InCells%Alloc))CALL MondoHalt(DRIV_ERROR,'Resetting of non allocated cell sets')
       ! OK, basis set is already allocated, stick with existing number
       ! of cells, but reset the existing displacements etc
       DO N=1,GM%InCells%NCells
          X=GM%InCells%CellCarts%D(1,N)
          Y=GM%InCells%CellCarts%D(2,N)
          Z=GM%InCells%CellCarts%D(3,N)
          I=X*OldInvBox(1,1)+Y*OldInvBox(1,2)+Z*OldInvBox(1,3)
          J=X*OldInvBox(2,1)+Y*OldInvBox(2,2)+Z*OldInvBox(2,3)
          K=X*OldInvBox(3,1)+Y*OldInvBox(3,2)+Z*OldInvBox(3,3)
          X=I*GM%PBC%BoxShape%D(1,1)+J*GM%PBC%BoxShape%D(1,2)+K*GM%PBC%BoxShape%D(1,3)
          Y=I*GM%PBC%BoxShape%D(2,1)+J*GM%PBC%BoxShape%D(2,2)+K*GM%PBC%BoxShape%D(2,3)
          Z=I*GM%PBC%BoxShape%D(3,1)+J*GM%PBC%BoxShape%D(3,2)+K*GM%PBC%BoxShape%D(3,3)
          GM%InCells%CellCarts%D(:,N)=(/X,Y,Z/)
       ENDDO
       DO N=1,GM%OvCells%NCells
          X=GM%OvCells%CellCarts%D(1,N)
          Y=GM%OvCells%CellCarts%D(2,N)
          Z=GM%OvCells%CellCarts%D(3,N)
          I=X*OldInvBox(1,1)+Y*OldInvBox(1,2)+Z*OldInvBox(1,3)
          J=X*OldInvBox(2,1)+Y*OldInvBox(2,2)+Z*OldInvBox(2,3)
          K=X*OldInvBox(3,1)+Y*OldInvBox(3,2)+Z*OldInvBox(3,3)
          X=I*GM%PBC%BoxShape%D(1,1)+J*GM%PBC%BoxShape%D(1,2)+K*GM%PBC%BoxShape%D(1,3)
          Y=I*GM%PBC%BoxShape%D(2,1)+J*GM%PBC%BoxShape%D(2,2)+K*GM%PBC%BoxShape%D(2,3)
          Z=I*GM%PBC%BoxShape%D(3,1)+J*GM%PBC%BoxShape%D(3,2)+K*GM%PBC%BoxShape%D(3,3)
          GM%OvCells%CellCarts%D(:,N)=(/X,Y,Z/)
       ENDDO
    ENDIF
!!$

!    CALL Print_CRDS(GM,PrintGeom_O='XYZ')
!    CALL Print_CRDS(GM,Unit_O=6,PrintGeom_O='XSF')!,CrdInAng_O=.FALSE.)
    !
  END SUBROUTINE MkGeomPeriodic
  !--------------------------------------------------------------------------
  ! Here is the new Cell Sets routine. This should become the only routine
  ! used for setting the CellSet data structure.
  !--------------------------------------------------------------------------
  SUBROUTINE SetCellSets(CS,AW,MAT,Prog,Option, &
       CellSetsThreshold_O,MinExpt_O,NElect_O,MaxEll_O,Rmin_O,HardCell_O)
    TYPE(CellSet)                        :: CS
    INTEGER,DIMENSION(3)                 :: AW
    REAL(DOUBLE),DIMENSION(3,3)          :: MAT
    REAL(DOUBLE),DIMENSION(3)            :: UpperUC,LowerUC
    REAL(DOUBLE)                         :: CellSetsThreshold,MinExpt,HardCell,T,NElect,R2Min, &
         MaxT,TDa,TDb,TDc
    REAL(DOUBLE),OPTIONAL                :: Rmin_O,CellSetsThreshold_O,MinExpt_O
    INTEGER,OPTIONAL                     :: NElect_O,HardCell_O,MaxEll_O
    CHARACTER(LEN=*)                     :: Prog
    CHARACTER(LEN=*)                     :: Option
    CHARACTER(LEN=132)                   :: Mssg
    !
    INTEGER                              :: I,J,K,L,M,II,JJ,KK,MaxEll,CheckNCell
    INTEGER                              :: IXM,IYM,IZM,PMx,PMy,PMz,Ipm,Jpm,Kpm,NCELL
    REAL(DOUBLE)                         :: X,Y,Z,PQx,PQy,PQz,PQ,PQ2,PQMinSep,PQ2MinSep, &
         R,R2,R2Max,MaxTranslation
    LOGICAL                              :: InCell,MAC,PAC
    !
    HardCell=20
    !
    IF(Option=="Radial")THEN
       IF(PRESENT(RMin_O))THEN
          R2Min=RMin_O**2
       ELSE
          CALL Halt(' Optional argument RMin_O not present with Option=Radial in SetCellSets ')
       ENDIF
    ELSEIF(Option=='HardCell')THEN
       IF(PRESENT(HardCell_O))THEN
          HardCell=HardCell_O
       ELSE
          CALL Halt(' Optional argument HardCell_O not present with Option=HardCell in SetCellSets ')
       ENDIF
    ELSEIF(Option=='Penetration'.OR.Option=='Overlap')THEN
       IF(PRESENT(MinExpt_O))THEN
          MinExpt=MinExpt_O
       ELSE
          CALL Halt(' Optional argument MinExpt_O not present in SetCellSets ')
       ENDIF
       IF(PRESENT(CellSetsThreshold_O))THEN
          CellSetsThreshold=CellSetsThreshold_O
       ELSE
          CALL Halt(' Optional argument CellSetsThreshold_O not present in SetCellSets ')
       ENDIF
       IF(Option=='Penetration')THEN
          IF(PRESENT(MaxEll_O))THEN
             MaxEll=MaxEll_O
          ELSE
             CALL Halt(' Optional argument MaxEll_O not present in SetCellSets ')
          ENDIF
       ENDIF
       IF(PRESENT(NElect_O))THEN
          NElect=DBLE(NElect_O)
       ELSE
          CALL Halt(' Optional argument NElect_O not present in SetCellSets ')
       ENDIF
    ENDIF
    !
    IXM=0; PMX=0
    IYM=0; PMY=0
    IZM=0; PMZ=0
    IF(AW(1)==1)THEN
       IXM=HardCell
       PMX=1
    ENDIF
    IF(AW(2)==1)THEN
       IYM=HardCell
       PMY=1
    ENDIF
    IF(AW(3)==1)THEN
       IZM=HardCell
       PMZ=1
    ENDIF
    !
    IF(Option=='Penetration')THEN
       MaxTranslation=Zero
       DO II=0,1
          DO JJ=0,1
             DO KK=0,1
                MaxTranslation=MAX(MaxTranslation,VABS((/II*MAT(:,1)+JJ*MAT(:,2)+KK*MAT(:,3)/)))
             ENDDO
          ENDDO
       ENDDO
    ENDIF

    IF(Option=='HardCell')THEN
       !----------------------------------------------------------------------------------------
       ! HardCell sets a fixed cell set with rectalinear integer shape based on input parameters
       ! >>This option will give horrible results for strongly non-orthorombic cells<<
       !-----------------------------------------------------------------------------------------
       NCell=(2*IXM+1)*(2*IYM+1)*(2*IZM+1)
       CALL New(CS,NCell)
       NCell = 0
       DO I=-IXM,IXM
          DO J=-IYM,IYM
             DO K=-IZM,IZM
                X=I*MAT(1,1)+J*MAT(1,2)+K*MAT(1,3)
                Y=I*MAT(2,1)+J*MAT(2,2)+K*MAT(2,3)
                Z=I*MAT(3,1)+J*MAT(3,2)+K*MAT(3,3)
                NCell = NCell+1
                CS%CellCarts%D(:,NCell)=(/X,Y,Z/)
             ENDDO
          ENDDO
       ENDDO
    ELSEIF(Option=="Radial")THEN
       !------------------------------------------------------------------------------------
       ! Radial sets a fixed cell set with spherical shape based on input parameters
       !-----------------------------------------------------------------------------------
       NCell=0
       DO I=-IXM,IXM
          DO J=-IYM,IYM
             DO K=-IZM,IZM
                X=I*MAT(1,1)+J*MAT(1,2)+K*MAT(1,3)
                Y=I*MAT(2,1)+J*MAT(2,2)+K*MAT(2,3)
                Z=I*MAT(3,1)+J*MAT(3,2)+K*MAT(3,3)
                R2=X*X+Y*Y+Z*Z
                IF(R2<=R2Min)THEN
                   NCell = NCell+1
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       CALL New_CellSet(CS,NCell)
       CS%NCells=NCell
       NCell = 0
       DO I=-IXM,IXM
          DO J=-IYM,IYM
             DO K=-IZM,IZM
                X=I*MAT(1,1)+J*MAT(1,2)+K*MAT(1,3)
                Y=I*MAT(2,1)+J*MAT(2,2)+K*MAT(2,3)
                Z=I*MAT(3,1)+J*MAT(3,2)+K*MAT(3,3)
                R2=X*X+Y*Y+Z*Z
                IF(R2<=R2Min)THEN
                   NCell = NCell+1
                   CS%CellCarts%D(:,NCell)=(/X,Y,Z/)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ELSEIF(Option=='Penetration'.OR.Option=='Overlap')THEN
       !
       NCell=0
       R2Min=0
       DO I=-IXM,IXM
          DO J=-IYM,IYM
             DO K=-IZM,IZM
                !
                PQx=DBLE(I)*MAT(1,1)+DBLE(J)*MAT(1,2)+DBLE(K)*MAT(1,3)
                PQy=DBLE(I)*MAT(2,1)+DBLE(J)*MAT(2,2)+DBLE(K)*MAT(2,3)
                PQz=DBLE(I)*MAT(3,1)+DBLE(J)*MAT(3,2)+DBLE(K)*MAT(3,3)                !
                !-----------------------------------------------------------------------------------
                ! PQ is the cell-cell distance measured center to center, which is used
                ! in the multipole error criteria, while PQMinSep is the min cell-cell distance measured
                ! between closest corners, used for the penetration error criteria and the overlap
                !-----------------------------------------------------------------------------------
                PQ2=PQx**2+PQy**2+PQz**2
                PQ2MinSep=PQ2
                DO Ipm=-1,1
                   DO Jpm=-1,1
                      DO Kpm=-1,1
                         PQ2MinSep=MIN(PQ2MinSep, &
                               (PQx-(DBLE(Ipm)*MAT(1,1)+DBLE(Jpm)*MAT(1,2)+DBLE(Kpm)*MAT(1,3)))**2 &
                              +(PQy-(DBLE(Ipm)*MAT(2,1)+DBLE(Jpm)*MAT(2,2)+DBLE(Kpm)*MAT(2,3)))**2 &
                              +(PQz-(DBLE(Ipm)*MAT(3,1)+DBLE(Jpm)*MAT(3,2)+DBLE(Kpm)*MAT(3,3)))**2)
                      ENDDO
                   ENDDO
                ENDDO
                PQ=SQRT(PQ2)
                PQ2MinSep=PQ2MinSep+1D-10
                PQMinSep=SQRT(PQ2MinSep)

                InCell=.FALSE.
                IF(I==0.AND.J==0.AND.K==0)THEN
                   InCell=.TRUE.
                ELSEIF(Option=='Penetration')THEN
                   IF(PQ>Two*MaxTranslation)THEN

                      MAC=(NElect/(PQ-MaxTranslation))*(MaxTranslation/PQ)**(MaxEll) > CellSetsThreshold *1D5

                   ELSE
                      MAC=.TRUE.
                   ENDIF
                   !-----------------------------------------------------------------------------------
                   ! The penetration error, due to overlaping charge distributions must also be
                   ! considered.  Here, a simple s-s Coulomp error is used, with furthest corner
                   ! to closest corner distance
                   !-----------------------------------------------------------------------------------
                   PAC=NElect*EXP(-MinExpt*PQ2MinSep)/PQ2MinSep>CellSetsThreshold
                   !-----------------------------------------------------------------------------------
                   ! A false MAC or PAC adds this cell.
                   ! Note, this is the converse of MAC and PAC used in QCTC/TreeWalk.
                   !-----------------------------------------------------------------------------------
                   InCell=MAC.OR.PAC
                   IF(InCell)R2Min=MAX(PQ2,R2Min)

               ELSEIF(Option=='Overlap')THEN
                   !-----------------------------------------------------------------------------------
                   ! For the Overlap, PQ is the min separation that can occur between atom pairs A and B.
                   ! Note that due to periodic translation, this is not quite as simple as it
                   ! seems, since what we really need to wory about is the min value that A-(B+R)
                   ! can take.  Well, we know that the most A-B can be is +/- R.  Using PQ2MinSep accounts
                   ! for this.
                   !
                   ! And here is simply the overlap criterial, with Min(XiAB)=MinExpt/2
                   !-----------------------------------------------------------------------------------
                   T=Half*MinExpt*PQ2MinSep
                   InCell=NElect*EXP(-T)>CellSetsThreshold
                   IF(InCell)R2Min=MAX(PQ2MinSep,R2Min)
                ELSE
                   CALL Halt('A Unknown option to SetCellSets, OPTION=<'//Option//'>')
                ENDIF
                !-------------------------------------------------------------------------------------
                ! Here we are just counting how many cells.  See below for coments re double counting
                !-------------------------------------------------------------------------------------
                IF(InCell)THEN
                   NCell=NCell+1
                ENDIF
                !
             ENDDO
          ENDDO
       ENDDO
       !-------------------------------------------------------------------------------------
       ! Now allocate the CellSet data structure that was just counted out.
       !-------------------------------------------------------------------------------------
       CALL New_CellSet(CS,NCell)
       CS%NCells=NCell
       CheckNCell=NCell
       !
       NCell = 0
       DO I=-IXM,IXM
          DO J=-IYM,IYM
             DO K=-IZM,IZM
                !
                PQx=DBLE(I)*MAT(1,1)+DBLE(J)*MAT(1,2)+DBLE(K)*MAT(1,3)
                PQy=DBLE(I)*MAT(2,1)+DBLE(J)*MAT(2,2)+DBLE(K)*MAT(2,3)
                PQz=DBLE(I)*MAT(3,1)+DBLE(J)*MAT(3,2)+DBLE(K)*MAT(3,3)                !
                !----------------------------------------------------------------------
                ! Here we have same logic as above, but instead of just counting
                ! we now fill in the cell sets data structure. Any changes to the logic
                ! must be propigated both above and below, or there will be a mis-counting.
                ! Certainly, this double counting could be avoided with a more
                ! advanced data structure, like a linked list!  Cost for this double count
                ! should be quite small tho...
                !----------------------------------------------------------------------
                PQ2=PQx**2+PQy**2+PQz**2
                PQ2MinSep=PQ2
                DO Ipm=-1,1
                   DO Jpm=-1,1
                      DO Kpm=-1,1
                         PQ2MinSep=MIN(PQ2MinSep, &
                               (PQx-(DBLE(Ipm)*MAT(1,1)+DBLE(Jpm)*MAT(1,2)+DBLE(Kpm)*MAT(1,3)))**2 &
                              +(PQy-(DBLE(Ipm)*MAT(2,1)+DBLE(Jpm)*MAT(2,2)+DBLE(Kpm)*MAT(2,3)))**2 &
                              +(PQz-(DBLE(Ipm)*MAT(3,1)+DBLE(Jpm)*MAT(3,2)+DBLE(Kpm)*MAT(3,3)))**2)
                      ENDDO
                   ENDDO
                ENDDO
                PQ=SQRT(PQ2)
                PQ2MinSep=PQ2MinSep+1D-10
                PQMinSep=SQRT(PQ2MinSep)

                InCell=.FALSE.
                IF(I==0.AND.J==0.AND.K==0)THEN
                   InCell=.TRUE.
                ELSEIF(Option=='Penetration')THEN
                   IF(PQ>Two*MaxTranslation)THEN

                      MAC=(NElect/(PQ-MaxTranslation))*(MaxTranslation/PQ)**(MaxEll) > CellSetsThreshold *1D5

                   ELSE
                      MAC=.TRUE.
                   ENDIF
                   !-----------------------------------------------------------------------------------
                   ! The penetration error, due to overlaping charge distributions must also be
                   ! considered.  Here, a simple s-s Coulomp error is used, with furthest corner
                   ! to closest corner distance
                   !-----------------------------------------------------------------------------------
                   PAC=NElect*EXP(-MinExpt*PQ2MinSep)/PQ2MinSep>CellSetsThreshold
                   !-----------------------------------------------------------------------------------
                   ! A false MAC or PAC adds this cell.
                   ! Note, this is the converse of MAC and PAC used in QCTC/TreeWalk.
                   !-----------------------------------------------------------------------------------
                   InCell=MAC.OR.PAC
                   IF(InCell)R2Min=MAX(PQ2,R2Min)
               ELSEIF(Option=='Overlap')THEN
                   !-----------------------------------------------------------------------------------
                   ! For the Overlap, PQ is the min separation that can occur between atom pairs A and B.
                   ! Note that due to periodic translation, this is not quite as simple as it
                   ! seems, since what we really need to wory about is the min value that A-(B+R)
                   ! can take.  Well, we know that the most A-B can be is +/- R.  Using PQ2MinSep accounts
                   ! for this.
                   !
                   ! And here is simply the overlap criterial, with Min(XiAB)=MinExpt/2
                   !-----------------------------------------------------------------------------------
                   T=Half*MinExpt*PQ2MinSep
                   InCell=NElect*EXP(-T)>CellSetsThreshold
                ELSE
                   CALL Halt('B Unknown option to SetCellSets, OPTION=<'//Option//'>')
                ENDIF
                !----------------------------------------------------------------------
                ! Here is the filling in of the cell sets array
                !----------------------------------------------------------------------
                IF(InCell)THEN
                   NCell = NCell+1
                   X=DBLE(I)*MAT(1,1)+DBLE(J)*MAT(1,2)+DBLE(K)*MAT(1,3)
                   Y=DBLE(I)*MAT(2,1)+DBLE(J)*MAT(2,2)+DBLE(K)*MAT(2,3)
                   Z=DBLE(I)*MAT(3,1)+DBLE(J)*MAT(3,2)+DBLE(K)*MAT(3,3)
                   CS%CellCarts%D(:,NCell)=(/X,Y,Z/)
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       IF(NCell.NE.CheckNCell)THEN
          WRITE(*,*)' Option = ',Option
          WRITE(*,*)' NCell counting= ',CheckNCell
          WRITE(*,*)' NCell filling = ',NCell
          CALL Halt(' Counting and filling NCell values off in SetCellSets ')
       ENDIF

    ELSE
       CALL Halt('C Unknown option to SetCellSets, OPTION=<'//Option//'>')
    ENDIF
    !
    CS%Radius=SQRT(R2Min)

    !-------------------------------------------------------------------------------
    ! Sort so that we start furthest away, adding in the smallest contributions
    ! of whatever first (summing from smallest to largest)
    !-------------------------------------------------------------------------------
    CALL Sort_CellSet(CS,Order_O=2)

    CALL MondoLog(DEBUG_MEDIUM,Prog, &
         'Cells = <'//TRIM(IntToChar(CS%NCells))//'>, Radius = <'//TRIM(DblToShrtChar(CS%Radius))//'>' &
         ,TRIM(Option))

    !     CALL PPrint_CellSet(CS,Option,Unit_O=6)
  END SUBROUTINE SetCellSets
  !--------------------------------------------------------------------------
  ! Sortharthe Cells From Large R to Small R
  !--------------------------------------------------------------------------
  SUBROUTINE Sort_CellSet(CS,Order_O)
    TYPE(CellSet)                        :: CS
    TYPE(INT_VECT)                       :: IPnt
    TYPE(DBL_VECT)                       :: Vec
    TYPE(DBL_RNK2)                       :: CCarts
    INTEGER                              :: NC
    INTEGER, OPTIONAL                    :: Order_O
    INTEGER                              :: Order
    !

    IF(PRESENT(Order_O))THEN
       Order=Order_O
       IF(Order/=2)  &
            CALL Halt(' Bad logic in Sort_CellSet. Cells MUST be sorted in descending order! ')
    ELSE
       Order=2
    ENDIF

    CALL New(IPnt  ,CS%NCells)
    CALL New(Vec   ,CS%NCells)
    CALL New(CCarts,(/3,CS%NCells/))
    DO NC=1,CS%NCells
       IPnt%I(NC)= NC
       Vec%D(NC) = SQRT(CS%CellCarts%D(1,NC)**2+CS%CellCarts%D(2,NC)**2+CS%CellCarts%D(3,NC)**2)
    ENDDO
    !
    CALL Sort(Vec,IPnt,CS%NCells,Order)
    !
    DO NC=1,CS%NCells
       CCarts%D(1,NC) = CS%CellCarts%D(1,IPnt%I(NC))
       CCarts%D(2,NC) = CS%CellCarts%D(2,IPnt%I(NC))
       CCarts%D(3,NC) = CS%CellCarts%D(3,IPnt%I(NC))
    ENDDO
    !
    DO NC=1,CS%NCells
       CS%CellCarts%D(1,CS%NCells-NC+1) = CCarts%D(1,NC)
       CS%CellCarts%D(2,CS%NCells-NC+1) = CCarts%D(2,NC)
       CS%CellCarts%D(3,CS%NCells-NC+1) = CCarts%D(3,NC)
    ENDDO
    CALL Delete(IPnt)
    CALL Delete(Vec)
    CALL Delete(CCarts)
    !
  END SUBROUTINE Sort_CellSet
  !--------------------------------------------------------------------------
  ! Test if a Cell in CS has Coordinates (x,y,z)
  !--------------------------------------------------------------------------
  FUNCTION InCell_CellSet(CS,X,Y,Z)
    TYPE(CellSet)      :: CS
    INTEGER            :: NC
    REAL(DOUBLE)       :: X,Y,Z
    LOGICAL            :: InCell_CellSet
    !
    InCell_CellSet = .FALSE.
    DO NC = 1,CS%NCells
       IF(ABS(CS%CellCarts%D(1,NC)-X) < 1.0D-10 .AND. &
            ABS(CS%CellCarts%D(2,NC)-Y) < 1.0D-10 .AND. &
            ABS(CS%CellCarts%D(3,NC)-Z) < 1.0D-10) THEN
          InCell_CellSet = .TRUE.
          RETURN
       ENDIF
    ENDDO
    !
  END FUNCTION InCell_CellSet
  !-------------------------------------------------------------------------------
  ! Convert from Atomic Coordinates  to Fractional Coordinates
  !-------------------------------------------------------------------------------
  FUNCTION AtomToFrac(GM,VecA) RESULT(VecF)
    TYPE(CRDS)                 :: GM
    REAL(DOUBLE),DIMENSION(3)  :: VecA,VecF
    VecF=MATMUL(GM%PBC%InvBoxSh%D,VecA)
  END FUNCTION AtomToFrac
  !-------------------------------------------------------------------------------
  ! Convert from Fractional Coordinates to Atomic Coordinates
  !-------------------------------------------------------------------------------
  FUNCTION FracToAtom(GM,VecF) RESULT(VecA)
    TYPE(CRDS)                 :: GM
    REAL(DOUBLE),DIMENSION(3)  :: VecA,VecF
!!$!    VecA=MATMUL(GM%PBC%BoxShape%D,VecF)
!!$    VecA=MATMUL(VecF,GM%PBC%BoxShape%D)



    VecA(1) = VecF(1)*GM%PBC%BoxShape%D(1,1) + VecF(2)*GM%PBC%BoxShape%D(1,2) + VecF(3)*GM%PBC%BoxShape%D(1,3)
    VecA(2) = VecF(1)*GM%PBC%BoxShape%D(2,1) + VecF(2)*GM%PBC%BoxShape%D(2,2) + VecF(3)*GM%PBC%BoxShape%D(2,3)
    VecA(3) = VecF(1)*GM%PBC%BoxShape%D(3,1) + VecF(2)*GM%PBC%BoxShape%D(3,2) + VecF(3)*GM%PBC%BoxShape%D(3,3)
  END FUNCTION FracToAtom
  !-----------------------------------------------------------------------------------
  ! Subroutine PWrap wraps and cell-centers a the position of a primitive distribution
  !-----------------------------------------------------------------------------------
  SUBROUTINE PWrap(GM,P,Wrap)
    TYPE(CRDS)     :: GM
    TYPE(PrimPair) :: P
    LOGICAL        :: Wrap
    REAL(DOUBLE),DIMENSION(3):: PT1,PT2,PT3,PT4
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    IF(GM%PBC%Dimen>0.AND.Wrap)THEN
       ! Fractional coordinates:
       P%Pw=AtomToFrac(GM,P%P)
       ! Wrapped fractional coordinates:
       CALL FracCyclic(GM,P%Pw)
       ! Convert back to wrapped atomic coordinates
       P%Pw=FracToAtom(GM,P%Pw)
       PwMax=MAX(PwMax,SQRT(DOT_PRODUCT(P%PW,P%PW) ))
       ! DOUBLE CHECK: WRAPPED VECTORS SHOULD NOW BE STATIONARY WRT TO WRAPPING
       ! SMALL ROUND OFF ERRORS CAN SUBVERT THIS, ADDRESSED IN FracCyclic.
       PT1=P%PW
       PT2=AtomToFrac(GM,PT1)
       PT3=PT2
       CALL FracCyclic(GM,PT3)
       PT4=FracToAtom(GM,PT3)
       IF(DOT_PRODUCT(PT4-PT1,PT4-PT1).GT.1D-8)THEN
          WRITE(*,*)'- - - - - - - - - - - - - - - '
          WRITE(*,*)' P  = ',P%P
          WRITE(*,*)' Pw = ',P%Pw
          WRITE(*,*)' '
          WRITE(*,*)' P1 = ',PT1
          WRITE(*,*)' P2 = ',PT2
          WRITE(*,*)' P3 = ',PT3
          WRITE(*,*)' P4 = ',PT4
          STOP
       ENDIF
       ! Cell centering: Note, cell centered primitives are implicit throughout
       ! the new QCTC, and changing this policy will break the code badly!
       P%Pw=P%Pw-GM%PBC%CellCenter%D
    ELSE
       ! Nothing to do here, just fill the vector.
       P%Pw=P%P
    ENDIF
  END SUBROUTINE PWrap
  !-------------------------------------------------------------------------------
  ! In Perodic Systems, Cyclically put positions back in the Fracional Box (in Frac Coord)
  !-------------------------------------------------------------------------------
  SUBROUTINE FracCyclic(GM,VecF)
    TYPE(CRDS)                 :: GM
    REAL(DOUBLE),DIMENSION(3)  :: VecF
    INTEGER                    :: I
    REAL(DOUBLE),PARAMETER     :: Delta=1D-15
    DO I=1,3
       IF(GM%PBC%AutoW%I(I)==1)VecF(I)=MODULO(VecF(I)+Delta,ONE)
    ENDDO
  END SUBROUTINE FracCyclic
  !-------------------------------------------------------------------------------
  ! In Perodic Systems, Cyclically put positions back in the Atomic Box (in Atom Coord)
  !-------------------------------------------------------------------------------
  SUBROUTINE AtomCyclic(GM,VecA)
    TYPE(CRDS)                 :: GM
    REAL(DOUBLE),DIMENSION(3)  :: VecA,VecF
    VecF = AtomToFrac(GM,VecA)
    CALL FracCyclic(GM,VecF)
    VecA = FracToAtom(GM,VecF)
  END SUBROUTINE AtomCyclic
  !===================================================================================
  !
  !===================================================================================
  SUBROUTINE  WrapAtoms(G)
    TYPE(CRDS)     :: G
    INTEGER        :: I
    !   Check for no wrapping at all
    IF(G%PBC%Dimen==0) RETURN
    !   Wrap Fractional coordinates
    DO I=1,G%NAtms
       CALL FracCyclic(G,G%BoxCarts%D(:,I))
    ENDDO
    !   Wrap Atomic  coordinates
    DO I=1,G%NAtms
       CALL AtomCyclic(G,G%Carts%D(:,I))
    ENDDO
    !
  END SUBROUTINE WrapAtoms
  !===================================================================================
  !
  !===================================================================================
  SUBROUTINE  CalFracCarts(GM)
    TYPE(CRDS)                 :: GM
    INTEGER                    :: I
    DO I=1,GM%NAtms
       !      Compute fracts from Carts:
       GM%BoxCarts%D(:,I) = AtomToFrac(GM,GM%Carts%D(:,I))
    ENDDO
  END SUBROUTINE CalFracCarts
  !===================================================================================
  !
  !===================================================================================
  SUBROUTINE  CalAtomCarts(GM)
    TYPE(CRDS)                 :: GM
    INTEGER                    :: I
    DO I=1,GM%NAtms
       !      Only wrap carts
       GM%Carts%D(:,I)   = FracToAtom(GM,GM%BoxCarts%D(:,I))
    ENDDO
  END SUBROUTINE CalAtomCarts
  !===================================================================================
  !
  !===================================================================================
  SUBROUTINE  Translate(GM,ATvec)
    TYPE(CRDS)                 :: GM
    REAL(DOUBLE),DIMENSION(3)  :: ATvec,FTvec
    INTEGER                    :: I
    !
    FTvec(:) = AtomToFrac(GM,ATvec(:))
    !
    !   Tranaslate The Atoms
    !
    DO I=1,GM%NAtms
       GM%Carts%D(:,I)      = GM%Carts%D(:,I)    + ATvec(:)
       GM%Carts%D(:,I)    = GM%Carts%D(:,I)  + ATvec(:)
       GM%BoxCarts%D(:,I)   = GM%BoxCarts%D(:,I) + FTvec(:)
    ENDDO
    !
  END SUBROUTINE Translate
  !-------------------------------------------------------------------------------
  ! Test to See if in the Box (Fractional Coordinates)
  !-------------------------------------------------------------------------------
  FUNCTION InFracBox(GM,VecF)
    TYPE(CRDS)                 :: GM
    LOGICAL                    :: InFracBox
    INTEGER                    :: I
    REAL(DOUBLE),DIMENSION(3)  :: VecF

    InFracBox = .TRUE.
    DO I=1,3
       IF(GM%PBC%AutoW%I(I)==1) THEN
          IF(VecF(I) < zero .OR. VecF(I) > one) THEN
             InFracBox = .FALSE.
             RETURN
          ENDIF
       ENDIF
    ENDDO

  END FUNCTION InFracBox
  !-------------------------------------------------------------------------------
  ! Test to See if in the Box (Atomic Coordinates)
  !-------------------------------------------------------------------------------
  FUNCTION InAtomBox(GM,VecA)
    TYPE(CRDS)                 :: GM
    LOGICAL                    :: InAtomBox
    REAL(DOUBLE),DIMENSION(3)  :: VecA,VecF

    VecF = AtomToFrac(GM,VecA)
    InAtomBox = InFracBox(GM,VecF)

  END FUNCTION InAtomBox
  !-------------------------------------------------------------------------------
  FUNCTION CellCenter(PBC)
    TYPE(PBCInfo)                     :: PBC
    REAL(DOUBLE),DIMENSION(3)   :: CellCenter
    INTEGER                     :: I,J
    !   Find the center of the cell
    CellCenter=0D0
    DO I=1,3
       IF(PBC%AutoW%I(I)==1)THEN
          DO J=1,3
             IF(PBC%AutoW%I(J)==1)THEN
                CellCenter(I)=CellCenter(I)+Half*PBC%BoxShape%D(I,J)
             ENDIF
          ENDDO
       ENDIF
    ENDDO
  END FUNCTION CellCenter
  !-------------------------------------------------------------------------------
  FUNCTION CellVolume(BoxShape,AutoW)
    REAL(DOUBLE)                :: CellVolume
    REAL(DOUBLE),DIMENSION(3,3) :: BoxShape
    REAL(DOUBLE),DIMENSION(3)   :: AB
    INTEGER,DIMENSION(3)        :: AutoW
    INTEGER                     :: I,J,D

    D=0
    DO I=1,3
       D = D + AutoW(I)
    ENDDO

    CellVolume = -1.0D0

    IF(D==1) THEN
       IF(AutoW(1)==1) CellVolume=BoxShape(1,1)
       IF(AutoW(2)==1) CellVolume=BoxShape(2,2)
       IF(AutoW(3)==1) CellVolume=BoxShape(3,3)
    ELSEIF(D==2) THEN
       IF(AutoW(1)==0) THEN
          CellVolume = BoxShape(2,2)*BoxShape(3,3)-BoxShape(3,2)*BoxShape(2,3)
       ENDIF
       IF(AutoW(2)==0) THEN
          CellVolume = BoxShape(1,1)*BoxShape(3,3)-BoxShape(1,3)*BoxShape(3,1)
       ENDIF
       IF(AutoW(3)==0) THEN
          CellVolume = BoxShape(1,1)*BoxShape(2,2)-BoxShape(1,2)*BoxShape(2,1)
       ENDIF
    ELSEIF(D==3) THEN
       AB=CROSS_PRODUCT(BoxShape(:,1),BoxShape(:,2))
       CellVolume = DOT_PRODUCT(AB,BoxShape(:,3))
    ENDIF
    !
  END FUNCTION CellVolume
  !-------------------------------------------------------------------------------
  FUNCTION DivCellVolume(BoxShape,AutoW) RESULT(DivCV)
    REAL(DOUBLE),DIMENSION(3,3) :: BoxShape,DivCV,Temp,P,M
    INTEGER,DIMENSION(3)        :: AutoW
    REAL(DOUBLE)                :: Phase
    INTEGER                     :: I,J,D
    !
    D=0
    DO I=1,3
       D = D + AutoW(I)
    ENDDO
    !
    DivCV = Zero
    IF(D==1) THEN
       IF(AutoW(1)==1) DivCV(1,1) = One
       IF(AutoW(2)==1) DivCV(2,2) = One
       IF(AutoW(3)==1) DivCV(3,3) = One
    ELSEIF(D==2) THEN
       IF(AutoW(1)==0) THEN
          DivCV(2,2) =  BoxShape(3,3)
          DivCV(2,3) = -BoxShape(3,2)
          DivCV(3,2) = -BoxShape(2,3)
          DivCV(3,3) =  BoxShape(2,2)
       ENDIF
       IF(AutoW(2)==0) THEN
          DivCV(1,1) =  BoxShape(3,3)
          DivCV(1,3) = -BoxShape(3,1)
          DivCV(3,1) = -BoxShape(1,3)
          DivCV(3,3) =  BoxShape(1,1)
       ENDIF
       IF(AutoW(3)==0) THEN
          DivCV(1,1) =  BoxShape(2,2)
          DivCV(1,2) = -BoxShape(2,1)
          DivCV(2,1) = -BoxShape(1,2)
          DivCV(2,2) =  BoxShape(1,1)
       ENDIF
    ELSEIF(D==3) THEN

       DO I=1,3
          DO J=1,3
             Temp       = BoxShape
             Temp(:,J)  = Zero
             Temp(I,J)  = One
             DivCV(I,J) = CellVolume(Temp,AutoW)
          ENDDO
       ENDDO

!!$       DO I=1,3
!!$          DO J=1,3
!!$             P=BoxShape
!!$             M=BoxShape
!!$             P(I,J)=P(I,J)+1D-5
!!$             M(I,J)=M(I,J)-1D-5
!!$             DivCV(I,J) = (CellVolume(P,AutoW)+CellVolume(M,AutoW))/(2D-5)
!!$             WRITE(*,*)I,J,' div2 = ',DivCV(I,J)
!!$          ENDDO
!!$       ENDDO

    ENDIF

    Phase = CellVolume(BoxShape,AutoW)/ABS(CellVolume(BoxShape,AutoW))
    DivCV = Phase*DivCV
    !
  END FUNCTION DivCellVolume

  !--------------------------------------------------------------------------------
  SUBROUTINE BoxParsToCart(Vec,BoxShape)
    REAL(DOUBLE),DIMENSION(6)   :: Vec
    REAL(DOUBLE),DIMENSION(3,3) :: BoxShape
    BoxShape=Zero
    BoxShape(1,1)=Vec(1)
    BoxShape(1,2)=Vec(2)*COS(Vec(6))
    BoxShape(2,2)=Vec(2)*SIN(Vec(6))
    BoxShape(1,3)=Vec(3)*COS(Vec(5))
    BoxShape(2,3)=(Vec(2)*Vec(3)*COS(Vec(4)) &
         -BoxShape(1,2)*BoxShape(1,3))/BoxShape(2,2)
    BoxShape(3,3)=SQRT(Vec(3)**2-BoxShape(1,3)**2-BoxShape(2,3)**2)
  END SUBROUTINE BoxParsToCart
  !--------------------------------------------------------------------
  SUBROUTINE CalcBoxPars(Vec,BoxShape)
    REAL(DOUBLE),DIMENSION(6)  :: Vec
    REAL(DOUBLE),DIMENSION(3)  :: VecA,VecB,VecC,VecCr
    REAL(DOUBLE),DIMENSION(:,:):: BoxShape
    !
    Vec(1)=SQRT(DOT_PRODUCT(BoxShape(1:3,1),BoxShape(1:3,1)))
    Vec(2)=SQRT(DOT_PRODUCT(BoxShape(1:3,2),BoxShape(1:3,2)))
    Vec(3)=SQRT(DOT_PRODUCT(BoxShape(1:3,3),BoxShape(1:3,3)))
    VecA=BoxShape(1:3,1)/Vec(1)
    VecB=BoxShape(1:3,2)/Vec(2)
    VecC=BoxShape(1:3,3)/Vec(3)
    Vec(4)=DOT_PRODUCT(VecB,VecC)
    Vec(5)=DOT_PRODUCT(VecC,VecA)
    Vec(6)=DOT_PRODUCT(VecA,VecB)
    Vec(4)=ACOS(Vec(4))
    Vec(5)=ACOS(Vec(5))
    Vec(6)=ACOS(Vec(6))
  END SUBROUTINE CalcBoxPars
  !-------------------------------------------------------------------------------
  !  Calculate Inverse of the Box Shape (Transpose of the reciprocal lattice vectors)
  !  Note:  For 1-D and 2-D periodicity, there is an assumption of a 90 degree
  !  lattice vector, of unit 1 if not specified, so that the iverse box shape works
  !  out. See MondoSCF/ParsePeriodic for details.
  !-------------------------------------------------------------------------------
  FUNCTION InverseBoxShape(Mat,Dim) RESULT(InvMat)
    REAL(DOUBLE)                    :: Det,Norm
    REAL(DOUBLE),DIMENSION(3,3)     :: Mat,InvMat
    INTEGER                         :: Dim
!!$    IF(Dim==1)THEN
!!$       InvMat=0D0
!!$       InvMat(1,1)=One/Mat(1,1)  ! Periodic along x-axis
!!$    ELSEIF(Dim==2)THEN
!!$       InvMat=0D0
!!$       InvMat(1:2,1:2)=InverseMatrix2x2(Mat(1:2,1:2))  ! Periodic in x-y plane
!!$    ELSE
    InvMat=InverseMatrix3x3(Mat)
!!$    ENDIF
  END FUNCTION InverseBoxShape

  FUNCTION InverseMatrix2x2(Mat) RESULT(InvMat)
    REAL(DOUBLE)                    :: Det,Norm
    REAL(DOUBLE),DIMENSION(2,2)     :: Mat,InvMat
    !
    Det  = Mat(1,1)*Mat(2,2)-Mat(1,2)*Mat(2,1)
    Norm = One/Det
    InvMat(1,1) = Norm*Mat(2,2)
    InvMat(1,2) =-Norm*Mat(1,2)
    InvMat(2,1) =-Norm*Mat(2,1)
    InvMat(2,2) = Norm*Mat(1,1)
  END FUNCTION InverseMatrix2x2

  FUNCTION InverseMatrix3x3(Mat) RESULT(InvMat)
    REAL(DOUBLE)                    :: Det,Norm
    REAL(DOUBLE),DIMENSION(3,3)     :: Mat,InvMat
    !
    Det  = Mat(1,1)*Mat(2,2)*Mat(3,3) + Mat(1,2)*Mat(2,3)*Mat(3,1) + Mat(1,3)*Mat(2,1)*Mat(3,2) &
         - Mat(1,3)*Mat(2,2)*Mat(3,1) - Mat(1,1)*Mat(2,3)*Mat(3,2) - Mat(1,2)*Mat(2,1)*Mat(3,3)
    Norm = One/Det
    InvMat(1,1) = Norm*(Mat(2,2)*Mat(3,3) - Mat(2,3)*Mat(3,2))
    InvMat(1,2) = Norm*(Mat(1,3)*Mat(3,2) - Mat(1,2)*Mat(3,3))
    InvMat(1,3) = Norm*(Mat(1,2)*Mat(2,3) - Mat(1,3)*Mat(2,2))
    InvMat(2,1) = Norm*(Mat(2,3)*Mat(3,1) - Mat(2,1)*Mat(3,3))
    InvMat(2,2) = Norm*(Mat(1,1)*Mat(3,3) - Mat(1,3)*Mat(3,1))
    InvMat(2,3) = Norm*(Mat(1,3)*Mat(2,1) - Mat(1,1)*Mat(2,3))
    InvMat(3,1) = Norm*(Mat(2,1)*Mat(3,2) - Mat(2,2)*Mat(3,1))
    InvMat(3,2) = Norm*(Mat(1,2)*Mat(3,1) - Mat(1,1)*Mat(3,2))
    InvMat(3,3) = Norm*(Mat(1,1)*Mat(2,2) - Mat(1,2)*Mat(2,1))
    !
  END FUNCTION InverseMatrix3x3



  SUBROUTINE New_CellSet_Sphere(CS,AW,MAT,Radius,Rmin_O,MaxCell_O)
    TYPE(CellSet)                        :: CS
    INTEGER,DIMENSION(3)                 :: AW
    REAL(DOUBLE),DIMENSION(3,3)          :: MAT
    REAL(DOUBLE)                         :: Radius,Radius_min
    REAL(DOUBLE),OPTIONAL                :: Rmin_O
    INTEGER,OPTIONAL                     :: MaxCell_O

    INTEGER                              :: I,J,K
    INTEGER                              :: IXM,IYM,IZM,NCell
    REAL(DOUBLE)                         :: X,Y,Z,Rad,R


    ! Oh god, I hate this routine. It is not as bad as MakeDivTensor2D tho,
    ! which is the only routine that calls this one.  Its clear how to
    ! go forward, but it will take some effort.


    !
    IF(PRESENT(Rmin_O)) THEN
       Radius_min = Rmin_O
    ELSE
       Radius_min = Zero
    ENDIF
    !

    IXM = 0
    IYM = 0
    IZM = 0
    IF(AW(1)==1) IXM = 2*(1+INT(Radius/SQRT(MAT(1,1)**2+MAT(2,1)**2+MAT(3,1)**2)))
    IF(AW(2)==1) IYM = 2*(1+INT(Radius/SQRT(MAT(1,2)**2+MAT(2,2)**2+MAT(3,2)**2)))
    IF(AW(3)==1) IZM = 2*(1+INT(Radius/SQRT(MAT(1,3)**2+MAT(2,3)**2+MAT(3,3)**2)))


    !
    NCell = 0
    DO I=-IXM,IXM
       DO J=-IYM,IYM
          DO K=-IZM,IZM
             X  = I*MAT(1,1)+J*MAT(1,2)+K*MAT(1,3)
             Y  = I*MAT(2,1)+J*MAT(2,2)+K*MAT(2,3)
             Z  = I*MAT(3,1)+J*MAT(3,2)+K*MAT(3,3)
             R = SQRT(X*X+Y*Y+Z*Z)
             IF(R .GE. Radius_min .AND. R .LT. Radius) THEN
                NCell = NCell+1
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !
    !
    IF(PRESENT(MaxCell_O)) THEN
       IF(NCELL > MaxCell_O) THEN
          WRITE(*,*) 'NCell   = ',NCell
          WRITE(*,*) 'MaxCell = ',MaxCell_O
          CALL Halt('NCELL is Greater then MaxCell')
       ELSE
          CALL New(CS,MaxCell_O)
          CS%NCells=NCELL
       ENDIF
    ELSE
       CALL New(CS,NCELL)
       CS%NCells=NCELL
    ENDIF

    NCELL=0
    DO I=-IXM,IXM
       DO J=-IYM,IYM
          DO K=-IZM,IZM
             X  = I*MAT(1,1)+J*MAT(1,2)+K*MAT(1,3)
             Y  = I*MAT(2,1)+J*MAT(2,2)+K*MAT(2,3)
             Z  = I*MAT(3,1)+J*MAT(3,2)+K*MAT(3,3)
             R = SQRT(X*X+Y*Y+Z*Z)
             IF(R .GE. Radius_min .AND. R .LT. Radius) THEN
                NCELL = NCELL+1
                CS%CellCarts%D(1,NCELL)= X
                CS%CellCarts%D(2,NCELL)= Y
                CS%CellCarts%D(3,NCELL)= Z
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    !
  END SUBROUTINE New_CellSet_Sphere
!!$  !--------------------------------------------------------------------------
!!$  ! Set up the CellSet from {-N,N} minus the Cells in {-M,M}
!!$  !--------------------------------------------------------------------------
!!$  SUBROUTINE New_CellSet_Cube(CS,AW,MAT,N,MaxCell_O)
!!$    TYPE(CellSet)                        :: CS
!!$    INTEGER,DIMENSION(3)                 :: AW
!!$    REAL(DOUBLE),DIMENSION(3,3)          :: MAT
!!$    INTEGER,DIMENSION(3)                 :: N
!!$    INTEGER,OPTIONAL                     :: MaxCell_O
!!$    !
!!$    INTEGER                              :: I,J,K,NCELL
!!$    REAL(DOUBLE)                         :: X,Y,Z
!!$    !
!!$
!!$    WRITE(*,*)' OLD CELLSET CUBE HAS BEEN DEPRICATED, PLEASE GET RID OF THIS CALL!! '
!!$    WRITE(*,*)' OLD CELLSET CUBE HAS BEEN DEPRICATED, PLEASE GET RID OF THIS CALL!! '
!!$    WRITE(*,*)' OLD CELLSET CUBE HAS BEEN DEPRICATED, PLEASE GET RID OF THIS CALL!! '
!!$    WRITE(*,*)' OLD CELLSET CUBE HAS BEEN DEPRICATED, PLEASE GET RID OF THIS CALL!! '
!!$    WRITE(*,*)' OLD CELLSET CUBE HAS BEEN DEPRICATED, PLEASE GET RID OF THIS CALL!! '
!!$
!!$    IF(N(1) .LT. 0) CALL Halt('N(1) is Miss Dimensioning in New_CellSet_Box')
!!$    IF(N(2) .LT. 0) CALL Halt('N(2) is Miss Dimensioning in New_CellSet_Box')
!!$    IF(N(3) .LT. 0) CALL Halt('N(3) is Miss Dimensioning in New_CellSet_Box')
!!$    !
!!$    IF(AW(1)==0) N(1)=0
!!$    IF(AW(2)==0) N(2)=0
!!$    IF(AW(3)==0) N(3)=0
!!$    !
!!$    NCELL = 0
!!$    DO I = -N(1),N(1)
!!$       DO J = -N(2),N(2)
!!$          DO K = -N(3),N(3)
!!$             NCELL = NCELL+1
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
!!$    !
!!$    IF(PRESENT(MaxCell_O)) THEN
!!$       IF(NCELL > MaxCell_O) THEN
!!$          WRITE(*,*) 'NCELL   = ',NCELL
!!$          WRITE(*,*) 'MaxCell = ',MaxCell_O
!!$          CALL Halt('NCELL is Greater then MaxCell')
!!$       ELSE
!!$          CALL New_CellSet(CS,MaxCell_O)
!!$          CS%NCells=NCELL
!!$       ENDIF
!!$    ELSE
!!$       CALL New_CellSet(CS,NCELL)
!!$       CS%NCells=NCELL
!!$    ENDIF
!!$    !
!!$    NCELL = 0
!!$    DO I = -N(1),N(1)
!!$       DO J = -N(2),N(2)
!!$          DO K = -N(3),N(3)
!!$             NCELL = NCELL+1
!!$             X  = I*MAT(1,1)+J*MAT(1,2)+K*MAT(1,3)
!!$             Y  = I*MAT(2,1)+J*MAT(2,2)+K*MAT(2,3)
!!$             Z  = I*MAT(3,1)+J*MAT(3,2)+K*MAT(3,3)
!!$             CS%CellCarts%D(1,NCELL)= X
!!$             CS%CellCarts%D(2,NCELL)= Y
!!$             CS%CellCarts%D(3,NCELL)= Z
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
!!$    !
!!$  END SUBROUTINE New_CellSet_Cube

END MODULE PBC
