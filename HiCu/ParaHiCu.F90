! Parallel HiCu  --- author: Chee K. Gan (25 June 2002)
MODULE ParallelHiCu

#ifdef PARALLEL 
  USE DerivedTypes
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  USE Indexing
  USE InOut
  USE Macros
  USE CubeGrid
  USE BoundingBox
  USE RhoTree
  USE Functionals
  USE Thresholding
  USE HiCuThresholds
  USE AtomPairs
  USE SpecFun
  USE PrettyPrint
  USE CubeTree
  USE KxcGen
  IMPLICIT NONE
!===============================================================================
  TYPE(DBL_RNK2)::LCoor,RCoor       ! Coordinates of the bounding boxes
  TYPE(DBL_RNK2)::RepLCoor,RepRCoor ! Coordinates of the repartitioned boxes
  INTEGER::NVol                     ! number of the bounding boxes
  REAL(DOUBLE)::MyLeavesTm          ! time to create all leaves in a subvolume
  TYPE(DBL_VECT)::VolLeavesTm       ! the array to store MyLeavesTm
  REAL(DOUBLE)::VolTm               ! total time to work on a subvolume
  TYPE(DBL_VECT)::TmNode            ! the array to stor VolTm
  REAL(DOUBLE)::BoxPoint(6)
  INTEGER,PARAMETER::LCoor_T=26,RCoor_T=27,Time_T=28,LeavesTm_T=29
   
!===============================================================================
  CONTAINS
!===============================================================================

  SUBROUTINE GetBBox()
    LOGICAL::Exist

    NVol = NPrc
    WRITE(*,*) 'MyID = ', MyID, ', NVol = ', NVol
    CALL New(LCoor,(/3,NVol/))
    CALL New(RCoor,(/3,NVol/))
    INQUIRE(FILE='BegEnd.dat',EXIST=Exist)

    IF(Exist) THEN
      IF(MyID == 0) THEN
        WRITE(*,*) 'BegEnd.dat exists. It is going to be opened..'
        OPEN(UNIT=53,FILE='BegEnd.dat',FORM='formatted',Status='unknown')
        READ(53,*) NVol
        IF(NVol /= NPrc) THEN
          WRITE(*,*) 'NVol = ',NVol, ',  NPrc = ',NPrc
          STOP 'ERROR: BegEnd.dat -- NVol is not equal to NPrc!'
        ENDIF
        READ(53,*) LCoor%D(1:3,1:NVol)
        READ(53,*) RCoor%D(1:3,1:NVol)
        CLOSE(53)
      ENDIF
    ELSE
      ! assign the first box
      IF(MyID == 0) THEN
        WRITE(*,*) 'BegEnd.dat does not exists.'
        LCoor%D(1:3,1) = RhoRoot%Box%BndBox(1:3,1)
        RCoor%D(1:3,1) = RhoRoot%Box%BndBox(1:3,2)
        CALL Get_SubVol(LCoor,RCoor,NVol)
      ENDIF
    ENDIF

  END SUBROUTINE GetBBox

!===============================================================================
  SUBROUTINE SendBBox()
    INTEGER::I,IErr
    REAL(DOUBLE)::DblArr(3)
    INTEGER,DIMENSION(MPI_STATUS_SIZE)::Status
    IF(MyID == ROOT) THEN
      DO I = 1, NPrc-1
        DblArr(1:3) = LCoor%D(1:3,I+1)
        CALL MPI_Send(DblArr(1),3,MPI_DOUBLE_PRECISION,I,LCoor_T,MPI_COMM_WORLD,IErr)
        DblArr(1:3) = RCoor%D(1:3,I+1)
        CALL MPI_Send(DblArr(1),3,MPI_DOUBLE_PRECISION,I,RCoor_T,MPI_COMM_WORLD,IErr)
      ENDDO
    ELSE
      CALL MPI_Recv(DblArr(1),3,MPI_DOUBLE_PRECISION,0,LCoor_T,MPI_COMM_WORLD,Status,IErr)
      LCoor%D(1:3,1) = DblArr(1:3)
      CALL MPI_Recv(DblArr(1),3,MPI_DOUBLE_PRECISION,0,RCoor_T,MPI_COMM_WORLD,Status,IErr)
      RCoor%D(1:3,1) = DblArr(1:3)
    ENDIF
  END SUBROUTINE SendBBox

!===============================================================================
  SUBROUTINE WorkBBox(Kxc)
    TYPE(BBox) :: WBox
    REAL(DOUBLE)::TotRho,TotExc,SubVolRho,SubVolExc
    REAL(DOUBLE)::TmBegM,TmEndM
    TYPE(BCSR):: Kxc
    
    WBox%BndBox(1:3,1) = LCoor%D(1:3,1)
    WBox%BndBox(1:3,2) = RCoor%D(1:3,1)
    CALL CalCenterAndHalf(WBox)
    TmBegM = MPI_WTime()
    CALL GridGen(WBox,SubVolRho,SubVolExc)
    MyLeavesTm = LeavesTmCount(CubeRoot)
    CALL MakeKxc(Kxc,CubeRoot)
    TmEndM = MPI_WTime()
    VolTm = TmEndM-TmBegM
    TotRho = Reduce(SubVolRho)
    TotExc = Reduce(SubVolExc)
    IF(MyID == ROOT) THEN
      WRITE(*,*) 'TotRho = ',TotRho
      WRITE(*,*) 'TotExc = ',TotExc
    ENDIF
  END SUBROUTINE WorkBBox

!===============================================================================
  SUBROUTINE CollectLeavesTime()
    INTEGER::I,IErr
    INTEGER,DIMENSION(MPI_STATUS_SIZE)::Status

    IF(MyID == 0) THEN
      CALL New(VolLeavesTm,NPrc)
      VolLeavesTm%D(1) = MyLeavesTm
      DO I = 1, NPrc-1
        CALL MPI_Recv(VolLeavesTm%D(I+1),1,MPI_DOUBLE_PRECISION,I,LeavesTm_T,MPI_COMM_WORLD,Status,IErr)
      ENDDO
      OPEN(unit=53,file='VolLeavesTm.dat',status='unknown')
      DO I = 1, NPrc
        WRITE(53,*) I, VolLeavesTm%D(I)
      ENDDO
      CLOSE(53)
    ELSE
      CALL MPI_Send(MyLeavesTm,1,MPI_DOUBLE_PRECISION,0,LeavesTm_T,MPI_COMM_WORLD,IErr)
    ENDIF
  END SUBROUTINE CollectLeavesTime

!===============================================================================
  SUBROUTINE RepartitionVol()
    
    LOGICAL::Busy
    TYPE(DBL_VECT)::RepLeavesTm
    INTEGER::SmallN,I,J,Ierr,RootFindIter,Stage,PIndex,CIndex,DirInt
    REAL(DOUBLE)::x0,x1,f0,f1,x2,f2,oldf2,LinDim,MaxDim,OldVol,ThisVol,&
      AbsSum,LocalLeavesTm,LeavesTmInside,OrigLeafTm,NewVol
    REAL(DOUBLE)::Power2(0:31)
    REAL(DOUBLE),PARAMETER::RootTau=1.0D-4

    CALL New(RepLCoor,(/3,NVol/))
    CALL New(RepRCoor,(/3,NVol/))
    CALL New(RepLeavesTm,NVol)

    Busy = .TRUE.
    IF(MyID == 0) THEN
      RepLCoor%D(1:3,1) = RhoRoot%Box%BndBox(1:3,1)
      RepRCoor%D(1:3,1) = RhoRoot%Box%BndBox(1:3,2)
      RepLeavesTm%D(1) = Sum(VolLeavesTm%D(:))
      WRITE(*,*) 'The total leave times is ',RepLeavesTm%D(1)
      NewVol = 1.0D0
      DO I = 1, 3
        NewVol = NewVol*(RepRCoor%D(I,1)-RepLCoor%D(I,1))
      ENDDO
      WRITE(*,*) 'The initial volume to start with is ',NewVol
      
      Power2(0) = 1
      DO I = 1, 31
        Power2(I) = Power2(I-1)*2
      ENDDO
      SmallN = NINT(LOG(NVol*1.0D0)/LOG(2.0D0))
      IF(NVol /= Power2(SmallN)) THEN
        WRITE(*,*) 'NVol is ',NVol
        WRITE(*,*) 'NVol is not a power of 2!'
        STOP 
      ENDIF
  
      OPEN(unit=53,FILE='DirInt.dat',form='formatted',status='unknown')
      DO Stage = 1, SmallN
        WRITE(*,*) 'Stage = ',Stage
        DO PIndex = 1, Power2(Stage-1)
          CIndex = PIndex + Power2(Stage-1)
          !! copy the coordinates to the child box
          RepLCoor%D(1:3,CIndex) = RepLCoor%D(1:3,PIndex)
          RepRCoor%D(1:3,CIndex) = RepRCoor%D(1:3,PIndex)
  
          !! direction to split
          MaxDim = -1.0D0
          DO I = 1, 3
            LinDim = RepRCoor%D(I,PIndex)-RepLCoor%D(I,PIndex)
            IF(LinDim <= 0) STOP 'ERROR: LinDim is not positive in RepartionVol'
            IF(LinDim > MaxDim) THEN
              MaxDim = LinDim
              DirInt = I
            ENDIF
          ENDDO
          WRITE(53,*) DirInt
  
          OrigLeafTm = RepLeavesTm%D(PIndex)
          x0 = RepLCoor%D(DirInt,PIndex)
          x1 = RepRCoor%D(DirInt,PIndex)
          f0 = -0.50D0
          f1 =  0.50D0
  
          RootFindIter = 0
          f2 = 1000.0D0
          DO 
            RootFindIter = RootFindIter + 1
            IF(RootFindIter > 100) THEN 
              WRITE(*,*) 'Stage = ',Stage,', PIndex = ',PIndex,', CIndex =',CIndex
              WRITE(*,*) 'ERROR: Regula-falsi does not converged.'
            ENDIF
            x2 = x1 - f1*(x1-x0)/(f1-f0)
            !! calculate f2 now
            !! first define the box
            RepRCoor%D(DirInt,PIndex) = x2
            !! copy the coordinates of RepLCoor and RepRCoor to define a box
            BoxPoint(1:3) = RepLCoor%D(1:3,PIndex)
            BoxPoint(4:6) = RepRCoor%D(1:3,PIndex)
            CALL MPI_BCast(BoxPoint(1),6,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
            LocalLeavesTm = CalLeavesTm(BoxPoint)
            LeavesTmInside = Reduce(LocalLeavesTm)
            oldf2 = f2
            f2 = LeavesTmInside*1.0D0/(OrigLeafTm*1.0D0)-0.5D0
            IF(ABS(f2) < RootTau .OR. ABS(f2-oldf2) < 1.0D-5) THEN
              ! RepRCoor%D(1:3,PIndex) has been defined
              WRITE(*,*) 'RootFindIter = ',RootFindIter,', f2 = ',f2
              RepLCoor%D(DirInt,CIndex) = x2
              RepLeavesTm%D(PIndex) = LeavesTmInside
              RepLeavesTm%D(CIndex) = OrigLeafTm - LeavesTmInside
              EXIT
            ELSE
              !! convergence not achieved, going to the next iteration
              WRITE(*,*) 'RootFindIter = ',RootFindIter,', f2 = ',f2
              IF(f2*f1 < 0.0D0) THEN
                x0 = x1; f0 = f1
                x1 = x2; f1 = f2
              ELSE
                x1 = x2; f1 = f2
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      CLOSE(UNIT=53)
      
      OPEN(unit=53,FILE='LeafTmInVol.dat',form='formatted',status='unknown')
      DO I = 1, NVol
        WRITE(53,*) I, RepLeavesTm%D(I)
      ENDDO
      CLOSE(53)
      !! ask all processors to quit now!
      BoxPoint(1:6) = 0.0D0
      CALL MPI_BCast(BoxPoint(1),6,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
  
      !! check the volume before and after are the same
      OldVol = ZERO
      OPEN(unit=53,FILE='RepartionVol.dat',form='formatted',status='unknown')
      DO I = 1, NVol
        ThisVol = 1.0D0
        DO J = 1, 3
          ThisVol = ThisVol*(RCoor%D(J,I)-LCoor%D(J,I))
        ENDDO
        WRITE(53,*) I, ThisVol
        OldVol = OldVol + ThisVol
      ENDDO
      CLOSE(53)
  
      NewVol = ZERO
      DO I = 1, NVol
        ThisVol = 1.0D0
        DO J = 1, 3
          ThisVol = ThisVol*(RepRCoor%D(J,I)-RepLCoor%D(J,I))
        ENDDO
        NewVol = NewVol + ThisVol
      ENDDO
      WRITE(*,*) 'OldVol = ',OldVol, ', NewVol = ', NewVol
      WRITE(*,*) 'ABS(DIFF) = ',ABS(NewVol-OldVol)
      ! WRITE to a file: unit=53 is used.
      open(unit=53,file='BegEnd.dat',form='formatted',status='unknown')
      WRITE(53,*) NVol
    
      !! check for the positive definiteness 
      DO I = 1, NVol
        DO J = 1, 3
          LinDim = RepRCoor%D(J,I)-RepLCoor%D(J,I)
          IF(LinDim <= 0) THEN
            WRITE(*,*) 'LinDim is non-positive!'
            WRITE(*,*) 'LinDim = ',LinDim
            WRITE(*,*) 'Box = ',I, ',  Direction = ',J
            STOP
          ENDIF
        ENDDO
      ENDDO
      WRITE(53,*) RepLCoor%D(1:3,1:NVol)
      WRITE(53,*) RepRCoor%D(1:3,1:NVol)
      CLOSE(53)
    ELSE
      DO 
        CALL MPI_BCast(BoxPoint(1),6,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
        AbsSum = Sum( ABS(BoxPoint(:)))
        ! WRITE(*,*) 'MyID = ',MyID,', AbsSum = ', AbsSum
        IF(ABS(AbsSum) < 1.0D-14) THEN
          Busy = .FALSE.
        ENDIF
        IF(.NOT. Busy) EXIT
        LocalLeavesTm = CalLeavesTm(BoxPoint)
        LeavesTmInside = Reduce(LocalLeavesTm)
      ENDDO
    ENDIF
    CALL Delete(RepLCoor)
    CALL Delete(RepRCoor)
    CALL Delete(RepLeavesTm)
  END SUBROUTINE RepartitionVol
  
!===============================================================================
     FUNCTION CalLeavesTm(BoxPoint)
       REAL(DOUBLE)::BoxPoint(6),CalLeavesTm,LE(3),UE(3)
       INTEGER::I

       LE(1:3) = CubeRoot%Box%BndBox(1:3,1)
       UE(1:3) = CubeRoot%Box%BndBox(1:3,2)

       DO I = 1, 3
         IF(UE(I) <= LE(I)) STOP 'ERROR: In CalLeavesTm, UE <= LE !'
         IF(BoxPoint(I+3) <= BoxPoint(I)) STOP 'ERROR: BoxPoint(I+3) <= BoxPoint(I)'
       ENDDO
       IF(UE(1) <= BoxPoint(1) .OR. LE(1) >= BoxPoint(4) .OR. &
          UE(2) <= BoxPoint(2) .OR. LE(2) >= BoxPoint(5) .OR. &
          UE(3) <= BoxPoint(3) .OR. LE(3) >= BoxPoint(6)) THEN
         CalLeavesTm = 0
       ELSE IF(LE(1) >= BoxPoint(1) .AND. UE(1) <= BoxPoint(4) .AND. &
               LE(2) >= BoxPoint(2) .AND. UE(2) <= BoxPoint(5) .AND. &
               LE(3) >= BoxPoint(3) .AND. UE(3) <= BoxPoint(6)) THEN
         CalLeavesTm = MyLeavesTm
       ELSE
         CalLeavesTm = RecurCalLeavesTm(CubeRoot)
       ENDIF
     END FUNCTION CalLeavesTm
     
!===============================================================================
     RECURSIVE FUNCTION RecurCalLeavesTm(Cube)
       TYPE(CubeNode),Pointer :: Cube
       REAL(DOUBLE)::RecurCalLeavesTm,OverlapV,CubeV
       REAL(DOUBLE)::C1,C2,C3,NC1,NC2,NC3,LE(3),UE(3),P(3),Q(3),Ovlap(3)
       INTEGER::I

#undef CenterConsideration
#define FractionalConsideration
       
#ifdef CenterConsideration
       IF(Cube%Leaf) THEN
         C1 = Cube%Box%Center(1)
         C2 = Cube%Box%Center(2)
         C3 = Cube%Box%Center(3)

         IF(C1 >= BoxPoint(1) .AND. C1 < BoxPoint(4) .AND. &
            C2 >= BoxPoint(2) .AND. C2 < BoxPoint(5) .AND. &
            C3 >= BoxPoint(3) .AND. C3 < BoxPoint(6)) THEN
           RecurCalLeavesTm = Cube%LayGridCost
         ELSE
           RecurCalLeavesTm = 0.0D0
         ENDIF
#endif
#ifdef FractionalConsideration

       IF(Cube%Leaf) THEN
!! two possibilities: This leaf partially overlaps with BoxPoint or not!
         LE(1:3) = Cube%Box%BndBox(1:3,1)
         UE(1:3) = Cube%Box%BndBox(1:3,2)

         P(1:3) = BoxPoint(1:3)
         Q(1:3) = BoxPoint(4:6)

         IF(UE(1) <= P(1) .OR. LE(1) >= Q(1) .OR. &
            UE(2) <= P(2) .OR. LE(2) >= Q(2) .OR. &
            UE(3) <= P(3) .OR. LE(3) >= Q(3)) THEN
           RecurCalLeavesTm = 0.0D0
         ELSE 
           DO I = 1, 3
             !! I-direction Overlap
             IF(LE(I) < P(I)) THEN
               IF(UE(I) <= Q(I)) THEN
                 Ovlap(I) = UE(I) - P(I)
               ELSE
                 Ovlap(I) = Q(I) - P(I)
               ENDIF
             ELSE
               IF(UE(I) <= Q(I)) THEN
                 Ovlap(I) = UE(I) - LE(I)
               ELSE
                 Ovlap(I) = Q(I) - LE(I)
               ENDIF
             ENDIF
           ENDDO
           OverlapV = Ovlap(1)*Ovlap(2)*Ovlap(3)
           IF(OverlapV <= 0.0D0) STOP 'ERROR: OverlapV is negative!'
           CubeV = (UE(1)-LE(1))*(UE(2)-LE(2))*(UE(3)-LE(3))
           RecurCalLeavesTm = (Cube%LayGridCost)*OverlapV/CubeV
         ENDIF
#endif
       ELSE
         RecurCalLeavesTm = RecurCalLeavesTm(Cube%Descend)+&
                            RecurCalLeavesTm(Cube%Descend%Travrse)
       ENDIF
     END FUNCTION RecurCalLeavesTm

!===============================================================================
  RECURSIVE FUNCTION LeavesTmCount(Cube)
    TYPE(CubeNode), POINTER    :: Cube
    REAL(DOUBLE)               :: LeavesTmCount
    IF(Cube%Leaf)THEN
      LeavesTmCount=Cube%LayGridCost
    ELSE
      LeavesTmCount=LeavesTmCount(Cube%Descend)+LeavesTmCount(Cube%Descend%Travrse)
    ENDIF
  END FUNCTION LeavesTmCount

!===============================================================================
  SUBROUTINE CollectTime()
    INTEGER::IErr,I
    INTEGER,DIMENSION(MPI_STATUS_SIZE)::Status

    CALL New(TmNode,NPrc)
    IF(MyID == 0) THEN
      TmNode%D(1) = VolTm
      DO I = 1, NPrc-1
        CALL MPI_Recv(TmNode%D(I+1),1,MPI_DOUBLE_PRECISION,I,Time_T,MPI_COMM_WORLD,Status,IErr)
      ENDDO
      DO I = 0, NPrc-1
        WRITE(*,*) 'Time on node = '//TRIM(IntToChar(I))//' is ',TmNode%D(I+1)
      ENDDO
    ELSE
      CALL MPI_Send(VolTm,1,MPI_DOUBLE_PRECISION,0,Time_T,MPI_COMM_WORLD,Status,IErr)
    ENDIF
    
  END SUBROUTINE CollectTime
!===============================================================================
  SUBROUTINE CalImbalance
    REAL(DOUBLE)::TmMax,Imbalance,DevSum
    INTEGER::I
  
    TmMax = -100.0D0
    IF(MyID == 0) THEN
      DO I = 1, NPrc
        TmMax = Max(TmMax,TmNode%D(I))
      ENDDO
      WRITE(*,*) 'TmMax = ',TmMax
      DevSum = 0.0D0
      DO I = 1, NPrc
        DevSum = DevSum + (TmMax-TmNode%D(I))
      ENDDO
      Imbalance = DevSum/(NPrc*TmMax)
      WRITE(*,*) 'Imbalance = ',Imbalance
    ENDIF
  END SUBROUTINE CalImbalance

!===============================================================================
  FUNCTION AtomRad(Z)
    REAL(DOUBLE)::Z,AtomRad,TableAtomRad(200)

    !! Values are obtained from Slater, JCP 41 (1964) 3199
    TableAtomRad(1) = 0.25*AngstromsToAU
    TableAtomRad(8) = 0.60*AngstromsToAU
    IF(ABS(Z-1.0D0) > 1.0E-10 .AND. ABS(Z-8.0D0) > 1.0E-10) THEN
      WRITE(*,*) 'ERROR: Atom radius for this Z has not been set yet!'
      STOP
    ENDIF
    AtomRad = TableAtomRad(NINT(Z))
  END FUNCTION AtomRad

!===============================================================================
  SUBROUTINE Get_SubVol(LCoor,RCoor,N)
    TYPE(DBL_RNK2)::LCoor,RCoor
    INTEGER::N,I,CIndex,J,K,SmallN,PIndex,Stage,DirInt
    INTEGER::Power2(0:31)
    REAL(DOUBLE)::OvLen(3),XEff(3),Sum1,Sum2,CenterZ,&
      ZEff,MaxDim,LinDim,Rad,RangeL,RangeR,BoxL,BoxR
    TYPE(DBL_VECT)::DLCoor,DRCoor

    IF(MyID == 0) THEN
      WRITE(*,*) 'ThreeDFracZ_p2_partition is implemented...'
      CALL OpenASCII(OutFile,Out)
      WRITE(Out,*) 'ThreeDFracZ_p2_partition is implemented...'
      CLOSE(Out,STATUS='KEEP')
    ENDIF
    Power2(0) = 1
    DO I = 1,31
      Power2(I) = Power2(I-1)*2
    ENDDO
    ! check if N is a power of 2
    SmallN = Nint(Log(N*1.0)/Log(2.0D0))
    IF(N /= Power2(SmallN)) THEN
      WRITE(*,*) 'ERROR: N is ', N, ', N is not a power of 2, stop!!!'
      STOP
    ENDIF
  
    CALL New(DLCoor,3)
    CALL New(DRCoor,3)
    DO Stage = 1, SmallN
      DO PIndex = 1, Power2(Stage-1)
        ! WRITE(*,*) 'Stage = ', Stage, ', PIndex = ',PIndex
        !! CIndex is the index of the child subvolume
        CIndex = PIndex + Power2(Stage-1)
        !! Determine the direction to split
        MaxDim = -1.0D0
        DO J = 1, 3
          LinDim = RCoor%D(J,PIndex)-LCoor%D(J,PIndex)
          IF(LinDim > MaxDim) THEN
            DirInt = J
            MaxDim = LinDim
          ENDIF
        ENDDO
        !! The direction is DirInt

        !! now go through all atoms and check if each of them
        !! is inside an expanded volume.
        Sum1 = 0.0D0
        Sum2 = 0.0D0
        ! WRITE(*,*) 'NAtoms is ', NAtoms
        DO J = 1, NAtoms
          ! Rad = AtomRad(GM%AtNum%I(J))
          Rad = AtomRad(GM%AtNum%D(J))
          ! WRITE(*, *) 'Rad = ', Rad
          !! expand the volume 
          DLCoor%D(1:3) = LCoor%D(1:3,PIndex) - Rad
          DRCoor%D(1:3) = RCoor%D(1:3,PIndex) + Rad
          IF(GM%Carts%D(1,J) > DLCoor%D(1) .AND. &
             GM%Carts%D(2,J) > DLCoor%D(2) .AND. &
             GM%Carts%D(3,J) > DLCoor%D(3) .AND. &
             GM%Carts%D(1,J) < DRCoor%D(1) .AND. &
             GM%Carts%D(2,J) < DRCoor%D(2) .AND. &
             GM%Carts%D(3,J) < DRCoor%D(3)) THEN
            DO K = 1, 3
              RangeL = GM%Carts%D(K,J)-Rad
              RangeR = GM%Carts%D(K,J)+Rad
              BoxL = LCoor%D(K,PIndex)
              BoxR = RCoor%D(K,PIndex)
              !! now check 11 cases
              IF(RangeL < BoxL .AND. RangeR <= BoxL) THEN
                WRITE(*,*) 'ERROR: case 1 not possible!'
                STOP
              ELSE IF(RangeL<BoxL .AND. RangeR<=BoxR) THEN 
                OvLen(K) = RangeR - BoxL
                XEff(K) = BoxL+OvLen(K)/2.0D0
              ELSE IF(RangeL<BoxL .AND. RangeR>BoxR) THEN
                OvLen(K) = BoxR-BoxL
                XEff(K) = BoxL+OvLen(K)/2.0D0
              ELSE IF(RangeL==BoxL .AND. RangeR<=BoxR) THEN
                OvLen(K) = RangeR-BoxL
                XEff(K) = BoxL+OvLen(K)/2.0D0
              ELSE IF(RangeL==BoxL .AND. RangeR>BoxR) THEN
                OvLen(K) = BoxR-BoxL
                XEff(K) = BoxL+OvLen(K)/2.0D0
              ELSE IF(RangeL<BoxR .AND. RangeR<=BoxR) THEN
                OvLen(K) = RangeR-RangeL
                XEff(K) = RangeL+OvLen(K)/2.0D0
              ELSE IF(RangeL<BoxR .AND. RangeR>BoxR) THEN
                OvLen(K) = BoxR-RangeL
                XEff(K) = RangeL+OvLen(K)/2.0D0
              ELSE
                WRITE(*,*) 'ERROR: A missing case ??'
                STOP
              ENDIF
            ENDDO
            !! 4.18879D0 is 4Pi/3
            ZEff = GM%AtNum%D(J)*OvLen(1)*OvLen(2)*OvLen(3)/(4.18879D0*Rad**3.0)
            ! WRITE(*,*) 'ZEff: GM_AtNum_D = ',GM%AtNum%D(J)
            Sum1 = Sum1 + ZEff*XEff(DirInt)
            Sum2 = Sum2 + ZEff
          ENDIF
        ENDDO
        !! check the condition
        IF(Sum2 == 0.0D0) THEN
          STOP 'ERROR: The sum of effective charge in vol PIndex is zero!'
        ENDIF
        CenterZ = Sum1/Sum2
        !! split the box
        RCoor%D(1:3,CIndex) = RCoor%D(1:3,PIndex)
        LCoor%D(1:3,CIndex) = LCoor%D(1:3,PIndex)
        RCoor%D(DirInt,PIndex) = CenterZ
        LCoor%D(DirInt,CIndex) = CenterZ
      ENDDO
    ENDDO
    CALL Delete(DLCoor)
    CALL Delete(DRCoor)
    
  END SUBROUTINE Get_SubVol
#endif
END MODULE
