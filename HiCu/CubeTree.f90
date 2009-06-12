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

#undef USE_LEAF_COUNT

MODULE CubeTree
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
  USE MondoLogger

  IMPLICIT NONE
  !====================================================================================
  !  Hierarchical cubature node
  !====================================================================================
  TYPE CubeNode
#ifdef PARALLEL
    REAL(DOUBLE)                            :: LayGridCost
#endif
    LOGICAL                                 :: Leaf
    !     Intermediate values
    INTEGER                                 :: ISplit
    REAL(DOUBLE),DIMENSION(2)               :: ICube
    REAL(DOUBLE)                            :: ECube
    !     Bounding box
    TYPE(BBox)                              :: Box
    !     Links
    TYPE(CubeNode), POINTER                 :: Travrse  ! Next node in tree traversal
    TYPE(CubeNode), POINTER                 :: Descend  ! Next node in tree descent
    !     Cubature grid
    REAL(DOUBLE),   POINTER, DIMENSION(:,:) :: Grid     ! Transformed grid
    REAL(DOUBLE),   POINTER, DIMENSION(:,:) :: Vals     ! Values at each grid pt
    REAL(DOUBLE),   POINTER, DIMENSION(:)   :: Wght     ! Transformed Wght
  END TYPE CubeNode
  !----------------------------------------------------------------------------------
  !  Global variables for statistics accumulation
  INTEGER, PARAMETER             :: MaxTier=200     !  Maximum recursion
  INTEGER                        :: MaxLevel        !  Current maximum recursion depth
  INTEGER, DIMENSION(0:MaxTier)  :: GlobalCubes     !  Number of cubes per tier
  INTEGER                        :: CubeNodes       !  Global Cube counter
  REAL(DOUBLE)                   :: Exc             !  Exchange correlation energy
  REAL(DOUBLE),PARAMETER         :: ResSpan=1.D4    !  Controls span of accuracy
  INTEGER,     PARAMETER         :: MaxRes=12       !  Max intervals to divide accuracy over
  !  Global variables for current Cube
  TYPE(BBox)                     :: Box             !  Global BBox, set to current Cube%BBox
  REAL(DOUBLE),DIMENSION(NGrid,3):: Grid            !  Global Grid, set to current Cube%Grid
  REAL(DOUBLE),DIMENSION(NGrid * 3,4):: RhoV            !  Global Vals, set to current Cube%Vals  !<<< SPIN
  REAL(DOUBLE)                   :: Pop,dPopX,dPopY,dPopZ  !  Global Pop
  !  Global Objects
  TYPE(BSET)                     :: BS              !  Global basis set
  TYPE(CRDS)                     :: GM              !  Global molecular geometry
  TYPE(CubeNode), POINTER        :: CubeRoot

  INTEGER::NSDen
  !-----------!
CONTAINS !
  !================================================================================
  !     Grid generation routine
  !================================================================================
  SUBROUTINE GridGen(WBox,SubVolRho,SubVolExc)
    TYPE(BBox)                       :: WBox
    REAL(DOUBLE)                     :: SubVolRho,SubVolExc
    REAL(DOUBLE),   DIMENSION(2)     :: TotalError,LocalError,GlobalError, &
         RelativeError
    REAL(DOUBLE),   DIMENSION(2)     :: NewCubes,OldCubes
    REAL(DOUBLE)                     :: MaxError,BoxSep,Delta,TargtThresh,IXact
    REAL(DOUBLE), DIMENSION(0:MaxRes):: MaxRelError
    INTEGER                          :: I,J,K,NRes,ErrCount,PtsPerAtom,PU
    TYPE(BBox)                       :: CubeBox
    CHARACTER(LEN=DEFAULT_CHR_LEN)   :: Mssg
    INTEGER                            :: NC
    REAL(DOUBLE), DIMENSION(3)         :: BoxCenter,BoxBndLow,BoxBndHig
    !---------------------------------------------------------------------------------
    !        Initialize some global variables
    GlobalCubes=0
    GlobalError=Zero
    MaxLevel=0
    CubeNodes=0
    !        Initialize the CubeRoot and set thresholding
    CALL InitCubeRoot(CubeRoot,WBox)

    !        CubeRoot%Box%BndBox(1:3,1:2) = WBox%BndBox(1:3,1:2)
    !        CubeRoot%Box%Center(1:3) = (CubeRoot%Box%BndBox(1:3,1)+CubeRoot%Box%BndBox(1:3,2))*Half
    !        CubeRoot%Box%Half(1:3) = (CubeRoot%Box%BndBox(1:3,2)-CubeRoot%Box%BndBox(1:3,1))*Half
    !        Compute total electron population in this box
    CALL SetBBox(CubeRoot%Box,Box)
    IXact=Zero
    BoxCenter(:) = Box%Center(:)
    BoxBndLow(:) = Box%BndBox(:,1)
    BoxBndHig(:) = Box%BndBox(:,2)
    DO NC = 1,CS_OUT%NCells
      Box%Center(:)   = BoxCenter(:)+CS_OUT%CellCarts%D(:,NC)
      Box%BndBox(:,1) = BoxBndLow(:)+CS_OUT%CellCarts%D(:,NC)
      Box%BndBox(:,2) = BoxBndHig(:)+CS_OUT%CellCarts%D(:,NC)
      IXact=IXact+PopInBox(RhoRoot)
    ENDDO

    Box%Center(:)   = BoxCenter(:)
    Box%BndBox(:,1) = BoxBndLow(:)
    Box%BndBox(:,2) = BoxBndHig(:)
    !        Seting parameters for grid generation with variable accuracy
    MaxRelError=BIG_DBL
    Delta=1.D-1
    NRes=LOG10(1.D-1/TauRel)
    TauRel=1.D-1
    PU=OpenPU(); CALL PrintProtectL(PU); CLOSE(PU)
    !        Begin generation of the hierarchical grid
    DO J=1,NRes
      TauRel=TauRel*Delta          !  Current target error
      CALL GridRefine(CubeRoot)    !  Refine grid for this value of TauRel
      NewCubes=CubeWalk(CubeRoot)  !  Compute Exc and convergence parameters
      !           Compute and print convergence statistics
      PtsPerAtom=INT(DBLE(NGrid*LeafCount(CubeRoot))/DBLE(NAtoms))
      RelativeError(1)=ABS(IXact-NewCubes(1))/IXact
      Exc=NewCubes(2)
#ifdef PARALLEL
      !This breaks the load balancing algorithm!
      !Mssg=ProcessName('HiCu.GridGen')                      &
      !     //'Tau = ' //TRIM(DblToShrtChar(TauRel))         &
      !     //', <Rho> = '//TRIM(DblToMedmChar(Reduce(NewCubes(1)))) &
      !     //', <Exc> = '//TRIM(DblToMedmChar(Reduce(Exc)))         &
      !     //', Pts/Atom = '//TRIM(IntToChar(Reduce(PtsPerAtom)))
      !IF(MyID.EQ.0)THEN
      !   CALL OpenASCII(OutFile,Out)
      !   WRITE(Out,*)TRIM(Mssg)
      !   CLOSE(Out)
      !ENDIF
#else
      CALL MondoLog(DEBUG_MAXIMUM, 'HiCu', &
           'Tau = ' //TRIM(DblToShrtChar(TauRel))         &
           //', <Rho> = '//TRIM(DblToMedmChar(NewCubes(1))) &
           //', <Exc> = '//TRIM(DblToMedmChar(Exc))         &
           //', Pts/Atom = '//TRIM(IntToChar(PtsPerAtom)), "GridGen")
#endif
      OldCubes=NewCubes
    ENDDO
    SubVolRho = NewCubes(1)
    SubVolExc = Exc
#ifdef PARALLEL
#else
    CALL MondoLog(DEBUG_MEDIUM, 'HiCu', &
         'TauRel = '//TRIM(DblToShrtChar(TauRel))              &
         //', RhoErr = '//TRIM(DblToShrtChar(RelativeError(1)))  &
         //', <Exc> = '//TRIM(DblToMedmChar(Exc)), "GridGen")
#endif
    PU=OpenPU(); CALL PrintProtectR(PU); CLOSE(PU)
  END SUBROUTINE GridGen
  !================================================================================
  !
  !================================================================================
  RECURSIVE SUBROUTINE GridRefine(Cube)
    TYPE(CubeNode), POINTER :: Cube
    !------------------------------------------------------------------
    IF(Cube%Leaf)THEN
      IF(Cube%ECube>TauRel)THEN
        CALL SplitCube(Cube)
        CALL GridRefine(Cube%Descend)
        CALL GridRefine(Cube%Descend%Travrse)
      ENDIF
    ELSE
      CALL GridRefine(Cube%Descend)
      CALL GridRefine(Cube%Descend%Travrse)
    ENDIF
  END SUBROUTINE GridRefine
  !==============================================================================================
  !
  !==============================================================================================
  SUBROUTINE SplitCube(Cube)
    TYPE(CubeNode), POINTER  :: Cube
    TYPE(CubeNode), POINTER  :: Left,Right
    REAL(DOUBLE),DIMENSION(3):: DD,MaxVar
    REAL(DOUBLE)             :: Q,MaxDir
    INTEGER                  :: J,ISplit
    REAL(DOUBLE)             :: StartWTime,EndWTime
#ifdef PARALLEL
    REAL(DOUBLE),EXTERNAL    :: MondoTimer
#endif
    !------------------------------------------------------------------------
    !        Stupid should be painfull...
    IF(.NOT.Cube%Leaf)CALL Halt(' Logic error in SplitCube ')
    Cube%Leaf=.FALSE.
    !        Determine direction to split
    ISplit=MaxSplit(Cube)
    !        Free grid memory
    CALL DeleteCubeGrid(Cube)
    !        Allocate new cubes
    CALL NewCubeNode(Cube%Descend,Cube%Box%Tier+1)
    CALL NewCubeNode(Cube%Descend%Travrse,Cube%Box%Tier+1)
    !        Set links
    Cube%Descend%Travrse%Travrse=>Cube%Travrse
    Left =>Cube%Descend
    Right=>Cube%Descend%Travrse
    CALL SplitBox(Cube%Box,Left%Box,Right%Box,ISplit)
    !        Compute approximate and exact integrals of the density
#ifdef PARALLEL
    StartWTime = MondoTimer()
#endif
    CALL LayGrid(Left)
#ifdef PARALLEL
    EndWTime = MondoTimer()
    Left%LayGridCost = EndWTime - StartWTime
#ifdef USE_LEAF_COUNT
    Left%LayGridCost = 1.0D0
#endif
    StartWTime = MondoTimer()
#endif
    CALL LayGrid(Right)
#ifdef PARALLEL
    EndWTime = MondoTimer()
    Right%LayGridCost = EndWTime - StartWTime
#ifdef USE_LEAF_COUNT
    Right%LayGridCost = 1.0D0
#endif
#endif
    !        Compute the exact cubature error for the density
  END SUBROUTINE SplitCube
#ifdef NEWSPLIT
  !===============================================================================
  !    New max split that tries to follow the largest integration error
  !    if the cube aspect ratio isn't too bad
  !===============================================================================
  FUNCTION MaxSplit(Cube) RESULT(ISplit)
    TYPE(CubeNode), POINTER   :: Cube
    REAL(DOUBLE)              :: MaxDim,MinDim,LinDim
    REAL(DOUBLE),DIMENSION(3) :: DD,MaxVar
    REAL(DOUBLE)              :: Q,MaxDir
    INTEGER                   :: I,J,ISplit
    MaxDim = -1.0D0
    MinDim =  1D10
    DO I = 1, 3
      LinDim = Cube%Box%BndBox(I,2)-Cube%Box%BndBox(I,1)
      MinDim=MIN(MinDim,LinDim)
      IF(LinDim > MaxDim) THEN
        MaxDim = LinDim
        ISplit = I
      ENDIF
    ENDDO
    IF(MaxDim/MinDim>Three)RETURN
    DD=Zero
    DO J=1,NGrid
      DD(1)=MAX(DD(1),ABS(Cube%Vals(J,3)))
      DD(2)=MAX(DD(2),ABS(Cube%Vals(J,4)))
      DD(3)=MAX(DD(3),ABS(Cube%Vals(J,5)))
    ENDDO
    MaxDir=Zero
    DO J=1,3
      MaxVar(J)=DD(J)*Cube%Box%Half(J)**2
      MaxDir=MAX(MaxDir,MaxVar(J))
    ENDDO
    DO J=1,3
      IF(ABS(MaxDir-MaxVar(J))<1.D-8)THEN
        ISplit=J
        EXIT
      ENDIF
    ENDDO
  END FUNCTION MaxSplit
#else
  !===============================================================================
  !    Chee Kwans new max split that splits the largest dimension
  !===============================================================================
  FUNCTION MaxSplit(Cube) RESULT(ISplit)
    TYPE(CubeNode), POINTER   :: Cube
    REAL(DOUBLE)              :: MaxDim,LinDim
    INTEGER                   :: ISplit,I
    MaxDim = -1.0D0
    DO I = 1, 3
      LinDim = Cube%Box%BndBox(I,2)-Cube%Box%BndBox(I,1)
      IF(LinDim > MaxDim) THEN
        MaxDim = LinDim
        ISplit = I
      ENDIF
    ENDDO
  END FUNCTION MaxSplit
#endif
  !===============================================================================
  !     Lay the density out on the cubes grid
  !===============================================================================
  SUBROUTINE LayGrid(Cube)
    TYPE(CubeNode), POINTER            :: Cube
    REAL(DOUBLE),   DIMENSION(NGrid * 3):: Rho,AbsGradRho2         !<<< SPIN
    REAL(DOUBLE),   DIMENSION(NGrid * 3):: E,dEdRho,dEdAbsGradRho2 !<<< SPIN
    REAL(DOUBLE),   DIMENSION(NGrid * 6):: Buf
    REAL(DOUBLE),   DIMENSION(3)       :: GradRhoOnTheCube
    REAL(DOUBLE)                       :: EOnTheCube,dEd1OnTheCube, &
         dEd2OnTheCube,RhoOnTheCube, &
         DeltaRho,DeltaGrad
    INTEGER                            :: I,J,K,iSDen,iGBeg,iGEnd,I1,I2
    INTEGER                            :: NC
    REAL(DOUBLE), DIMENSION(3)         :: BoxCenter,BoxBndLow,BoxBndHig
    REAL(DOUBLE), DIMENSION(NGRID,3)   :: GridOld
    REAL(DOUBLE)                       :: Rsum,PopOld
    !--------------------------------------------------------------------------
    !        Transform cubature rule to this nodes bounding box
    CALL CubeRule(Cube)
    !--------------------------------------------------------------------------
    !        Initialize global gridding variables
    Pop=Zero
    dPopX=Zero
    dPopY=Zero
    dPopZ=Zero
    RhoV=Zero
    Grid=Cube%Grid
    CALL SetBBox(Cube%Box,Box)
    !-------------------------------------------------------------------------
    !        Lay the grid
    BoxCenter(:) = Box%Center(:)
    BoxBndLow(:) = Box%BndBox(:,1)
    BoxBndHig(:) = Box%BndBox(:,2)
    GridOld(:,:) = Grid(:,:)
    DO NC = 1,CS_OUT%NCells
      Box%Center(:)   = BoxCenter(:)+CS_OUT%CellCarts%D(:,NC)
      Box%BndBox(:,1) = BoxBndLow(:)+CS_OUT%CellCarts%D(:,NC)
      Box%BndBox(:,2) = BoxBndHig(:)+CS_OUT%CellCarts%D(:,NC)
      DO I=1,NGrid
        Grid(I,:) = GridOld(I,:)+CS_OUT%CellCarts%D(:,NC)
      ENDDO
      CALL RhoOnGrid(RhoRoot)
    ENDDO
    Box%Center(:)   = BoxCenter(:)
    Box%BndBox(:,1) = BoxBndLow(:)
    Box%BndBox(:,2) = BoxBndHig(:)
    Grid(:,:)       = GridOld(:,:)
    !----------------------------------------------------------------------------------
    ! Rho(2*NGrid)
    ! AbsGradRho2(3*NGrid)
    ! E(NGrid)
    ! dEdRho(2*NGrid)
    ! dEdAbsGradRho2(3*NGrid)
    ! Transfer global Cube values to local Cube
    IF(NSDen.EQ.1)THEN
      DO I=1,NGrid
        Rho(I)=RhoV(I,1)!Rho_tot
        AbsGradRho2(I)=RhoV(I,2)**2+RhoV(I,3)**2+RhoV(I,4)**2!DOT(GRho_tot,GRho_tot)
        Cube%Vals(I,3)=RhoV(I,2)!GRho_tot_x
        Cube%Vals(I,4)=RhoV(I,3)!GRho_tot_y
        Cube%Vals(I,5)=RhoV(I,4)!GRho_tot_z
      ENDDO
    ELSEIF(NSDen.EQ.3)THEN
      DO I=1,NGrid
        I1=I+  NGrid
        I2=I+2*NGrid
        Rho(I )=RhoV(I1,1)!Rho_a
        Rho(I1)=RhoV(I2,1)!Rho_b
        AbsGradRho2(I )=RhoV(I1,2)**2+RhoV(I1,3)**2+RhoV(I1,4)**2!DOT(GRho_a,GRho_a)
        AbsGradRho2(I1)=RhoV(I2,2)**2+RhoV(I2,3)**2+RhoV(I2,4)**2!DOT(GRho_b,GRho_b)
        AbsGradRho2(I2)=RhoV(I1,2)*RhoV(I2,2)+RhoV(I1,3)*RhoV(I2,3)+RhoV(I1,4)*RhoV(I2,4)!DOT(GRho_a,GRho_b)
        Cube%Vals(I ,3)=RhoV(I1,2)!GRho_a_x
        Cube%Vals(I ,4)=RhoV(I1,3)!GRho_a_y
        Cube%Vals(I ,5)=RhoV(I1,4)!GRho_a_z
        Cube%Vals(I1,3)=RhoV(I2,2)!GRho_b_x
        Cube%Vals(I1,4)=RhoV(I2,3)!GRho_b_y
        Cube%Vals(I1,5)=RhoV(I2,4)!GRho_b_z
      ENDDO
    ELSE
      CALL Halt('LayGrid: Wrong NSDen')
    ENDIF
    ! Evaluate Exc, dExcdRho, and dExcdAbsGradRho2 on the grid
    CALL ExcOnTheGrid_ClSh(NGrid,Rho(1),AbsGradRho2(1),E(1),dEdRho(1),dEdAbsGradRho2(1),Buf(1),NSDen)
    ! Transfer global values to the cube and compute approximate integrals
    ! of E, dE/dRho, dE/d(AbsGradRho^2), Rho and GradRho over the cube
    EOnTheCube=Zero
    dEd1OnTheCube=Zero
    dEd2OnTheCube=Zero
    RhoOnTheCube=Zero
    GradRhoOnTheCube=Zero

    IF(NSDen.EQ.1)THEN
      DO I=1,NGrid
        Cube%Vals(I,1)=dEdRho(I)
        Cube%Vals(I,2)=dEdAbsGradRho2(I)
        EOnTheCube=EOnTheCube+Cube%Wght(I)*E(I)                   ! E on the cube
        dEd1OnTheCube=dEd1OnTheCube+Cube%Wght(I)*dEdRho(I)        ! dEdRho on the cube
        dEd2OnTheCube=dEd2OnTheCube+Cube%Wght(I)*dEdAbsGradRho2(I)! dEdAbsGradRho on the cube
        RhoOnTheCube=RhoOnTheCube+Cube%Wght(I)*RhoV(I,1)          ! Rho on the cube
        GradRhoOnTheCube(1:3)=GradRhoOnTheCube(1:3) &             ! GradRho on the cube
             +Cube%Wght(I)*RhoV(I,2:4)
      ENDDO
    ELSEIF(NSDen.EQ.3)THEN
      CALL DCOPY(2*NGrid,dEdRho(1)        ,1,Cube%Vals(1,1),1)
      CALL DCOPY(3*NGrid,dEdAbsGradRho2(1),1,Cube%Vals(1,2),1)
      DO I=1,NGrid
        EOnTheCube=EOnTheCube+Cube%Wght(I)*E(I)                   ! E on the cube
        dEd1OnTheCube=dEd1OnTheCube+Cube%Wght(I)* &
             (dEdRho(I)+dEdRho(I+NGrid))                 ! dEdRho on the cube
        dEd2OnTheCube=dEd2OnTheCube+Cube%Wght(I)* &
             (dEdAbsGradRho2(I)+dEdAbsGradRho2(I+NGrid)) ! dEdAbsGradRho on the cube
        RhoOnTheCube=RhoOnTheCube+Cube%Wght(I)*RhoV(I,1)          ! Rho on the cube
        GradRhoOnTheCube(1:3)=GradRhoOnTheCube(1:3) &             ! GradRho on the cube
             +Cube%Wght(I)*RhoV(I,2:4)
      ENDDO
    ELSE
      CALL Halt('LayGrid: Wrong NSDen')
    ENDIF
    !----------------------------------------------------------------------------------
    Cube%ICube(1)=RhoOnTheCube
    Cube%ICube(2)=EOnTheCube
    !        Compute local error estimates
    DeltaRho=ABS(Pop-RhoOnTheCube)
    DeltaGrad=(dPopX-GradRhoOnTheCube(1))**2 &
         +(dPopY-GradRhoOnTheCube(2))**2 &
         +(dPopZ-GradRhoOnTheCube(3))**2
    Cube%ECube=DeltaRho+SQRT(DeltaGrad)/Three
    !        Determine optimal direction for spliting
  END SUBROUTINE LayGrid
  !=================================================================================
  !     Sums the significant contributions (leaves) to the density at NGrid points
  !     in a Cube, and evaluates its exact contribution to the total electron count
  !=================================================================================
  RECURSIVE SUBROUTINE RhoOnGrid(Node)
    TYPE(RhoNode), POINTER                     :: Node
    REAL(DOUBLE)                               :: Tx,Ty,Tz
    REAL(DOUBLE), DIMENSION(0:HGEll+1)         :: LLambdaX,LLambdaY,LLambdaZ, &
         ULambdaX,ULambdaY,ULambdaZ, &
         LambdaX,LambdaY,LambdaZ
    REAL(DOUBLE)                               :: RQx,RQy,RQz,RQ2,Z,X,W,Sgn,Xpt,Co,        &
         LQx,LQy,LQz,LQ2,UQx,UQy,UQz,UQ2,         &
         LXpt,UXpt,TwoZ,SqZ,CoFact,RL1,TmpX,TmpY, &
         LXptX,LXptY,LXptZ,UXptX,UXptY,UXptZ,RL2
    INTEGER                                    :: I,J,IQ,IC,JQ,JC,KQ,KC,L,Ell,L1,L2,M,N,LMN,IGrid
    INTEGER :: OffSDen,iSDen,I0
    !-------------------------------------------------------------------------------------------------------
    Tx=ABS(Node%Box%Center(1)-Box%Center(1))
    IF(Tx>Node%Box%Half(1)+Box%Half(1))RETURN
    Ty=ABS(Node%Box%Center(2)-Box%Center(2))
    IF(Ty>Node%Box%Half(2)+Box%Half(2))RETURN
    Tz=ABS(Node%Box%Center(3)-Box%Center(3))
    IF(Tz>Node%Box%Half(3)+Box%Half(3))RETURN
    !
    IF(Node%Leaf)THEN
      OffSDen=LHGTF(Node%Ell)
#ifdef EXPLICIT_SOURCE
      INCLUDE 'ExplicitLeafPopulation.Inc'
      INCLUDE 'ExplicitLeafContribution.Inc'
#else
      INCLUDE 'GeneralLeafPopulation.Inc'
      INCLUDE 'GeneralLeafContribution.Inc'
#endif
    ELSE
      CALL RhoOnGrid(Node%Descend)
      CALL RhoOnGrid(Node%Descend%Travrse)
    ENDIF
  END SUBROUTINE RhoOnGrid
  !=================================================================================
  !     Sums the significant contributions (leaves) to the density at NGrid points
  !     in a Cube, and evaluates its exact contribution to the total electron count
  !=================================================================================
  RECURSIVE FUNCTION PopInBox(Node) RESULT(EPop)
    TYPE(RhoNode), POINTER                     :: Node
    REAL(DOUBLE)                               :: Tx,Ty,Tz
    REAL(DOUBLE), DIMENSION(0:HGEll+1)         :: LLambdaX,LLambdaY,LLambdaZ, &
         ULambdaX,ULambdaY,ULambdaZ, &
         LambdaX,LambdaY,LambdaZ
    REAL(DOUBLE)                               :: RQx,RQy,RQz,RQ2,Z,X,W,Sgn,Xpt,Co,              &
         LQx,LQy,LQz,LQ2,UQx,UQy,UQz,UQ2,         &
         LXpt,UXpt,TwoZ,SqZ,CoFact,RL1,TmpX,TmpY, &
         LXptX,LXptY,LXptZ,UXptX,UXptY,UXptZ,EPop,RL2
    INTEGER                                    :: I,J,IQ,IC,JQ,JC,KQ,KC,L,Ell,L1,L2,M,N,LMN,GKount
    !INTEGER::OffSDen
    !-------------------------------------------------------------------------------------------------------
    EPop=Zero
    Tx=ABS(Node%Box%Center(1)-Box%Center(1))
    IF(Tx>Node%Box%Half(1)+Box%Half(1))RETURN
    Ty=ABS(Node%Box%Center(2)-Box%Center(2))
    IF(Ty>Node%Box%Half(2)+Box%Half(2))RETURN
    Tz=ABS(Node%Box%Center(3)-Box%Center(3))
    IF(Tz>Node%Box%Half(3)+Box%Half(3))RETURN
    IF(Node%Leaf)THEN
      !           Intermediates for computation and thresholding of electron count contributions
      Pop=Zero
      !OffSDen=0*LHGTF(Node%Ell)!if you use that be careful with RhoOnGrid!
#ifdef EXPLICIT_SOURCE
      INCLUDE 'ExplicitLeafPopulation.Inc'
#else
      INCLUDE 'GeneralLeafPopulation.Inc'
#endif
      EPop=Pop
    ELSE
      EPop=PopInBox(Node%Descend) &
           +PopInBox(Node%Descend%Travrse)
    ENDIF
  END FUNCTION PopInBox
  !=====================================================================
  !     Generate a cubature rule for the bounding box, performing
  !     affine transformations and possible non-linear coordinate
  !     mappings to improve convergence
  !=====================================================================
  SUBROUTINE CubeRule(Node)
    TYPE(CubeNode), POINTER          :: Node
    REAL(DOUBLE)                     :: Shift,Slope
    INTEGER                          :: I,J
    !----------------------------------------------------------------------
    !        Transform from the [-1,1]x[-1,1]x[-1,1] rule to a new box
    DO J=1,NGrid
      Node%Wght(J)=CubeRuleWght(J)
    ENDDO
    DO I=1,3
      Shift=Half*(Node%Box%BndBox(I,2)+Node%Box%BndBox(I,1))
      Slope=Half*(Node%Box%BndBox(I,2)-Node%Box%BndBox(I,1))
      DO J=1,NGrid
        Node%Grid(J,I)=Shift+CubeRuleGrid(I,J)*Slope
        Node%Wght(J)=Node%Wght(J)*Slope
      ENDDO
    ENDDO
  END SUBROUTINE CubeRule
  !=====================================================================
  !     Can now hopefully just use RhoTrees BBox, since this should now
  !     also change slowly and continously with geometry
  !=====================================================================
  SUBROUTINE InitCubeRoot(CubeRoot,WBox)
    TYPE(BBox)                       :: WBox
    TYPE(CubeNode),POINTER           :: CubeRoot
    REAL(DOUBLE)                     :: BoxSep,Delta,MidPop,TargetError, &
         REl,MidSep,BisSep,DelSep,FMid
    INTEGER                          :: I,J,K
    TYPE(BBox)                       :: CubeBox
    CHARACTER(LEN=DEFAULT_CHR_LEN)   :: Mssg
    !------------------------------------------------------------------------
    !        Allocate cube node
    CALL NewCubeNode(CubeRoot,0)
    !        Use the RhoTrees BBox
    !        CALL SetBBox(RhoRoot%Box,CubeRoot%Box)
    CALL SetBBox(WBox,CubeRoot%Box)
    ! moved to HiCu.F90 and XCForce.F90
    !        Set global Cube BBox
    !        CALL SetBBox(CubeRoot%Box,Box)
    CubeRoot%ISplit=1
    CubeRoot%Box%Tier=0
    CubeRoot%ECube=BIG_DBL
  END SUBROUTINE InitCubeRoot
  !==========================================================================
  !
  !==========================================================================
  SUBROUTINE NewCubeNode(Node,Level)
    TYPE(CubeNode),POINTER :: Node
    INTEGER                 :: Level
    INTEGER                 :: Status
    ALLOCATE(Node,STAT=Status)
    IF(Status/=SUCCEED) &
         CALL Halt(' ALLOCATE 1 FAILED IN NewCubeNode')
#ifdef PARALLEL
    Node%LayGridCost = 0.0D0
#endif
    Node%Box%Tier=Level
    Node%Box%Number=CubeNodes+1
    CubeNodes=CubeNodes+1
    GlobalCubes(Level)=GlobalCubes(Level)+1
    MaxLevel=MAX(MaxLevel,Level)
    Node%Leaf=.TRUE.
    Node%ECube=BIG_DBL
    NULLIFY(Node%Travrse)
    NULLIFY(Node%Descend)
    ALLOCATE(Node%Grid(NGrid,3),STAT=Status)
    CALL IncMem(Status,0,3*NGrid,'HiCu.CubeTree.Node%Grid')
    ALLOCATE(Node%Wght(NGrid),STAT=Status)
    CALL IncMem(Status,0,NGrid,'HiCu.CubeTree.Node%Wght')
    ALLOCATE(Node%Vals(NGrid*NSDen,5),STAT=Status)! <<< SPIN
    CALL IncMem(Status,0,5*NGrid,'HiCu.CubeTree.Node%Vals')
  END SUBROUTINE NewCubeNode
  !==========================================================================
  !
  !==========================================================================
  SUBROUTINE DeleteCubeGrid(Node)
    TYPE(CubeNode), POINTER  :: Node
    INTEGER                  :: Status
    DEALLOCATE(Node%Grid,STAT=Status)
    CALL DecMem(Status,0,3*NGrid)
    DEALLOCATE(Node%Wght,STAT=Status)
    CALL DecMem(Status,0,NGrid)
    DEALLOCATE(Node%Vals,STAT=Status)
    CALL DecMem(Status,0,5*NGrid)
    NULLIFY(Node%Grid)
    NULLIFY(Node%Wght)
    NULLIFY(Node%Vals)
  END SUBROUTINE DeleteCubeGrid
  !================================================================================
  !
  !================================================================================
  RECURSIVE FUNCTION CubeWalk(Cube) RESULT(ICube)
    TYPE(CubeNode), POINTER :: Cube
    REAL(DOUBLE),DIMENSION(2) :: ICube
    !--------------------------------------------------------------------------
    IF(Cube%Leaf)THEN
      ICube=Cube%ICube
    ELSE
      ICube=CubeWalk(Cube%Descend)+CubeWalk(Cube%Descend%Travrse)
    ENDIF
  END FUNCTION CubeWalk

  RECURSIVE FUNCTION ErrWalk(Cube) RESULT(ECube)
    TYPE(CubeNode), POINTER    :: Cube
    REAL(DOUBLE), DIMENSION(3) :: ECube
    !--------------------------------------------------------------------------
    IF(Cube%Leaf)THEN
      ECube=Cube%ECube
    ELSE
      ECube=ErrWalk(Cube%Descend)+ErrWalk(Cube%Descend%Travrse)
    ENDIF
  END FUNCTION ErrWalk
  !================================================================================
  !
  !================================================================================
  RECURSIVE FUNCTION LeafCount(Cube) RESULT(Count)
    TYPE(CubeNode), POINTER    :: Cube
    INTEGER                    :: Count
    INTEGER ::X
    !--------------------------------------------------------------------------
    IF(Cube%Leaf)THEN
      Count=1
    ELSE
      Count=LeafCount(Cube%Descend)+LeafCount(Cube%Descend%Travrse)
    ENDIF
  END FUNCTION LeafCount
  !================================================================================
  !================================================================================
  !
  !================================================================================
  RECURSIVE SUBROUTINE PrintCubeLeaves(Node,Level_O)
    TYPE(CubeNode), POINTER       :: Node
    INTEGER :: I
    INTEGER,DIMENSION(2),OPTIONAL :: Level_O
    IF(Node%Box%Tier==0)THEN
      CALL OpenASCII('Cubes.mma',Out,.TRUE.)
      WRITE(Out,*)'Needs["Graphics`Shapes`"];'
      !            WRITE(Out,*)'SetOptions[Graphics3D,Boxed->False]'
      WRITE(Out,*)'CubeList={}; '
      CLOSE(Out)
    ENDIF
    IF(Node%Leaf)THEN
      IF(PRESENT(Level_O))THEN
        IF(Node%Box%Tier<Level_O(1).OR. &
             Node%Box%Tier>Level_O(2))RETURN
      ENDIF
      CALL OpenASCII('Cubes.mma',Out)
      DO I=1,NGrid
        WRITE(Out,55)Node%Grid(I,1:3)

55      FORMAT('CubeList=Append[CubeList,Point[{', &
             F12.6,', ',F12.6,', ',F12.6,'}]]')
      ENDDO
      !            WRITE(Out,55)Node%Box%BndBox(1:3,1),Node%Box%BndBox(1:3,2)
      !            55 FORMAT('CubeList=Append[CubeList,Cuboid[{', &
      !                       F12.6,', ',F12.6,', ',F12.6,'},{',  &
      !                       F12.6,', ',F12.6,', '`,F12.6,'}]];')
      CLOSE(Out)
    ELSE
      CALL PrintCubeLeaves(Node%Descend,Level_O)
      CALL PrintCubeLeaves(Node%Descend%Travrse,Level_O)
    ENDIF
    IF(Node%Box%Tier==0)THEN
      CALL OpenASCII('Cubes.mma',Out)
      WRITE(Out,*)'Show[Graphics3D[{PointSize[0.001],CubeList}]];'
      !            WRITE(Out,*)'Show[WireFrame[Graphics3D[CubeList]]];'
      CLOSE(Out)
    ENDIF
  END SUBROUTINE PrintCubeLeaves
  !==========================================================================
  !
  !==========================================================================
  SUBROUTINE PrintCube(Node)
    TYPE(CubeNode), POINTER          :: Node
    INTEGER                          :: I,J
    Node%ICube=Zero
    DO J=1,NGrid
      Node%ICube=Node%ICube+Node%Wght(J)*Node%Vals(J,1)
    ENDDO
    !         WRITE(*,*)' ICube = ',Node%ICube
55  FORMAT('Rho[',D8.2,', ',D8.2,', ',D8.2,']= ',D12.6)
  END SUBROUTINE PrintCube
  !==========================================================================
  !
  !==========================================================================
  SUBROUTINE MakeBoxPeriodic(Box)
    TYPE(BBox)       :: Box
    INTEGER          :: I
    !
    DO I = 1,3
      IF(GM%PBC%AutoW%I(I)==1) THEN
        IF(Box%BndBox(I,1) < Zero) THEN
          Box%BndBox(I,1) = Zero
          Box%BndBox(I,2) = GM%PBC%BoxShape%D(I,I)
        ENDIF
        IF(Box%BndBox(I,2) > GM%PBC%BoxShape%D(I,I)) THEN
          Box%BndBox(I,1) = Zero
          Box%BndBox(I,2) = GM%PBC%BoxShape%D(I,I)
        ENDIF
      ENDIF
    ENDDO
    DO I = 1,3
      Box%Half(I)   = Half*(Box%BndBox(I,2)-Box%BndBox(I,1))
      Box%Center(I) = Half*(Box%BndBox(I,2)+Box%BndBox(I,1))
    ENDDO
  END SUBROUTINE MakeBoxPeriodic
  !==========================================================================
END MODULE CubeTree

