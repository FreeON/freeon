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
   USE AtomPairs
   IMPLICIT NONE
!====================================================================================
!  Hierarchical cubature node
!====================================================================================
   TYPE CubeNode
      LOGICAL                                 :: Leaf
!     Intermediate values
      INTEGER                                 :: ISplit
      REAL(DOUBLE)                            :: IXact
      REAL(DOUBLE),DIMENSION(3)               :: ICube    ! 
      REAL(DOUBLE),DIMENSION(3)               :: ECube    ! 
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
!  Global thresholding parameters
   REAL(DOUBLE)                   :: TauRel          !  Relative accuracy sought in the
                                                     !  total integrated electron density.
   REAL(DOUBLE)                   :: TauBox          !  Local, numerically integrated 
                                                     !  accuracy sought in electron density
   REAL(DOUBLE)                   :: TauRho          !  Local accuracy of the electron density 
                                                     !  evaluated at a point
   REAL(DOUBLE),PARAMETER         :: DeltaRel=0.2D0  !  Scales TauRel to give TauBox
   REAL(DOUBLE),PARAMETER         :: DeltaRho=1.D-2  !  Scales TauBox to give TauRho
   REAL(DOUBLE),PARAMETER         :: ResSpan=1.D4    !  Controls span of accuracy
   INTEGER,     PARAMETER         :: MaxRes=10       !  Intervals to divide accuracy over
!  Global variables for current Cube
   TYPE(BBox)                     :: Box             !  Global BBox, set to current Cube%BBox        
   REAL(DOUBLE),DIMENSION(NGrid,3):: Grid            !  Global Grid, set to current Cube%Grid
   REAL(DOUBLE),DIMENSION(NGrid,4):: RhoV            !  Global Vals, set to current Cube%Vals
   REAL(DOUBLE)                   :: Pop             !  Global Pop,  set to current Cube%IXact
!  Interpolation grid for fast evaluation of Exp[-x] and Erf[x]
   INCLUDE 'Exp.Inc'   
   INCLUDE 'Erf.Inc'
!-----------!
   CONTAINS !
!================================================================================
!     Grid generation routine     
!================================================================================
      SUBROUTINE GridGen(CubeRoot)
         TYPE(CubeNode), POINTER          :: CubeRoot
         REAL(DOUBLE),   DIMENSION(3)     :: TotalError,LocalError,GlobalError, &
                                             RelativeError,NewCubes,OldCubes
         REAL(DOUBLE)                     :: MaxError,BoxSep,Delta,TargtThresh
         REAL(DOUBLE), DIMENSION(0:MaxRes):: MaxRelError
         INTEGER                          :: I,J,K,ErrCount,PtsPerAtom
         TYPE(BBox)                       :: CubeBox    
         CHARACTER(LEN=DEFAULT_CHR_LEN)   :: Mssg     
!---------------------------------------------------------------------------------
!        Initialize some global variables
         GlobalCubes=0
         GlobalError=Zero
         MaxLevel=0       
!        Set thresholding
         TauRel=Thresholds%Cube
!        Initialize the CubeRoot and set thresholding
         CALL InitCubeRoot(CubeRoot)
!        Seting parameters for grid generation with variable accuracy
         MaxRelError=BIG_DBL
         TauRel=TauRel*ResSpan
         Delta=(One/ResSpan)**(One/DBLE(MaxRes))
!        Begin generation of the hierarchical grid
         DO J=1,MaxRes
!           Target relative error
            TauRel=TauRel*Delta
!           Accuracy of numerically integrated density in a Cube
            TauBox=TauRel*DeltaRel
!           Accuracy of density at a point           
            TauRho=TauBox*DeltaRho 
            CALL SetPenetrationThresh(TauRho,Exp_Switch)
            CALL ExpandBoxWalk(RhoRoot,One)
            CALL GridRefine(CubeRoot)
            CALL ExpandBoxWalk(RhoRoot,-One)
!           Compute Exc and convergence parameters
            NewCubes=CubeWalk(CubeRoot)
            PtsPerAtom=INT(DBLE(NGrid*LeafCount(CubeRoot))/DBLE(NAtoms))
            RelativeError(1)=ABS(CubeRoot%IXact-NewCubes(1))/CubeRoot%IXact
            Exc=NewCubes(2) 
            RelativeError(2)=ABS((NewCubes(2)-OldCubes(2))/NewCubes(2))               
            RelativeError(3)=ABS((NewCubes(3)-OldCubes(3))/NewCubes(3))               
            IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
               Mssg=ProcessName('HiCu.GridGen')                    &
                  //'TauBox = ' //TRIM(DblToShrtChar(TauBox))      &
                  //', <Rho> = '//TRIM(DblToMedmChar(NewCubes(1))) &
                  //', <Exc> = '//TRIM(DblToMedmChar(Exc))         &
                  //', Pts/Atom = '//TRIM(IntToChar(PtsPerAtom))
               CALL OpenASCII(OutFile,Out)         
               WRITE(*,*)TRIM(Mssg)
               WRITE(Out,*)TRIM(Mssg)
               CLOSE(Out)
            ENDIF
            OldCubes=NewCubes
         ENDDO
         IF(PrintFlags%Key>DEBUG_MEDIUM)THEN
            Mssg=ProcessName('HiCu.GridGen')                            &
                //'TauRel = '//TRIM(DblToShrtChar(TauRel))              &
                //', RhoErr = '//TRIM(DblToShrtChar(RelativeError(1)))  &
                //', <Exc> = '//TRIM(DblToMedmChar(Exc))         
            CALL OpenASCII(OutFile,Out)         
            WRITE(*,*)TRIM(Mssg)
            WRITE(Out,*)TRIM(Mssg)
            CLOSE(Out)
         ENDIF
     END SUBROUTINE GridGen  
!================================================================================
!          
!================================================================================
      RECURSIVE SUBROUTINE GridRefine(Cube)
         TYPE(CubeNode), POINTER :: Cube,Left,Right
         REAL(DOUBLE)            :: AbsErr,DHalf
         INTEGER                 :: ISplit
!------------------------------------------------------------------
         IF(Cube%Leaf)THEN
            IF(Cube%ECube(1)>TauBox)THEN 
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
         TYPE(CubeNode), POINTER :: Cube,Left,Right
!        Split this node
         IF(.NOT.Cube%Leaf)CALL Halt(' Logic error in SplitCube ')
         Cube%Leaf=.FALSE.
!        Free grid memory
         CALL DeleteCubeGrid(Cube)
!        Allocate new cubes
         CALL NewCubeNode(Cube%Descend,Cube%Box%Tier+1)
         CALL NewCubeNode(Cube%Descend%Travrse,Cube%Box%Tier+1)
!        Set links
         Cube%Descend%Travrse%Travrse=>Cube%Travrse
         Left =>Cube%Descend
         Right=>Cube%Descend%Travrse
         CALL SplitBox(Cube%Box,Left%Box,Right%Box,Cube%ISplit)
!        Compute approximate and exact integrals of the density
         CALL LayGrid(Left)
         CALL LayGrid(Right)
!        Compute the exact cubature error for the density 
      END SUBROUTINE SplitCube
!===============================================================================
!     Lay the density out on the cubes grid
!===============================================================================
      SUBROUTINE LayGrid(Cube)
         TYPE(CubeNode), POINTER            :: Cube
         REAL(DOUBLE),   DIMENSION(NGrid)   :: Rho,AbsGradRho2
         REAL(DOUBLE),   DIMENSION(NGrid)   :: E,dEdRho,dEdAbsGradRho2
         REAL(DOUBLE),   DIMENSION(3)       :: MaxGrad
         REAL(DOUBLE)                       :: MaxDir
         INTEGER                            :: I,J
#ifdef PERIODIC 
         INTEGER                            :: NC
         REAL(DOUBLE), DIMENSION(3)         :: BoxCenter,BoxBndLow,BoxBndHig
         REAL(DOUBLE), DIMENSION(3,NGRID)   :: GridOld
         REAL(DOUBLE)                       :: Rsum,PopOld
#endif
!--------------------------------------------------------------------------
!        Transform cubature rule to this nodes bounding box
         CALL CubeRule(Cube)
!        Lay the grid
         MaxGrad=Zero
!--------------------------------------------------------------------------
!        Initialize global griding variables
         Pop=Zero
         RhoV=Zero
         Grid=Cube%Grid
         CALL SetBBox(Cube%Box,Box)
!-------------------------------------------------------------------------
!        Lay grid
#ifdef PERIODIC
         BoxCenter(:) = Box%Center(:)
         BoxBndLow(:) = Box%BndBox(:,1)
         BoxBndHig(:) = Box%BndBox(:,2)
         GridOld(:,:) = Grid(:,:)
         DO NC = 1,CS%NCells
            Box%Center(:)   = BoxCenter(:)+CS%CellCarts%D(:,NC)
            Box%BndBox(:,1) = BoxBndLow(:)+CS%CellCarts%D(:,NC)
            Box%BndBox(:,2) = BoxBndHig(:)+CS%CellCarts%D(:,NC)        
            DO I=1,NGrid
               Grid(I,:) = GridOld(I,:)+CS%CellCarts%D(NC,:)
            ENDDO
            CALL RhoOnGrid(RhoRoot)
         ENDDO
         Box%Center(:)   = BoxCenter(:)
         Box%BndBox(:,1) = BoxBndLow(:)
         Box%BndBox(:,2) = BoxBndHig(:)  
         Grid(:,:)       = GridOld(:,:)
#else
         CALL RhoOnGrid(RhoRoot)
#endif
!        Set grid variables
         Cube%IXact=Pop
         DO I=1,NGrid             
            Rho(I)        =RhoV(I,1)
            AbsGradRho2(I)=RhoV(I,2)**2+RhoV(I,3)**2+RhoV(I,4)**2
            Cube%Vals(I,3)=RhoV(I,2)
            Cube%Vals(I,4)=RhoV(I,3)
            Cube%Vals(I,5)=RhoV(I,4)
            DO J=1,3; MaxGrad(J)=Max(MaxGrad(J),ABS(RhoV(I,J+1))); ENDDO
         ENDDO
!        Evaluate Exc, dExcdRho, and dExcdAbsGradRho2 on the grid
         CALL ExcOnTheGrid(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2) 
!        Cubature of Exc, dExcdRho, and dExcdAbsGradRho2
         Cube%ICube=Zero
         DO I=1,NGrid         
            Cube%Vals(I,1)=dEdRho(I)
            Cube%Vals(I,2)=dEdAbsGradRho2(I)
            Cube%ICube(1)=Cube%ICube(1)+Cube%Wght(I)*Rho(I)
            Cube%ICube(2)=Cube%ICube(2)+Cube%Wght(I)*E(I)
            Cube%ICube(3)=Cube%ICube(3)+Cube%Wght(I)*(dEdRho(I)*Rho(I)+dEdAbsGradRho2(I)*AbsGradRho2(I) )
         ENDDO
         Cube%ECube(1)=ABS(Cube%IXact-Cube%ICube(1)) 
!        Determine the optimal ordinate for bisection of this node
         IF(Cube%Box%Tier<10)THEN
            Cube%ISplit=MOD(Cube%Box%Tier,3)+1 
         ELSE
            MaxGrad=Zero
            DO I=1,NGrid            
               DO J=1,3; MaxGrad(J)=Max(MaxGrad(J),ABS(Cube%Vals(I,J+2))); ENDDO
            ENDDO
            MaxDir=Zero
            DO J=1,3
               MaxGrad(J)=MaxGrad(J)*Cube%Box%Half(J)**2
               MaxDir=MAX(MaxDir,MaxGrad(J))
            ENDDO
            Cube%ISplit=0       
            DO J=1,3         
               IF(ABS(MaxDir-MaxGrad(J))<1.D-8)THEN
                  Cube%ISplit=J
                  EXIT
               ENDIF
            ENDDO
            IF(Cube%ISplit==0)CALL Halt('split logic hosed ')
         ENDIF
     END SUBROUTINE LayGrid
!=================================================================================
!     Sums the significant contributions (leaves) to the density at NGrid points
!     in a Cube, and evaluates its exact contribution to the total electron count 
!=================================================================================
      RECURSIVE SUBROUTINE RhoOnGrid(Node)
         TYPE(RhoNode), POINTER                     :: Node
         REAL(DOUBLE)                               :: Tx,Ty,Tz
         REAL(DOUBLE), DIMENSION(0:MaxEll+1)        :: LLambdaX,LLambdaY,LLambdaZ, &
                                                       ULambdaX,ULambdaY,ULambdaZ, &
                                                       LambdaX,LambdaY,LambdaZ
         REAL(DOUBLE)                               :: RQx,RQy,RQz,RQ2,Z,X,W,Sgn,Xpt,Co,        &
                                                       LQx,LQy,LQz,LQ2,UQx,UQy,UQz,UQ2,         &
                                                       LXpt,UXpt,TwoZ,SqZ,CoFact,RL1,TmpX,TmpY, &
                                                       LXptX,LXptY,LXptZ,UXptX,UXptY,UXptZ, pop2
         INTEGER                                    :: I,J,IQ,IC,JQ,JC,KQ,KC,L,Ell,L1,L2,M,N,LMN,IGrid
!-------------------------------------------------------------------------------------------------------
         Tx=ABS(Node%Box%Center(1)-Box%Center(1))
         IF(Tx>Node%Box%Half(1)+Box%Half(1))RETURN
         Ty=ABS(Node%Box%Center(2)-Box%Center(2))
         IF(Ty>Node%Box%Half(2)+Box%Half(2))RETURN
         Tz=ABS(Node%Box%Center(3)-Box%Center(3))
         IF(Tz>Node%Box%Half(3)+Box%Half(3))RETURN
         IF(Node%MaxAmp+PenetratDistanceThreshold<Zero)RETURN
         IF(Node%Leaf)THEN            
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
         REAL(DOUBLE), DIMENSION(0:MaxEll+1)        :: LLambdaX,LLambdaY,LLambdaZ, &
                                                       ULambdaX,ULambdaY,ULambdaZ, &
                                                       LambdaX,LambdaY,LambdaZ
         REAL(DOUBLE)                               :: RQx,RQy,RQz,RQ2,Z,X,W,Sgn,Xpt,Co,              &
                                                       LQx,LQy,LQz,LQ2,UQx,UQy,UQz,UQ2,         &
                                                       LXpt,UXpt,TwoZ,SqZ,CoFact,RL1,TmpX,TmpY, &
                                                       LXptX,LXptY,LXptZ,UXptX,UXptY,UXptZ,EPop
         INTEGER                                    :: I,J,IQ,IC,JQ,JC,KQ,KC,L,Ell,L1,L2,M,N,LMN,GKount
!-------------------------------------------------------------------------------------------------------
         EPop=Zero
         Tx=ABS(Node%Box%Center(1)-Box%Center(1))
         IF(Tx>Node%Box%Half(1)+Box%Half(1))RETURN
         Ty=ABS(Node%Box%Center(2)-Box%Center(2))
         IF(Ty>Node%Box%Half(2)+Box%Half(2))RETURN
         Tz=ABS(Node%Box%Center(3)-Box%Center(3))
         IF(Tz>Node%Box%Half(3)+Box%Half(3))RETURN         
         IF(Node%MaxAmp+PenetratDistanceThreshold<Zero)RETURN
         IF(Node%Leaf)THEN            
!           Intermediates for computation and thresholding of electron count contributions
            Pop=Zero
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
!     Find a BBox for the CubeRoot that achieves a given percent
!     accuracy for the integrated electron density, and which changes
!     smoothly with changes in molecular geometry.
!=====================================================================
      SUBROUTINE InitCubeRoot(CubeRoot)
         TYPE(CubeNode), POINTER          :: CubeRoot
         REAL(DOUBLE)                     :: BoxSep,Delta,MidPop,TargetError, &
                                             REl,MidSep,BisSep,DelSep,FMid
         INTEGER                          :: I,J,K
         TYPE(BBox)                       :: CubeBox    
         CHARACTER(LEN=DEFAULT_CHR_LEN)   :: Mssg     
!------------------------------------------------------------------------         
!        Allocate cube node
         CALL NewCubeNode(CubeRoot,0) 
!        A minimal BBox that just encloses nuclear centers
         CALL SetBBox(RhoRoot%Box,CubeBox)
!        Exact electron count for closed shell density
         REl=Half*DBLE(NEl)
!        Recursive bisection to determine largest BBox for CubeRoot that integrates
!        the density to within TauRel
         BisSep=Zero
         DelSep=20.0D0
         DO J=1,100
!           Half the step size
            DelSep=Half*DelSep
!           New midpoint
            MidSep=BisSep+DelSep
!           Set the BBox over which to integrate the density,
            CubeRoot%Box=ExpandBox(CubeBox,MidSep)
#ifdef PERIODIC
            CALL MakeBoxPeriodic(CubeRoot%Box)
#endif
!           and the cooresponding penetration distance threshold
            PenetratDistanceThreshold=MidSep**2
!           Expand native Rho BBox by this amount
            CALL SetBBox(CubeRoot%Box,Box)
!           Find the exact electron density in the expanded BBox
            MidPop=PopInBox(RhoRoot)
!           Achieved-target errors         
            FMid=ABS((MidPop-REl)/REl)-TauRel
!           Convergence test
            IF(DelSep<TauRel*1.D1)EXIT
!           If still to the left, increment bisection point
            IF(FMid>Zero)BisSep=MidSep
WRITE(*,*)' MidPop = ',MidPop
         ENDDO
!        MidPop is now the most accurate evaluation of the integrated density
!        that HiCu can achieve with current thresholds.
         CubeRoot%IXact=MidPop
         CubeRoot%ISplit=1
         CubeRoot%Box%Tier=0
         CubeRoot%ECube=BIG_DBL
      END SUBROUTINE InitCubeRoot
!==========================================================================
!        
!==========================================================================
      SUBROUTINE NewCubeNode(Node,Level)
         TYPE(CubeNode), POINTER :: Node
         INTEGER                 :: Level
         INTEGER                 :: Status        
         ALLOCATE(Node,STAT=Status)
         IF(Status/=SUCCEED) &
            CALL Halt(' ALLOCATE 1 FAILED IN NewCubeNode')
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
         ALLOCATE(Node%Vals(NGrid,5),STAT=Status)
         CALL IncMem(Status,0,5*NGrid,'HiCu.CubeTree.Node%Vals')
      END SUBROUTINE NewCubeNode
!==========================================================================
!
!==========================================================================
      SUBROUTINE DeleteCubeGrid(Node)
         TYPE(CubeNode), POINTER :: Node
         INTEGER                 :: Status        
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
!========================================================================================
!     Compute the erf function
!========================================================================================
      FUNCTION ERF(W)
         REAL(DOUBLE),INTENT(IN) :: W
         REAL(DOUBLE)            :: ERF,X,Sgn
         INTEGER                 :: I,J
         Sgn=One; IF(W<0.0D0)Sgn=-One
         X=Sgn*W
         IF(X>Erf_Switch)THEN
            Erf=Sgn*1.0D0
         ELSE
            J=AINT(X*Erf_Grid)
            Erf=Sgn*(Erf_0(J)+X*(Erf_1(J)+X*(Erf_2(J)+X*Erf_3(J))))
         ENDIF
      END FUNCTION ERF
!========================================================================================
!     Compute the inverse exponential function, EXP(-X)
!========================================================================================
      FUNCTION ExpInv(X)
         REAL(DOUBLE), INTENT(IN) :: X
         REAL(DOUBLE)             :: EXPInv
         INTEGER                  :: J ,I
         IF(X.GE.Exp_Switch)THEN
            Expinv=0.0D0
         ELSE
            J=AINT(X*Exp_Grid)
            Expinv=Exp_0(J)+X*(Exp_1(J)+X*(Exp_2(J)+X*(Exp_3(J)+X*Exp_4(J))))
         ENDIF
      END FUNCTION ExpInv
!================================================================================
!     
!================================================================================
      RECURSIVE FUNCTION CubeWalk(Cube) RESULT(ICube)
         TYPE(CubeNode), POINTER    :: Cube
         REAL(DOUBLE), DIMENSION(3) :: ICube
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
      RECURSIVE FUNCTION LeafCount(Cube) 
         TYPE(CubeNode), POINTER    :: Cube
         INTEGER                    :: LeafCount
         INTEGER ::X
!--------------------------------------------------------------------------
         IF(Cube%Leaf)THEN
            LeafCount=1
         ELSE
            LeafCount=LeafCount(Cube%Descend)+LeafCount(Cube%Descend%Travrse)
         ENDIF
       END FUNCTION LeafCount
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

               55 FORMAT('CubeList=Append[CubeList,Point[{', &
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
      55 FORMAT('Rho[',D8.2,', ',D8.2,', ',D8.2,']= ',D12.6)
      END SUBROUTINE PrintCube
#ifdef PERIODIC
!==========================================================================
!
!==========================================================================
      SUBROUTINE MakeBoxPeriodic(Box)
        TYPE(BBox)       :: Box
        INTEGER          :: I
!
        DO I = 1,3
           IF(GM%AutoW(I)) THEN
              IF(Box%BndBox(I,1) < 1.D-9) THEN
                 Box%BndBox(I,1) = Zero
              ENDIF
              IF(Box%BndBox(I,2) > GM%BoxShape%D(I,I)-1.D-9) THEN
                 Box%BndBox(I,2) = GM%BoxShape%D(I,I)
              ENDIF
           ENDIF
        ENDDO
        DO I = 1,3
           Box%Half(I)   = Half*(Box%BndBox(I,2)-Box%BndBox(I,1))
           Box%Center(I) = Half*(Box%BndBox(I,2)+Box%BndBox(I,1))
        ENDDO
     END SUBROUTINE MakeBoxPeriodic
#endif
!
!
!
END MODULE

