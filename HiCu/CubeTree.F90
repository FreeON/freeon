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
      LOGICAL                                     :: Leaf
!     Intermediate values
      INTEGER                                     :: ISplit
      REAL(DOUBLE)                                :: IXact
      REAL(DOUBLE),DIMENSION(3)                   :: ICube    ! 
      REAL(DOUBLE),DIMENSION(3)                   :: ECube    ! 
!     Bounding box 
      TYPE(BBox)                                  :: Box          
!     Links
      TYPE(CubeNode), POINTER                     :: Travrse  ! Next node in tree traversal
      TYPE(CubeNode), POINTER                     :: Descend  ! Next node in tree descent
!     Cubature grid
      REAL(DOUBLE),   POINTER, DIMENSION(:,:) :: Grid     ! Transformed grid
      REAL(DOUBLE),   POINTER, DIMENSION(:)   :: Wght     ! Transformed Wght
      REAL(DOUBLE),   POINTER, DIMENSION(:,:) :: Vals     ! Values at each grid pt
   END TYPE CubeNode
!----------------------------------------------------------------------------------
   INTEGER, PARAMETER :: MaxTier=100
!  Global scalars 
   INTEGER            :: KOpt
   INTEGER            :: CubeNodes,CubeLevel
   REAL(DOUBLE)       :: LocalThresh
   REAL(DOUBLE)       :: Exc 
!  Global arrays
   INTEGER, DIMENSION(0:MaxTier) :: GlobalCubes 
!  Global Types 
   TYPE(BSET)                    :: BS 
   TYPE(CRDS)                    :: GM 
   TYPE(DBL_RNK4)                :: MD
!-----------!
   CONTAINS !
!================================================================================
!          
!================================================================================
      SUBROUTINE GridGen(CubeRoot)
         TYPE(CubeNode), POINTER          :: CubeRoot
         REAL(DOUBLE),PARAMETER           :: ErrorDeAmplification=2.D-1
         REAL(DOUBLE),   DIMENSION(3)     :: TotalError,LocalError,GlobalError, &
                                             RelativeError,NewCubes,OldCubes
         REAL(DOUBLE)                     :: MaxError,BoxSep,FullPop,OldPop, &
                                             Delta,TargtThresh
         INTEGER, PARAMETER               :: MaxIt=50
         REAL(DOUBLE), DIMENSION(0:MaxIt) :: MaxRelError
         INTEGER                          :: I,J,K,ErrCount,PtsPerAtom
         TYPE(BBox)                       :: CubeBox    
         CHARACTER(LEN=DEFAULT_CHR_LEN)   :: Mssg     
!---------------------------------------------------------------------------------
!        Initialize global variables
         GlobalCubes=0
         GlobalError=Zero       
!        Initialize root cube
         CALL NewCubeNode(CubeRoot,0) 
!        Find the CubeRoots bounding box
         FullPop=SetCubeRootBBox(CubeRoot)
         CubeRoot%ISplit=1
         CubeRoot%Box%Tier=0
         CALL ExpandBoxWalk(RhoRoot,-One)
         CubeRoot%ECube=BIG_DBL
!        Begin generation of the hierarchical grid
         MaxRelError=BIG_DBL
         LocalThresh=Thresholds%Cube*1.D2
         TargtThresh=Thresholds%Cube*ErrorDeAmplification
         Delta=(TargtThresh/LocalThresh)**(1.0D0/DBLE(MaxIt))
         DO J=1,MaxIt
            LocalThresh=LocalThresh*Delta
            CALL SetPenetrationThresh(LocalThresh*1.D-2,Exp_Switch)
            CALL ExpandBoxWalk(RhoRoot,One)
            CALL GridRefine(CubeRoot)
            CALL ExpandBoxWalk(RhoRoot,-One)
            NewCubes=CubeWalk(CubeRoot)
            PtsPerAtom=INT(DBLE(NGrid*LeafCount(CubeRoot))/DBLE(NAtoms))
            RelativeError(1)=ABS(FullPop-NewCubes(1))/DBLE(FullPop)
            IF(J==1)THEN
               OldCubes=NewCubes
            ELSE
               Exc=NewCubes(2) 
               RelativeError(2)=ABS((NewCubes(2)-OldCubes(2))/NewCubes(2))               
               RelativeError(3)=ABS((NewCubes(3)-OldCubes(3))/NewCubes(3))               
!               WRITE(*,9) LocalThresh,PtsPerAtom,RelativeError(1:3) 
             9 FORMAT(' LocThr = ',D8.2,' Pts/Atom= ',I5,' ERRORS: Rho=',D8.2, &
                      ' Exc=',D8.2,' Gxc=',D8.2)
            ENDIF
            OldCubes=NewCubes
         ENDDO
         IF(PrintFlags%Key>DEBUG_MEDIUM)THEN
            CALL OpenASCII(OutFile,Out)         
            WRITE(  *,10)PtsPerAtom,LocalThresh,RelativeError(1) 
!            WRITE(Out,10)PtsPerAtom,LocalThresh,RelativeError(1) 
            CLOSE(Out)
         10 FORMAT(' HiCu.GridGen     :: Pts/Atom= ',I6,' LocThr=',D8.2,', RhoErr=',D15.9)
         ENDIF
     END SUBROUTINE GridGen  
!================================================================================
!
!=================================================================================
     FUNCTION SetCubeRootBBox(CubeRoot) RESULT(FullPop)
         TYPE(CubeNode), POINTER          :: CubeRoot
         REAL(DOUBLE)                     :: BoxSep,Delta,NewPop,FullPop, &
                                             RelativeError
         REAL(DOUBLE), PARAMETER          :: PcntPop=1.D-2
         INTEGER                          :: I,J,K,KDig
         TYPE(BBox)                       :: CubeBox    
         CHARACTER(LEN=DEFAULT_CHR_LEN)   :: Mssg     
         
         CALL SetBBox(RhoRoot%Box,CubeBox)
         CALL SetPenetrationThresh(1.D-15)
         CALL ExpandBoxWalk(RhoRoot,One)
         CubeRoot%Box=ExpandBox(CubeBox,30.0D0)
#ifdef PERIODIC
         CALL MakeBoxPeriodic(CubeRoot%Box)
#endif
         CALL LayGrid(CubeRoot)
         CubeRoot%ISplit=1
         FullPop=CubeRoot%IXact
!
         BoxSep=1.0D0
         Delta=1.0D0
         KDig=0
         DO J=1,1000
          2 BoxSep=BoxSep+Delta
            CubeRoot%Box=ExpandBox(CubeBox,BoxSep)
#ifdef PERIODIC
            CALL MakeBoxPeriodic(CubeRoot%Box)
#endif
            CALL LayGrid(CubeRoot)
            NewPop=CubeRoot%IXact
            RelativeError=ABS((NewPop-FullPop)/FullPop)
            IF(RelativeError<Thresholds%Cube*PcntPop)THEN
               IF(KDig>8)GOTO 1
               KDig=KDig+1
               BoxSep=BoxSep-Delta
               Delta=Delta*1.D-1
            ENDIF
         ENDDO
         Mssg='Failed to find bounding box for CubeRoot! '//Rtrn &
            //'Integrated density in largest BBox = '            &
            //TRIM(DblToChar(Two*CubeRoot%IXact))//Rtrn          &
            //' Relative error = '//TRIM(DblToShrtChar(RelativeError))
         CALL Halt(Mssg)          
       1 CONTINUE
!
!         WRITE(*,*) ' Full Pop #1 = ',FullPop
!         WRITE(*,*) ' Box Size:'
!         WRITE(*,*) ' BndBox(xl) =',CubeRoot%Box%BndBox(1,1), &
!                    ' BndBox(xh) =',CubeRoot%Box%BndBox(1,2)
!         WRITE(*,*) ' BndBox(yl) =',CubeRoot%Box%BndBox(2,1), &
!                    ' BndBox(yh) =',CubeRoot%Box%BndBox(2,2)
!         WRITE(*,*) ' BndBox(zl) =',CubeRoot%Box%BndBox(3,1), &
!                    ' BndBox(zh) =',CubeRoot%Box%BndBox(3,2)
!
      END FUNCTION SetCubeRootBBox
!================================================================================
!          
!================================================================================
      RECURSIVE SUBROUTINE GridRefine(Cube)
         TYPE(CubeNode), POINTER :: Cube,Left,Right
         REAL(DOUBLE)            :: AbsErr,DHalf
         INTEGER                 :: ISplit
!------------------------------------------------------------------
         IF(Cube%Leaf)THEN
!            IF(Cube%IXact<LocalThresh)THEN
!               RETURN
            IF(Cube%ECube(1)>LocalThresh)THEN 
               CALL SplitCube(Cube)
               CALL GridRefine(Cube%Descend)
               CALL GridRefine(Cube%Descend%Travrse)
            ENDIF           
         ELSEIF(Branch(Cube))THEN
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
               Grid(:,I) = GridOld(:,I)+CS%CellCarts%D(:,NC)
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
            Rho(I)        =RhoV(1,I)
            AbsGradRho2(I)=RhoV(2,I)**2+RhoV(3,I)**2+RhoV(4,I)**2
            Cube%Vals(I,3)=RhoV(2,I)
            Cube%Vals(I,4)=RhoV(3,I)
            Cube%Vals(I,5)=RhoV(4,I)
            DO J=1,3; MaxGrad(J)=Max(MaxGrad(J),ABS(RhoV(J+1,I))); ENDDO
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
!            Cube%ICube(4)=Cube%ICube(4)+Cube%Wght(I)*dEdAbsGradRho2(I)*AbsGradRho2(I)
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
      FUNCTION Branch(Node)
         TYPE(CubeNode) :: Node
         LOGICAL        :: Branch
         Branch=.FALSE.
         IF(ASSOCIATED(Node%Descend)) Branch=.TRUE.
      END FUNCTION Branch
!================================================================================
!     
!================================================================================
      RECURSIVE FUNCTION LeafCount(Cube) 
         TYPE(CubeNode), POINTER    :: Cube
         INTEGER                    :: LeafCount
!--------------------------------------------------------------------------
         IF(Cube%Leaf)THEN
            LeafCount=1
         ELSE
            LeafCount=LeafCount(Cube%Descend)+LeafCount(Cube%Descend%Travrse)
         ENDIF
       END FUNCTION LeafCount
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
               Node%Grid(I,J)=Shift+CubeRuleGrid(I,J)*Slope
               Node%Wght(J)=Node%Wght(J)*Slope               
            ENDDO
         ENDDO
      END SUBROUTINE CubeRule
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
         Node%Leaf=.TRUE.
         Node%ECube=BIG_DBL
         NULLIFY(Node%Travrse)
         NULLIFY(Node%Descend)
         ALLOCATE(Node%Grid(3,NGrid),STAT=Status)
         CALL IncMem(Status,0,3*NGrid,'HiCu.CubeTree.Node%Grid')
         ALLOCATE(Node%Wght(NGrid),STAT=Status)
         CALL IncMem(Status,0,NGrid,'HiCu.CubeTree.Node%Wght')
         ALLOCATE(Node%Vals(NGrid,5),STAT=Status)
         CALL IncMem(Status,0,5*NGrid,'HiCu.CubeTree.Node%Vals')
      END SUBROUTINE NewCubeNode

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
!==========================================================================
!
!==========================================================================
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
               WRITE(Out,55)Node%Grid(1:3,I)

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
              IF(Box%BndBox(I,1) < 1.D-6) THEN
                 Box%BndBox(I,1) = Zero
              ENDIF
              IF(Box%BndBox(I,2) > GM%BoxShape%D(I,I)-1.D-6) THEN
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

