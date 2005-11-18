MODULE DIPMWTree
   USE DerivedTypes
   USE GlobalScalars   
   USE GlobalObjects
   USE ProcessControl
   USE Indexing
   USE InOut
   USE Macros
   USE BoundingBox
   USE RhoTree
   USE Functionals
   USE Thresholding
   USE AtomPairs
   USE SpecFun  
   USE DIPMWThresholds
   IMPLICIT NONE
!=================================================================
!  PIG Wavelet Node Types 
!=================================================================
   TYPE LinkArray
      TYPE(PIGWNode),POINTER             :: Link
   END TYPE LinkArray
   TYPE PIGWNode
      INTEGER                            :: Level
      INTEGER                            :: IType
      CHARACTER(LEN=8)                   :: LeafType
!     Links
      INTEGER                            :: NLinks
      TYPE(LinkArray),ALLOCATABLE        :: Links(:)
!     Bounding box 
      TYPE(BBox)                         :: Box       
!     The WCoefs    
      REAL(DOUBLE),ALLOCATABLE           :: WCoef(:)
!     For Integration
      REAL(DOUBLE),ALLOCATABLE           :: IntFactor(:)
   END TYPE PIGWNode
!==================================================================
!  Globals
!==================================================================
   TYPE(PIGWNode), POINTER         :: PIGWRoot 
   INTEGER                         :: CurrentLevel,NCallToRho
   INTEGER,PARAMETER               :: BigJ=20
   INTEGER,DIMENSION(0:BigJ)       :: PIGWNodes
   REAL(DOUBLE),DIMENSION(1:3)     :: X0,Length
   REAL(DOUBLE),DIMENSION(0:BigJ)  :: JFactor,InvJFactor,ScaleX,ScaleY,ScaleZ
!------------------------------------------------------------------
! Points to check around Wavelet Node
   INTEGER,PARAMETER               :: NCheck=26
   REAL(DOUBLE),DIMENSION(NCheck)  :: XCheck=(/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, &
                                                0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1/)
   REAL(DOUBLE),DIMENSION(NCheck)  :: YCheck=(/ 1, 1, 1, 0, 0, 0,-1,-1,-1, 1, 1, 1, 0, &
                                                0,-1,-1,-1, 1, 1, 1, 0, 0, 0,-1,-1,-1/)
   REAL(DOUBLE),DIMENSION(NCheck)  :: ZCheck=(/-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, &
                                                1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1/)
!------------------------------------------------------------------
!  Wavelet Order Stuff
!
!
   INTEGER,PARAMETER                 :: DOrder=3
   INTEGER,PARAMETER                 :: NOrder=0
   INTEGER,PARAMETER                 :: DMax=DOrder+1
   INTEGER,PARAMETER                 :: LenDOrder=(DOrder+1)**3
   INTEGER,PARAMETER                 :: LenNOrder=2*NOrder+2
   REAL(DOUBLE),DIMENSION(LenDOrder) :: RhoTmp
   REAL(DOUBLE),DIMENSION(LenNOrder,DMax)      :: WvltInt
   REAL(DOUBLE),DIMENSION(LenNOrder,DMax,DMax) :: WFilter
!------------------------------------------------------------------
!  RhoAtPoint and HasNode Stuff
!
   LOGICAL                           :: IsTrue
   REAL(DOUBLE),DIMENSION(LenDOrder) :: RhoSum
   REAL(DOUBLE),DIMENSION(1:3)       :: Xpos
   TYPE(BBox)                        :: GBox       
!------------------------------------------------------------------
!  Exc Energy
!
   REAL(DOUBLE)                      :: XCEnergy
   CONTAINS 
!==================================================================
!  Generate the Wavelet Rep of the XC potential
!==================================================================
   SUBROUTINE  DIPMWTree(MaxJ)
     INTEGER                                :: MaxJ,J,I,K,TotWav,NN
     REAL(DOUBLE),DIMENSION(3)              :: X
     REAL(DOUBLE)                           :: TauDIPMW_old,DX,DY,DZ,FFac,SUM,RhoInt
     TYPE(TIME)                             :: TimeForLevel
!------------------------------------------------------------------
!    Error Check
     IF(MaxJ > BigJ) THEN
        CALL Halt("MaxJ > BigJ")
     ENDIF
!    Initialize counters
     TauDIPMW_old = TauDIPMW
     PIGWNodes(:) = 0
!    Intialize Filters
     CALL InitializeWFilters()
!    Initialize the lengths  
     Length(:)        = ABS(RhoRoot%Box%BndBox(:,2)-RhoRoot%Box%BndBox(:,1))
     X0(:)            = RhoRoot%Box%BndBox(:,1)+Half*Length(:)
!    Global Bounding Box
     GBox%Center(:)   = RhoRoot%Box%Center(:)
     GBox%Half(:)     = RhoRoot%Box%Half(:)
     GBox%BndBox(:,1) = RhoRoot%Box%BndBox(:,1)-1.D-10
     GBox%BndBox(:,2) = RhoRoot%Box%BndBox(:,2)+1.D-10    
!!$!
!!$     Length(1) = 2.D0
!!$     Length(2) = 2.D0
!!$     Length(3) = 2.D0
!!$     X0(1)     = 0.D0
!!$     X0(2)     = 0.D0
!!$     X0(3)     = 0.D0
!!$!
!!$     GBox%Center(:)   = X0(:)
!!$     GBox%Half(:)     = Half*Length(:)
!!$     GBox%BndBox(:,1) = X0(:)-GBox%Half(:)-1.D-10
!!$     GBox%BndBox(:,2) = X0(:)+GBox%Half(:)+1.D-10
!
!    Initialize JFactor
     DO J=0,BigJ
        JFactor(J)    = Two**(J+1)
        InvJFactor(J) = Half**(J+1)  
        ScaleX(J)     = InvJFactor(J)*Length(1)
        ScaleY(J)     = InvJFactor(J)*Length(2)
        ScaleZ(J)     = InvJFactor(J)*Length(3)  
     ENDDO
!
     NCallToRho=0
!    Initialize the Root Node 
     CALL InitializePIGWRoot(X0,0)
!
!!$     WRITE(*,*) "Length"
!!$     WRITE(*,*)  Length(1:3)
!!$     WRITE(*,*) "Box"
!!$     WRITE(*,*)  RhoRoot%Box%BndBox(:,1)
!!$     WRITE(*,*)  RhoRoot%Box%BndBox(:,2) 
!!$     WRITE(*,*) "X0"
!!$     WRITE(*,*) X0(1:3)
!!$     WRITE(*,*)
!!$     CALL PrintNode(PIGWRoot)
!
!    Recursively generate the Wavelet Rep of the Density, XC Energy and XC potential
     TotWav  = 0
     DO CurrentLevel=1,MaxJ
        WRITE(*,*) 'CurrentLevel = ',CurrentLevel,TauDIPMW
!
        CALL Elapsed_Time(TimeForLevel,'Init')
        CALL MakePIGWTree(PIGWRoot)
        CALL Elapsed_Time(TimeForLevel,'Accum')
!  
        TotWav = TotWav+PIGWNodes(CurrentLevel)
        WRITE(*,*)  "   ",PIGWNodes(CurrentLevel)
        WRITE(*,*)  "   ",TimeForLevel%CPUS,TimeForLevel%CPUS/PIGWNodes(CurrentLevel-1)
        IF(PIGWNodes(CurrentLevel) == 0) EXIT
     ENDDO
!
     WRITE(*,*) "Total Wavelets        = ",TotWav
     WRITE(*,*) "Total Call to RhoTree = ",NCallToRho
     WRITE(*,*) "Wavelets per Atom     = ",TotWav/GM%Natms
     XCEnergy = Zero
     CALL ComputeExc(PIGWRoot)
     WRITE(*,*) "XCEnergy=",XCEnergy 
!
   END SUBROUTINE DIPMWTree
!==================================================================
!  Initialize the Wavelet Filters
!==================================================================
   SUBROUTINE InitializeWFilters()
!
     WvltInt = Zero
     SELECT CASE(DOrder)
     CASE(0)
        SELECT CASE(NOrder)
        CASE(0)
           WvltInt(1,1) = 0.0D0
           WvltInt(2,1) = 5.0D-1
!
           WFilter(1,1,1)=-5.0D-1
           WFilter(2,1,1)=-5.0D-1
        CASE(1)
           WvltInt(1,1) = 0.0D0
           WvltInt(2,1) = 13.0D0/24.0D0
           WvltInt(3,1) = 5.0D-1
           WvltInt(4,1) = 5.0D-1
!
           WFilter(1,1,1)= 6.25D-2
           WFilter(2,1,1)=-5.625D-1
           WFilter(3,1,1)=-5.625D-1
           WFilter(4,1,1)= 6.25D-2
        CASE(2)
           WvltInt(1,1) = 0.0D0
           WvltInt(2,1) = 0.0D0
           WvltInt(3,1) = 0.0D0
           WvltInt(4,1) = 0.0D0
           WvltInt(5,1) = 0.0D0
           WvltInt(6,1) = 5.0D-1
!
           WFilter(1,1,1)=-1.171875D-2
           WFilter(2,1,1)= 9.765625D-2
           WFilter(3,1,1)=-5.859375D-1
           WFilter(4,1,1)=-5.859375D-1
           WFilter(5,1,1)= 9.765625D-2
           WFilter(6,1,1)=-1.171875D-2
        CASE(3)
           WvltInt(1,1) = 0.0D0
           WvltInt(2,1) = 0.0D0
           WvltInt(3,1) = 0.0D0
           WvltInt(4,1) = 0.0D0
           WvltInt(5,1) = 0.0D0
           WvltInt(6,1) = 0.0D0
           WvltInt(7,1) = 0.0D0
           WvltInt(8,1) = 5.0D-1
!
           WFilter(1,1,1)= 2.4414062500D-3
           WFilter(2,1,1)=-2.3925781250D-2
           WFilter(3,1,1)= 1.1962890625D-1
           WFilter(4,1,1)=-5.9814453125D-1
           WFilter(5,1,1)=-5.9814453125D-1
           WFilter(6,1,1)= 1.1962890625D-1
           WFilter(7,1,1)=-2.3925781250D-2
           WFilter(8,1,1)= 2.4414062500D-3
        CASE DEFAULT
           CALL Halt('WFilter Size for DOrder==0 not defined')
        END SELECT
     CASE(1)
        SELECT CASE(NOrder)
        CASE(0)
           WvltInt(1,1) = 0.0000000000000000000D0
           WvltInt(2,1) = 5.0000000000000000000D-1
           WvltInt(1,2) = 0.0000000000000000000D0
           WvltInt(2,2) = 8.3333333333333333333D-2
!
           WFilter(1,1,1) = -5.0D-1 
           WFilter(2,1,1) = -5.0D-1
           WFilter(1,2,1) = -2.5D-1
           WFilter(2,2,1) =  2.5D-1
!
           WFilter(1,1,2) =  7.5D-1
           WFilter(2,1,2) = -7.5D-1
           WFilter(1,2,2) =  2.5D-1
           WFilter(2,2,2) =  2.5D-1
        CASE DEFAULT
           CALL Halt('WFilter Size for DOrder==1 not defined')
        END SELECT
     CASE(2)
        SELECT CASE(NOrder)
        CASE(0)
           WvltInt(1,1) = 0.0000000000000000000D0
           WvltInt(2,1) = 5.0000000000000000000D-1
           WvltInt(1,2) = 0.0000000000000000000D0
           WvltInt(2,2) = 1.0000000000000000000D-1
           WvltInt(1,3) = 0.0000000000000000000D0
           WvltInt(2,3) = 8.3333333333333333333D-3
!
           WFilter(1,1,1) = -5.000D-1
           WFilter(2,1,1) = -5.000D-1
           WFilter(1,2,1) = -3.125D-1
           WFilter(2,2,1) =  3.125D-1
           WFilter(1,3,1) = -6.250D-2
           WFilter(2,3,1) = -6.250D-2
!
           WFilter(1,1,2) =  9.375D-1
           WFilter(2,1,2) = -9.375D-1
           WFilter(1,2,2) =  4.375D-1
           WFilter(2,2,2) =  4.375D-1
           WFilter(1,3,2) =  6.250D-2
           WFilter(2,3,2) = -6.250D-2
!
           WFilter(1,1,3) =  0.000D0
           WFilter(2,1,3) =  0.000D0
           WFilter(1,2,3) =  7.500D-1
           WFilter(2,2,3) = -7.500D-1
           WFilter(1,3,3) =  2.500D-1
           WFilter(2,3,3) =  2.500D-1
        CASE DEFAULT
           CALL Halt('WFilter Size for DOrder==2 not defined')
        END SELECT
     CASE(3)
        SELECT CASE(NOrder)
        CASE(0)
           WvltInt(1,1) = 0.0000000000000000000D0
           WvltInt(2,1) = 5.0000000000000000000D-1
           WvltInt(1,2) = 0.0000000000000000000D0
           WvltInt(2,2) = 1.0714285714285714286D-1
           WvltInt(1,3) = 0.0000000000000000000D0
           WvltInt(2,3) = 1.1904761904761904762D-2
           WvltInt(1,4) = 0.0000000000000000000D0
           WvltInt(2,4) = 5.9523809523809523809D-4
!
           WFilter(1,1,1) = -5.0000D-1
           WFilter(2,1,1) = -5.0000D-1
           WFilter(1,2,1) = -3.4375D-1
           WFilter(2,2,1) =  3.4375D-1
           WFilter(1,3,1) = -9.3750D-2
           WFilter(2,3,1) = -9.3750D-2
           WFilter(1,4,1) = -1.041666666666666666D-2
           WFilter(2,4,1) =  1.041666666666666666D-2
!
           WFilter(1,1,2) =  1.09375D0
           WFilter(2,1,2) = -1.09375D0
           WFilter(1,2,2) =  5.93750D-1
           WFilter(2,2,2) =  5.93750D-1
           WFilter(1,3,2) =  1.25000D-1
           WFilter(2,3,2) = -1.25000D-1
           WFilter(1,4,2) =  1.041666666666666666D-2
           WFilter(2,4,2) =  1.041666666666666666D-2
!
           WFilter(1,1,3) =  0.000D0
           WFilter(2,1,3) =  0.000D0
           WFilter(1,2,3) =  9.375D-1
           WFilter(2,2,3) = -9.375D-1
           WFilter(1,3,3) =  4.375D-1
           WFilter(2,3,3) =  4.375D-1
           WFilter(1,4,3) =  6.250D-2
           WFilter(2,4,3) = -6.250D-2
!
           WFilter(1,1,4) = -6.5625D0
           WFilter(2,1,4) =  6.5625D0
           WFilter(1,2,4) = -6.5625D0
           WFilter(2,2,4) = -6.5625D0
           WFilter(1,3,4) = -1.8750D0
           WFilter(2,3,4) =  1.8750D0
           WFilter(1,4,4) = -1.8750D-1
           WFilter(2,4,4) = -1.8750D-1
        CASE DEFAULT
           CALL Halt('WFilter Size for DOrder==3 not defined')
        END SELECT
     CASE DEFAULT
        CALL Halt('WFilter DOrder>3 not defined')
     END SELECT
!
   END SUBROUTINE InitializeWFilters
!==================================================================
!  Initialize PIGWTree
!==================================================================
   SUBROUTINE InitializePIGWRoot(X,ILevel)
     INTEGER                            :: ILevel,Status,ILink,IX,IY,IZ,JLink
     REAL(DOUBLE),DIMENSION(LenDOrder)  :: RhoX,ExcX,VxcX
     REAL(DOUBLE),DIMENSION(3)          :: X,XX,DX
!        
!     PIGWRoot
!
     CALL NewPIGWNode(PIGWRoot,ILevel,1000,X)
!    Compute the coefs of 000
     RhoX = RhoAtX(X)
     CALL ExcVxc(RhoX,ExcX,VxcX)
     CALL NormalizeWCoefs(ExcX,ILevel,1000)
!    Store
     PIGWRoot%WCoef   = ExcX
     PIGWRoot%NLinks  = 26
     PIGWRoot%LeafType='NodeNode'
!    Allocate the Links
     ALLOCATE(PIGWRoot%Links(PIGWRoot%NLinks),STAT=Status)
     IF(Status/=SUCCEED)CALL Halt('Link ALLOCATE failed in InitializePIGWRoot ')
!
!    Surface Scaling Functions 
!
     DX(1)    = ScaleX(0)
     DX(2)    = ScaleY(0)
     DX(3)    = ScaleZ(0)
     DO ILink = 1, PIGWRoot%NLinks
        IX = XCheck(ILink)
        IY = YCheck(ILink) 
        IZ = ZCheck(ILink)
        XX(1) = PIGWRoot%Box%Center(1)+IX*DX(1)
        XX(2) = PIGWRoot%Box%Center(2)+IY*DX(2)
        XX(3) = PIGWRoot%Box%Center(3)+IZ*DX(3)
        CALL NewPIGWNode(PIGWRoot%Links(ILink)%Link,ILevel,1000,XX)
        RhoX = RhoAtX(XX)
        CALL ExcVxc(RhoX,ExcX,VxcX)
        CALL NormalizeWCoefs(ExcX,ILevel,1000)
        PIGWRoot%Links(ILink)%Link%WCoef = ExcX
     ENDDO
!
   END SUBROUTINE InitializePIGWRoot
!==================================================================
!  Print a Node
!==================================================================
   SUBROUTINE PrintNode(Node)
     TYPE(PIGWNode), POINTER         :: Node
     INTEGER                         :: I,J,K,IJK
!
     WRITE(*,*) "----------------------------"
     WRITE(*,*) 'Level    = ',Node%Level
     WRITE(*,*) 'IType    = ',Node%IType
     WRITE(*,*) 'LeafType = ',Node%LeafType
     WRITE(*,*) 'NLinks   = ',Node%NLinks
     WRITE(*,*) 'Box'
     WRITE(*,*) Node%Box%Center(1:3)
     WRITE(*,*) Node%Box%Half(1:3)
     WRITE(*,*) Node%Box%BndBox(1:3,1)
     WRITE(*,*) Node%Box%BndBox(1:3,2)
     WRITE(*,*) 'WCoef'
     DO K=0,DOrder
        DO J=0,DOrder
           IJK = 1+DMax*J+DMax*DMax*K
           WRITE(*,'3(5X,I1,I1,I1,1X,D16.8)') (I,J,K,Node%WCoef(IJK+I),I=0,DOrder) 
        ENDDO
     ENDDO
     WRITE(*,*) 'Integration  Factor'
     DO K=0,DOrder
        DO J=0,DOrder
           IJK = 1+DMax*J+DMax*DMax*K
           WRITE(*,'3(5X,I1,I1,I1,1X,D16.8)') (I,J,K,Node%IntFactor(IJK+I),I=0,DOrder) 
        ENDDO
     ENDDO
     WRITE(*,*) "----------------------------"
!
   END SUBROUTINE PrintNode
!==================================================================
!  Allocate Node
!==================================================================
   SUBROUTINE NewPIGWNode(Node,ILevel,IType,X)
     TYPE(PIGWNode), POINTER         :: Node
     INTEGER                         :: ILevel,Status
     INTEGER                         :: IType
     REAL(DOUBLE),DIMENSION(3)       :: X
!
     ALLOCATE(Node,STAT=Status)
     IF(Status/=SUCCEED)CALL Halt(' Node ALLOCATE failed in NewPIGWNode ')
!
     Node%LeafType = 'LeafLeaf'
     Node%Level    = ILevel
     Node%IType    = IType
     Node%NLinks   = 0
!    Wavelet Coefs
     ALLOCATE(Node%WCoef(1:LenDOrder),STAT=Status)
     IF(Status/=SUCCEED)CALL Halt(' Node%WCoef ALLOCATE failed in NewPIGWNode ')
!    Integration Factor
     ALLOCATE(Node%IntFactor(1:LenDOrder),STAT=Status)
     IF(Status/=SUCCEED)CALL Halt(' Node%IntFactor ALLOCATE failed in NewPIGWNode ')
!    Intitialize WCoefs
     Node%WCoef(:) = Zero
!    Initialize the Bounding Box
     CALL NodeBoundingBox(Node%Level,Node%IType,X,Node%Box)
!    Integration Factor
     CALL WWeight(Node%Level,Node%IType,Node%Box,Node%IntFactor)
!    Globals
     PIGWNodes(ILevel) = PIGWNodes(ILevel)+1
!
   END SUBROUTINE NewPIGWNode
!======================================================================================
!
!======================================================================================
   SUBROUTINE NodeBoundingBox(Level,IType,XX,Box)
     TYPE(BBox)                      :: Box     
     INTEGER                         :: Level,IType
     REAL(DOUBLE),DIMENSION(3)       :: XX
     REAL(DOUBLE)                    :: BoxFactor

!    Initialize Center
     Box%Center(:) = XX(:)
!    Make Bounding Box
     BoxFactor=DBLE(2*NOrder+1)
     IF(IType==1000) THEN
        Box%Half(1)     = 1.0D0*BoxFactor*ScaleX(Level)
        Box%Half(2)     = 1.0D0*BoxFactor*ScaleY(Level)
        Box%Half(3)     = 1.0D0*BoxFactor*ScaleZ(Level)
     ELSEIF(IType==1100) THEN
        Box%Half(1)     = 1.0D0*BoxFactor*ScaleX(Level)
        Box%Half(2)     = 2.0D0*BoxFactor*ScaleY(Level)
        Box%Half(3)     = 2.0D0*BoxFactor*ScaleZ(Level)
     ELSEIF(IType==1010) THEN
        Box%Half(1)     = 1.0D0*BoxFactor*ScaleX(Level)
        Box%Half(2)     = 1.0D0*BoxFactor*ScaleY(Level)
        Box%Half(3)     = 2.0D0*BoxFactor*ScaleZ(Level)
     ELSEIF(IType==1001) THEN
        Box%Half(1)     = 1.0D0*BoxFactor*ScaleX(Level)
        Box%Half(2)     = 1.0D0*BoxFactor*ScaleY(Level)
        Box%Half(3)     = 1.0D0*BoxFactor*ScaleZ(Level)
     ENDIF
!    Boundding Box
     Box%BndBox(:,1) = Box%Center(:)-Box%Half(:)
     Box%BndBox(:,2) = Box%Center(:)+Box%Half(:)
     Box%BndBox(:,1) = MAX(Box%BndBox(:,1),GBox%BndBox(:,1))
     Box%BndBox(:,2) = MIN(Box%BndBox(:,2),GBox%BndBox(:,2))
!
   END SUBROUTINE NodeBoundingBox
!======================================================================================
!
!======================================================================================
   SUBROUTINE WWeight(Level,IType,Box,IWeight)
     TYPE(BBox)                        :: Box     
     INTEGER                           :: Level,IType,I,J,K,IJK,NLx,NLy,NLz,NRx,NRy,NRz
     REAL(DOUBLE)                      :: SX,SY,SZ,SXYZ,Lx,Ly,Lz,Rx,Ry,Rz,BoxFactor
     REAL(DOUBLE)                      :: WX,WY,WZ
     REAL(DOUBLE),DIMENSION(LenDOrder) :: IWeight
!
     BoxFactor = DBLE(2*NOrder+1)
     Lx = (Box%Center(1)-Box%BndBox(1,1))/Box%Half(1)
     Ly = (Box%Center(2)-Box%BndBox(2,1))/Box%Half(2)
     Lz = (Box%Center(3)-Box%BndBox(3,1))/Box%Half(3)
     Rx = (Box%BndBox(1,2)-Box%Center(1))/Box%Half(1)
     Ry = (Box%BndBox(2,2)-Box%Center(2))/Box%Half(2)
     Rz = (Box%BndBox(3,2)-Box%Center(3))/Box%Half(3)
     NLx = INT(BoxFactor*(Lx+1.D-10))+1
     NLy = INT(BoxFactor*(Ly+1.D-10))+1
     NLz = INT(BoxFactor*(Lz+1.D-10))+1
     NRx = INT(BoxFactor*(Rx+1.D-10))+1
     NRy = INT(BoxFactor*(Ry+1.D-10))+1
     NRz = INT(BoxFactor*(Rz+1.D-10))+1
!
     SXYZ     = (Box%Half(1)/BoxFactor)*(Box%Half(2)/BoxFactor)*(Box%Half(3)/BoxFactor)
     IWeight  = Zero
     DO K=0,DOrder
        DO J=0,DOrder
           DO I=0,DOrder
              WX  = ((-One)**I)*WvltInt(NLx,I+1)+WvltInt(NRx,I+1)
              WY  = ((-One)**J)*WvltInt(NLy,J+1)+WvltInt(NRy,J+1)
              WZ  = ((-One)**K)*WvltInt(NLz,K+1)+WvltInt(NRz,K+1)
              IJK = 1+I+DMax*J+DMax*DMax*K
              IWeight(IJK) = WX*WY*WZ*SXYZ
           ENDDO 
        ENDDO
     ENDDO
!
   END SUBROUTINE WWeight
!==================================================================
!  Recursively Genergate the Wavelet Rep of the XC potential
!==================================================================
   RECURSIVE SUBROUTINE MakePIGWTree(Node)
     TYPE(PIGWNode), POINTER               :: Node
     TYPE(BBox)                            :: NewBox     
     REAL(DOUBLE)                          :: TotalWeight
     REAL(DOUBLE),DIMENSION(3)             :: DX,XX
     INTEGER                               :: NewLevel,ICk,IX,IY,IZ,NewIType,ILink,Status,NLinks,I,J,K
     REAL(DOUBLE),DIMENSION(LenDOrder)     :: WCoef,IntValue
!
     INTEGER,DIMENSION(26)                 :: ITypeArray
     INTEGER,DIMENSION(26)                 :: LeafTypeArray
     REAL(DOUBLE),DIMENSION(26,3)          :: XposArray
     REAL(DOUBLE), DIMENSION(26,LenDOrder) :: WCoefArray,IntVArray
!
     IF(Node%LeafType=='LeafLeaf') THEN 
        NewLevel = Node%Level+1
        DX(1)    = ScaleX(NewLevel)
        DX(2)    = ScaleY(NewLevel)
        DX(3)    = ScaleZ(NewLevel)
!       Determine the Sub-Wavelets that will exist
        NLinks = 0
        DO ICk=1,NCheck
           IX = XCheck(ICk)
           IY = YCheck(ICk) 
           IZ = ZCheck(ICk)
           IF(IX .NE. 0) NewIType = 1100  
           IF(IY .NE. 0) NewIType = 1010  
           IF(IZ .NE. 0) NewIType = 1001  
           XX(1) = Node%Box%Center(1)+IX*DX(1)
           XX(2) = Node%Box%Center(2)+IY*DX(2)
           XX(3) = Node%Box%Center(3)+IZ*DX(3)
           IF(.NOT. HasNode(XX)) THEN
              CALL ComputeWCoef(XX,NewLevel,NewIType,WCoef)
              CALL NodeBoundingBox(NewLevel,NewIType,XX,NewBox)
              CALL WWeight(NewLevel,NewIType,NewBox,IntValue)
              TotalWeight = Zero
              DO I=1,LenDOrder
                 TotalWeight = TotalWeight+IntValue(I)*ABS(WCoef(I))
              ENDDO
              NLinks = NLinks+1
              ITypeArray(NLinks)              = NewIType
              IF(ABS(TotalWeight) > TauDIPMW) THEN
                 LeafTypeArray(NLinks)        = 0
              ELSE
                 LeafTypeArray(NLinks)        = 1
              ENDIF
              XposArray(NLinks,1:3)          = XX(1:3)
              WCoefArray(NLinks,1:LenDOrder) = WCoef(1:LenDOrder)
              IntVArray(NLinks,1:LenDOrder)  = IntValue(1:LenDOrder)
           ENDIF
        ENDDO 
!       Allocate the Link Node
        Node%NLinks=NLinks
        IF(Node%NLinks == 0) THEN
           Node%LeafType='EnddLeaf'
        ELSE
           Node%LeafType='NodeNode'
        ENDIF
!       Allocate the Links
        ALLOCATE(Node%Links(Node%NLinks),STAT=Status)
        IF(Status/=SUCCEED)CALL Halt('Link ALLOCATE failed in MakePIGWTree ')
!       Create the New SubWavelets Nodes and Store the Relavent Info
        DO ILink = 1, Node%NLinks
           CALL  NewPIGWNode(Node%Links(ILink)%Link,NewLevel,ITypeArray(ILink),XposArray(ILink,1:3))
           Node%Links(ILink)%Link%WCoef     = WCoefArray(ILink,1:LenDOrder)
           Node%Links(ILink)%Link%IntFactor = IntVArray(ILink,1:LenDOrder)
           IF(LeafTypeArray(ILink)==0) THEN   
              Node%Links(ILink)%Link%LeafType="LeafLeaf"
           ELSE
!              Node%Links(ILink)%Link%LeafType="ContLeaf" 
              Node%Links(ILink)%Link%LeafType="EnddLeaf" 
           ENDIF
        ENDDO
!       We are Done, Return
        RETURN
     ELSEIF(Node%LeafType=='NodeNode') THEN
        DO ILink=1,Node%NLinks
           CALL MakePIGWTree(Node%Links(ILink)%Link)
        ENDDO
     ENDIF
     RETURN
!
   END SUBROUTINE MakePIGWTree
!==================================================================
! Normailze the WCoefs
!==================================================================
   SUBROUTINE NormalizeWCoefs(WCoef,ILevel,IType)
     INTEGER                           :: ILevel,IType,I,J,K,IJK
     REAL(DOUBLE),DIMENSION(LenDOrder) :: WCoef
     REAL(DOUBLE)                      :: DX,DY,DZ       
!
     SELECT CASE(IType)
     CASE(1000)
        DX = One*ScaleX(ILevel)
        DY = One*ScaleY(ILevel)
        DZ = One*ScaleZ(ILevel)
     CASE(1100)
        DX = One*ScaleX(ILevel)
        DY = Two*ScaleY(ILevel)
        DZ = Two*ScaleZ(ILevel)
     CASE(1010)
        DX = One*ScaleX(ILevel)
        DY = One*ScaleY(ILevel)
        DZ = Two*ScaleZ(ILevel)
     CASE(1001)
        DX = One*ScaleX(ILevel)
        DY = One*ScaleY(ILevel)
        DZ = One*ScaleZ(ILevel)
     END SELECT
     DO K=0,DOrder
        DO J=0,DOrder
           DO I=0,DOrder
              IJK = 1+I+DMax*J+DMax*DMax*K
              WCoef(IJK) = WCoef(IJK)*(DX**I)*(DY**J)*(DZ**K)
           ENDDO
        ENDDO
     ENDDO
!
   END SUBROUTINE NormalizeWCoefs
!==================================================================
!  Calculate Coefs from the Density Rho
!==================================================================
   SUBROUTINE ComputeWCoef(XX,NewLevel,NewIType,WCoef)
       REAL(DOUBLE),DIMENSION(3)         :: XX,XF
       REAL(DOUBLE),DIMENSION(3)         :: DX,DX2,DX3,XD,XD2,XD3
       INTEGER                           :: NewLevel,NewIType,NN,I,J,K
       INTEGER                           :: IJK0,IJK1,IJK2,IJK3
       REAL(DOUBLE),DIMENSION(LenDOrder) :: RhoX,WCoef,Tmp1,Tmp2
       REAL(DOUBLE)                      :: SX,SY,SZ
!
       DX(1) = ScaleX(NewLevel)
       DX(2) = ScaleY(NewLevel)
       DX(3) = ScaleZ(NewLevel)
       XD(1) = One/DX(1)
       XD(2) = One/DX(2)
       XD(3) = One/DX(3)
!
       DX2(1) = DX(1)*DX(1) 
       DX2(2) = DX(2)*DX(2) 
       DX2(3) = DX(3)*DX(3)
       XD2(1) = XD(1)*XD(1)
       XD2(2) = XD(2)*XD(2)
       XD2(3) = XD(3)*XD(3)
!
       DX3(1) = DX2(1)*DX(1) 
       DX3(2) = DX2(2)*DX(2) 
       DX3(3) = DX2(3)*DX(3)
       XD3(1) = XD2(1)*XD(1)
       XD3(2) = XD2(2)*XD(2)
       XD3(3) = XD2(3)*XD(3)
!
       XF    = XX
       RhoX  = RhoAtX(XX)
       CALL ExcVxc(RhoX,WCoef,Tmp1)
       SELECT CASE(DOrder)
       CASE(0)
          SELECT CASE(NewIType)
          CASE(1100)
             DO NN=1,LenNOrder
                XF(1) = XX(1)+DX(1)*(2*NN-1-LenNOrder)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp2,Tmp1)
                WCoef(1) = WCoef(1) + WFilter(NN,1,1)*Tmp2(1)
             ENDDO
          CASE(1010)
             DO NN=1,LenNOrder
                XF(2) = XX(2)+DX(2)*(2*NN-1-LenNOrder)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp2,Tmp1)
                WCoef(1) = WCoef(1) + WFilter(NN,1,1)*Tmp2(1)
             ENDDO
          CASE(1001)
             DO NN=1,LenNOrder
                XF(3) = XX(3)+DX(3)*(2*NN-1-LenNOrder)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp2,Tmp1)
                WCoef(1) = WCoef(1) + WFilter(NN,1,1)*Tmp2(1)
             ENDDO
          END SELECT
       CASE(1)
          SELECT CASE(NewIType)
          CASE(1100)
             DO NN=1,LenNOrder
                XF(1) = XX(1)+DX(1)*(2*NN-1-LenNOrder)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp2,Tmp1)
                DO K=0,DOrder
                   DO J=0,DOrder
                      IJK0 = 1+DMax*J+DMax*DMax*K
                      IJK1 = IJK0+1
                      WCoef(IJK0) = WCoef(IJK0) +       WFilter(NN,1,1)*Tmp2(IJK0)   &
                                                + DX(1)*WFilter(NN,2,1)*Tmp2(IJK1)
                      WCoef(IJK1) = WCoef(IJK1) + XD(1)*WFilter(NN,1,2)*Tmp2(IJK0)   &
                                                +       WFilter(NN,2,2)*Tmp2(IJK1)
                   ENDDO
                ENDDO
             ENDDO
          CASE(1010)
             DO NN=1,LenNOrder
                XF(2) = XX(2)+DX(2)*(2*NN-1-LenNOrder)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp2,Tmp1)
                DO K=0,DOrder
                   DO I=0,DOrder
                      IJK0 = 1+I+DMax*DMax*K
                      IJK1 = IJK0+DMax
                      WCoef(IJK0) = WCoef(IJK0) +       WFilter(NN,1,1)*Tmp2(IJK0)    &
                                                + DX(2)*WFilter(NN,2,1)*Tmp2(IJK1)
                      WCoef(IJK1) = WCoef(IJK1) + XD(2)*WFilter(NN,1,2)*Tmp2(IJK0)    &
                                                +       WFilter(NN,2,2)*Tmp2(IJK1)
                   ENDDO
                ENDDO
             ENDDO
          CASE(1001)
            DO NN=1,LenNOrder
                XF(3) = XX(3)+DX(3)*(2*NN-1-LenNOrder)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp2,Tmp1)
                DO J=0,DOrder
                   DO I=0,DOrder
                      IJK0 = 1+I+DMax*J
                      IJK1 = IJK0+DMax*DMax
                      WCoef(IJK0) = WCoef(IJK0) +       WFilter(NN,1,1)*Tmp2(IJK0)    &
                                                + DX(3)*WFilter(NN,2,1)*Tmp2(IJK1)
                      WCoef(IJK1) = WCoef(IJK1) + XD(3)*WFilter(NN,1,2)*Tmp2(IJK0)    &
                                                +       WFilter(NN,2,2)*Tmp2(IJK1)
                   ENDDO
                ENDDO
             ENDDO
          END SELECT
       CASE(2)
          SELECT CASE(NewIType)
          CASE(1100)
             DO NN=1,LenNOrder
                XF(1) = XX(1)+DX(1)*(2*NN-1-LenNOrder)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp2,Tmp1)
                DO K=0,DOrder
                   DO J=0,DOrder
                      IJK0 = 1+DMax*J+DMax*DMax*K
                      IJK1 = IJK0+1
                      IJK2 = IJK0+2
                      WCoef(IJK0) = WCoef(IJK0)  +        WFilter(NN,1,1)*Tmp2(IJK0)   &
                                                 +  DX(1)*WFilter(NN,2,1)*Tmp2(IJK1)   &
                                                 + DX2(1)*WFilter(NN,3,1)*Tmp2(IJK2)  
!
                      WCoef(IJK1) = WCoef(IJK1)  +  XD(1)*WFilter(NN,1,2)*Tmp2(IJK0)   &
                                                 +        WFilter(NN,2,2)*Tmp2(IJK1)   &
                                                 +  DX(1)*WFilter(NN,3,2)*Tmp2(IJK2)
!
                      WCoef(IJK2) = WCoef(IJK2)  + XD2(1)*WFilter(NN,1,3)*Tmp2(IJK0)   &
                                                 +  XD(1)*WFilter(NN,2,3)*Tmp2(IJK1)   &
                                                 +        WFilter(NN,3,3)*Tmp2(IJK2)
                   ENDDO
                ENDDO
             ENDDO
          CASE(1010)
             DO NN=1,LenNOrder
                XF(2) = XX(2)+DX(2)*(2*NN-1-LenNOrder)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp2,Tmp1)
                DO K=0,DOrder
                   DO I=0,DOrder
                      IJK0  = 1+I+DMax*DMax*K
                      IJK1  = IJK0+1*DMax
                      IJK2  = IJK0+2*DMax
                      WCoef(IJK0) = WCoef(IJK0) +        WFilter(NN,1,1)*Tmp2(IJK0)   &
                                                +  DX(2)*WFilter(NN,2,1)*Tmp2(IJK1)   &
                                                + DX2(2)*WFilter(NN,3,1)*Tmp2(IJK2)  
!
                      WCoef(IJK1) = WCoef(IJK1) +  XD(2)*WFilter(NN,1,2)*Tmp2(IJK0)   &
                                                +        WFilter(NN,2,2)*Tmp2(IJK1)   &
                                                +  DX(2)*WFilter(NN,3,2)*Tmp2(IJK2)
!
                      WCoef(IJK2) = WCoef(IJK2) + XD2(2)*WFilter(NN,1,3)*Tmp2(IJK0)   &
                                                +  XD(2)*WFilter(NN,2,3)*Tmp2(IJK1)   &
                                                +        WFilter(NN,3,3)*Tmp2(IJK2)
                   ENDDO
                ENDDO
             ENDDO
          CASE(1001)
            DO NN=1,LenNOrder
                XF(3) = XX(3)+DX(3)*(2*NN-1-LenNOrder)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp2,Tmp1)
                DO J=0,DOrder
                   DO I=0,DOrder
                      IJK0  = 1+I+DMax*J
                      IJK1  = IJK0+1*DMax*DMax
                      IJK2  = IJK0+2*DMax*DMax
                      WCoef(IJK0) = WCoef(IJK0) +        WFilter(NN,1,1)*Tmp2(IJK0)   &
                                                +  DX(3)*WFilter(NN,2,1)*Tmp2(IJK1)   &
                                                + DX2(3)*WFilter(NN,3,1)*Tmp2(IJK2)  
!
                      WCoef(IJK1) = WCoef(IJK1) +  XD(3)*WFilter(NN,1,2)*Tmp2(IJK0)   &
                                                +        WFilter(NN,2,2)*Tmp2(IJK1)   &
                                                +  DX(3)*WFilter(NN,3,2)*Tmp2(IJK2)
!
                      WCoef(IJK2) = WCoef(IJK2) + XD2(3)*WFilter(NN,1,3)*Tmp2(IJK0)   &
                                                +  XD(3)*WFilter(NN,2,3)*Tmp2(IJK1)   &
                                                +        WFilter(NN,3,3)*Tmp2(IJK2)
                   ENDDO
                ENDDO
             ENDDO
          END SELECT
       CASE(3)
          SELECT CASE(NewIType)
          CASE(1100)
             DO NN=1,LenNOrder
                XF(1) = XX(1)+DX(1)*(2*NN-1-LenNOrder)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp2,Tmp1)
                DO K=0,DOrder
                   DO J=0,DOrder
                      IJK0 = 1+DMax*J+DMax*DMax*K
                      IJK1 = IJK0+1
                      IJK2 = IJK0+2
                      IJK3 = IJK0+3
                      WCoef(IJK0) = WCoef(IJK0) +        WFilter(NN,1,1)*Tmp2(IJK0)  &
                                                +  DX(1)*WFilter(NN,2,1)*Tmp2(IJK1)  &
                                                + DX2(1)*WFilter(NN,3,1)*Tmp2(IJK2)  &
                                                + DX3(1)*WFilter(NN,4,1)*Tmp2(IJK3) 
!
                      WCoef(IJK1) = WCoef(IJK1) +  XD(1)*WFilter(NN,1,2)*Tmp2(IJK0)  &
                                                +        WFilter(NN,2,2)*Tmp2(IJK1)  &
                                                +  DX(1)*WFilter(NN,3,2)*Tmp2(IJK2)  &
                                                + DX2(1)*WFilter(NN,4,2)*Tmp2(IJK3)
!
                      WCoef(IJK2) = WCoef(IJK2) + XD2(1)*WFilter(NN,1,3)*Tmp2(IJK0)  &
                                                +  XD(1)*WFilter(NN,2,3)*Tmp2(IJK1)  &
                                                +        WFilter(NN,3,3)*Tmp2(IJK2)  &
                                                +  DX(1)*WFilter(NN,4,3)*Tmp2(IJK3)
!
                      WCoef(IJK3) = WCoef(IJK3) + XD3(1)*WFilter(NN,1,4)*Tmp2(IJK0)  &
                                                + XD2(1)*WFilter(NN,2,4)*Tmp2(IJK1)  &
                                                +  XD(1)*WFilter(NN,3,4)*Tmp2(IJK2)  &
                                                +        WFilter(NN,4,4)*Tmp2(IJK3)
                   ENDDO
                ENDDO
             ENDDO
          CASE(1010)
             DO NN=1,LenNOrder
                XF(2) = XX(2)+DX(2)*(2*NN-1-LenNOrder)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp2,Tmp1)
                DO K=0,DOrder
                   DO I=0,DOrder
                      IJK0 = 1+I+DMax*DMax*K
                      IJK1 = IJK0+1*DMax
                      IJK2 = IJK0+2*DMax
                      IJK3 = IJK0+3*DMax
                      WCoef(IJK0) = WCoef(IJK0) +        WFilter(NN,1,1)*Tmp2(IJK0)  &
                                                +  DX(2)*WFilter(NN,2,1)*Tmp2(IJK1)  &
                                                + DX2(2)*WFilter(NN,3,1)*Tmp2(IJK2)  &
                                                + DX3(2)*WFilter(NN,4,1)*Tmp2(IJK3) 
!
                      WCoef(IJK1) = WCoef(IJK1) +  XD(2)*WFilter(NN,1,2)*Tmp2(IJK0)  &
                                                +        WFilter(NN,2,2)*Tmp2(IJK1)  &
                                                +  DX(2)*WFilter(NN,3,2)*Tmp2(IJK2)  &
                                                + DX2(2)*WFilter(NN,4,2)*Tmp2(IJK3)
!
                      WCoef(IJK2) = WCoef(IJK2) + XD2(2)*WFilter(NN,1,3)*Tmp2(IJK0)  &
                                                +  XD(2)*WFilter(NN,2,3)*Tmp2(IJK1)  &
                                                +        WFilter(NN,3,3)*Tmp2(IJK2)  &
                                                +  DX(2)*WFilter(NN,4,3)*Tmp2(IJK3)
!
                      WCoef(IJK3) = WCoef(IJK3) + XD3(2)*WFilter(NN,1,4)*Tmp2(IJK0)  &
                                                + XD2(2)*WFilter(NN,2,4)*Tmp2(IJK1)  &
                                                +  XD(2)*WFilter(NN,3,4)*Tmp2(IJK2)  &
                                                +        WFilter(NN,4,4)*Tmp2(IJK3)
                   ENDDO
                ENDDO
             ENDDO
          CASE(1001)
             DO NN=1,LenNOrder
                XF(3) = XX(3)+DX(3)*(2*NN-1-LenNOrder)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp2,Tmp1)
                DO J=0,DOrder
                   DO I=0,DOrder
                      IJK0 = 1+I+DMax*J
                      IJK1 = IJK0+1*DMax*DMax
                      IJK2 = IJK0+2*DMax*DMax
                      IJK3 = IJK0+3*DMax*DMax
                      WCoef(IJK0) = WCoef(IJK0) +        WFilter(NN,1,1)*Tmp2(IJK0)  &
                                                +  DX(3)*WFilter(NN,2,1)*Tmp2(IJK1)  &
                                                + DX2(3)*WFilter(NN,3,1)*Tmp2(IJK2)  &
                                                + DX3(3)*WFilter(NN,4,1)*Tmp2(IJK3) 
!
                      WCoef(IJK1) = WCoef(IJK1) +  XD(3)*WFilter(NN,1,2)*Tmp2(IJK0)  &
                                                +        WFilter(NN,2,2)*Tmp2(IJK1)  &
                                                +  DX(3)*WFilter(NN,3,2)*Tmp2(IJK2)  &
                                                + DX2(3)*WFilter(NN,4,2)*Tmp2(IJK3)
!
                      WCoef(IJK2) = WCoef(IJK2) + XD2(3)*WFilter(NN,1,3)*Tmp2(IJK0)  &
                                                +  XD(3)*WFilter(NN,2,3)*Tmp2(IJK1)  &
                                                +        WFilter(NN,3,3)*Tmp2(IJK2)  &
                                                +  DX(3)*WFilter(NN,4,3)*Tmp2(IJK3)
!
                      WCoef(IJK3) = WCoef(IJK3) + XD3(3)*WFilter(NN,1,4)*Tmp2(IJK0)  &
                                                + XD2(3)*WFilter(NN,2,4)*Tmp2(IJK1)  &
                                                +  XD(3)*WFilter(NN,3,4)*Tmp2(IJK2)  &
                                                +        WFilter(NN,4,4)*Tmp2(IJK3)
                   ENDDO
                ENDDO
             ENDDO
          END SELECT
       END SELECT
       CALL NormalizeWCoefs(WCoef,NewLevel,NewIType)
!               
   END SUBROUTINE ComputeWCoef
!==================================================================
!  Wrapper For HasThisNode
!==================================================================
   FUNCTION HasNode(XX)
     LOGICAL                    :: HasNode
     REAL(DOUBLE),DIMENSION(3)  :: XX
!
     IF(ABS(XX(1)-GBox%Center(1))>GBox%Half(1)) THEN
        HasNode = .TRUE.
        RhoSum  = Zero
        RETURN
     ELSEIF(ABS(XX(2)-GBox%Center(2))>GBox%Half(2)) THEN
        HasNode = .TRUE.
        RhoSum  = Zero
        RETURN
     ELSEIF(ABS(XX(3)-GBox%Center(3))>GBox%Half(3)) THEN
        HasNode = .TRUE.
        RhoSum  = Zero
        RETURN
     ELSE
        IsTrue = .FALSE.
        Xpos   = XX
        RhoSum = Zero
        CALL HasThisNode(PIGWRoot)
        HasNode= IsTrue
     ENDIF
!
   END FUNCTION HasNode
!==================================================================
!  Determine if PIGWTree has this Node
!==================================================================
   RECURSIVE SUBROUTINE HasThisNode(Node)
     TYPE(PIGWNode), POINTER      :: Node
     INTEGER                      :: ILink
     REAL(DOUBLE)                 :: RQx,RQy,RQz
     REAL(DOUBLE),PARAMETER       :: DTol=1.D-10
!
     IF(IsTrue) RETURN
!
     RQx=Node%Box%Center(1)-Xpos(1)
     IF(ABS(RQx)>Node%Box%Half(1)) RETURN
     RQy=Node%Box%Center(2)-Xpos(2)
     IF(ABS(RQy)>Node%Box%Half(2)) RETURN
     RQz=Node%Box%Center(3)-Xpos(3)
     IF(ABS(RQz)>Node%Box%Half(3)) RETURN
!
     IF(ABS(RQx) < DTol .AND. ABS(RQy) < DTol .AND. ABS(RQz) < DTol) THEN
        IsTrue = .TRUE.
        RhoSum = Node%WCoef
        RETURN
     ENDIF
!
     IF(Node%NLinks==0) RETURN
!
     DO ILink=1,Node%NLinks
        CALL HasThisNode(Node%Links(ILink)%Link)
        IF(IsTrue) EXIT
     ENDDO
!
   END SUBROUTINE HasThisNode

!=================================================================================
!  Wrapper for RhoAtPoint
!=================================================================================
   FUNCTION RhoAtX(XX) 
     REAL(DOUBLE),DIMENSION(LenDOrder) :: RhoAtX
     REAL(DOUBLE),DIMENSION(3)         :: XX
     REAL(DOUBLE)                      :: Bet,Fun,D1X,D1Y,D1Z,D2X,D2Y,D2Z
!
     IF(ABS(XX(1)-GBox%Center(1))>GBox%Half(1)) THEN
        RhoAtX = Zero
        RETURN
     ELSEIF(ABS(XX(2)-GBox%Center(2))>GBox%Half(2)) THEN
        RhoAtX = Zero
        RETURN
     ELSEIF(ABS(XX(3)-GBox%Center(3))>GBox%Half(3)) THEN
        RhoAtX = Zero
        RETURN
     ELSE
        NCallToRho = NCallToRho+1
!
        RhoSum = Zero        
        Xpos   = XX
!        CALL TestRho()
        CALL RhoAtPointX(RhoRoot)
        RhoAtX = RhoSum
     ENDIF

   END FUNCTION RhoAtX
!=================================================================================
!  Test Function Integration
!=================================================================================
   SUBROUTINE TestRho()
     REAL(DOUBLE)                    :: R2,Expt,Beta
     INTEGER                         :: I,J,K,IJK
     REAL(DOUBLE),DIMENSION(0:3)     :: DivX,DivY,DivZ
!
     Beta = One
     R2   = Xpos(1)*Xpos(1)+Xpos(2)*Xpos(2)+Xpos(3)*Xpos(3)
     Expt = EXP(-Beta*R2)
!
     DivX(0) = One
     DivY(0) = One
     DivZ(0) = One
!     
     DivX(1) = -Two*Beta*Xpos(1)
     DivY(1) = -Two*Beta*Xpos(2)
     DivZ(1) = -Two*Beta*Xpos(3)
!     
     DivX(2) = -Two*Beta+(2*Beta*Xpos(1))**2
     DivY(2) = -Two*Beta+(2*Beta*Xpos(2))**2
     DivZ(2) = -Two*Beta+(2*Beta*Xpos(3))**2
!     
     DivX(3) = 12.D0*Beta*Beta*Xpos(1)-(Two*Beta*Xpos(1))**3
     DivY(3) = 12.D0*Beta*Beta*Xpos(2)-(Two*Beta*Xpos(2))**3
     DivZ(3) = 12.D0*Beta*Beta*Xpos(3)-(Two*Beta*Xpos(3))**3
!
     DO K=0,DOrder
        DO J=0,DOrder
           DO I=0,DOrder
              IJK = 1+I+DMax*J+DMax*DMax*K
              RhoSum(IJK) = RhoSum(IJK)+DivX(I)*DivY(J)*DivZ(K)*Expt
           ENDDO
        ENDDO
     ENDDO
!
   END SUBROUTINE TestRho
!=================================================================================
!  Sums the significant contributions (leaves) to the density at a point
!=================================================================================
   RECURSIVE SUBROUTINE RhoAtPointX(Node)
     TYPE(RhoNode), POINTER                     :: Node
     INTEGER                                    :: L,M,N,LMN,L1,L2,J
     REAL(DOUBLE)                               :: RQx,RQy,RQz,RQ2
     REAL(DOUBLE)                               :: Xpt,Zeta,TwoZeta,X
     REAL(DOUBLE),DIMENSION(1:142)              :: W
!----------------------------------------------------------------------------------
     RQx=Xpos(1)-Node%Box%Center(1)
     IF(ABS(RQx)>Node%Box%Half(1)) RETURN
     RQy=Xpos(2)-Node%Box%Center(2)
     IF(ABS(RQy)>Node%Box%Half(2)) RETURN
     RQz=Xpos(3)-Node%Box%Center(3)
     IF(ABS(RQz)>Node%Box%Half(3)) RETURN
!    
     IF(Node%Leaf)THEN   
        RQ2     = RQx*RQx+RQy*RQy+RQz*RQz
        Zeta    = Node%Zeta
        TwoZeta = Two*Zeta
        X       = Node%Zeta*RQ2
        IF(X<Exp_Switch)THEN
           SELECT CASE(DOrder)
           CASE(0)
              J=AINT(X*Exp_Grid)
              Xpt=Exp_0(J)+X*(Exp_1(J)+X*(Exp_2(J)+X*(Exp_3(J)+X*Exp_4(J))))
              SELECT CASE(Node%Ell)
              CASE(0)
!                 INCLUDE "MMA/DivRho_0_0.Inc"    
              CASE(1)
!                 INCLUDE "MMA/DivRho_1_0.Inc"
              CASE(2)
!                 INCLUDE "MMA/DivRho_2_0.Inc"
              CASE(3)
!                 INCLUDE "MMA/DivRho_3_0.Inc"
              CASE(4)
!                 INCLUDE "MMA/DivRho_4_0.Inc"
              END SELECT
           CASE(1)
              J=AINT(X*Exp_Grid)
              Xpt=Exp_0(J)+X*(Exp_1(J)+X*(Exp_2(J)+X*(Exp_3(J)+X*Exp_4(J))))
              SELECT CASE(Node%Ell)
              CASE(0)
!                 INCLUDE "MMA/DivRho_0_1.Inc"
              CASE(1)
!                 INCLUDE "MMA/DivRho_1_1.Inc"
              CASE(2)
!                 INCLUDE "MMA/DivRho_2_1.Inc"
              END SELECT
           CASE(2)
              J=AINT(X*Exp_Grid)
              Xpt=Exp_0(J)+X*(Exp_1(J)+X*(Exp_2(J)+X*(Exp_3(J)+X*Exp_4(J))))
              SELECT CASE(Node%Ell)
              CASE(0)
!                 INCLUDE "MMA/DivRho_0_2.Inc"
              CASE(1)
!                 INCLUDE "MMA/DivRho_1_2.Inc"
              CASE(2)
!                 INCLUDE "MMA/DivRho_2_2.Inc"
              END SELECT
           CASE(3)
              J=AINT(X*Exp_Grid)
              Xpt=Exp_0(J)+X*(Exp_1(J)+X*(Exp_2(J)+X*(Exp_3(J)+X*Exp_4(J))))
              SELECT CASE(Node%Ell)
              CASE(0)
!                 INCLUDE "MMA/DivRho_0_3.Inc"
              CASE(1)
!                 INCLUDE "MMA/DivRho_1_3.Inc"
              CASE(2)
!                 INCLUDE "MMA/DivRho_2_3.Inc"
              END SELECT
           END SELECT
        ENDIF
     ELSE
        CALL RhoAtPointX(Node%Descend)
        CALL RhoAtPointX(Node%Descend%Travrse)
     ENDIF
   END SUBROUTINE RhoAtPointX
!====================================================================================
!  Calculate Exc 
!====================================================================================
   SUBROUTINE ExcVxc(RhoXX,Exc,Vxc)
     REAL(DOUBLE),DIMENSION(LenDOrder) :: RhoXX,Exc,Vxc
     REAL(DOUBLE)                      :: InvRho,d0Exc,d1Exc,d2Exc,d3Exc,d10Exc
     REAL(DOUBLE)                      :: d4Exc,d5Exc,d6Exc,d7Exc,d8Exc,d9Exc
     REAL(DOUBLE)                      :: d0Vxc,d1Vxc,d2Vxc,d3Vxc
     REAL(DOUBLE)                      :: d4Vxc,d5Vxc,d6Vxc,d7Vxc,d8Vxc,d9Vxc
     REAL(DOUBLE)                      :: RhoTmp1,RhoTmp2,RhoTmp3,RhoTmp4,RhoTmp5
     REAL(DOUBLE)                      :: RhoTmp6,RhoTmp7,RhoTmp8,RhoTmp9
     REAL(DOUBLE),DIMENSION(1:179)     :: W
     REAL(DOUBLE),PARAMETER            :: D1A  =  4.D0/3.D0
     REAL(DOUBLE),PARAMETER            :: D2A  =  1.D0/3.D0
     REAL(DOUBLE),PARAMETER            :: D3A  = -2.D0/3.D0
     REAL(DOUBLE),PARAMETER            :: D4A  = -5.D0/3.D0
     REAL(DOUBLE),PARAMETER            :: D5A  = -8.D0/3.D0
     REAL(DOUBLE),PARAMETER            :: D6A  = -11.D0/3.D0
     REAL(DOUBLE),PARAMETER            :: D7A  = -14.D0/3.D0
     REAL(DOUBLE),PARAMETER            :: D8A  = -17.D0/3.D0
     REAL(DOUBLE),PARAMETER            :: D9A  = -20.D0/3.D0
     REAL(DOUBLE),PARAMETER            :: D10A = -23.D0/3.D0
     CHARACTER(LEN=20)                 :: ExchangeType
!
!!$     Exc = RhoXX
!!$     Vxc = Zero
!!$     RETURN
!
     ExchangeType = 'SlaterDirac'
     SELECT CASE(ExchangeType) 
!    Slater-Dirac exchange
     CASE('SlaterDirac')
        SELECT CASE(DOrder)
        CASE(0)
           IF(RhoXX(1) .LE. Zero) THEN
              Exc = Zero
              Vxc = Zero
              RETURN
           ELSE
              InvRho =  One/RhoXX(1)
              d0Exc  = -1.8610514726982D0*RhoXX(1)**(4.D0/3.D0)
              d1Exc  =  D1A*d0Exc*InvRho
              d0Vxc  =  d1Exc
           ENDIF
           Exc(1) = d0Exc
           Vxc(1) = d0Vxc
        CASE(1)
           IF(RhoXX(1) .LE. Zero) THEN
              Exc = Zero
              Vxc = Zero
              RETURN
           ELSE
              InvRho =  One/RhoXX(1)
              d0Exc  = -1.8610514726982D0*RhoXX(1)**(4.D0/3.D0) 
              d1Exc  =  D1A*d0Exc*InvRho
              d2Exc  =  D2A*d1Exc*InvRho 
              d3Exc  =  D3A*d2Exc*InvRho 
              d4Exc  =  D4A*d3Exc*InvRho 
              d0Vxc  =  d1Exc
              d1Vxc  =  d2Exc
              d2Vxc  =  d3Exc
              d3Vxc  =  d4Exc
           ENDIF
           INCLUDE "MMA/RhoToExc_1.Inc"
        CASE(2) 
           IF(RhoXX(1) .LE. Zero) THEN
              Exc = Zero
              Vxc = Zero
              RETURN
           ELSE
              InvRho =  One/RhoXX(1)
              d0Exc  = -1.8610514726982D0*RhoXX(1)**(4.D0/3.D0)
              d1Exc  =  D1A*d0Exc*InvRho
              d2Exc  =  D2A*d1Exc*InvRho 
              d3Exc  =  D3A*d2Exc*InvRho 
              d4Exc  =  D4A*d3Exc*InvRho 
              d5Exc  =  D5A*d4Exc*InvRho 
              d6Exc  =  D6A*d5Exc*InvRho 
              d7Exc  =  D7A*d6Exc*InvRho 
              d0Vxc  =  d1Exc
              d1Vxc  =  d2Exc
              d2Vxc  =  d3Exc
              d3Vxc  =  d4Exc
              d4Vxc  =  d5Exc
              d5Vxc  =  d6Exc
              d6Vxc  =  d7Exc
           ENDIF
           INCLUDE "MMA/RhoToExc_2.Inc"
        CASE(3) 
           IF(RhoXX(1) .LE. Zero) THEN
              Exc = Zero
              Vxc = Zero
              RETURN
           ELSE
              InvRho =  One/RhoXX(1)
              d0Exc  = -1.8610514726982D0*RhoXX(1)**(4.D0/3.D0)
              d1Exc  =  D1A*d0Exc*InvRho
              d2Exc  =  D2A*d1Exc*InvRho 
              d3Exc  =  D3A*d2Exc*InvRho 
              d4Exc  =  D4A*d3Exc*InvRho 
              d5Exc  =  D5A*d4Exc*InvRho 
              d6Exc  =  D6A*d5Exc*InvRho 
              d7Exc  =  D7A*d6Exc*InvRho 
              d8Exc  =  D8A*d7Exc*InvRho 
              d9Exc  =  D9A*d8Exc*InvRho 
              d10Exc =  D10A*d9Exc*InvRho 
              d0Vxc  =  d1Exc
              d1Vxc  =  d2Exc
              d2Vxc  =  d3Exc
              d3Vxc  =  d4Exc
              d4Vxc  =  d5Exc
              d5Vxc  =  d6Exc
              d6Vxc  =  d7Exc
              d7Vxc  =  d8Exc
              d8Vxc  =  d9Exc
              d9Vxc  =  d10Exc
           ENDIF
!           INCLUDE "MMA/RhoToExc_3.Inc"
        END SELECT
     END SELECT
!
   END SUBROUTINE ExcVxc
!====================================================================================
!  Calculate Exc and Vxc
!====================================================================================
   RECURSIVE SUBROUTINE ComputeExc(Node)
     TYPE(PIGWNode), POINTER                    :: Node
     REAL(DOUBLE)                               :: NodeEnergy
     INTEGER                                    :: ILink,I
!
     NodeEnergy = Zero
     DO I=1,LenDOrder
        NodeEnergy = NodeEnergy + Node%WCoef(I)*Node%IntFactor(I)
     ENDDO
     XCEnergy = XCEnergy + NodeEnergy
!
!     WRITE(*,*) Node%Box%Center(1:3)
!     WRITE(*,*) NodeEnergy
!
     IF(Node%LeafType=='EnddLeaf' .OR. Node%LeafType=='LeafLeaf' ) RETURN
     DO ILink=1,Node%NLinks 
        CALL ComputeExc(Node%Links(ILink)%Link)
     ENDDO
!
   END SUBROUTINE ComputeExc
!
END MODULE DIPMWTree
