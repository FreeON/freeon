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
!     Potistion and Dimenstion of Wavelte
      REAL(DOUBLE), DIMENSION(3)         :: X
      REAL(DOUBLE), DIMENSION(3)         :: DX
!     Bounding box of All and Below
      TYPE(BBox)                         :: Box    
!     The WCoefs    
      REAL(DOUBLE),ALLOCATABLE           :: WCoef(:)
!     For Integration
      REAL(DOUBLE),ALLOCATABLE           :: IntFactor(:)
   END TYPE PIGWNode
!==================================================================
!  Globals
!==================================================================
   TYPE(PIGWNode), POINTER               :: PIGWRoot 
   INTEGER                               :: CurrentLevel,NCallToRho
   INTEGER,PARAMETER                     :: BigJ=20
   INTEGER,DIMENSION(0:BigJ)             :: PIGWNodes
   REAL(DOUBLE),DIMENSION(1:3)           :: X0,Length
   REAL(DOUBLE),DIMENSION(0:BigJ)        :: JFactor,InvJFactor,ScaleX,ScaleY,ScaleZ
!------------------------------------------------------------------
!  Wavelet Order Stuff
!
   INTEGER,PARAMETER                     :: DOrder=3
   INTEGER,PARAMETER                     :: DMax=DOrder+1
   INTEGER,PARAMETER                     :: LenDOrder=(DOrder+1)**3
   REAL(DOUBLE),DIMENSION(LenDOrder)     :: RhoTmp
   REAL(DOUBLE),DIMENSION(DMax)          :: WvltInt
   REAL(DOUBLE),DIMENSION(2,DMax,DMax)   :: WFilter
!------------------------------------------------------------------
!  RhoAtPoint and HasNode Stuff
!
   LOGICAL                               :: IsTrue
   REAL(DOUBLE)                          :: VxcSum 
   REAL(DOUBLE),DIMENSION(LenDOrder)     :: RhoSum
   REAL(DOUBLE),DIMENSION(1:3)           :: Xpos
   TYPE(BBox)                            :: RhoBox,XCBox       
!------------------------------------------------------------------
!  Exc Energy
!
   REAL(DOUBLE)                          :: XCEnergy
   CONTAINS 
!==================================================================
!  Generate the Wavelet Rep of the XC potential
!==================================================================
   SUBROUTINE DIPMWTree(MaxJ)
     INTEGER                             :: MaxJ,J,I,K,TotWav,NN,MaxCurrent
     REAL(DOUBLE),DIMENSION(3)           :: X,XX
     REAL(DOUBLE)                        :: TauDIPMW_old,MinDelta,DX,DY,DZ
     REAL(DOUBLE)                        :: FFac,SUM,RhoInt,Tmp1,Tmp2
     CHARACTER(LEN=DEFAULT_CHR_LEN)      :: Mssg 
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
!    Rho Bounding Box
     RhoBox%Center(:)   = RhoRoot%Box%Center(:)
     RhoBox%Half(:)     = RhoRoot%Box%Half(:)
     RhoBox%BndBox(:,1) = RhoRoot%Box%BndBox(:,1)+1.D-12
     RhoBox%BndBox(:,2) = RhoRoot%Box%BndBox(:,2)-1.D-12
!    ExchangeCorrelation Regions Bounding Box
     XCBox%Center(:)    = RhoRoot%Box%Center(:)
     XCBox%Half(:)      = RhoRoot%Box%Half(:)
     XCBox%BndBox(:,1)  = RhoRoot%Box%BndBox(:,1)+1.D-12
     XCBox%BndBox(:,2)  = RhoRoot%Box%BndBox(:,2)-1.D-12
     CALL MakeBoxPeriodic(XCBox)
!    Initialize the lengths  
     Length(:)        = ABS(XCBox%BndBox(:,2)-XCBox%BndBox(:,1))
     X0(:)            = Half*(XCBox%BndBox(:,2)+XCBox%BndBox(:,1))
!    Initialize JFactor
     DO J=0,BigJ
        JFactor(J)    = Two**(J+1)
        InvJFactor(J) = Half**(J+1)  
        ScaleX(J)     = InvJFactor(J)*Length(1)
        ScaleY(J)     = InvJFactor(J)*Length(2)
        ScaleZ(J)     = InvJFactor(J)*Length(3)  
     ENDDO
!    Initialize the Root Node and XC Energy
     NCallToRho= 0
     XCEnergy  = Zero
     CALL InitializePIGWRoot(X0,0)
!    Recursively generate the Wavelet Rep of the Density, XC Energy and XC potential
     TotWav  = PIGWNodes(0)
     DO CurrentLevel=1,MaxJ
        CALL MakePIGWTree(PIGWRoot)
        IF(PIGWNodes(CurrentLevel) == 0) EXIT
        TotWav = TotWav+PIGWNodes(CurrentLevel)
!
        Mssg=ProcessName('DIPMW.TreeGen')                             &
                  // 'Level = ' // TRIM(IntToChar(CurrentLevel))      &
                  //', Tau = ' //TRIM(DblToShrtChar(TauDIPMW))        &
                  //', <Exc> = '//TRIM(DblToMedmChar(XCEnergy))       &
                  //', DipMW/Atom  = '//TRIM(IntToChar(TotWav/GM%Natms)) &
                  //', DipMW/Level = '//TRIM(IntToChar(PIGWNodes(CurrentLevel)/GM%Natms)) 
        CALL OpenASCII(OutFile,Out)         
        WRITE(*,*)TRIM(Mssg)
        WRITE(Out,*)TRIM(Mssg)
        CLOSE(Out)
     ENDDO
!
     WRITE(*,*) "Total Wavelets        = ",TotWav
     WRITE(*,*) "Total Call to RhoTree = ",NCallToRho
!
!    Compute Pointwise the Vxc for reference
!!$     DO I=0,256
!!$        X(1) =  XCBox%BndBox(1,1)+I*Length(1)/256.D0
!!$        X(2) =  XCBox%BndBox(2,1)+Length(2)/2.0D0
!!$        X(3) =  XCBox%BndBox(3,1)+Length(3)/2.0D0
!!$        RhoTmp = RhoAtX(X)
!!$        Tmp1   = VxcAtX(X)
!!$        Tmp2   = -1.8610514726982D0*(2.D0/3.D0)*(RhoTmp(1))**(1.D0/3.D0)
!!$        WRITE(99,'4(D16.8,3x)') X(1),Tmp1,Tmp2,Tmp2-Tmp1
!!$     ENDDO
!!$
!!$     IF(.TRUE.) STOP
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
        WvltInt(1) = 5.0D-1
!
        WFilter(1,1,1) = -5.0D-1
        WFilter(2,1,1) = -5.0D-1
     CASE(1)
        WvltInt(1) = 5.0000000000000000000D-1
        WvltInt(2) = 8.3333333333333333333D-2
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
     CASE(2)
 
        WvltInt(1) = 5.0000000000000000000D-1
        WvltInt(2) = 1.0000000000000000000D-1
        WvltInt(3) = 8.3333333333333333333D-3
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
     CASE(3)
        WvltInt(1) = 5.0000000000000000000D-1
        WvltInt(2) = 1.0714285714285714286D-1
        WvltInt(3) = 1.1904761904761904762D-2
        WvltInt(4) = 5.9523809523809523809D-4
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
        CALL Halt('WFilter DOrder>3 not defined')
     END SELECT
!
   END SUBROUTINE InitializeWFilters
!==================================================================
!  Initialize PIGWTree
!==================================================================
   SUBROUTINE InitializePIGWRoot(X,ILevel)
     INTEGER                            :: ILevel,Status,ILink,IX,IY,IZ,JLink,I
     REAL(DOUBLE),DIMENSION(LenDOrder)  :: RhoX,ExcX,VxcX
     REAL(DOUBLE),DIMENSION(3)          :: X,XX,DX,XL,XH
     LOGICAL                            :: DoNode
!        
!    PIGWRoot
!
!    Allocate the Root Node: There is no Inforamtion on this Node
     CALL NewPIGWNode(PIGWRoot,ILevel,1000,X)
     PIGWRoot%LeafType ='NodeNode'
     DX                = PIGWRoot%DX
!    Bounding Box of Root Node
     XL = X-DX
     XH = X+DX
     CALL NodeBoundingBox(XL,XH,PIGWRoot%Box)
!    Allocate the Links
     PIGWRoot%NLinks  = 3*3*3
     ALLOCATE(PIGWRoot%Links(PIGWRoot%NLinks),STAT=Status)
     IF(Status/=SUCCEED)CALL Halt('Link ALLOCATE failed in InitializePIGWRoot ')
!
!    Surface Scaling Functions 
!
     JLink = 0
     DO IX =-1,1
        DO IY = -1,1
           DO IZ =-1,1
              JLink = JLink+1
              XX(1) = X(1)+IX*DX(1)
              XX(2) = X(2)+IY*DX(2)
              XX(3) = X(3)+IZ*DX(3)
!             Allocate Node
              CALL NewPIGWNode(PIGWRoot%Links(JLink)%Link,ILevel,1000,XX)
!             Bounding Box of Root Node
              XL = XX-DX
              XH = XX+DX
              CALL NodeBoundingBox(XL,XH,PIGWRoot%Links(JLink)%Link%Box)
!             Integration Weight
              CALL WWeight(XX,DX,PIGWRoot%Links(JLink)%Link%Box,PIGWRoot%Links(JLink)%Link%IntFactor)
              RhoX = RhoAtX(XX)
              CALL ExcVxc(RhoX,ExcX,VxcX)
              CALL NormalizeWCoefs(VxcX,DX(1),DX(2),DX(3))
              CALL NormalizeWCoefs(ExcX,DX(1),DX(2),DX(3))
              PIGWRoot%Links(JLink)%Link%WCoef = VxcX
!             Compute Initial Exc
              DO I=1,LenDOrder
                 XCEnergy = XCEnergy + PIGWRoot%Links(JLink)%Link%IntFactor(I)*ExcX(I)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
!
   END SUBROUTINE InitializePIGWRoot
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
!    Intialize Node Info
     Node%LeafType        = 'LeafLeaf'
     Node%Level           = ILevel
     Node%IType           = IType
     Node%NLinks          = 0
!    Potition of Wavelet and Dimension
     Node%X               = X
     CALL ComputeDX(ILevel,IType,Node%DX)
!    Allocate Wavelet Coefs
     ALLOCATE(Node%WCoef(1:LenDOrder),STAT=Status)
     IF(Status/=SUCCEED)CALL Halt(' Node%WCoef ALLOCATE failed in NewPIGWNode ')
!    Allocate Integration Factor
     ALLOCATE(Node%IntFactor(1:LenDOrder),STAT=Status)
     IF(Status/=SUCCEED)CALL Halt(' Node%IntFactor ALLOCATE failed in NewPIGWNode ')
!    Intitialize WCoefs and IntFactor
     Node%WCoef(:)        = Zero
     Node%IntFactor(:)    = Zero
!    Initialize Box
     Node%Box%Center      = Zero
     Node%Box%Half        = Zero
     Node%Box%BndBox(:,1) = Zero
     Node%Box%BndBox(:,2) = Zero
!    Bump Globals Counter
     PIGWNodes(ILevel) = PIGWNodes(ILevel)+1
!
   END SUBROUTINE NewPIGWNode
!======================================================================================
!
!======================================================================================
   SUBROUTINE ComputeDX(ILevel,IType,DX)
     INTEGER                     :: ILevel,IType
     REAL(DOUBLE),DIMENSION(3)   :: DX
!    Intialize Dimensions
     IF(IType==1000) THEN
        DX(1)     = 1.0D0*ScaleX(ILevel)
        DX(2)     = 1.0D0*ScaleY(ILevel)
        DX(3)     = 1.0D0*ScaleZ(ILevel)
     ELSEIF(IType==1100) THEN
        DX(1)     = 1.0D0*ScaleX(ILevel)
        DX(2)     = 2.0D0*ScaleY(ILevel)
        DX(3)     = 2.0D0*ScaleZ(ILevel)
     ELSEIF(IType==1010) THEN
        DX(1)     = 1.0D0*ScaleX(ILevel)
        DX(2)     = 1.0D0*ScaleY(ILevel)
        DX(3)     = 2.0D0*ScaleZ(ILevel)
     ELSEIF(IType==1001) THEN
        DX(1)     = 1.0D0*ScaleX(ILevel)
        DX(2)     = 1.0D0*ScaleY(ILevel)
        DX(3)     = 1.0D0*ScaleZ(ILevel)
     ENDIF
   END SUBROUTINE ComputeDX
!======================================================================================
!
!======================================================================================
   SUBROUTINE NodeBoundingBox(XL,XH,Box)
     TYPE(BBox)                      :: Box     
     INTEGER                         :: Level,IType
     REAL(DOUBLE),DIMENSION(3)       :: XL,XH
!    Bounding Box
     Box%BndBox(:,1) = XL
     Box%BndBox(:,2) = XH
!    Box Extents
     Box%BndBox(1,1) = MAX(Box%BndBox(1,1),XCBox%BndBox(1,1))
     Box%BndBox(1,2) = MIN(Box%BndBox(1,2),XCBox%BndBox(1,2))
!
     Box%BndBox(2,1) = MAX(Box%BndBox(2,1),XCBox%BndBox(2,1))
     Box%BndBox(2,2) = MIN(Box%BndBox(2,2),XCBox%BndBox(2,2))
!
     Box%BndBox(3,1) = MAX(Box%BndBox(3,1),XCBox%BndBox(3,1))
     Box%BndBox(3,2) = MIN(Box%BndBox(3,2),XCBox%BndBox(3,2))
!    Center and Half
     Box%Center(:) = Half*(Box%BndBox(:,2)+Box%BndBox(:,1))
     Box%Half(:)   = Half*(Box%BndBox(:,2)-Box%BndBox(:,1))
!
   END SUBROUTINE NodeBoundingBox
!======================================================================================
!
!======================================================================================
   SUBROUTINE WWeight(X,DX,Box,IWeight)
     TYPE(BBox)                        :: Box   
     REAL(DOUBLE),DIMENSION(3)         :: X,DX  
     INTEGER                           :: I,J,K,IJK
     REAL(DOUBLE)                      :: Lx,Ly,Lz,Rx,Ry,Rz
     REAL(DOUBLE)                      :: WX,WY,WZ
     REAL(DOUBLE),DIMENSION(LenDOrder) :: IWeight

!
     Lx = MIN(X(1)-Box%BndBox(1,1),DX(1))
     Ly = MIN(X(2)-Box%BndBox(2,1),DX(2))
     Lz = MIN(X(3)-Box%BndBox(3,1),DX(3))
     Rx = MIN(Box%BndBox(1,2)-X(1),DX(1))
     Ry = MIN(Box%BndBox(2,2)-X(2),DX(2))
     Rz = MIN(Box%BndBox(3,2)-X(3),DX(3))
!
     IWeight  = Zero
     DO K=0,DOrder
        DO J=0,DOrder
           DO I=0,DOrder
              WX  = ((-One)**I)*Lx*WvltInt(I+1)+Rx*WvltInt(I+1)
              WY  = ((-One)**J)*Ly*WvltInt(J+1)+Ry*WvltInt(J+1)
              WZ  = ((-One)**K)*Lz*WvltInt(K+1)+Rz*WvltInt(K+1)
              IJK = 1+I+DMax*J+DMax*DMax*K
              IWeight(IJK) = WX*WY*WZ
           ENDDO 
        ENDDO
     ENDDO
!
   END SUBROUTINE WWeight
!======================================================================================
!
!======================================================================================
   SUBROUTINE MaxWWeight(X,DX,Box,IWeight)
     TYPE(BBox)                        :: Box   
     REAL(DOUBLE),DIMENSION(3)         :: X,DX  
     INTEGER                           :: I,J,K,IJK
     REAL(DOUBLE)                      :: Lx,Ly,Lz,Rx,Ry,Rz
     REAL(DOUBLE)                      :: WX,WY,WZ
     REAL(DOUBLE),DIMENSION(LenDOrder) :: IWeight

!
     Lx = MIN(X(1)-Box%BndBox(1,1),DX(1))
     Ly = MIN(X(2)-Box%BndBox(2,1),DX(2))
     Lz = MIN(X(3)-Box%BndBox(3,1),DX(3))
     Rx = MIN(Box%BndBox(1,2)-X(1),DX(1))
     Ry = MIN(Box%BndBox(2,2)-X(2),DX(2))
     Rz = MIN(Box%BndBox(3,2)-X(3),DX(3))
!
     IWeight  = Zero
     DO K=0,DOrder
        DO J=0,DOrder
           DO I=0,DOrder
              WX  = Lx*WvltInt(I+1)+Rx*WvltInt(I+1)
              WY  = Ly*WvltInt(J+1)+Ry*WvltInt(J+1)
              WZ  = Lz*WvltInt(K+1)+Rz*WvltInt(K+1)
              IJK = 1+I+DMax*J+DMax*DMax*K
              IWeight(IJK) = WX*WY*WZ
           ENDDO 
        ENDDO
     ENDDO
!
   END SUBROUTINE MaxWWeight
!==================================================================
!  Recursively Genergate the Wavelet Rep of the XC potential
!==================================================================
   RECURSIVE SUBROUTINE MakePIGWTree(Node)
     TYPE(PIGWNode), POINTER               :: Node
     TYPE(BBox)                            :: NewBox     
     REAL(DOUBLE)                          :: TotalWeight
     REAL(DOUBLE),DIMENSION(3)             :: DX,XX,DXX,XL,XH
     INTEGER                               :: NewLevel,IX,IY,IZ,NewIType,ILink,Status,NLinks,I,J,K
     REAL(DOUBLE),DIMENSION(LenDOrder)     :: WECoef,WVCoef,IntValue
!
     INTEGER,DIMENSION(26)                 :: ITypeArray,ILeafType
     REAL(DOUBLE),DIMENSION(26,3)          :: XposArray
     REAL(DOUBLE),DIMENSION(26,LenDOrder)  :: WCoefArray,IntVArray
     REAL(DOUBLE),PARAMETER                :: Mix=0.01D0
!
     IF(Node%LeafType=='LeafLeaf') THEN 
        NewLevel = Node%Level+1
        DX(1)    = ScaleX(NewLevel)
        DX(2)    = ScaleY(NewLevel)
        DX(3)    = ScaleZ(NewLevel)
!       Determine the Sub-Wavelets that will exist
        NLinks = 0
        DO IX=-1,1
        DO IY=-1,1
        DO IZ=-1,1
           IF(IX .NE. 0) NewIType = 1100  
           IF(IY .NE. 0) NewIType = 1010  
           IF(IZ .NE. 0) NewIType = 1001  
           XX(1) = Node%X(1)+IX*DX(1)
           XX(2) = Node%X(2)+IY*DX(2)
           XX(3) = Node%X(3)+IZ*DX(3)
           IF(.NOT. HasNode(XX)) THEN
!             Compute the Wavelet Coefs
              CALL ComputeWCoef(XX,NewLevel,NewIType,WECoef,WVCoef)
              CALL ComputeDX(NewLevel,NewIType,DXX)
!             Compute the Bounding Box
              XL = XX-DXX
              XH = XX+DXX
              CALL NodeBoundingBox(XL,XH,NewBox)
!             Compute The MaxWeight 
              CALL MaxWWeight(XX,DXX,NewBox,IntValue)
              TotalWeight = Zero
              DO I=1,LenDOrder
                 TotalWeight = TotalWeight+IntValue(I)*((1.D0-Mix)*ABS(WECoef(I))+Mix*ABS(WVCoef(I)))
              ENDDO
!             Compute the Energy Contribution to this node
              CALL WWeight(XX,DXX,NewBox,IntValue)
              DO I=1,LenDOrder
                 XCEnergy = XCEnergy + IntValue(I)*WECoef(I)
              ENDDO
!             Do the Links
              IF(TotalWeight > TauDIPMW) THEN
                 NLinks = NLinks+1
                 ILeafType(NLinks)              = 1
                 ITypeArray(NLinks)             = NewIType
                 XposArray(NLinks,1:3)          = XX(1:3)
                 WCoefArray(NLinks,1:LenDOrder) = WVCoef(1:LenDOrder)
                 IntVArray(NLinks,1:LenDOrder)  = IntValue(1:LenDOrder)
              ELSE
                 NLinks = NLinks+1
                 ILeafType(NLinks)              = 0
                 ITypeArray(NLinks)             = NewIType
                 XposArray(NLinks,1:3)          = XX(1:3)
                 WCoefArray(NLinks,1:LenDOrder) = WVCoef(1:LenDOrder)
                 IntVArray(NLinks,1:LenDOrder)  = IntValue(1:LenDOrder)
              ENDIF
           ENDIF
        ENDDO
        ENDDO
        ENDDO 
!       Allocate the Link Node, Test if it is an EnddLeaf
        Node%NLinks=NLinks
        Node%LeafType='NodeNode'
        IF(Node%NLinks == 0) THEN
           Node%LeafType='EnddLeaf'
           RETURN
        ENDIF
!       Allocate the Links
        ALLOCATE(Node%Links(Node%NLinks),STAT=Status)
        IF(Status/=SUCCEED)CALL Halt('Link ALLOCATE failed in MakePIGWTree ')
!       Create the New SubWavelets Nodes and Store the Relavent Info
        DO ILink = 1, Node%NLinks
!          Allocate Node
           CALL  NewPIGWNode(Node%Links(ILink)%Link,NewLevel,ITypeArray(ILink),XposArray(ILink,1:3))
!          Compute Bounding Box
           XL = Node%Links(ILink)%Link%X-Node%Links(ILink)%Link%DX
           XH = Node%Links(ILink)%Link%X+Node%Links(ILink)%Link%DX
           CALL NodeBoundingBox(XL,XH,Node%Links(ILink)%Link%Box)
!          Store the Coefs
           Node%Links(ILink)%Link%WCoef     = WCoefArray(ILink,1:LenDOrder)
           Node%Links(ILink)%Link%IntFactor = IntVArray(ILink,1:LenDOrder)
           IF(ILeafType(ILink)==1) THEN
               Node%Links(ILink)%Link%LeafType='LeafLeaf'
           ELSE
               Node%Links(ILink)%Link%LeafType='EnddLeaf'
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
   SUBROUTINE NormalizeWCoefs(WCoef,DX,DY,DZ)
     INTEGER                           :: I,J,K,IJK
     REAL(DOUBLE),DIMENSION(LenDOrder) :: WCoef
     REAL(DOUBLE)                      :: DX,DY,DZ  
!
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
   SUBROUTINE ComputeWCoef(XX,NewLevel,NewIType,WECoef,WVCoef)
       REAL(DOUBLE),DIMENSION(3)         :: XX,XF
       REAL(DOUBLE),DIMENSION(3)         :: DX,DX2,DX3,XD,XD2,XD3
       INTEGER                           :: NewLevel,NewIType,NN,I,J,K
       INTEGER                           :: IJK0,IJK1,IJK2,IJK3
       REAL(DOUBLE),DIMENSION(LenDOrder) :: RhoX,WECoef,WVCoef,Tmp1,Tmp2
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
       CALL ExcVxc(RhoX,Tmp1,Tmp2)
       WECoef = Tmp1
       WVCoef = Tmp2
       SELECT CASE(DOrder)
       CASE(0)
          SELECT CASE(NewIType)
          CASE(1100)
             DO NN=1,2
                XF(1) = XX(1)+DX(1)*(2*NN-3)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp1,Tmp2)
                WECoef(1) = WECoef(1) + WFilter(NN,1,1)*Tmp1(1)
                WVCoef(1) = WVCoef(1) + WFilter(NN,1,1)*Tmp2(1)
             ENDDO
          CASE(1010)
             DO NN=1,2
                XF(2) = XX(2)+DX(2)*(2*NN-3)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp1,Tmp2)
                WECoef(1) = WECoef(1) + WFilter(NN,1,1)*Tmp1(1)
                WVCoef(1) = WVCoef(1) + WFilter(NN,1,1)*Tmp2(1)
             ENDDO
          CASE(1001)
             DO NN=1,2
                XF(3) = XX(3)+DX(3)*(2*NN-3)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp1,Tmp2)
                WECoef(1) = WECoef(1) + WFilter(NN,1,1)*Tmp1(1)
                WVCoef(1) = WVCoef(1) + WFilter(NN,1,1)*Tmp2(1)
             ENDDO
          END SELECT
       CASE(1)
          SELECT CASE(NewIType)
          CASE(1100)
             DO NN=1,2
                XF(1) = XX(1)+DX(1)*(2*NN-3)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp1,Tmp2)
                DO K=0,DOrder
                   DO J=0,DOrder
                      IJK0 = 1+DMax*J+DMax*DMax*K
                      IJK1 = IJK0+1
                      WECoef(IJK0) = WECoef(IJK0) +       WFilter(NN,1,1)*Tmp1(IJK0)   &
                                                  + DX(1)*WFilter(NN,2,1)*Tmp1(IJK1)
                      WECoef(IJK1) = WECoef(IJK1) + XD(1)*WFilter(NN,1,2)*Tmp1(IJK0)   &
                                                  +       WFilter(NN,2,2)*Tmp1(IJK1)
!
                      WVCoef(IJK0) = WVCoef(IJK0) +       WFilter(NN,1,1)*Tmp2(IJK0)   &
                                                  + DX(1)*WFilter(NN,2,1)*Tmp2(IJK1)
                      WVCoef(IJK1) = WVCoef(IJK1) + XD(1)*WFilter(NN,1,2)*Tmp2(IJK0)   &
                                                  +       WFilter(NN,2,2)*Tmp2(IJK1)
                   ENDDO
                ENDDO
             ENDDO
             CALL NormalizeWCoefs(WECoef,DX(1),Two*DX(2),Two*DX(3))
             CALL NormalizeWCoefs(WVCoef,DX(1),Two*DX(2),Two*DX(3))
          CASE(1010)
             DO NN=1,2
                XF(2) = XX(2)+DX(2)*(2*NN-3)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp1,Tmp2)
                DO K=0,DOrder
                   DO I=0,DOrder
                      IJK0 = 1+I+DMax*DMax*K
                      IJK1 = IJK0+DMax
                      WECoef(IJK0) = WECoef(IJK0) +       WFilter(NN,1,1)*Tmp1(IJK0)    &
                                                  + DX(2)*WFilter(NN,2,1)*Tmp1(IJK1)
                      WECoef(IJK1) = WECoef(IJK1) + XD(2)*WFilter(NN,1,2)*Tmp1(IJK0)    &
                                                  +       WFilter(NN,2,2)*Tmp1(IJK1)
!
                      WVCoef(IJK0) = WVCoef(IJK0) +       WFilter(NN,1,1)*Tmp2(IJK0)    &
                                                  + DX(2)*WFilter(NN,2,1)*Tmp2(IJK1)
                      WVCoef(IJK1) = WVCoef(IJK1) + XD(2)*WFilter(NN,1,2)*Tmp2(IJK0)    &
                                                  +       WFilter(NN,2,2)*Tmp2(IJK1)
                   ENDDO
                ENDDO
             ENDDO
             CALL NormalizeWCoefs(WECoef,DX(1),DX(2),Two*DX(3))
             CALL NormalizeWCoefs(WVCoef,DX(1),DX(2),Two*DX(3))
          CASE(1001)
            DO NN=1,2
                XF(3) = XX(3)+DX(3)*(2*NN-3)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp1,Tmp2)
                DO J=0,DOrder
                   DO I=0,DOrder
                      IJK0 = 1+I+DMax*J
                      IJK1 = IJK0+DMax*DMax
                      WECoef(IJK0) = WECoef(IJK0) +       WFilter(NN,1,1)*Tmp1(IJK0)    &
                                                  + DX(3)*WFilter(NN,2,1)*Tmp1(IJK1)
                      WECoef(IJK1) = WECoef(IJK1) + XD(3)*WFilter(NN,1,2)*Tmp1(IJK0)    &
                                                  +       WFilter(NN,2,2)*Tmp1(IJK1)
!
                      WVCoef(IJK0) = WVCoef(IJK0) +       WFilter(NN,1,1)*Tmp2(IJK0)    &
                                                  + DX(3)*WFilter(NN,2,1)*Tmp2(IJK1)
                      WVCoef(IJK1) = WVCoef(IJK1) + XD(3)*WFilter(NN,1,2)*Tmp2(IJK0)    &
                                                  +       WFilter(NN,2,2)*Tmp2(IJK1)
                   ENDDO
                ENDDO
             ENDDO
             CALL NormalizeWCoefs(WECoef,DX(1),DX(2),DX(3))
             CALL NormalizeWCoefs(WVCoef,DX(1),DX(2),DX(3))
          END SELECT
       CASE(2)
          SELECT CASE(NewIType)
          CASE(1100)
             DO NN=1,2
                XF(1) = XX(1)+DX(1)*(2*NN-3)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp1,Tmp2)
                DO K=0,DOrder
                   DO J=0,DOrder
                      IJK0 = 1+DMax*J+DMax*DMax*K
                      IJK1 = IJK0+1
                      IJK2 = IJK0+2
                      WECoef(IJK0) = WECoef(IJK0)  +        WFilter(NN,1,1)*Tmp1(IJK0)   &
                                                   +  DX(1)*WFilter(NN,2,1)*Tmp1(IJK1)   &
                                                   + DX2(1)*WFilter(NN,3,1)*Tmp1(IJK2)  
                      WECoef(IJK1) = WECoef(IJK1)  +  XD(1)*WFilter(NN,1,2)*Tmp1(IJK0)   &
                                                   +        WFilter(NN,2,2)*Tmp1(IJK1)   &
                                                   +  DX(1)*WFilter(NN,3,2)*Tmp1(IJK2)
                      WECoef(IJK2) = WECoef(IJK2)  + XD2(1)*WFilter(NN,1,3)*Tmp1(IJK0)   &
                                                   +  XD(1)*WFilter(NN,2,3)*Tmp1(IJK1)   &
                                                   +        WFilter(NN,3,3)*Tmp1(IJK2)
!
                      WVCoef(IJK0) = WVCoef(IJK0)  +        WFilter(NN,1,1)*Tmp2(IJK0)   &
                                                   +  DX(1)*WFilter(NN,2,1)*Tmp2(IJK1)   &
                                                   + DX2(1)*WFilter(NN,3,1)*Tmp2(IJK2)  
                      WVCoef(IJK1) = WVCoef(IJK1)  +  XD(1)*WFilter(NN,1,2)*Tmp2(IJK0)   &
                                                   +        WFilter(NN,2,2)*Tmp2(IJK1)   &
                                                   +  DX(1)*WFilter(NN,3,2)*Tmp2(IJK2)
                      WVCoef(IJK2) = WVCoef(IJK2)  + XD2(1)*WFilter(NN,1,3)*Tmp2(IJK0)   &
                                                   +  XD(1)*WFilter(NN,2,3)*Tmp2(IJK1)   &
                                                   +        WFilter(NN,3,3)*Tmp2(IJK2)
                   ENDDO
                ENDDO
             ENDDO
             CALL NormalizeWCoefs(WECoef,DX(1),Two*DX(2),Two*DX(3))
             CALL NormalizeWCoefs(WVCoef,DX(1),Two*DX(2),Two*DX(3))
          CASE(1010)
             DO NN=1,2
                XF(2) = XX(2)+DX(2)*(2*NN-3)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp1,Tmp2)
                DO K=0,DOrder
                   DO I=0,DOrder
                      IJK0  = 1+I+DMax*DMax*K
                      IJK1  = IJK0+1*DMax
                      IJK2  = IJK0+2*DMax
                      WECoef(IJK0) = WECoef(IJK0) +        WFilter(NN,1,1)*Tmp1(IJK0)   &
                                                  +  DX(2)*WFilter(NN,2,1)*Tmp1(IJK1)   &
                                                  + DX2(2)*WFilter(NN,3,1)*Tmp1(IJK2)  
                      WECoef(IJK1) = WECoef(IJK1) +  XD(2)*WFilter(NN,1,2)*Tmp1(IJK0)   &
                                                  +        WFilter(NN,2,2)*Tmp1(IJK1)   &
                                                  +  DX(2)*WFilter(NN,3,2)*Tmp1(IJK2)
                      WECoef(IJK2) = WECoef(IJK2) + XD2(2)*WFilter(NN,1,3)*Tmp1(IJK0)   &
                                                  +  XD(2)*WFilter(NN,2,3)*Tmp1(IJK1)   &
                                                  +        WFilter(NN,3,3)*Tmp1(IJK2)
!
                      WVCoef(IJK0) = WVCoef(IJK0) +        WFilter(NN,1,1)*Tmp2(IJK0)   &
                                                  +  DX(2)*WFilter(NN,2,1)*Tmp2(IJK1)   &
                                                  + DX2(2)*WFilter(NN,3,1)*Tmp2(IJK2)  
                      WVCoef(IJK1) = WVCoef(IJK1) +  XD(2)*WFilter(NN,1,2)*Tmp2(IJK0)   &
                                                  +        WFilter(NN,2,2)*Tmp2(IJK1)   &
                                                  +  DX(2)*WFilter(NN,3,2)*Tmp2(IJK2)
                      WVCoef(IJK2) = WVCoef(IJK2) + XD2(2)*WFilter(NN,1,3)*Tmp2(IJK0)   &
                                                  +  XD(2)*WFilter(NN,2,3)*Tmp2(IJK1)   &
                                                  +        WFilter(NN,3,3)*Tmp2(IJK2)
                   ENDDO
                ENDDO
             ENDDO
             CALL NormalizeWCoefs(WECoef,DX(1),DX(2),Two*DX(3))
             CALL NormalizeWCoefs(WVCoef,DX(1),DX(2),Two*DX(3))
          CASE(1001)
            DO NN=1,2
                XF(3) = XX(3)+DX(3)*(2*NN-3)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp1,Tmp2)
                DO J=0,DOrder
                   DO I=0,DOrder
                      IJK0  = 1+I+DMax*J
                      IJK1  = IJK0+1*DMax*DMax
                      IJK2  = IJK0+2*DMax*DMax
                      WECoef(IJK0) = WECoef(IJK0) +        WFilter(NN,1,1)*Tmp1(IJK0)   &
                                                  +  DX(3)*WFilter(NN,2,1)*Tmp1(IJK1)   &
                                                  + DX2(3)*WFilter(NN,3,1)*Tmp1(IJK2)  
                      WECoef(IJK1) = WECoef(IJK1) +  XD(3)*WFilter(NN,1,2)*Tmp1(IJK0)   &
                                                  +        WFilter(NN,2,2)*Tmp1(IJK1)   &
                                                  +  DX(3)*WFilter(NN,3,2)*Tmp1(IJK2)
                      WECoef(IJK2) = WECoef(IJK2) + XD2(3)*WFilter(NN,1,3)*Tmp1(IJK0)   &
                                                  +  XD(3)*WFilter(NN,2,3)*Tmp1(IJK1)   &
                                                  +        WFilter(NN,3,3)*Tmp1(IJK2)
!
                      WVCoef(IJK0) = WVCoef(IJK0) +        WFilter(NN,1,1)*Tmp2(IJK0)   &
                                                  +  DX(3)*WFilter(NN,2,1)*Tmp2(IJK1)   &
                                                  + DX2(3)*WFilter(NN,3,1)*Tmp2(IJK2)  
                      WVCoef(IJK1) = WVCoef(IJK1) +  XD(3)*WFilter(NN,1,2)*Tmp2(IJK0)   &
                                                  +        WFilter(NN,2,2)*Tmp2(IJK1)   &
                                                  +  DX(3)*WFilter(NN,3,2)*Tmp2(IJK2)
                      WVCoef(IJK2) = WVCoef(IJK2) + XD2(3)*WFilter(NN,1,3)*Tmp2(IJK0)   &
                                                  +  XD(3)*WFilter(NN,2,3)*Tmp2(IJK1)   &
                                                  +        WFilter(NN,3,3)*Tmp2(IJK2)
                   ENDDO
                ENDDO
             ENDDO
             CALL NormalizeWCoefs(WECoef,DX(1),DX(2),DX(3))
             CALL NormalizeWCoefs(WVCoef,DX(1),DX(2),DX(3))
          END SELECT
       CASE(3)
          SELECT CASE(NewIType)
          CASE(1100)
             DO NN=1,2
                XF(1) = XX(1)+DX(1)*(2*NN-3)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp1,Tmp2)
                DO K=0,DOrder
                   DO J=0,DOrder
                      IJK0 = 1+DMax*J+DMax*DMax*K
                      IJK1 = IJK0+1
                      IJK2 = IJK0+2
                      IJK3 = IJK0+3
                      WECoef(IJK0) = WECoef(IJK0) +        WFilter(NN,1,1)*Tmp1(IJK0)  &
                                                  +  DX(1)*WFilter(NN,2,1)*Tmp1(IJK1)  &
                                                  + DX2(1)*WFilter(NN,3,1)*Tmp1(IJK2)  &
                                                  + DX3(1)*WFilter(NN,4,1)*Tmp1(IJK3) 
                      WECoef(IJK1) = WECoef(IJK1) +  XD(1)*WFilter(NN,1,2)*Tmp1(IJK0)  &
                                                  +        WFilter(NN,2,2)*Tmp1(IJK1)  &
                                                  +  DX(1)*WFilter(NN,3,2)*Tmp1(IJK2)  &
                                                  + DX2(1)*WFilter(NN,4,2)*Tmp1(IJK3)
                      WECoef(IJK2) = WECoef(IJK2) + XD2(1)*WFilter(NN,1,3)*Tmp1(IJK0)  &
                                                  +  XD(1)*WFilter(NN,2,3)*Tmp1(IJK1)  &
                                                  +        WFilter(NN,3,3)*Tmp1(IJK2)  &
                                                  +  DX(1)*WFilter(NN,4,3)*Tmp1(IJK3)
                      WECoef(IJK3) = WECoef(IJK3) + XD3(1)*WFilter(NN,1,4)*Tmp1(IJK0)  &
                                                  + XD2(1)*WFilter(NN,2,4)*Tmp1(IJK1)  &
                                                  +  XD(1)*WFilter(NN,3,4)*Tmp1(IJK2)  &
                                                  +        WFilter(NN,4,4)*Tmp1(IJK3)
!
                      WVCoef(IJK0) = WVCoef(IJK0) +        WFilter(NN,1,1)*Tmp2(IJK0)  &
                                                  +  DX(1)*WFilter(NN,2,1)*Tmp2(IJK1)  &
                                                  + DX2(1)*WFilter(NN,3,1)*Tmp2(IJK2)  &
                                                  + DX3(1)*WFilter(NN,4,1)*Tmp2(IJK3) 
                      WVCoef(IJK1) = WVCoef(IJK1) +  XD(1)*WFilter(NN,1,2)*Tmp2(IJK0)  &
                                                  +        WFilter(NN,2,2)*Tmp2(IJK1)  &
                                                  +  DX(1)*WFilter(NN,3,2)*Tmp2(IJK2)  &
                                                  + DX2(1)*WFilter(NN,4,2)*Tmp2(IJK3)
                      WVCoef(IJK2) = WVCoef(IJK2) + XD2(1)*WFilter(NN,1,3)*Tmp2(IJK0)  &
                                                  +  XD(1)*WFilter(NN,2,3)*Tmp2(IJK1)  &
                                                  +        WFilter(NN,3,3)*Tmp2(IJK2)  &
                                                  +  DX(1)*WFilter(NN,4,3)*Tmp2(IJK3)
                      WVCoef(IJK3) = WVCoef(IJK3) + XD3(1)*WFilter(NN,1,4)*Tmp2(IJK0)  &
                                                  + XD2(1)*WFilter(NN,2,4)*Tmp2(IJK1)  &
                                                  +  XD(1)*WFilter(NN,3,4)*Tmp2(IJK2)  &
                                                  +        WFilter(NN,4,4)*Tmp2(IJK3)
                   ENDDO
                ENDDO
             ENDDO
             CALL NormalizeWCoefs(WECoef,DX(1),Two*DX(2),Two*DX(3))
             CALL NormalizeWCoefs(WVCoef,DX(1),Two*DX(2),Two*DX(3))
          CASE(1010)
             DO NN=1,2
                XF(2) = XX(2)+DX(2)*(2*NN-3)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp1,Tmp2)
                DO K=0,DOrder
                   DO I=0,DOrder
                      IJK0 = 1+I+DMax*DMax*K
                      IJK1 = IJK0+1*DMax
                      IJK2 = IJK0+2*DMax
                      IJK3 = IJK0+3*DMax
                      WECoef(IJK0) = WECoef(IJK0) +        WFilter(NN,1,1)*Tmp1(IJK0)  &
                                                  +  DX(2)*WFilter(NN,2,1)*Tmp1(IJK1)  &
                                                  + DX2(2)*WFilter(NN,3,1)*Tmp1(IJK2)  &
                                                  + DX3(2)*WFilter(NN,4,1)*Tmp1(IJK3) 
                      WECoef(IJK1) = WECoef(IJK1) +  XD(2)*WFilter(NN,1,2)*Tmp1(IJK0)  &
                                                  +        WFilter(NN,2,2)*Tmp1(IJK1)  &
                                                  +  DX(2)*WFilter(NN,3,2)*Tmp1(IJK2)  &
                                                  + DX2(2)*WFilter(NN,4,2)*Tmp1(IJK3)
                      WECoef(IJK2) = WECoef(IJK2) + XD2(2)*WFilter(NN,1,3)*Tmp1(IJK0)  &
                                                  +  XD(2)*WFilter(NN,2,3)*Tmp1(IJK1)  &
                                                  +        WFilter(NN,3,3)*Tmp1(IJK2)  &
                                                  +  DX(2)*WFilter(NN,4,3)*Tmp1(IJK3)
                      WECoef(IJK3) = WECoef(IJK3) + XD3(2)*WFilter(NN,1,4)*Tmp1(IJK0)  &
                                                  + XD2(2)*WFilter(NN,2,4)*Tmp1(IJK1)  &
                                                  +  XD(2)*WFilter(NN,3,4)*Tmp1(IJK2)  &
                                                  +        WFilter(NN,4,4)*Tmp1(IJK3)
!
                      WVCoef(IJK0) = WVCoef(IJK0) +        WFilter(NN,1,1)*Tmp2(IJK0)  &
                                                  +  DX(2)*WFilter(NN,2,1)*Tmp2(IJK1)  &
                                                  + DX2(2)*WFilter(NN,3,1)*Tmp2(IJK2)  &
                                                  + DX3(2)*WFilter(NN,4,1)*Tmp2(IJK3) 
                      WVCoef(IJK1) = WVCoef(IJK1) +  XD(2)*WFilter(NN,1,2)*Tmp2(IJK0)  &
                                                  +        WFilter(NN,2,2)*Tmp2(IJK1)  &
                                                  +  DX(2)*WFilter(NN,3,2)*Tmp2(IJK2)  &
                                                  + DX2(2)*WFilter(NN,4,2)*Tmp2(IJK3)
                      WVCoef(IJK2) = WVCoef(IJK2) + XD2(2)*WFilter(NN,1,3)*Tmp2(IJK0)  &
                                                  +  XD(2)*WFilter(NN,2,3)*Tmp2(IJK1)  &
                                                  +        WFilter(NN,3,3)*Tmp2(IJK2)  &
                                                  +  DX(2)*WFilter(NN,4,3)*Tmp2(IJK3)
                      WVCoef(IJK3) = WVCoef(IJK3) + XD3(2)*WFilter(NN,1,4)*Tmp2(IJK0)  &
                                                  + XD2(2)*WFilter(NN,2,4)*Tmp2(IJK1)  &
                                                  +  XD(2)*WFilter(NN,3,4)*Tmp2(IJK2)  &
                                                  +        WFilter(NN,4,4)*Tmp2(IJK3)
                   ENDDO
                ENDDO
             ENDDO
             CALL NormalizeWCoefs(WECoef,DX(1),DX(2),Two*DX(3))
             CALL NormalizeWCoefs(WVCoef,DX(1),DX(2),Two*DX(3))
          CASE(1001)
             DO NN=1,2
                XF(3) = XX(3)+DX(3)*(2*NN-3)
                RhoX  = RhoAtX(XF)
                CALL ExcVxc(RhoX,Tmp1,Tmp2)
                DO J=0,DOrder
                   DO I=0,DOrder
                      IJK0 = 1+I+DMax*J
                      IJK1 = IJK0+1*DMax*DMax
                      IJK2 = IJK0+2*DMax*DMax
                      IJK3 = IJK0+3*DMax*DMax
                      WECoef(IJK0) = WECoef(IJK0) +        WFilter(NN,1,1)*Tmp1(IJK0)  &
                                                  +  DX(3)*WFilter(NN,2,1)*Tmp1(IJK1)  &
                                                  + DX2(3)*WFilter(NN,3,1)*Tmp1(IJK2)  &
                                                 + DX3(3)*WFilter(NN,4,1)*Tmp1(IJK3) 
                      WECoef(IJK1) = WECoef(IJK1) +  XD(3)*WFilter(NN,1,2)*Tmp1(IJK0)  &
                                                  +        WFilter(NN,2,2)*Tmp1(IJK1)  &
                                                  +  DX(3)*WFilter(NN,3,2)*Tmp1(IJK2)  &
                                                  + DX2(3)*WFilter(NN,4,2)*Tmp1(IJK3)
                      WECoef(IJK2) = WECoef(IJK2) + XD2(3)*WFilter(NN,1,3)*Tmp1(IJK0)  &
                                                  +  XD(3)*WFilter(NN,2,3)*Tmp1(IJK1)  &
                                                  +        WFilter(NN,3,3)*Tmp1(IJK2)  &
                                                  +  DX(3)*WFilter(NN,4,3)*Tmp1(IJK3)
                      WECoef(IJK3) = WECoef(IJK3) + XD3(3)*WFilter(NN,1,4)*Tmp1(IJK0)  &
                                                  + XD2(3)*WFilter(NN,2,4)*Tmp1(IJK1)  &
                                                  +  XD(3)*WFilter(NN,3,4)*Tmp1(IJK2)  &
                                                  +        WFilter(NN,4,4)*Tmp1(IJK3)
!
                      WVCoef(IJK0) = WVCoef(IJK0) +        WFilter(NN,1,1)*Tmp2(IJK0)  &
                                                +  DX(3)*WFilter(NN,2,1)*Tmp2(IJK1)  &
                                                + DX2(3)*WFilter(NN,3,1)*Tmp2(IJK2)  &
                                                + DX3(3)*WFilter(NN,4,1)*Tmp2(IJK3) 
                      WVCoef(IJK1) = WVCoef(IJK1) +  XD(3)*WFilter(NN,1,2)*Tmp2(IJK0)  &
                                                +        WFilter(NN,2,2)*Tmp2(IJK1)  &
                                                +  DX(3)*WFilter(NN,3,2)*Tmp2(IJK2)  &
                                                + DX2(3)*WFilter(NN,4,2)*Tmp2(IJK3)
                      WVCoef(IJK2) = WVCoef(IJK2) + XD2(3)*WFilter(NN,1,3)*Tmp2(IJK0)  &
                                                +  XD(3)*WFilter(NN,2,3)*Tmp2(IJK1)  &
                                                +        WFilter(NN,3,3)*Tmp2(IJK2)  &
                                                +  DX(3)*WFilter(NN,4,3)*Tmp2(IJK3)
                      WVCoef(IJK3) = WVCoef(IJK3) + XD3(3)*WFilter(NN,1,4)*Tmp2(IJK0)  &
                                                + XD2(3)*WFilter(NN,2,4)*Tmp2(IJK1)  &
                                                +  XD(3)*WFilter(NN,3,4)*Tmp2(IJK2)  &
                                                +        WFilter(NN,4,4)*Tmp2(IJK3)
                   ENDDO
                ENDDO
             ENDDO
             CALL NormalizeWCoefs(WECoef,DX(1),DX(2),DX(3))
             CALL NormalizeWCoefs(WVCoef,DX(1),DX(2),DX(3))
          END SELECT
       END SELECT
 
!               
   END SUBROUTINE ComputeWCoef
!==================================================================
!  Wrapper For HasThisNode
!==================================================================
   FUNCTION HasNode(XX)
     LOGICAL                    :: HasNode
     REAL(DOUBLE),DIMENSION(3)  :: XX
!
     IF(ABS(XX(1)-XCBox%Center(1))>XCBox%Half(1)) THEN
        HasNode = .TRUE.
        RETURN
     ELSEIF(ABS(XX(2)-XCBox%Center(2))>XCBox%Half(2)) THEN
        HasNode = .TRUE.
        RETURN
     ELSEIF(ABS(XX(3)-XCBox%Center(3))>XCBox%Half(3)) THEN
        HasNode = .TRUE.
        RETURN
     ELSE
        IsTrue = .FALSE.
        Xpos   = XX
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
!    Determine if in in Bounding Box 
     IF(ABS(Node%Box%Center(1)-Xpos(1)) > Node%Box%Half(1)) RETURN
     IF(ABS(Node%Box%Center(2)-Xpos(2)) > Node%Box%Half(2)) RETURN
     IF(ABS(Node%Box%Center(3)-Xpos(3)) > Node%Box%Half(3)) RETURN
!    Determine is the X positions are te Same
     IF(ABS(Node%X(1)-Xpos(1)) < DTol .AND. &
        ABS(Node%X(2)-Xpos(2)) < DTol .AND. &
        ABS(Node%X(3)-Xpos(3)) < DTol) THEN
        IsTrue = .TRUE.
        RETURN
     ENDIF
!    If No links, return
     IF(Node%NLinks==0) RETURN
!    Reccur Down
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
     INTEGER                           :: NC
!
     RhoAtX(:)  = Zero
     NCallToRho = NCallToRho+1
     DO NC=1,CS_IN%NCells
        Xpos(1:3) = XX(1:3)+CS_IN%CellCarts%D(1:3,NC)
        RhoSum(:)=Zero
        CALL RhoAtPointX(RhoRoot)
        RhoAtX(:) = RhoAtX(:)+RhoSum(:)
     ENDDO
!
   END FUNCTION RhoAtX
!=================================================================================
!  Sums the significant contributions (leaves) to the density at a point
!=================================================================================
   RECURSIVE SUBROUTINE RhoAtPointX(Node)
     TYPE(RhoNode), POINTER                     :: Node
     INTEGER                                    :: L,M,N,LMN,L1,L2,J
     REAL(DOUBLE)                               :: RQx,RQy,RQz,RQ2,RL1
     REAL(DOUBLE)                               :: Xpt,Zeta,TwoZ,X
     REAL(DOUBLE),DIMENSION(1:43)               :: W
     REAL(DOUBLE),DIMENSION(0:7)                :: LamX,LamY,LamZ
!----------------------------------------------------------------------------------
     RQx=Node%Box%Center(1)-Xpos(1)
     IF(ABS(RQx)>Node%Box%Half(1)) RETURN
     RQy=Node%Box%Center(2)-Xpos(2)
     IF(ABS(RQy)>Node%Box%Half(2)) RETURN
     RQz=Node%Box%Center(3)-Xpos(3)
     IF(ABS(RQz)>Node%Box%Half(3)) RETURN
!    
     IF(Node%Leaf)THEN   
        RQ2  = RQx*RQx+RQy*RQy+RQz*RQz
        Zeta = Node%Zeta
        TwoZ = Two*Zeta
        X    = Node%Zeta*RQ2
        IF(X<Exp_Switch)THEN
           J=AINT(X*Exp_Grid)
           LamX(0) = One
           LamY(0) = One
           LamZ(0) = Exp_0(J)+X*(Exp_1(J)+X*(Exp_2(J)+X*(Exp_3(J)+X*Exp_4(J))))
           LamX(1) = TwoZ*RQx
           LamY(1) = TwoZ*RQy
           LamZ(1) = TwoZ*RQz*LamZ(0)
           DO L=2,DOrder+Node%Ell
              L1=L-1
              L2=L-2
              RL1=DBLE(L1)
              LamX(L)=TwoZ*(RQx*LamX(L1)-RL1*LamX(L2))
              LamY(L)=TwoZ*(RQy*LamY(L1)-RL1*LamY(L2))
              LamZ(L)=TwoZ*(RQz*LamZ(L1)-RL1*LamZ(L2))
           ENDDO
           SELECT CASE(DOrder)
           CASE(0)
              SELECT CASE(Node%Ell)
              CASE(0)
                 INCLUDE "MMA/DivRho_0_0.Inc"
              CASE(1)
                 INCLUDE "MMA/DivRho_1_0.Inc"
              CASE(2)
                 INCLUDE "MMA/DivRho_2_0.Inc"
              CASE(3)
!                 INCLUDE "MMA/DivRho_3_0.Inc"
              CASE(4)
!                 INCLUDE "MMA/DivRho_4_0.Inc"
              END SELECT
           CASE(1)
              SELECT CASE(Node%Ell)
              CASE(0)
                 INCLUDE "MMA/DivRho_0_1.Inc"
              CASE(1)
                 INCLUDE "MMA/DivRho_1_1.Inc"
              CASE(2)
                 INCLUDE "MMA/DivRho_2_1.Inc"
              CASE(3)
!                 INCLUDE "MMA/DivRho_3_1.Inc"
              CASE(4)
!                 INCLUDE "MMA/DivRho_4_1.Inc"
              END SELECT
           CASE(2)
              SELECT CASE(Node%Ell)
              CASE(0)
                 INCLUDE "MMA/DivRho_0_2.Inc"
              CASE(1)
                 INCLUDE "MMA/DivRho_1_2.Inc"
              CASE(2)
                 INCLUDE "MMA/DivRho_2_2.Inc"
              CASE(3)
!                 INCLUDE "MMA/DivRho_3_2.Inc"
              CASE(4)
!                 INCLUDE "MMA/DivRho_4_2.Inc"
              END SELECT
           CASE(3)
              SELECT CASE(Node%Ell)
              CASE(0)
                 INCLUDE "MMA/DivRho_0_3.Inc"
              CASE(1)
                 INCLUDE "MMA/DivRho_1_3.Inc"
              CASE(2)
                 INCLUDE "MMA/DivRho_2_3.Inc"
              CASE(3)
!                 INCLUDE "MMA/DivRho_3_3.Inc"
              CASE(4)
!                 INCLUDE "MMA/DivRho_4_3.Inc"
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
!     Exc = RhoXX
!     Vxc = RhoXX
!     RETURN
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
              d0Vxc  =  Half*d1Exc
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
              d0Vxc  =  Half*d1Exc
              d1Vxc  =  Half*d2Exc
              d2Vxc  =  Half*d3Exc
              d3Vxc  =  Half*d4Exc
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
              d0Vxc  =  Half*d1Exc
              d1Vxc  =  Half*d2Exc
              d2Vxc  =  Half*d3Exc
              d3Vxc  =  Half*d4Exc
              d4Vxc  =  Half*d5Exc
              d5Vxc  =  Half*d6Exc
              d6Vxc  =  Half*d7Exc
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
              d0Vxc  =  Half*d1Exc
              d1Vxc  =  Half*d2Exc
              d2Vxc  =  Half*d3Exc
              d3Vxc  =  Half*d4Exc
              d4Vxc  =  Half*d5Exc
              d5Vxc  =  Half*d6Exc
              d6Vxc  =  Half*d7Exc
              d7Vxc  =  Half*d8Exc
              d8Vxc  =  Half*d9Exc
              d9Vxc  =  Half*d10Exc
           ENDIF
           INCLUDE "MMA/RhoToExc_3.Inc"
        END SELECT
     END SELECT
!
   END SUBROUTINE ExcVxc
!==========================================================================
!  Make the Bounding Box Periodic
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
!  Delete the Tree
!==========================================================================
   RECURSIVE SUBROUTINE DeleteDIPMWTree(Node)
     TYPE(PIGWNode), POINTER     :: Node
     INTEGER                     :: ILink,Status
!     
     DEALLOCATE(Node%WCoef,STAT=Status)
     IF(Status/=SUCCEED) CALL Halt('Node%WCoef DEALLOCATE failed in DeleteDIPMWTree ')          
     DEALLOCATE(Node%IntFactor,STAT=Status)
     IF(Status/=SUCCEED) CALL Halt('Node%IntFactor DEALLOCATE failed in DeleteDIPMWTree ') 
     IF(Node%LeafType=='EnddLeaf' .OR. Node%LeafType=='LeafLeaf' ) THEN
        NULLIFY(Node)
        RETURN
     ELSE
        DO ILink=1,Node%NLinks
           CALL DeleteDIPMWTree(Node%Links(ILink)%Link)
           NULLIFY(Node%Links(ILink)%Link)
        ENDDO
        DEALLOCATE(Node%Links,STAT=Status)
        IF(Status/=SUCCEED) CALL Halt('Node%Links DEALLOCATE failed in DeleteDIPMWTree ') 
     ENDIF

   END SUBROUTINE DeleteDIPMWTree
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
     IF(Node%LeafType=='EnddLeaf' .OR. Node%LeafType=='LeafLeaf' ) RETURN
     DO ILink=1,Node%NLinks 
        CALL ComputeExc(Node%Links(ILink)%Link)
     ENDDO
!
   END SUBROUTINE ComputeExc
!=================================================================================
!  Wrapper for VxcAtPoint
!=================================================================================
   FUNCTION VxcAtX(XX) 
     REAL(DOUBLE)                      :: VxcAtX
     REAL(DOUBLE),DIMENSION(3)         :: XX
     INTEGER                           :: NC
!
     VxcSum = Zero  
     Xpos   = XX
     CALL VxcAtPointX(PIGWRoot)
     VxcAtX = VxcSum
!
   END FUNCTION VxcAtX
!=================================================================================
!  Sums the significant contributions (leaves) to the density at a point
!=================================================================================
   RECURSIVE SUBROUTINE VxcAtPointX(Node)
     TYPE(PIGWNode), POINTER                    :: Node
     INTEGER                                    :: I,J,K,IJK,ILink
     REAL(DOUBLE)                               :: RQx,RQy,RQz,XR,YR,ZR,FacX,FacY,FacZ
     REAL(DOUBLE),DIMENSION(0:DOrder)           :: WVx,WVy,WVz
!----------------------------------------------------------------------------------
!    See if point overlaps with bounding box 
!    Test to See If We Overlap Wavelet in the Box
     RQx=Xpos(1)-Node%Box%Center(1)
     IF(ABS(RQx)>Node%Box%Half(1)) RETURN
     RQy=Xpos(2)-Node%Box%Center(2)
     IF(ABS(RQy)>Node%Box%Half(2)) RETURN
     RQz=Xpos(3)-Node%Box%Center(3)
     IF(ABS(RQz)>Node%Box%Half(3)) RETURN
!    Compute the Cintribution From the X Z and Z wavelets
     XR = -(Node%X(1)-Xpos(1) )/Node%DX(1)
     YR = -(Node%X(2)-Xpos(2) )/Node%DX(2)
     ZR = -(Node%X(3)-Xpos(3) )/Node%DX(3)
!
     WVx = 0.D0
     WVy = 0.D0
     WVz = 0.D0
     SELECT CASE (DOrder)
     CASE(0)
        WVx(0) = 1.D0-ABS(XR)
        WVy(0) = 1.D0-ABS(YR)
        WVz(0) = 1.D0-ABS(ZR)
     CASE(1)
        WVx(0) = 1.D0-3.D0*XR*XR+2.D0*ABS(XR)**3
        WVy(0) = 1.D0-3.D0*YR*YR+2.D0*ABS(YR)**3
        WVz(0) = 1.D0-3.D0*ZR*ZR+2.D0*ABS(ZR)**3
        !
        WVx(1) = XR*(1.D0-2.D0*ABS(XR)+XR*XR)
        WVy(1) = YR*(1.D0-2.D0*ABS(YR)+YR*YR)
        WVz(1) = ZR*(1.D0-2.D0*ABS(ZR)+ZR*ZR)
     CASE(2) 
        WVx(0) = 1.D0-10.D0*ABS(XR)**3+15.D0*XR**4-6.D0*ABS(XR)**5
        WVy(0) = 1.D0-10.D0*ABS(YR)**3+15.D0*YR**4-6.D0*ABS(YR)**5
        WVz(0) = 1.D0-10.D0*ABS(ZR)**3+15.D0*ZR**4-6.D0*ABS(ZR)**5
        !
        WVx(1) = XR*(1.D0-6.D0*XR**2+8.D0*ABS(XR)**3-3.D0*XR**4)
        WVy(1) = YR*(1.D0-6.D0*YR**2+8.D0*ABS(YR)**3-3.D0*YR**4)
        WVz(1) = ZR*(1.D0-6.D0*ZR**2+8.D0*ABS(ZR)**3-3.D0*ZR**4)
        !
        WVx(2) = XR*XR*(0.5D0-1.5D0*ABS(XR)+1.5D0*XR**2-0.5D0*ABS(XR)**3)
        WVy(2) = YR*YR*(0.5D0-1.5D0*ABS(YR)+1.5D0*YR**2-0.5D0*ABS(YR)**3)
        WVz(2) = ZR*ZR*(0.5D0-1.5D0*ABS(ZR)+1.5D0*ZR**2-0.5D0*ABS(ZR)**3)
     CASE(3)
        WVx(0) = 1.D0-35.D0*XR**4+84.D0*ABS(XR)**5-70.D0*XR**6+20.D0*ABS(XR)**7
        WVy(0) = 1.D0-35.D0*YR**4+84.D0*ABS(YR)**5-70.D0*YR**6+20.D0*ABS(YR)**7
        WVz(0) = 1.D0-35.D0*ZR**4+84.D0*ABS(ZR)**5-70.D0*ZR**6+20.D0*ABS(ZR)**7
        !
        WVx(1) = XR*(1.D0-20.D0*ABS(XR)**3+45.D0*XR**4-36.D0*ABS(XR)**5+10.D0*XR**6)
        WVy(1) = YR*(1.D0-20.D0*ABS(YR)**3+45.D0*YR**4-36.D0*ABS(YR)**5+10.D0*YR**6)
        WVz(1) = ZR*(1.D0-20.D0*ABS(ZR)**3+45.D0*ZR**4-36.D0*ABS(ZR)**5+10.D0*ZR**6)
        !
        WVx(2) = XR*XR*(0.5D0-5.D0*XR**2+10.D0*ABS(XR)**3-7.5D0*XR**4+2.D0*ABS(XR)**5)
        WVy(2) = YR*YR*(0.5D0-5.D0*YR**2+10.D0*ABS(YR)**3-7.5D0*YR**4+2.D0*ABS(YR)**5)
         WVz(2) = ZR*ZR*(0.5D0-5.D0*ZR**2+10.D0*ABS(ZR)**3-7.5D0*ZR**4+2.D0*ABS(ZR)**5)
        !
        WVx(3) = XR*XR*XR*((1.D0/6.D0)-(2.D0/3.D0)*ABS(XR)+XR**2-(2.D0/3.D0)*ABS(XR)**3+(1.D0/6.D0)*ABS(XR)**4)
        WVy(3) = YR*YR*YR*((1.D0/6.D0)-(2.D0/3.D0)*ABS(YR)+YR**2-(2.D0/3.D0)*ABS(YR)**3+(1.D0/6.D0)*ABS(YR)**4)
        WVz(3) = ZR*ZR*ZR*((1.D0/6.D0)-(2.D0/3.D0)*ABS(ZR)+ZR**2-(2.D0/3.D0)*ABS(ZR)**3+(1.D0/6.D0)*ABS(ZR)**4)
     END SELECT 
!    Sum it Up
     DO K=0,DOrder
        DO J=0,DOrder
           DO I=0,DOrder  
              IJK = 1+I+DMax*J+DMax*DMax*K
              VxcSum = VxcSum+Node%WCoef(IJK)*WVx(I)*WVy(J)*WVz(K)
           ENDDO
        ENDDO
     ENDDO
!    Return if endleaf
     IF(Node%LeafType=='EnddLeaf' .OR. Node%LeafType=='LeafLeaf' ) RETURN
!    Recure Down tree
     DO ILink=1,Node%NLinks 
        CALL VxcAtPointX(Node%Links(ILink)%Link)
     ENDDO

   END SUBROUTINE VxcAtPointX
!==================================================================
!  Print a Node
!==================================================================
   SUBROUTINE PrintNode(Node)
     TYPE(PIGWNode), POINTER            :: Node
     INTEGER                            :: I,J,K,IJK
     REAL(DOUBLE),DIMENSION(LenDOrder)  :: IntValue
     REAL(DOUBLE)                       :: MaxWeight
!
     WRITE(*,*) "----------------------------"
     WRITE(*,*) 'Level    = ',Node%Level
     WRITE(*,*) 'IType    = ',Node%IType
     WRITE(*,*) 'LeafType = ',Node%LeafType
     WRITE(*,*) 'NLinks   = ',Node%NLinks
     WRITE(*,'(A5,1X,F20.10,1X,F20.10,1X,F20.10)') 'X   :',Node%X(1:3)
     WRITE(*,'(A5,1X,F20.10,1X,F20.10,1X,F20.10)') 'DX  :',Node%DX(1:3)
     WRITE(*,*) 'Box'
     WRITE(*,'(A5,1X,F20.10,1X,F20.10,1X,F20.10)') 'Cent:',Node%Box%Center(1:3)
     WRITE(*,'(A5,1X,F20.10,1X,F20.10,1X,F20.10)') 'Half:',Node%Box%Half(1:3)
     WRITE(*,'(A5,1X,F20.10,1X,F20.10,1X,F20.10)') 'Xlow:',Node%Box%BndBox(1:3,1)
     WRITE(*,'(A5,1X,F20.10,1X,F20.10,1X,F20.10)') 'Xhig:',Node%Box%BndBox(1:3,2)
     WRITE(*,*) 'WCoef'
     DO K=0,DOrder
        DO J=0,DOrder
           IJK = 1+DMax*J+DMax*DMax*K
           WRITE(*,'4(5X,I1,I1,I1,1X,D16.8)') (I,J,K,Node%WCoef(IJK+I),I=0,DOrder) 
        ENDDO
     ENDDO
     WRITE(*,*) 'Integration  Factor'
     DO K=0,DOrder
        DO J=0,DOrder
           IJK = 1+DMax*J+DMax*DMax*K
           WRITE(*,'4(5X,I1,I1,I1,1X,D16.8)') (I,J,K,Node%IntFactor(IJK+I),I=0,DOrder) 
        ENDDO
     ENDDO
     CALL MaxWWeight(Node%X,Node%DX,Node%Box,IntValue)
     MaxWeight = Zero
     DO I=1,LenDOrder
        MaxWeight = MaxWeight+ABS(IntValue(I))*ABS(Node%WCoef(I))
     ENDDO
     WRITE(*,*) 'Max Weight = ',MaxWeight

     WRITE(*,*) "----------------------------"
!
   END SUBROUTINE PrintNode
!
!
!
END MODULE DIPMWTree
!!$        CASE(1)
!!$           WvltInt(1,1) = 0.0D0
!!$           WvltInt(2,1) = 13.0D0/24.0D0
!!$           WvltInt(3,1) = 5.0D-1
!!$           WvltInt(4,1) = 5.0D-1
!!$!
!!$           WFilter(1,1,1)= 6.25D-2
!!$           WFilter(2,1,1)=-5.625D-1
!!$           WFilter(3,1,1)=-5.625D-1
!!$           WFilter(4,1,1)= 6.25D-2
!!$        CASE(2)
!!$           WvltInt(1,1) = 0.0D0
!!$           WvltInt(2,1) = 0.0D0
!!$           WvltInt(3,1) = 0.0D0
!!$           WvltInt(4,1) = 0.0D0
!!$           WvltInt(5,1) = 0.0D0
!!$           WvltInt(6,1) = 5.0D-1
!!$!
!!$           WFilter(1,1,1)=-1.171875D-2
!!$           WFilter(2,1,1)= 9.765625D-2
!!$           WFilter(3,1,1)=-5.859375D-1
!!$           WFilter(4,1,1)=-5.859375D-1
!!$           WFilter(5,1,1)= 9.765625D-2
!!$           WFilter(6,1,1)=-1.171875D-2
!!$        CASE(3)
!!$           WvltInt(1,1) = 0.0D0
!!$           WvltInt(2,1) = 0.0D0
!!$           WvltInt(3,1) = 0.0D0
!!$           WvltInt(4,1) = 0.0D0
!!$           WvltInt(5,1) = 0.0D0
!!$           WvltInt(6,1) = 0.0D0
!!$           WvltInt(7,1) = 0.0D0
!!$           WvltInt(8,1) = 5.0D-1
!!$!
!!$           WFilter(1,1,1)= 2.4414062500D-3
!!$           WFilter(2,1,1)=-2.3925781250D-2
!!$           WFilter(3,1,1)= 1.1962890625D-1
!!$           WFilter(4,1,1)=-5.9814453125D-1
!!$           WFilter(5,1,1)=-5.9814453125D-1
!!$           WFilter(6,1,1)= 1.1962890625D-1
!!$           WFilter(7,1,1)=-2.3925781250D-2
!!$           WFilter(8,1,1)= 2.4414062500D-3
!!$!=================================================================================
!!$!  Test Function Integration
!!$!=================================================================================
!!$   SUBROUTINE TestRho()
!!$     REAL(DOUBLE)                    :: R2,Expt,Beta
!!$     INTEGER                         :: I,J,K,IJK
!!$     REAL(DOUBLE),DIMENSION(0:3)     :: DivX,DivY,DivZ
!!$!
!!$     Beta  = 2.D0
!!$!
!!$     R2    = (Xpos(1)-1.D0)**2+Xpos(2)**2+Xpos(3)**2
!!$     Expt  = EXP(-Beta*R2)
!!$!
!!$     DivX(0) = One
!!$     DivY(0) = One
!!$     DivZ(0) = One
!!$!     
!!$     DivX(1) = -Two*Beta*(Xpos(1)-1.D0)
!!$     DivY(1) = -Two*Beta*Xpos(2)
!!$     DivZ(1) = -Two*Beta*Xpos(3)
!!$!     
!!$     DivX(2) = -Two*Beta+(2*Beta*(Xpos(1)-1.D0))**2
!!$     DivY(2) = -Two*Beta+(2*Beta*Xpos(2))**2
!!$     DivZ(2) = -Two*Beta+(2*Beta*Xpos(3))**2
!!$!     
!!$     DivX(3) = 12.D0*Beta*Beta*(Xpos(1)-1.D0)-(Two*Beta*(Xpos(1)-1.D0))**3
!!$     DivY(3) = 12.D0*Beta*Beta*Xpos(2)-(Two*Beta*Xpos(2))**3
!!$     DivZ(3) = 12.D0*Beta*Beta*Xpos(3)-(Two*Beta*Xpos(3))**3
!!$!
!!$     DO K=0,DOrder
!!$        DO J=0,DOrder
!!$           DO I=0,DOrder
!!$              IJK = 1+I+DMax*J+DMax*DMax*K
!!$              RhoSum(IJK) = RhoSum(IJK)+DivX(I)*DivY(J)*DivZ(K)*Expt
!!$           ENDDO
!!$        ENDDO
!!$     ENDDO
!!$!
!!$     R2    = (Xpos(1)+1.D0)**2+Xpos(2)**2+Xpos(3)**2
!!$     Expt  = EXP(-Beta*R2)
!!$!
!!$     DivX(0) = One
!!$     DivY(0) = One
!!$     DivZ(0) = One
!!$!     
!!$     DivX(1) = -Two*Beta*(Xpos(1)+1.D0)
!!$     DivY(1) = -Two*Beta*Xpos(2)
!!$     DivZ(1) = -Two*Beta*Xpos(3)
!!$!     
!!$     DivX(2) = -Two*Beta+(2*Beta*(Xpos(1)+1.D0))**2
!!$     DivY(2) = -Two*Beta+(2*Beta*Xpos(2))**2
!!$     DivZ(2) = -Two*Beta+(2*Beta*Xpos(3))**2
!!$!     
!!$     DivX(3) = 12.D0*Beta*Beta*(Xpos(1)+1.D0)-(Two*Beta*(Xpos(1)+1.D0))**3
!!$     DivY(3) = 12.D0*Beta*Beta*Xpos(2)-(Two*Beta*Xpos(2))**3
!!$     DivZ(3) = 12.D0*Beta*Beta*Xpos(3)-(Two*Beta*Xpos(3))**3
!!$!
!!$     DO K=0,DOrder
!!$        DO J=0,DOrder
!!$           DO I=0,DOrder
!!$              IJK = 1+I+DMax*J+DMax*DMax*K
!!$              RhoSum(IJK) = RhoSum(IJK)+DivX(I)*DivY(J)*DivZ(K)*Expt
!!$           ENDDO
!!$        ENDDO
!!$     ENDDO
!!$!
!!$   END SUBROUTINE TestRho
