MODULE RhoTree
   USE Derivedtypes
   USE GlobalScalars   
   USE GlobalObjects
   USE ProcessControl
   USE Indexing
   USE InOut
#ifdef USE_SPECFUN
   USE SpecFun
#endif
   USE Macros
   USE Thresholding
   USE BoundingBox
   USE CubeGrid
   IMPLICIT NONE
!=================================================================================================
!  Hierarchical density node
!=================================================================================================
   TYPE RhoNode
      Logical                           :: Leaf     ! Is this a data containing node?
!     Indexes USEd for tree building
      INTEGER                           :: Bdex     ! Begining index of ORB list for this node
      INTEGER                           :: Edex     ! ENDign index of ORB list for this node
      INTEGER                           :: NQ       ! Number of centers
!     Bounding box 
      REAL(DOUBLE)                      :: ZetaMin  ! Minimum exponent in this node
      REAL(DOUBLE)                      :: MaxAmp   ! Max amplitude in this node
      REAL(DOUBLE)                      :: Extent   ! Penetration extent
      TYPE(BBox)                        :: Box          
!     Links
      TYPE(RhoNode),            POINTER :: Travrse  ! Next node in tree traversal
      TYPE(RhoNode),            POINTER :: Descend  ! Next node in tree descent
!     Density 
      INTEGER          :: Ell                       ! Ell type
      REAL(DOUBLE)     :: Zeta
      REAL(DOUBLE)     :: Qx
      REAL(DOUBLE)     :: Qy
      REAL(DOUBLE)     :: Qz
      REAL(DOUBLE),DIMENSION(:),POINTER :: Co       ! Coefficients of the HGTF density
   END TYPE                                      
!----------------------------------------------------------------------------------
!  Global parameters
!
   INCLUDE 'Exp.Inc'
   INCLUDE 'Erf.Inc'
   INTEGER, Parameter :: DistPerBox=1
!  Global scalars
   INTEGER            :: RhoNodes,RhoLevel
   REAL(DOUBLE)       :: MaxAmp,MiniumExp
!----------------------------------------------------------------------------------
!  Global density in array form
   TYPE(HGRho)                           :: Rho
   INTEGER,     DIMENSION(:),Allocatable :: Qdex      
   INTEGER,     DIMENSION(:),Allocatable :: Cdex
   INTEGER,     DIMENSION(:),Allocatable :: Ldex
   REAL(DOUBLE),DIMENSION(:),Allocatable :: RList      
   REAL(DOUBLE),DIMENSION(:),Allocatable :: Zeta
   REAL(DOUBLE),DIMENSION(:),Allocatable :: Amp
!----------------------------------------------------------------------------------
!  Global Rho trees
   TYPE(RhoNode), POINTER :: RhoRoot  ! Root of the tree 
   TYPE(RhoNode), POINTER :: RhoSubTree
!----------------------------------------------------------------------------------
!  Global griding variables
!
   TYPE(BBox)                       :: Box
   REAL(DOUBLE), DIMENSION(3,NGrid) :: Grid
   REAL(DOUBLE), DIMENSION(4,NGrid) :: RhoV 
   REAL(DOUBLE)                     :: Pop
   INTEGER                          :: IGrid
!-----------!
   CONTAINS !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!=================================================================================
!     
!=================================================================================
      SUBROUTINE RhoToTree(Args)
         TYPE(ARGMT)               :: Args
         INTEGER                   :: Status,K,I
!-------------------------------------------------------------------
!        Initialize the density
         CALL InitRho(Args)
!        Possibly print the density out
!         CALL Print_HGRho
!        Initialize counters
         RhoNodes=0
         RhoLevel=0
!        Initialize the root node
         CALL NewRhoNode(RhoRoot,0)
         CALL InitRhoRoot
!        Convert the density into a 4-D BinTree
         CALL SplitRho(RhoRoot)
!        Delete the density in the array structure 
         CALL DeleteRho
      END SUBROUTINE RhoToTree
!===============================================================================
!     
!=======================================================================
      SUBROUTINE InitRhoRoot
         TYPE(DBL_RNK2) :: BndBox
         RhoRoot%Bdex=1
         RhoRoot%Edex=Rho%NDist
         RhoRoot%NQ=Rho%NDist
!        Get the nuclear bounding box
         CALL New(BndBox,(/3,2/))
         CALL Get(BndBox,'boundingbox',CurGeom)
         RhoRoot%Box%BndBox(1:3,1:2)=BndBox%D(1:3,1:2)
         CALL Delete(BndBox)
!        Set the largest range (smallest exponent) of the roots BB
         RhoRoot%ZetaMin=Rho%Expt%D(1)
         MiniumExp=Rho%Expt%D(1)
!        Set the largest amplituded (max norm of the rho coefficients)
         RhoRoot%MaxAmp=MaxAmp
!        Set the center and width of the Cartesian part
         RhoRoot%Box%Half=Half*(RhoRoot%Box%BndBox(1:3,2)-RhoRoot%Box%BndBox(1:3,1))
         RhoRoot%Box%Center=Half*(RhoRoot%Box%BndBox(1:3,2)+RhoRoot%Box%BndBox(1:3,1))
     END SUBROUTINE InitRhoRoot
!===================================================================
!
!===================================================================
      RECURSIVE SUBROUTINE SplitRho(Node)
         TYPE(RhoNode), POINTER :: Node,Left,Right
!--------------------------------------------------------------
!         CALL PrintBBox(Node%Box,Node%ZetaMin,Node%MaxAmp)
         IF(Node%NQ<=DistPerBox)THEN
            CALL FillRhoLeaf(Node)
!            CALL CheckBounds(Node)
         ELSE 
!           Allocate new children 
            CALL NewRhoNode(Node%Descend,Node%Box%Tier+1)
            CALL NewRhoNode(Node%Descend%Travrse,Node%Box%Tier+1)
!           Set links
            Node%Descend%Travrse%Travrse=>Node%Travrse
            Left=>Node%Descend
            Right=>Node%Descend%Travrse
            CALL SplitRhoBox(Node,Left,Right)
!           Recur
            CALL SplitRho(Left)
            CALL SplitRho(Right)
         ENDIF
       END SUBROUTINE SplitRho
!================================================================================
!     Bisection   
!================================================================================
      SUBROUTINE SplitRhoBox(Node,Left,Right)
         TYPE(RhoNode), POINTER :: Node,Left,Right
         REAL(DOUBLE)           :: Section
         INTEGER                :: B,E,N,ISplit,Split,I,J,k
!--------------------------------------------------------------
!        Indexing
         J=0
         B=Node%Bdex
         E=Node%Edex
         N=E-B+1
!        Orthogonal direction
         ISplit=Mod(Node%Box%Tier,5)+1
         IF(ISplit==1)THEN
            DO I=B,E
               J=J+1;K=Qdex(I)
               RList(J)=Zeta(K)
            ENDDO
         ELSEIF(ISplit==2)THEN
            DO I=B,E
               J=J+1;K=Qdex(I)
               RList(J)=Amp(K)
            ENDDO
         ELSEIF(ISplit==3)THEN
            DO I=B,E
               J=J+1;K=Qdex(I)
               RList(J)=Rho%Qx%D(K)
            ENDDO
         ELSEIF(ISplit==4)THEN
            DO I=B,E
               J=J+1;K=Qdex(I)
               RList(J)=Rho%Qy%D(K)
            ENDDO
         ELSEIF(ISplit==5)THEN
            DO I=B,E
               J=J+1;K=Qdex(I)
               RList(J)=Rho%Qz%D(K)
            ENDDO
         ENDIF
!        Sort
         CALL DblIntSort77(N,RList,Qdex(B:E),2)
!        Orthogonal RECURSIVE bisection (ORB)
         Section   =RList(1)+Half*(RList(N)-RList(1))
!         Split     =BinarySearch(N,RList,Section)
         Split     =N/2
         Left%Bdex =B
         Right%Edex=E
         IF(Split==1)THEN 
            Split=2
            Left%Edex =Min(B,B+Split-2)
            Right%Bdex=Min(B+Split-1,E)
         ELSE
            Left%Edex =B+Split-2 
            Right%Bdex=B+Split-1 
         ENDIF
!        Counters
         Left%NQ=Left%Edex-Left%Bdex+1
         Right%NQ=Right%Edex-Right%Bdex+1
!        Find optimal boxes
         CALL NewRhoBox(Left)
         CALL NewRhoBox(Right)
      END SUBROUTINE SplitRhoBox
!================================================================================
!     New RhoBox
!================================================================================
      SUBROUTINE NewRhoBox(Node)
         TYPE(RhoNode), POINTER :: Node
         INTEGER                :: J,I,IQ
         REAL(DOUBLE)           :: Test
         Node%ZetaMin=1.D12
         Node%MaxAmp=-1.D12
         DO I=1,3; Node%Box%BndBox(I,1)=1.D12; ENDDO
         DO I=1,3; Node%Box%BndBox(I,2)=-1.D12; ENDDO
         DO J=Node%Bdex,Node%Edex 
            IQ=Qdex(J)
            Node%ZetaMin=Min(Node%ZetaMin,Zeta(IQ))
            Node%MaxAmp =Max(Node%MaxAmp ,Amp(IQ))
            Node%Box%BndBox(1,1)=Min(Node%Box%BndBox(1,1),Rho%Qx%D(IQ))
            Node%Box%BndBox(1,2)=Max(Node%Box%BndBox(1,2),Rho%Qx%D(IQ))
            Node%Box%BndBox(2,1)=Min(Node%Box%BndBox(2,1),Rho%Qy%D(IQ))
            Node%Box%BndBox(2,2)=Max(Node%Box%BndBox(2,2),Rho%Qy%D(IQ))
            Node%Box%BndBox(3,1)=Min(Node%Box%BndBox(3,1),Rho%Qz%D(IQ))
            Node%Box%BndBox(3,2)=Max(Node%Box%BndBox(3,2),Rho%Qz%D(IQ))
         ENDDO
!        New sides
         Node%Box%Half(1:3)=Half*(Node%Box%BndBox(:,2)-Node%Box%BndBox(:,1))
!        New center
         Node%Box%Center(1:3)=Half*(Node%Box%BndBox(:,2)+Node%Box%BndBox(:,1))
      END SUBROUTINE NewRhoBox
!================================================================================
!     Binary search
!================================================================================
      FUNCTION BinarySearch(M,List,Key) Result(I)
         INTEGER :: I,K,M,High,Low      
         REAL(DOUBLE)                 :: Key
         REAL(DOUBLE), DIMENSION(1:M) :: List 
         IF(List(1)>Key)THEN
            CALL Halt(' Logical error 1 in BinarySearch ')
         ELSEIF(List(M)<Key)THEN
            CALL Halt(' Logical error 2 in BinarySearch ')
         ELSE
            High=M
            Low=0
            DO K=1,M
              I=(High+Low)/2
              IF(Key<List(I))THEN
                 High=I
              ELSE
                 Low=I
              ENDIF
              IF(High-Low<=1)Return
           ENDDO
           CALL Halt(' Logical error 3 in BinarySearch ')
         ENDIF
      END FUNCTION BinarySearch
!================================================================================
!     Compute distance criteria for controling penetration errors
!================================================================================
      RECURSIVE SUBROUTINE ExpandBoxWalk(Node,Sign)
         TYPE(RhoNode), POINTER :: Node
         REAL(DOUBLE)           :: Extent,Sign
!--------------------------------------------------------------------------
         Extent=GaussianExtent(Node%ZetaMin,Node%MaxAmp)        
         Node%Box=ExpandBox(Node%Box,Sign*Extent)
         IF(Node%Leaf)RETURN
         CALL ExpandBoxWalk(Node%Descend,Sign)
         CALL ExpandBoxWalk(Node%Descend%Travrse,Sign)
      END SUBROUTINE ExpandBoxWalk
!===============================================================================
!     Sums the significant contributions (leaves) to the density at point R
!===============================================================================
      RECURSIVE SUBROUTINE RhoOnGrid(Node)
         TYPE(RhoNode)                              :: Node
         REAL(DOUBLE)                               :: Tx,Ty,Tz
         REAL(DOUBLE), DIMENSION(0:MaxEll)          :: LLambdaX,LLambdaY,LLambdaZ, &
                                                       ULambdaX,ULambdaY,ULambdaZ, &
                                                       LambdaX,LambdaY,LambdaZ
         REAL(DOUBLE)                               :: LQx,LQy,LQz,LQ2,UQx,UQy,UQz,UQ2, &
                                                       LXpt,UXpt,TwoZ,SqZ,CoFact,RL1,TmpX,TmpY, &
                                                       LXptX,LXptY,LXptZ,UXptX,UXptY,UXptZ,Z
         INTEGER                                    :: IQ,IC,JQ,JC,KQ,KC,L,Ell,L1,L2,M,N,LMN
!-------------------------------------------------------------------------------
         Tx=ABS(Node%Box%Center(1)-Box%Center(1))
         IF(Tx>Node%Box%Half(1)+Box%Half(1))RETURN
         Ty=ABS(Node%Box%Center(2)-Box%Center(2))
         IF(Ty>Node%Box%Half(2)+Box%Half(2))RETURN
         Tz=ABS(Node%Box%Center(3)-Box%Center(3))
         IF(Tz>Node%Box%Half(3)+Box%Half(3))RETURN
         IF(Node%Leaf)THEN
            DO IGrid=1,NGrid
               CALL RhoAtAPoint(Node)
            ENDDO
            INCLUDE 'ExplicitLeafPopulation.Inc'         
            RETURN
!
             Z=Node%Zeta
             LQx=Box%BndBox(1,1)-Node%Qx
             LQy=Box%BndBox(2,1)-Node%Qy
             LQz=Box%BndBox(3,1)-Node%Qz
             UQx=Box%BndBox(1,2)-Node%Qx
             UQy=Box%BndBox(2,2)-Node%Qy
             UQz=Box%BndBox(3,2)-Node%Qz
             LQ2=LQx*LQx+LQy*LQy+LQz*LQz
             UQ2=UQx*UQx+UQy*UQy+UQz*UQz
             TwoZ=Two*Z
             SqZ=SQRT(Z)
             CoFact=SqrtPi/(Two*SqZ)
             LLambdaX(0)=-CoFact*ERF(SqZ*LQx)
             LLambdaY(0)=-CoFact*ERF(SqZ*LQy)
             LLambdaZ(0)=-CoFact*ERF(SqZ*LQz)
             ULambdaX(0)=-CoFact*ERF(SqZ*UQx)
             ULambdaY(0)=-CoFact*ERF(SqZ*UQy)
             ULambdaZ(0)=-CoFact*ERF(SqZ*UQz)
             LLambdaX(1)=TwoZ*LQx
             LLambdaY(1)=TwoZ*LQy
             LLambdaZ(1)=TwoZ*LQz
             ULambdaX(1)=TwoZ*UQx
             ULambdaY(1)=TwoZ*UQy
             ULambdaZ(1)=TwoZ*UQz
             DO L=2,Node%Ell
                L1=L-1
                L2=L-2
                RL1=DBLE(L1)
                LLambdaX(L)=TwoZ*(LQx*LLambdaX(L1)-RL1*LLambdaX(L2))
                LLambdaY(L)=TwoZ*(LQy*LLambdaY(L1)-RL1*LLambdaY(L2))
                LLambdaZ(L)=TwoZ*(LQz*LLambdaZ(L1)-RL1*LLambdaZ(L2))
                ULambdaX(L)=TwoZ*(UQx*ULambdaX(L1)-RL1*ULambdaX(L2))
                ULambdaY(L)=TwoZ*(UQy*ULambdaY(L1)-RL1*ULambdaY(L2))
                ULambdaZ(L)=TwoZ*(UQz*ULambdaZ(L1)-RL1*ULambdaZ(L2))
             ENDDO
             DO L=Node%Ell,1,-1
                LLambdaX(L)=LLambdaX(L-1)*LXptX
                LLambdaY(L)=LLambdaY(L-1)*LXptY
                LLambdaZ(L)=LLambdaZ(L-1)*LXptZ
                ULambdaX(L)=ULambdaX(L-1)*UXptX
                ULambdaY(L)=ULambdaY(L-1)*UXptY
                ULambdaZ(L)=ULambdaZ(L-1)*UXptZ
             ENDDO
             SqZ=SQRT(Z)
             CoFact=SqrtPi/(Two*SqZ)
             LLambdaX(0)=-CoFact*ERF(SqZ*LQx)
             LLambdaY(0)=-CoFact*ERF(SqZ*LQy)
             LLambdaZ(0)=-CoFact*ERF(SqZ*LQz)
             ULambdaX(0)=-CoFact*ERF(SqZ*UQx)
             ULambdaY(0)=-CoFact*ERF(SqZ*UQy)
             ULambdaZ(0)=-CoFact*ERF(SqZ*UQz)
             DO L=0,Node%Ell
                TmpX=LLambdaX(L)-ULambdaX(L)
                DO M=0,Node%Ell-L
                   TmpY=TmpX*(LLambdaY(M)-ULambdaY(M))
                   DO N=0,Node%Ell-M-L
                      LMN=LMNDex(L,M,N)
                      Pop=Pop+TmpY*(LLambdaZ(N)-ULambdaZ(N))*Node%Co(LMN)
                   ENDDO
                ENDDO
             ENDDO
         ELSE
            CALL RhoOnGrid(Node%Descend)
            CALL RhoOnGrid(Node%Descend%Travrse)
         ENDIF
      END SUBROUTINE RhoOnGrid
!===============================================================================
!     Sums the significant contributions (leaves) to the density at point R
!===============================================================================
      RECURSIVE SUBROUTINE RhoAtAPoint(Node)
         TYPE(RhoNode),         INTENT(IN)    :: Node
         REAL(DOUBLE), DIMENSION(0:MaxEll+1)  :: LambdaX,LambdaY,LambdaZ
         REAL(DOUBLE)                         :: Tx,Ty,Tz
         REAL(DOUBLE)                         :: Z,RQx,RQy,RQz,RQ2,Xpt,TwoZ,Co,RL1,X
         INTEGER                              :: IC,IQ,JC,JQ,KC,KQ,NQ,Ell, &
                                                 I,L,M,N,L1,L2,LMN,LMNLen,J
!-------------------------------------------------------------------------------
         IF(Grid(1,IGrid)<Node%Box%BndBox(1,1).OR. &
            Grid(2,IGrid)<Node%Box%BndBox(2,1).OR. &
            Grid(3,IGrid)<Node%Box%BndBox(3,1))RETURN
         IF(Grid(1,IGrid)>Node%Box%BndBox(1,2).OR. &
            Grid(2,IGrid)>Node%Box%BndBox(2,2).OR. &
            Grid(3,IGrid)>Node%Box%BndBox(3,2))RETURN
         IF(Node%Leaf)THEN            
#ifdef EXPLICIT_SOURCE
            INCLUDE 'ExplicitLeafContribution.Inc'
#else
            RQx=Grid(1,IGrid)-Node%Qx            
            RQy=Grid(2,IGrid)-Node%Qy
            RQz=Grid(3,IGrid)-Node%Qz
            RQ2=RQx*RQx+RQy*RQy+RQz*RQz
            Xpt=EXP(-Node%Zeta*RQ2)
            IF(Xpt>1.D-16)THEN
               TwoZ=Two*Node%Zeta
               LambdaX(0)=Xpt
               LambdaY(0)=One
               LambdaZ(0)=One
               LambdaX(1)=TwoZ*RQx*Xpt
               Lambday(1)=TwoZ*RQy
               LambdaZ(1)=TwoZ*RQz
               DO L=2,Node%Ell+1
                  L1=L-1
                  L2=L-2
                  RL1=DBLE(L1)
                  LambdaX(L)=TwoZ*(RQx*LambdaX(L1)-RL1*LambdaX(L2))
                  LambdaY(L)=TwoZ*(RQy*LambdaY(L1)-RL1*LambdaY(L2))
                  LambdaZ(L)=TwoZ*(RQz*LambdaZ(L1)-RL1*LambdaZ(L2))
               ENDDO
               DO L=0,Node%Ell
                  DO M=0,Node%Ell-L
                     DO N=0,Node%Ell-M-L
                        LMN=LMNDex(L,M,N)
                        RhoV(1,IGrid)=RhoV(1,IGrid)+LambdaX(L  )*LambdaY(M  )*LambdaZ(N  )*Node%Co(LMN)
                        RhoV(2,IGrid)=RhoV(2,IGrid)-LambdaX(L+1)*LambdaY(M  )*LambdaZ(N  )*Node%Co(LMN)
                        RhoV(3,IGrid)=RhoV(3,IGrid)-LambdaX(L  )*LambdaY(M+1)*LambdaZ(N  )*Node%Co(LMN)
                        RhoV(4,IGrid)=RhoV(4,IGrid)-LambdaX(L  )*LambdaY(M  )*LambdaZ(N+1)*Node%Co(LMN)
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
#endif
         ELSE
            CALL RhoAtAPoint(Node%Descend)
            CALL RhoAtAPoint(Node%Descend%Travrse)
         ENDIF
      END SUBROUTINE RhoAtAPoint
!=====================================================================================
!
!=====================================================================================
      SUBROUTINE FillRhoLeaf(Node)
         TYPE(RhoNode), POINTER :: Node
         INTEGER                :: I,IQ,IC,J,K,KQ,KC,L,B,E,N,NQ,NC,LMNLen,LTmp,Status        
         REAL(DOUBLE) :: RhoSum
         Interface 
            SUBROUTINE DblIntSort77(N,X,Y,Ordr)
               USE DerivedTYPEs
               INTEGER,                  Intent(IN)    :: N,Ordr
               REAL(DOUBLE),DIMENSION(N),Intent(INOUT) :: X
               INTEGER,     DIMENSION(N),Intent(INOUT) :: Y
            END SUBROUTINE
         END Interface
!-------------------------------------------------------------------------------------
!        Set leaf logical
         Node%Leaf=.True.
!        POINTERs to boundaries in the ordered lists
         B=Node%Bdex
         E=Node%Edex
         NQ=E-B+1
         IF(NQ/=1.OR.B/=E)  &
            CALL Halt('Bad Logic in FillRhoLeaf ')
         KQ=Qdex(B)
         KC=Cdex(KQ)
         Node%Ell=Ldex(KQ)
         LMNLen=LHGTF(Node%Ell)
         Node%Qx=Rho%Qx%D(KQ)
         Node%Qy=Rho%Qy%D(KQ)
         Node%Qz=Rho%Qz%D(KQ)
         Node%Zeta=Zeta(KQ)
         ALLOCATE(Node%Co(1:LMNLen),STAT=Status)
         CALL IncMem(Status,0,LMNLen)
         Node%Co(1:LMNLen)=Rho%Co%D(KC:KC+LMNLen-1)
      END SUBROUTINE FillRhoLeaf
!==========================================================================
!     Initialize a new RhoNode
!==========================================================================
      SUBROUTINE NewRhoNode(Node,Level)
         TYPE(RhoNode), POINTER   :: Node
         INTEGER                  :: Level,I,Status        
         ALLOCATE(Node,STAT=Status)
         IF(Status/=SUCCEED)  &
            CALL Halt(' Node ALLOCATE failed in NewRhoNode ')
         Node%Leaf=.False.
         Node%Box%Tier=Level
         Node%Box%Number=RhoNodes
         RhoNodes=RhoNodes+1
         NULLIFY(Node%Travrse)
         NULLIFY(Node%Descend)
         NULLIFY(Node%Co)
      END SUBROUTINE NewRhoNode
!==========================================================================
!     Recusrively delete the density tree
!==========================================================================
      RECURSIVE SUBROUTINE DeleteRhoTree(Node)
         TYPE(RhoNode), POINTER   :: Node
         INTEGER                  :: Level,I,LMNLen,Status        
         IF(Node%Leaf)THEN
            DEALLOCATE(Node%Co,STAT=Status)
            CALL DecMem(Status,LHGTF(Node%Ell),0)
         ELSEIF(Associated(Node%Descend))THEN
            CALL DeleteRhoTree(Node%Descend%Travrse)
            NULLIFY(Node%Descend%Travrse) 
            CALL DeleteRhoTree(Node%Descend)
            NULLIFY(Node%Descend)
         ENDIF
         DEALLOCATE(Node,STAT=Status)
         IF(Status/=SUCCEED)  &
            CALL Halt(' NODE DEALLOCATION FAILED IN DeleteRhoTree')
      END SUBROUTINE DeleteRhoTree
!========================================================================================
!     ALLOCATE and read in the density, initalize global lists 
!========================================================================================
      SUBROUTINE InitRho(Args)
         TYPE(ARGMT)  :: Args
         INTEGER      :: z,oq,or,iq,NQ,Q,Ell,Status,I,IOS,LMNLen
         REAL(DOUBLE) :: Dummy
!----------------------------------------------------------------------------------------
!        Get the density
         Open(UNIT=Seq,FILE=TrixFile('Rho',Args,0),STATUS='OLD', &
              FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
         Read(UNIT=Seq,Err=202,IOSTAT=IOS)Rho%NExpt,Rho%NDist,Rho%NCoef
         CALL New(Rho%NQ  ,Rho%NExpt)
         CALL New(Rho%OffQ,Rho%NExpt)
         CALL New(Rho%OffR,Rho%NExpt)
         CALL New(Rho%Lndx,Rho%NExpt)
         CALL New(Rho%Expt,Rho%NExpt)
         CALL New(Rho%Qx,  Rho%NDist)
         CALL New(Rho%Qy,  Rho%NDist)
         CALL New(Rho%Qz,  Rho%NDist)
         CALL New(Rho%Co,  Rho%NCoef)
         Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%NQ%I  (i),i=1,Rho%NExpt)
         Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%OffQ%I(i),i=1,Rho%NExpt)
         Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%OffR%I(i),i=1,Rho%NExpt)
         Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%Lndx%I(i),i=1,Rho%NExpt)
         Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%Expt%D(i),i=1,Rho%NExpt)
         Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%Qx%D  (i),i=1,Rho%NDist)
         Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%Qy%D  (i),i=1,Rho%NDist)
         Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%Qz%D  (i),i=1,Rho%NDist)
         Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Dummy        ,i=1,Rho%NDist)
         Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%Co%D  (i),i=1,Rho%NCoef)
         Close(UNIT=Seq,STATUS='KEEP')
!-------------------------------------------------------------
!        ALLOCATE global lists
         ALLOCATE(Qdex(1:Rho%NDist),STAT=Status)
         CALL IncMem(Status,Rho%NDist,0)
         ALLOCATE(Cdex(1:Rho%NDist),STAT=Status)
         CALL IncMem(Status,Rho%NDist,0)
         ALLOCATE(Ldex(1:Rho%NDist),STAT=Status)
         CALL IncMem(Status,Rho%NDist,0)
         ALLOCATE(RList(1:Rho%NDist),STAT=Status)
         CALL IncMem(Status,0,Rho%NDist)
         ALLOCATE(Zeta(1:Rho%NDist),STAT=Status)
         CALL IncMem(Status,0,Rho%NDist)
         ALLOCATE(Amp(1:Rho%NDist),STAT=Status)
         CALL IncMem(Status,0,Rho%NDist)
!
!        Fill the global indeces
!
         IQ=1
         MaxAmp=Zero
         DO z=1,Rho%NExpt-1
            oq  =Rho%OffQ%I(z)   
            or  =Rho%OffR%I(z)   
            Ell =Rho%Lndx%I(z)   
            LMNLen=LHGTF(Ell)
            NQ  =Rho%NQ%I(z)
            DO Q=1,NQ
               Zeta(IQ)=Rho%Expt%D(z)
               Qdex(IQ)=oq+Q
               Cdex(IQ)=or+(Q-1)*LMNLen +1
               Ldex(IQ)=Ell
               Amp(IQ)=Zero
               DO I=Cdex(IQ),Cdex(IQ)+LMNLen-1
                  Amp(IQ)=Amp(IQ)+Rho%Co%D(I)**2
               ENDDO
               Amp(IQ)=Log(Sqrt(Amp(IQ))+1.0D-100)
               MaxAmp=Max(MaxAmp,Amp(IQ))
               IQ=IQ+1
            ENDDO
         ENDDO        
!        Redefine NDist to exclude nuclear charges
         Rho%NDist=IQ-1
!        Later 
         Return            
!        Bomb on IO error
    202  CALL Halt('Died in PutRho, IOSTAT = '//Trim(IntToChar(IOS)))
      END SUBROUTINE InitRho
!========================================================================================
!     Delete globals associated with the array representation of the density
!========================================================================================
      SUBROUTINE DeleteRho
         INTEGER :: Status
         CALL Delete(Rho%NQ)
         CALL Delete(Rho%OffQ)
         CALL Delete(Rho%OffR)
         CALL Delete(Rho%Lndx)
         CALL Delete(Rho%Expt)
         CALL Delete(Rho%Qx)
         CALL Delete(Rho%Qy)
         CALL Delete(Rho%Qz)
         CALL Delete(Rho%Co)
!        DEALLOCATE global allocatables
         DEALLOCATE(Qdex,STAT=Status)
         CALL DecMem(Status,Rho%NDist,0)
         DEALLOCATE(Cdex,STAT=Status)
         CALL DecMem(Status,Rho%NDist,0)
         DEALLOCATE(Ldex,STAT=Status)
         CALL DecMem(Status,Rho%NDist,0)
         DEALLOCATE(Zeta,STAT=Status)
         CALL DecMem(Status,0,Rho%NDist)
         DEALLOCATE(Amp,STAT=Status)
         CALL DecMem(Status,0,Rho%NDist)
         DEALLOCATE(RList,STAT=Status) 
         CALL DecMem(Status,0,Rho%NDist)
      END SUBROUTINE DeleteRho
#ifdef MULTIPLE_DIST
!===================================================================
!     Check RhoBox agains contents
!===================================================================
      SUBROUTINE CheckBounds(Node)
         TYPE(RhoNode), POINTER :: Node
         REAL(DOUBLE)           :: ZMin,MaxAm,AM
         REAL(DOUBLE), Parameter :: Eps=1.D-8
         INTEGER :: I,IQ,JQ,JC,KQ,KC,L,LQ,LMNLen
         CALL PrintBBox(Node%Box,Node%ZetaMin,Node%MaxAmp)
         IF(Node%Leaf)THEN
         ZMin=1.D12
         MaxAm=-1.D12
         DO L=0,MaxEll
           JQ=Node%Qdex(L)
           JC=Node%Cdex(L)
           LMNLen=LHGTF(L)
           DO IQ=1,Node%NEll(L)                     
              KQ=JQ+IQ
              KC=JC+(IQ-1)*LMNLen                   
              ZMin=Min(ZMin,Node%Zeta(KQ))
              Am=Zero
              DO I=KC+1,KC+LMNLen
                 Am=Am+Node%Co(I)**2
              ENDDO
              Am=Log(Sqrt(Am))
              MaxAm=Max(MaxAm,Am)
              Write(22,33)Node%Zeta(KQ),Am,Node%Qx(KQ),Node%Qy(KQ),Node%Qz(KQ)
!             Tests
           33 Format(' Z = ',F10.4,', Amp = ',F10.4,', Q = ',F10.7,', ',F10.7,', ',F10.7)
              IF(Node%Qx(KQ)+Eps<Node%Box%BndBox(1,1).Or.Node%Qx(KQ)-Eps>Node%Box%BndBox(1,2))THEN
                 Write(*,*)Node%Qx(KQ)<Node%Box%BndBox(1,1),Node%Qx(KQ)>Node%Box%BndBox(1,2),&
                           Node%Box%BndBox(1,1),Node%Qx(KQ),Node%Box%BndBox(1,2)
                 Stop 'Bounding box hosed !'
              ENDIF
              IF(Node%Qy(KQ)+Eps<Node%Box%BndBox(2,1).Or.Node%Qy(KQ)-Eps>Node%Box%BndBox(2,2))THEN
                 Write(*,*)Node%Qy(KQ)<Node%Box%BndBox(2,1),Node%Qy(KQ)>Node%Box%BndBox(2,2), &
                           Node%Box%BndBox(2,1),Node%Qy(KQ),Node%Box%BndBox(2,2)
                 Stop 'Bounding box hosed !'
              ENDIF
              IF(Node%Qz(KQ)+Eps<Node%Box%BndBox(3,1).Or.Node%Qz(KQ)-Eps>Node%Box%BndBox(3,2))THEN
                 Write(*,*)Node%Qz(KQ)<Node%Box%BndBox(3,1),Node%Qz(KQ)>Node%Box%BndBox(3,2),&
                           Node%Box%BndBox(3,1),Node%Qz(KQ),Node%Box%BndBox(3,2)
                 Stop 'Bounding box hosed !'
              ENDIF 
           ENDDO
         ENDDO
         IF(Node%ZetaMin>ZMin)THEN
            Write(22,*)Node%ZetaMin,ZMin,Node%ZetaMin>ZMin
            Stop ' zeta '
          ENDIF
         IF(Node%MaxAmp+Eps<MaxAm)THEN
            Write(22,*)Node%MaxAmp,MaxAm,Node%MaxAmp<MaxAm
            Stop ' amp '
         ENDIF
         ENDIF
      END SUBROUTINE CheckBounds  
#endif
!========================================================================================
!     Compute the erf function
!========================================================================================
      FUNCTION ERF(W)
         REAL(DOUBLE), INTENT(IN) :: W
         REAL(DOUBLE)             :: ERF,X,Sgn
         INTEGER                  :: I,J
         IF(W<0.0D0)THEN
            Sgn=-1.0D0
         ELSE
            Sgn=1.0D0 
         ENDIF        
         X=Sgn*W
         IF(X.GE.Erf_Switch)THEN
            Erf=1.0D0
         ELSE
            J=AINT(X*Erf_Grid)
            Erf=Erf_0(J)+X*(Erf_1(J)+X*(Erf_2(J)+X*Erf_3(J))) 
         ENDIF
         Erf=Sgn*Erf
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
END Module
