MODULE PoleTree
   USE Derivedtypes
   USE GlobalScalars   
   USE GlobalObjects
   USE ProcessControl
   USE Indexing
   USE Parse
   USE InOut
   USE Macros
   USE Thresholding
   USE BoundingBox
   USE MondoPoles
   IMPLICIT NONE
!----------------------------------------------------------------------------------
!  Global parameters
!
!  Global scalars
   INTEGER            :: PoleNodes,RhoLevel,MaxTier,NTier
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
   TYPE(PoleNode), POINTER :: PoleRoot  ! Root of the tree 
!-----------!
   CONTAINS !
!=================================================================================
!     
!=================================================================================
      SUBROUTINE RhoToPoleTree(Args)
         TYPE(ARGMT)               :: Args
         INTEGER                   :: Status,K,I
!-------------------------------------------------------------------
!        Initialize the density
         CALL InitRho(Args)
!        Initialize counters
         PoleNodes=0
         RhoLevel=0
!        Initialize the root node
         CALL NewPoleNode(PoleRoot,0)
         CALL InitPoleRoot
!        Convert the density into a 3-D BinTree
         CALL SplitPole(PoleRoot)
!        Delete the global array based density
         CALL DeleteRho
!        Make PoleTree tier by tier, recuring up from the bottom
         DO NTier=MaxTier,0,-1         
            CALL MakePoleTree(PoleRoot) 
         ENDDO 
        CALL Print_PoleNode(PoleRoot,'Root')
      END SUBROUTINE RhoToPoleTree
!==========================================================================
!
!==========================================================================
      RECURSIVE SUBROUTINE MakePoleTree(Node)
         TYPE(PoleNode) :: Node
!--------------------------------------------------------------
         IF(Node%Box%Tier<NTier.AND.(.NOT.Node%Leaf))THEN
            CALL MakePoleTree(Node%Descend)
            CALL MakePoleTree(Node%Descend%Travrse)
         ELSEIF(Node%Box%Tier==NTier)THEN
            IF(Node%Leaf)THEN
               Node%C=Zero; Node%S=Zero
               CALL HGToSP(Node)
            ELSE            
               Node%C=Zero; Node%S=Zero
!              Translation from left Q to P center
               CALL XLate(Node%Descend,Node)
!              Translation from right Q to P center
               CALL XLate(Node%Descend%Travrse,Node)
            ENDIF
         ENDIF
     END SUBROUTINE MakePoleTree
!===============================================================================
!     
!=======================================================================
      SUBROUTINE InitPoleRoot
         TYPE(DBL_RNK2) :: BndBox
         PoleRoot%Bdex=1
         PoleRoot%Edex=Rho%NDist
         PoleRoot%NQ=Rho%NDist
!        Get the nuclear bounding box
         CALL New(BndBox,(/3,2/))
         CALL Get(BndBox,'boundingbox',CurGeom)
         PoleRoot%Box%BndBox(1:3,1:2)=BndBox%D(1:3,1:2)
         CALL Delete(BndBox)
!        Set the largest range (smallest exponent) of the roots BB
         PoleRoot%Zeta=Rho%Expt%D(1)
         MiniumExp=Rho%Expt%D(1)
!        Set the center and width of the Cartesian part
         PoleRoot%Box%Half  =Half*(PoleRoot%Box%BndBox(:,2)-PoleRoot%Box%BndBox(:,1))
         PoleRoot%Box%Center=Half*(PoleRoot%Box%BndBox(:,2)+PoleRoot%Box%BndBox(:,1))
     END SUBROUTINE InitPoleRoot
!===================================================================
!
!===================================================================
      RECURSIVE SUBROUTINE SplitPole(Node)
         TYPE(PoleNode), POINTER :: Node,Left,Right
!--------------------------------------------------------------
         IF(Node%NQ==1)THEN
            CALL FillRhoLeaf(Node)
         ELSE 
!           Allocate new children 
            CALL NewPoleNode(Node%Descend,Node%Box%Tier+1)
            CALL NewPoleNode(Node%Descend%Travrse,Node%Box%Tier+1)
!           Set links
            Node%Descend%Travrse%Travrse=>Node%Travrse
            Left=>Node%Descend
            Right=>Node%Descend%Travrse
            CALL SplitPoleBox(Node,Left,Right)
!           Recur
            CALL SplitPole(Left)
            CALL SplitPole(Right)
!           Min exponent
            Node%Zeta=MIN(Left%Zeta,Right%Zeta)
!           Max distance to the BBox surface squared
            Node%D2=Node%Box%Half(1)**2+Node%Box%Half(2)**2+Node%Box%Half(3)**2+1.D-16
         ENDIF
       END SUBROUTINE SplitPole
!=====================================================================================
!
!=====================================================================================
      SUBROUTINE FillRhoLeaf(Node)
         TYPE(PoleNode), POINTER :: Node
         INTEGER                 :: I,IQ,IC,J,K,KQ,KC,L,B,E,N,NQ,NC,LMNLen,Status        
         REAL(DOUBLE)            :: RhoSum
!-------------------------------------------------------------------------------------
!        Set leaf logical
         Node%Leaf=.True.
!        POINTERs to boundaries in the ordered lists
         B=Node%Bdex
         E=Node%Edex
         NQ=E-B+1
         IF(NQ/=1.OR.B/=E)CALL Halt('Bad Logic in FillRhoLeaf ')
         KQ=Qdex(B)
         KC=Cdex(KQ)
         Node%Box%Center(1)=Rho%Qx%D(KQ)
         Node%Box%Center(2)=Rho%Qy%D(KQ)
         Node%Box%Center(3)=Rho%Qz%D(KQ)
         Node%Zeta=Zeta(KQ)
         Node%Ell=Ldex(KQ)
         LMNLen=LHGTF(Node%Ell)
         ALLOCATE(Node%Co(1:LMNLen),STAT=Status)
         CALL IncMem(Status,0,LMNLen)
         Node%Co(1:LMNLen)=Rho%Co%D(KC:KC+LMNLen-1)
         Node%Bdex=MAX(0,KQ-SUM(Rho%NQ%I(1:Rho%NExpt-1)))
!         IF(Node%Bdex>0) &
!         WRITE(*,*)SUM(Rho%NQ%I(1:Rho%NExpt-1)),KQ,' B = ',Node%Bdex,' Co = ',Node%Co(1)
      END SUBROUTINE FillRhoLeaf
!==========================================================================
!     Initialize a new PoleNode
!==========================================================================
      SUBROUTINE NewPoleNode(Node,Level)
         TYPE(PoleNode), POINTER   :: Node
         INTEGER                  :: Level,I,Status        
         ALLOCATE(Node,STAT=Status)
         IF(Status/=SUCCEED)CALL Halt(' Node ALLOCATE failed in NewPoleNode ')
         Node%Leaf=.False.
         Node%Box%Tier=Level
         Node%Box%Number=PoleNodes
         Node%Ell=SPEll
         MaxTier=MAX(MaxTier,Level)
         PoleNodes=PoleNodes+1
         NULLIFY(Node%Travrse)
         NULLIFY(Node%Descend)
         ALLOCATE(Node%S(0:SPLen),STAT=Status)
         CALL IncMem(Status,0,SPLen+1)
         ALLOCATE(Node%C(0:SPLen),STAT=Status)
         CALL IncMem(Status,0,SPLen+1)
         Node%S(0:SPLen)=Zero
         Node%C(0:SPLen)=Zero
      END SUBROUTINE NewPoleNode
!================================================================================
!     New RhoBox
!================================================================================
      SUBROUTINE NewRhoBox(Node)
         TYPE(PoleNode), POINTER :: Node
         INTEGER                :: J,I,IQ
         REAL(DOUBLE)           :: Test
         Node%Zeta=1.D12
         Node%Box%BndBox(:,1)= 1.D12
         Node%Box%BndBox(:,2)=-1.D12
         DO J=Node%Bdex,Node%Edex 
            IQ=Qdex(J)
            Node%Box%BndBox(1,1)=Min(Node%Box%BndBox(1,1),Rho%Qx%D(IQ))
            Node%Box%BndBox(1,2)=Max(Node%Box%BndBox(1,2),Rho%Qx%D(IQ))
            Node%Box%BndBox(2,1)=Min(Node%Box%BndBox(2,1),Rho%Qy%D(IQ))
            Node%Box%BndBox(2,2)=Max(Node%Box%BndBox(2,2),Rho%Qy%D(IQ))
            Node%Box%BndBox(3,1)=Min(Node%Box%BndBox(3,1),Rho%Qz%D(IQ))
            Node%Box%BndBox(3,2)=Max(Node%Box%BndBox(3,2),Rho%Qz%D(IQ))
         ENDDO
!        New sides
         Node%Box%Half(:)=Half*(Node%Box%BndBox(:,2)-Node%Box%BndBox(:,1))
!        New center
         Node%Box%Center(:)=Half*(Node%Box%BndBox(:,2)+Node%Box%BndBox(:,1))
      END SUBROUTINE NewRhoBox
!========================================================================================
!     ALLOCATE and read in the density, initalize global lists 
!========================================================================================
      SUBROUTINE InitRho(Args)
         TYPE(ARGMT)  :: Args
         INTEGER      :: z,oq,or,iq,jq,NQ,Q,Ell,Status,I,IOS,LMNLen
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
!        Fill the global indeces
         IQ=1
         JQ=1
         MaxAmp=Zero
         DO z=1,Rho%NExpt!,Rho%NExpt
            oq  =Rho%OffQ%I(z)   
            or  =Rho%OffR%I(z)   
            Ell =Rho%Lndx%I(z)   
            LMNLen=LHGTF(Ell)
            NQ  =Rho%NQ%I(z)
            DO Q=1,NQ
               Qdex(IQ)=oq+Q
               Cdex(oq+Q)=or+(Q-1)*LMNLen +1 
               Zeta(oq+Q)=Rho%Expt%D(z)
               Ldex(oq+Q)=Ell
               IQ=IQ+1
               JQ=JQ+LMNLen
            ENDDO
         ENDDO     
!        Redefine NDist to exclude nuclear charges
         Rho%NDist=IQ-1
!        Later 
         Return            
!        Bomb on IO error
    202  CALL Halt('Died in PutRho, IOSTAT = '//Trim(IntToChar(IOS)))
      END SUBROUTINE InitRho
!================================================================================
!     Bisection   
!================================================================================
      SUBROUTINE SplitPoleBox(Node,Left,Right)
         TYPE(PoleNode), POINTER :: Node,Left,Right
         REAL(DOUBLE)            :: Section
         INTEGER                 :: B,E,N,ISplit,Split,I,J,k
!--------------------------------------------------------------
!        Indexing
         J=0
         B=Node%Bdex
         E=Node%Edex
         N=E-B+1
!        Orthogonal direction
         ISplit=Mod(Node%Box%Tier,3)+1
         IF(ISplit==1)THEN
            DO I=B,E
               J=J+1;K=Qdex(I)
               RList(J)=Rho%Qx%D(K)
            ENDDO
         ELSEIF(ISplit==2)THEN
            DO I=B,E
               J=J+1;K=Qdex(I)
               RList(J)=Rho%Qy%D(K)
            ENDDO
         ELSEIF(ISplit==3)THEN
            DO I=B,E
               J=J+1;K=Qdex(I)
               RList(J)=Rho%Qz%D(K)
            ENDDO
         ENDIF
!        Sort
         CALL DblIntSort77(N,RList,Qdex(B:E),2)
!        Orthogonal RECURSIVE bisection (ORB)
         Section   =RList(1)+Half*(RList(N)-RList(1))
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
      END SUBROUTINE SplitPoleBox
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
!
END MODULE
