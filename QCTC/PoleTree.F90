MODULE PoleTree
   USE Derivedtypes
   USE GlobalScalars   
   USE GlobalObjects
   USE ProcessControl
   USE Indexing
   USE Parse
   USE InOut
   USE Macros
   USE QCTCThresholds
   USE BoundingBox
   USE MondoPoles
   USE Globals
   IMPLICIT NONE
!----------------------------------------------------------------------------------
!  Globals
   TYPE(PoleNode), POINTER               :: PoleRoot  ! Root of the pole tree 
   INTEGER                               :: PoleNodes
   INTEGER                               :: RhoLevel
   INTEGER                               :: CurrentTier
   INTEGER                               :: MaxTier
   INTEGER                               :: NElecDist
   INTEGER,     DIMENSION(:),Allocatable :: Qdex      
   INTEGER,     DIMENSION(:),Allocatable :: Cdex
   INTEGER,     DIMENSION(:),Allocatable :: Ldex
   REAL(DOUBLE),DIMENSION(:),Allocatable :: RList      
   REAL(DOUBLE),DIMENSION(:),Allocatable :: Zeta
   REAL(DOUBLE),DIMENSION(:),Allocatable :: Ext
!-----------!
   CONTAINS !
!=================================================================================
!     
!=================================================================================
      SUBROUTINE RhoToPoleTree
        INTEGER                   :: Status,K,I
!-------------------------------------------------------------------
!       Initialize counters
        PoleNodes=0
        RhoLevel=0
!       Initialize the root node
        CALL NewPoleNode(PoleRoot,0)
        CALL NewSPArrays(PoleRoot)
        PoleRoot%Bdex=1
        PoleRoot%Edex=Rho%NDist
        PoleRoot%NQ=Rho%NDist
!       Convert the density into a 3-D BinTree
        CALL SplitPole(PoleRoot)
!       Make PoleTree tier by tier, recuring up from the bottom
        DO CurrentTier=MaxTier,0,-1         
           CALL MakePoleTree(PoleRoot) 
        ENDDO
!       Reset Ell of PoleRoot
        PoleRoot%Ell=SPEll
!        CALL Print_PoleNode(PoleRoot,'Root')
      END SUBROUTINE RhoToPoleTree
!=====================================================================================
!
!=====================================================================================
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
!           Allocate multipole arrays
            CALL NewSPArrays(Left)
            CALL NewSPArrays(Right)
         ENDIF
       END SUBROUTINE SplitPole
!=================================================================================
!     
!=================================================================================
      RECURSIVE SUBROUTINE MakePoleTree(P)
         TYPE(PoleNode)                   :: P
         TYPE(PoleNode),POINTER           :: LeftQ,RightQ
         INTEGER                          :: K
         REAL(DOUBLE)                     :: RightDist,LeftDist
         REAL(DOUBLE),PARAMETER           :: SPEllPlus2=SPEll+2
         REAL(DOUBLE),PARAMETER           :: UnsoldExp=Two/SPEllPlus2
!--------------------------------------------------------------
         IF(P%Box%Tier==CurrentTier.AND.P%Leaf)THEN
            P%C=Zero
            P%S=Zero
            P%Strength=Zero
            CALL HGToSP(P)
         ELSEIF(P%Leaf)THEN
            RETURN
         ELSEIF(P%Box%Tier==CurrentTier)THEN 
            P%C=Zero 
            P%S=Zero
            LeftQ=>P%Descend
            RightQ=>P%Descend%Travrse
!           Compute new nodes BBox
            DO K=1,3
               P%Box%BndBox(K,1)=MIN(LeftQ%Box%BndBox(K,1),RightQ%Box%BndBox(K,1))
               P%Box%BndBox(K,2)=MAX(LeftQ%Box%BndBox(K,2),RightQ%Box%BndBox(K,2))
            ENDDO
            P%Box%Half(:)  = Half*(P%Box%BndBox(:,2)-P%Box%BndBox(:,1))
            P%Box%Center(:)= Half*(P%Box%BndBox(:,2)+P%Box%BndBox(:,1))
!           Compute new nodes DBox
            DO K=1,3
               P%DBox%BndBox(K,1)=MIN(LeftQ%DBox%BndBox(K,1),RightQ%DBox%BndBox(K,1))
               P%DBox%BndBox(K,2)=MAX(LeftQ%DBox%BndBox(K,2),RightQ%DBox%BndBox(K,2))
            ENDDO
            P%DBox%Half(:)  = Half*(P%DBox%BndBox(:,2)-P%DBox%BndBox(:,1))
            P%DBox%Center(:)= Half*(P%DBox%BndBox(:,2)+P%DBox%BndBox(:,1))
!           Compute DMax
            P%DMax2 = Zero
            DO K=1,3
               P%DMax2 = P%DMax2+(P%DBox%Center(K)-P%DBox%BndBox(K,2))**2
            ENDDO
!           Compute Zeta
            P%Zeta=MIN(LeftQ%Zeta,RightQ%Zeta)
!-----------------------------------------------------------------------------
!           Translate Left and Right with SPEll+1
!           Translate LeftQ-> P
            CALL XLate(LeftQ,P)
!           Translate RightQ-> P
            CALL XLate(RightQ,P)
!           Reset Ell of Right and Left nodes if not leafs
            IF(.NOT.LeftQ%Leaf)   LeftQ%Ell=SPell
            IF(.NOT.RightQ%Leaf) RightQ%Ell=SPell
!           Compute the multipole strength [O^P_(L+1)]^(2/(2+L))
            P%Strength=Unsold1(SPEll+1,P%C,P%S)**UnsoldExp
         ELSE
!           Keep on truckin ...
            CALL MakePoleTree(P%Descend)
            CALL MakePoleTree(P%Descend%Travrse)
         ENDIF
     END SUBROUTINE MakePoleTree
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
         Node%Zeta=Zeta(KQ)
!        Reset leaf nodes ell to min value: its untranslated!
         Node%Ell=Ldex(KQ)
!        Set and inflate this nodes BBox
         Node%Box%BndBox(1,:)=Rho%Qx%D(KQ)
         Node%Box%BndBox(2,:)=Rho%Qy%D(KQ)
         Node%Box%BndBox(3,:)=Rho%Qz%D(KQ)
         Node%Box=ExpandBox(Node%Box,Ext(KQ))
!        Set this nodes DBBox
         Node%DBox%BndBox(1,:)=Rho%Qx%D(KQ)
         Node%DBox%BndBox(2,:)=Rho%Qy%D(KQ)
         Node%DBox%BndBox(3,:)=Rho%Qz%D(KQ)
         Node%DBox%Half(:)    =Zero
         Node%DBox%Center(:)  =Node%DBox%BndBox(:,1)     
         Node%DMax2           =Zero
!        Allocate and fill HGTF coefficients
         LMNLen=LHGTF(Node%Ell)
         ALLOCATE(Node%Co(1:LMNLen),STAT=Status)
         CALL IncMem(Status,0,LMNLen)
         KC=Cdex(KQ)
         Node%Co(1:LMNLen)=Rho%Co%D(KC:KC+LMNLen-1)
!        Mark node to identify nuclear self interaction
         Node%Bdex=MAX(0,KQ-NElecDist)
      END SUBROUTINE FillRhoLeaf
!==========================================================================
!     Initialize a new PoleNode
!==========================================================================
      SUBROUTINE NewPoleNode(Node,Level)
         TYPE(PoleNode), POINTER   :: Node
         INTEGER                   :: Level,I,Status        
         ALLOCATE(Node,STAT=Status)
         IF(Status/=SUCCEED)CALL Halt(' Node ALLOCATE failed in NewPoleNode ')
         Node%Leaf=.False.
         Node%Box%Tier=Level
         Node%Box%Number=PoleNodes
         Node%Ell=SPEll+1
         MaxTier=MAX(MaxTier,Level)
         PoleNodes=PoleNodes+1
         NULLIFY(Node%Travrse)
         NULLIFY(Node%Descend)
      END SUBROUTINE NewPoleNode
!==========================================================================
!     Initialize a new PoleNodes Array
!==========================================================================
      SUBROUTINE NewSPArrays(Node)
         TYPE(PoleNode), POINTER   :: Node
         INTEGER                   :: LenSP,Status        
         LenSP=LSP(Node%Ell)
         ALLOCATE(Node%S(0:LenSP),STAT=Status)
         CALL IncMem(Status,0,LenSP+1)
         ALLOCATE(Node%C(0:LenSP),STAT=Status)
         CALL IncMem(Status,0,LenSP+1)
         Node%S=Zero
         Node%C=Zero
      END SUBROUTINE NewSPArrays
!================================================================================
!     Bisection   
!================================================================================
      SUBROUTINE SplitPoleBox(Node,Left,Right)
         TYPE(PoleNode), POINTER :: Node,Left,Right
         REAL(DOUBLE)            :: Section,MaxBox
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
      END SUBROUTINE SplitPoleBox
!========================================================================================
!     ALLOCATE and read in the density, initalize global lists 
!========================================================================================
      SUBROUTINE InitRhoAux
         INTEGER      :: z,oq,or,iq,jq,NQ,Q,Ell,Status,I,IOS,LMNLen,CD,QD
         REAL(DOUBLE) :: ZE,EX
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
         ALLOCATE(Ext(1:Rho%NDist),STAT=Status)
         CALL IncMem(Status,0,Rho%NDist)
!        Fill the global indeces
         IQ=1
         DO z=1,Rho%NExpt
            oq =Rho%OffQ%I(z)   
            or =Rho%OffR%I(z)   
            Ell=Rho%Lndx%I(z)   
            ZE =Rho%Expt%D(z)
            LMNLen=LHGTF(Ell)
            DO Q=1,Rho%NQ%I(z)
               QD=oq+Q
               CD=or+(Q-1)*LMNLen+1
               EX=Extent(Ell,ZE,Rho%Co%D(CD:CD+LMNLen-1),Tau_O=TauPAC,Potential_O=.TRUE.,ExtraEll_O=0)
!              Threshold out distributions with zero extent 
               IF(EX>Zero)THEN
                  Qdex(IQ)=QD
                  Cdex(QD)=CD
                  Ext( QD)=EX
                  Zeta(QD)=ZE
                  Ldex(QD)=Ell
                  IQ=IQ+1
               ENDIF
            ENDDO
         ENDDO        
!        Recompute the number of distributions
         Rho%NDist=IQ-1
!        Number of distributions excluding nuclei (used for
!        identifiying nuclear selfinteraction, note that nuclei
!        must be in correct order for this to work...)
         NElecDist=SUM(Rho%NQ%I(1:Rho%NExpt-1))
       END SUBROUTINE InitRhoAux
!========================================================================================
!     Delete globals associated with the array representation of the density
!========================================================================================
       SUBROUTINE DeleteRhoAux
         INTEGER :: Status
!        DEALLOCATE global allocatables
         DEALLOCATE(Qdex,STAT=Status)
         CALL DecMem(Status,Rho%NDist,0)
         DEALLOCATE(Cdex,STAT=Status)
         CALL DecMem(Status,Rho%NDist,0)
         DEALLOCATE(Ldex,STAT=Status)
         CALL DecMem(Status,Rho%NDist,0)
         DEALLOCATE(Zeta,STAT=Status)
         CALL DecMem(Status,0,Rho%NDist)
         DEALLOCATE(Ext,STAT=Status)
         CALL DecMem(Status,0,Rho%NDist)
         DEALLOCATE(RList,STAT=Status) 
         CALL DecMem(Status,0,Rho%NDist)
       END SUBROUTINE DeleteRhoAux
!
END MODULE
