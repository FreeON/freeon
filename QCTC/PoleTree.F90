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
   TYPE(PoleNode), POINTER               :: PoleRoot ! Root of the pole tree 
   TYPE(PoleNode), POINTER               :: PR1
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
        MaxTier=0
!       Initialize the root node   
        CALL NewPoleNode(PoleRoot,0)
!        CALL NewSPArrays(PoleRoot) 
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
        PoleRoot%Ell     = SPEll
!
      END SUBROUTINE RhoToPoleTree
!=====================================================================================
!     
!=====================================================================================
      RECURSIVE SUBROUTINE CheckNodes(Node)
        TYPE(PoleNode), POINTER :: Node
!
        IF(Node%Leaf)THEN
           WRITE(*,*) 'Leaf Node',Node%Ell,Node%EllCD,Node%Box%Tier
           WRITE(*,*) Node%DMax2
           WRITE(*,*) Node%WCoef
           WRITE(*,*) Node%Box%BndBox(1,1:2)
           RETURN   
        ELSE
           WRITE(*,*) 'Pole Node',Node%Ell,Node%EllCD,Node%Box%Tier
           WRITE(*,*) Node%DMax2
           WRITE(*,*) Node%WCoef
           WRITE(*,*) Node%Box%BndBox(1,1:2)
!           CALL CheckNodes(Node%Descend)
           CALL CheckNodes(Node%Descend%Travrse)
        ENDIF
!
      END SUBROUTINE CheckNodes
!=====================================================================================
!
!=====================================================================================
      RECURSIVE SUBROUTINE SplitPole(Node)
         TYPE(PoleNode), POINTER :: Node,Left,Right
!--------------------------------------------------------------
         IF(Node%NQ==1)THEN
            CALL FillRhoLeaf(Node)
            CALL NewSPArrays(Node)
         ELSE 
!           Allocate new children 
            CALL NewPoleNode(Node%Descend,Node%Box%Tier+1)
            CALL NewPoleNode(Node%Descend%Travrse,Node%Box%Tier+1)
!           Set links
            Node%Descend%Travrse%Travrse=>Node%Travrse
            Left=>Node%Descend
            Right=>Node%Descend%Travrse
!           Compute the Bounding Box
!            CALL ComputeBoundingBox(Node)
!           Split the Box
            CALL SplitPoleBox(Node,Left,Right)
!           Recur
            CALL SplitPole(Left)
            CALL SplitPole(Right)
!           Allocate multipole arrays
            ! CALL NewSPArrays(Left)
            ! CALL NewSPArrays(Right)
            CALL NewSPArrays(Node)
         ENDIF
       END SUBROUTINE SplitPole
!=================================================================================
!     
!=================================================================================
      RECURSIVE SUBROUTINE MakePoleTree(P)
         TYPE(PoleNode)                        :: P
         TYPE(PoleNode),POINTER                :: LeftQ,RightQ
         INTEGER                               :: I,J,K,L,LL,LM
         REAL(DOUBLE)                          :: MaxUnsold,PMax
         REAL(DOUBLE)                          :: M1,M2,M3,CQ,Del
!--------------------------------------------------------------
         IF(P%Box%Tier==CurrentTier.AND.P%Leaf)THEN
            ! P%C=Zero
            ! P%S=Zero
            P%Strength=Zero
            CALL HGToSP(P)
         ELSEIF(P%Leaf)THEN
            RETURN
         ELSEIF(P%Box%Tier==CurrentTier)THEN 
            ! P%C=Zero 
            ! P%S=Zero
            LeftQ=>P%Descend
            RightQ=>P%Descend%Travrse
!           Compute new nodes BBox
            DO K=1,3
               P%Box%BndBox(K,1)=MIN(LeftQ%Box%BndBox(K,1),RightQ%Box%BndBox(K,1))
               P%Box%BndBox(K,2)=MAX(LeftQ%Box%BndBox(K,2),RightQ%Box%BndBox(K,2))
            ENDDO
            P%Box%Half(1)  = Half*(P%Box%BndBox(1,2)-P%Box%BndBox(1,1))
            P%Box%Half(2)  = Half*(P%Box%BndBox(2,2)-P%Box%BndBox(2,1))
            P%Box%Half(3)  = Half*(P%Box%BndBox(3,2)-P%Box%BndBox(3,1))
            P%Box%Center(1)= Half*(P%Box%BndBox(1,2)+P%Box%BndBox(1,1))
            P%Box%Center(2)= Half*(P%Box%BndBox(2,2)+P%Box%BndBox(2,1))
            P%Box%Center(3)= Half*(P%Box%BndBox(3,2)+P%Box%BndBox(3,1))
!           Compute DMax2
            P%DMax2 = Zero
            DO I=P%Bdex,P%Edex
               J=Qdex(I)
               PMax = (Rho%Qx%D(J)-P%Box%Center(1))**2+(Rho%Qy%D(J)-P%Box%Center(2))**2+(Rho%Qz%D(J)-P%Box%Center(3))**2
               P%DMax2=MAX(P%DMax2,PMax)
            ENDDO
!           Compute Zeta
            P%Zeta=MIN(LeftQ%Zeta,RightQ%Zeta)
!           Set the Ells and the Strengths
            IF(SQRT(P%DMax2) < 1.D-12) THEN
!              Reset Ell of Right and Left nodes if not leafs
               IF(.NOT.LeftQ%Leaf) THEN
                  LeftQ%Ell     = MIN(MAX(LeftQ%Descend%Ell,LeftQ%Descend%Travrse%Ell),SPEll)
               ENDIF
               IF(.NOT.RightQ%Leaf) THEN 
                  RightQ%Ell    = MIN(MAX(RightQ%Descend%Ell,RightQ%Descend%Travrse%Ell),SPEll)
               ENDIF
!              Accumulate the Multipoles in P
               DO LM=0,LSP(LeftQ%Ell)
                  P%C(LM) = LeftQ%C(LM)
                  P%S(LM) = LeftQ%S(LM)
               ENDDO
               DO LM=0,LSP(RightQ%Ell)
                  P%C(LM) = P%C(LM)+RightQ%C(LM)
                  P%S(LM) = P%S(LM)+RightQ%S(LM)
               ENDDO
!              Compute the multipole strength  MAX{[O^P_(L+1)]^(2/(2+L))}
               P%Strength=Zero
            ELSE
!              Translate the Nodes
!              Translate LeftQ-> P
               CALL XLate(LeftQ,P)
!              Translate RightQ-> P
               CALL XLate(RightQ,P)
!              Reset Ell of Right and Left nodes if not leafs
               IF(.NOT.LeftQ%Leaf)  LeftQ%Ell  = SPEll
               IF(.NOT.RightQ%Leaf) RightQ%Ell = SPEll 
!              Compute the multipole strength MAX{[O^P_(L+1)]^(2/(2+L))}
               CQ=Zero
               DO L=SPEll+1,SPEll+MaxUEll
                  CQ = MAX(CQ,Unsold0(L,P%C,P%S)/(P%DMax2**(Half*DBLE(L))))
               ENDDO
               MaxUnsold  = CQ*(P%DMax2**(Half*DBLE(SPEll+1)))
               P%Strength = MaxUnsold**(Two/DBLE(SPEll+2))
            ENDIF
#ifdef NewPAC
!           Accumulate Info for PAC from Left and Right Nodes
            P%WCoef = LeftQ%WCoef+RightQ%WCoef
            P%EllCD = MAX(LeftQ%EllCD,RightQ%EllCD)
#endif
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
         Node%Ell   = Ldex(KQ)
!        Set and inflate this nodes BBox
         Node%Box%BndBox(1,:)= Rho%Qx%D(KQ)
         Node%Box%BndBox(2,:)= Rho%Qy%D(KQ)
         Node%Box%BndBox(3,:)= Rho%Qz%D(KQ)
#ifdef NewPAC
         Node%Box=ExpandBox(Node%Box,1.D-16)
#else
         Node%Box=ExpandBox(Node%Box,Ext(KQ))
#endif
!        Allocate and fill HGTF coefficients
         LMNLen=LHGTF(Node%Ell)
         ALLOCATE(Node%Co(1:LMNLen),STAT=Status)
         CALL IncMem(Status,0,LMNLen)
         KC=Cdex(KQ)
         Node%Co(1:LMNLen)=Rho%Co%D(KC:KC+LMNLen-1)
!        Initialize the MAC Stuff to be Zero
         Node%DMax2    = Zero
         Node%Strength = Zero
#ifdef NewPAC
!        Accumulate Info for the New PAC
         Node%WCoef  = NodeWeight(Node%Ell,Node%Zeta,Node%Co)
         Node%EllCD  = Node%Ell
#endif
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
         Node%Ell  =SPEll+MaxUEll
#ifdef NewPAC
         Node%EllCD=0
#endif
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
         ! Node%S=Zero
         ! Node%C=Zero
         CALL DBL_VECT_EQ_DBL_SCLR(LenSP+1,Node%S(0),Zero)
         CALL DBL_VECT_EQ_DBL_SCLR(LenSP+1,Node%C(0),Zero)
      END SUBROUTINE NewSPArrays
!================================================================================
!     Compute the Bounding Box and 
!     The maxium distance from a dist to the box center
!================================================================================
      SUBROUTINE ComputeBoundingBox(Node)
         TYPE(PoleNode), POINTER     :: Node
         INTEGER                     :: I,J
         REAL(DOUBLE)                :: PMAx
!
         J=Qdex(Node%Bdex)
         Node%Box%BndBox(1,1) = Rho%Qx%D(J)-Ext(J)
         Node%Box%BndBox(1,2) = Rho%Qx%D(J)+Ext(J)
         Node%Box%BndBox(2,1) = Rho%Qy%D(J)-Ext(J)
         Node%Box%BndBox(2,2) = Rho%Qy%D(J)+Ext(J)
         Node%Box%BndBox(3,1) = Rho%Qz%D(J)-Ext(J)
         Node%Box%BndBox(3,2) = Rho%Qz%D(J)+Ext(J)
         DO I=Node%Bdex+1,Node%Edex
            J=Qdex(I)
            Node%Box%BndBox(1,1) = MIN(Node%Box%BndBox(1,1),Rho%Qx%D(J)-Ext(J))
            Node%Box%BndBox(1,2) = MAX(Node%Box%BndBox(1,1),Rho%Qx%D(J)+Ext(J))
            Node%Box%BndBox(2,1) = MIN(Node%Box%BndBox(2,1),Rho%Qy%D(J)-Ext(J))
            Node%Box%BndBox(2,2) = MAX(Node%Box%BndBox(2,1),Rho%Qy%D(J)+Ext(J))
            Node%Box%BndBox(3,1) = MIN(Node%Box%BndBox(3,1),Rho%Qz%D(J)-Ext(J))
            Node%Box%BndBox(3,2) = MAX(Node%Box%BndBox(3,1),Rho%Qz%D(J)+Ext(J))
         ENDDO
         Node%Box%Half(1)  = Half*(Node%Box%BndBox(1,2)-Node%Box%BndBox(1,1))
         Node%Box%Half(2)  = Half*(Node%Box%BndBox(2,2)-Node%Box%BndBox(2,1))
         Node%Box%Half(3)  = Half*(Node%Box%BndBox(3,2)-Node%Box%BndBox(3,1))
         Node%Box%Center(1)= Half*(Node%Box%BndBox(1,2)+Node%Box%BndBox(1,1))
         Node%Box%Center(2)= Half*(Node%Box%BndBox(2,2)+Node%Box%BndBox(2,1))
         Node%Box%Center(3)= Half*(Node%Box%BndBox(3,2)+Node%Box%BndBox(3,1))
!
         Node%DMax2 = Zero
         DO I=Node%Bdex,Node%Edex
            J=Qdex(I)
            PMax = (Rho%Qx%D(J)-Node%Box%Center(1))**2 &
                  +(Rho%Qy%D(J)-Node%Box%Center(2))**2 &
                  +(Rho%Qz%D(J)-Node%Box%Center(3))**2
            Node%DMax2=MAX(Node%DMax2,PMax)
         ENDDO
!
       END SUBROUTINE ComputeBoundingBox
!================================================================================
!     Bisection   
!================================================================================
      SUBROUTINE SplitPoleBox(Node,Left,Right)
         TYPE(PoleNode), POINTER     :: Node,Left,Right
         REAL(DOUBLE)                :: Section,MaxBox,Extent,MaxExt
         REAL(DOUBLE),DIMENSION(3)   :: MaxQL,MaxQH
         INTEGER                     :: B,E,N,ISplit,Split,I,J,k,IS
!--------------------------------------------------------------
!        Indexing
         B=Node%Bdex
         E=Node%Edex
         N=E-B+1
!        Determine the Split in the Largest Direction 
         ISplit = Mod(Node%Box%Tier,3)+1
!!$         MaxBox = Zero
!!$         DO I=1,3
!!$            IF(MaxBox<Node%Box%Half(I))THEN
!!$               MaxBox=Node%Box%Half(I)
!!$               ISplit=I
!!$            ENDIF
!!$         ENDDO
!!$!        Find the max extent, if > MaxBox*Two, split under extents
!!$         MaxExt=Zero
!!$         DO I=B,E
!!$            K=Qdex(I)
!!$            MaxExt=MAX(MaxExt,Ext(K))
!!$         ENDDO
!!$         IF(MaxExt>Two*MaxBox) ISplit=4
!        Split the Distibutuins
         J=0
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
         ELSEIF(ISplit==4)THEN
            DO I=B,E
               J=J+1;K=Qdex(I)
               RList(J)=Ext(K)
            ENDDO
         ENDIF
!        Sort
         CALL DblIntSort77(N,RList,Qdex(B:E),2)
!        Orthogonal RECURSIVE bisection (ORB)
         Left%Bdex =B
         Right%Edex=E
         Section   =RList(1)+Half*(RList(N)-RList(1))
         Split     =N/2
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
