!------------------------------------------------------------------------------
!--  This code is part of the MondoSCF suite of programs for linear scaling 
!    electronic structure theory and ab initio molecular dynamics.
!
!--  Copyright (c) 2001, the Regents of the University of California.  
!    This SOFTWARE has been authored by an employee or employees of the 
!    University of California, operator of the Los Alamos National Laboratory 
!    under Contract No. W-7405-ENG-36 with the U.S. Department of Energy.  
!    The U.S. Government has rights to use, reproduce, and distribute this 
!    SOFTWARE.  The public may copy, distribute, prepare derivative works 
!    and publicly display this SOFTWARE without charge, provided that this 
!    Notice and any statement of authorship are reproduced on all copies.  
!    Neither the Government nor the University makes any warranty, express 
!    or implied, or assumes any liability or responsibility for the use of 
!    this SOFTWARE.  If SOFTWARE is modified to produce derivative works, 
!    such modified SOFTWARE should be clearly marked, so as not to confuse 
!    it with the version available from LANL.  The return of derivative works
!    to the primary author for integration and general release is encouraged. 
!    The first publication realized with the use of MondoSCF shall be
!    considered a joint work.  Publication of the results will appear
!    under the joint authorship of the researchers nominated by their
!    respective institutions. In future publications of work performed
!    with MondoSCF, the use of the software shall be properly acknowledged,
!    e.g. in the form "These calculations have been performed using MondoSCF, 
!    a suite of programs for linear scaling electronic structure theory and
!    ab initio molecular dynamics", and given appropriate citation.  
!------------------------------------------------------------------------------
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
!
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
      SUBROUTINE RhoToPoleTree
        INTEGER                   :: Status,K,I
!-------------------------------------------------------------------
!       Initialize counters
        PoleNodes=0
        RhoLevel=0
!       Initialize the root node
        CALL NewPoleNode(PoleRoot,0)
        CALL InitPoleRoot
!       Convert the density into a 3-D BinTree
        CALL SplitPole(PoleRoot)
!       Make PoleTree tier by tier, recuring up from the bottom
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
!===============================================================================
      SUBROUTINE InitPoleRoot
         PoleRoot%Bdex=1
         PoleRoot%Edex=Rho%NDist
         PoleRoot%NQ=Rho%NDist
!        Get the Bounding Box for PoleRoot
         CALL NewRhoBox(PoleRoot)
!        Set the largest range (smallest exponent) of the roots BB
         PoleRoot%Zeta=Rho%Expt%D(1)
         MiniumExp=Rho%Expt%D(1)
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


!IF(Node%Zeta>100.D0)Node%Co=Zero

!WRITE(77,*)'test2['//TRIM(IntToChar(B))//']={qx->'//TRIM(DblToMMAChar(Node%Box%Center(1))) &
!                                         //',qy->'//TRIM(DblToMMAChar(Node%Box%Center(2))) &
!                                         //',qz->'//TRIM(DblToMMAChar(Node%Box%Center(3))) &
!                                         //',zq->'//TRIM(DblToMMAChar(Node%Zeta))          &
!                                         //',qco->'//TRIM(DblToMMAChar(Node%Co(1)))//'};'


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
!     ALLOCATE and read in the density, initalize global lists 
!========================================================================================
      SUBROUTINE InitRhoAux
         INTEGER      :: z,oq,or,iq,jq,NQ,Q,Ell,Status,I,IOS,LMNLen
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
         DEALLOCATE(Amp,STAT=Status)
         CALL DecMem(Status,0,Rho%NDist)
         DEALLOCATE(RList,STAT=Status) 
         CALL DecMem(Status,0,Rho%NDist)
       END SUBROUTINE DeleteRhoAux
!
END MODULE

