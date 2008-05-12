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
MODULE PoleTree
  USE Derivedtypes
  USE GlobalScalars   
  USE GlobalObjects
  USE PoleNodeType
  USE ProcessControl
  USE Indexing
  USE Parse
  USE InOut
  USE Macros
  USE QCTCThresholds
  USE QCTCIndexing
  USE BoundingBox
  USE MondoPoles
  USE Globals
  USE Density
  USE Order
  USE Clock
  USE MatFunk
  IMPLICIT NONE
  !  Globals
  TYPE(PoleNode), POINTER               :: PoleRoot ! Root of the pole tree 
  TYPE(PoleNode), POINTER               :: PR1
  INTEGER                               :: PoleNodes
  INTEGER                               :: RhoLevel
  INTEGER                               :: CurrentTier
  INTEGER                               :: MaxTier
  INTEGER                               :: NElecDist
  INTEGER,     DIMENSION(:),Allocatable :: Qdex,Ndex,Cdex,Ldex      
  REAL(DOUBLE),DIMENSION(:),Allocatable :: RList,Zeta,Ext      
  TYPE(INT_RNK2)                        :: Qndx,Cndx
  LOGICAL PPPRINT
  INTEGER                               :: MaxRhoEll,NGaussTotal
  !
  INTEGER, PARAMETER                    :: EllFit=1
  INTEGER, PARAMETER                    :: LenFit=(EllFit+1)*(EllFit+2)*(EllFit+3)/6    
  ! 
  INTEGER                               :: iHGStack
  REAL(DOUBLE),ALLOCATABLE,DIMENSION(:) :: HGStack
  !                               
  REAL(DOUBLE) :: PACFudgeFactor=1.0D0 ! If FF= 1/2, rigorous Cramer's bound
  REAL(DOUBLE),DIMENSION(0:12),PARAMETER :: Fact=(/1D0,1D0,2D0,6D0,24D0,120D0,      &
       720D0,5040D0,40320D0,362880D0,3628800D0,39916800D0,479001600D0/)
  !
!!$  REAL(DOUBLE)                          :: NodeWghts,MaxEx
!  REAL(DOUBLE), EXTERNAL                :: MondoTimer
  REAL(DOUBLE) :: Decompose_Time_Start,Decompose_Time, &
       TreeMake_Time_Start,TreeMake_Time
CONTAINS 
  !=================================================================================
  !     
  !=================================================================================
  SUBROUTINE RhoToPoleTree
    INTEGER                   :: Status,K,I,J
    !-------------------------------------------------------------------
    iHGStack=1
    ! 
    CALL New(Qndx,(/ClusterSize,MaxRhoEll/),(/1,0/))
    CALL New(Cndx,(/ClusterSize,MaxRhoEll/),(/1,0/))
    CALL New(CMTmp,(/ClusterSize,SPLen/),(/1,0/))
    CALL New(SMTmp,(/ClusterSize,SPLen/),(/1,0/))
    CALL SetDSYEVWork(LenFit)
    !
    PoleRoot%BdexE=1
    PoleRoot%EdexE=Rho%NDist
    PoleRoot%BdexN=1
    PoleRoot%EdexN=GM%NAtms
    !
    PoleNodes=2
    PoleRoot%Box%Tier=0
    PoleRoot%NQ=Rho%NDist
    PoleRoot%Box%Number=1
    PoleRoot%NAtms=GM%NAtms
    !
    PoleRoot%Leaf=.False.

    NULLIFY(PoleRoot%Next)
    ! Recursively split the density into a k-d tree
    Decompose_Time_Start=MTimer()
    NGaussTotal=0;
    CALL SplitPole(PoleRoot,PoleRoot%Next)
    WRITE(*,*)' NGaussTotal = ',NGaussTotal
    
    TreeMake_Time_Start=MTimer()
    Decompose_Time=TreeMake_Time_Start-Decompose_Time_Start
    ! Make PoleTree tier by tier, recuring up from the bottom
    CALL MakePoleTree(PoleRoot) 
    !  CALL Print_POLEROOT_PAC(PoleRoot)
    !  STOP
    TreeMake_Time=MTimer()-TreeMake_Time_Start
    !
    ! Delete global arrays etc
    CALL Delete(Qndx)
    CALL Delete(Cndx)
    CALL Delete(CMTmp)
    CALL Delete(SMTmp)
    CALL UnSetDSYEVWork
    !
    ALLOCATE(HGStack(iHGStack),STAT=Status)
    CALL IncMem(Status,0,iHGStack)
    J=iHGStack
    iHGStack=1
    CALL PackStack(PoleRoot)    

  END SUBROUTINE RhoToPoleTree
  !=====================================================================================
  !     
  !=====================================================================================
  RECURSIVE SUBROUTINE SplitPole(Node,Next)
    TYPE(PoleNode), POINTER :: Node,Next
    !--------------------------------------------------------------
! IF(Node%NAtms<10)THEN
    IF(Node%NAtms==1)THEN
       ! We are at a leaf
       NULLIFY(Node%Left)
       NULLIFY(Node%Right)
       NULLIFY(Node%Next)
       CALL FillRhoLeaf(Node)
       ! Set link
       Node%Next=>Next
    ELSE 
       CALL NewPoleNode(Node%Left,Node%Box%Tier+1)
       CALL NewPoleNode(Node%Right,Node%Box%Tier+1)
       ! Split the Box
       CALL SplitPoleBox(Node,Node%Left,Node%Right)
       ! Recur
       Node%Next=>Next
       CALL SplitPole(Node%Left,Node%Right)
       CALL SplitPole(Node%Right,Next)
    ENDIF
  END SUBROUTINE SplitPole
  !
  RECURSIVE SUBROUTINE MakePoleTree(P)
    TYPE(PoleNode),POINTER  :: P
    REAL(DOUBLE)    :: PMAX
    INTEGER         :: I,J,Status
    !--------------------------------------------------------------
    IF(.NOT.P%Leaf)THEN
       CALL MakePoleTree(P%Left)
       CALL MakePoleTree(P%Right)
       ! Walking backwards up the tree, merge left and right nodes:
       ! Merge of Hermite Gaussian exponents
       P%Herm%Ell=MAX(P%Left%Herm%Ell,P%Right%Herm%Ell)
       ! Merge the integral estimates ...
       P%IHalf=P%Left%IHalf+P%Right%IHalf
       ! ... the boxes ...
       CALL BoxMerge(P%Left%Box,P%Right%Box,P%Box)
       P%Pole%Center=P%Box%Center  
       ! ... and merge the multipole expansions ...
       CALL PoleMerge(P%Left%Pole,P%Right%Pole,P%Pole)
       CALL SetSignedPAC(P%PAC,P)
       !       ALLOCATE(P%MAC%O(0:P%HERM%Ell),STAT=Status)
       !       CALL IncMem(Status,0,P%HERM%Ell+1)       
	CALL SetMAC(P)
    ENDIF
  END SUBROUTINE MakePoleTree

  SUBROUTINE SplitPoleBox(Node,Left,Right)
    TYPE(PoleNode), POINTER     :: Node,Left,Right
    REAL(DOUBLE)                :: Section,MaxBox,Extent,MaxExt,IHalfHalf,BalanceQ,TotCh2,ChHalf2
    REAL(DOUBLE),DIMENSION(3)   :: MaxQL,MaxQH
    INTEGER                     :: Ne,Nn,Be,Ee,Bn,En,Je,Jn,ISplit,Split,I,J,k,IS,NewSplit
    INTEGER                     :: MaxBoxedEll,NBoxedEll,EllSplit,Cd,SplitI,SplitJ
    REAL(DOUBLE)                :: Volume,MaxZeta,MinZeta,ZetaHalf,ChHalf,TotCh,PiZ
    REAL(DOUBLE)                :: BEBalance,ChBalance,XBalance
    INTEGER,DIMENSION(:), ALLOCATABLE :: IList,JList
    REAL(DOUBLE),DIMENSION(:), ALLOCATABLE :: XList,CList
    LOGICAL :: Tru
    !--------------------------------------------------------------
    ! Indexing to start

!    WRITE(*,*)'===================================================================='
!    WRITE(*,*)Node%Box%Tier,Node%Box%Number

    ! Begin and end index for charges (electrons and nuclei)
    Be=Node%BdexE
    Ee=Node%EdexE
    ! Begin and end index for nuclei only
    Bn=Node%BdexN
    En=Node%EdexN
    ! Number of each in this box
    Ne=Ee-Be+1
    Nn=En-Bn+1
    !
    IF(Nn<1)CALL Halt(' Logic error in SplitPoleBox ')
    ! Determine the longest coordinate ...
    MaxBox=Zero
    DO I=1,3
       IF(MaxBox<Node%Box%Half(I))THEN
          MaxBox=Node%Box%Half(I)
          ISplit=I
       ENDIF
    ENDDO
    !
    IF(ISplit==1)THEN
       CALL FairSplit(Ne,Qdex(Be:Ee),Rho%Qx%D,Je,Nn,ISplit,NDex(Bn:En),GM%Carts%D,Jn)
    ELSEIF(ISplit==2)THEN
       CALL FairSplit(Ne,Qdex(Be:Ee),Rho%Qy%D,Je,Nn,ISplit,NDex(Bn:En),GM%Carts%D,Jn)
    ELSEIF(ISplit==3)THEN
       CALL FairSplit(Ne,Qdex(Be:Ee),Rho%Qz%D,Je,Nn,ISplit,NDex(Bn:En),GM%Carts%D,Jn)
    ENDIF
    ! Split the charges
    Split=Be+Je-1
    Left%BdexE=Be
    Left%EdexE=Split
    Right%EdexE=Ee
    Right%BdexE=Split+1
    IF(Left%EdexE<Left%BdexE.OR.Right%EdexE<Right%BdexE)THEN
       CALL Halt(' Charge indexing hosed in SplitPoleBox ' )
    ENDIF
    ! ... L&R counters ...
    Left%NQ=Left%EdexE-Left%BdexE+1
    Right%NQ=Right%EdexE-Right%BdexE+1
    ! Split the nuclei
    Split=Bn+Jn-1
    Left%BdexN=Bn
    Left%EdexN=Split
    Right%EdexN=En 
    Right%BdexN=Split+1
    IF(Left%EdexN<Left%BdexN.OR.Right%EdexN<Right%BdexN)THEN
       CALL Halt(' Nuclei indexing hosed in SplitPoleBox ' )
    ENDIF
    ! ... L&R counters ...
    Left%NAtms=Left%EdexN-Left%BdexN+1
    Right%NAtms=Right%EdexN-Right%BdexN+1
    ! Ok, I hate this part but there seems to be no escaping 
    ! it for now.  Grit our teeth and do a top down sintering of the BBoxes
    Left%Box%BndBox(:,1)=1D20
    Left%Box%BndBox(:,2)=-1D20
    Left%IHMin=1D20
    Left%IHalf=-1D20
    DO I=Left%BdexE,Left%EdexE
       K=Qdex(I)
       Left%Box%BndBox(1,1)=MIN(Left%Box%BndBox(1,1),Rho%Qx%D(K))
       Left%Box%BndBox(1,2)=MAX(Left%Box%BndBox(1,2),Rho%Qx%D(K))
       Left%Box%BndBox(2,1)=MIN(Left%Box%BndBox(2,1),Rho%Qy%D(K))
       Left%Box%BndBox(2,2)=MAX(Left%Box%BndBox(2,2),Rho%Qy%D(K))
       Left%Box%BndBox(3,1)=MIN(Left%Box%BndBox(3,1),Rho%Qz%D(K))
       Left%Box%BndBox(3,2)=MAX(Left%Box%BndBox(3,2),Rho%Qz%D(K))
       Left%IHMin=MIN(Left%IHMin,Ext(K))
       Left%IHalf=MAX(Left%IHalf,Ext(K))
    ENDDO
    Left%Box%Half(:)=Half*(Left%Box%BndBox(:,2)-Left%Box%BndBox(:,1))
    Left%Box%Center(:)=Left%Box%BndBox(:,1)+Left%Box%Half(:)
    !
    Right%Box%BndBox(:,1)=1D20
    Right%Box%BndBox(:,2)=-1D20
    Right%IHMin=1D20
    Right%IHalf=-1D20
    DO I=Right%BdexE,Right%EdexE
       K=Qdex(I)
       Right%Box%BndBox(1,1)=MIN(Right%Box%BndBox(1,1),Rho%Qx%D(K))
       Right%Box%BndBox(1,2)=MAX(Right%Box%BndBox(1,2),Rho%Qx%D(K))
       Right%Box%BndBox(2,1)=MIN(Right%Box%BndBox(2,1),Rho%Qy%D(K))
       Right%Box%BndBox(2,2)=MAX(Right%Box%BndBox(2,2),Rho%Qy%D(K))
       Right%Box%BndBox(3,1)=MIN(Right%Box%BndBox(3,1),Rho%Qz%D(K))
       Right%Box%BndBox(3,2)=MAX(Right%Box%BndBox(3,2),Rho%Qz%D(K))
       Right%IHMin=MIN(Right%IHMin,Ext(K))
       Right%IHalf=MAX(Right%IHalf,Ext(K))
    ENDDO
    Right%Box%Half(:)=Half*(Right%Box%BndBox(:,2)-Right%Box%BndBox(:,1))
    Right%Box%Center(:)=Right%Box%BndBox(:,1)+Right%Box%Half(:)


!    WRITE(*,55)ISplit,Left%Box%BndBox(ISplit,1),Left%Box%BndBox(ISplit,2), &
!                  Right%Box%BndBox(ISplit,1),Right%Box%BndBox(ISplit,2)
!
!55  FORMAT("  Split = ",I2," [",D12.6,", ",D12.6,"] [",D12.6,", ",D12.6," ] ")

!!$    WRITE(*,55)ISplit,Node%Box%Number,Left%BdexE,Left%EDexE,Left%EdexE-Left%BDexE, &
!!$                                      Right%BdexE,Right%EDexE,Right%EdexE-Right%BDexE, &
!!$                                      Left%BdexN,Left%EDexN,Right%BdexN,Right%EDexN
!!$
!!$55  FORMAT("  Split = ",I2,", Node = ",I4," [",I6,", ",I6,";",I6,"] [",I6,", ",I6,";",I6," ] // <",I2,", ",I2,"><",I2,", ",I2,"> ")


  END SUBROUTINE SplitPoleBox
  !
  SUBROUTINE FairSplit(Ne,Qe,RhoX,Je,Nn,ISplit,Qn,NucXYZ,Jn)
    INTEGER :: Ne,Nn,Je,Jn,J,ISplit
    REAL(DOUBLE) :: Section
    INTEGER,DIMENSION(1:Ne) :: Qe
    INTEGER,DIMENSION(1:Nn) :: Qn
    REAL(DOUBLE),DIMENSION(:)   :: RhoX
    REAL(DOUBLE),DIMENSION(:,:) :: NucXYZ
    REAL(DOUBLE),DIMENSION(1:Ne) :: X
    !
    DO J=1,Nn
       X(J)=NucXYZ(ISplit,Qn(J))
    ENDDO
    CALL DblIntSort77(Nn,X,Qn,2)             
    !
    Jn=Nn/2
    Section=Half*(X(Jn)+X(Jn+1))
    !
    DO J=1,Ne
       X(J)=RhoX(Qe(J))
    ENDDO
    CALL DblIntSort77(Ne,X,Qe,2)             
    !
    DO J=1,Ne
       IF(X(J)>Section)THEN
          Je=J-1
          EXIT
       ENDIF
    ENDDO
    !
  END SUBROUTINE FairSplit
  !
  RECURSIVE SUBROUTINE PackStack(Q)
    TYPE(PoleNode),POINTER  :: Q
    INTEGER :: I,J,Nq,EllQ,LenQ
    !--------------------------------------------------------------
    IF(.NOT.Q%Leaf)THEN
       CALL PackStack(Q%Left)
       CALL PackStack(Q%Right)
    ELSE
       Q%Herm%Stack=iHGStack
!       WRITE(*,*)'====================================================='
!       WRITE(*,*)'   Node Number = ',Q%BOX%Number
!       WRITE(*,*)'   Stack     # = ',iHGStack,Q%Herm%Stack
       HGStack(iHGStack)=Q%HERM%Ell
!       WRITE(*,*)iHGStack,' Ell = ',Q%HERM%Ell
       iHGStack=iHGStack+1       

       DO EllQ=0,Q%HERM%Ell
          LenQ=LHGTF(EllQ)
          Nq=Q%HERM%Nq(EllQ)
          HGStack(iHGStack)=Nq
          iHGStack=iHGStack+1
          DO I=1,Nq
             HGStack(iHGStack)=Q%HERM%IHlf(EllQ)%D(I)
             iHGStack=iHGStack+1
          ENDDO
          DO I=1,Nq
             HGStack(iHGStack)=Q%HERM%Zeta(EllQ)%D(I)
             iHGStack=iHGStack+1
          ENDDO
          DO I=1,Nq
             HGStack(iHGStack)=Q%HERM%Cent(EllQ)%D(1,I)
             iHGStack=iHGStack+1
             HGStack(iHGStack)=Q%HERM%Cent(EllQ)%D(2,I)
             iHGStack=iHGStack+1
             HGStack(iHGStack)=Q%HERM%Cent(EllQ)%D(3,I)
             iHGStack=iHGStack+1
          ENDDO
          DO I=1,Nq
             DO J=1,LenQ
                HGStack(iHGStack)=Q%HERM%Coef(EllQ)%D(J,I)
                iHGStack=iHGStack+1
             ENDDO
          ENDDO
       ENDDO
       DO EllQ=0,Q%HERM%Ell
	  IF(Q%HERM%NQ(EllQ).NE.0)THEN
!!$          CALL Delete(Q%HERM%IHlf(EllQ))
!!$	     CALL Delete(Q%HERM%Zeta(EllQ))
!!$	     CALL Delete(Q%HERM%Cent(EllQ))
!!$	     CALL Delete(Q%HERM%Coef(EllQ))
          ENDIF
       ENDDO	
!!$       DEALLOCATE(Q%HERM%NQ)
!!$       DEALLOCATE(Q%HERM%Coef)
!!$       DEALLOCATE(Q%HERM%Cent)
!!$       DEALLOCATE(Q%HERM%Zeta)
!!$       DEALLOCATE(Q%HERM%IHlf)
    ENDIF
  END SUBROUTINE PackStack
  !=====================================================================================
  !
  !=====================================================================================
  SUBROUTINE FillRhoLeaf(Node)
    TYPE(POLENode), POINTER :: Node
    INTEGER                 :: I,IQ,IC,J,K,KQ,KC,L,B,E,N,NQ,NC,LMNLen,Status,Ell,Qd,Cd        
    REAL(DOUBLE)            :: RhoSum,ZE,EX,delta
    !-------------------------------------------------------------------------------------
    ! Set leaf logical
    Node%Leaf=.True.
    !
    ALLOCATE(Node%HERM%NQ(0:MaxRhoEll))
    ALLOCATE(Node%HERM%Coef(0:MaxRhoEll))
    ALLOCATE(Node%HERM%Cent(0:MaxRhoEll))
    ALLOCATE(Node%HERM%Zeta(0:MaxRhoEll))
    ALLOCATE(Node%HERM%IHlf(0:MaxRhoEll))
    Node%HERM%NQ=Zero
    Node%IHalf=Zero
    !
    Cndx%I=Zero
    Qndx%I=Zero
    Node%HERM%Ell=0
    Node%BOX%BndBOX(:,1)=1D20
    Node%BOX%BndBOX(:,2)=-1D20
    ! 
    ! Order by integral magnitude
    J=0
    B=Node%BdexE
    E=Node%EdexE
    N=E-B+1
    DO I=B,E
       J=J+1
       K=Qdex(I)
       RList(J)=Ext(K)
    ENDDO
    CALL DblIntSort77(N,RList,Qdex(B:E),-2)
    !
    ! This for MaxEll
    iHGStack=iHGStack+1
    !
    delta=0D0
    Node%POLE%Charge=GM%AtNum%D(NDex(Node%BdexN))
    Node%POLE%Center=GM%Carts%D(:,NDex(Node%BdexN))

    DO I=Node%BdexE,Node%EdexE

       Qd=Qdex(I)
       !
       Ell=LDex(Qd)
       Node%HERM%Ell=MAX(Node%HERM%Ell,Ell)
       Node%HERM%NQ(Ell)=Node%HERM%NQ(Ell)+1
       N=Node%HERM%NQ(Ell)       
       !
       Qndx%I(N,Ell)=Qd
       Cndx%I(N,Ell)=Cdex(Qd)
       !
       ZE=Zeta(Qd)
       Cd=Cdex(Qd)
       LMNLen=LHGTF(Ell)
       ! Set and this nodes BBOX
       Node%BOX%BndBOX(1,1)=MIN(Node%BOX%BndBOX(1,1),Rho%Qx%D(Qd))
       Node%BOX%BndBOX(1,2)=MAX(Node%BOX%BndBOX(1,2),Rho%Qx%D(Qd))
       Node%BOX%BndBOX(2,1)=MIN(Node%BOX%BndBOX(2,1),Rho%Qy%D(Qd))
       Node%BOX%BndBOX(2,2)=MAX(Node%BOX%BndBOX(2,2),Rho%Qy%D(Qd))
       Node%BOX%BndBOX(3,1)=MIN(Node%BOX%BndBOX(3,1),Rho%Qz%D(Qd))
       Node%BOX%BndBOX(3,2)=MAX(Node%BOX%BndBOX(3,2),Rho%Qz%D(Qd))


       delta=MAX(delta,SQRT(DOT_PRODUCT( (/Rho%Qx%D(Qd)-NODE%Pole%Center(1),Rho%Qy%D(Qd)-NODE%Pole%Center(2),Rho%Qz%D(Qd)-NODE%Pole%Center(3)/) , &
                                         (/Rho%Qx%D(Qd)-NODE%Pole%Center(1),Rho%Qy%D(Qd)-NODE%Pole%Center(2),Rho%Qz%D(Qd)-NODE%Pole%Center(3)/) )))

!!$!        Set and inflate this nodes BBox
!!$         Node%Box%BndBox(1,:)=Rho%Qx%D(KQ)
!!$         Node%Box%BndBox(2,:)=Rho%Qy%D(KQ)
!!$         Node%Box%BndBox(3,:)=Rho%Qz%D(KQ)
!!$         Node%Box=ExpandBox(Node%Box,Ext(KQ))


    ENDDO
    Node%POLE%delta=delta

    ! The leaf multipole center should always be at the atomic center
    Node%BOX%Center(:)=Node%BOX%BndBOX(:,1)+Node%Box%Half    
    Node%BOX%Half(:)=Half*(Node%BOX%BndBOX(:,2)-Node%BOX%BndBOX(:,1))

    DO Ell=0,Node%HERM%Ell
       N=Node%HERM%NQ(Ell)
       iHGStack=iHGStack+1
       IF(N.NE.0)THEN
          LMNLen=LHGTF(Ell)
          CALL New(Node%HERM%Coef(Ell),(/LMNLen,N/))
          CALL New(Node%HERM%Cent(Ell),(/3,N/))
          CALL New(Node%HERM%Zeta(Ell),N)          
          CALL New(Node%HERM%IHlf(Ell),N)
          iHGStack=iHGStack+LMNLen*N+5*N
          DO I=1,N
             Qd=Qndx%I(I,Ell)
             Cd=Cndx%I(I,Ell)
             Node%HERM%Coef(Ell)%D(1:LMNLen,I)=Rho%Co%D(Cd:Cd+LMNLen-1)
             Node%HERM%Cent(Ell)%D(:,I)=(/Rho%Qx%D(Qd),Rho%Qy%D(Qd),Rho%Qz%D(Qd)/)
             Node%HERM%Zeta(Ell)%D(I)=Zeta(Qd)
             Node%HERM%IHlf(Ell)%D(I)=EXT(Qd)
             Node%IHalf=Node%IHalf+EXT(Qd)
          ENDDO
       ENDIF
    ENDDO
    ! Fill in the PAC and the MAC
    CALL SetSignedPAC(Node%PAC,Node)
    ! Fill in the multiPOLE 
    Node%POLE%Ell=MaxPoleEll
    CALL AllocSP(Node%POLE)
    ! Double check we have a valid expansion center (ie in the box!)
    IF(PointOutSideBox(Node%POLE%Center,Node%BOX))THEN
       WRITE(*,32)Node%POLE%Center
32     FORMAT(' POLE%Center = ',3(D12.6,', '))
       CALL PrintBBox(Node%BOX)
       CALL Halt(' In FillRhoLeaf: Multipole center outside of BBox ')
    ENDIF
    !
    CALL HGToSP_POLENODE(Node%HERM,Node%POLE)
    CALL SetMAC(Node)

    NGaussTotal=NGaussTotal+(Node%EdexE-Node%BdexE)
    !
    WRITE(*,33)Node%BOX%Number,Node%POLE%Charge,Node%EdexE-Node%BdexE,NODE%POLE%Delta
33  FORMAT(' Node = ',I4,' Z = ',F4.1,' NGauss = ',I4,' Delta = ',D8.2)
    !
  END SUBROUTINE FillRhoLeaf
  !=================================================================================
  !     
  !=================================================================================     
  SUBROUTINE PoleMerge(LPole,RPole,PPole)
    TYPE(Pole)                        :: LPole,RPole,PPole       
    !------------------------------------------------------------------------------------       
    PPole%Ell=MAX(LPole%Ell,RPole%Ell)
    PPole%Center=(LPole%Charge*LPole%Center+RPole%Charge*RPole%Center)/(LPole%Charge+RPole%Charge)
     CALL AllocSP(PPole)
    ! Move the left center to the new midpoint
    CALL XLate(LPole,PPole)
    ! Move the right center to the new midpoint
    CALL XLate(RPole,PPole)
  END SUBROUTINE PoleMerge
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
    MaxTier=MAX(MaxTier,Level)
    PoleNodes=PoleNodes+1
    NULLIFY(Node%Left)
    NULLIFY(Node%Right)
  END SUBROUTINE NewPoleNode
  !==========================================================================
  !     Initialize a new PoleNodes Array
  !==========================================================================
  SUBROUTINE AllocSP(Node)
    TYPE(Pole)  :: Node
    INTEGER                   :: LenSP,Status        
    LenSP=LSP(Node%Ell)
    ALLOCATE(Node%S(0:LenSP),STAT=Status)
    CALL IncMem(Status,0,LenSP+1)
    ALLOCATE(Node%C(0:LenSP),STAT=Status)
    CALL IncMem(Status,0,LenSP+1)
    CALL DBL_VECT_EQ_DBL_SCLR(LenSP+1,Node%S(0),Zero)
    CALL DBL_VECT_EQ_DBL_SCLR(LenSP+1,Node%C(0),Zero)
  END SUBROUTINE AllocSP

  !========================================================================================
  !     ALLOCATE and read in the density, initalize global lists 
  !========================================================================================
  SUBROUTINE InitRhoAux
    INTEGER      :: z,oq,or,iq,jq,NQ,Q,Ell,Status,I,IOS,LMNLen,CD,QD,K
    REAL(DOUBLE) :: ZE,EX,CheckChg
    TYPE(DBL_VECT) :: Est
    !-------------------------------------------------------------
    !  ALLOCATE global lists
    ALLOCATE(Qdex(1:Rho%NDist),STAT=Status)
    CALL IncMem(Status,Rho%NDist,0)
    ALLOCATE(Ndex(1:GM%NAtms),STAT=Status)
    CALL IncMem(Status,GM%NAtms,0)
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
    !
    PoleRoot%Box%BndBox(:,1)=1D20
    PoleRoot%Box%BndBox(:,2)=-1D20
    PoleRoot%IHMin=1D20
    PoleRoot%IHalf=-1D20
    !
    I=0
    IQ=1
    MaxRhoEll=0
    DO z=1,Rho%NExpt
       oq =Rho%OffQ%I(z)   
       or =Rho%OffR%I(z)   
       Ell=Rho%Lndx%I(z)   
       MaxRhoEll=MAX(MaxRhoEll,Ell)
       ZE =Rho%Expt%D(z)
       LMNLen=LHGTF(Ell)
       !
       CALL New(Est,Rho%NQ%I(z))
       !
       CALL IBounds(Ell,LMNLen,Rho%NQ%I(z),ZE,Rho%Co%D(or+1),  &
            PLMNx(Ell,Ell)%I(1),QLMNx(Ell,Ell)%I(1),   &
            PQLMNx(Ell,Ell)%I(1),Est%D(1))
       !
       DO Q=1,Rho%NQ%I(z)          
          QD=oq+Q
          CD=or+(Q-1)*LMNLen+1


!!$          EX=Extent(Ell,ZE,Rho%Co%D(CD:CD+LMNLen-1),Tau_O=TauPAC,Potential_O=.TRUE.,ExtraEll_O=0)
!!$          IF(EX>Zero)THEN


          EX=Est%D(Q)  
          IF(EX>TauTwo*1D-5)THEN
             ! Intialize PoleRoot boundaries
             PoleRoot%Box%BndBox(1,1)=MIN(PoleRoot%Box%BndBox(1,1),Rho%Qx%D(IQ))
             PoleRoot%Box%BndBox(1,2)=MAX(PoleRoot%Box%BndBox(1,1),Rho%Qx%D(IQ))
             PoleRoot%Box%BndBox(2,1)=MIN(PoleRoot%Box%BndBox(2,1),Rho%Qy%D(IQ))
             PoleRoot%Box%BndBox(2,2)=MAX(PoleRoot%Box%BndBox(2,1),Rho%Qy%D(IQ))
             PoleRoot%Box%BndBox(3,1)=MIN(PoleRoot%Box%BndBox(3,1),Rho%Qz%D(IQ))
             PoleRoot%Box%BndBox(3,2)=MAX(PoleRoot%Box%BndBox(3,1),Rho%Qz%D(IQ))
             PoleRoot%IHalf=MAX(PoleRoot%IHalf,EX)
             PoleRoot%IHMin=MIN(PoleRoot%IHMin,EX)
             ! Here are the rest of the indexing variables
             Qdex(IQ)=QD
             Cdex(QD)=CD
             Ext( QD)=EX
             Zeta(QD)=ZE
             Ldex(QD)=Ell
             IQ=IQ+1
          ENDIF
       ENDDO
       CALL Delete(Est)
    ENDDO
    !
    Rho%NDist=IQ-1
    ! Number of distributions excluding nuclei (used for
    ! identifiying nuclear selfinteraction, note that nuclei
    ! must be in correct order for this to work...)
    NElecDist=SUM(Rho%NQ%I(1:Rho%NExpt-1))
    ! Here is indexing of the nuclear centers only
    DO I=1,GM%NAtms
       NDex(I)=I
    ENDDO
    !
    PoleRoot%Box=ExpandBox(PoleRoot%Box,1D-20)
  END SUBROUTINE InitRhoAux
  !========================================================================================
  ! Delete globals associated with the density indexing
  !========================================================================================
  SUBROUTINE DeleteRhoAux
    INTEGER :: Status
    DEALLOCATE(Qdex,STAT=Status)
    CALL DecMem(Status,Rho%NDist,0)
    DEALLOCATE(Cdex,STAT=Status)
    CALL DecMem(Status,Rho%NDist,0)
    DEALLOCATE(Ndex,STAT=Status)
    CALL DecMem(Status,GM%NAtms,0)
    DEALLOCATE(Ldex,STAT=Status)
    CALL DecMem(Status,Rho%NDist,0)
    DEALLOCATE(Zeta,STAT=Status)
    CALL DecMem(Status,0,Rho%NDist)
    DEALLOCATE(Ext,STAT=Status)
    CALL DecMem(Status,0,Rho%NDist)
    DEALLOCATE(RList,STAT=Status) 
    CALL DecMem(Status,0,Rho%NDist)
  END SUBROUTINE DeleteRhoAux
  !=================================================================================
  !     
  !=================================================================================
  SUBROUTINE SetSignedPAC(PC,Node)
    TYPE(PoleNode),POINTER    :: Node,Q
    TYPE(PAC)                 :: PC
    INTEGER                          :: Ell,L,M,N,LMN,I,J,ICrnr,LP,MP,NP,LQ,MQ,NQ,Pdex,Qdex
    REAL(DOUBLE)                     :: RTE,RPE,T,Omega,Upq
    REAL(DOUBLE)                     :: NodeWeightC,Zeta,Zeta2,PiZ,MinMax,ZetaFac,PCFit
    REAL(DOUBLE)                     :: CramCo,MixMax,ScaledTau,HGInEq,ZZ,CO,PZ,MaxExtent,Ch,WMin,WMax
    REAL(DOUBLE),DIMENSION(0:EllFit)  :: Wght
    REAL(DOUBLE),DIMENSION(3,14)      :: Centers
    REAL(DOUBLE),DIMENSION(LenFit,14) :: ChProj
    REAL(DOUBLE),DIMENSION(LenFit)    :: FitMax,Fit
    REAL(DOUBLE),DIMENSION(NAtoms)    :: NukeMask
    REAL(DOUBLE), DIMENSION(0:2*EllFit) :: AuxR
    REAL(DOUBLE), DIMENSION(LenFit,LenFit) :: FitTrix,InvTrix
    REAL(DOUBLE),DIMENSION(0:2*EllFit,0:2*EllFit,0:2*EllFit,0:2*EllFit) :: MDR
    LOGICAL PrintFlag
    !-------------------------------------------------------------------------------    
    ! H_n(t) < K 2^(n/2) SQRT(n!) EXP(t^2/2), with K=1.09
    !-------------------------------------------------------------------------------    
    PC%Zeta=PACZetaMin(Node)

!    IF(PC%Zeta==NuclearExpnt)PC%Zeta=PC%Zeta-1D10


    ! Corners
    Centers(:,1)=(/Node%Box%BndBox(1,1),Node%Box%BndBox(2,1),Node%Box%BndBox(3,1)/)
    Centers(:,2)=(/Node%Box%BndBox(1,1),Node%Box%BndBox(2,1),Node%Box%BndBox(3,2)/)
    Centers(:,3)=(/Node%Box%BndBox(1,1),Node%Box%BndBox(2,2),Node%Box%BndBox(3,1)/)
    Centers(:,4)=(/Node%Box%BndBox(1,2),Node%Box%BndBox(2,1),Node%Box%BndBox(3,1)/)
    Centers(:,5)=(/Node%Box%BndBox(1,2),Node%Box%BndBox(2,2),Node%Box%BndBox(3,1)/)
    Centers(:,6)=(/Node%Box%BndBox(1,2),Node%Box%BndBox(2,1),Node%Box%BndBox(3,2)/)
    Centers(:,7)=(/Node%Box%BndBox(1,1),Node%Box%BndBox(2,2),Node%Box%BndBox(3,2)/)
    Centers(:,8)=(/Node%Box%BndBox(1,2),Node%Box%BndBox(2,2),Node%Box%BndBox(3,2)/) 
    ! Sides
    Centers(:,9)=(/Node%Box%Center(1),Node%Box%BndBox(2,1),Node%Box%BndBox(3,1)/)
    Centers(:,10)=(/Node%Box%BndBox(1,1),Node%Box%Center(2),Node%Box%BndBox(3,1)/)
    Centers(:,11)=(/Node%Box%BndBox(1,1),Node%Box%BndBox(2,1),Node%Box%Center(3)/)
    Centers(:,12)=(/Node%Box%Center(1),Node%Box%BndBox(2,2),Node%Box%BndBox(3,2)/)
    Centers(:,13)=(/Node%Box%BndBox(1,2),Node%Box%Center(2),Node%Box%BndBox(3,2)/)
    Centers(:,14)=(/Node%Box%BndBox(1,2),Node%Box%BndBox(2,2),Node%Box%Center(3)/)

    ChProj=PACChargeProject(Node,PC%Zeta,Centers)
!    WRITE(*,*)'  ChProj = ',ChProj

    RTE=PC%Zeta*PC%Zeta
    RPE=PC%Zeta+PC%Zeta
    T=0D0
    Omega=RTE/RPE 
#ifdef OVERLAP_PROJECTION    
    Upq=(Pi/RPE)**1.5D0
    CALL OvrInts(2*EllFit,2*EllFit,AuxR,Omega,T)
#else
    Upq=TwoPi5x2/(RTE*SQRT(RPE))
    CALL AuxInts(2*EllFit,2*EllFit,AuxR,Omega,T)
#endif
    CALL MD3TRR(2*EllFit,2*EllFit,MDR,AuxR,Upq,0D0,0D0,0D0)
    DO LP=0,EllFit
       DO MP=0,EllFit-LP
          DO NP=0,EllFit-LP-MP 
             PDex=LMNx(LP,MP,NP)                   
             DO LQ=0,EllFit
                DO MQ=0,EllFit-LQ
                   DO NQ=0,EllFit-LQ-MQ
                      QDex=LMNx(LQ,MQ,NQ)
                      FitTrix(PDex,QDex)=MDR(LP+LQ,MP+MQ,NP+NQ,0)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    CALL FunkOnSqMat(LenFit,Inverse,FitTrix,InvTrix,PosDefMat_O=.FALSE.,EigenThresh_O=0D0)

    PC%Wght=0D0
    DO J=1,14
       Fit=MATMUL(InvTrix,ChProj(:,J))
       Wght=Zero
       DO L=0,EllFit
          DO M=0,EllFit-L
             DO N=0,EllFit-L-M
                LMN=LMNx(L,M,N)
                MixMax=Fact(L)*Fact(M)*Fact(N)
                HGInEq=SQRT(MixMax*(Two*PC%Zeta)**(L+M+N))*Fit(LMN)
                Wght(L+M+N)=MAX(Wght(L+M+N),ABS(HGInEq))
             ENDDO
          ENDDO
       ENDDO
       PC%Wght=MAX(PC%Wght,SUM(Wght(0:EllFit)))
    ENDDO
    PiZ=(Pi/PC%Zeta)**(ThreeHalves)
    PC%Wght=PC%Wght*PiZ
    IF(EllFit.NE.0)PC%Zeta=PC%Zeta*PACFudgeFactor

!!$    IF(Node%Leaf.and.PC%Wght<1D-40)THEN
!!$       WRITE(*,*)' Leaf? = ',Node%Leaf,' Wght = ',PC%Wght
!!$       WRITE(*,*)' Node center = ',Node%HERM%Cent(0)%D(:,1)
!!$       WRITE(*,*)' ChProj  = ',ChProj
!!$       WRITE(*,*)' FitTrix = ',FitTrix
!!$       DO J=1,14
!!$          Fit=MATMUL(InvTrix,ChProj(:,J))
!!$          WRITE(*,*)' InVTrix = ',InvTrix
!!$          WRITE(*,*)' Fit =',Fit
!!$          Wght=Zero
!!$          DO L=0,EllFit
!!$             DO M=0,EllFit-L
!!$                DO N=0,EllFit-L-M
!!$                   LMN=LMNx(L,M,N)
!!$                   MixMax=Fact(L)*Fact(M)*Fact(N)
!!$                   HGInEq=SQRT(MixMax*(Two*PC%Zeta)**(L+M+N))*Fit(LMN)
!!$                   Wght(L+M+N)=MAX(Wght(L+M+N),ABS(HGInEq))
!!$                ENDDO
!!$             ENDDO
!!$          ENDDO
!!$          PC%Wght=MAX(PC%Wght,SUM(Wght(0:EllFit)))
!!$       ENDDO
!!$       
!!$
!!$
!!$!       DO I=1,14
!!$!          WRITE(*,*)' Centers     = ',Centers(:,I)
!!$!       ENDDO
!!$    ENDIF
  END SUBROUTINE SetSignedPAC

  RECURSIVE FUNCTION PACChargeProject(Q,PCZeta,Corners) RESULT(ChFit)
    TYPE(PoleNode)                 :: Q
    REAL(DOUBLE)                   :: PCZeta,PiZ,Ch
    REAL(DOUBLE),DIMENSION(3,14)    :: Corners
    REAL(DOUBLE),DIMENSION(LenFit,14)      :: ChFit
    INTEGER                        :: I,EllQ,Nuke,LenP,PQEll,EllP,PQLen,ICrnr
    REAL(DOUBLE),DIMENSION(NAtoms) :: NukeMask

    IF(Q%Leaf)THEN
       EllP=EllFit
       LenP=LenFit
       ChFit=0D0
       DO EllQ=0,Q%HERM%Ell
          IF(Q%HERM%Nq(EllQ).NE.0)THEN 
             PQEll=EllP+EllQ
             PQLen=LHGTF(PQEll)
             Nuke=0
             DO I=1,Q%HERM%Nq(EllQ)
!!$                IF(Q%HERM%Zeta(EllQ)%D(I)==NuclearExpnt)THEN
!!$                   Nuke=Nuke+1
!!$                   NukeMask(Nuke)=Q%HERM%Coef(EllQ)%D(1,I)
!!$                   Q%HERM%Coef(EllQ)%D(1,I)=0D0
!!$                ELSE

                   PiZ=(Pi/Q%HERM%Zeta(EllQ)%D(I))**(3D0/2D0)

                   !                   ChFit(1)=ChFit(1)+Q%HERM%Coef(EllQ)%D(1,I)*PiZ
                   !                   Q%HERM%Coef(EllQ)%D(:,I)=Q%HERM%Coef(EllQ)%D(:,I)*PiZ
!!$                ENDIF
             ENDDO

!             WRITE(*,*)' QHERM + ',Q%HERM%Coef(EllQ)%D(1,1)

             DO ICrnr=1,14
!                WRITE(*,*)' QCent = ',Q%Herm%Cent(EllQ)%D(1,1),' Corner = ',Corners(1,ICrnr)
#ifdef OVERLAP_PROJECTION    
                CALL AAOverlap(HGEll4,Q%HERM%Nq(EllQ),              &
                     EllP,EllQ,LenP,LHGTF(EllQ),PQEll,PQLen,               &
                     NuclearExpnt,1D-100,Q%HERM%IHlf(EllQ)%D(1),             &
                     PCZeta,Q%HERM%Zeta(EllQ)%D(1),Corners(:,ICrnr),       &
                     Q%HERM%Cent(EllQ)%D(1,1),Q%HERM%Coef(EllQ)%D(1,1),ChFit(:,ICrnr))             
#else
                CALL RAhmadiJAlmlof95(HGEll4,Q%HERM%Nq(EllQ),              &
                     EllP,EllQ,LenP,LHGTF(EllQ),PQEll,PQLen,               &
                     NuclearExpnt,1D-100,Q%HERM%IHlf(EllQ)%D(1),             &
                     PCZeta,Q%HERM%Zeta(EllQ)%D(1),Corners(:,ICrnr),       &
                     Q%HERM%Cent(EllQ)%D(1,1),Q%HERM%Coef(EllQ)%D(1,1),ChFit(:,ICrnr))             
!                WRITE(*,*)' COUL = ',ChFit(:,ICrnr)             
#endif

             ENDDO
             Nuke=0
             DO I=1,Q%HERM%Nq(EllQ)
!!$                IF(Q%HERM%Zeta(EllQ)%D(I)==NuclearExpnt)THEN
!!$                   Nuke=Nuke+1
!!$                   Q%HERM%Coef(EllQ)%D(1,I)=NukeMask(Nuke)
!!$                ELSE
                   PiZ=(Pi/Q%HERM%Zeta(EllQ)%D(I))**(-3D0/2D0)
                   !                   Q%HERM%Coef(EllQ)%D(:,I)=Q%HERM%Coef(EllQ)%D(:,I)*PiZ
!!$                ENDIF
             ENDDO
          ENDIF
       ENDDO
    ELSEIF(ASSOCIATED(Q%Left).AND.ASSOCIATED(Q%Right))THEN
       ChFit=PACChargeProject(Q%Left,PCZeta,Corners)+PACChargeProject(Q%Right,PCZeta,Corners)
    ELSEIF(ASSOCIATED(Q%Left))THEN
       ChFit=PACChargeProject(Q%Left,PCZeta,Corners)
    ELSEIF(ASSOCIATED(Q%Right))THEN
       ChFit=PACChargeProject(Q%Right,PCZeta,Corners)
    ENDIF
  END FUNCTION PACChargeProject

  RECURSIVE FUNCTION PACZetaMin(Q) RESULT(ZetaMin)
    TYPE(PoleNode)   :: Q
    REAL(DOUBLE)     :: ZetaMin
    INTEGER          :: I,EllQ
    IF(Q%Leaf)THEN
       ZetaMin=1D20
       DO EllQ=0,Q%HERM%Ell
          DO I=1,Q%HERM%Nq(EllQ)
             ZetaMin=MIN(ZetaMin,Q%HERM%Zeta(EllQ)%D(I))
          ENDDO
       ENDDO
    ELSEIF(ASSOCIATED(Q%Left).AND.ASSOCIATED(Q%Right))THEN
       ZetaMin=MIN(PACZetaMin(Q%Left),PACZetaMin(Q%Right))
    ELSEIF(ASSOCIATED(Q%Left))THEN
       ZetaMin=PACZetaMin(Q%Left)
    ELSEIF(ASSOCIATED(Q%Right))THEN
       ZetaMin=PACZetaMin(Q%Right)
    ENDIF
  END FUNCTION PACZetaMin
  !=================================================================================
  !     
  !=================================================================================
  SUBROUTINE SetSerialPAC(PC,SHG)
    TYPE(ScalarHerm) :: SHG
    TYPE(PAC)                 :: PC
    INTEGER                          :: Ell,L,M,N,LMN,I,Nq,NTot,EllQ,PQLen,PQEll,EllP,LenP
    REAL(DOUBLE)                     :: NodeWeightC,Zeta,Zeta2,PiZ,MinMax,ZetaFac,PCFit
    REAL(DOUBLE)                     :: CramCo,MixMax,ScaledTau,HGInEq,ZZ,CO,PZ,MaxExtent
    REAL(DOUBLE),PARAMETER           :: K3=1.09D0**3
    REAL(DOUBLE),DIMENSION(0:SHG%Ell) :: Wght
    !-------------------------------------------------------------------------------    
    ! H_n(t) < K 2^(n/2) SQRT(n!) EXP(t^2/2), with K=1.09
    !-------------------------------------------------------------------------------    
    PC%Zeta=1D20   
    PC%Wght=Zero
    PC%Zeta=SHG%Zeta
    PiZ = (Pi/PC%Zeta)**(ThreeHalves)
    Wght = Zero
    DO L=0,SHG%Ell
       DO M=0,SHG%Ell-L
          DO N=0,SHG%Ell-L-M
             LMN=LMNDex(L,M,N)
             MixMax=Fact(L)*Fact(M)*Fact(N)
             HGInEq=SQRT(MixMax*(Two*PC%Zeta)**(L+M+N))*SHG%Coef(LMN)
             Wght(L+M+N)=MAX(Wght(L+M+N),ABS(HGInEq))
          ENDDO
       ENDDO
    ENDDO
    PC%Wght=SUM(Wght(0:SHG%Ell))*PiZ
    IF(SHG%Ell.NE.0)PC%Zeta=PC%Zeta*PACFudgeFactor
  END SUBROUTINE SetSerialPAC
  !=================================================================================
  !     
  !=================================================================================
  SUBROUTINE SetSerialMAC(MC,SHG)
    TYPE(ScalarHerm)          :: SHG
    TYPE(MAC)                 :: MC
    INTEGER                   :: Ell,L,M,N,LMN,I,Nq,Ell2
    REAL(DOUBLE)              :: PiZ
    REAL(DOUBLE),DIMENSION(0:HGEll) :: OTmp
    !-------------------------------------------------------------------------------    
    MC%Delta=Zero
    PiZ=(Pi/SHG%Zeta)**(ThreeHalves)
    MC%O=Zero
    DO L=0,SHG%Ell
       DO M=0,SHG%Ell-L
          DO N=0,SHG%Ell-L-M
             Ell=L+M+N
             LMN=LMNDex(L,M,N)
             MC%O(Ell)=MAX(MC%O(Ell),SHG%Coef(LMN)*SHG%Coef(LMN))
          ENDDO
       ENDDO
    ENDDO
    DO Ell=0,SHG%Ell
       MC%O(Ell)=SQRT(MC%O(Ell))*PiZ
    ENDDO
  END SUBROUTINE SetSerialMAC

  SUBROUTINE SetMAC(Q)
    TYPE(PoleNode)        :: Q
    TYPE(MAC)             :: MC
    INTEGER               :: L
    CALL MacProject(Q,Q%POLE%Center,MC)
    Q%MAC%Delta=MC%Delta
    IF(MC%Delta==Zero)THEN	
       Q%MAC%O=Zero
    ELSE
       DO L=0,4
          Q%MAC%O(L)=ABS(MC%O(L)/Q%MAC%Delta**DBLE(MaxPoleEll+1)) 
       ENDDO
    ENDIF
  END SUBROUTINE SetMAC

  RECURSIVE SUBROUTINE MACProject(Q,Center,MC)
    TYPE(PoleNode)        :: Q
    TYPE(MAC)             :: MC,MCLeft,MCRight
    REAL(DOUBLE),DIMENSION(3) :: Center
    REAL(DOUBLE)          :: PiZ,Delta,DeltaEll    
    INTEGER               :: I,Ell,L,M,N,LMN,Nq
    REAL(DOUBLE),DIMENSION(0:4) :: OTmp

    IF(Q%Leaf)THEN
       MC%O=Zero
       MC%Delta=Zero
       DO Ell=0,Q%HERM%Ell
          Nq=Q%HERM%NQ(Ell)
          DO I=1,Nq
             OTmp=0D0
             PiZ=(Pi/Q%HERM%Zeta(Ell)%D(I))**(1.5D0)
             Delta=SQRT(DOT_PRODUCT(Center-Q%HERM%Cent(Ell)%D(:,I),Center-Q%HERM%Cent(Ell)%D(:,I)))
             MC%Delta=MAX(MC%Delta,Delta)
             DO L=0,Ell
                DO M=0,Ell-L
                   DO N=0,Ell-L-M
                      LMN=LMNx(L,M,N)
                      OTmp(L+M+N)=OTmp(L+M+N)+Q%HERM%Coef(Ell)%D(LMN,I)
                   ENDDO
                ENDDO
             ENDDO
             DO L=0,Ell
                DeltaEll=Delta**DBLE(MaxPoleEll)
                MC%O(L)=MC%O(L)+DeltaEll*ABS(OTmp(L))*PiZ
             ENDDO
          ENDDO
       ENDDO
    ELSEIF(ASSOCIATED(Q%Left).AND.ASSOCIATED(Q%Right))THEN
       CALL MACProject(Q%Left,Center,MCLeft)
       CALL MACProject(Q%Right,Center,MCRight)
       MC%Delta=MAX(MCLeft%Delta,MCRight%Delta)
       MC%O=MCLeft%O+MCRight%O
    ELSEIF(ASSOCIATED(Q%Left))THEN
       CALL MACProject(Q%Left,Center,MC)
    ELSEIF(ASSOCIATED(Q%Right))THEN
       CALL MACProject(Q%Right,Center,MC)
    ENDIF
  END SUBROUTINE MACProject

  RECURSIVE FUNCTION DeltaInBox(Q,Center) RESULT(Delta)
    TYPE(PoleNode) :: Q
    REAL(DOUBLE), DIMENSION(3) :: Center,CQ
    REAL(DOUBLE) :: Delta
    INTEGER :: I,Ell
    Delta=0D0
    IF(Q%Leaf)THEN
       DO Ell=0,MaxRhoEll
          DO I=1,Q%Herm%NQ(Ell)
             CQ=Center-Q%Herm%Cent(Ell)%D(:,I)
             Delta=MAX(Delta,SQRT(DOT_PRODUCT(CQ,CQ)))
          END DO
       END DO
    ELSEIF(ASSOCIATED(Q%Left).AND.ASSOCIATED(Q%Right))THEN
       Delta=MAX(DeltaInBox(Q%Left,Center),DeltaInBox(Q%Right,Center))
    ELSEIF(ASSOCIATED(Q%Left))THEN
       Delta=DeltaInBox(Q%Left,Center)
    ELSEIF(ASSOCIATED(Q%Right))THEN
       Delta=DeltaInBox(Q%Right,Center)
    ENDIF
  END FUNCTION DeltaInBox

  !
  FUNCTION NodeWeightC(Ell,Zeta,HGCo) 
    INTEGER                          :: Ell,L,M,N,LMN
    REAL(DOUBLE)                     :: NodeWeightC,Zeta,Zeta2,PiZ,MinMax,ZetaFac
    REAL(DOUBLE), DIMENSION(1:)      :: HGCo
    REAL(DOUBLE)                    :: CramCo,MixMax,ScaledTau,HGInEq
    REAL(DOUBLE),PARAMETER          :: K3=1.09D0**3
    REAL(DOUBLE),DIMENSION(0:12),PARAMETER :: Fact=(/1D0,1D0,2D0,6D0,24D0,120D0,      &
         720D0,5040D0,40320D0,362880D0,   &
         3628800D0,39916800D0,479001600D0/)
    !
    ! H_n(t) < K 2^(n/2) SQRT(n!) EXP(t^2/2), with K=1.09
    PiZ = (Pi/Zeta)**(ThreeHalves)
    NodeWeightC = Zero
    DO L=0,Ell
       DO M=0,Ell-L
          DO N=0,Ell-L-M
             LMN=LMNDex(L,M,N)
             MixMax=Fact(L)*Fact(M)*Fact(N)
             HGInEq=SQRT(MixMax*(Two*Zeta)**(L+M+N))*HGCo(LMN)
             NodeWeightC = NodeWeightC + ABS(HGInEq)
          ENDDO
       ENDDO
    ENDDO
    NodeWeightC = NodeWeightC*PiZ
  END FUNCTION NodeWeightC
  !======================================================================
  !
  !======================================================================
  RECURSIVE SUBROUTINE Print_POLEROOT_LINKS(A)
    TYPE(PoleNode),POINTER :: A
    INTEGER :: NextNumber
    !----------------------------------------------------------------------
    IF(ASSOCIATED(A%Left))CALL Print_POLEROOT_LINKS(A%Left)
    IF(ASSOCIATED(A%Right))CALL Print_POLEROOT_LINKS(A%Right)      
    IF(ASSOCIATED(A%Next))THEN
       NextNumber=A%Next%Box%Number
    ELSE
       NextNumber=-100
    ENDIF
    IF(ASSOCIATED(A%Left).AND.ASSOCIATED(A%Right))THEN
       IF(PPPRINT)WRITE(*,73)A%Box%Tier,A%Box%Number,A%BdexE,A%EdexE,A%Left%Box%Number,A%Right%Box%Number,NextNumber,A%Pole%C(0)
    ELSEIF(ASSOCIATED(A%Left))THEN
       IF(PPPRINT)WRITE(*,74)A%Box%Tier,A%Box%Number,A%BdexE,A%EdexE,A%Left%Box%Number,NextNumber,A%Pole%C(0)
    ELSEIF(ASSOCIATED(A%Right))THEN
       IF(PPPRINT)WRITE(*,77)A%Box%Tier,A%Box%Number,A%BdexE,A%EdexE,A%Right%Box%Number,NextNumber,A%Pole%C(0)
    ELSE
       IF(PPPRINT)WRITE(*,75)A%Box%Tier,A%Box%Number,A%BdexE,A%EdexE,NextNumber,A%Pole%C(0)
    ENDIF
73  FORMAT('xT=',I4,', Num=',I4,', [',I4,',',I4,'], L#=',I4,' R#=',I4,' N#=',I4,' C = ',D12.6)           
77  FORMAT('xT=',I4,', Num=',I4,', [',I4,',',I4,'], R#=',I4,' N#=',I4,' C = ',D12.6)                   
74  FORMAT('xT=',I4,', Num=',I4,', [',I4,',',I4,'], L#=',I4,' N#=',I4,' C = ',D12.6)                   
75  FORMAT('xT=',I4,', Num=',I4,', [',I4,',',I4,'], N#=',I4,' C = ',D12.6)                              
  END SUBROUTINE Print_POLEROOT_LINKS
  !======================================================================
  !
  !======================================================================
  RECURSIVE SUBROUTINE Print_POLEROOT_BOUNDS(A)
    TYPE(PoleNode),POINTER :: A
    INTEGER :: NextNumber,L
    REAL(DOUBLE)           :: Vol,IMx,H1,H2,H3
    !----------------------------------------------------------------------
    IF(ASSOCIATED(A%Left))CALL Print_POLEROOT_BOUNDS(A%Left)
    IF(ASSOCIATED(A%Right))CALL Print_POLEROOT_BOUNDS(A%Right)      
    IF(A%Box%Half(1)==0D0)THEN
       H1=1D0
    ELSE
       H1=A%Box%Half(1)
    ENDIF
    IF(A%Box%Half(2)==0D0)THEN
       H2=1D0
    ELSE
       H2=A%Box%Half(2)
    ENDIF
    IF(A%Box%Half(3)==0D0)THEN
       H3=1D0
    ELSE
       H3=A%Box%Half(3)
    ENDIF
    Vol=H1*H2*H3

    IMx=A%IHalf
    L=A%Herm%Ell
    IF(A%Box%Tier<10)THEN
       IF(ASSOCIATED(A%Left).AND.ASSOCIATED(A%Right))THEN
          WRITE(*,73)A%Box%Tier,A%Box%Number,A%BdexE,A%EdexE,A%Left%Box%Number,A%Right%Box%Number,Vol,L,IMx
       ELSEIF(ASSOCIATED(A%Left))THEN
          WRITE(*,74)A%Box%Tier,A%Box%Number,A%BdexE,A%EdexE,A%Left%Box%Number,Vol,L,IMx
       ELSEIF(ASSOCIATED(A%Right))THEN
          WRITE(*,77)A%Box%Tier,A%Box%Number,A%BdexE,A%EdexE,A%Right%Box%Number,Vol,L,IMx
       ELSE
          WRITE(*,75)A%Box%Tier,A%Box%Number,A%BdexE,A%EdexE,Vol,L,IMx
       ENDIF
    ENDIF
73  FORMAT('xT=',I4,', Num=',I4,', [',I6,',',I6,'], L#=',I4,' R#=',I4,' Vol=',D12.6,', L = ',I2,' I = ',D12.6)           
77  FORMAT('xT=',I4,', Num=',I4,', [',I6,',',I6,'], R#=',I4,' Vol=',D12.6,', L = ',I2,' I  = ',D12.6)                   
74  FORMAT('xT=',I4,', Num=',I4,', [',I6,',',I6,'], L#=',I4,' Vol=',D12.6,', L = ',I2,' I = ',D12.6)                   
75  FORMAT('xT=',I4,', Num=',I4,', [',I6,',',I6,'], Vol=',D12.6,', L = ',I2,' I = ',D12.6)                              
  END SUBROUTINE Print_POLEROOT_BOUNDS





  RECURSIVE SUBROUTINE Print_POLEROOT_LEAF(A)
    TYPE(PoleNode),POINTER :: A
    INTEGER :: NextNumber
    !----------------------------------------------------------------------
    IF(ASSOCIATED(A%Left))CALL Print_POLEROOT_LEAF(A%Left)
    IF(ASSOCIATED(A%Right))CALL Print_POLEROOT_LEAF(A%Right)      
    IF(A%LEAF)THEN
       WRITE(*,73)A%Box%Tier,A%Box%Number,A%PAC%Zeta,A%PAC%Wght,SUM(A%HERM%NQ),A%HERM%NQ
    ENDIF
73  FORMAT('T=',I4,', #=',I4,', Z=',D8.3,', P=',D8.3,', NT = ',I4,', N=',8(I2,","))
  END SUBROUTINE Print_POLEROOT_LEAF




  RECURSIVE SUBROUTINE Print_POLEROOT_PAC(A)
    TYPE(PoleNode),POINTER :: A
    INTEGER :: NextNumber
    !----------------------------------------------------------------------
    IF(ASSOCIATED(A%Left))CALL Print_POLEROOT_PAC(A%Left)
    IF(ASSOCIATED(A%Right))CALL Print_POLEROOT_PAC(A%Right)      
    IF(A%Box%Tier.LT.6)THEN
       IF(ASSOCIATED(A%Left).AND.ASSOCIATED(A%Right))THEN
          WRITE(*,73)A%Box%Number,A%BdexE,A%EdexE,A%Left%Box%Number,A%Right%Box%Number,A%PAC%Zeta,A%PAC%Wght
       ELSEIF(ASSOCIATED(A%Left))THEN
          WRITE(*,74)A%Box%Number,A%BdexE,A%EdexE,A%Left%Box%Number,A%PAC%Zeta,A%PAC%Wght
       ELSEIF(ASSOCIATED(A%Right))THEN
          WRITE(*,77)A%Box%Number,A%BdexE,A%EdexE,A%Right%Box%Number,A%PAC%Zeta,A%PAC%Wght
       ELSE
          WRITE(*,75)A%Box%Number,A%BdexE,A%EdexE,A%PAC%Zeta,A%PAC%Wght
       ENDIF
    ENDIF
73  FORMAT('#=',I4,', [',I8,',',I8,'], L#=',I4,' R#=',I4,' Z=',D12.6,' W = ',D12.6)
77  FORMAT('#=',I4,', [',I8,',',I8,'], R#=',I4,' Z=',D12.6,' W = ',D12.6)
74  FORMAT('#=',I4,', [',I8,',',I8,'], L#=',I4,' Z=',D12.6,' W = ',D12.6)
75  FORMAT('#=',I4,', [',I8,',',I8,'], Z=',D12.6,' W = ',D12.6)
  END SUBROUTINE Print_POLEROOT_PAC

  RECURSIVE SUBROUTINE Print_POLEROOT_MAC(A)
    TYPE(PoleNode),POINTER :: A
    INTEGER :: NextNumber
    !----------------------------------------------------------------------
    IF(ASSOCIATED(A%Left))CALL Print_POLEROOT_MAC(A%Left)
    IF(ASSOCIATED(A%Right))CALL Print_POLEROOT_MAC(A%Right)      
    IF(A%Box%Tier.LT.6)THEN
       IF(ASSOCIATED(A%Left).AND.ASSOCIATED(A%Right))THEN
          WRITE(*,73)A%Box%Number,A%BdexE,A%EdexE,A%Left%Box%Number,A%Right%Box%Number,A%MAC%Delta,A%MAC%O
       ELSEIF(ASSOCIATED(A%Left))THEN
          WRITE(*,74)A%Box%Number,A%BdexE,A%EdexE,A%Left%Box%Number,A%MAC%Delta,A%MAC%O
       ELSEIF(ASSOCIATED(A%Right))THEN
          WRITE(*,77)A%Box%Number,A%BdexE,A%EdexE,A%Right%Box%Number,A%MAC%Delta,A%MAC%O
       ELSE
          WRITE(*,75)A%Box%Number,A%BdexE,A%EdexE,A%MAC%Delta,A%MAC%O
       ENDIF
    ENDIF
73  FORMAT('#=',I4,', [',I8,',',I8,'], L#=',I4,' R#=',I4,' D=',D12.6,' O = ',6(D12.6,", "))
77  FORMAT('#=',I4,', [',I8,',',I8,'], R#=',I4,' D=',D12.6,' O = ',6(D12.6,", "))
74  FORMAT('#=',I4,', [',I8,',',I8,'], L#=',I4,' D=',D12.6,' O = ',6(D12.6,", "))
75  FORMAT('#=',I4,', [',I8,',',I8,'], D=',D12.6,' O = ',6(D12.6,", "))
  END SUBROUTINE Print_POLEROOT_MAC
  !--------------------------------------------------------------
  ! Generate The Estimate
  !--------------------------------------------------------------
  FUNCTION Estimate(LKet,Zt,Coef)
    INTEGER                        :: LKet,LenKet,LKet2
    INTEGER                        :: I,J,K,L,M,N,IJK,LMN,IL,JM,KN
    REAL(DOUBLE)                   :: Coef1,Coef2,Estimate,Zt
    REAL(DOUBLE),DIMENSION(:)      :: Coef
    REAL(DOUBLE),DIMENSION(0:2*LKet)                      :: AuxR
    REAL(DOUBLE),DIMENSION(-1:2*LKet,-1:2*LKet,-1:2*LKet) :: BigR
    CALL SetSameCenterR(Zt,LKet,AuxR,BigR)
    LKet2    = LKet
    Estimate = Zero
    DO K=0,LKet2
       DO J=0,LKet2-K
          DO I=0,LKet2-K-J
             IJK=LMNdex(I,J,K)
             Coef1 = ((-One)**(I+J+K))*Coef(IJK)
             DO N=0,LKet2
                KN = K+N
                DO M=0,LKet2-N
                   JM = J+M
                   DO L=0,LKet2-N-M
                      IL = I+L
                      LMN=LMNdex(L,M,N)
                      Coef2 = Coef(LMN)
                      Estimate = Estimate + Coef1*Coef2*BigR(IL,JM,KN)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    Estimate = SQRT(ABS(Estimate))*Sqrt2Pi5x2*Zt**(-FiveFourths)
  END FUNCTION Estimate
  !--------------------------------------------------------------
  ! Generate R and AuxRfun
  !--------------------------------------------------------------
  SUBROUTINE SetSameCenterR(Zt,LKet,AuxR,BigR)
    INTEGER            :: I0,I1,I2,L,LKet,LKet2
    REAL(DOUBLE)       :: Zt,Omega,o1,o2
    REAL(DOUBLE),DIMENSION(0:2*LKet)                      :: AuxR
    REAL(DOUBLE),DIMENSION(-1:2*LKet,-1:2*LKet,-1:2*LKet) :: BigR
    BigR=Zero
    AuxR=Zero
    Omega  = Half*Zt
    LKet2  = 2*LKet
    o1     = One
    o2     = -Two*Omega
    DO L=0,LKet2
       AuxR(L)=o1/DBLE(2*L+1)
       o1=o1*o2
    ENDDO
    DO L=LKet2,0,-1
       DO I0=LKet2-L,1,-1
          BigR(I0,0,0)=DBLE(I0-1)*BigR(I0-2,0,0)
          DO I1=LKet2-L-I0,1,-1
             BigR(I1,I0,0)=DBLE(I1-1)*BigR(I1-2,I0,0)
             DO I2=LKet2-L-I0-I1,1,-1
                BigR(I2,I1,I0)=DBLE(I2-1)*BigR(I2-2,I1,I0)
             ENDDO
             BigR(0,I1,I0)=DBLE(I1-1)*BigR(0,I1-2,I0)
             BigR(I1,0,I0)=DBLE(I1-1)*BigR(I1-2,0,I0)
          ENDDO
          BigR(0,I0,0)=DBLE(I0-1)*BigR(0,I0-2,0)
          BigR(0,0,I0)=DBLE(I0-1)*BigR(0,0,I0-2)
       ENDDO
       BigR(0,0,0)=AuxR(L)
    ENDDO
  END SUBROUTINE SetSameCenterR


  FUNCTION Estimate2(LTot,Zeta,BCo) RESULT(Est)

    INTEGER :: LTot
    REAL(DOUBLE),DIMENSION(0:LTot*2,0:LTot*2,0:LTot*2)  ::  R
    REAL(DOUBLE),DIMENSION(0:LTot*2)   ::  AuxR
    REAL(DOUBLE),DIMENSION(:) :: BCo
    REAL(DOUBLE) :: Omega,Zeta,Upq,O1,O2,Est,RInt,TTMP
    INTEGER :: LTot2,I,J,K,I0,I1,I2,L,M,N,LMN,Phase,IJK

    Omega=0.5D0*Zeta
    Upq=2.473942945119315D1*Zeta**(-2.5D0)
    LTot2=LTot*2
    o1=1.0D0
    o2=-2.0D0*Omega
    DO J=0,LTot2
       AuxR(J)=o1/DBLE(2*J+1)
       o1=o1*o2
    ENDDO
    DO J=LTot2,0,-1
       DO I0=LTot2-J,1,-1
          IF(I0-1.LE.0)THEN
             R(I0,0,0)=0.0D0
          ELSE
             R(I0,0,0)=DBLE(I0-1)*R(I0-2,0,0)
          ENDIF
          DO I1=LTot2-J-I0,1,-1
             IF(I1-1.LE.0)THEN
                R(I1,I0,0)=0.0D0
             ELSE
                R(I1,I0,0)=DBLE(I1-1)*R(I1-2,I0,0)
             ENDIF
             DO I2=LTot2-J-I0-I1,1,-1
                IF(I2-1.LE.0)THEN
                   R(I2,I1,I0)=0.0D0
                ELSE
                   R(I2,I1,I0)=DBLE(I2-1)*R(I2-2,I1,I0)
                ENDIF
             ENDDO
             IF(I1-1.LE.0)THEN
                R(0,I1,I0)=0.0D0
                R(I1,0,I0)=0.0D0
             ELSE
                R(0,I1,I0)=DBLE(I1-1)*R(0,I1-2,I0)
                R(I1,0,I0)=DBLE(I1-1)*R(I1-2,0,I0)
             ENDIF
          ENDDO
          IF(I0-1.LE.0)THEN
             R(0,I0,0)=0.0D0
             R(0,0,I0)=0.0D0
          ELSE
             R(0,I0,0)=DBLE(I0-1)*R(0,I0-2,0)
             R(0,0,I0)=DBLE(I0-1)*R(0,0,I0-2)
          ENDIF
       ENDDO
       R(0,0,0)=AuxR(J)
    ENDDO

    !      DO LMN=1,LHGTF(LTot)
    !         WRITE(*,2)LMN,BCo(LMN)
    !2        FORMAT(I2,", B R = ",D12.6)
    !      ENDDO

    RInt=0.0D0
    DO K=0,LTot
       DO J=0,LTot-k
          DO I=0,LTot-K-J
             IJK=LMNdex(I,J,K)
             Phase=(-1.0D0)**(I+J+K)
             TTmp=Phase*BCo(IJK)
             DO N=0,LTot
                DO M=0,LTot-N
                   DO L=0,LTot-N-M
                      LMN=LMNdex(L,M,N)
                      RInt=RInt+TTMP*BCo(LMN)*R(I+L,J+M,K+N)

                      !     WRITE(*,3)IJK,LMN,0,BCo(IJK),BCo(LMN),R(I+L,J+M,K+N)
                      ! 3   FORMAT(' A ',3(I3,', '),3(D12.6,', '))



                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    Est=DSQRT(Upq*RInt)

    !   IF(LTot>0)STOP


  END FUNCTION Estimate2



  SUBROUTINE BBoxCheck(P)
    TYPE(PoleNode) :: P
    REAL(DOUBLE) :: Delta
    LOGICAL      :: InBox
    Delta=DeltaInBox(P,P%Pole%Center)
    IF(Delta>P%MAC%Delta)THEN
       STOP 'Fucked delta'
    ENDIF
    InBox=DistInBox(P,P%Box)
    IF(.NOT.InBox)THEN
       STOP 'Fucked box'
    ENDIF
  END SUBROUTINE BBoxCheck

  RECURSIVE FUNCTION DistInBox(Q,Box) RESULT(InBox)
    TYPE(PoleNode) :: Q
    TYPE(BBox)     :: Box
    INTEGER :: I,Ell
    LOGICAL :: InBox
    InBox=.TRUE.
    IF(Q%Leaf)THEN
       DO Ell=0,MaxRhoEll
          DO I=1,Q%Herm%NQ(Ell)
             InBox=InBox.AND.(.NOT.PointOutSideBox(Q%Herm%Cent(Ell)%D(:,I),Box))
             IF(.NOT.InBox)THEN
                WRITE(*,*)' Q%Center = ',Q%Herm%Cent(Ell)%D(:,I)
                WRITE(*,*)' BndBox 1 = ',Box%BndBox(:,1)
                WRITE(*,*)' BndBox 2 = ',Box%BndBox(:,2)
             ENDIF
          END DO
       END DO
    ELSEIF(ASSOCIATED(Q%Left).AND.ASSOCIATED(Q%Right))THEN
       InBox=DistInBox(Q%Left,Box).AND.DistInBox(Q%Right,Box)
    ELSEIF(ASSOCIATED(Q%Left))THEN
       InBox=DistInBox(Q%Left,Box)
    ELSEIF(ASSOCIATED(Q%Right))THEN
       InBox=DistInBox(Q%Right,Box)
    ENDIF
  END FUNCTION DistInBox


!!$
!!$
!!$  SUBROUTINE ChargeSplit(N,Q,C,Ze,ExEst,JS) 
!!$    INTEGER :: N,JS,I,J,K,NAts,NBnd,L,JP,JM
!!$    REAL(DOUBLE) :: IBM,IBP,XBM,XBP
!!$    REAL(DOUBLE) :: PiZ,IB,ChB,XB,Ch,ChLeft,ChRight,ChHalf,Section
!!$    REAL(DOUBLE) :: ChPLeft,ChPRight,ChMLeft,ChMRight
!!$    INTEGER,DIMENSION(:)   :: C
!!$    INTEGER,DIMENSION(1:N) :: Q
!!$    REAL(DOUBLE),DIMENSION(:) :: Ze,ExEst
!!$    REAL(DOUBLE),DIMENSION(1:N) :: Charge,X
!!$    LOGICAL OnAtom,Partitioned,FromPos,POutOfBounds,MOutOfBounds
!!$    REAL(DOUBLE),PARAMETER :: IMax=10D0,XMax=3D0
!!$
!!$    DO J=1,N
!!$       Charge(J)=ExEst(Q(J))
!!$    ENDDO
!!$    CALL DblIntSort77(N,Charge,Q,2)
!!$    ! Split on order of magnitude 
!!$    Section=10D0**(LOG10(Charge(1))+Half*(LOG10(Charge(N))-LOG10(Charge(1))))
!!$    DO J=1,N
!!$       IF(Charge(J)>Section)THEN
!!$          JS=J
!!$          EXIT
!!$       ENDIF
!!$    ENDDO
!!$    IB=DBLE(N-JS)/DBLE(JS)
!!$    IF(IB<0.5D0.OR.IB>2.0D0)THEN
!!$       JS=N/2
!!$    ENDIF
!!$  END SUBROUTINE ChargeSplit
!!$
!!$  FUNCTION BalancedSplit(N,Q,C,Ze,RhoX,Co,JS) RESULT(Partitioned)
!!$    INTEGER :: N,JS,I,J,K,NAts,NBnd,L,JP,JM
!!$    REAL(DOUBLE) :: IBM,IBP,XBM,XBP
!!$    REAL(DOUBLE) :: PiZ,IB,ChB,XB,Ch,ChLeft,ChRight,ChHalf,Section
!!$    REAL(DOUBLE) :: ChPLeft,ChPRight,ChMLeft,ChMRight
!!$    INTEGER,DIMENSION(:)   :: C
!!$    INTEGER,DIMENSION(1:N) :: Q,NewQ
!!$    REAL(DOUBLE),DIMENSION(:) :: RhoX,Ze,Co
!!$    REAL(DOUBLE),DIMENSION(1:N) :: X,CL,Charge
!!$    LOGICAL OnAtom,Partitioned,FromPos,POutOfBounds,MOutOfBounds
!!$    REAL(DOUBLE),PARAMETER :: IMax=20D0,XMax=3D0 !2D0
!!$
!!$    Partitioned=.FALSE.
!!$
!!$    Ch=0D0
!!$    DO J=1,N
!!$       NewQ(J)=Q(J)
!!$       X(J)=RhoX(Q(J))       
!!$       IF(Ze(Q(J)).NE.NuclearExpnt) &
!!$            Ch=Ch+Co(C(Q(J)))*(Pi/Ze(Q(J)))**(3D0/2D0)
!!$    ENDDO
!!$    !   WRITE(*,*)' CHARGE = ',Ch
!!$
!!$    CALL DblIntSort77(N,X,Q,2)
!!$
!!$    ! Eliminate redundancies
!!$
!!$    J=1
!!$    NBnd=0
!!$    Charge=0D0
!!$    DO WHILE(J.LE.N)
!!$       IF(Ze(Q(J)).NE.NuclearExpnt)THEN
!!$          NBnd=NBnd+1
!!$          X(NBnd)=RhoX(Q(J))
!!$          Charge(NBnd)=Co(C(Q(J)))*(Pi/Ze(Q(J)))**(3D0/2D0)
!!$          !          Ch=Ch+Co(C(Q(J)))*(Pi/Ze(Q(J)))**(3D0/2D0)
!!$          DO K=J+1,N
!!$             IF(X(NBnd)==RhoX(Q(K)))THEN
!!$                IF(Ze(Q(K)).NE.NuclearExpnt)THEN
!!$                   Charge(NBnd)=Charge(NBnd)+Co(C(Q(K)))*(Pi/Ze(Q(K)))**(3D0/2D0)
!!$                ENDIF
!!$             ELSE
!!$                EXIT
!!$             ENDIF
!!$             J=J+1
!!$          ENDDO
!!$       ENDIF
!!$       J=J+1
!!$    ENDDO
!!$    IF(ABS((Ch-SUM(Charge(1:NBnd)))/Ch)>1D-10)THEN       
!!$       WRITE(*,*)' - - - - - - - - - - - - - - - - - - - - - '
!!$       WRITE(*,*)' NEl = ',NEl/2
!!$       WRITE(*,*)'2 CHARGE = ',Ch,SUM(Charge(1:NBnd))
!!$       STOP
!!$    ENDIF
!!$    IF(Ch<0D0.OR.Ch>NEl)THEN
!!$       Q=NewQ
!!$       RETURN
!!$    ENDIF
!!$
!!$    Section=X(1)+Half*(X(NBnd)-X(1))
!!$    FromPos=.TRUE.
!!$
!!$    DO J=1,NBnd
!!$       IF(X(J)>Section)THEN
!!$          JS=J-1
!!$          EXIT
!!$       ENDIF
!!$    ENDDO
!!$
!!$    ChLeft=SUM(Charge(1:JS))
!!$    ChRight=SUM(Charge(JS+1:NBnd))
!!$    Ch=ChLeft+ChRight
!!$    ChB=ChRight/ChLeft   
!!$!    WRITE(*,111)JS,ChLeft,ChRight,Ch,ChLeft/ChRight
!!$
!!$    IF(.NOT.(ChRight>0D0.AND.ChLeft>0D0  &   
!!$         .AND.ChRight<Ch .AND.ChLeft<Ch) )THEN
!!$
!!$       JM=JS       
!!$       DO WHILE(JM>0)
!!$          JM=JM-1
!!$          IBM=DBLE(JM)/DBLE(N-JM)
!!$          IF(IBM<1D0/IMax)THEN
!!$             WRITE(*,*)' SUCK 1 MINUS XB = ',XBM,' IBM = ',IBM,IBM<1D0/IMax,1D0/IMax
!!$             EXIT
!!$          ENDIF
!!$          XBM=(X(JM)-X(1))/(X(NBnd)-X(JM))
!!$          IF(XBM<1D0/XMax.OR.XBM>XMax)THEN
!!$             WRITE(*,*)' SUCK 2 MINUS XB = ',XBM,' IBM = ',IBM
!!$             EXIT
!!$          ENDIF
!!$          ChMLeft=SUM(Charge(1:JM))
!!$          ChMRight=SUM(Charge(JM+1:NBnd))
!!$          IF(ChMRight>0D0.AND.ChMLeft>0D0)EXIT
!!$       ENDDO
!!$
!!$       JP=JS
!!$       DO WHILE(JP<NBnd)
!!$          JP=JP+1
!!$          IBP=DBLE(JP)/DBLE(N-JP)
!!$          IF(IBP>IMax)THEN
!!$             WRITE(*,*)' SUCK 1 PLUS IB = ',XBP,' IBP = ',IBP
!!$             EXIT
!!$          ENDIF
!!$          XBP=(X(JP)-X(1))/(X(NBnd)-X(JP))
!!$          IF(XBP<1D0/XMax.OR.XBP>XMax)THEN
!!$             WRITE(*,*)' SUCK 2 PLUS XB = ',XBP,' IBP = ',IBP
!!$             EXIT 
!!$          ENDIF
!!$          ChPLeft=SUM(Charge(1:JP))
!!$          ChPRight=SUM(Charge(JP+1:NBnd))
!!$          IF(ChPRight>0D0.AND.ChPLeft>0D0)EXIT
!!$       ENDDO
!!$
!!$       MOutOfBounds=XBM<=1D0/XMax.OR.XBM>=XMax.OR.IBM<=1D0/IMax.OR.IBM>=IMax
!!$       POutOfBounds=XBP<=1D0/XMax.OR.XBP>=XMax.OR.IBP<=1D0/IMax.OR.IBP>=IMax       
!!$
!!$       WRITE(*,*)' XBM = ',XBM,' IBM = ',IBM
!!$       WRITE(*,*)' MOUNT = ',XBM<1D0/XMax,XBM>XMax,IBM<1D0/IMax,IBM>IMax
!!$       WRITE(*,*)' MOut = ',MOutOfBounds,' POut = ',POutOfBounds
!!$
!!$       IF(MOutOfBounds.AND.POutOfBounds)THEN
!!$          Q=NewQ
!!$          RETURN
!!$       ELSEIF(.NOT.MOutOfBounds.AND..NOT.POutOfBounds)THEN
!!$          IF(JP-JS>JS-JM)THEN
!!$             JS=JM
!!$             FromPos=.FALSE.
!!$!             WRITE(*,*)' --- '
!!$!             WRITE(*,112)JM,ChLeft,ChRight,Ch,ChLeft/ChRight
!!$          ELSE
!!$             JS=JP
!!$             FromPos=.TRUE.
!!$!             WRITE(*,*)' ++++ '
!!$!             WRITE(*,112)JP,ChLeft,ChRight,Ch,ChLeft/ChRight
!!$          ENDIF
!!$       ELSEIF(.NOT.MOutOfBounds)THEN
!!$          JS=JM
!!$          FromPos=.FALSE.
!!$!          WRITE(*,*)' --- '
!!$!          WRITE(*,111)JM,ChLeft,ChRight,Ch,ChLeft/ChRight
!!$       ELSEIF(.NOT.POutOfBounds)THEN
!!$          JS=JP
!!$          FromPos=.TRUE.
!!$!          WRITE(*,*)' ++++ '
!!$!          WRITE(*,111)JP,ChLeft,ChRight,Ch,ChLeft/ChRight
!!$       ENDIF
!!$
!!$    ENDIF
!!$
!!$    Section=X(JS)
!!$
!!$    ChLeft=0D0
!!$    DO J=1,N
!!$       IF(RhoX(Q(J)).LE.Section.AND.Ze(Q(J)).NE.NuclearExpnt)THEN
!!$          ChLeft=ChLeft+Co(C(Q(J)))*(Pi/Ze(Q(J)))**(3D0/2D0)!!$    NAts=0
!!$       ENDIF
!!$    END DO
!!$
!!$    ChRight=0D0
!!$    DO J=1,N
!!$       IF(RhoX(Q(J)).GT.Section.AND.Ze(Q(J)).NE.NuclearExpnt)THEN
!!$          ChRight=ChRight+Co(C(Q(J)))*(Pi/Ze(Q(J)))**(3D0/2D0)!!$    NAts=0
!!$       ENDIF
!!$    END DO
!!$
!!$!    WRITE(*,111)JS,ChLeft,ChRight,ChLeft+ChRight,ChLeft/ChRight
!!$
!!$    DO J=1,N
!!$       X(J)=RhoX(Q(J))
!!$    ENDDO
!!$
!!$    CALL DblIntSort77(N,X,Q,2)       
!!$
!!$    DO J=1,N
!!$       IF(Ze(Q(J)).NE.NuclearExpnt)THEN
!!$          Charge(J)=Co(C(Q(J)))*(Pi/Ze(Q(J)))**(3D0/2D0)
!!$       ENDIF
!!$    END DO
!!$
!!$    IF(FromPos)THEN
!!$       DO J=N,1,-1
!!$          IF(X(J)==Section)THEN
!!$             JS=J
!!$             IBP=DBLE(JS)/DBLE(N-JS)
!!$             IF(IBP<1D0/IMax.OR.IBP>IMax)THEN
!!$                Q=NewQ
!!$                RETURN
!!$             ENDIF
!!$             EXIT
!!$          ENDIF
!!$       ENDDO
!!$       ChLeft=SUM(Charge(1:JS))
!!$       ChRight=SUM(Charge(JS+1:N))
!!$       !       WRITE(*,111)JS,ChLeft,ChRight,Ch,ChLeft/ChRight
!!$    ELSE
!!$       DO J=1,N
!!$          IF(X(J)==Section)THEN
!!$             JS=J
!!$             IBM=DBLE(JS)/DBLE(N-JS)
!!$             IF(IBM<1D0/IMax.OR.IBM>IMax)THEN
!!$                Q=NewQ
!!$                RETURN
!!$             ENDIF
!!$             EXIT
!!$          ENDIF
!!$       ENDDO
!!$       ChLeft=SUM(Charge(1:JS))
!!$       ChRight=SUM(Charge(JS+1:N))
!!$       !       WRITE(*,111)JS,ChLeft,ChRight,Ch,ChLeft/ChRight
!!$    ENDIF
!!$
!!$    ChLeft=0D0
!!$    DO J=1,JS
!!$       IF(Ze(Q(J)).NE.NuclearExpnt)THEN
!!$          ChLeft=ChLeft+Co(C(Q(J)))*(Pi/Ze(Q(J)))**(3D0/2D0)!!$    NAts=0
!!$       ENDIF
!!$    END DO
!!$
!!$    ChRight=0D0
!!$    DO J=JS+1,N
!!$       IF(Ze(Q(J)).NE.NuclearExpnt)THEN
!!$          ChRight=ChRight+Co(C(Q(J)))*(Pi/Ze(Q(J)))**(3D0/2D0)!!$    NAts=0
!!$       ENDIF
!!$    END DO
!!$
!!$!    WRITE(*,112)JS,ChLeft,ChRight,ChLeft+ChRight,ChLeft/ChRight
!!$
!!$    !    WRITE(*,55)RhoX(Q(1)),RhoX(Q(JS)),RhoX(Q(JS+1)),RhoX(Q(N))
!!$    !55  FORMAT(" Split [",D12.6,", ",D12.6,"] [",D12.6,", ",D12.6," ]")
!!$
!!$    Partitioned=.TRUE.                     ! found a good partition
!!$
!!$111 FORMAT('  TRIAL:',I5,' ChLeft = ',D12.6,' ChRight = ',D12.6,' ChTot = ',D12.6,' Balance = ',F6.3)
!!$112 FORMAT('  PASSD:',I5,' ChLeft = ',D12.6,' ChRight = ',D12.6,' ChTot = ',D12.6,' Balance = ',F6.3)
!!$
!!$  END FUNCTION BalancedSplit

END MODULE PoleTree
