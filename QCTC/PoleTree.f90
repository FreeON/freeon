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
    CALL New(Qndx,(/MaxCluster,MaxRhoEll/),(/1,0/))
    CALL New(Cndx,(/MaxCluster,MaxRhoEll/),(/1,0/))
    CALL New(CMTmp,(/MaxCluster,SPLen/),(/1,0/))
    CALL New(SMTmp,(/MaxCluster,SPLen/),(/1,0/))
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
!    IF(Node%NQ<=MaxCluster)THEN
    IF(Node%Leaf)THEN
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
       CALL BoxMerge(P%Left%Box,P%Right%Box,P%Box)
       CALL SetMAC(P)
    ENDIF
  END SUBROUTINE MakePoleTree

  SUBROUTINE SplitPoleBox(Node,Left,Right)
    TYPE(PoleNode), POINTER     :: Node,Left,Right
    REAL(DOUBLE)                :: Section,MaxBox,Extent,MaxExt,IHalfHalf,BalanceQ,TotCh2,ChHalf2,Ex
    REAL(DOUBLE),DIMENSION(3)   :: MaxQL,MaxQH
    INTEGER                     :: Ne,Nn,Be,Ee,Bn,En,Je,Jn,ISplit,Split,I,J,k,IS,NL,NR
    INTEGER                     :: MaxBoxedEll,NBoxedEll,EllSplit,Cd,SplitI,SplitJ
    REAL(DOUBLE)                :: Volume,MaxZeta,MinZeta,ZetaHalf,ChHalf,TotCh,PiZ,RL,RR,RN
    REAL(DOUBLE)                :: BEBalance,ChBalance,XBalance
    INTEGER,DIMENSION(:), ALLOCATABLE :: IList,JList
    REAL(DOUBLE),DIMENSION(:), ALLOCATABLE :: XList,CList
    TYPE(BBox)                  :: BB
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
    !
    IF(Nn==1)THEN
       ! Sort on PAC + HGTF distance from the nuclear center, with the goal of 
       ! creating smaller and smaller boxes that may satisfy both MAC and PAC
       CALL PACSplit(Ne,Qdex(Be:Ee),Ext,Rho%Qx%D,Rho%Qy%D,Rho%Qz%D,GM%Carts%D(:,NDex(Bn)), &
                     MaxCluster,Je)       
       Jn=1
       iSplit=4
    ELSE
       CALL AxisSplit(Ne,Nn,Qdex(Be:Ee),NDex(Bn:En),Ext,Rho%Qx%D,Rho%Qy%D,Rho%Qz%D, &
                      GM%Carts%D,GM%PBC%CellCenter%D,GM%AtNum%D,iSplit,Je,Jn)


!!$       WRITE(*,*)' Jn = ',Jn
!!$       IF(Jn==1)THEN
!!$          CALL PACSplit(Ne,Qdex(Be:Ee),Ext,Rho%Qx%D,Rho%Qy%D,Rho%Qz%D,GM%Carts%D(:,NDex(Bn)), &
!!$                        MaxCluster,Je)       
!!$          iSplit=4 
!!$      ELSE
       IF(iSplit==4)THEN
          CALL ExtSplit(Ne,Qdex(Be:Ee),Ext,MaxCluster,Je)
       ENDIF
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
    ! ... L&R charge count ...
    Left%NQ=Left%EdexE-Left%BdexE+1
    Right%NQ=Right%EdexE-Right%BdexE+1

    ! Split the nuclei 
    IF(ISplit<4)THEN
       ! If we are still doing ORB, then fill in pointer (NDex) locations
       Split=Bn+Jn-1
       Left%BdexN=Bn
       Left%EdexN=Split      
       Right%EdexN=En 
       Right%BdexN=Split+1
    ELSE
       Left%BdexN=Bn
       Left%EdexN=En
       Right%BdexN=Bn
       Right%EdexN=En
    ENDIF
    !
    IF(ISplit==4)THEN
       Left%Leaf=.TRUE.
    ELSEIF(Left%NQ<=MinCluster)THEN
       Left%Leaf=.TRUE.
    ENDIF
    IF(Right%NQ<=MinCluster)Right%Leaf=.TRUE.
    ! ... L&R atom count ...
    Left%NAtms=Left%EdexN-Left%BdexN+1
    Right%NAtms=Right%EdexN-Right%BdexN+1
    !
    K=Qdex(Left%BdexE)
    Left%Box=ExpandPoint((/Rho%Qx%D(K),Rho%Qy%D(K),Rho%Qz%D(K)/),Ext(K) ,Left%Box)
    DO I=Left%BdexE+1,Left%EdexE
       K=Qdex(I)
       BB=ExpandPoint((/Rho%Qx%D(K),Rho%Qy%D(K),Rho%Qz%D(K)/),Ext(K))
       CALL BoxMerge(Left%Box,BB,Left%Box)
    ENDDO

    K=Qdex(Right%BdexE)
    Right%Box=ExpandPoint((/Rho%Qx%D(K),Rho%Qy%D(K),Rho%Qz%D(K)/),Ext(K) ,Right%Box)
    DO I=Right%BdexE+1,Right%EdexE
       K=Qdex(I)
       BB=ExpandPoint((/Rho%Qx%D(K),Rho%Qy%D(K),Rho%Qz%D(K)/),Ext(K))
       CALL BoxMerge(Right%Box,BB,Right%Box)
    ENDDO

    IF(ISplit<4)THEN
       WRITE(*,55)ISplit,Node%Box%Number,Left%BdexE,Left%EDexE,Left%EDexE-Left%BdexE, &
                                         Right%BdexE,Right%EDexE,Right%EDexE-Right%BdexE, &
                                         SUM(GM%AtNum%D(NDex(Left%BdexN:Left%EdexN))),   &
                                         SUM(GM%AtNum%D(NDex(Right%BdexN:Right%EdexN)))
       
55     FORMAT("  Split = ",I2,", Node = ",I4," [",I6,", ",I6,"; ",I3,"] [",I6,", ",I6,"; ",I3," ] // <",F10.4,"><",F10.4,"> ")
    ELSE

!      RL=FurthestPointInBox(GM%Carts%D(:,NDex(Left%BDexN)),Left%Box)
!      RR=FurthestPointInBox(GM%Carts%D(:,NDex(Right%BDexN)),Right%Box)


!      CALL PrintBBox(Right%Box,6)

       RL=BoxMargin(Left%Box)
       RR=BoxMargin(Right%Box)

       NL=Left%EDexE-Left%BdexE+1
       NR=Right%EDexE-Right%BdexE+1
       WRITE(*,56)ISplit,Node%BdexN,GM%AtNum%D(NDex(Node%BdexN)),Left%BdexE,Left%EDexE,NL, &
                                      Right%BdexE,Right%EDexE,NR,RL,RR

56     FORMAT("  Split = ",I2," At = ",I2,", Z = ",F4.1," [",I5,", ",I5,";N = ",I3,"] [",I5,", ",I5,";N = ",I3, &
              "] // <",D10.5,"><",D10.5,"> ")

      IF(RR>RL)STOP

    ENDIF

!    IF(Node%Box%Number>2)STOP

  END SUBROUTINE SplitPoleBox
  !

  SUBROUTINE AxisSplit(Ne,Nn,Qe,Qn,Ex,Qx,Qy,Qz,Nuc,Cntr,Chg,Axis,Je,Jn)
    INTEGER                       :: Ne,Nn,Je,Jn,J,K,Axis,iDim
    INTEGER,DIMENSION(1:Ne)       :: Qe,Qej,QeTmp
    INTEGER,DIMENSION(1:Nn)       :: Qn,Qnj,QnTmp
    REAL(DOUBLE),DIMENSION(:)     :: Ex,Qx,Qy,Qz,Chg
    REAL(DOUBLE),DIMENSION(:,:)   :: Nuc
    REAL(DOUBLE),DIMENSION(3)     :: Cntr
    INTEGER,DIMENSION(1:Ne)       :: Ue
    INTEGER,DIMENSION(1:Nn)       :: Un
    REAL(DOUBLE),DIMENSION(1:Ne)  :: Xe
    REAL(DOUBLE),DIMENSION(1:Nn)  :: Xn
    REAL(DOUBLE)                  :: LpRVol,LpRMrg,LpRSid,NodeSid
    REAL(DOUBLE)                  :: XnLeft,XnRight,XSplit,NodeVol,NodeMrg,MarginMin,MaxDim,MaxExt
    TYPE(BBox)                    :: NodeBox
    TYPE(BBox)                    :: LeftBox,RightBox
    TYPE(BBox),DIMENSION(1:Ne)    :: Bb
    !-------------------------------------------------------------------
    MaxExt = -BIG_DBL
    DO J=1,Ne
       K=Qe(J)
       Bb(J)%BndBox(:,1)=(/Qx(K),Qy(K),Qz(K)/)
       Bb(J)%BndBox(:,2)=(/Qx(K),Qy(K),Qz(K)/)
       Bb(J)=ExpandBox(Bb(J),Ex(K))
       MaxExt=MAX(MaxExt,Ex(K))
    ENDDO

    NodeBox=Bb(1)
    DO J=2,Ne
       CALL BoxMerge(NodeBox,Bb(J),NodeBox)
    ENDDO
    NodeVol=BoxVolume(NodeBox)
    NodeMrg=BoxMargin(NodeBox)

    MaxDim=-BIG_DBL
    DO iDim=1,3
       MaxDim=MAX(MaxDim,NodeBox%BndBox(iDim,2)-NodeBox%BndBox(iDim,1))
!       WRITE(*,*)' BoxSz = ',NodeBox%BndBox(iDim,2)-NodeBox%BndBox(iDim,1)
    ENDDO
    DO iDim=1,3
!       WRITE(*,*)iDIM,MaxDim==(NodeBox%BndBox(iDim,2)-NodeBox%BndBox(iDim,1))
       IF(ABS(MaxDim-(NodeBox%BndBox(iDim,2)-NodeBox%BndBox(iDim,1)))<1D-10)THEN
          Axis=iDim
          Exit
       ENDIF
    ENDDO
    
    XSplit=Zero
    DO J=1,Nn
       K=Qn(J)
       Qnj(J)=J
       Xn(J)=Nuc(Axis,K)-Cntr(Axis)
       XSplit=XSplit+Xn(J)*Chg(K)
    ENDDO
    XSplit=XSplit/SUM(Chg(Qn(1:Nn)))
    !
    CALL DblIntSort77(Nn,Xn,Qnj,2)                    
    !
    DO J=1,Nn
       IF(Xn(J)>XSplit)THEN
          Jn=J-1
          EXIT
       ENDIF
    ENDDO
    IF(ABS(Xn(Nn)-Xn(1))<1D-5)Jn=Nn/2

!!$

    !
    DO J=1,Ne
       Qej(J)=J
       Xe(J)=Bb(J)%Center(Axis)
    ENDDO
    CALL DblIntSort77(Ne,Xe,Qej,2)                    

    DO J=1,Ne
       IF(Xe(J)>XSplit)THEN
          Je=J-1
          EXIT
       ENDIF
    ENDDO    

!!$
!!$    WRITE(*,*)' Jn= ',Jn,' Nn= ',NN
!!$    WRITE(*,*)' Je= ',Je,' Ne= ',Ne
!!$    WRITE(*,*)' XNSplit = ',Xn(1),XSplit,Xn(Nn)
!!$    WRITE(*,*)' XESplit = ',Xe(1),XSplit,Xe(Ne)
    !
    LeftBox=Bb(Qej(1))
    DO J=1,Je
       CALL BoxMerge(LeftBox,Bb(Qej(J)),LeftBox)
    ENDDO

    RightBox=Bb(Qej(Je+1))
    DO J=Je+1,Ne
       CALL BoxMerge(RightBox,Bb(Qej(J)),RightBox)
    ENDDO
!    WRITE(*,*)'Left = ',LeftBox%BndBox(Axis,1),LeftBox%BndBox(Axis,2)
!    WRITE(*,*)'Right = ',RightBox%BndBox(Axis,1),RightBox%BndBox(Axis,2)
!    WRITE(*,*)' Node = ',NodeBox%BndBox(Axis,1),NodeBox%BndBox(Axis,2)

!!$
!!$    WRITE(*,*)' Node = ',NodeBox%BndBox(Axis,2)-NodeBox%BndBox(Axis,1), &
!!$         LeftBox%BndBox(Axis,2) - LeftBox%BndBox(Axis,1) +              &
!!$         RightBox%BndBox(Axis,2)-RightBox%BndBox(Axis,1)
!!$
!!$    WRITE(*,*)' Left = ',LeftBox%BndBox(Axis,2)- LeftBox%BndBox(Axis,1)
!!$    WRITE(*,*)' Right= ',RightBox%BndBox(Axis,2)-RightBox%BndBox(Axis,1)

!!$    WRITE(*,*)' % = ',One/((NodeBox%BndBox(Axis,2)-NodeBox%BndBox(Axis,1))/( &
!!$         (LeftBox%BndBox(Axis,2)- LeftBox%BndBox(Axis,1))+  &
!!$         (RightBox%BndBox(Axis,2)-RightBox%BndBox(Axis,1))))

    NodeSid=NodeBox%BndBox(Axis,2)-NodeBox%BndBox(Axis,1)
    LpRSid=LeftBox%BndBox(Axis,2)- LeftBox%BndBox(Axis,1)+  &
          RightBox%BndBox(Axis,2)-RightBox%BndBox(Axis,1)
    LpRMrg=BoxMargin(LeftBox)+BoxMargin(RightBox)
    LpRVol=BoxVolume(LeftBox)+BoxVolume(RightBox)



!    WRITE(*,33)Axis,LeftBox%BndBox(Axis,2) - LeftBox%BndBox(Axis,1),     &
!                   RightBox%BndBox(Axis,2)-RightBox%BndBox(Axis,1),      &
!                   NodeSid,LpRSid/NodeSid

!    WRITE(*,33)Axis,BoxMargin(LeftBox),BoxMargin(RightBox),BoxMargin(NodeBox),LpRMrg/NodeMrg
33  FORMAT(I3,", LMarg = ",F12.6,' RMarg = ',F12.6," NodeMarg = ",F12.6," % = ",F12.6)
    !----------------------------------------------------
    ! Lots left to do here.  For example, can used a sorted
    ! quadratic search to find min overlap section within
    ! [JeLeft,JeRight], or could do a locally exponential search.
    ! For now, just split the difference.
    !
!    STOP

    IF(LpRMrg>3D0*NodeMrg)THEN

!    WRITE(*,*) LpRSid/NodeSid,LpRSid<1.3D0*NodeSid,Axis,MaxExt
!    WRITE(*,*)'- - - - - - - - - - - - - - - - - - - - - - - - - - - -'
    !
!   IF(LpRSid>1.7D0*NodeSid)THEN
       Axis=4
       RETURN
    ENDIF
    !
    DO J=1,Ne
       QeTmp(J)=Qe(Qej(J))
    ENDDO
    Qe(1:Ne)=QeTmp(1:Ne)
    !
    DO J=1,Nn
       QnTmp(J)=Qn(Qnj(J))
    ENDDO
    Qn(1:Nn)=QnTmp(1:Nn)
    !
  END SUBROUTINE AxisSplit


  SUBROUTINE AxisSplit2(Ne,Nn,Qe,Qn,Ex,Qx,Qy,Qz,Nuc,Chg,Axis,Je,Jn)
    INTEGER                       :: Ne,Nn,Je,Jn,J,K,Axis,iDim,MinDim
    INTEGER,DIMENSION(1:Ne)       :: Qe
    INTEGER,DIMENSION(1:Nn)       :: Qn
    REAL(DOUBLE),DIMENSION(:)     :: Ex,Qx,Qy,Qz,Chg
    REAL(DOUBLE),DIMENSION(:,:)   :: Nuc
    INTEGER,DIMENSION(1:Ne,3)     :: Ue
    INTEGER,DIMENSION(1:Nn,3)     :: Un
    REAL(DOUBLE),DIMENSION(1:Ne)  :: Xe
    REAL(DOUBLE),DIMENSION(1:Nn)  :: Xn
    INTEGER,DIMENSION(3)          :: JnLeft,JnRight,JeLeft,JeRight
    REAL(DOUBLE),DIMENSION(3)     :: LpRVol,LpRMrg
    REAL(DOUBLE)                  :: XnLeft,XnRight,XSplit,NodeVol,NodeMrg,MarginMin

    TYPE(BBox)                    :: NodeBox
    TYPE(BBox),DIMENSION(3)       :: LeftBox,RightBox
    TYPE(BBox),DIMENSION(1:Ne)    :: Bb
    !-------------------------------------------------------------------

    WRITE(*,*)' Ne = ',NE

    DO J=1,Ne
       K=Qe(J)
       Ue(J,:)=J          

       WRITE(*,*)J,' Q = ',(/Qx(K),Qy(K),Qz(K)/)


       Bb(J)%BndBox(:,1)=(/Qx(K),Qy(K),Qz(K)/)
       Bb(J)%BndBox(:,2)=(/Qx(K),Qy(K),Qz(K)/)
       Bb(J)=ExpandBox(Bb(J),Ex(K))    
    ENDDO
    !
    NodeBox=Bb(1)
    DO J=2,Ne
       CALL BoxMerge(NodeBox,Bb(J),NodeBox)
    ENDDO
    NodeVol=BoxVolume(NodeBox)
    NodeMrg=BoxMargin(NodeBox)
!    WRITE(*,*)' NodeVolume = ',NodeVol



    XSplit=0D0
!    WRITE(*,*)' NodeMargin = ',NodeMrg
    !---------------------------------
    DO iDim=1,3
       !---------------------------------
       DO J=1,Ne
          Ue(J,iDim)=J
          Xe(J)=Bb(J)%Center(iDim)
       ENDDO
       CALL DblIntSort77(Ne,Xe,Ue(:,iDim),2)                    
       !---------------------------------
       XSplit=0D0
       DO J=1,Nn
          K=Qn(J)
          Un(J,iDim)=J
          Xn(J)=Nuc(iDim,K)
          XSplit=XSplit+Nuc(iDim,K)*Chg(K)
       ENDDO
       XSplit=XSplit/SUM(Chg(Qn(1:Nn)))
       CALL DblIntSort77(Nn,Xn,Un(:,iDim),2)                    
       DO J=1,Nn
          IF(Xn(J)>XSplit)THEN
             Jn=J-1
             EXIT
          ENDIF
       ENDDO
       ! Left and right should be the split interval 
       ! spanning the sep between two atoms along iDim

!       write(*,*)' XSplit= ',XSplit
!       WRITE(*,*)' Xn = ',Xn
!       WRITE(*,*)' Chg = ',Chg
!       WRITE(*,*)' Jn = ',Jn

       JnLeft(iDim)=Jn
       XnLeft=Xn(Jn)
       XnRight=Xn(Jn+1)

!       WRITE(*,*)' Jn = ',Jn,' Xn(Jn) = ',XnLeft,XSplit,XnRight

       !       WRITE(*,*)' Charge = ',SUM(Chg(Un(1:Jn,iDim))),SUM(Chg(Un(Jn+1:Nn,iDim)))

       !---------------------------------
       ! Find the electron charge sitting between XnLeft and XnRight
       ! for sure, a binary search here is a good idea
       DO J=1,Ne
          IF(Xe(J)>XnLeft)THEN
             JeLeft(iDim)=J
             EXIT
          ENDIF
       ENDDO
       DO J=Ne,1,-1
          IF(Xe(J)<XnRight)THEN
             JeRight(iDim)=J
             EXIT
          ENDIF
       ENDDO

       JeLeft(iDim)=JeLeft(iDim)+(JeRight(iDim)-JeLeft(iDim))/2
       JeRight(iDim)=JeLeft(iDim)+1

       !---------------------------------
       ! Fill in the left and right BBoxes, not including
       ! the charges sitting in the gap.  That inter-gap
       ! split is to be determined later, if we can make this split.
       !
       LeftBox(iDim)=Bb(Ue(1,iDim))
       DO J=1,JeLeft(iDim)
          CALL BoxMerge(LeftBox(iDim),Bb(Ue(J,iDim)),LeftBox(iDim))
       ENDDO
       RightBox(iDim)=Bb(Ue(Ne,iDim))
       DO J=Ne,JeRight(iDim),-1
          CALL BoxMerge(RightBox(iDim),Bb(Ue(J,iDim)),RightBox(iDim))
       ENDDO
       LpRMrg(iDim)=BoxMargin(LeftBox(iDim))+BoxMargin(RightBox(iDim))
       LpRVol(iDim)=BoxVolume(LeftBox(iDim))+BoxVolume(RightBox(iDim))

!       WRITE(*,*)' L  = ',BoxMargin(LeftBox(iDim)),' R = ',BoxMargin(RightBox(iDim))
!       WRITE(*,*)' L+R Margin = ',LpRMrg(iDim),' (L+R)/Node = ',LpRMrg(iDim)/NodeMrg

       !
    ENDDO
    !
    MinDim=0
    MarginMin=BIG_DBL
    DO iDim=1,3
       IF(LpRMrg(iDim)<1.85D0*NodeMrg)THEN
          MarginMin=MIN(MarginMin,LpRMrg(iDim))
          IF(iDim==1)THEN
             MinDim=1
          ELSEIF(LpRMrg(iDim)<LpRMrg(iDim-1))THEN
             MinDim=iDim
          ENDIF
       ENDIF
    ENDDO

!    IF(MinDim.NE.0)WRITE(*,33)MinDim,LpRMrg(MinDim),BoxMargin(LeftBox(MinDim)),BoxMargin(RightBox(MinDim)),BoxMargin(NodeBox)
33  FORMAT(I3,", ",6(D12.6,", "))
    !----------------------------------------------------
    ! Lots left to do here.  For example, can used a sorted
    ! quadratic search to find min overlap section within
    ! [JeLeft,JeRight], or could do a locally exponential search.
    ! For now, just split the difference.
    !
    IF(MinDim==0)THEN
       Axis=4
       RETURN
    ELSE
       ! We are splitting along this Cartesian axis.
       ! Sort the nuclear and electronic indexes according
       ! to the best split, using idim+1 or idim-1 as scratch      
       Axis=MinDim
       IF(MinDim==1)THEN
          DO J=1,Ne
             Ue(J,MinDim+1)=Qe(Ue(J,MinDim))
          ENDDO
          Qe(:)=Ue(:,MinDim+1)
          DO J=1,Nn
             Un(J,MinDim+1)=Qn(Un(J,MinDim))
          ENDDO
          Qn(:)=Un(:,MinDim+1)
       ELSE
          DO J=1,Ne
             Ue(J,MinDim-1)=Qe(Ue(J,MinDim))
          ENDDO
          Qe(:)=Ue(:,MinDim-1)
          DO J=1,Nn
             Un(J,MinDim-1)=Qn(Un(J,MinDim))
          ENDDO
          Qn(:)=Un(:,MinDim-1)
       ENDIF
       !
       Jn=JnLeft(MinDim)
       Je=JeLeft(MinDim)!+(JeRight(MinDim)-JeLeft(MinDim))/2
    ENDIF
    !
  END SUBROUTINE AxisSplit2


  SUBROUTINE ExtSplit(Ne,Qe,Ex,MxCluster,Je)
    INTEGER :: Ne,Nn,Je,Jn,J,ISplit,K,MxCluster
    REAL(DOUBLE) :: Section
    INTEGER,DIMENSION(1:Ne) :: Qe
    REAL(DOUBLE),DIMENSION(:) :: Ex
    REAL(DOUBLE),DIMENSION(1:Ne) :: X
    !
    DO J=1,Ne
       K=Qe(J)
       X(J)=Ex(K)
    ENDDO
    ! Sort with increasing first, so that we terminate the tree as early as possible with
    ! as many small distributions as possible taken care of using the multipole approximation
    CALL DblIntSort77(Ne,X,Qe,-2)             
    !    
    Section=Half*X(1)
    !
    DO J=1,Ne
       IF(X(J)<Section)THEN
          Je=J-1
          EXIT
       ENDIF
    ENDDO

!    WRITE(*,*)' X1 = ',X(1),X(Je),' X(Ne) = ',X(Ne)

    ! Here we may overide the above bisection,  just pulling off the CluserSize 
    ! large distance distributions from the bottom of the list.  Thus, the left
    ! node generated from this split must be a leaf:
    IF(Je>MxCluster)THEN       
       Je=MxCluster
!       WRITE(*,*)' X1 = ',X(1),X(Je),' X(Ne) = ',X(Ne)
    ENDIF
    !
  END SUBROUTINE ExtSplit



  SUBROUTINE PACSplit(Ne,Qe,Ex,Qx,Qy,Qz,Nuc,MxCluster,Je)
    INTEGER :: Ne,Nn,Je,Jn,J,ISplit,K,MxCluster
    REAL(DOUBLE) :: Section
    INTEGER,DIMENSION(1:Ne) :: Qe
    REAL(DOUBLE),DIMENSION(:) :: Ex,Qx,Qy,Qz
    REAL(DOUBLE),DIMENSION(3) :: Nuc
    REAL(DOUBLE),DIMENSION(1:Ne) :: X
    TYPE(BBox)                   :: Box
    !
    DO J=1,Ne
       K=Qe(J)
       Box%BndBox(:,1)=Nuc
       Box%BndBox(:,2)=Nuc
       CALL BoxMerge( Box, ExpandPoint((/Qx(K),Qy(K),Qz(K)/),Ext(K)) , Box )
       X(J)=FurthestPointInBox(Nuc,Box)       
!       X(J)=BoxMargin(Box)       
    ENDDO
    ! Sort with increasing first, so that we terminate the tree as early as possible with
    ! as many small distributions as possible taken care of using the multipole approximation
    CALL DblIntSort77(Ne,X,Qe,-2)             
    !
    Box%BndBox(:,1)=Nuc
    Box%BndBox(:,2)=Nuc
    DO J=Ne,1,-1
       K=Qe(J)
       CALL BoxMerge( Box, ExpandPoint((/Qx(K),Qy(K),Qz(K)/),Ext(K)) , Box )
!       X(J)=BoxVolume(Box)       
       X(J)=BoxMargin(Box)       
!       X(J)=FurthestPointInBox(Nuc,Box)       
    ENDDO
    ! Sort with increasing first, so that we terminate the tree as early as possible with
    ! as many small distributions as possible taken care of using the multipole approximation
    CALL DblIntSort77(Ne,X,Qe,-2)             
    !    
    Section=X(Ne)+Half*(X(1)-X(Ne))
    !
    DO J=1,Ne
       IF(X(J)<Section)THEN
          Je=J-1
          EXIT
       ENDIF
    ENDDO

!    WRITE(*,*)' Section = ',Section

!    WRITE(*,*)' X = ',X(1),X(Je),X(Je+1),X(Ne)


    ! Here we may overide the above bisection,  just pulling off the CluserSize 
    ! large distance distributions from the bottom of the list.  Thus, the left
    ! node generated from this split must be a leaf:
    IF(Je>MxCluster)Je=MxCluster

!!$
!!$    Box%BndBox(:,1)=Nuc
!!$    Box%BndBox(:,2)=Nuc
!!$    DO J=Ne,Je+1,-1
!!$       K=Qe(J)
!!$       CALL BoxMerge( Box, ExpandPoint((/Qx(K),Qy(K),Qz(K)/),Ext(K)) , Box )
!!$       WRITE(*,*)J,X(J),FurthestPointInBox(Nuc,Box)       
!!$    ENDDO
!!$


!    CALL PrintBBox(Box,6)

    !
  END SUBROUTINE PACSplit

  SUBROUTINE ChargeSplit(Ne,Qe,RhoX,Je,Nn,ISplit,Qn,NucXYZ,Charge,Jn)
    INTEGER :: Ne,Nn,Je,Jn,J,ISplit
    REAL(DOUBLE) :: Section
    INTEGER,DIMENSION(1:Ne) :: Qe
    INTEGER,DIMENSION(1:Nn) :: Qn
    REAL(DOUBLE),DIMENSION(:)   :: RhoX,Charge
    REAL(DOUBLE),DIMENSION(:,:) :: NucXYZ
    REAL(DOUBLE),DIMENSION(1:Ne) :: X
    !

    DO J=1,Nn
       X(J)=NucXYZ(ISplit,Qn(J))
    ENDDO
    CALL DblIntSort77(Nn,X,Qn,2)             
    !
    Section=Zero
    DO J=1,Nn
       Section=Section+Charge(Qn(J))*X(J)       
    ENDDO
    Section=Section/SUM(Charge(Qn(1:Nn)))
    !
    DO J=1,Nn
       IF(X(J)>Section)THEN
          Jn=J-1
          EXIT
       ENDIF
    ENDDO
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
  END SUBROUTINE ChargeSplit


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
    !
    ! This for MaxEll
    iHGStack=iHGStack+1
    !
    delta=0D0
    Node%POLE%Charge=GM%AtNum%D(NDex(Node%BdexN))
    Node%POLE%Center=GM%Carts%D(:,NDex(Node%BdexN))
    Node%BOX%BndBOX(:,1)=Node%POLE%Center
    Node%BOX%BndBOX(:,2)=Node%POLE%Center
    ! 
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
       !
       Ex=Ext(Qd)
       Node%BOX%BndBOX(1,1)=MIN(Node%BOX%BndBOX(1,1),Rho%Qx%D(Qd)-Ex)
       Node%BOX%BndBOX(1,2)=MAX(Node%BOX%BndBOX(1,2),Rho%Qx%D(Qd)+Ex)
       Node%BOX%BndBOX(2,1)=MIN(Node%BOX%BndBOX(2,1),Rho%Qy%D(Qd)-Ex)
       Node%BOX%BndBOX(2,2)=MAX(Node%BOX%BndBOX(2,2),Rho%Qy%D(Qd)+Ex)
       Node%BOX%BndBOX(3,1)=MIN(Node%BOX%BndBOX(3,1),Rho%Qz%D(Qd)-Ex)
       Node%BOX%BndBOX(3,2)=MAX(Node%BOX%BndBOX(3,2),Rho%Qz%D(Qd)+Ex)
       !
       delta=MAX(delta,SQRT(DOT_PRODUCT( (/Rho%Qx%D(Qd)-NODE%Pole%Center(1),Rho%Qy%D(Qd) &
                               -NODE%Pole%Center(2),Rho%Qz%D(Qd)-NODE%Pole%Center(3)/) , &
                                         (/Rho%Qx%D(Qd)-NODE%Pole%Center(1),Rho%Qy%D(Qd) &
                                -NODE%Pole%Center(2),Rho%Qz%D(Qd)-NODE%Pole%Center(3)/) )))
       !
    ENDDO
    Node%POLE%delta=delta
    Node%BOX%Half(:)=Half*(Node%BOX%BndBOX(:,2)-Node%BOX%BndBOX(:,1))
    Node%BOX%Center(:)=Node%BOX%BndBOX(:,1)+Node%Box%Half    
    !
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

!!$
!!$
!!$
!!$       Ex=Ext(Qd)
!!$       
!!$       Ex=Extent(Ell,Node%Herm%Zeta(Ell)%D(I),Node%Herm%Coef(Ell)%D(:,I),Tau_O=TauPAC,Potential_O=.TRUE.,ExtraEll_O=0)
!!$
!!$       IF(Ex.NE.Ext(Qd))THEN
!!$          WRITE(*,*)Ex,Ext(Qd)
!!$          STOP
!!$       ENDIF





          ENDDO
       ENDIF
    ENDDO
    ! Fill in the multiPOLE 
    Node%POLE%Ell=MaxPoleEll
    CALL AllocSP(Node%POLE)
    ! Double check we have a valid expansion center (ie in the box!)
    IF(PointOutSideBox(Node%POLE%Center,Node%BOX))THEN
       WRITE(*,32)Node%POLE%Center
32     FORMAT(' POLE%Center = ',3(D12.6,', '))
       CALL PrintBBox(Node%BOX,6)
       CALL Halt(' In FillRhoLeaf: Multipole center outside of BBox ')
    ENDIF
    !
    CALL HGToSP_POLENODE(Node%HERM,Node%POLE)
    CALL SetMAC(Node)

    NGaussTotal=NGaussTotal+(Node%EdexE-Node%BdexE+1)
    !
!    WRITE(*,33)Node%BOX%Number,Node%POLE%Charge,Node%EdexE-Node%BdexE+1,NODE%POLE%Delta, &
!               FurthestPointInBox(GM%Carts%D(:,NDex(Node%BDexN)),Node%Box)

!    WRITE(*,44)Node%Box%BndBox(:,1),Node%Box%BndBox(:,2),
!44     format('Box = /',3(D10.4,', '),'/,/',3(D10.4,', '),' C= ',3(D10.4,', ') )

33  FORMAT(' Node = ',I4,' Z = ',F4.1,' NGauss = ',I4,' Delta = ',D8.2,' BoxSz = ',D8.2)
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
          EX=Extent(Ell,ZE,Rho%Co%D(CD:CD+LMNLen-1),Tau_O=TauPAC,Potential_O=.TRUE.,ExtraEll_O=0)
          IF(EX>Zero)THEN
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

END MODULE PoleTree
