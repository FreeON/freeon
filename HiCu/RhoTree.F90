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
!    BUILD A HIERARCHICAL REPRESENTATION OF THE ELECTRON DENSITY USING
!    THE K-D TREE DATA STRUCTURE TO ENABLE EFFICIENT RANGE QUERRIES TO
!    ACCESS MINIMALLY LOCAL PORTIONS OF THE DENSITY
!
!    Author: Matt Challacombe
!------------------------------------------------------------------------------
MODULE RhoTree
  USE Derivedtypes
  USE GlobalScalars
  USE GlobalObjects
  USE ProcessControl
  USE Indexing
  USE InOut
  USE Macros
  USE Thresholding
  USE HiCuThresholds
  USE BoundingBox
  USE CubeGrid
  USE MondoLogger

  IMPLICIT NONE

  !=================================================================================================
  !  Density node
  !=================================================================================================
  TYPE RhoNode
    Logical                           :: Leaf         ! Is this a data containing node?
    INTEGER                           :: Bdex         ! Begining index of ORB list for this node
    INTEGER                           :: Edex         ! ENDign index of ORB list for this node
    INTEGER                           :: NQ           ! Number of centers
    REAL(DOUBLE)                      :: Extent       ! Penetration extent
    TYPE(BBox)                        :: Box          ! Bounding box
    INTEGER                           :: Ell          ! Ell type
    REAL(DOUBLE)                      :: Zeta         ! Exponent
    REAL(DOUBLE)                      :: Qx,Qy,Qz     ! Position
#ifdef POINTERS_IN_DERIVED_TYPES
    REAL(DOUBLE),DIMENSION(:),POINTER :: Co           ! Coefficients of the HGTF density
#else
    REAL(DOUBLE),DIMENSION(:), &
         ALLOCATABLE :: Co           ! Coefficients of the HGTF density
#endif
    TYPE(RhoNode),            POINTER :: Travrse      ! Next node in tree traversal
    TYPE(RhoNode),            POINTER :: Descend      ! Next node in tree descent
  END TYPE RhoNode
  !=================================================================================================
  !  Globals
  !=================================================================================================
  TYPE(RhoNode), POINTER                :: RhoRoot     ! Root of the tree
  INTEGER                               :: RhoNodes    ! Number of nodes in the RhoTree
  INTEGER                               :: RhoLevel    ! Number of tiers in the RhoTree
  INTEGER                               :: CurrentTier ! Current tier for level by level addressing
  TYPE(HGRho)                           :: Rho         ! Density
  INTEGER,     DIMENSION(:),Allocatable :: Qdex        ! Distribution pointer
  INTEGER,     DIMENSION(:),Allocatable :: Cdex        ! Coefficient pointer
  INTEGER,     DIMENSION(:),Allocatable :: Ldex        ! Ell index
  REAL(DOUBLE),DIMENSION(:),Allocatable :: RList       ! Bisection array
  REAL(DOUBLE),DIMENSION(:),Allocatable :: Zeta        ! Exponent list
  REAL(DOUBLE),DIMENSION(:),Allocatable :: Ext         ! Extent array
  !-----------!
CONTAINS !
  !=================================================================================
  !     Build the RhoTree from Rho in HGTF array form
  !=================================================================================
  SUBROUTINE RhoToTree(Args)
    TYPE(ARGMT)               :: Args
    INTEGER                   :: NSDen,Status,K,I
    !--------------------------------------------------------------------------------
    !        Initialize the density
    CALL InitRho(Args)
    NSDen=Rho%NSDen
    !        Initialize counters
    RhoNodes=0
    RhoLevel=0
    !        Initialize the root node
    CALL NewRhoNode(RhoRoot,0)
    CALL InitRhoRoot
    !        Recursively partition the density into a 4-D BinTree
    CALL SplitRho(RhoRoot)
    !        Recursively construct BBoxs for the RhoTree
    DO CurrentTier=RhoLevel,0,-1
      CALL MergeRhoBBox(RhoRoot)
    ENDDO
    !        Delete the arrayed version of the density
    CALL DeleteRho
  END SUBROUTINE RhoToTree
  !===============================================================================
  !     Initialize the root of the RhoTree
  !===============================================================================
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
    !        Set the center and width of the Cartesian part
    RhoRoot%Box%Half=Half*(RhoRoot%Box%BndBox(1:3,2)-RhoRoot%Box%BndBox(1:3,1))
    RhoRoot%Box%Center=Half*(RhoRoot%Box%BndBox(1:3,2)+RhoRoot%Box%BndBox(1:3,1))
  END SUBROUTINE InitRhoRoot
  !===================================================================
  !     Recursively partition the density into a 4-D BinTree
  !===================================================================
  RECURSIVE SUBROUTINE SplitRho(Node)
    TYPE(RhoNode), POINTER :: Node,Left,Right
    !--------------------------------------------------------------
    IF(Node%NQ==1)THEN
      CALL FillRhoLeaf(Node)
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
  !     Orthogonal Recusive Bisection
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
    !        Choose direction to section
    ISplit=Mod(Node%Box%Tier,4)+1
    IF(ISplit==1)THEN
      !           Split X dir
      DO I=B,E
        J=J+1
        K=Qdex(I)
        RList(J)=Rho%Qx%D(K)
      ENDDO
    ELSEIF(ISplit==2)THEN
      !           Split on Y dir
      DO I=B,E
        J=J+1
        K=Qdex(I)
        RList(J)=Rho%Qy%D(K)
      ENDDO
    ELSEIF(ISplit==3)THEN
      !           Split on Z dir
      DO I=B,E
        J=J+1
        K=Qdex(I)
        RList(J)=Rho%Qz%D(K)
      ENDDO
    ELSEIF(ISplit==4)THEN
      !           Split on box size
      DO I=B,E
        J=J+1
        K=Qdex(I)
        RList(J)=Ext(K)
      ENDDO
    ENDIF
    !        Sort
    CALL DblIntSort77(N,RList,Qdex(B:E),2)
    !        Orthogonal recursive bisection (ORB)
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
    !        Set counters
    Left%NQ=Left%Edex-Left%Bdex+1
    Right%NQ=Right%Edex-Right%Bdex+1
  END SUBROUTINE SplitRhoBox
  !================================================================================
  !     Recursively build up BBoxes for each node in the RhoTree
  !================================================================================
  RECURSIVE SUBROUTINE MergeRhoBBox(Node)
    TYPE(RhoNode)         :: Node
    TYPE(RhoNode),POINTER :: Left,Right
    REAL(DOUBLE)          :: Ex
    INTEGER               :: K
    !--------------------------------------------------------------------------
    IF(Node%Box%Tier==CurrentTier.AND.Node%Leaf)THEN
      Node%Box%BndBox(:,1)=(/Node%Qx,Node%Qy,Node%Qz/)
      Node%Box%BndBox(:,2)=(/Node%Qx,Node%Qy,Node%Qz/)
      Node%Box=ExpandBox(Node%Box,Node%Extent)
    ELSEIF(Node%Leaf)THEN
      RETURN
    ELSEIF(Node%Box%Tier==CurrentTier)THEN
      Left=>Node%Descend
      Right=>Node%Descend%Travrse
      DO K=1,3
        Node%Box%BndBox(K,1)=MIN(Left%Box%BndBox(K,1),Right%Box%BndBox(K,1))
        Node%Box%BndBox(K,2)=MAX(Left%Box%BndBox(K,2),Right%Box%BndBox(K,2))
      ENDDO
      Node%Box%Half(:)=Half*(Node%Box%BndBox(:,2)-Node%Box%BndBox(:,1))
      Node%Box%Center(:)=Half*(Node%Box%BndBox(:,2)+Node%Box%BndBox(:,1))
    ELSE
      CALL MergeRhoBBox(Node%Descend)
      CALL MergeRhoBBox(Node%Descend%Travrse)
    ENDIF
  END SUBROUTINE MergeRhoBBox
  !=====================================================================================
  !     Fill leaf nodes with data
  !=====================================================================================
  SUBROUTINE FillRhoLeaf(Node)
    TYPE(RhoNode), POINTER :: Node
    INTEGER                :: I,IQ,IC,J,K,KQ,KC,L,B,E,N,NQ,NC,LMNLen,LTmp,Status,iSDen
    REAL(DOUBLE) :: RhoSum
    Interface
      SUBROUTINE DblIntSort77(N,X,Y,Ordr)
        USE DerivedTYPEs
        INTEGER,                  Intent(IN)    :: N,Ordr
        REAL(DOUBLE),DIMENSION(N),Intent(INOUT) :: X
        INTEGER,     DIMENSION(N),Intent(INOUT) :: Y
      END SUBROUTINE DblIntSort77
    END Interface
    !-------------------------------------------------------------------------------------
    !        Set leaf logical
    Node%Leaf=.True.
    !        Boundaries in the ordered lists
    B=Node%Bdex
    E=Node%Edex
    NQ=E-B+1
    !        Stupid check
    IF(NQ/=1.OR.B/=E)  &
         CALL Halt('Bad Logic in FillRhoLeaf ')
    !        Filler up...
    KQ=Qdex(B)
    KC=Cdex(KQ)
    Node%Ell=Ldex(KQ)
    Node%Zeta=Zeta(KQ)
    Node%Extent=Ext(KQ)
    Node%Qx=Rho%Qx%D(KQ)
    Node%Qy=Rho%Qy%D(KQ)
    Node%Qz=Rho%Qz%D(KQ)
    !        Allocate HGTF coefficients array
    LMNLen=LHGTF(Node%Ell)
    ALLOCATE(Node%Co(1:LMNLen*Rho%NSDen),STAT=Status)!<<< SPIN
    CALL IncMem(Status,0,LMNLen*Rho%NSDen)
    !        Transfer data in
    DO iSDen=1,Rho%NSDen
      CALL DCOPY(LMNLen,Rho%Co%D(KC+(iSDen-1)*Rho%NCoef),1,Node%Co(1+(iSDen-1)*LMNLen),1)
      !old            Node%Co(1:LMNLen)=Rho%Co%D(KC:KC+LMNLen-1)
    ENDDO
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
    RhoLevel=MAX(RhoLevel,Level)
    Node%Box%Number=RhoNodes
    RhoNodes=RhoNodes+1
    NULLIFY(Node%Travrse)
    NULLIFY(Node%Descend)
#ifdef POINTERS_IN_DERIVED_TYPES
    NULLIFY(Node%Co)
#endif
  END SUBROUTINE NewRhoNode
  !==========================================================================
  !     Recusrively delete RhoTree
  !==========================================================================
  RECURSIVE SUBROUTINE DeleteRhoTree(Node)
    TYPE(RhoNode), POINTER   :: Node
    INTEGER :: Status
    IF(Node%Leaf)THEN
      DEALLOCATE(Node%Co,STAT=Status)
      CALL DecMem(Status,LHGTF(Node%Ell),0)
      IF(Status/=SUCCEED) &
           CALL Halt(' Leaf Node DEALLOCATE failed in DeleteRhoNode ')
      NULLIFY(Node)
    ELSEIF(ASSOCIATED(Node%Descend%Travrse))THEN
      CALL DeleteRhoTree(Node%Descend%Travrse)
      NULLIFY(Node%Descend%Travrse)
    ELSEIF(ASSOCIATED(Node%Descend))THEN
      CALL DeleteRhoTree(Node%Descend)
      NULLIFY(Node%Descend)
    ENDIF
  END SUBROUTINE DeleteRhoTree
  !========================================================================================
  !     ALLOCATE and read in the density, initalize global lists
  !========================================================================================
  SUBROUTINE InitRho(Args)
    TYPE(ARGMT)  :: Args
    INTEGER      :: z,oq,or,iq,NQ,Q,Ell,Status,I,IOS,LMNLen,QD,CD
    REAL(DOUBLE) :: ZE,EX,Dummy
    !----------------------------------------------------------------------------------------
    !        Get the density

    ! Index check...
    IF(SCFActn=='ForceEvaluation')THEN
#ifdef PARALLEL
      Open(UNIT=Seq,FILE=TrixFile('Rho'//IntToChar(MyID),Args,1),STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
#else
      CALL MondoLog(DEBUG_MEDIUM, "InitRho", "Index Check: opening "//TRIM(TrixFile('Rho',Args,1)))
      Open(UNIT=Seq,FILE=TrixFile('Rho',Args,1),STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
#endif
    ELSE
#ifdef PARALLEL
      Open(UNIT=Seq,FILE=TrixFile('Rho'//IntToChar(MyID),Args,0),STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
#else
      CALL MondoLog(DEBUG_MEDIUM, "InitRho", "Index Check: opening "//TRIM(TrixFile('Rho',Args,0)))
      Open(UNIT=Seq,FILE=TrixFile('Rho',Args,0),STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
#endif
    ENDIF
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)Rho%NSDen,Rho%NExpt,Rho%NDist,Rho%NCoef!<<< SPIN
    CALL MondoLog(DEBUG_MAXIMUM, "InitRho", "Rho%NSDen = "//TRIM(IntToChar(Rho%NSDen))//", " &
      //"Rho%NExpt = "//TRIM(IntToChar(Rho%NExpt))//", " &
      //"Rho%NDist = "//TRIM(IntToChar(Rho%NDist))//", " &
      //"Rho%NCoef = "//TRIM(IntToChar(Rho%NCoef))//", " &
      //"MyID = "//TRIM(IntToChar(MyID)))
    CALL New(Rho%NQ  ,Rho%NExpt)
    CALL New(Rho%OffQ,Rho%NExpt)
    CALL New(Rho%OffR,Rho%NExpt)
    CALL New(Rho%Lndx,Rho%NExpt)
    CALL New(Rho%Expt,Rho%NExpt)
    CALL New(Rho%Qx,  Rho%NDist)
    CALL New(Rho%Qy,  Rho%NDist)
    CALL New(Rho%Qz,  Rho%NDist)
    CALL New(Rho%Co,  Rho%NCoef*Rho%NSDen)!<<< SPIN
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%NQ%I  (i),i=1,Rho%NExpt)
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%OffQ%I(i),i=1,Rho%NExpt)
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%OffR%I(i),i=1,Rho%NExpt)
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%Lndx%I(i),i=1,Rho%NExpt)
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%Expt%D(i),i=1,Rho%NExpt)
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%Qx%D  (i),i=1,Rho%NDist)
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%Qy%D  (i),i=1,Rho%NDist)
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%Qz%D  (i),i=1,Rho%NDist)
    Read(UNIT=Seq,Err=202,IOSTAT=IOS)(Rho%Co%D  (i),i=1,Rho%NCoef*Rho%NSDen)!<<< SPIN
    Close(UNIT=Seq,STATUS='KEEP')
    !        ALLOCATE global arrays
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
    DO z=1,Rho%NExpt-1
      oq =Rho%OffQ%I(z)
      or =Rho%OffR%I(z)
      Ell=Rho%Lndx%I(z)
      ZE =Rho%Expt%D(z)
      LMNLen=LHGTF(Ell)
      DO Q=1,Rho%NQ%I(z)
        QD=oq+Q
        CD=or+(Q-1)*LMNLen+1
        EX=Extent(Ell,ZE,Rho%Co%D(CD:CD+LMNLen-1),TauRho,ExtraEll_O=1)
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
    !        Redefine NDist
    Rho%NDist=IQ-1
    !        Normal exit
    RETURN
    !        Hurl on IO error
202 CALL Halt('Died in PutRho, IOSTAT = '//Trim(IntToChar(IOS)))
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
    DEALLOCATE(Ext,STAT=Status)
    CALL DecMem(Status,0,Rho%NDist)
    DEALLOCATE(RList,STAT=Status)
    CALL DecMem(Status,0,Rho%NDist)
  END SUBROUTINE DeleteRho
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
END MODULE RhoTree
