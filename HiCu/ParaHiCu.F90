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
! Parallel HiCu  --- author: Chee K. Gan (2002 June 25)

#include "MondoConfig.h"

MODULE ParallelHiCu
#ifdef PARALLEL
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalObjects
  USE ProcessControl
  USE Indexing
  USE InOut
  USE Macros
  USE CubeGrid
  USE BoundingBox
  USE RhoTree
  USE Functionals
  USE Thresholding
  USE HiCuThresholds
  USE AtomPairs
  USE SpecFun
  USE PrettyPrint
  USE CubeTree
  USE KxcGen
  IMPLICIT NONE
  !===============================================================================
  TYPE(DBL_RNK2)::LCoor,RCoor       ! Coordinates of the bounding boxes
  TYPE(DBL_RNK2)::RepLCoor,RepRCoor ! Coordinates of the repartitioned boxes
  INTEGER::NVol                     ! number of the bounding boxes
  REAL(DOUBLE)::MyLeavesTm          ! time to create all leaves in a subvolume
  REAL(DOUBLE)::BoxPoint(6)
  REAL(DOUBLE)::TotExc

  TYPE VaryLengLs
     INTEGER,POINTER,DIMENSION(:) :: Ls
  END TYPE VaryLengLs
  INTEGER :: NumOfDist
  INTEGER,ALLOCATABLE :: NewQdex(:),NewCdex(:),NewLdex(:)
  REAL(DOUBLE),ALLOCATABLE :: NewRList(:),NewExt(:),NewZeta(:),NewX(:),&
       NewY(:),NewZ(:),NewCo(:)
  TYPE(BBox) :: GRhoBBox,LocalRhoBBox
  REAL(DOUBLE) :: GRhoBBoxVol
  !===============================================================================
CONTAINS
  !===============================================================================

  !===============================================================================
  SUBROUTINE ParaInitRho(Args)
    TYPE(ARGMT)               :: Args
    INTEGER :: IErr,I

    CALL InitRho(Args)
    CALL GetLocalBoundingBox()
    !! now get the Global bounding box
    CALL MPI_AllReduce(LocalRhoBBox%BndBox(1,1),GRhoBBox%BndBox(1,1),&
         3,MPI_DOUBLE_PRECISION,MPI_MIN,MONDO_COMM,IErr)
    CALL MPI_AllReduce(LocalRhoBBox%BndBox(1,2),GRhoBBox%BndBox(1,2),&
         3,MPI_DOUBLE_PRECISION,MPI_MAX,MONDO_COMM,IErr)
    IF(MyID == 0) THEN
       CALL MakeBoxPeriodic(GRhoBBox)
    ENDIF
    CALL MPI_Bcast(GRhoBBox%BndBox(1,1),3,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
    CALL MPI_Bcast(GRhoBBox%BndBox(1,2),3,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
    GRhoBBoxVol = 1.0D0
    DO I = 1, 3
       GRhoBBoxVol = GRhoBBoxVol*(GRhoBBox%BndBox(I,2)-GRhoBBox%BndBox(I,1))
    ENDDO
    ! Messy output is unacceptable and uneccesary:  Protect with ifdef if you want to keep it
    !
    !    IF(MyID == 0) THEN
    !      CALL OpenASCII(OutFile,Out)
    !      WRITE(Out,*) 'MyID = 0, GRhoBBox = ',GRhoBBox%BndBox(1:3,1),GRhoBBox%BndBox(1:3,2), ', GRhoBBoxVol = ',GRhoBBoxVol
    !      CLOSE(Out,STATUS='KEEP')
    !ENDIF
    CALL AlignNodes()
  END SUBROUTINE ParaInitRho

  RECURSIVE FUNCTION CountRhoLeafNode(RhoTreeNode) RESULT(LeafNode)
    TYPE(RhoNode), POINTER :: RhoTreeNode
    INTEGER                :: LeafNode

    LeafNode = 0
    IF(RhoTreeNode%Leaf) THEN
       LeafNode = 1
    ELSE
       LeafNode = CountRhoLeafNode(RhoTreeNode%Descend) + &
                  CountRhoLeafNode(RhoTreeNode%Descend%Travrse)
    ENDIF
  END FUNCTION CountRhoLeafNode

  SUBROUTINE GetBBox()
    INTEGER        :: Power2(0:31),SmallN,PIndex,CIndex,Stage,DirInt,I,J,Ind,LineLocExist
    REAL(DOUBLE)   :: x2,LinDim
    TYPE(INT_VECT) :: ETDirArr
    TYPE(DBL_VECT) :: ETRootArr

    NVol = NPrc
    CALL New(LCoor,(/3,NVol/))
    CALL New(RCoor,(/3,NVol/))

    CALL Get(LineLocExist,'LineLocExist')
    IF(LineLocExist < 0) THEN
      STOP 'ERR: LineLocExist must be non-negative!'
    ELSEIF(LineLocExist > 0) THEN
       CALL New(ETDirArr,NPrc-1)
       CALL New(ETRootArr,NPrc-1)
       CALL Get(ETDirArr,'ETDirArr')
       CALL Get(ETRootArr,'ETRootArr')
       IF(MyID == 0) THEN
          NVol = LineLocExist
          IF(NVol /= NPrc) THEN
             WRITE(*,*) 'NVol = ',NVol, ',  NPrc = ',NPrc
             STOP 'ERROR: LineLoc.dat -- NVol is not equal to NPrc!'
          ENDIF
          LCoor%D(1:3,1) = GRhoBBox%BndBox(1:3,1)
          RCoor%D(1:3,1) = GRhoBBox%BndBox(1:3,2)

          Power2(0) = 1
          DO I = 1, 30 ! 31--this is an overflow on tru64
             Power2(I) = Power2(I-1)*2
          ENDDO
          SmallN = NInt( Log(NVol*1.0D0)/Log(2.0D0) )
          IF(NVol /= Power2(SmallN) ) THEN
             WRITE(*,*) 'NVol is ', NVol
             WRITE(*,*) 'NVol is not a power of 2!'
             STOP
          ENDIF

          Ind = 1
          DO Stage = 1, SmallN
             DO PIndex = 1, Power2(Stage-1)
                CIndex = PIndex + Power2(Stage-1)
                LCoor%D(1:3,CIndex) = LCoor%D(1:3,PIndex)
                RCoor%D(1:3,CIndex) = RCoor%D(1:3,PIndex)

                ! READ(54,*) DirInt
                ! READ(54,*) x2
                DirInt = ETDirArr%I(Ind)
                x2 = ETRootArr%D(Ind)
                Ind = Ind + 1
                RCoor%D(DirInt,PIndex) = x2
                LCoor%D(DirInt,CIndex) = x2
             ENDDO
          ENDDO

          !! check the volumes' dimensions are positive.
          DO I = 1, NVol
             DO J = 1, 3
                LinDim = RCoor%D(J,I)-LCoor%D(J,I)
                IF(LinDim <= 0) THEN
                   WRITE(*,*) 'ERROR: Volume assignment. LinDim is not positive!'
                   STOP
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       CALL Delete(ETDirArr)
       CALL Delete(ETRootArr)
    ELSE
       ! assign the first box
       IF(MyID == 0) THEN
          CALL OpenASCII(OutFile,Out)
          WRITE(Out,*) 'LineLoc is not found.'
          CLOSE(Out,STATUS='KEEP')

          ! use the GRhoBBox information rather than RhoRoot
          ! because RhoRoot is not assigned yet
          LCoor%D(1:3,1) = GRhoBBox%BndBox(1:3,1)
          RCoor%D(1:3,1) = GRhoBBox%BndBox(1:3,2)
          CALL Get_SubVol(LCoor,RCoor,NVol)
       ENDIF
    ENDIF

  END SUBROUTINE GetBBox

  SUBROUTINE GetNewLocalRhoBoundingBox()
    INTEGER :: I,KQ,J
    TYPE(BBox) :: NodeBox
    DO I = 1, NumOfDist
       KQ = NewQdex(I)
       NodeBox%BndBox(1:3,1) = (/NewX(KQ),NewY(KQ),NewZ(KQ)/)
       NodeBox%BndBox(1:3,2) = (/NewX(KQ),NewY(KQ),NewZ(KQ)/)
       NodeBox=ExpandBox(NodeBox,NewExt(KQ))
       IF(I == 1) THEN
          LocalRhoBBox%BndBox(1:3,1:2) = NodeBox%BndBox(1:3,1:2)
       ELSE
          DO J = 1, 3
             LocalRhoBBox%BndBox(J,1) = Min(LocalRhoBBox%BndBox(J,1),NodeBox%BndBox(J,1))
             LocalRhoBBox%BndBox(J,2) = Max(LocalRhoBBox%BndBox(J,2),NodeBox%BndBox(J,2))
          ENDDO
       ENDIF
    ENDDO
    LocalRhoBBox%Half(1:3) = (LocalRhoBBox%BndBox(1:3,2)-LocalRhoBBox%BndBox(1:3,1))*Half
    LocalRhoBBox%Center(1:3) = (LocalRhoBBox%BndBox(1:3,2)+LocalRhoBBox%BndBox(1:3,1))*Half
    CALL AlignNodes()
    !! WRITE(*,*) 'Finding out what new rho bounding box is : MyID = ', MyID, 'LocalRhoBBox = ',LocalRhoBBox%BndBox(1:3,1),LocalRhoBBox%BndBox(1:3,2)
    CALL AlignNodes()

  END SUBROUTINE GetNewLocalRhoBoundingBox

  SUBROUTINE GetLocalBoundingBox()
    INTEGER :: NDist,I,KQ,J
    TYPE(BBox) :: NodeBox

    NDist = Rho%NDist
    DO I = 1, NDist !! run through all distributions
       KQ = QDex(I) !! KQ is the index of a distribution
       NodeBox%BndBox(1:3,1) = (/Rho%Qx%D(KQ),Rho%Qy%D(KQ),Rho%Qz%D(KQ)/)
       NodeBox%BndBox(1:3,2) = (/Rho%Qx%D(KQ),Rho%Qy%D(KQ),Rho%Qz%D(KQ)/)
       NodeBox=ExpandBox(NodeBox,Ext(KQ))
       IF(I == 1) THEN
          LocalRhoBBox%BndBox(1:3,1:2) = NodeBox%BndBox(1:3,1:2)
       ELSE
          DO J = 1, 3
             LocalRhoBBox%BndBox(J,1) = Min(LocalRhoBBox%BndBox(J,1),NodeBox%BndBox(J,1))
             LocalRhoBBox%BndBox(J,2) = Max(LocalRhoBBox%BndBox(J,2),NodeBox%BndBox(J,2))
          ENDDO
       ENDIF
    ENDDO
    LocalRhoBBox%Half(1:3) = (LocalRhoBBox%BndBox(1:3,2)-LocalRhoBBox%BndBox(1:3,1))*Half
    LocalRhoBBox%Center(1:3) = (LocalRhoBBox%BndBox(1:3,2)+LocalRhoBBox%BndBox(1:3,1))*Half
    !! CALL AlignNodes()
    !! WRITE(*,*) 'MyID = ', MyID, 'LocalRhoBBox = ',LocalRhoBBox%BndBox(1:3,1),LocalRhoBBox%BndBox(1:3,2)
    CALL AlignNodes()
  END SUBROUTINE GetLocalBoundingBox

  SUBROUTINE SendBBox()
    INTEGER::I,IErr,NumDbl
    REAL(DOUBLE)::DblArr(3)
    INTEGER,DIMENSION(MPI_STATUS_SIZE)::Status

    NumDbl = 3*NPrc
    CALL MPI_BCast(LCoor%D(1,1),NumDbl,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
    CALL MPI_BCast(RCoor%D(1,1),NumDbl,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)

  END SUBROUTINE SendBBox

  SUBROUTINE DistDist()
    IMPLICIT NONE
    INTEGER :: NC,Tot,FloatDblIndex,ActIntRecvAmt,ActDblRecvAmt,&
         DblSendBufIndex,IntSendBufIndex,IntAmtRecv,DblAmtRecv,&
         SendTo,RecvFrom,CD,DistIndex,EllFillIndex,CoFillIndex,NumOfCo,&
         GIntSendMaxSize,GDblSendMaxSize,IntSendToOthers,DblSendToOthers,&
         I,KQ,NDist,J,Ell,LMNLen,Amt,LocalIntSendMaxSize,LocalDblSendMaxSize,IErr
    TYPE(BBox) :: NodeBox
    REAL(DOUBLE) :: BoxL(3),BoxU(3),QBL(3),QBU(3)
    TYPE(INT_VECT) :: CoToSend,IntToSend,IndexToSend,DblToSend,&
         TotIntSendFromProc,TotDblSendFromProc
    TYPE(VaryLengLs),ALLOCATABLE :: HeadArr(:)
    INTEGER,ALLOCATABLE :: IntSendBuf(:),IntRecvBuf(:),SF(:),Dest(:)
    REAL(DOUBLE),ALLOCATABLE :: DblSendBuf(:),DblRecvBuf(:)
    INTEGER,DIMENSION(MPI_STATUS_SIZE) :: IntStatus
    INTEGER,DIMENSION(MPI_STATUS_SIZE) :: DblStatus
    REAL(DOUBLE) :: StartTm,EndTm,TotTm
    REAL(DOUBLE),EXTERNAL    :: MondoTimer
    INTEGER::iSDen,NCoef,CD1

write(*,*)'Enter DistDist',MyID
call mpi_barrier(mondo_comm,ierr)
call mpi_barrier(mondo_comm,ierr)
call mpi_barrier(mondo_comm,ierr)

    StartTm = MondoTimer()
    !    IF(MyID == 0) THEN
    !      CALL OpenASCII(OutFile,Out)
    !      WRITE(Out,*) 'DistDist is entered...'
    !      CLOSE(Out,STATUS='KEEP')
    !    ENDIF

    IF(MyID == 0) THEN
       !! WRITE(*,*) 'Cube Box info:'
       DO I = 1, NPrc
          !! WRITE(*,*) 'Proc = ', I, ', LCoor%D(1:3,I) = ',LCoor%D(1:3,I),', RCoor%D(1:3,I) = ', RCoor%D(1:3,I)
       ENDDO
    ENDIF

    CALL New(IndexToSend,NPrc,0)
    CALL New(CoToSend,NPrc,0)
    CALL New(IntToSend,NPrc,0)
    CALL New(DblToSend,NPrc,0)

    IndexToSend%I(0:NPrc-1) = 0
    CoToSend%I(0:NPrc-1) = 0
    IntToSend%I(0:NPrc-1) = 0
    DblToSend%I(0:NPrc-1) = 0

    NDist = Rho%NDist
    NCoef = Rho%NCoef      ! <<< SPIN
    NSDen = Rho%NSDen      ! <<< SPIN

    write(*,*) 'DistDist NSDen',NSDen,' NCoef',NCoef,'NDist',NDist,MyID

    DO I = 1, NDist
       KQ = QDex(I)
       NodeBox%BndBox(1:3,1) = (/Rho%Qx%D(KQ),Rho%Qy%D(KQ),Rho%Qz%D(KQ)/)
       NodeBox%BndBox(1:3,2) = (/Rho%Qx%D(KQ),Rho%Qy%D(KQ),Rho%Qz%D(KQ)/)
       NodeBox=ExpandBox(NodeBox,Ext(KQ))
       !! Distribution box obtained
       !! WRITE(*,*) 'Distribution with index KQ = ',KQ, ' is to be inserted.'
       DO J = 0, NPrc-1
          DO NC = 1, CS_OUT%NCells
             BoxL(1:3) = LCoor%D(1:3,J+1)+CS_OUT%CellCarts%D(:,NC)
             BoxU(1:3) = RCoor%D(1:3,J+1)+CS_OUT%CellCarts%D(:,NC)
             QBL(1:3) = NodeBox%BndBox(1:3,1)
             QBU(1:3) = NodeBox%BndBox(1:3,2)
             IF(QBU(1) <= BoxL(1) .OR. QBL(1) >= BoxU(1) .OR. &
                  QBU(2) <= BoxL(2) .OR. QBL(2) >= BoxU(2) .OR. &
                  QBU(3) <= BoxL(3) .OR. QBL(3) >= BoxU(3)) THEN
                !! do nothing, no intersection
             ELSE
                !! Put KQ into the list
                !! first go the pointer that has is pointing to J
                !! WRITE(*,*) 'Distribution with index KQ = ',KQ, ' is to be inserted.'
                !! J is the proc to send
                IndexToSend%I(J) = IndexToSend%I(J) + 1
                IntToSend%I(J) = IntToSend%I(J) + 1 !! Ell
                Ell = Ldex(KQ)
                LMNLen = LHGTF(Ell)
                CoToSend%I(J) = CoToSend%I(J) + LMNLen*NSDen           ! <<< SPIN
                DblToSend%I(J) = DblToSend%I(J) + 5 + LMNLen*NSDen     ! <<< SPIN !! x,y,z,zeta,extent,co
                EXIT ! crucial exit
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !! allocate memory
    ALLOCATE(HeadArr(0:NPrc-1))
    CALL AlignNodes()
    DO I = 0, NPrc-1
       Amt = IndexToSend%I(I) + 1
       ALLOCATE(HeadArr(I)%LS(Amt))
    ENDDO
    CALL AlignNodes()
    IndexToSend%I(0:NPrc-1) = 0
    DO I = 1, NDist
       KQ = QDex(I)
       NodeBox%BndBox(1:3,1) = (/Rho%Qx%D(KQ),Rho%Qy%D(KQ),Rho%Qz%D(KQ)/)
       NodeBox%BndBox(1:3,2) = (/Rho%Qx%D(KQ),Rho%Qy%D(KQ),Rho%Qz%D(KQ)/)
       NodeBox=ExpandBox(NodeBox,Ext(KQ))
       !! Distribution box obtained
       DO J = 0, NPrc-1
          DO NC = 1, CS_OUT%NCells
             BoxL(1:3) = LCoor%D(1:3,J+1)+CS_OUT%CellCarts%D(:,NC)
             BoxU(1:3) = RCoor%D(1:3,J+1)+CS_OUT%CellCarts%D(:,NC)
             QBL(1:3) = NodeBox%BndBox(1:3,1)
             QBU(1:3) = NodeBox%BndBox(1:3,2)
             IF(QBU(1) <= BoxL(1) .OR. QBL(1) >= BoxU(1) .OR. &
                  QBU(2) <= BoxL(2) .OR. QBL(2) >= BoxU(2) .OR. &
                  QBU(3) <= BoxL(3) .OR. QBL(3) >= BoxU(3)) THEN
                !! do nothing, no intersection
             ELSE
                !! Put KQ into the list
                !! first go the pointer that has is pointing to J
                !! WRITE(*,*) 'Distribution with index KQ = ',KQ, ' is to be inserted.'
                !! J is the proc to send
                IndexToSend%I(J) = IndexToSend%I(J) + 1
                HeadArr(J)%LS(IndexToSend%I(J)) = I
                EXIT ! crucial exit
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    CALL AlignNodes()
    IF(MyID == 0) THEN
       !! WRITE(*,*) 'MyID = ',MyID, ' done with storing indices'
    ENDIF
    CALL AlignNodes()

    !! LocalIntSendMaxSize and LocalDblSendMaxSize determine the sizes of the send buffers
    LocalIntSendMaxSize = 0
    LocalDblSendMaxSize = 0
    DO I = 0, NPrc-1
       IF(I == MyID) Cycle
       LocalIntSendMaxSize = Max(LocalIntSendMaxSize,IntToSend%I(I)+1) !! 1 for dummy
       LocalDblSendMaxSize = Max(LocalDblSendMaxSize,DblToSend%I(I)+1) !! 1 for dummy
    ENDDO
    CALL AlignNodes()
    !! WRITE(*,*) 'MyID = ',MyID, ', LocalIntSendMaxSize=',LocalIntSendMaxSize,' ,LocalDblSendMaxSize = ',LocalDblSendMaxSize
    CALL AlignNodes()
    !! GIntSendMaxSize and GDblSendMaxSize determine the sizes of the receive buffers. Of course
    !! this may be smaller on some processors, but for simplicity we assume it is the same for all.
    CALL MPI_AllReduce(LocalIntSendMaxSize,GIntSendMaxSize,1,MPI_INTEGER,MPI_MAX,MONDO_COMM,IErr)
    CALL MPI_AllReduce(LocalDblSendMaxSize,GDblSendMaxSize,1,MPI_INTEGER,MPI_MAX,MONDO_COMM,IErr)
    !!
    !    IF(MyID == 0) THEN
    !      CALL OpenASCII(OutFile,Out)
    !      WRITE(Out,*) 'MyID=0, GIntSendMaxSize=',GIntSendMaxSize,' ,GDblSendMaxSize = ',GDblSendMaxSize
    !      CLOSE(Out,STATUS='KEEP')
    !    ENDIF
    IntSendToOthers = 0
    DblSendToOthers = 0
    DO I = 0, NPrc-1
       IF(I == MyID) Cycle
       IntSendToOthers = IntSendToOthers + IntToSend%I(I)
       DblSendToOthers = DblSendToOthers + DblToSend%I(I)
    ENDDO
    !! WRITE(*,*) 'MYID = ',MyID,', IntSendToOthers = ',IntSendToOthers,', DblSendToOthers = ',DblSendToOthers
    CALL New(TotIntSendFromProc,NPrc,0)
    CALL New(TotDblSendFromProc,NPrc,0)
    CALL MPI_GATHER(IntSendToOthers,1,MPI_INTEGER,TotIntSendFromProc%I(0),1,MPI_INTEGER,0,MONDO_COMM,IErr)
    CALL MPI_GATHER(DblSendToOthers,1,MPI_INTEGER,TotDblSendFromProc%I(0),1,MPI_INTEGER,0,MONDO_COMM,IErr)
    IF(MyID == 0) THEN
       !! WRITE(*,*) 'Gather on MyID = 0, TotIntSendFromProc%I(0:NPrc-1) = ',TotIntSendFromProc%I(0:NPrc-1)
       !! WRITE(*,*) 'Gather on MyID = 0, TotDblSendFromProc%I(0:NPrc-1) = ',TotDblSendFromProc%I(0:NPrc-1)
       Tot = 0
       DO I = 0, NPrc-1
          Tot = Tot + TotIntSendFromProc%I(I)
       ENDDO
       !      CALL OpenASCII(OutFile,Out)
       !      WRITE(Out,*) 'DistDist: Gather on ROOT, Tot Int to pass around is ',Tot
       Tot = 0
       DO I = 0, NPrc-1
          Tot = Tot + TotDblSendFromProc%I(I)
       ENDDO
       !      WRITE(Out,*) 'DistDist: Gather on ROOT, Tot Dbl to pass around is ',Tot
       !      CLOSE(Out,STATUS='KEEP')
    ENDIF

    !! allocate memory enough for packing and sending
    ALLOCATE(IntSendBuf(LocalIntSendMaxSize),DblSendBuf(LocalDblSendMaxSize))
    ALLOCATE(IntRecvBuf(GIntSendMaxSize),DblRecvBuf(GDblSendMaxSize))

    !! allocate all arrays for unpacking.
    DO I = 0, NPrc-1
       CALL MPI_Reduce(IndexToSend%I(I),NumOfDist,1,MPI_INTEGER,MPI_SUM,I,MONDO_COMM,IErr)
       CALL MPI_Reduce(CoToSend%I(I),NumOfCo,1,MPI_INTEGER,MPI_SUM,I,MONDO_COMM,IErr)
    ENDDO
    !! WRITE(*,*) 'MyID = ', MyID, ', NumOfDist = ',NumOfDist, ', NumOfCo = ',NumOfCo
    ALLOCATE(NewQdex(NumOfDist))
    ALLOCATE(NewCdex(NumOfDist))
    ALLOCATE(NewLdex(NumOfDist))

    ALLOCATE(NewExt(NumOfDist))
    ALLOCATE(NewZeta(NumOfDist))
    ALLOCATE(NewX(NumOfDist))
    ALLOCATE(NewY(NumOfDist))
    ALLOCATE(NewZ(NumOfDist))
    ALLOCATE(NewCo(NumOfCo))

    !! Fill in with the local data
    EllFillIndex = 0
    CoFillIndex = 0
    DO I = 1, IndexToSend%I(MyID)
       EllFillIndex = EllFillIndex + 1
       DistIndex = HeadArr(MyID)%Ls(EllFillIndex)
       KQ = Qdex(DistIndex) !! absolute distribution index in Rho
       Ell = Ldex(KQ)
       NewLdex(EllFillIndex) = Ell
       NewX(EllFillIndex) = Rho%Qx%D(KQ)
       NewY(EllFillIndex) = Rho%Qy%D(KQ)
       NewZ(EllFillIndex) = Rho%Qz%D(KQ)
       NewExt(EllFillIndex) = Ext(KQ)
       NewZeta(EllFillIndex) = Zeta(KQ)
       LMNLen = LHGTF(Ell)
       CD = Cdex(KQ)
       NewCdex(EllFillIndex) = CoFillIndex+1
       DO iSDen=1,NSDen                                                        ! <<< SPIN
          CD1=CD+(iSDen-1)*NCoef                                               ! <<< SPIN
          NewCo(CoFillIndex+1:CoFillIndex+LMNLen) = Rho%Co%D(CD1:CD1+LMNLen-1) ! <<< SPIN
          CoFillIndex = CoFillIndex + LMNLen                                   ! <<< SPIN
       ENDDO                                                                   ! <<< SPIN
    ENDDO
    !! WRITE(*,*) 'MyID = ',MyID, ', EllFillIndex = ',EllFillIndex,', CoFillIndex = ',CoFillIndex, ', IndexToSend%I(MyID) = ', IndexToSend%I(MyID), ', CoToSend%I(MyID) = ', CoToSend%I(MyID)

    IntAmtRecv = GIntSendMaxSize
    DblAmtRecv = GDblSendMaxSize

    ALLOCATE(SF(0:NPrc-1),Dest(0:NPrc-1))
    DO I = 1, NPrc-1
       DO J = 0, NPrc-1
          SendTo = Modulo(J+I,NPrc)
          Dest(J) = SendTo
          SF(J) = 0
       ENDDO
       DO J = 0, NPrc-1
          IF(SF(J) == 0 .AND. SF(Dest(J)) == 0) THEN
             SF(J) = 1
             SF(Dest(J)) = 2
          ENDIF
       ENDDO
       SendTo = Modulo(MyID+I,NPrc)
       RecvFrom = Modulo(MyID-I,NPrc)
       IF(SF(MyID) == 1) THEN
          !! packing int
          IntSendBufIndex = 0
          DblSendBufIndex = 0
          DO J = 1, IndexToSend%I(SendTo)
             IntSendBufIndex = IntSendBufIndex + 1
             DistIndex = HeadArr(SendTo)%Ls(IntSendBufIndex)
             KQ = Qdex(DistIndex)
             Ell = Ldex(KQ)
             IntSendBuf(IntSendBufIndex) = Ell
             DblSendBuf(DblSendBufIndex+1:DblSendBufIndex+3) = (/Rho%Qx%D(KQ),Rho%Qy%D(KQ),Rho%Qz%D(KQ)/)
             DblSendBufIndex = DblSendBufIndex + 3 + 1
             DblSendBuf(DblSendBufIndex) = Zeta(KQ)
             DblSendBufIndex = DblSendBufIndex+1
             DblSendBuf(DblSendBufIndex) = Ext(KQ)
             LMNLen = LHGTF(Ell)
             CD = Cdex(KQ)
             DO iSDen=1,NSDen                                                                     ! <<< SPIN
                CD1=CD+(iSDen-1)*NCoef                                                            ! <<< SPIN
                DblSendBuf(DblSendBufIndex+1:DblSendBufIndex+LMNLen) = Rho%Co%D(CD1:CD1+LMNLen-1) ! <<< SPIN
                DblSendBufIndex = DblSendBufIndex + LMNLen                                        ! <<< SPIN
             ENDDO                                                                                ! <<< SPIN
          ENDDO
          IF(IntToSend%I(SendTo) /= IntSendBufIndex &
               .OR. DblToSend%I(SendTo) /= DblSendBufIndex) THEN
             WRITE(*,*) 'MyID = ',MyID, ', IntToSend%I(SendTo) = ',IntToSend%I(SendTo), ', IntSendBufIndex = ',IntSendBufIndex
             WRITE(*,*) 'DblToSend%I(SendTo) = ', DblToSend%I(SendTo),', DblSendBufIndex = ',DblSendBufIndex
             STOP 'ERR: The numbers should the same!'
          ENDIF
          IntSendBufIndex = IntSendBufIndex + 1
          DblSendBufIndex = DblSendBufIndex + 1
          IntSendBuf(IntSendBufIndex) = 100000
          DblSendBuf(DblSendBufIndex) = 100000.0D0
          CALL MPI_Send(IntSendBuf(1),IntSendBufIndex,MPI_INTEGER,SendTo,MyID,MONDO_COMM,IErr)
          CALL MPI_Send(DblSendBuf(1),DblSendBufIndex,MPI_DOUBLE_PRECISION,SendTo,MyID,MONDO_COMM,IErr)

          !! receiving
          CALL MPI_Recv(IntRecvBuf(1),GIntSendMaxSize,MPI_INTEGER,RecvFrom,RecvFrom,MONDO_COMM,IntStatus,IErr)
          CALL MPI_Get_Count(IntStatus,MPI_INTEGER,ActIntRecvAmt,IErr)
          CALL MPI_Recv(DblRecvBuf(1),GDblSendMaxSize,MPI_DOUBLE_PRECISION,RecvFrom,RecvFrom,MONDO_COMM,DblStatus,IErr)
          CALL MPI_Get_Count(DblStatus,MPI_DOUBLE_PRECISION,ActDblRecvAmt,IErr)
          !! unpacking.
          FloatDblIndex = 0
          DO J = 1, ActIntRecvAmt-1
             EllFillIndex = EllFillIndex+1
             Ell = IntRecvBuf(J)
             LMNLen = LHGTF(Ell)
             NewLdex(EllFillIndex) = Ell
             FloatDblIndex = FloatDblIndex+1
             NewX(EllFillIndex) = DblRecvBuf(FloatDblIndex)
             FloatDblIndex = FloatDblIndex+1
             NewY(EllFillIndex) = DblRecvBuf(FloatDblIndex)
             FloatDblIndex = FloatDblIndex+1
             NewZ(EllFillIndex) = DblRecvBuf(FloatDblIndex)
             FloatDblIndex = FloatDblIndex+1
             NewZeta(EllFillIndex) = DblRecvBuf(FloatDblIndex)
             FloatDblIndex = FloatDblIndex+1
             NewExt(EllFillIndex) = DblRecvBuf(FloatDblIndex)
             NewCdex(EllFillIndex) = CoFillIndex+1
             DO iSDen=1,NSDen                                                                   ! <<< SPIN
                NewCo(CoFillIndex+1:CoFillIndex+LMNLen) = DblRecvBuf(FloatDblIndex+1:FloatDblIndex+LMNLen)
                FloatDblIndex = FloatDblIndex + LMNLen                                          ! <<< SPIN
                CoFillIndex = CoFillIndex + LMNLen                                              ! <<< SPIN
             ENDDO                                                                              ! <<< SPIN
          ENDDO

       ELSE
          CALL MPI_Recv(IntRecvBuf(1),GIntSendMaxSize,MPI_INTEGER,RecvFrom,RecvFrom,MONDO_COMM,IntStatus,IErr)
          CALL MPI_Get_Count(IntStatus,MPI_INTEGER,ActIntRecvAmt,IErr)
          CALL MPI_Recv(DblRecvBuf(1),GDblSendMaxSize,MPI_DOUBLE_PRECISION,RecvFrom,RecvFrom,MONDO_COMM,DblStatus,IErr)
          CALL MPI_Get_Count(DblStatus,MPI_DOUBLE_PRECISION,ActDblRecvAmt,IErr)
          !! unpacking.
          FloatDblIndex = 0
          DO J = 1, ActIntRecvAmt-1
             EllFillIndex = EllFillIndex+1
             Ell = IntRecvBuf(J)
             LMNLen = LHGTF(Ell)
             NewLdex(EllFillIndex) = Ell
             FloatDblIndex = FloatDblIndex+1
             NewX(EllFillIndex) = DblRecvBuf(FloatDblIndex)
             FloatDblIndex = FloatDblIndex+1
             NewY(EllFillIndex) = DblRecvBuf(FloatDblIndex)
             FloatDblIndex = FloatDblIndex+1
             NewZ(EllFillIndex) = DblRecvBuf(FloatDblIndex)
             FloatDblIndex = FloatDblIndex+1
             NewZeta(EllFillIndex) = DblRecvBuf(FloatDblIndex)
             FloatDblIndex = FloatDblIndex+1
             NewExt(EllFillIndex) = DblRecvBuf(FloatDblIndex)
             NewCdex(EllFillIndex) = CoFillIndex+1
             DO iSDen=1,NSDen                                                                   ! <<< SPIN
                NewCo(CoFillIndex+1:CoFillIndex+LMNLen) = DblRecvBuf(FloatDblIndex+1:FloatDblIndex+LMNLen)
                FloatDblIndex = FloatDblIndex + LMNLen                                          ! <<< SPIN
                CoFillIndex = CoFillIndex + LMNLen                                              ! <<< SPIN
             ENDDO                                                                              ! <<< SPIN

          ENDDO

          !! packing int
          IntSendBufIndex = 0
          DblSendBufIndex = 0
          DO J = 1, IndexToSend%I(SendTo)
             IntSendBufIndex = IntSendBufIndex + 1
             DistIndex = HeadArr(SendTo)%Ls(IntSendBufIndex)
             KQ = Qdex(DistIndex)
             Ell = Ldex(KQ)
             IntSendBuf(IntSendBufIndex) = Ell
             DblSendBuf(DblSendBufIndex+1:DblSendBufIndex+3) = (/Rho%Qx%D(KQ),Rho%Qy%D(KQ),Rho%Qz%D(KQ)/)
             DblSendBufIndex = DblSendBufIndex + 4
             DblSendBuf(DblSendBufIndex) = Zeta(KQ)
             DblSendBufIndex = DblSendBufIndex+1
             DblSendBuf(DblSendBufIndex) = Ext(KQ)
             LMNLen = LHGTF(Ell)
             CD = Cdex(KQ)
             DO iSDen=1,NSDen                                                       ! <<< SPIN
                CD1=CD+(iSDen-1)*NCoef                                              ! <<< SPIN
                DblSendBuf(DblSendBufIndex+1:DblSendBufIndex+LMNLen) = Rho%Co%D(CD1:CD1+LMNLen-1)
                DblSendBufIndex = DblSendBufIndex + LMNLen                          ! <<< SPIN
             ENDDO                                                                  ! <<< SPIN

          ENDDO
          IF(IntToSend%I(SendTo) /= IntSendBufIndex &
               .OR. DblToSend%I(SendTo) /= DblSendBufIndex) THEN
             WRITE(*,*) 'MyID = ',MyID, ', IntToSend%I(SendTo) = ',IntToSend%I(SendTo), ', IntSendBufIndex = ',IntSendBufIndex
             WRITE(*,*) 'DblToSend%I(SendTo) = ', DblToSend%I(SendTo),', DblSendBufIndex = ',DblSendBufIndex
             STOP 'ERR: The numbers should the same!'
          ENDIF
          IntSendBufIndex = IntSendBufIndex + 1
          DblSendBufIndex = DblSendBufIndex + 1
          IntSendBuf(IntSendBufIndex) = 100000
          DblSendBuf(DblSendBufIndex) = 100000.0D0
          CALL MPI_Send(IntSendBuf(1),IntSendBufIndex,MPI_INTEGER,SendTo,MyID,MONDO_COMM,IErr)
          CALL MPI_Send(DblSendBuf(1),DblSendBufIndex,MPI_DOUBLE_PRECISION,SendTo,MyID,MONDO_COMM,IErr)

       ENDIF
    ENDDO
    DO I = 1, NumOfDist
       NewQdex(I) = I
    ENDDO

    !! WRITE(*,*) 'MyID = ',MyID, ' NumOfDist = ', NumOfDist, ', EllFillIndex = ',EllFillIndex, ', NumOfCo = ', NumOfCo,', CoFillIndex = ',CoFillIndex
    IF(NumOfDist /= EllFillIndex) STOP 'ERR: NumOfDist /= EllFillIndex!'
    IF(NumOfCo /= CoFillIndex) STOP 'ERR: NumOfCo /= CoFillIndex!'
    CALL AlignNodes()
    !    IF(MyID == 0) THEN
    !      CALL OpenASCII(OutFile,Out)
    !      WRITE(Out,*) 'MyID=0, end of DistDist.'
    !      CLOSE(Out,STATUS='KEEP')
    !    ENDIF
    CALL Delete(IndexToSend)
    CALL Delete(CoToSend)
    CALL Delete(IntToSend)
    CALL Delete(DblToSend)
    CALL Delete(TotIntSendFromProc)
    CALL Delete(TotDblSendFromProc)
    DEALLOCATE(IntSendBuf,DblSendBuf)
    DEALLOCATE(IntRecvBuf,DblRecvBuf)
    DEALLOCATE(SF,Dest)
    DO I = 0, NPrc-1
       DEALLOCATE(HeadArr(I)%LS)
    ENDDO
    DEALLOCATE(HeadArr)
    CALL DeleteRho
    CALL AlignNodes()
    EndTm = MondoTimer()
    TotTm = EndTm - StartTm

write(*,*)'Exit DistDist',MyID
call mpi_barrier(mondo_comm,ierr)
call mpi_barrier(mondo_comm,ierr)
call mpi_barrier(mondo_comm,ierr)

    !    IF(MyID == ROOT) THEN
    !      CALL OpenASCII(OutFile,Out)
    !      WRITE(Out,*) 'Total time to do DistDist is ', TotTm
    !      CLOSE(Out,STATUS='KEEP')
    !    ENDIF
  END SUBROUTINE DistDist

  SUBROUTINE ParaRhoToTree()
    INTEGER :: I

    RhoNodes = 0
    RhoLevel = 0
    CALL NewRhoNode(RhoRoot,0)
    CALL ParaInitRhoRoot()
    ALLOCATE(NewRList(NumOfDist))
    CALL ParaSplitRho(RhoRoot)

    CALL GetNewLocalRhoBoundingBox()

    DEALLOCATE(NewX,NewY,NewZ,NewLdex,NewQdex,NewExt,NewZeta,NewCo,NewCdex,NewRList)
    DO CurrentTier = RhoLevel,0,-1
       CALL MergeRhoBBox(RhoRoot)
    ENDDO
    !! WRITE(*,*) 'After merging: MyID = ',MyID, 'RhoRoot%Box%BndBox(1:3,1) = ',RhoRoot%Box%BndBox(1:3,1), ' RhoRoot%Box%BndBox(1:3,2) = ', RhoRoot%Box%BndBox(1:3,2)
    DO I = 1, 3
       IF(ABS(RhoRoot%Box%BndBox(I,1)-LocalRhoBBox%BndBox(I,1)) > 1.0D-10) THEN
          WRITE(*,*) 'L: Rho Bounding Box dimension is not consistent!'
          STOP
       ENDIF
       IF(ABS(RhoRoot%Box%BndBox(I,2)-LocalRhoBBox%BndBox(I,2)) > 1.0D-10) THEN
          WRITE(*,*) 'R: Rho Bounding Box dimension is not consistent!'
          STOP
       ENDIF
    ENDDO

  END SUBROUTINE ParaRhoToTree

  !===================================================================
  !     Recursively partition the density into a 4-D BinTree
  !===================================================================
  RECURSIVE SUBROUTINE ParaSplitRho(Node)
    TYPE(RhoNode), POINTER :: Node,Left,Right
    IF(Node%NQ==1)THEN
       CALL ParaFillRhoLeaf(Node)
    ELSE
       ! Allocate new children
       CALL NewRhoNode(Node%Descend,Node%Box%Tier+1)
       CALL NewRhoNode(Node%Descend%Travrse,Node%Box%Tier+1)
       ! Set links
       Node%Descend%Travrse%Travrse=>Node%Travrse
       Left=>Node%Descend
       Right=>Node%Descend%Travrse
       CALL ParaSplitRhoBox(Node,Left,Right)
       ! Recur
       CALL ParaSplitRho(Left)
       CALL ParaSplitRho(Right)
    ENDIF
  END SUBROUTINE ParaSplitRho

  ! Orthogonal Recusive Bisection
  SUBROUTINE ParaSplitRhoBox(Node,Left,Right)
    TYPE(RhoNode), POINTER :: Node,Left,Right
    REAL(DOUBLE)           :: Section
    INTEGER                :: B,E,N,ISplit,Split,I,J,k
    ! Indexing
    J=0
    B=Node%Bdex
    E=Node%Edex
    N=E-B+1
    ! Choose direction to section
    ISplit=Mod(Node%Box%Tier,4)+1
    IF(ISplit==1)THEN
       ! Split X dir
       DO I=B,E
          J=J+1
          K=NewQdex(I)
          NewRList(J)=NewX(K)
       ENDDO
    ELSEIF(ISplit==2)THEN
       ! Split on Y dir
       DO I=B,E
          J=J+1
          K=NewQdex(I)
          NewRList(J)=NewY(K)
       ENDDO
    ELSEIF(ISplit==3)THEN
       ! Split on Z dir
       DO I=B,E
          J=J+1
          K=NewQdex(I)
          NewRList(J)=NewZ(K)
       ENDDO
    ELSEIF(ISplit==4)THEN
       ! Split on box size
       DO I=B,E
          J=J+1
          K=NewQdex(I)
          NewRList(J)=NewExt(K)
       ENDDO
    ENDIF
    ! Sort
    CALL DblIntSort77(N,NewRList,NewQdex(B:E),2)
    ! Orthogonal recursive bisection (ORB)
    Section   =NewRList(1)+Half*(NewRList(N)-NewRList(1))
    ! Split     =BinarySearch(N,NewRList,Section)
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
    ! Set counters
    Left%NQ=Left%Edex-Left%Bdex+1
    Right%NQ=Right%Edex-Right%Bdex+1
  END SUBROUTINE ParaSplitRhoBox

  !===============================================================================
  !     Fill leaf nodes with data
  !===============================================================================
  SUBROUTINE ParaFillRhoLeaf(Node)
    IMPLICIT NONE
    TYPE(RhoNode), POINTER :: Node
    INTEGER   :: I,IQ,IC,J,K,KQ,KC,L,B,E,N,NQ,NC,LMNLen,LTmp,Status
    REAL(DOUBLE) :: RhoSum
    ! Set leaf logical
    Node%Leaf=.True.
    ! Boundaries in the ordered lists
    B=Node%Bdex
    E=Node%Edex
    NQ=E-B+1
    IF(NQ/=1.OR.B/=E)  &
         CALL Halt('Bad Logic in FillRhoLeaf ')
    ! Filler up...
    KQ=NewQdex(B)
    KC=NewCdex(KQ)
    Node%Ell=NewLdex(KQ)
    Node%Zeta=NewZeta(KQ)
    Node%Extent=NewExt(KQ)
    Node%Qx=NewX(KQ)
    Node%Qy=NewY(KQ)
    Node%Qz=NewZ(KQ)
    ! Allocate HGTF coefficients array
    LMNLen=LHGTF(Node%Ell)
    ALLOCATE(Node%Co(1:LMNLen*NSDen),STAT=Status)   !<<< SPIN
    CALL IncMem(Status,0,LMNLen*NSDen)              !<<< SPIN
    ! Transfer data in
    Node%Co(1:LMNLen*NSDen)=NewCo(KC:KC+LMNLen*NSDen-1)  !<<< SPIN
  END SUBROUTINE ParaFillRhoLeaf

  SUBROUTINE ParaInitRhoRoot()
    RhoRoot%Bdex = 1
    RhoRoot%Edex = NumOfDist
    RhoRoot%NQ = NumOfDist
  END SUBROUTINE ParaInitRhoRoot

  SUBROUTINE ParaGridGen()
    TYPE(BBox) :: WBox
    REAL(DOUBLE)::TotRho,SubVolRho,SubVolExc

    WBox%BndBox(1:3,1) = LCoor%D(1:3,MyID+1)
    WBox%BndBox(1:3,2) = RCoor%D(1:3,MyID+1)

    CALL CalCenterAndHalf(WBox)

    !    IF(MyID == 0) THEN
    !      CALL OpenASCII(OutFile,Out)
    !      WRITE(Out,*) 'Calling GridGen...'
    !      CLOSE(Out,STATUS='KEEP')
    !    ENDIF
    CALL GridGen(WBox,SubVolRho,SubVolExc)
    TotRho = Reduce(SubVolRho)
    TotExc = Reduce(SubVolExc)

if(myid==0)write(*,*)'ParaGridGen: TotRho=',TotRho,' TotExc=',TotExc,MyID

    !    IF(MyID == ROOT) THEN
    !      CALL OpenASCII(OutFile,Out)
    !      WRITE(Out,*) 'ParaGridGen: TotRho = ',TotRho, ', TotExc = ',TotExc
    !      WRITE(Out,*) 'ParaGridGen is done.'
    !      CLOSE(Out,STATUS='KEEP')
    !    ENDIF

  END SUBROUTINE ParaGridGen

  SUBROUTINE WorkBBox(Kxc)
    TYPE(BBox) :: WBox
    REAL(DOUBLE)::TotRho,SubVolRho,SubVolExc
    REAL(DOUBLE)::HiCuBegTm,HiCuEndTm,HiCuTm
    TYPE(DBL_VECT) :: TmHiCuArr
    TYPE(FastMat),POINTER  :: Kxc
    INTEGER :: IErr
    REAL(DOUBLE),EXTERNAL    :: MondoTimer

    WBox%BndBox(1:3,1) = LCoor%D(1:3,MyID+1)
    WBox%BndBox(1:3,2) = RCoor%D(1:3,MyID+1)

    CALL CalCenterAndHalf(WBox)
    HiCuBegTm = MondoTimer()
    !    IF(MyID == 0) THEN
    !      CALL OpenASCII(OutFile,Out)
    !      WRITE(Out,*) 'Calling GridGen...'
    !      CLOSE(Out,STATUS='KEEP')
    !    ENDIF
    CALL GridGen(WBox,SubVolRho,SubVolExc)
    TotRho = Reduce(SubVolRho)
    TotExc = Reduce(SubVolExc)

if(myid==0)write(*,*)'WorkBBox: TotRho=',TotRho,' TotExc=',TotExc,MyID

    MyLeavesTm = LeavesTmCount(CubeRoot)
    !    IF(MyID == ROOT) THEN
    !      CALL OpenASCII(OutFile,Out)
    !      WRITE(Out,*) 'GridGen: TotRho = ',TotRho,', TotExc = ',TotExc
    !      WRITE(Out,*) 'GridGen is done. Calling MakeKxc...'
    !      CLOSE(Out,STATUS='KEEP')
    !    ENDIF
    CALL MakeKxc(Kxc,CubeRoot)
    HiCuEndTm = MondoTimer()
    ! VolTm = TmEndM-TmBegM
    HiCuTm = HiCuEndTm-HiCuBegTm
    CALL New(TmHiCuArr,NPrc)
    CALL MPI_Gather(HiCuTm,1,MPI_DOUBLE_PRECISION,TmHiCuArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
    !    IF(MyID == ROOT) THEN
    !      CALL PImbalance(TmHiCuArr,NPrc,Prog_O='HiCu')
    !    ENDIF
    CALL Delete(TmHiCuArr)

  END SUBROUTINE WorkBBox

  SUBROUTINE RepartitionVol()

    LOGICAL::Busy
    TYPE(DBL_VECT)::RepLeavesTm
    INTEGER::SmallN,I,J,Ierr,RootFindIter,Stage,PIndex,CIndex,DirInt
    REAL(DOUBLE)::DIFF,x0,x1,f0,f1,x2,f2,oldf2,LinDim,MaxDim,OldVol,ThisVol,&
         AbsSum,LocalLeavesTm,LeavesTmInside,OrigLeafTm,NewVol
    INTEGER::Power2(0:31)
    REAL(DOUBLE),PARAMETER::RootTau=1.0D-4
    REAL(DOUBLE)::StartTm,EndTm,TotTm,AllLeavesTm
    REAL(DOUBLE),PARAMETER::TauBeg=1.0D-5,TauEnd=1.0D-4
    REAL(DOUBLE) :: lnTauBeg,lnTauEnd,lnTauH,lnTau
    TYPE(DBL_VECT) :: TauArr,ETRootArr
    TYPE(INT_VECT) :: ETDirArr
    INTEGER :: RootNum,RootIndex
    REAL(DOUBLE),EXTERNAL    :: MondoTimer


    StartTm = MondoTimer()

    CALL New(RepLCoor,(/3,NVol/))
    CALL New(RepRCoor,(/3,NVol/))
    CALL New(RepLeavesTm,NVol)

    AllLeavesTm = Reduce(MyLeavesTm)

    Busy = .TRUE.
    IF(MyID == 0) THEN
       ! Messy output is unacceptable and uneccesary:  Protect with ifdef if you want to keep it
       !CALL OpenASCII(OutFile,Out)
       !WRITE(Out,*) 'RepartitionVol: RhoRoot Box  ',RhoRoot%Box%BndBox(1:3,1),RhoRoot%Box%BndBox(1:3,2)
       !WRITE(Out,*) 'RepartitionVol: GRhoBBox   Box  ',GRhoBBox%BndBox(1:3,1),GRhoBBox%BndBox(1:3,2)
       !CLOSE(Out,STATUS='KEEP')
       RepLCoor%D(1:3,1) = GRhoBBox%BndBox(1:3,1)
       RepRCoor%D(1:3,1) = GRhoBBox%BndBox(1:3,2)
       !! RepLeavesTm%D(1) = Sum(VolLeavesTm%D(:))
       RepLeavesTm%D(1) = AllLeavesTm
       NewVol = 1.0D0
       DO I = 1, 3
          NewVol = NewVol*(RepRCoor%D(I,1)-RepLCoor%D(I,1))
       ENDDO
       ! Messy output is unacceptable and uneccesary:  Protect with ifdef if you want to keep it
       !CALL OpenASCII(OutFile,Out)
       !WRITE(Out,*) 'The total leave times is ',RepLeavesTm%D(1)
       !WRITE(Out,*) 'The initial volume to start with is ',NewVol
       !CLOSE(Out,STATUS='KEEP')

       Power2(0) = 1
       DO I = 1, 30 ! 31--overflow
          Power2(I) = Power2(I-1)*2
       ENDDO
       SmallN = NINT(LOG(NVol*1.0D0)/LOG(2.0D0))
       IF(NVol /= Power2(SmallN)) THEN
          WRITE(*,*) 'NVol is ',NVol
          WRITE(*,*) 'NVol is not a power of 2!'
          STOP
       ENDIF

#undef REGULA_FALSI
#ifdef REGULA_FALSI
       ! Messy output is unacceptable and uneccesary:  Protect with ifdef if you want to keep it

       !      CALL OpenASCII(OutFile,Out)
       !      WRITE(Out,*) 'REGULA_FALSI is defined!'
       !      CLOSE(Out,STATUS='KEEP')
#else
       ! Messy output is unacceptable and uneccesary:  Protect with ifdef if you want to keep it
       !CALL OpenASCII(OutFile,Out)
       !WRITE(Out,*) 'REGULA_FALSI is not defined. Bisection search is used.'
       !WRITE(Out,*) 'TauBeg = ',TauBeg,' ,TauEnd = ', TauEnd
       !CLOSE(Out,STATUS='KEEP')
#endif

#ifdef REGULA_FALSI
#else
       CALL New(TauArr,SmallN)
       IF(SmallN == 1) THEN
          TauArr%D(1) = TauBeg
       ELSE
          lnTauBeg = Log(TauBeg)
          lnTauEnd = Log(TauEnd)
          lnTauH = (lnTauEnd - lnTauBeg)/(SmallN-1.0D0)
          DO I  = 1, SmallN
             lnTau = lnTauBeg + (I-1)*lnTauH
             TauArr%D(I) = EXP(lnTau)
             IF(I == 1 .OR. I == SmallN) THEN
                ! Messy output is unacceptable and uneccesary:  Protect with ifdef if you want to keep it
                !            CALL OpenASCII(OutFile,Out)
                !            WRITE(Out,*) 'TauArr%D(Stage=',I,')=',TauArr%D(I)
                !            CLOSE(Out,STATUS='KEEP')
             ENDIF
          ENDDO
       ENDIF
#endif

       RootNum = Power2(SmallN)-1
       CALL New(ETRootArr,RootNum)
       CALL New(ETDirArr,RootNum)
       RootIndex = 0
       DO Stage = 1, SmallN
          !! WRITE(*,*) 'Stage = ',Stage
          DO PIndex = 1, Power2(Stage-1)
             CIndex = PIndex + Power2(Stage-1)
             !! copy the coordinates to the child box
             RepLCoor%D(1:3,CIndex) = RepLCoor%D(1:3,PIndex)
             RepRCoor%D(1:3,CIndex) = RepRCoor%D(1:3,PIndex)

             !! direction to split
             MaxDim = -1.0D0
             DO I = 1, 3
                LinDim = RepRCoor%D(I,PIndex)-RepLCoor%D(I,PIndex)
                IF(LinDim <= 0) STOP 'ERROR: LinDim is not positive in RepartionVol'
                IF(LinDim > MaxDim) THEN
                   MaxDim = LinDim
                   DirInt = I
                ENDIF
             ENDDO
             ! WRITE(54,*) DirInt
             RootIndex = RootIndex + 1
             ETDirArr%I(RootIndex) = DirInt

             OrigLeafTm = RepLeavesTm%D(PIndex)
             x0 = RepLCoor%D(DirInt,PIndex)
             x1 = RepRCoor%D(DirInt,PIndex)
             f0 = -0.50D0
             f1 =  0.50D0

             RootFindIter = 0


#ifdef REGULA_FALSI
             f2 = 1000.0D0
             DO
                RootFindIter = RootFindIter + 1
                IF(RootFindIter > 100) THEN
                   WRITE(*,*) 'Stage = ',Stage,', PIndex = ',PIndex,', CIndex =',CIndex
                   WRITE(*,*) 'ERROR: Regula-falsi does not converged.'
                   STOP
                ENDIF
                x2 = x1 - f1*(x1-x0)/(f1-f0)
                !! calculate f2 now
                !! first define the box
                RepRCoor%D(DirInt,PIndex) = x2
                !! copy the coordinates of RepLCoor and RepRCoor to define a box
                BoxPoint(1:3) = RepLCoor%D(1:3,PIndex)
                BoxPoint(4:6) = RepRCoor%D(1:3,PIndex)
                CALL MPI_BCast(BoxPoint(1),6,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
                LocalLeavesTm = CalLeavesTm(BoxPoint)
                LeavesTmInside = Reduce(LocalLeavesTm)
                oldf2 = f2
                f2 = LeavesTmInside*1.0D0/(OrigLeafTm*1.0D0)-0.5D0
                IF(ABS(f2) < RootTau .OR. ABS(f2-oldf2) < 1.0D-5) THEN
                   ! RepRCoor%D(1:3,PIndex) has been defined
                   ! Messy output is unacceptable and uneccesary:  Protect with ifdef if you want to keep it
                   ! WRITE(*,*) 'RootFindIter = ',RootFindIter,', f2 = ',f2
                   ! WRITE(54,*) x2
                   ETRootArr%D(RootIndex) = x2
                   RepLCoor%D(DirInt,CIndex) = x2
                   RepLeavesTm%D(PIndex) = LeavesTmInside
                   RepLeavesTm%D(CIndex) = OrigLeafTm - LeavesTmInside
                   EXIT
                ELSE
                   !! convergence not achieved, going to the next iteration
                   ! WRITE(*,*) 'RootFindIter = ',RootFindIter,', f2 = ',f2
                   IF(f2*f1 < 0.0D0) THEN
                      x0 = x1; f0 = f1
                      x1 = x2; f1 = f2
                   ELSE
                      x1 = x2; f1 = f2
                   ENDIF
                ENDIF
             ENDDO
#else
             !! bisection
             DO
                RootFindIter = RootFindIter + 1
                IF(RootFindIter > 100) THEN
                   WRITE(*,*) 'Stage = ',Stage,', PIndex = ',PIndex,', CIndex =',CIndex
                   WRITE(*,*) 'ERROR: Regula-falsi does not converged.'
                   STOP
                ENDIF
                x2 = (x0+x1)/2.0D0
                RepRCoor%D(DirInt,PIndex) = x2
                !! copy the coordinates of RepLCoor and RepRCoor to define a box
                BoxPoint(1:3) = RepLCoor%D(1:3,PIndex)
                BoxPoint(4:6) = RepRCoor%D(1:3,PIndex)
                CALL MPI_BCast(BoxPoint(1),6,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
                LocalLeavesTm = CalLeavesTm(BoxPoint)
                LeavesTmInside = Reduce(LocalLeavesTm)
                f2 = LeavesTmInside*1.0D0/(OrigLeafTm*1.0D0)-0.5D0
                IF(ABS(f2) < TauArr%D(Stage)) THEN
                   !! WRITE(*,*) 'RootFindIter = ',RootFindIter,', f2 = ',f2
                   ! WRITE(54,*) x2
                   ETRootArr%D(RootIndex) = x2
                   RepLCoor%D(DirInt,CIndex) = x2
                   RepLeavesTm%D(PIndex) = LeavesTmInside
                   RepLeavesTm%D(CIndex) = OrigLeafTm - LeavesTmInside
                   EXIT
                ELSE
                   !! WRITE(*,*) 'RootFindIter = ',RootFindIter,', f2 = ',f2
                   IF(f2*f0 < 0.0D0) THEN
                      x1 = x2; f1 = f2
                   ELSE
                      x0 = x2; f0 = f2
                   ENDIF
                ENDIF
             ENDDO

#endif

          ENDDO
       ENDDO
       IF(RootIndex /= RootNum) THEN
          WRITE(*,*) 'RootIndex = ',RootIndex, ', RootNum = ',RootNum
          STOP 'RootIndex not equal to RootNum'
       ENDIF
       CALL Put(NPrc,'LineLocExist')
       CALL Put(ETDirArr,'ETDirArr')
       CALL Put(ETRootArr,'ETRootArr')

       ! ask all processors to quit now!
       BoxPoint(1:6) = 0.0D0
       CALL MPI_BCast(BoxPoint(1),6,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)

       !! check the volume before and after are the same
       OldVol = ZERO
       DO I = 1, NVol
          ThisVol = 1.0D0
          DO J = 1, 3
             ThisVol = ThisVol*(RCoor%D(J,I)-LCoor%D(J,I))
          ENDDO
          OldVol = OldVol + ThisVol
       ENDDO

       NewVol = ZERO
       DO I = 1, NVol
          ThisVol = 1.0D0
          DO J = 1, 3
             LinDim = RepRCoor%D(J,I)-RepLCoor%D(J,I)
             !! check for the positive definiteness
             IF(LinDim <= 0) THEN
                WRITE(*,*) 'LinDim is non-positive!'
                WRITE(*,*) 'LinDim = ',LinDim
                WRITE(*,*) 'Box = ',I, ',  Direction = ',J
                STOP
             ENDIF
             ThisVol = ThisVol*LinDim
          ENDDO
          NewVol = NewVol + ThisVol
       ENDDO
       DIFF = NewVol-OldVol
       ! Messy output is unacceptable and uneccesary:  Protect with ifdef if you want to keep it
       !CALL OpenASCII(OutFile,Out)
       !WRITE(Out,*) 'OldVol = ',OldVol, ', NewVol = ', NewVol
       !WRITE(Out,*) 'DIFF = ',DIFF
       !CLOSE(Out,STATUS='KEEP')
       IF(ABS(DIFF) > 1.0D-2) THEN
          STOP 'ERR: DIFF in Vol (RepartionVol) is wrong!'
       ENDIF
    ELSE
       DO
          CALL MPI_BCast(BoxPoint(1),6,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
          AbsSum = Sum( ABS(BoxPoint(:)))
          ! WRITE(*,*) 'MyID = ',MyID,', AbsSum = ', AbsSum
          IF(ABS(AbsSum) < 1.0D-14) THEN
             Busy = .FALSE.
          ENDIF
          IF(.NOT. Busy) EXIT
          LocalLeavesTm = CalLeavesTm(BoxPoint)
          LeavesTmInside = Reduce(LocalLeavesTm)
       ENDDO
    ENDIF
#ifdef REGULA_FALSI
#else
    IF(MyID == 0) THEN
       CALL Delete(TauArr)
    ENDIF
#endif
    CALL Delete(RepLCoor)
    CALL Delete(RepRCoor)
    CALL Delete(RepLeavesTm)
    EndTm = MondoTimer()
    TotTm = EndTm - StartTm
    ! Messy output is unacceptable and uneccesary:  Protect with ifdef if you want to keep it
    !IF(MyID == 0) THEN
    ! CALL OpenASCII(OutFile,Out)
    !WRITE(Out,*) 'Time to repartition the total volume is ',TotTm
    !CLOSE(Out,STATUS='KEEP')
    !ENDIF
  END SUBROUTINE RepartitionVol

  FUNCTION CalLeavesTm(BoxPoint)
    REAL(DOUBLE)::BoxPoint(6),CalLeavesTm,LE(3),UE(3)
    INTEGER::I

    LE(1:3) = CubeRoot%Box%BndBox(1:3,1)
    UE(1:3) = CubeRoot%Box%BndBox(1:3,2)

    DO I = 1, 3
       IF(UE(I) <= LE(I)) STOP 'ERROR: In CalLeavesTm, UE <= LE !'
       IF(BoxPoint(I+3) <= BoxPoint(I)) STOP 'ERROR: BoxPoint(I+3) <= BoxPoint(I)'
    ENDDO
    IF(UE(1) <= BoxPoint(1) .OR. LE(1) >= BoxPoint(4) .OR. &
         UE(2) <= BoxPoint(2) .OR. LE(2) >= BoxPoint(5) .OR. &
         UE(3) <= BoxPoint(3) .OR. LE(3) >= BoxPoint(6)) THEN
       CalLeavesTm = 0
    ELSE IF(LE(1) >= BoxPoint(1) .AND. UE(1) <= BoxPoint(4) .AND. &
         LE(2) >= BoxPoint(2) .AND. UE(2) <= BoxPoint(5) .AND. &
         LE(3) >= BoxPoint(3) .AND. UE(3) <= BoxPoint(6)) THEN
       CalLeavesTm = MyLeavesTm
    ELSE
       CalLeavesTm = RecurCalLeavesTm(CubeRoot)
    ENDIF
  END FUNCTION CalLeavesTm
  !===============================================================================
  RECURSIVE FUNCTION RecurCalLeavesTm(Cube) RESULT(LeavesTm)
    TYPE(CubeNode),Pointer :: Cube
    REAL(DOUBLE)::LeavesTm,OverlapV,CubeV
    REAL(DOUBLE)::C1,C2,C3,NC1,NC2,NC3,LE(3),UE(3),P(3),Q(3),Ovlap(3)
    INTEGER::I

#undef CenterConsideration
#define FractionalConsideration

#ifdef CenterConsideration
    IF(Cube%Leaf) THEN
       C1 = Cube%Box%Center(1)
       C2 = Cube%Box%Center(2)
       C3 = Cube%Box%Center(3)

       IF(C1 >= BoxPoint(1) .AND. C1 < BoxPoint(4) .AND. &
            C2 >= BoxPoint(2) .AND. C2 < BoxPoint(5) .AND. &
            C3 >= BoxPoint(3) .AND. C3 < BoxPoint(6)) THEN
          LeavesTm = Cube%LayGridCost
       ELSE
          LeavesTm = 0.0D0
       ENDIF
#endif
#ifdef FractionalConsideration
       IF(Cube%Leaf) THEN
          !! two possibilities: This leaf partially overlaps with BoxPoint or not!
          LE(1:3) = Cube%Box%BndBox(1:3,1)
          UE(1:3) = Cube%Box%BndBox(1:3,2)

          P(1:3) = BoxPoint(1:3)
          Q(1:3) = BoxPoint(4:6)

          IF(UE(1) <= P(1) .OR. LE(1) >= Q(1) .OR. &
               UE(2) <= P(2) .OR. LE(2) >= Q(2) .OR. &
               UE(3) <= P(3) .OR. LE(3) >= Q(3)) THEN
             LeavesTm = 0.0D0
          ELSE
             DO I = 1, 3
                !! I-direction Overlap
                IF(LE(I) < P(I)) THEN
                   IF(UE(I) <= Q(I)) THEN
                      Ovlap(I) = UE(I) - P(I)
                   ELSE
                      Ovlap(I) = Q(I) - P(I)
                   ENDIF
                ELSE
                   IF(UE(I) <= Q(I)) THEN
                      Ovlap(I) = UE(I) - LE(I)
                   ELSE
                      Ovlap(I) = Q(I) - LE(I)
                   ENDIF
                ENDIF
             ENDDO
             OverlapV = Ovlap(1)*Ovlap(2)*Ovlap(3)
             IF(OverlapV <= 0.0D0) STOP 'ERROR: OverlapV is negative!'
             CubeV = (UE(1)-LE(1))*(UE(2)-LE(2))*(UE(3)-LE(3))
             LeavesTm = (Cube%LayGridCost)*OverlapV/CubeV
          ENDIF
#endif
       ELSE
          LeavesTm = RecurCalLeavesTm(Cube%Descend)+&
                     RecurCalLeavesTm(Cube%Descend%Travrse)
       ENDIF
     END FUNCTION RecurCalLeavesTm

     !===============================================================================
     RECURSIVE FUNCTION LeavesTmCount(Cube) RESULT(TmCount)
       TYPE(CubeNode), POINTER    :: Cube
       REAL(DOUBLE)               :: TmCount
       IF(Cube%Leaf)THEN
          TmCount=Cube%LayGridCost
       ELSE
          TmCount=LeavesTmCount(Cube%Descend)+LeavesTmCount(Cube%Descend%Travrse)
       ENDIF
     END FUNCTION LeavesTmCount

     !===============================================================================
     FUNCTION AtomRad(Z)
       REAL(DOUBLE)::Z,AtomRad
       AtomRad = SLRadii(NINT(Z))*AngstromsToAU
     END FUNCTION AtomRad

     !===============================================================================
     SUBROUTINE Get_SubVol(LCoor,RCoor,N)
       TYPE(DBL_RNK2)::LCoor,RCoor
       INTEGER::N,I,CIndex,J,K,SmallN,PIndex,Stage,DirInt
       INTEGER::Power2(0:31)
       REAL(DOUBLE)::OvLen(3),XEff(3),Sum1,Sum2,CenterZ,&
            ZEff,MaxDim,LinDim,Rad,RangeL,RangeR,BoxL,BoxR
       TYPE(DBL_VECT)::DLCoor,DRCoor

       ! Messy output is unacceptable and uneccesary:  Protect with ifdef if you want to keep it
       !IF(MyID == 0) THEN
       ! WRITE(*,*) 'Center of Mass Partition is implemented...'
       !CALL OpenASCII(OutFile,Out)
       !WRITE(Out,*) 'Center of Mass Partition is implemented...'
       !CLOSE(Out,STATUS='KEEP')
       !ENDIF
       Power2(0) = 1
       DO I = 1,30 !31 is an overflow!
          Power2(I) = Power2(I-1)*2
       ENDDO
       ! check if N is a power of 2
       SmallN = Nint(Log(N*1.0)/Log(2.0D0))
       IF(N /= Power2(SmallN)) THEN
          WRITE(*,*) 'ERROR: N is ', N, ', N is not a power of 2, stop!!!'
          STOP
       ENDIF

       CALL New(DLCoor,3)
       CALL New(DRCoor,3)
       DO Stage = 1, SmallN
          DO PIndex = 1, Power2(Stage-1)
             ! WRITE(*,*) 'Stage = ', Stage, ', PIndex = ',PIndex
             !! CIndex is the index of the child subvolume
             CIndex = PIndex + Power2(Stage-1)
             !! Determine the direction to split
             MaxDim = -1.0D0
             DO J = 1, 3
                LinDim = RCoor%D(J,PIndex)-LCoor%D(J,PIndex)
                IF(LinDim > MaxDim) THEN
                   DirInt = J
                   MaxDim = LinDim
                ENDIF
             ENDDO
             !! The direction is DirInt

             !! now go through all atoms and check if each of them
             !! is inside an expanded volume.
             Sum1 = 0.0D0
             Sum2 = 0.0D0
             ! WRITE(*,*) 'NAtoms is ', NAtoms
             DO J = 1, NAtoms
                ! Rad = AtomRad(GM%AtNum%I(J))
                Rad = AtomRad(GM%AtNum%D(J))
                ! WRITE(*, *) 'Rad = ', Rad
                !! expand the volume
                DLCoor%D(1:3) = LCoor%D(1:3,PIndex) - Rad
                DRCoor%D(1:3) = RCoor%D(1:3,PIndex) + Rad
                IF(GM%Carts%D(1,J) > DLCoor%D(1) .AND. &
                     GM%Carts%D(2,J) > DLCoor%D(2) .AND. &
                     GM%Carts%D(3,J) > DLCoor%D(3) .AND. &
                     GM%Carts%D(1,J) < DRCoor%D(1) .AND. &
                     GM%Carts%D(2,J) < DRCoor%D(2) .AND. &
                     GM%Carts%D(3,J) < DRCoor%D(3)) THEN
                   DO K = 1, 3
                      RangeL = GM%Carts%D(K,J)-Rad
                      RangeR = GM%Carts%D(K,J)+Rad
                      BoxL = LCoor%D(K,PIndex)
                      BoxR = RCoor%D(K,PIndex)
                      !! now check 11 cases
                      IF(RangeL < BoxL .AND. RangeR <= BoxL) THEN
                         WRITE(*,*) 'ERROR: case 1 not possible!'
                         STOP
                      ELSE IF(RangeL<BoxL .AND. RangeR<=BoxR) THEN
                         OvLen(K) = RangeR - BoxL
                         XEff(K) = BoxL+OvLen(K)/2.0D0
                      ELSE IF(RangeL<BoxL .AND. RangeR>BoxR) THEN
                         OvLen(K) = BoxR-BoxL
                         XEff(K) = BoxL+OvLen(K)/2.0D0
                      ELSE IF(RangeL==BoxL .AND. RangeR<=BoxR) THEN
                         OvLen(K) = RangeR-BoxL
                         XEff(K) = BoxL+OvLen(K)/2.0D0
                      ELSE IF(RangeL==BoxL .AND. RangeR>BoxR) THEN
                         OvLen(K) = BoxR-BoxL
                         XEff(K) = BoxL+OvLen(K)/2.0D0
                      ELSE IF(RangeL<BoxR .AND. RangeR<=BoxR) THEN
                         OvLen(K) = RangeR-RangeL
                         XEff(K) = RangeL+OvLen(K)/2.0D0
                      ELSE IF(RangeL<BoxR .AND. RangeR>BoxR) THEN
                         OvLen(K) = BoxR-RangeL
                         XEff(K) = RangeL+OvLen(K)/2.0D0
                      ELSE
                         WRITE(*,*) 'ERROR: A missing case ??'
                         STOP
                      ENDIF
                   ENDDO
                   !! 4.18879D0 is 4Pi/3
                   ZEff = GM%AtNum%D(J)*OvLen(1)*OvLen(2)*OvLen(3)/(4.18879D0*Rad**3.0)
                   ! WRITE(*,*) 'ZEff: GM_AtNum_D = ',GM%AtNum%D(J)
                   Sum1 = Sum1 + ZEff*XEff(DirInt)
                   Sum2 = Sum2 + ZEff
                ENDIF
             ENDDO
             !! check the condition
             IF(Sum2 == 0.0D0) THEN
                STOP 'ERROR: The sum of effective charge in vol PIndex is zero!'
             ENDIF
             CenterZ = Sum1/Sum2
             !! split the box
             RCoor%D(1:3,CIndex) = RCoor%D(1:3,PIndex)
             LCoor%D(1:3,CIndex) = LCoor%D(1:3,PIndex)
             RCoor%D(DirInt,PIndex) = CenterZ
             LCoor%D(DirInt,CIndex) = CenterZ
          ENDDO
       ENDDO
       CALL Delete(DLCoor)
       CALL Delete(DRCoor)

     END SUBROUTINE Get_SubVol
#endif
   END MODULE
