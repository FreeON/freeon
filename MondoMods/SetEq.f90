MODULE SetXYZ
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalObjects
   USE MemMan

#ifdef PARALLEL
   USE MondoMPI
#endif

   IMPLICIT NONE

   INTERFACE SetEq
      MODULE PROCEDURE Set_TIME_EQ_TIME,          &
                       Set_BCSR_EQ_BCSR,          &
                       Set_RNK2_EQ_BCSR,          &
                       Set_BCSR_EQ_RNK2,          &
                       Set_DBL_VECT_EQ_BCSRColVect,&
                       Set_INTC_EQ_INTC,          &
                       Set_AtmB_EQ_AtmB,          &
                       Set_BMATR_EQ_BMATR,        &
                       Set_Chol_EQ_Chol,          &
                       Set_BONDDATA_EQ_BONDDATA,  &
                       Set_CRDS_EQ_CRDS,          &
                       Set_PBCInfo_EQ_PBCInfo,    &
#ifdef PARALLEL
                       Set_DBCSR_EQ_BCSR,         &
                       Set_BCSR_EQ_DBCSR,         &
                       Set_RNK2_EQ_DBCSR,         &
                       Set_DBCSR_EQ_DBCSR,        &
#endif
                       Set_INT_VECT_EQ_INT_SCLR,  &
                       Set_DBL_VECT_EQ_DBL_SCLR,  &
                       Set_INT_VECT_EQ_INT_VECT,  &
                       Set_DBL_VECT_EQ_DBL_VECT,  &
                       Set_DBL_RNK2_EQ_DBL_RNK2,  &
                       Set_CHR10_VECT_EQ_CHR10_VECT, &
                       Set_CellSet_EQ_CellSet
   END INTERFACE

   EXTERNAL bcsr_to_dens

   CONTAINS

!======================================================================
!    Set Vector to
!======================================================================
     SUBROUTINE Set_PBCInfo_EQ_PBCInfo(PBCNew,PBCOld)
       TYPE(PBCInfo) :: PBCNew,PBCOld

       IF(.NOT.AllocQ(PBCNew%Alloc)) CALL Halt('PBCNew needs to be allocated.')
       PBCNew%Dimen=PBCOld%Dimen
       PBCNew%PFFMaxEll=PBCOld%PFFMaxEll
       PBCNew%PFFWelSep=PBCOld%PFFWelSep
       ! Depricated
!!       PBCNew%AtomW=PBCOld%AtomW
       PBCNew%InVecForm=PBCOld%InVecForm
       PBCNew%InAtomCrd=PBCOld%InAtomCrd
       PBCNew%Translate=PBCOld%Translate
       PBCNew%CellVolume=PBCOld%CellVolume
       PBCNew%Epsilon=PBCOld%Epsilon
       PBCNew%DipoleFAC=PBCOld%DipoleFAC
       PBCNew%QupoleFAC=PBCOld%QupoleFAC
       PBCNew%AutoW%I=PBCOld%AutoW%I
       PBCNew%CellCenter%D=PBCOld%CellCenter%D
       PBCNew%TransVec%D=PBCOld%TransVec%D
       PBCNew%BoxShape%D=PBCOld%BoxShape%D
       PBCNew%InvBoxSh%D=PBCOld%InvBoxSh%D
       PBCNew%LatFrc%D=PBCOld%LatFrc%D
     END SUBROUTINE Set_PBCInfo_EQ_PBCInfo

     SUBROUTINE VecToAng(PBC,A,B,C,Alpha,Beta,Gamma)
        TYPE(PBCInfo)               :: PBC
        REAL(DOUBLE)                :: A,B,C,Alpha,Beta,Gamma
        REAL(DOUBLE),PARAMETER      :: DegToRad =  1.745329251994329576923D-2
        A = SQRT(PBC%BoxShape%D(1,1)**2 + PBC%BoxShape%D(2,1)**2+ PBC%BoxShape%D(3,1)**2)
        B = SQRT(PBC%BoxShape%D(1,2)**2 + PBC%BoxShape%D(2,2)**2+ PBC%BoxShape%D(3,2)**2)
        C = SQRT(PBC%BoxShape%D(1,3)**2 + PBC%BoxShape%D(2,3)**2+ PBC%BoxShape%D(3,3)**2)
        Gamma = ACOS((PBC%BoxShape%D(1,1)*PBC%BoxShape%D(1,2))/(A*B))/DegToRad
        Beta  = ACOS((PBC%BoxShape%D(1,1)*PBC%BoxShape%D(1,3))/(A*C))/DegToRad
        Alpha = PBC%BoxShape%D(1,2)*PBC%BoxShape%D(1,3)+PBC%BoxShape%D(2,2)*PBC%BoxShape%D(2,3)
        Alpha = ACOS(Alpha/(B*C))/DegToRad
      END SUBROUTINE VecToAng

!======================================================================
!     Copy a Timer
!======================================================================
      SUBROUTINE Set_TIME_EQ_TIME(T2,T1)
         TYPE(TIME), INTENT(IN)  :: T1
         TYPE(TIME), INTENT(OUT) :: T2
         T2%CPUS=T1%CPUS
         T2%Wall=T1%Wall
         T2%FLOP=T1%FLOP
      END SUBROUTINE Set_TIME_EQ_TIME

!======================================================================
! Set a Column-vector from BCSR form to DBL_VECT form
!
      SUBROUTINE Set_DBL_VECT_EQ_BCSRColVect(A,B)
        TYPE(BCSR) :: B
        TYPE(DBL_VECT) :: A
        INTEGER :: I,J,K,L,II,KK,M

        DO II=1,B%NAtms
          DO I=B%RowPt%I(II),B%RowPt%I(II+1)-1
            IF(B%ColPt%I(I)==1) THEN
              J=II
              K=OffS%I(J)
              L=BSiz%I(J)
              KK=B%BlkPt%I(I)
              DO M=1,L
                A%D(K+M-1)=B%MTrix%D(KK+M-1)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      END SUBROUTINE Set_DBL_VECT_EQ_BCSRColVect

!======================================================================
!     Copy a BCSR matrix
!======================================================================
      SUBROUTINE Set_BCSR_EQ_BCSR(B,A)
         TYPE(BCSR), INTENT(INOUT) :: A
         TYPE(BCSR), INTENT(INOUT) :: B
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            IF(AllocQ(B%Alloc).AND. &
                 & (B%NSMat.NE.A%NSMat.OR.B%NAtms<A%NAtms.OR.B%NBlks<A%NBlks.OR.B%NNon0<A%NNon0) )THEN
               CALL Delete(B)
               CALL New(B,NSMat_O=A%NSMat)
            ELSE
               IF(.NOT.AllocQ(B%Alloc))CALL New(B,NSMat_O=A%NSMat)
            ENDIF
            B%NSMat=A%NSMat; B%NAtms=A%NAtms; B%NBlks=A%NBlks; B%NNon0=A%NNon0
            B%RowPt%I(1:A%NAtms+1)=A%RowPt%I(1:A%NAtms+1)
            B%ColPt%I(1:A%NBlks)  =A%ColPt%I(1:A%NBlks)
            !B%BlkPt%I(1:A%NBlks*A%NSMat)=A%BlkPt%I(1:A%NBlks*A%NSMat)
            B%BlkPt%I(1:A%NBlks)  =A%BlkPt%I(1:A%NBlks)
            B%MTrix%D(1:A%NNon0)  =A%MTrix%D(1:A%NNon0)
#ifdef PARALLEL
         ENDIF
#endif
      END SUBROUTINE Set_BCSR_EQ_BCSR
#ifdef PARALLEL
!======================================================================
!     Copy a DBCSR matrix
!======================================================================
      SUBROUTINE Set_DBCSR_EQ_DBCSR(B,A)
         TYPE(DBCSR), INTENT(INOUT) :: A
         TYPE(DBCSR), INTENT(INOUT) :: B
         LOGICAL                    :: LimitsQ
         IF(AllocQ(B%Alloc).AND. &
              & (B%NSMat.NE.A%NSMat.OR.B%NAtms<A%NAtms.OR.B%NBlks<A%NBlks.OR.B%NNon0<A%NNon0) )THEN
            CALL Delete(B)
            CALL New(B,NSMat_O=A%NSMat)
         ELSE
            IF(.NOT.AllocQ(B%Alloc))CALL New(B,NSMat_O=A%NSMat)
         ENDIF
!        Local
         B%NSMat=A%NSMat
         B%NAtms=A%NAtms
         B%NBlks=A%NBlks
         B%NNon0=A%NNon0
         B%Node=A%Node
         CALL SetEq(B%RowPt,A%RowPt,A%NAtms+1)
         CALL SetEq(B%ColPt,A%ColPt,A%NBlks)
         CALL SetEq(B%BlkPt,A%BlkPt,A%NBlks)
         CALL SetEq(B%MTrix,A%MTrix,A%NNon0)
!        Global
         B%GUpDate=A%GUpDate
         CALL SetEq(B%GRwPt,A%GRwPt,NAtoms+1)
         CALL SetEq(B%GClPt,A%GClPt) ! ??
      END SUBROUTINE Set_DBCSR_EQ_DBCSR
#endif
!======================================================================
!     Transform a BCSR matrix into a dense matrix (Rank 2 array)
!======================================================================
      SUBROUTINE Set_RNK2_EQ_BCSR(B,A)
         TYPE(BCSR),    INTENT(IN)    :: A
         TYPE(DBL_RNK2),INTENT(INOUT) :: B
         INTEGER                      :: NRow,NCol
         INTERFACE
            SUBROUTINE BCSR_TO_DENS(NRow,NCol,NBasF,NSMat,NAtoms,NBlks,NNon0,BSiz,OffS, &
                                    A,MTrix,RowPt,ColPt,BlkPt)
              INTEGER                            ,INTENT(IN)  :: NRow,NCol,NBasF,NSMat,NAtoms, &
                                                                   NBlks,NNon0
              INTEGER, PARAMETER                              :: DOUBLE=KIND(0.D0)
              REAL(DOUBLE),DIMENSION(NRow,NCol),  INTENT(OUT) :: A
              REAL(DOUBLE),DIMENSION(NNon0),      INTENT(IN)  :: MTrix
              INTEGER     ,DIMENSION(NAtoms+1),   INTENT(IN)  :: RowPt
              INTEGER     ,DIMENSION(NBlks),      INTENT(IN)  :: ColPt
              INTEGER     ,DIMENSION(NBlks*NSMat),INTENT(IN)  :: BlkPt
              INTEGER     ,DIMENSION(NAtoms),     INTENT(IN)  :: BSiz,OffS
            END SUBROUTINE BCSR_TO_DENS
         END INTERFACE
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            SELECT CASE(A%NSMat)
            CASE(1);NRow=  NBasF;NCol=  NBasF
            CASE(2);NRow=  NBasF;NCol=2*NBasF
            CASE(4);NRow=2*NBasF;NCol=2*NBasF
            CASE DEFAULT;CALL Halt(' Set_RNK2_EQ_BCSR: A%NSMat doesn''t have an expected value! ')
            END SELECT

            IF(.NOT.AllocQ(B%Alloc))THEN
               CALL New(B,(/NRow,NCol/))
            ELSEIF(SIZE(B%D,1)/=NRow.OR.SIZE(B%D,2)/=NCol)THEN
               CALL Delete(B)
               CALL New(B,(/NRow,NCol/))
            ENDIF


!     SUBROUTINE BCSR_TO_DENS(NRow,NCol,NBasF,NSMat,NAtoms,NBlks,NNon0,
!     >                        MSiz,OffS,A,MTrix,RowPt,ColPt,BlkPt)

!!$write(*,*) 'NRow',NRow
!!$write(*,*) 'NCol',NCol
!!$write(*,*) 'NBasF',NBasF
!!$write(*,*) 'A%NSMat',A%NSMat
!!$write(*,*) 'NAtoms',NAtoms
!!$write(*,*) 'A%NBlks',A%NBlks
!!$write(*,*) 'A%NNon0',A%NNon0
!!$write(*,*) 'BSiz%I',BSiz%I
!!$write(*,*) 'OffS%I',OffS%I
!!$!!write(*,*) 'B%D',B%D
!!$write(*,*) 'A%MTrix%D',A%MTrix%D
!!$write(*,*) 'A%RowPt%I',A%RowPt%I
!!$write(*,*) 'A%ColPt%I',A%ColPt%I
!!$write(*,*) 'A%BlkPt%I',A%BlkPt%I

            CALL BCSR_To_DENS(NRow,NCol,NBasF,A%NSMat,NAtoms,A%NBlks,A%NNon0,BSiz%I,OffS%I, &
                 &            B%D,A%MTrix%D,A%RowPt%I,A%ColPt%I,A%BlkPt%I)

#ifdef PARALLEL
         ENDIF
#endif
      END SUBROUTINE Set_RNK2_EQ_BCSR
!
!---------------------------------------------------------------------
!
      SUBROUTINE Set_RNK2_EQ_BCSR_Dim(B,A,Dim)
        TYPE(BCSR) :: A
        TYPE(DBL_RNK2) :: B
        INTEGER :: Dim,I,J,K,L,KK,NDimS,M,N,MM,NN,II
! Dim: # of real atoms
!
        IF(AllocQ(B%Alloc).AND.SIZE(B%D,1)/=Dim) CALL Delete(B)
        IF(.NOT.AllocQ(B%Alloc)) CALL New(B,(/Dim,Dim/))
        B%D=Zero
!
        DO I=1,A%NAtms
              IF(OffS%I(I)>Dim) EXIT
          DO J=A%RowPt%I(I),A%RowPt%I(I+1)-1
            L=A%ColPt%I(J)
            IF(OffS%I(L)>Dim) CYCLE
            KK=A%BlkPt%I(J)
            II=0
            DO M=1,BSiz%I(L)
              MM=OffS%I(L)-1+M
              DO N=1,BSiz%I(I)
                NN=OffS%I(I)-1+N
                II=II+1
                B%D(MM,NN)=A%MTrix%D(KK-1+II)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      END SUBROUTINE Set_RNK2_EQ_BCSR_Dim
!
!---------------------------------------------------------------------
!
!======================================================================
!     Transform a dense matrix (Rank 2 array) into a BCSR matrix
!======================================================================
      SUBROUTINE Set_BCSR_EQ_RNK2(B,A,NSMat_O)
         TYPE(DBL_RNK2), INTENT(INOUT) :: A
         TYPE(BCSR),     INTENT(INOUT) :: B
         INTEGER,OPTIONAL              :: NSMat_O
         INTEGER                       :: I,J,P,Q,OI,OJ,OI0,OJ0,MA,NA,NSMat,iSMat
         NSMat=1
         IF(PRESENT(NSMat_O))NSMat=NSMat_O
         IF(AllocQ(B%Alloc))  &
         CALL Delete(B)
         CALL New(B,NSMat_O=NSMat)
         P=1
         Q=1
         OI0=0
         B%RowPt%I(1)=1
         DO I=1,NAtoms
            OJ0=0
            MA=BSiz%I(I)
            DO J=1,NAtoms
               NA=BSiz%I(J)
               B%BlkPt%I(P)=Q
               DO iSMat=1,NSMat
                  IF(iSMat.EQ.1)THEN
                     OI=OI0
                     OJ=OJ0
                  ELSEIF(iSMat.EQ.2)THEN
                     OI=OI0
                     OJ=OJ0+NBasF
                  ELSEIF(iSMat.EQ.3)THEN
                     OI=OI0+NBasF
                     OJ=OJ0
                  ELSEIF(iSMat.EQ.4)THEN
                     OI=OI0+NBasF
                     OJ=OJ0+NBasF
                  ELSE
                     STOP 'Set_BCSR_EQ_RNK2'
                  ENDIF
                  !CALL BlockToBlock(MA,NA,OI,OJ,A%D,B%MTrix%D(Q:))
                  CALL BlockToBlock_(MA,NA,OI,OJ,A,B%MTrix%D(Q))
                  Q=Q+MA*NA
               ENDDO
               B%ColPt%I(P)=J
               P=P+1
               B%RowPt%I(I+1)=P
               OJ0=OJ0+NA
            ENDDO
            OI0=OI0+MA
         ENDDO
         B%NAtms=NAtoms
         B%NBlks=P-1
         B%NNon0=Q-1
      END SUBROUTINE Set_BCSR_EQ_RNK2

      SUBROUTINE BlockToBlock_(M,N,OI,OJ,A,B)
        IMPLICIT NONE
        INTEGER :: M,N,OI,OJ
        TYPE(DBL_RNK2), INTENT(INOUT) :: A
        REAL(DOUBLE),DIMENSION(M,N) :: B
        INTEGER :: I,J
        DO J=1,N
           DO I=1,M
              B(I,J)=A%D(OI+I,OJ+J)
           ENDDO
        ENDDO
      END SUBROUTINE BlockToBlock_

      SUBROUTINE BlockToBlock(M,N,OI,OJ,A,B)
         IMPLICIT NONE
         INTEGER,                     INTENT(IN)    :: M,N,OI,OJ
         REAL(DOUBLE),DIMENSION(:,:), INTENT(IN)    :: A
         REAL(DOUBLE),DIMENSION(M,N), INTENT(INOUT) :: B
         INTEGER                                    :: I,J

         DO J=1,N
            DO I=1,M
               B(I,J)=A(OI+I,OJ+J)
            ENDDO
         ENDDO
      END SUBROUTINE BlockToBlock
!
      SUBROUTINE BlockToBlock2(M,N,A,B)
         IMPLICIT NONE
         INTEGER,                     INTENT(IN)    :: M,N
         REAL(DOUBLE),DIMENSION(:,:), INTENT(IN)    :: A
         REAL(DOUBLE),DIMENSION(M,N), INTENT(INOUT) :: B
         INTEGER                                    :: I,J

           DO J=1,N
              DO I=1,M
                 B(I,J)=A(I,J)
              ENDDO
           ENDDO
      END SUBROUTINE BlockToBlock2
!
      SUBROUTINE Set_BCSR_EQ_VECT(B,A)
         ! turn column vector into sparse matrix (BCSR) form
         TYPE(DBL_VECT), INTENT(INOUT) :: A
         TYPE(BCSR),     INTENT(INOUT) :: B
         INTEGER                       :: I,J,P,Q,OI,OJ,MA,NA
         IF(AllocQ(B%Alloc))  &
         CALL Delete(B)
         CALL New(B)
         P=1
         Q=1
         B%RowPt%I(1)=1
         DO I=1,NAtoms
            OI=OffS%I(I)
            MA=BSiz%I(I)
            NA=BSiz%I(1)
            CALL VectToBlock2(MA,NA,A%D(OI:),B%MTrix%D(Q:))
            B%BlkPt%I(P)=Q
            B%ColPt%I(P)=1
            Q=Q+MA*NA
            P=P+1
            B%RowPt%I(I+1)=P
         ENDDO
         B%NAtms=NAtoms
         B%NBlks=P-1
         B%NNon0=Q-1
      END SUBROUTINE Set_BCSR_EQ_VECT

      SUBROUTINE VectToBlock2(M,N,A,B)
         IMPLICIT NONE
         INTEGER,                     INTENT(IN)    :: M,N
         REAL(DOUBLE),DIMENSION(:),   INTENT(IN)    :: A
         REAL(DOUBLE),DIMENSION(M,N), INTENT(INOUT) :: B
         INTEGER                                    :: I,J
         B=Zero
         DO I=1,M
            B(I,1)=A(I)
         ENDDO
      END SUBROUTINE VectToBlock2
!
      SUBROUTINE Set_VECT_EQ_BCSR(B,A)
         TYPE(BCSR),     INTENT(INOUT) :: A
         TYPE(DBL_VECT), INTENT(INOUT) :: B
         INTEGER                       :: I,J,P,Q,OI,JP,MA,NA
         DO I=1,NAtoms
            OI=OffS%I(I)
            MA=BSiz%I(I)
            NA=BSiz%I(1)
            DO JP=A%RowPt%I(I),A%RowPt%I(I+1)-1
               J=A%ColPt%I(JP)
               IF(J==1)THEN
                  Q=A%BlkPt%I(JP)
                  CALL BlockToVect2(MA,NA,B%D(OI:),A%MTrix%D(Q:))
               ENDIF
            ENDDO
         ENDDO
      END SUBROUTINE Set_VECT_EQ_BCSR

      SUBROUTINE BlockToVect2(M,N,A,B)
         IMPLICIT NONE
         INTEGER,                     INTENT(IN)    :: M,N
         REAL(DOUBLE),DIMENSION(:),   INTENT(INOUT) :: A
         REAL(DOUBLE),DIMENSION(M,N), INTENT(IN)    :: B
         INTEGER                                    :: I,J
         DO I=1,M
            A(I)=B(I,1)
         ENDDO
      END SUBROUTINE BlockToVect2

#ifdef PARALLEL
!===================================================================
!     Load the Global Row Pointer from the local Row Pointer
!===================================================================
      SUBROUTINE Load(A,RowPt,Id_O)
         TYPE(DBCSR),     INTENT(IN)    :: A
         TYPE(INT_VECT),  INTENT(INOUT) :: RowPt
         INTEGER,OPTIONAL,INTENT(IN)    :: Id_O
         INTEGER                        :: I,K,Id,I1,I2
!--------------------------------------------------------------------
         Id=A%Node; IF(PRESENT(Id_O))Id=Id_O
         I1=OffSt%I(Id)+1
         I2=I1+A%NAtms
         IF(I2>SIZE(RowPt%I))CALL Halt(' Indexing hosed in Load ')
         K=1
         DO I=I1,I2
            RowPt%I(I)=A%RowPt%I(K)
            K=K+1
         ENDDO
      END SUBROUTINE Load
!===================================================================
!     Zero out the Global Row Pointer
!===================================================================
      SUBROUTINE Clear(A,RowPt,Id_O)
         TYPE(DBCSR),     INTENT(IN)    :: A
         TYPE(INT_VECT),  INTENT(INOUT) :: RowPt
         INTEGER,OPTIONAL,INTENT(IN)    :: Id_O
         INTEGER                        :: I,Id,I1,I2
!--------------------------------------------------------------------
         Id=A%Node; IF(PRESENT(Id_O))Id=Id_O
         I1=OffSt%I(Id)+1
         I2=I1+A%NAtms
         DO I=I1,I2
            RowPt%I(I)=0
         ENDDO
      END SUBROUTINE Clear
!======================================================================
!     Transform a DBCSR matrix into locally a dense matrices,
!     one for each node.  Usefull only for printing.
!======================================================================
      SUBROUTINE Set_RNK2_EQ_DBCSR(B,A,Id_O)
         TYPE(DBCSR),    INTENT(INOUT)       :: A
         TYPE(DBL_RNK2),INTENT(INOUT)     :: B
         INTEGER,OPTIONAL                 :: Id_O
         INTEGER                          :: K,L,JP,P,II,JJ,MA,MA1,NA, &
                                             MOff,KOff,IStrt,IStop
         INTEGER                          :: IL,IG,JG,Id
         IF(B%Alloc/=ALLOCATED_TRUE.AND.B%Alloc/=ALLOCATED_TEMP)THEN
            CALL New(B,(/NBasF,NBasF/))
         ELSEIF(SIZE(B%D,1)/=NBasF.OR.SIZE(B%D,2)/=NBasF)THEN
            CALL Delete(B)
            CALL New(B,(/NBasF,NBasF/))
         ENDIF
         Id=A%Node; IF(PRESENT(Id_O))Id=Id_O
         B%D(1:NBasF,1:NBasF)=Zero
         DO IL=1,A%NAtms
            IG=OffSt%I(Id)+IL
            MA=BSiz%I(IG)
            II=OffS%I(IG)
            MA1=MA-1
            IStrt=A%RowPt%I(IL)
            IStop=A%RowPt%I(IL+1)-1
            IF(IStrt/=0.AND.IStop/=0)THEN
            DO JP=IStrt,IStop
               JG=A%ColPt%I(JP)
               P =A%BlkPt%I(JP)
               NA=BSiz%I(JG)
               JJ=OffS%I(JG)
               DO K=0,NA-1
                  KOff=JJ+K
                  MOff=P+K*MA
                  DO L=0,MA1
                     B%D(II+L,KOff)=A%MTrix%D(MOff+L)
                  ENDDO
               ENDDO
            ENDDO
            ENDIF
         ENDDO
      END SUBROUTINE Set_RNK2_EQ_DBCSR
!======================================================================
!     Gather a distributed BCSR matrix to a serial BCSR matrix on ROOT
!======================================================================
      SUBROUTINE Set_BCSR_EQ_DBCSR(B,A)
         TYPE(DBCSR), INTENT(INOUT) :: A
         TYPE(BCSR),  INTENT(INOUT) :: B
         INTEGER                    :: I,K,NAtms
         TYPE(INT_VECT)             :: MA,MB,MN,NA,NB,NN
!-----------------------------------------------------------------------
!        Allocate indecies
!
         IF(MyId==ROOT)THEN
            CALL New(MA,NPrc,M_O=0)
            CALL New(MB,NPrc,M_O=0)
            CALL New(MN,NPrc,M_O=0)
            CALL New(NA,NPrc,M_O=0)
            CALL New(NB,NPrc,M_O=0)
            CALL New(NN,NPrc,M_O=0)
         ENDIF
!----------------------------------------------------------
!        Allocate the serial matrix
!
         IF(AllocQ(B%Alloc))CALL Delete(B)
         B%NAtms=NAtoms
         B%NBlks=Reduce(A%NBlks)
         B%NNon0=Reduce(A%NNon0)
!-1
!MINUS
         CALL New(B,(/B%NAtms,B%NBlks,B%NNon0/))
         B%NSMat=A%NSMat
!------------------------------------
!        Number of atoms per node
!
         NAtms=End%I(MyId)-Beg%I(MyId)+1
         IF(MyID==NPrc-1)NAtms=NAtms+1
!------------------------------------
!        Gather the indecies to root
!
         CALL Gather(NAtms  ,MA)
         CALL Gather(A%NBlks,MB)
         CALL Gather(A%NNon0,MN)
         IF(MyID==ROOT)THEN
!---------------------------------------------------
!           Calculate index offsets (displacements)
!
            NA%I(0)=  NAtms
            NB%I(0)=A%NBlks
            NN%I(0)=A%NNon0
            DO I=1,NPrc-1
               NA%I(I)=MA%I(I)+NA%I(I-1)
               NB%I(I)=MB%I(I)+NB%I(I-1)
               NN%I(I)=MN%I(I)+NN%I(I-1)
            ENDDO
            DO I=NPrc,1,-1
               NA%I(I)=NA%I(I-1)
               NB%I(I)=NB%I(I-1)
               NN%I(I)=NN%I(I-1)
            ENDDO
            NA%I(0)=0
            NB%I(0)=0
            NN%I(0)=0
!---------------------------------------------
!           Copy local portions of the matirx
!
            CALL SetEq(B%RowPt,A%RowPt,  NAtms)
            CALL SetEq(B%ColPt,A%ColPt,A%NBlks)
            CALL SetEq(B%BlkPt,A%BlkPt,A%NBlks)
            CALL SetEq(B%MTrix,A%MTrix,A%NNon0)
         ENDIF
!----------------------------------
!        Gather the matrix to ROOT
!
         CALL Gather(A%RowPt,B%RowPt,  NAtms,MA,NA)
         CALL Gather(A%ColPt,B%ColPt,A%NBlks,MB,NB)
         CALL Gather(A%BlkPt,B%BlkPt,A%NBlks,MB,NB)
         CALL Gather(A%MTrix,B%MTrix,A%NNon0,MN,NN)
!--------------------------------------------------
!        Add offsets to achieve correct indexing
!
         IF(MyId==ROOT)THEN
            DO I=1,NPrc-1
               DO K=NA%I(I)+1,NA%I(I+1)
                  B%RowPt%I(K)=B%RowPt%I(K)+NB%I(I)
               ENDDO
               DO K=NB%I(I)+1,NB%I(I+1)
                  B%BlkPt%I(K)=B%BlkPt%I(K)+NN%I(I)
               ENDDO
            ENDDO
         ENDIF
!------------------
!        Tidy up
!
         IF(MyId==ROOT)THEN
            CALL Delete(MA)
            CALL Delete(MB)
            CALL Delete(MN)
            CALL Delete(NA)
            CALL Delete(NB)
            CALL Delete(NN)
         ENDIF
      END SUBROUTINE Set_BCSR_EQ_DBCSR
!============================================================================
!     Scatter a serial BCSR matrix from ROOT to a distributed BCSR matrix
!============================================================================
      SUBROUTINE Set_DBCSR_EQ_BCSR(B,A)
         TYPE(BCSR),  INTENT(INOUT)  :: A
         TYPE(DBCSR), INTENT(INOUT)  :: B
         INTEGER                     :: I,Id,J,JG,M,MN,MN1,P,  &
                                        NAtms,NBlks,NNon0
         LOGICAL                     :: ReAllocate
!-----------------------------------------------------------------------
!        Allocate if required
!
         CALL BCast(A%NSMat)
         IF(.NOT.AllocQ(B%Alloc))CALL New(B,NSMat_O=A%NSMat)
         IF(B%NSMat.NE.A%NSMat) THEN
            CALL Delete(B)
            CALL New(B,NSMat_O=A%NSMat)
         ENDIF
!------------------------------------------------
!        Distribute to each processor
!
         DO Id=NPrc-1,0,-1
            IF(MyId==ROOT)THEN
               B%NAtms=0
               B%NBlks=1
               B%NNon0=1
               B%RowPt%I(1)=1
               DO I=Beg%I(Id),End%I(Id)
                  M=BSiz%I(I)
                  B%NAtms=B%NAtms+1
                  DO J=A%RowPt%I(I),A%RowPt%I(I+1)-1
                     JG=A%ColPt%I(J)
                     B%ColPt%I(B%NBlks)=JG
                     B%BlkPt%I(B%NBlks)=B%NNon0
                     B%NBlks=B%NBlks+1
                     MN=M*BSiz%I(JG)     *A%NSMat !<<< SPIN
                     MN1=MN-1
                     P=A%BlkPt%I(J)
                     B%MTrix%D(B%NNon0:B%NNon0+MN1)=A%MTrix%D(P:P+MN1)
                     B%NNon0=B%NNon0+MN
                  ENDDO
                  B%RowPt%I(B%NAtms+1)=B%NBlks
               ENDDO
               B%NBlks=B%NBlks-1
               B%NNon0=B%NNon0-1
! MINUS
               IF(Id/=ROOT)THEN
                  CALL Send(B%NAtms,Id,1)
                  CALL Send(B%NBlks,Id,2)
                  CALL Send(B%NNon0,Id,3)
                  CALL Send(B%RowPt,B%NAtms+1,Id,4)
                  CALL Send(B%ColPt,B%NBlks,Id,5)
                  CALL Send(B%BlkPt,B%NBlks,Id,6)
                  CALL Send(B%MTrix,B%NNon0,Id,7)
               ENDIF
            ELSEIF(MyId==Id)THEN
               CALL Recv(B%NAtms,ROOT,1)
               CALL Recv(B%NBlks,ROOT,2)
               CALL Recv(B%NNon0,ROOT,3)
               ReAllocate=(SIZE(B%RowPt%I)<B%NAtms+1).OR. &
                          (SIZE(B%ColPt%I)<B%NBlks)  .OR. &
                          (SIZE(B%BlkPt%I)<B%NBlks)  .OR. &
                          (SIZE(B%MTrix%D)<B%NNon0)
               IF(ReAllocate)THEN
                  NAtms=B%NAtms; NBlks=B%NBlks; NNon0=B%NNon0
                  CALL Delete(B)
                  CALL New(B,(/NAtms,NBlks,NNon0/))
               ENDIF
               CALL Recv(B%RowPt,B%NAtms+1,ROOT,4)
               CALL Recv(B%ColPt,B%NBlks,ROOT,5)
               CALL Recv(B%BlkPt,B%NBlks,ROOT,6)
               CALL Recv(B%MTrix,B%NNon0,ROOT,7)
            ENDIF
         ENDDO
!-------------------------------------------------------
         CALL BCast(A%NBlks)
         IF(MyId==ROOT)THEN
           CALL SetEq(B%GRwPt,A%RowPt,NAtoms+1)
           CALL SetEq(B%GClPt,A%ColPt,A%NBlks)
         ENDIF
         CALL BCast(B%GRwPt,NAtoms+1)
         CALL BCast(B%GClPt,A%NBlks)
         B%Node=MyId
         B%GUpDate=STATUS_TRUE
!
      END SUBROUTINE Set_DBCSR_EQ_BCSR
!===================================================================
!     Load the Global Row Pointer from the local Row Pointer
!===================================================================
      SUBROUTINE Load_DBCSR_GRwPt(A,Id_O)
         TYPE(DBCSR),     INTENT(INOUT) :: A
         INTEGER,OPTIONAL,INTENT(IN)    :: Id_O
         INTEGER                        :: I,K,Id,I1,I2
!--------------------------------------------------------------------
         IF(PRESENT(Id_O))THEN
            Id=Id_O
         ELSE
            Id=MyId
         ENDIF
         K=1
         I1=OffSt%I(Id)+1
         I2=I1+A%NAtms
         IF(I2>NAtoms+1)CALL Halt(' Indexing hosed in Load_DBCSR_GRwPt ')
         DO I=I1,I2
            A%GRwPt%I(I)=A%RowPt%I(K)
            K=K+1
         ENDDO
         CALL AlignNodes()
      END SUBROUTINE Load_DBCSR_GRwPt
!===================================================================
!     Zero out the Global Row Pointer
!===================================================================
      SUBROUTINE Clear_DBCSR_GRwPt(A,Id_O)
         TYPE(DBCSR),     INTENT(INOUT) :: A
         INTEGER,OPTIONAL,INTENT(IN)    :: Id_O
         INTEGER                        :: I,Id,I1,I2
!--------------------------------------------------------------------
         IF(PRESENT(Id_O))THEN
            Id=Id_O
         ELSE
            Id=MyId
         ENDIF
         I1=OffSt%I(Id)+1
         I2=I1+A%NAtms
         DO I=I1,I2
            A%GRwPt%I(I)=0
         ENDDO
      END SUBROUTINE Clear_DBCSR_GRwPt
#endif
!===============================================================
!     Copy a double vector to a double vector: Y(1:N)=X(1:N)
!===============================================================
      SUBROUTINE Set_DBL_VECT_EQ_DBL_VECT(Y,X,N_O)
         TYPE(DBL_VECT), INTENT(IN)  :: X
         TYPE(DBL_VECT), INTENT(INOUT) :: Y
         INTEGER,OPTIONAL,INTENT(IN) :: N_O
         INTEGER                     :: N
         INTERFACE
            SUBROUTINE DBL_VECT_EQ_DBL_VECT(N,Y,X)
               INTEGER, INTENT(IN)                   :: N
               INTEGER, PARAMETER                    :: DOUBLE=KIND(0.D0)
               REAL(DOUBLE),DIMENSION(N),INTENT(IN)  :: X
               REAL(DOUBLE),DIMENSION(N),INTENT(OUT) :: Y
            END SUBROUTINE DBL_VECT_EQ_DBL_VECT
         END INTERFACE
         N=SIZE(Y%D); IF(PRESENT(N_O))N=N_O
         IF(SIZE(X%D)<N.OR.SIZE(Y%D)<N) &
            CALL Halt(' Dimensions off in Set_DBL_VECT_DBL_VECT')
         CALL DBL_VECT_EQ_DBL_VECT(N,Y%D,X%D)
      END SUBROUTINE Set_DBL_VECT_EQ_DBL_VECT
!===============================================================
!     Copy an integer vector to an integer vector: Y(1:N)=X(1:N)
!===============================================================
      SUBROUTINE Set_INT_VECT_EQ_INT_VECT(Y,X,N_O)
         TYPE(INT_VECT), INTENT(IN)    :: X
         TYPE(INT_VECT), INTENT(INOUT) :: Y
         INTEGER,OPTIONAL,INTENT(IN)   :: N_O
         INTEGER                       :: N
         INTERFACE
            SUBROUTINE INT_VECT_EQ_INT_VECT(N,Y,X)
               INTEGER, INTENT(IN)              :: N
               INTEGER,DIMENSION(N),INTENT(IN)  :: X
               INTEGER,DIMENSION(N),INTENT(OUT) :: Y
            END SUBROUTINE INT_VECT_EQ_INT_VECT
         END INTERFACE
         N=SIZE(Y%I); IF(PRESENT(N_O))N=N_O
         IF(SIZE(X%I)<N.OR.SIZE(Y%I)<N) &
            CALL Halt(' Dimensions off in Set_INT_VECT_INT_VECT')
         CALL INT_VECT_EQ_INT_VECT(N,Y%I,X%I)
      END SUBROUTINE Set_INT_VECT_EQ_INT_VECT
!===============================================================
!     Set an integer scalar to an integer vector: Y(1:N)=X
!===============================================================
      SUBROUTINE Set_INT_VECT_EQ_INT_SCLR(Y,X,N_O)
         INTEGER,        INTENT(IN)    :: X
         TYPE(INT_VECT), INTENT(INOUT) :: Y
         INTEGER,OPTIONAL,INTENT(IN)   :: N_O
         INTEGER                       :: N
         INTERFACE
            SUBROUTINE INT_VECT_EQ_INT_SCLR(N,Y,X)
               INTEGER,INTENT(IN)               :: N
               INTEGER,INTENT(IN)               :: X
               INTEGER,DIMENSION(N),INTENT(OUT) :: Y
            END SUBROUTINE INT_VECT_EQ_INT_SCLR
         END INTERFACE
         N=SIZE(Y%I); IF(PRESENT(N_O))N=N_O
         IF(SIZE(Y%I)<N) &
            CALL Halt(' Dimensions off in Set_INT_SCLR_INT_VECT')
         CALL INT_VECT_EQ_INT_SCLR(N,Y%I,X)
      END SUBROUTINE Set_INT_VECT_EQ_INT_SCLR
!===============================================================
!     Set a double scalar to a double vector: Y(1:N)=X
!===============================================================
      SUBROUTINE Set_DBL_VECT_EQ_DBL_SCLR(Y,X,N_O)
         REAL(DOUBLE),   INTENT(IN)    :: X
         TYPE(DBL_VECT), INTENT(INOUT) :: Y
         INTEGER,OPTIONAL,INTENT(IN)   :: N_O
         INTEGER                       :: N
         INTERFACE
            SUBROUTINE DBL_VECT_EQ_DBL_SCLR(N,Y,X)
               INTEGER,INTENT(IN)                    :: N
               INTEGER, PARAMETER                    :: DOUBLE=KIND(0.D0)
               REAL(DOUBLE),INTENT(IN)               :: X
               REAL(DOUBLE),DIMENSION(N),INTENT(OUT) :: Y
            END SUBROUTINE DBL_VECT_EQ_DBL_SCLR
         END INTERFACE
         N=SIZE(Y%D); IF(PRESENT(N_O))N=N_O
         IF(SIZE(Y%D)<N) &
            CALL Halt(' Dimensions off in Set_DBL_SCLR_DBL_VECT')
         CALL DBL_VECT_EQ_DBL_SCLR(N,Y%D,X)
      END SUBROUTINE Set_DBL_VECT_EQ_DBL_SCLR
!
!===============================================================
!
      SUBROUTINE Set_AtmB_EQ_AtmB(A,B,N2_O)
        TYPE(ATOMBONDS)  :: A,B
        INTEGER          :: N2,I,J
        INTEGER,OPTIONAL :: N2_O
        !
        N2=B%N2
        IF(PRESENT(N2_O)) N2=N2_O
        IF(.NOT.AllocQ(A%Alloc)) THEN
          CALL New(A,B%N1,N2)
        ELSE
          CALL Delete(A)
          CALL New(A,B%N1,N2)
        ENDIF
        A%Count%I=B%Count%I
        DO I=1,B%N1
          DO J=1,B%N2
            A%Bonds%I(I,J)=B%Bonds%I(I,J)
            A%Atoms%I(I,J)=B%Atoms%I(I,J)
          ENDDO
        ENDDO
      END SUBROUTINE Set_AtmB_EQ_AtmB
!
!===============================================================
!
      SUBROUTINE Set_BONDDATA_EQ_BONDDATA(A,B,NewDim_O,OldDim_O)
        TYPE(BONDDATA)          :: A,B
        INTEGER,OPTIONAL        :: NewDim_O,OldDim_O
        INTEGER                 :: I,NewDim,OldDim

        ! In case we have not allocated B, we don't have to worry about copying
        ! any bond information from B to A.
        IF(.NOT.AllocQ(B%Alloc) .AND. B%N /= 0) THEN
          WRITE(*,*) "[Set_BONDDATA_EQ_BONDDATA] B is not allocated and B%N /= 0!"
          CALL Halt("Fatal")
          OldDim = 0
        ELSE
          OldDim=B%N
        ENDIF
        NewDim=OldDim
        IF(PRESENT(NewDim_O)) NewDim=NewDim_O
        IF(PRESENT(OldDim_O) .AND. AllocQ(B%Alloc)) OldDim=OldDim_O

        IF(.NOT.AllocQ(A%Alloc)) THEN
          CALL New(A,NewDim)
        ELSE
          CALL Delete(A)
          CALL New(A,NewDim)
        ENDIF
        A%N=NewDim
        DO I=1,OldDim
          CALL Set_Bond_EQ_Bond(A,I,B,I)
        ENDDO
      END SUBROUTINE Set_BONDDATA_EQ_BONDDATA

      !===============================================================================
      ! Make a deep copy of the CRDS structure
      ! Also figure out PBC issue.
      !
      ! G1 = G2
      !===============================================================================
      SUBROUTINE Set_CRDS_EQ_CRDS(G1,G2)
        TYPE(CRDS) :: G1,G2

        IF(.NOT. AllocQ(G1%Alloc)) THEN
          CALL Halt("[Set_CRDS_EQ_CRDS] G1 is not allocated")
        ENDIF

        G1%InAU=G2%InAU
        G1%NElec=G2%NElec
        G1%Ordrd=G2%Ordrd
        G1%Confg=G2%Confg
        G1%NElec=G2%NElec
        G1%Multp=G2%Multp
        G1%TotCh=G2%TotCh
        G1%NAlph=G2%NAlph
        G1%NBeta=G2%NBeta
        G1%ETotal=G2%ETotal
        CALL SetEq(G1%ETotalPerSCF, G2%ETotalPerSCF)
        G1%GradRMS=G2%GradRMS
        G1%GradMax=G2%GradMax
        G1%Unstable=G2%Unstable
        CALL SetEq(G1%BndBox, G2%BndBox)
        CALL SetEq(G1%PBC, G2%PBC)
        ! A freshly New'ed CRDS will not have properly allocated OvCells and
        ! InCells. We need to be careful.
        IF(AllocQ(G2%OvCells%Alloc) .AND. .NOT. AllocQ(G1%OvCells%Alloc)) THEN
          CALL New(G1%OvCells, G2%OvCells%NCells)
        ENDIF
        IF(AllocQ(G2%InCells%Alloc) .AND. .NOT. AllocQ(G1%InCells%Alloc)) THEN
          CALL New(G1%InCells, G2%InCells%NCells)
        ENDIF
        CALL SetEq(G1%OvCells, G2%OvCells)
        CALL SetEq(G1%InCells, G2%InCells)
        G1%NAtms=G2%NAtms
        G1%Nkind=G2%Nkind
        CALL SetEq(G1%AtNum, G2%AtNum)
        CALL SetEq(G1%AtTyp, G2%AtTyp)
        CALL SetEq(G1%AtNam, G2%AtNam)
        CALL SetEq(G1%AtMss, G2%AtMss)
        CALL SetEq(G1%CConstrain, G2%CConstrain)
        CALL SetEq(G1%DoFreq, G2%DoFreq)
        CALL SetEq(G1%Carts, G2%Carts)
        CALL SetEq(G1%BoxCarts, G2%BoxCarts)
        CALL SetEq(G1%Velocity, G2%Velocity)
        CALL SetEq(G1%Gradients, G2%Gradients)
        CALL SetEq(G1%Fext, G2%Fext)
        CALL SetEq(G1%Displ, G2%Displ)
        CALL SetEq(G1%PBCDispl, G2%PBCDispl)
        G1%LatticeOnly=G2%LatticeOnly
        G1%AltCount=G2%AltCount
      END SUBROUTINE Set_CRDS_EQ_CRDS

      SUBROUTINE Set_Bond_EQ_Bond(A,IA,B,IB)
        TYPE(BONDDATA) :: A,B
        INTEGER        :: IA,IB

        A%IJ%I(1:2,IA)=B%IJ%I(1:2,IB)
        A%Length%D(IA)=B%Length%D(IB)
        A%Type%C(IA)=B%Type%C(IB)
      END SUBROUTINE Set_Bond_EQ_Bond

      SUBROUTINE  Set_BMATR_EQ_BMATR(A,B)
        TYPE(BMATR) :: A,B
        INTEGER     :: NIntC

        NIntC=SIZE(B%IB%I,1)

        IF(.NOT.AllocQ(A%Alloc)) THEN
          CALL New(A,NIntC)
        ELSE
          IF(SIZE(A%IB%I,1)/=NIntC) THEN
            CALL Delete(A)
            CALL New(A,NIntC)
          ENDIF
        ENDIF
!
        A%IB%I =B%IB%I
        A%B%D  =B%B%D
        A%BLI%I =B%BLI%I
        A%BL%D =B%BL%D
      END SUBROUTINE  Set_BMATR_EQ_BMATR
!
!===============================================================
!
      SUBROUTINE Set_Chol_EQ_Chol(A,B)
        TYPE(Cholesky) :: A,B
        INTEGER :: NCart,ChNon0
        NCart=SIZE(B%Perm%I)
        ChNon0=B%ChRowPt%I(NCart+1)-1
        IF(.NOT.AllocQ(A%Alloc)) CALL New(A,NCart,ChNon0)
        A%GcScale%D=B%GcScale%D
        A%ChFact%D=B%ChFact%D
        A%ChDiag%D=B%ChDiag%D
        A%Perm%I=B%Perm%I
        A%IPerm%I=B%IPerm%I
        A%ChRowPt%I=B%ChRowPt%I
        A%ChColPt%I=B%ChColPt%I
      END SUBROUTINE Set_Chol_EQ_Chol
!
!===============================================================
!
      SUBROUTINE Set_INTC_EQ_INTC(IntCs,IntCsCopy,From,To,Start)
!
        TYPE(INTC) :: IntCs,IntCsCopy
        INTEGER    :: I,J,From,To,Start,NSize,To2

        NSize=IntCsCopy%N
        IF(NSize<Start+(To-From)) THEN
          CALL Halt('Dimensionality error in IntCCopy')
        ENDIF

        To2=Start+(To-From)

        IF(To2 < Start) THEN
          CALL MondoLog(DEBUG_MAXIMUM, "Set_INTC_EQ_INTC", "To2 < Start")
          RETURN
        ELSE IF(To < From) THEN
          CALL MondoLog(DEBUG_MAXIMUM, "Set_INTC_EQ_INTC", "To < From")
          RETURN
        ENDIF

        IntCsCopy%Def%C(Start:To2)        =IntCs%Def%C(From:To)
        IntCsCopy%Atoms%I(Start:To2,1:4)  =IntCs%Atoms%I(From:To,1:4)
        IntCsCopy%Cells%I(Start:To2,1:12) =IntCs%Cells%I(From:To,1:12)
        IntCsCopy%Value%D(Start:To2)      =IntCs%Value%D(From:To)
        IntCsCopy%Constraint%L(Start:To2) =IntCs%Constraint%L(From:To)
        IntCsCopy%ConstrValue%D(Start:To2)=IntCs%ConstrValue%D(From:To)
        IntCsCopy%Active%L(Start:To2)     =IntCs%Active%L(From:To)
        IntCsCopy%PredVal%D(Start:To2)    =IntCs%PredVal%D(From:To)
        IntCsCopy%PredGrad%D(Start:To2)   =IntCs%PredGrad%D(From:To)
!
      END SUBROUTINE Set_INTC_EQ_INTC
!===============================================================
   SUBROUTINE CartRNK1ToCartRNK2(VectCart,ActCarts,Add_O)
     REAL(DOUBLE),DIMENSION(:)   :: VectCart
     REAL(DOUBLE),DIMENSION(:,:) :: ActCarts(:,:)
     INTEGER :: I,J,NatmsLoc
     LOGICAL,OPTIONAL :: Add_O
     !
     NatmsLoc=SIZE(ActCarts,2)
     IF(PRESENT(Add_O)) THEN
       IF(Add_O) THEN
         DO I=1,NatmsLoc
           J=3*(I-1)
           ActCarts(1,I)=ActCarts(1,I)+VectCart(J+1)
           ActCarts(2,I)=ActCarts(2,I)+VectCart(J+2)
           ActCarts(3,I)=ActCarts(3,I)+VectCart(J+3)
         ENDDO
       ENDIF
     ELSE
       DO I=1,NatmsLoc
         J=3*(I-1)
         ActCarts(1,I)=VectCart(J+1)
         ActCarts(2,I)=VectCart(J+2)
         ActCarts(3,I)=VectCart(J+3)
       ENDDO
     ENDIF
   END SUBROUTINE CartRNK1ToCartRNK2
!===============================================================
   SUBROUTINE CartRNK2ToCartRNK1(VectCart,ActCarts,Add_O)
     !
     REAL(DOUBLE),DIMENSION(:)   :: VectCart
     REAL(DOUBLE),DIMENSION(:,:) :: ActCarts(:,:)
     INTEGER :: I,J,NatmsLoc
     LOGICAL,OPTIONAL :: Add_O
     !
     NatmsLoc=SIZE(ActCarts,2)
     !
     IF(PRESENT(Add_O)) THEN
       IF(Add_O) THEN
         DO I=1,NatmsLoc
           J=3*(I-1)
           VectCart(J+1)=VectCart(J+1)+ActCarts(1,I)
           VectCart(J+2)=VectCart(J+2)+ActCarts(2,I)
           VectCart(J+3)=VectCart(J+3)+ActCarts(3,I)
         ENDDO
       ENDIF
     ELSE
       DO I=1,NatmsLoc
         J=3*(I-1)
         VectCart(J+1)=ActCarts(1,I)
         VectCart(J+2)=ActCarts(2,I)
         VectCart(J+3)=ActCarts(3,I)
       ENDDO
     ENDIF
   END SUBROUTINE CartRNK2ToCartRNK1

   SUBROUTINE Set_DBL_RNK2_EQ_DBL_RNK2(A, B)
    TYPE(DBL_RNK2), INTENT(IN)    :: B
    TYPE(DBL_RNK2), INTENT(INOUT) :: A

    IF(SIZE(A%D, 1) /= SIZE(B%D, 1)) THEN
      WRITE(*,*) "SIZE(A%D, 1) = ", SIZE(A%D, 1)
      WRITE(*,*) "SIZE(B%D, 1) = ", SIZE(B%D, 1)
      CALL Halt("size mismatch in dimension 1")
    ENDIF

    IF(SIZE(A%D, 2) /= SIZE(B%D, 2)) THEN
      WRITE(*,*) "SIZE(A%D, 2) = ", SIZE(A%D, 2)
      WRITE(*,*) "SIZE(B%D, 2) = ", SIZE(B%D, 2)
      CALL Halt("size mismatch in dimension 2")
    ENDIF

    A%D = B%D
   END SUBROUTINE Set_DBL_RNK2_EQ_DBL_RNK2

   SUBROUTINE Set_CHR10_VECT_EQ_CHR10_VECT(A, B)
    TYPE(CHR10_VECT), INTENT(IN)    :: B
    TYPE(CHR10_VECT), INTENT(INOUT) :: A

    A%C = B%C
   END SUBROUTINE Set_CHR10_VECT_EQ_CHR10_VECT

   SUBROUTINE Set_CellSet_EQ_CellSet(A, B)
    TYPE(CellSet), INTENT(IN)    :: B
    TYPE(CellSet), INTENT(INOUT) :: A

    IF(AllocQ(B%Alloc)) THEN
      A%NCells = B%NCells
      A%Radius = B%Radius
      CALL SetEq(A%CellCarts, B%CellCarts)
    ENDIF
   END SUBROUTINE Set_CellSet_EQ_CellSet

END MODULE
