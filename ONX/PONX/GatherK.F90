#ifdef PARALLEL
SUBROUTINE GatherK(BS,GM,NameBuf,KTotal,K)
  USE DerivedTypes
  USE GlobalScalars
  USE MemMan
  USE MondoMPI
  USE ONXParameters
  IMPLICIT NONE
  TYPE(BSET),INTENT(IN)     :: BS
  TYPE(CRDS),INTENT(IN)     :: GM
  TYPE(INT_VECT),INTENT(IN) :: NameBuf
  TYPE(BCSR),INTENT(INOUT)  :: KTotal
  TYPE(DBCSR),INTENT(IN)    :: K
  TYPE(INT_VECT)            :: AInfo
  TYPE(INT_RNK2)            :: ABuff
  TYPE(DBL_VECT)            :: MTrixB
  INTEGER                   :: I,J,NBlocks,Id,N
  INTEGER                   :: ri,ci,iPtrD,iPtrO,Ind
  INTEGER                   :: AtA,KA,NBFA
  INTEGER                   :: AtB,KB,NBFB,N2
  CALL New(ABuff,(/4,NElem/))
  CALL New(AInfo,3)
!--------------------------------------------------------------------------------
! Do the local copy first
!--------------------------------------------------------------------------------
  IF (MyID==ROOT) THEN
    DO ri=1,NRows
      AtA=NameBuf%I(ri)
      KA=GM%AtTyp%I(AtA)
      NBFA=BS%BfKnd%I(KA)
      DO ci=K%RowPt%I(ri),K%RowPt%I(ri+1)-1
        AtB=K%ColPt%I(ci)
        iPtrO=K%BlkPt%I(ci)
        KB=GM%AtTyp%I(AtB)
        NBFB=BS%BfKnd%I(KB)
        N2=NBFA*NBFB
        CALL GetAdrB(AtA,AtB,Ind,0,KTotal%ColPt%I,KTotal%RowPt%I)
        iPtrD=KTotal%BlkPt%I(Ind)
        CALL DBL_VECT_EQ_DBL_VECT(N2,KTotal%MTrix%D(iPtrD),K%MTrix%D(iPtrO))
      END DO ! ci
    END DO ! ri
  END IF 
!--------------------------------------------------------------------------------
! Communicate with the nodes
!--------------------------------------------------------------------------------
  IF (MyID==ROOT) THEN
    CALL New(MTrixB,MaxN2)
    DO I=1,NPrc-1
      CALL Recv(NBlocks,I,0)
      DO J=1,NBlocks
        Id=J+NBlocks
        CALL Recv(AInfo,3,I,J)
        AtA=AInfo%I(1)
        AtB=AInfo%I(2)
        N2 =AInfo%I(3)
        CALL GetAdrB(AtA,AtB,Ind,0,KTotal%ColPt%I,KTotal%RowPt%I)
        iPtrD=KTotal%BlkPt%I(Ind)
        CALL Recv(MTrixB,N2,I,Id)
        CALL DBL_VECT_PLS_DBL_VECT(N2,KTotal%MTrix%D(iPtrD),MTrixB%D)
      END DO ! J
    END DO ! I
    CALL Delete(MTrixB)
  ELSE
    NBlocks=0
    DO ri=1,NRows
      AtA=NameBuf%I(ri)
      KA=GM%AtTyp%I(AtA)
      NBFA=BS%BfKnd%I(KA)
      DO ci=K%RowPt%I(ri),K%RowPt%I(ri+1)-1
        AtB=K%ColPt%I(ci)
        iPtrO=K%BlkPt%I(ci)
        KB=GM%AtTyp%I(AtB)
        NBFB=BS%BfKnd%I(KB)
        N2=NBFA*NBFB
        NBlocks=NBlocks+1
        ABuff%I(1,NBlocks)=AtA
        ABuff%I(2,NBlocks)=AtB
        ABuff%I(3,NBlocks)=N2
        ABuff%I(4,NBlocks)=iPtrO
      END DO ! ci
    END DO ! ri
    CALL Send(NBlocks,ROOT,0)
    DO J=1,NBlocks
      Id=J+NBlocks
      AInfo%I(1)=ABuff%I(1,J)
      AInfo%I(2)=ABuff%I(2,J)
      AInfo%I(3)=ABuff%I(3,J)
      iPtrO=ABuff%I(4,J)
      N=AInfo%I(3)+iPtrO-1
      CALL Send(AInfo,3,ROOT,J)
      CALL Send(K%MTrix,N,ROOT,Id,iPtrO)
    END DO ! J
  END IF ! (MyID==ROOT)
  CALL Delete(AInfo)
  CALL Delete(ABuff)
END SUBROUTINE GatherK
#endif
! Added so there is some source in here, otherwise will
! break some compilers (eg. SGI).  Probably should put 
! in a module...
SUBROUTINE DUMMY1
END SUBROUTINE DUMMY1
