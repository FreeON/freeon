MODULE MatFilter
  USE DerivedTypes
  USE GlobalScalars
  USE MemMan
  USE LinAlg
  USE ONXParameters
  INTERFACE ONXFilter
#ifdef PARALLEL_ONX
    MODULE PROCEDURE ONXFilter_DBCSR 
#endif
    MODULE PROCEDURE ONXFilter_BCSR
  END INTERFACE 

  CONTAINS

#ifdef PARALLEL_ONX
SUBROUTINE ONXFilter_DBCSR(BS,GM,A,NameBuf,Tol_O)
  IMPLICIT NONE
  TYPE(BSET),INTENT(IN)             :: BS
  TYPE(CRDS),INTENT(IN)             :: GM
  TYPE(DBCSR),INTENT(INOUT)         :: A
  TYPE(INT_VECT),INTENT(IN)         :: NameBuf
  REAL(DOUBLE),OPTIONAL,INTENT(IN)  :: Tol_O
  INTEGER                           :: AtA,KA,NBFA
  INTEGER                           :: AtB,KB,NBFB,N2
  INTEGER                           :: K,K1,K2,Ind1,Ind2,iPtr,I
  REAL(DOUBLE)                      :: FN,Tol
  REAL(DOUBLE),EXTERNAL             :: DBL_Dot
  IF (PRESENT(Tol_O)) THEN
    Tol=Tol_O
  ELSE
    Tol=Thresholds%Trix
  END IF
  Ind1=1
  Ind2=1
  DO I=1,NRows
    AtA=NameBuf%I(I)
    KA=GM%AtTyp%I(AtA)
    NBFA=BS%BfKnd%I(KA)
    K1=A%RowPt%I(AtA)
    K2=A%RowPt%I(AtA+1)-1
    A%RowPt%I(AtA)=Ind2
    DO K=K1,K2
      AtB=A%ColPt%I(K)
      KB=GM%AtTyp%I(AtB)
      NBFB=BS%BfKnd%I(KB)
      N2=NBFA*NBFB
      iPtr=A%BlkPt%I(K)
      FN=SQRT(DBL_Dot(N2,A%MTrix%D(iPtr:iPtr+N2-1),A%MTrix%D(iPtr:iPtr+N2-1)))
!      FN=FNorm(N2,A%MTrix%D(iPtr:iPtr+N2-1))
      IF (FN.GT.Tol) THEN
        CALL DBL_VECT_EQ_DBL_VECT(N2,A%MTrix%D(Ind1),A%MTrix%D(iPtr))
        A%ColPt%I(Ind2)=AtB
        A%BlkPt%I(Ind2)=Ind1
        Ind1=Ind1+N2
        Ind2=Ind2+1
      END IF
    END DO ! K
  END DO ! NAtoms
  A%RowPt%I(NAtoms+1)=Ind2
  A%NBlks=Ind2-1
  A%NNon0=Ind1-1
END SUBROUTINE ONXFilter_DBCSR
#endif

SUBROUTINE ONXFilter_BCSR(BS,GM,A,NameBuf,Tol_O)
  IMPLICIT NONE
  TYPE(BSET),INTENT(IN)             :: BS
  TYPE(CRDS),INTENT(IN)             :: GM
  TYPE(BCSR),INTENT(INOUT)          :: A
  TYPE(INT_VECT),INTENT(IN)         :: NameBuf
  REAL(DOUBLE),OPTIONAL,INTENT(IN)  :: Tol_O
  INTEGER                           :: AtA,KA,NBFA
  INTEGER                           :: AtB,KB,NBFB,N2
  INTEGER                           :: K,K1,K2,Ind1,Ind2,iPtr
  REAL(DOUBLE)                      :: FN,Tol
  REAL(DOUBLE),EXTERNAL             :: DBL_Dot
  IF (PRESENT(Tol_O)) THEN
    Tol=Tol_O
  ELSE
    Tol=Thresholds%Trix
  END IF
  Ind1=1
  Ind2=1
  DO AtA=1,NAtoms
    KA=GM%AtTyp%I(AtA)
    NBFA=BS%BfKnd%I(KA)
    K1=A%RowPt%I(AtA)
    K2=A%RowPt%I(AtA+1)-1
    A%RowPt%I(AtA)=Ind2
    DO K=K1,K2
      AtB=A%ColPt%I(K)
      KB=GM%AtTyp%I(AtB)
      NBFB=BS%BfKnd%I(KB)
      N2=NBFA*NBFB
      iPtr=A%BlkPt%I(K)
      FN=SQRT(DBL_Dot(N2,A%MTrix%D(iPtr:iPtr+N2-1),A%MTrix%D(iPtr:iPtr+N2-1)))
!      FN=FNorm(N2,A%MTrix%D(iPtr:iPtr+N2-1))
      IF (FN.GT.Tol) THEN
        CALL DBL_VECT_EQ_DBL_VECT(N2,A%MTrix%D(Ind1),A%MTrix%D(iPtr))
        A%ColPt%I(Ind2)=AtB
        A%BlkPt%I(Ind2)=Ind1
        Ind1=Ind1+N2
        Ind2=Ind2+1
      END IF
    END DO ! K
  END DO ! NAtoms
  A%RowPt%I(NAtoms+1)=Ind2
  A%NBlks=Ind2-1
  A%NNon0=Ind1-1
END SUBROUTINE ONXFilter_BCSR

END MODULE MatFilter
