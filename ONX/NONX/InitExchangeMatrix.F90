MODULE InitExchangeMatrix
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE PrettyPrint
  USE ONXParameters
  IMPLICIT NONE

#ifdef PARALLEL
  INTERFACE InitK
    MODULE PROCEDURE InitK_DBCSR, InitK_BCSR
  END INTERFACE
#else
  INTERFACE InitK
    MODULE PROCEDURE InitK_BCSR
  END INTERFACE
#endif

  CONTAINS

SUBROUTINE InitK_BCSR(BS,GM,K,NameBuf)
  TYPE(BSET),INTENT(IN)        :: BS
  TYPE(CRDS),INTENT(IN)        :: GM
  TYPE(BCSR),INTENT(INOUT)     :: K
  TYPE(INT_VECT),INTENT(IN)    :: NameBuf
  INTEGER                      :: AtA,KA,NBFA,CFA
  INTEGER                      :: AtB,KB,NBFB
  INTEGER                      :: StartLA,StopLA,StrideA
  INTEGER                      :: IndexA1,IndexA2
  REAL(DOUBLE)                 :: Ax,Ay,Az,AB2
  REAL(DOUBLE)                 :: Bx,By,Bz
  INTEGER                      :: J,iPnt,N2,NBFM

  K%MTrix%D(1:NElem)=0.0D0
  j=1; iPnt=1
  DO AtA=1,NAtoms
    KA=GM%AtTyp%I(AtA)
    NBFA=BS%BfKnd%I(KA)
    Ax=GM%Carts%D(1,AtA)
    Ay=GM%Carts%D(2,AtA)
    Az=GM%Carts%D(3,AtA)
    K%RowPt%I(AtA)=j
    DO AtB=1,NAtoms
      KB=GM%AtTyp%I(AtB)
      NBFB=BS%BfKnd%I(KB)
      Bx=GM%Carts%D(1,AtB)
      By=GM%Carts%D(2,AtB)
      Bz=GM%Carts%D(3,AtB)
      AB2=(Ax-Bx)*(Ax-Bx) + &
          (Ay-By)*(Ay-By) + &
          (Az-Bz)*(Az-Bz)
      IF (SQRT(AB2).LE.ONXRange) THEN
        K%ColPt%I(j)=AtB
        K%BlkPt%I(j)=iPnt
        j=j+1
        iPnt=iPnt+NBFA*NBFB
      END IF
    END DO ! ci
  END DO ! ri
  K%RowPt%I(NRows+1)=j

END SUBROUTINE InitK_BCSR


#ifdef PARALLEL

SUBROUTINE InitK_DBCSR(BS,GM,K,NameBuf)
  TYPE(BSET),INTENT(IN)        :: BS
  TYPE(CRDS),INTENT(IN)        :: GM
  TYPE(DBCSR),INTENT(INOUT)    :: K
  TYPE(INT_VECT),INTENT(IN)    :: NameBuf
  INTEGER                      :: j,iPnt,ri,ci
  INTEGER                      :: IndexA1,IndexA2
  INTEGER                      :: AtA,KA,NBFA,CFA,StartLA,StopLA,StrideA
  INTEGER                      :: AtB,KB,NBFB
  REAL(DOUBLE)                 :: Ax,Ay,Az,Bx,By,Bz,AB2

  K%MTrix%D(1:NElem)=0.0D0
  j=1; iPnt=1
  DO ri=1,NRows
    AtA=NameBuf%I(ri)
    KA=GM%AtTyp%I(AtA)
    NBFA=BS%BfKnd%I(KA)
    Ax=GM%Carts%D(1,AtA)
    Ay=GM%Carts%D(2,AtA)
    Az=GM%Carts%D(3,AtA)
    K%RowPt%I(ri)=j
    K%GRwPt%I(AtA)=j
    DO ci=1,NRows
      AtB=NameBuf%I(ci)
      KB=GM%AtTyp%I(AtB)
      NBFB=BS%BfKnd%I(KB)
      Bx=GM%Carts%D(1,AtB)
      By=GM%Carts%D(2,AtB)
      Bz=GM%Carts%D(3,AtB)
      AB2=(Ax-Bx)*(Ax-Bx) + &
          (Ay-By)*(Ay-By) + &
          (Az-Bz)*(Az-Bz)
      IF (SQRT(AB2).LE.ONXRange) THEN
        K%ColPt%I(j)=AtB
        K%BlkPt%I(j)=iPnt
        j=j+1
        iPnt=iPnt+NBFA*NBFB
      END IF
    END DO ! ci
    K%GRwPt%I(AtA+1)=j
  END DO ! ri
  K%RowPt%I(NRows+1)=j  
  
END SUBROUTINE InitK_DBCSR
#endif

END MODULE
