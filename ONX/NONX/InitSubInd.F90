  SUBROUTINE InitSubInd(BS,GM,SubInd)
  USE DerivedTypes
  IMPLICIT NONE
  TYPE(BSET),INTENT(IN)        :: BS
  TYPE(CRDS),INTENT(IN)        :: GM
  TYPE(INT_RNK2),INTENT(INOUT) :: SubInd
  INTEGER                      :: IndexA1,IndexA2
  INTEGER                      :: AtA,KA,NBFA,CFA
  INTEGER                      :: StartLA,StopLA,StrideA
  IndexA1=0
  DO AtA=1,NAtoms
    KA=GM%AtTyp%I(AtA)
    NBFA=BS%BfKnd%I(KA)
    IndexA2=0
    DO CFA=1,BS%NCFnc%I(KA)
      IndexA1=IndexA1+1
      StartLA=BS%LStrt%I(CFA,KA)
      StopLA=BS%LStop%I(CFA,KA)
      StrideA=StopLA-StartLA+1
      SubInd%I(1,IndexA1)=AtA
      SubInd%I(2,IndexA1)=NBFA
      SubInd%I(3,IndexA1)=IndexA2+1
      IndexA2=IndexA2+StrideA
    END DO 
  END DO 
END SUBROUTINE InitSubInd
