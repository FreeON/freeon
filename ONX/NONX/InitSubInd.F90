SUBROUTINE InitSubInd(BS,GM,SubInd)
  USE DerivedTypes
!#ifdef PARALLEL_ONX          !vw
!  USE MondoMPI               !vw
!#endif                       !vw
  IMPLICIT NONE
  TYPE(BSET),INTENT(IN)        :: BS
  TYPE(CRDS),INTENT(IN)        :: GM
  TYPE(INT_RNK2),INTENT(INOUT) :: SubInd
  INTEGER                      :: IndexA1,IndexA2
  INTEGER                      :: AtA,KA,NBFA,CFA
  INTEGER                      :: StartLA,StopLA,StrideA
!#ifdef PARALLEL_ONX          !vw
!   INTEGER :: ri             !vw
!#endif                       !vw
!#ifdef PARALLEL_ONX             !vw
! DO AtA=Beg%I(MyID),End%I(MyID) !vw THESE WILL GIVE LOCAL ARRAYS
!    ri=AtA-Beg%I(MyID)+1        !vw THESE WILL GIVE LOCAL ARRAYS
!#else                           !vw

  IndexA1=0
  DO AtA=1,NAtoms
!#endif                          !vw
     KA=GM%AtTyp%I(AtA)
     NBFA=BS%BfKnd%I(KA)
     IndexA2=0
     DO CFA=1,BS%NCFnc%I(KA)
        IndexA1=IndexA1+1
        StartLA=BS%LStrt%I(CFA,KA)
        StopLA=BS%LStop%I(CFA,KA)
        StrideA=StopLA-StartLA+1
!#ifdef PARALLEL_ONX            !vw
!        SubInd%I(1,IndexA1)=ri !vw THESE WILL GIVE LOCAL ARRAYS
!#else                          !vw
        SubInd%I(1,IndexA1)=AtA              !vw THESE ARE GLOBAL ARRAYS
!#endif                         !vw
        SubInd%I(2,IndexA1)=NBFA             !vw THESE ARE GLOBAL ARRAYS
        SubInd%I(3,IndexA1)=IndexA2+1        !vw THESE ARE GLOBAL ARRAYS
        !write(*,*) 'MyID',MyID,SubInd%I(1,IndexA1),SubInd%I(2,IndexA1),SubInd%I(3,IndexA1)
        IndexA2=IndexA2+StrideA
     END DO
  END DO
END SUBROUTINE InitSubInd
