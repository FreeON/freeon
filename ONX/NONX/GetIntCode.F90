  SUBROUTINE GetIntCode(LTot,TBra,TKet,IntCodeV,IntCodeC,Explicit)
    USE DerivedTypes
    USE PrettyPrint
    IMPLICIT NONE
    INTEGER           :: LTot,TBra,TKet,IntCodeV,IntCodeC
    INTEGER           :: TB,TK
    LOGICAL           :: Explicit

    IF(LTot.GT.8) THEN
      Explicit=.FALSE.
      RETURN
    ELSE
      Explicit=.TRUE.
    ENDIF
    TB=TBra
    TK=TKet

    IF(TBra.eq.0301) TB=0201
    IF(TBra.eq.0303) TB=0202
    IF(TBra.eq.0302) TB=0202
    IF(TBra.eq.0603) TB=0602
    IF(TKet.eq.0301) TK=0201
    IF(TKet.eq.0303) TK=0202
    IF(TKet.eq.0302) TK=0202
    IF(TKet.eq.0603) TK=0602

    IntCodeC = TBra*10000+TKet
    IntCodeV = TB*10000+TK

!
! L=1 exceptions
!
!    IF (IntCodeC.eq.03010101) Explicit=.FALSE.
!    IF (IntCodeC.eq.01010301) Explicit=.FALSE.
!
! L=2 exceptions
!
!    IF (IntCodeC.eq.02010301) Explicit=.FALSE.
!    IF (IntCodeC.eq.03010201) Explicit=.FALSE.
!    IF (IntCodeC.eq.03010301) Explicit=.FALSE.
!    IF (IntCodeC.eq.06010101) Explicit=.FALSE.
!    IF (IntCodeC.eq.01010601) Explicit=.FALSE.
!
! L=3 exceptions
!
!    IF (IntCodeC.eq.02020201) Explicit=.FALSE.
!    IF (IntCodeC.eq.02010202) Explicit=.FALSE.  
!    IF (IntCodeC.eq.06020101) Explicit=.FALSE.
!    IF (IntCodeC.eq.01010602) Explicit=.FALSE.
!    IF (IntCodeC.eq.06010201) Explicit=.FALSE.
!    IF (IntCodeC.eq.02010601) Explicit=.FALSE.
!    IF (IntCodeC.eq.03020301) Explicit=.FALSE.  ! 02020201
!    IF (IntCodeC.eq.03030301) Explicit=.FALSE.  ! 02020201
!    IF (IntCodeC.eq.03020201) Explicit=.FALSE.  ! 02020201
!    IF (IntCodeC.eq.03030201) Explicit=.FALSE.  ! 02020201
!    IF (IntCodeC.eq.02020301) Explicit=.FALSE.  ! 02020201
!    IF (IntCodeC.eq.03010302) Explicit=.FALSE.  ! 02010202
!    IF (IntCodeC.eq.02010303) Explicit=.FALSE.  ! 02010202
!    IF (IntCodeC.eq.03010303) Explicit=.FALSE.  ! 02010202
!    IF (IntCodeC.eq.03010202) Explicit=.FALSE.  ! 02010202
!    IF (IntCodeC.eq.02010302) Explicit=.FALSE.  ! 02010202
!    IF (IntCodeC.eq.06030101) Explicit=.FALSE.  ! 06020101
!    IF (IntCodeC.eq.01010603) Explicit=.FALSE.  ! 01010602
!    IF (IntCodeC.eq.06010301) Explicit=.FALSE.  ! 06010201
!    IF (IntCodeC.eq.03010601) Explicit=.FALSE.  ! 02010601 
!
! L=4 exceptions
!
!    IF (IntCodeC.eq.02020202) Explicit=.FALSE.
!    IF (IntCodeC.eq.06060101) Explicit=.FALSE.
!    IF (IntCodeC.eq.01010606) Explicit=.FALSE.
!    IF (IntCodeC.eq.06010601) Explicit=.FALSE.
!    IF (IntCodeC.eq.06020201) Explicit=.FALSE.
!    IF (IntCodeC.eq.02010602) Explicit=.FALSE.
!    IF (IntCodeC.eq.06010202) Explicit=.FALSE.
!    IF (IntCodeC.eq.02020601) Explicit=.FALSE.
!    IF (IntCodeC.eq.03020202) Explicit=.FALSE.  ! 02020202
!    IF (IntCodeC.eq.02020302) Explicit=.FALSE.  ! 02020202
!    IF (IntCodeC.eq.03030202) Explicit=.FALSE.  ! 02020202
!    IF (IntCodeC.eq.02020303) Explicit=.FALSE.  ! 02020202
!    IF (IntCodeC.eq.03020302) Explicit=.FALSE.  ! 02020202
!    IF (IntCodeC.eq.03030302) Explicit=.FALSE.  ! 02020202
!    IF (IntCodeC.eq.03020303) Explicit=.FALSE.  ! 02020202
!    IF (IntCodeC.eq.03030303) Explicit=.FALSE.  ! 02020202
!    IF (IntCodeC.eq.03010602) Explicit=.FALSE.  ! 02010602
!    IF (IntCodeC.eq.02010603) Explicit=.FALSE.  ! 02010602
!    IF (IntCodeC.eq.03010603) Explicit=.FALSE.  ! 02010602
!    IF (IntCodeC.eq.06020301) Explicit=.FALSE.  ! 06020201
!    IF (IntCodeC.eq.06030201) Explicit=.FALSE.  ! 06020201
!    IF (IntCodeC.eq.06030301) Explicit=.FALSE.  ! 06020201
!    IF (IntCodeC.eq.06010302) Explicit=.FALSE.  ! 06010202
!    IF (IntCodeC.eq.06010303) Explicit=.FALSE.  ! 06010202
!    IF (IntCodeC.eq.03020601) Explicit=.FALSE.  ! 02020601
!    IF (IntCodeC.eq.03030601) Explicit=.FALSE.  ! 02020601
  END SUBROUTINE GetIntCode
  
