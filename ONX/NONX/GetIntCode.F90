  SUBROUTINE GetIntCode(LTot,TBra,TKet,IntCodeV,IntCodeC,Explicit)
    USE DerivedTypes
    IMPLICIT NONE
    INTEGER           :: LTot,TBra,TKet,IntCodeV,IntCodeC
    INTEGER           :: TB,TK
    LOGICAL           :: Explicit

    IF(LTot.GT.2) THEN
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
    IF(TKet.eq.0301) TK=0201
    IF(TKet.eq.0303) TK=0202
    IF(TKet.eq.0302) TK=0202
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
    IF (IntCodeC.eq.02020201) Explicit=.FALSE.
    IF (IntCodeC.eq.02010202) Explicit=.FALSE.
!
! L=4 exceptions
!
    IF (IntCodeC.eq.02020202) Explicit=.FALSE.
    IF (IntCodeC.eq.03020202) Explicit=.FALSE.
    IF (IntCodeC.eq.02020302) Explicit=.FALSE.
    IF (IntCodeC.eq.03020302) Explicit=.FALSE.
    IF (IntCodeC.eq.03030202) Explicit=.FALSE.
    IF (IntCodeC.eq.02020303) Explicit=.FALSE.
    IF (IntCodeC.eq.03030302) Explicit=.FALSE.
    IF (IntCodeC.eq.03020303) Explicit=.FALSE.
    IF (IntCodeC.eq.03030303) Explicit=.FALSE.

  END SUBROUTINE GetIntCode
  
