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
    
    IF (IntCodeC.eq.02010301) Explicit=.FALSE.
    IF (IntCodeC.eq.03010201) Explicit=.FALSE.
    IF (IntCodeC.eq.03010301) Explicit=.FALSE.
    IF (IntCodeC.eq.01010601) Explicit=.FALSE.
    IF (IntCodeC.eq.06010101) Explicit=.FALSE.

  END SUBROUTINE GetIntCode
  
