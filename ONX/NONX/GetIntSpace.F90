  SUBROUTINE GetIntSpace(TBra,TKet,LBra,LKet,IS)
    USE DerivedTypes
    USE PrettyPrint
    IMPLICIT NONE
    INTEGER           :: TBra,TKet
    INTEGER           :: LBra,LKet
    TYPE(ISpc)        :: IS
    INTEGER           :: IndV,IndAC,IndBD
    INTEGER           :: ITypeA,ITypeB,ITypeC,ITypeD
    INTEGER,PARAMETER :: Space1(0:4) = (/1,4,17,24,36/)
    INTEGER,PARAMETER :: Space2(36)  = (/ 1,   4,    3,    0,    0,    6,  &
                                          4,  17,   12,    0,    0,   22,  &
                                          3,  12,    9,    0,    0,   16,  &
                                          0,   0,    0,    0,    0,    0,  &
                                          0,   0,    0,    0,    0,    0,  &
                                          6,  22,   16,    0,    0,   31   /)
    INTEGER,PARAMETER :: SpaceV(25)  = (/ 1,   4,   17,   24,   48, &
                                          4,  16,   68,   97,  197, &
                                         17,  68,  289,  408,  816, &
                                         24,  88,  408,  526,  983, &
                                         48, 140,  816,  863, 1708  /)
    INTEGER,PARAMETER :: SpaceF(28)  = (/ 1,                  &
                                          4, 3,               &
                                         10, 9, 6,            &
                                         20,19,16,10,         &
                                         35,34,31,25,15,      &
                                         56,55,52,46,36,21,   &
                                         84,83,80,74,64,49,28 /)

    IF(LBra>4.OR.LKet>4) THEN
      CALL Halt('Illegal LBra or LKet in GetIntSpace')
    ENDIF

    ITypeA = Mod(TBra,100)
    ITypeC = (TBra-ITypeA)/100
    ITypeB=MOD(TKet,100)
    ITypeD=(TKet-ITypeB)/100
    IndV = (LBra+1)+LKet*5
    IndAC = ITypeA+(ITypeC-1)*6
    IndBD = ITypeB+(ITypeD-1)*6

    IS%NB1  = Space1(LBra)
    IS%NK1  = Space1(LKet)
    IS%NB2  = Space2(IndAC)
    IS%NK2  = Space2(IndBD)
    IS%L1   = SpaceF(ITypeC)
    IS%L2   = SpaceF(ITypeA)
    IS%L3   = SpaceF(ITypeD)
    IS%L4   = SpaceF(ITypeB)
    IS%NVRR = SpaceV(IndV)

  END SUBROUTINE GetIntSpace
  
