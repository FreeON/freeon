  FUNCTION NFinal(IType)
    USE InOut
    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: IType
    INTEGER             :: NFinal
    INTEGER, PARAMETER  :: NT=28
    INTEGER, PARAMETER  :: T2Size(NT) = (/ 1,                  &
                                           4, 3,               &
                                          10, 9, 6,            &
                                          20,19,16,10,         &
                                          35,34,31,25,15,      &
                                          56,55,52,46,36,21,   &
                                          84,83,80,74,64,49,28 /)
    IF (IType>0.AND.IType<=NT) THEN
      NFinal=T2Size(IType)
    ELSE
      CALL Halt(' Illegal IType in NFinal')
    END IF
  END FUNCTION NFinal

  FUNCTION iT(I,J)
    USE InOut
    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: I,J
    INTEGER             :: iT,Ind
    INTEGER, PARAMETER  :: NT=6
    INTEGER, PARAMETER  :: iTD(NT*NT) = (/ 1, 0, 0, 0, 0, 0,  &
                                           2, 3, 0, 0, 0, 0,  &
                                           4, 5, 6, 0, 0, 0,  &
                                           0, 0, 0, 0, 0, 0,  &
                                           0, 0, 0, 0, 0, 0,  &
                                           7, 8, 9, 0, 0, 10  /)
    Ind=I+(J-1)*NT
    IF (Ind>0.AND.Ind<=NT*NT) THEN
      iT=iTD(Ind)
    ELSE
      CALL Halt(' Illegal Ind in iT')
    END IF
  END FUNCTION iT

  FUNCTION LTotal(I)
    USE InOut
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: I
    INTEGER              :: LTotal
    INTEGER,PARAMETER    :: Lang(21) = (/ 0,            &
                                          1,1,          &
                                          2,2,2,        &
                                          3,3,3,3,      &
                                          4,4,4,4,4,    &
                                          5,5,5,5,5,5   /)
    IF (I<=21.AND.I>0) THEN
      LTotal = Lang(I)
    ELSE
      WRITE(*,*) ' iType = ',I
      CALL Halt('Illegal iType in Ltotal')
    END IF
  END FUNCTION LTotal

  FUNCTION VRRSpace(LBra,LKet)
    USE InOut
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: LBra,LKet
    INTEGER            :: VRRSpace,Ind
    INTEGER,PARAMETER  :: VSpace(25) = (/   1,   4,   17,   24,   48, &
                                            4,  16,   68,   97,  197, &
                                           17,  68,  289,  408,  816, &
                                           24,  80,  408,  493,  983, &
                                           48, 140,  816,  863, 1708  /)
    IF(LBra>4.OR.LKet>4) THEN
      CALL Halt('Illegal LBra or LKet in VRRSpace')
    ENDIF
    Ind = (LBra+1)+LKet*5
    VRRSpace=VSpace(Ind)
  END FUNCTION VRRSpace
