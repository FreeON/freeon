  SUBROUTINE VRRs(LBra,LKet,Drv)
    USE DerivedTypes
    USE InOut
    IMPLICIT NONE
    TYPE(IDrv),INTENT(INOUT)   :: Drv
    INTEGER,INTENT(IN)         :: LBra,LKet
    INTEGER                    :: VRRSpace
    INTEGER                    :: LMAX=4      ! Maximum distribution angular
                                              ! symmetry that is included in 
                                              ! the VRR driver files.
    IF (LBra>LMAX.OR.LKet>LMAX) THEN
      CALL Halt(' LBra or LKet too large in ONX:' &
              //' LBra = '//TRIM(IntToChar(LBra)) &
              //' LKet = '//TRIM(IntToChar(LKet)))
    END IF
    Drv%id = LBra+1+LKet*(LMAX+1)
    Drv%is = Drv%SLoc%I(1,Drv%id)
    Drv%nr = Drv%SLoc%I(2,Drv%id)
    Drv%ns = Drv%SLoc%I(3,Drv%id)
  END SUBROUTINE VRRs
