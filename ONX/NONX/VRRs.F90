  SUBROUTINE VRRs(LBra,LKet,id,is,nr,ns,SLoc)
    USE DerivedTypes
    USE InOut
    IMPLICIT NONE
    TYPE(INT_RNK2),INTENT(IN)  :: SLoc
    INTEGER,INTENT(IN)         :: LBra,LKet
    INTEGER,INTENT(OUT)        :: id,is,nr,ns
    INTEGER                    :: LMAX=4      ! Maximum distribution angular
                                              ! symmetry that is included in 
                                              ! the VRR driver files.
    IF (LBra>LMAX.OR.LKet>LMAX) THEN
      CALL Halt(' LBra or LKet too large in ONX')
    END IF
    id = LBra+1+LKet*(LMAX+1)
    is = SLoc%I(1,id)
    nr = SLoc%I(2,id)
    ns = SLoc%I(3,id)
  END SUBROUTINE VRRs
