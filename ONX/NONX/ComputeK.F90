SUBROUTINE ComputeK(BSc,GMc,BSp,GMp,D,K,DB,IB,MB,Drv)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE Thresholding
  USE ONXParameters
  USE ONXMemory
  USE Stats
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(DBCSR)         :: D
  TYPE(DBCSR)         :: K
  TYPE(BCSR)          :: KTotal
#else
  TYPE(BCSR)          :: D
  TYPE(BCSR)          :: K
  INTEGER             :: MyID=0
#endif
  TYPE(BSET),INTENT(IN)    :: BSc,BSp   ! basis set info
  TYPE(CRDS),INTENT(IN)    :: GMc,GMp   ! geometry info
  TYPE(DBuf)               :: DB        ! ONX distribution buffers
  TYPE(IBuf)               :: IB        ! ONX 2-e eval buffers
  TYPE(DML)                :: MB        ! ONX distribution pointers
  TYPE(IDrv)               :: Drv       ! VRR/contraction drivers

  CALL GetGammaTable(0,IB)
  CALL GetExpTable(IB)

  CALL Halt('enough for now')

END SUBROUTINE ComputeK
