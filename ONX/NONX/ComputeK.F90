SUBROUTINE ComputeK(BSc,GMc,BSp,GMp,D,K,DB,IB,MB,Drv)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE Thresholding
  USE ONXParameters
  USE ONXMemory
  USE Stats
  USE GetTables
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(DBCSR)           :: D
  TYPE(DBCSR)           :: K
  TYPE(BCSR)            :: KTotal
#else
  TYPE(BCSR)            :: D
  TYPE(BCSR)            :: K
  INTEGER               :: MyID=0
#endif
  TYPE(BSET),INTENT(IN) :: BSc,BSp   ! basis set info
  TYPE(CRDS),INTENT(IN) :: GMc,GMp   ! geometry info
  TYPE(DBuf)            :: DB        ! ONX distribution buffers
  TYPE(IBuf)            :: IB        ! ONX 2-e eval buffers
  TYPE(DML)             :: MB        ! ONX distribution pointers
  TYPE(IDrv)            :: Drv       ! VRR/contraction drivers
!-------------------------------------------------------------------
! Misc. internal variables
!-------------------------------------------------------------------
  INTEGER               :: iTBra,TBra,LBra,ITypeC,ITypeA,NA,NC
  INTEGER               :: iTKet,TKet,LKet,ITypeD,ITypeB,NB,ND
!-------------------------------------------------------------------
! Function calls
!-------------------------------------------------------------------
  INTEGER               :: LTotal,NFinal

!  CALL GetGammaTable(0,IB)
  CALL GetExpTable(IB)


  DO iTBra=1,DB%NTypes       ! Loop over angular symmetry types on the Bra
    TBra=DB%TCode%I(iTBra)
    write(*,*) "TBra=",TBra
    ITypeC=MOD(TBra,100)
    ITypeA=(TBra-ITypeC)/100
    write(*,*) "ITypeA=",ITypeA
    write(*,*) "ITypeC=",ITypeC
    write(*,*) "--------------"
    LBra=LTotal(ITypeA)+LTotal(ITypeC)
    NA=NFinal(ITypeA)
    NC=NFinal(ITypeC)

    DO iTKet=1,DB%NTypes     ! Loop over angular symmetry types on the Ket


    END DO ! iTKet
  END DO ! iTBra

  CALL Halt('enough for now')

END SUBROUTINE ComputeK
