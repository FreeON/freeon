MODULE Conflicted
  USE DerivedTypes
  USE ControlStructures
  IMPLICIT NONE
CONTAINS 
  SUBROUTINE ConflictCheck()
    WRITE(*,*)' NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW '
    WRITE(*,*) 
    WRITE(*,*)' My name is ConflictCheck, I am a new subroutine to check logical conficts '
    WRITE(*,*)' after parsing!  I am in MondoSCF2/ConflictDetect.F90. Add to me as you    '
    WRITE(*,*)' encounter conflicts that should be resolved at startup.'
    WRITE(*,*) 
    WRITE(*,*)' NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW '

!    CALL ConflictCheck1()
!    CALL ConflictCheck2()
!    CALL ConflictCheck3()
!    CALL ConflictCheck4()

  END SUBROUTINE ConflictCheck

#ifdef LJDFLSJDFLSDJFLSKJFDLSDJF
  SUBROUTINE ConflictCheck1(O,D)
    ! Is MondoSCF compiled correctly for parallel replicas?
    IF(O%GradOpt==GRAD_DYNAMICS.AND.D%MDAlgorithm==MD_PARALLEL_REP)THEN
#ifdef !defined(PARALLEL)
       CALL MondoHalt(PRSE_ERROR,' MondoSCF must be compiled in parallel for replica exchange to be active.')
#endif
    ENDIF
  END SUBROUTINE ConflictCheck1

  SUBROUTINE ConflictCheck2(O,D,B)
    ! Check that the number of options in the progression of models/methods/accuracies/basissets 
    ! are consistent with each other 
    N=0
    N=MAX(N,O%NModls)
    N=MAX(N,O%NMthds)
    N=MAX(N,O%NThrsh)
    N=MAX(N,B%NBSets)
    IF(N/=O%NModls) &
       CALL MondoHalt(PRSE_ERROR,' Model chemistries in sequence is short.')
    IF(N/=O%NMthds) &
       CALL MondoHalt(PRSE_ERROR,' SCF methods in sequence is short.')
    IF(N/=O%NThrsh) &
       CALL MondoHalt(PRSE_ERROR,' Accuracies in sequence is short.')
    IF(N/=B%NBSets) &
       CALL MondoHalt(PRSE_ERROR,' Basis sets in sequence is short.')
  END SUBROUTINE ConflictCheck2
#endif

END MODULE Conflicted
