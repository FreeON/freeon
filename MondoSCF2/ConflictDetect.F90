MODULE Conflicted
  USE DerivedTypes
  USE Parse
  USE ProcessControl
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






  SUBROUTINE GeoConflictCheck(G)
!H---------------------------------------------------------------------------------
!H SUBROUTINE GeoConflictCheck(G)
!H  Main driver for geometry conflict.
!H---------------------------------------------------------------------------------
    TYPE(Geometries) :: G
    !-------------------------------------------------------------------
    !
    CALL GeoConflictCheck1(G)
    !CALL GeoConflictCheck2()
    !CALL GeoConflictCheck3()
    !CALL GeoConflictCheck4()
    !
  END SUBROUTINE GeoConflictCheck
  !
  !
  SUBROUTINE GeoConflictCheck1(G)
!H---------------------------------------------------------------------------------
!H SUBROUTINE GeoConflictCheck1(G)
!H  This routine checks if a pair of atoms are too close from each other.
!H---------------------------------------------------------------------------------
    TYPE(Geometries)               :: G
    !-------------------------------------------------------------------
    INTEGER                        :: i,iClone,AtA,AtB,MaxClone
    INTEGER, DIMENSION(2)          :: NClone
    REAL(DOUBLE)                   :: Ax,Ay,Az,Bx,By,Bz,Dist,t1,t2
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: Text
    !-------------------------------------------------------------------
    REAL(DOUBLE), PARAMETER        :: MinDist = 0.5D+00
    !-------------------------------------------------------------------
    !
    !Check for clones.
    IF(G%Clones.EQ.1) THEN
       ! We have only clone.
       MaxClone=1
       NClone(1)=1
       NClone(2)=1
    ELSE
       ! We have several clones, check pairs 
       ! only for reactant and product.
       MaxClone=2
       NClone(1)=0             ! reactant.
       NClone(2)=G%Clones+1    ! product.
    ENDIF
    !
    DO i=1,MaxClone
       iClone=NClone(i)
       DO AtA=2,G%Clone(iClone)%NAtms
          Ax=G%Clone(iClone)%AbCarts%D(1,AtA)
          Ay=G%Clone(iClone)%AbCarts%D(2,AtA)
          Az=G%Clone(iClone)%AbCarts%D(3,AtA)
          DO AtB=1,AtA-1
             Bx=G%Clone(iClone)%AbCarts%D(1,AtB)
             By=G%Clone(iClone)%AbCarts%D(2,AtB)
             Bz=G%Clone(iClone)%AbCarts%D(3,AtB)
             Dist=SQRT((Ax-Bx)**2+(Ay-By)**2+(Az-Bz)**2)
             !write(*,*) 'Dist',Dist,Ax,Ay,Az,Bx,By,Bz
             IF(Dist.LE.MinDist) THEN
                Text=' The distance between the atoms '// &
                     & TRIM(IntToChar(AtB))//' '//TRIM(G%Clone(iClone)%AtNam%C(AtB))//' and '// &
                     & TRIM(IntToChar(AtA))//' '//TRIM(G%Clone(iClone)%AtNam%C(AtA))// &
                     & ' in the clone '//TRIM(IntToChar(iClone)) //' is '// &
                     & TRIM(DblToShrtChar(Dist))//' !'
                CALL MondoHalt(PRSE_ERROR,Text)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    !
  END SUBROUTINE GeoConflictCheck1


END MODULE Conflicted
