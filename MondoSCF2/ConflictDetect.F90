MODULE Conflicted
  USE DerivedTypes
  USE Parse
  USE ProcessControl
  USE ControlStructures
  IMPLICIT NONE
  !
  PRIVATE :: GlbConflictCheck1
  PRIVATE :: GeoConflictCheck1
  !
CONTAINS 
  !
  !
  SUBROUTINE OptConflictCheck(O)
!H---------------------------------------------------------------------------------
!H SUBROUTINE OptConflictCheck(O)
!H  Main driver for options conflict.
!H---------------------------------------------------------------------------------
    TYPE(Options)   :: O
    !-------------------------------------------------------------------
    !
    WRITE(*,*)' NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW '
    WRITE(*,*) 
    WRITE(*,*)' My name is ConflictCheck, I am a new subroutine to check logical conficts '
    WRITE(*,*)' after parsing!  I am in MondoSCF2/ConflictDetect.F90. Add to me as you    '
    WRITE(*,*)' encounter conflicts that should be resolved at startup.'
    WRITE(*,*) 
    WRITE(*,*)' NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW '
    ! CALL OptConflictCheck1()
    ! CALL OptConflictCheck2()
    ! CALL OptConflictCheck3()
    ! CALL OptConflictCheck4()
    !
  END SUBROUTINE OptConflictCheck
  !
  !
  SUBROUTINE GlbConflictCheck(C)
!H---------------------------------------------------------------------------------
!H SUBROUTINE GlbConflictCheck(C)
!H  Main driver for options conflict.
!H---------------------------------------------------------------------------------
    TYPE(Controls) :: C
    !-------------------------------------------------------------------
    !
    CALL GlbConflictCheck1(C)
    !CALL GlbConflictCheck2()
    !CALL GlbConflictCheck3()
    !CALL GlbConflictCheck4()
    !
  END SUBROUTINE GlbConflictCheck
  !
  !
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
  SUBROUTINE GlbConflictCheck1(C)
!H---------------------------------------------------------------------------------
!H SUBROUTINE GlbConflictCheck1(C)
!H  Check that the number of options in the progression of 
!H  models/methods/accuracies/basissets are consistent with each other 
!H---------------------------------------------------------------------------------
    TYPE(Controls) :: C
    !-------------------------------------------------------------------
    INTEGER        :: N
    !-------------------------------------------------------------------
    !
    N=0
    N=MAX(N,C%Opts%NModls)
    N=MAX(N,C%Opts%NMthds)
    N=MAX(N,C%Opts%NThrsh)
    N=MAX(N,C%Sets%NBSets)
    IF(N/=C%Opts%NModls) &
         & CALL MondoHalt(PRSE_ERROR,' Model chemistries in sequence is short.')
    IF(N/=C%Opts%NMthds) &
         & CALL MondoHalt(PRSE_ERROR,' SCF methods in sequence is short.')
    !
    !vw PROBLEM WITH THIS TEST, CAUSE WE NEED THE THRESHOLD ARRAY TO SET THE
    !vw BSCR DIMENSIONS IN LoadBasisSets! THE TEST MAY BE DONE IN LoadBasisSets
    !vw OR WE CAN SET A DEFAULT THRESHOLD ARRAY AND THEN DO THE TEST HERE...
    IF(N/=C%Opts%NThrsh) &
         & CALL MondoHalt(PRSE_ERROR,' Accuracies in sequence is short.') 
    IF(N/=C%Sets%NBSets) &
         & CALL MondoHalt(PRSE_ERROR,' Basis sets in sequence is short.')
    !
  END SUBROUTINE GlbConflictCheck1
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
  !
  !



#ifdef LJDFLSJDFLSDJFLSKJFDLSDJF
  SUBROUTINE ConflictCheck1(O,D)
    ! Is MondoSCF compiled correctly for parallel replicas?
    IF(O%GradOpt==GRAD_DYNAMICS.AND.D%MDAlgorithm==MD_PARALLEL_REP)THEN
#ifdef !defined(PARALLEL)
       CALL MondoHalt(PRSE_ERROR,' MondoSCF must be compiled in parallel for replica exchange to be active.')
#endif
    ENDIF
  END SUBROUTINE ConflictCheck1
#endif


END MODULE Conflicted
