MODULE Conflicted
  USE Parse
  USE Functionals
  USE DerivedTypes
  USE ProcessControl
  USE ControlStructures
  IMPLICIT NONE
  !
  PRIVATE :: GlbConflictCheck1
  PRIVATE :: GeoConflictCheck1
  PRIVATE :: PBCConflictCheck1
  PRIVATE :: BStConflictCheck1
  !
CONTAINS 
!H---------------------------------------------------------------------------------
!H SUBROUTINE ConflictCheck1(C)
!H  Checking of conflicts for options, PBCs, geometries, etc etc.
!H  >>>>> TO BE DONE ONLY AFTER ALL THE PARSING AND MASSAGING <<<<<<
!H---------------------------------------------------------------------------------
  SUBROUTINE ConflictCheck(C)
    TYPE(Controls) :: C
    CALL GlbConflictCheck(C)
    CALL OptConflictCheck(C%Opts)
    CALL GeoConflictCheck(C%Geos)
    CALL PBCConflictCheck(C%Geos)
    CALL BStConflictCheck(C)
  END SUBROUTINE ConflictCheck

  SUBROUTINE OptConflictCheck(O)
    TYPE(Options) :: O
    IF(O%NModls==0) &
         CALL MondoHalt(PRSE_ERROR,'Option '//MODEL_OPTION//' not set in input.'//RTRN   &
         //'Options include '//        &      
         MODEL_ExactX//', '//          &      
         MODEL_SD//', '//              &       
         MODEL_XA//', '//              &      
         MODEL_B88x//', '//            &      
         MODEL_PBEx//', '//            &      
         MODEL_PW91x//', '//           &      
         MODEL_VWN3//', '//            &      
         MODEL_VWN5//', '//            &      
         MODEL_PW91PW91//', '//        &      
         MODEL_PW91LYP//', '//         &      
         MODEL_BLYP//', '//            &      
         MODEL_PBEPBE//', '//          &      
         MODEL_B3LYP_VWN3//', '//      &      
         MODEL_B3LYP_VWN5//', '//      &      
         MODEL_PBE0//', '//            &      
         MODEL_X3LYP)      
  END SUBROUTINE OptConflictCheck
  !
  SUBROUTINE GlbConflictCheck(C)
    TYPE(Controls) :: C
    CALL GlbConflictCheck1(C)
  END SUBROUTINE GlbConflictCheck
!H---------------------------------------------------------------------------------
!H SUBROUTINE GlbConflictCheck1(C)
!H  Check that the number of options in the progression of 
!H  models/methods/accuracies/basissets are consistent with each other 
!H---------------------------------------------------------------------------------
  SUBROUTINE GlbConflictCheck1(C)
    TYPE(Controls) :: C
    INTEGER        ::   N=0
    N=MAX(N,C%Opts%NModls)
    N=MAX(N,C%Opts%NMthds)
    N=MAX(N,C%Opts%NThrsh)
    N=MAX(N,C%Sets%NBSets)
    IF(N/=C%Opts%NModls)                                                          &
       CALL MondoHalt(PRSE_ERROR,' Model chemistries in sequence is short.'//RTRN &
       //' NModls = '//IntToChar(C%Opts%NModls)//RTRN                             &
       //' NMthds = '//IntToChar(C%Opts%NMthds)//RTRN                             &
       //' NThrsh = '//IntToChar(C%Opts%NThrsh)//RTRN                             &
       //' NBSets = '//IntToChar(C%Sets%NBSets))
    IF(N/=C%Opts%NMthds) &
       CALL MondoHalt(PRSE_ERROR,' SCF methods in sequence is short.')
    !
    !vw PROBLEM WITH THIS TEST, CAUSE WE NEED THE THRESHOLD ARRAY TO SET THE
    !vw BSCR DIMENSIONS IN LoadBasisSets! THE TEST MAY BE DONE IN LoadBasisSets
    !vw OR WE CAN SET A DEFAULT THRESHOLD ARRAY AND THEN DO THE TEST HERE...
    IF(N/=C%Opts%NThrsh) &
         & CALL MondoHalt(PRSE_ERROR,' Accuracies in sequence is short.') 
    IF(N/=C%Sets%NBSets) &
         & CALL MondoHalt(PRSE_ERROR,' Basis sets in sequence is short.')
  END SUBROUTINE GlbConflictCheck1
  !
  SUBROUTINE GeoConflictCheck(G)
    TYPE(Geometries) :: G
!    CALL GeoConflictCheck1(G)
  END SUBROUTINE GeoConflictCheck
!H---------------------------------------------------------------------------------
!H SUBROUTINE GeoConflictCheck1(G)
!H  This routine checks if a pair of atoms are too close from each other.
!H---------------------------------------------------------------------------------
  SUBROUTINE GeoConflictCheck1(G)
    TYPE(Geometries)               :: G
    INTEGER                        :: i,iClone,AtA,AtB,MaxClone
    INTEGER, DIMENSION(2)          :: NClone
    REAL(DOUBLE)                   :: Ax,Ay,Az,Bx,By,Bz,Dist,t1,t2
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: Text
    REAL(DOUBLE), PARAMETER        :: MinDist = 0.5D+00
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

  SUBROUTINE PBCConflictCheck(G)
    TYPE(Geometries) :: G
    CALL PBCConflictCheck1(G)
  END SUBROUTINE PBCConflictCheck

  SUBROUTINE PBCConflictCheck1(G)
    TYPE(Geometries)  :: G
    CHARACTER(LEN=20) :: ErrString
    INTEGER :: I,J,K
    DO I=1,G%Clones
       ErrString=""
       DO K=1,3
          IF(G%Clone(I)%PBC%AutoW%I(K)==0)THEN
             ErrString=TRIM(ErrString)//'F'
          ELSEIF(G%Clone(I)%PBC%AutoW%I(K)==1)THEN
             ErrString=TRIM(ErrString)//'T'
          ELSE
             ErrString=TRIM(ErrString)//'?'
          ENDIF
       ENDDO
       DO J=1,3
          IF(G%Clone(I)%PBC%AutoW%I(J)==BIG_INT)THEN
             CALL MondoHalt(PRSE_ERROR,'Periodic invocation hosed: PBC = '//TRIM(ErrString))
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE PBCConflictCheck1
  !
  !
  SUBROUTINE BStConflictCheck(C)
    TYPE(Controls) :: C
    CALL BStConflictCheck1(C%Sets,C%Geos)
  END SUBROUTINE BStConflictCheck
  !
  !
  SUBROUTINE BStConflictCheck1(B,G)
    TYPE(BasisSets)  :: B
    TYPE(Geometries) :: G
    !
    IF(B%NBSets.EQ.0) CALL MondoHalt(PRSE_ERROR,'The Number of Basis Set is ZERO!')
    !
  END SUBROUTINE BStConflictCheck1
  !
  !

END MODULE Conflicted
