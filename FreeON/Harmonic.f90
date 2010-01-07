MODULE Harmonic

  USE ControlStructures
  USE MondoLogger
  USE Parse

  IMPLICIT NONE

CONTAINS

  SUBROUTINE HarmonicPotential(O, G)
    TYPE(Options)    :: O
    TYPE(Geometries) :: G
    INTEGER          :: iCLONE, i, j
    REAL(DOUBLE)     :: R, R2

    DO iCLONE = 1, G%Clones
      G%Clone(iCLONE)%ETotal = Zero
      DO i = 1, G%Clone(iCLONE)%NAtms
        DO j = i+1, G%Clone(iCLONE)%NAtms
          R2 = (G%Clone(iCLONE)%Carts%D(1, j)-G%Clone(iCLONE)%Carts%D(1, i))**2 &
             + (G%Clone(iCLONE)%Carts%D(2, j)-G%Clone(iCLONE)%Carts%D(2, i))**2 &
             + (G%Clone(iCLONE)%Carts%D(3, j)-G%Clone(iCLONE)%Carts%D(3, i))**2
          R = SQRT(R2)

          CALL MondoLog(DEBUG_NONE, "Harmonic", "R = "//TRIM(DblToChar(R*AUToAngstroms))//" A", &
            "i = "//TRIM(IntToChar(i))//", j = "//TRIM(IntToChar(j)))

          G%Clone(iCLONE)%ETotal = G%Clone(iCLONE)%ETotal + O%LennardJonesEpsilon*(R-O%LennardJonesR0)**2
        ENDDO
      ENDDO
      !CALL MondoLog(DEBUG_NONE, "Harmonic", "ETotal = " &
      !  //TRIM(FltToChar(G%Clone(iCLONE)%ETotal/O%LennardJonesEpsilon))//" epsilon = " &
      !  //TRIM(DblToChar(G%Clone(iCLONE)%ETotal*au2eV))//" eV", "Clone "//TRIM(IntToChar(iCLONE)))
      CALL MondoLog(DEBUG_NONE, "Harmonic", "ETotal = "//TRIM(DblToChar(G%Clone(iCLONE)%ETotal)), "Clone "//TRIM(IntToChar(iCLONE)))
    ENDDO
  END SUBROUTINE HarmonicPotential

  SUBROUTINE HarmonicForce(O, G)
    TYPE(Options)    :: O
    TYPE(Geometries) :: G
    INTEGER          :: iCLONE, i, j
    REAL(DOUBLE)     :: R, R2, FMagnitude

    DO iCLONE = 1, G%Clones
      G%Clone(iCLONE)%Gradients%D = Zero
      G%Clone(iCLONE)%GradMax = Zero
      G%Clone(iCLONE)%GradRMS = Zero
      DO i = 1, G%Clone(iCLONE)%NAtms
        DO j = 1, G%Clone(iCLONE)%NAtms
          IF(i == j) THEN
            CYCLE
          ENDIF

          R2 = (G%Clone(iCLONE)%Carts%D(1, i)-G%Clone(iCLONE)%Carts%D(1, j))**2 &
             + (G%Clone(iCLONE)%Carts%D(2, i)-G%Clone(iCLONE)%Carts%D(2, j))**2 &
             + (G%Clone(iCLONE)%Carts%D(3, i)-G%Clone(iCLONE)%Carts%D(3, j))**2
          R = SQRT(R2)

          FMagnitude = 2*O%LennardJonesEpsilon*(R-O%LennardJonesR0)/R
          G%Clone(iCLONE)%Gradients%D(1, i) = G%Clone(iCLONE)%Gradients%D(1, i) + FMagnitude*(G%Clone(iCLONE)%Carts%D(1, i)-G%Clone(iCLONE)%Carts%D(1, j))
          G%Clone(iCLONE)%Gradients%D(2, i) = G%Clone(iCLONE)%Gradients%D(2, i) + FMagnitude*(G%Clone(iCLONE)%Carts%D(2, i)-G%Clone(iCLONE)%Carts%D(2, j))
          G%Clone(iCLONE)%Gradients%D(3, i) = G%Clone(iCLONE)%Gradients%D(3, i) + FMagnitude*(G%Clone(iCLONE)%Carts%D(3, i)-G%Clone(iCLONE)%Carts%D(3, j))
        ENDDO

        G%Clone(iCLONE)%GradMax = MAX(G%Clone(iCLONE)%GradMax, ABS(G%Clone(iCLONE)%Gradients%D(1, i)))
        G%Clone(iCLONE)%GradMax = MAX(G%Clone(iCLONE)%GradMax, ABS(G%Clone(iCLONE)%Gradients%D(2, i)))
        G%Clone(iCLONE)%GradMax = MAX(G%Clone(iCLONE)%GradMax, ABS(G%Clone(iCLONE)%Gradients%D(3, i)))

        G%Clone(iCLONE)%GradRMS = G%Clone(iCLONE)%GradRMS &
          + G%Clone(iCLONE)%Gradients%D(1, i)**2 &
          + G%Clone(iCLONE)%Gradients%D(2, i)**2 &
          + G%Clone(iCLONE)%Gradients%D(3, i)**2
      ENDDO
      G%Clone(iCLONE)%GradRMS = SQRT(G%Clone(iCLONE)%GradRMS/DBLE(3*G%Clone(iCLONE)%NAtms))

      CALL MondoLog(DEBUG_NONE, "Harmonic", "Gradient", "Clone "//TRIM(IntToChar(iCLONE)))
      DO i = 1, G%Clone(iCLONE)%NAtms
        CALL MondoLog(DEBUG_NONE, "Harmonic", "Grad["//TRIM(IntToChar(i))//"] = " &
          //TRIM(FltToChar(G%Clone(iCLONE)%Gradients%D(1, i)))//" " &
          //TRIM(FltToChar(G%Clone(iCLONE)%Gradients%D(2, i)))//" " &
          //TRIM(FltToChar(G%Clone(iCLONE)%Gradients%D(3, i))), &
          "Clone "//TRIM(IntToChar(iCLONE)))
      ENDDO
    ENDDO
  END SUBROUTINE HarmonicForce

END MODULE Harmonic
