      FUNCTION iT(I,J)
      IMPLICIT NONE
      INTEGER I,J,iT,iTD(6,6)
      INTEGER T2Size(28)
      SAVE iTD
      DATA iTD/ 1, 0, 0, 0, 0, 0,
     $          2, 3, 0, 0, 0, 0,
     $          4, 5, 6, 0, 0, 0,
     $          0, 0, 0, 0, 0, 0,
     $          0, 0, 0, 0, 0, 0,
     $          7, 8, 9, 0, 0, 10 /
      iT=iTD(I,J)
      RETURN
      END

