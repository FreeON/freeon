      LOGICAL FUNCTION SISNAN(SIN)
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      REAL SIN
*     ..
*
*  Purpose
*  =======
*
*  SISNAN returns .TRUE. if its argument is NaN, and .FALSE.
*  otherwise.  To be replaced by the Fortran 2003 intrinsic in the
*  future.
*
*  Arguments
*  =========
*
*  SIN      (input) REAL
*          Input to test for NaN.
*
*  =====================================================================
*
*  .. External Functions ..
      LOGICAL SLAISNAN
      EXTERNAL SLAISNAN
*  ..
*  .. Executable Statements ..
      SISNAN = SLAISNAN(SIN,SIN)
      RETURN
      END
