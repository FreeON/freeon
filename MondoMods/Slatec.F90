MODULE SLATEC
USE GlobalScalars
USE ProcessControl
IMPLICIT REAL(DOUBLE) (A-H,O-Z)
IMPLICIT INTEGER (I-N)
CONTAINS
!
!--------------------------------------------------------------------
!
!DECK FDUMP
      SUBROUTINE FDUMP
!***BEGIN PROLOGUE  FDUMP
!***PURPOSE  Symbolic dump (should be locally written).
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3
!***TYPE      ALL (FDUMP-A)
!***KEYWORDS  ERROR, XERMSG
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!        ***Note*** Machine Dependent Routine
!        FDUMP is intended to be replaced by a locally written
!        version which produces a symbolic dump.  Failing this,
!        it should be replaced by a version which prints the
!        subprogram nesting list.  Note that this dump must be
!        printed on each of up to five files, as indicated by the
!        XGETUA routine.  See XSETUA and XGETUA for details.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  FDUMP
!***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END SUBROUTINE
!*DECK I1MACH
       INTEGER FUNCTION I1MACH (I)
!C***BEGIN PROLOGUE  I1MACH
!C***PURPOSE  Return integer machine dependent constants.
!C***LIBRARY   SLATEC
!C***CATEGORY  R1
!C***TYPE      INTEGER (I1MACH-I)
!C***KEYWORDS  MACHINE CONSTANTS
!C***AUTHOR  Fox, P. A., (Bell Labs)
!C           Hall, A. D., (Bell Labs)
!C           Schryer, N. L., (Bell Labs)
!C***DESCRIPTION
!C
!C   I1MACH can be used to obtain machine-dependent parameters for the
!C   local machine environment.  It is a function subprogram with one
!C   (input) argument and can be referenced as follows:
!C
!C        K = I1MACH(I)
!C
!C   where I=1,...,16.  The (output) value of K above is determined by
!C   the (input) value of I.  The results for various values of I are
!C   discussed below.
!C
!C   I/O unit numbers:
!C     I1MACH( 1) = the standard input unit.
!C     I1MACH( 2) = the standard output unit.
!C     I1MACH( 3) = the standard punch unit.
!C     I1MACH( 4) = the standard error message unit.
!C
!C   Words:
!C     I1MACH( 5) = the number of bits per integer storage unit.
!C     I1MACH( 6) = the number of characters per integer storage unit.
!C
!C   Integers:
!C     assume integers are represented in the S-digit, base-A form
!C
!C                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!C
!C                where 0 .LE. X(I) .LT. A for I=0,...,S-1.
!C     I1MACH( 7) = A, the base.
!C     I1MACH( 8) = S, the number of base-A digits.
!C     I1MACH( 9) = A**S - 1, the largest magnitude.
!C
!C   Floating-Point Numbers:
!C     Assume floating-point numbers are represented in the T-digit,
!C     base-B form
!C                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!C
!C                where 0 .LE. X(I) .LT. B for I=1,...,T,
!C                0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
!C     I1MACH(10) = B, the base.
!C
!C   Single-Precision:
!C     I1MACH(11) = T, the number of base-B digits.
!C     I1MACH(12) = EMIN, the smallest exponent E.
!C     I1MACH(13) = EMAX, the largest exponent E.
!C
!C   Double-Precision:
!C     I1MACH(14) = T, the number of base-B digits.
!C     I1MACH(15) = EMIN, the smallest exponent E.
!C     I1MACH(16) = EMAX, the largest exponent E.
!C
!C   To alter this function for a particular environment, the desired
!C   set of DATA statements should be activated by removing the C from
!C   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
!C   checked for consistency with the local operating system.
!C
!C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!C                 a portable library, ACM Transactions on Mathematical
!C                 Software 4, 2 (June 1978), pp. 177-188.
!C***ROUTINES CALLED  (NONE)
!C***REVISION HISTORY  (YYMMDD)
!C   750101  DATE WRITTEN
!C   891012  Added VAX G-floating constants.  (WRB)
!C   891012  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900618  Added DEC RISC constants.  (WRB)
!C   900723  Added IBM RS 6000 constants.  (WRB)
!C   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16.
!C           (RWC)
!C   910710  Added HP 730 constants.  (SMR)
!C   911114  Added Convex IEEE constants.  (WRB)
!C   920121  Added SUN -r8 compiler option constants.  (WRB)
!C   920229  Added Touchstone Delta i860 constants.  (WRB)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C   920625  Added Convex -p8 and -pd8 compiler option constants.
!C           (BKS, WRB)
!C   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
!C   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler
!C           options.  (DWL, RWC and WRB).
!C***END PROLOGUE  I1MACH
!C
       INTEGER IMACH(16),OUTPUT
       SAVE IMACH
       EQUIVALENCE (IMACH(4),OUTPUT)
!C
!C     MACHINE CONSTANTS FOR THE AMIGA
!C     ABSOFT COMPILER
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          5 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -126 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1022 /
!C     DATA IMACH(16) /       1023 /
!C
!C     MACHINE CONSTANTS FOR THE APOLLO
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          6 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -125 /
!C     DATA IMACH(13) /        129 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1021 /
!C     DATA IMACH(16) /       1025 /
!C
!C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
!C
!C     DATA IMACH( 1) /          7 /
!C     DATA IMACH( 2) /          2 /
!C     DATA IMACH( 3) /          2 /
!C     DATA IMACH( 4) /          2 /
!C     DATA IMACH( 5) /         36 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         33 /
!C     DATA IMACH( 9) / Z1FFFFFFFF /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -256 /
!C     DATA IMACH(13) /        255 /
!C     DATA IMACH(14) /         60 /
!C     DATA IMACH(15) /       -256 /
!C     DATA IMACH(16) /        255 /
!C
!C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          7 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         48 /
!C     DATA IMACH( 6) /          6 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         39 /
!C     DATA IMACH( 9) / O0007777777777777 /
!C     DATA IMACH(10) /          8 /
!C     DATA IMACH(11) /         13 /
!C     DATA IMACH(12) /        -50 /
!C     DATA IMACH(13) /         76 /
!C     DATA IMACH(14) /         26 /
!C     DATA IMACH(15) /        -50 /
!C     DATA IMACH(16) /         76 /
!C
!C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          7 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         48 /
!C     DATA IMACH( 6) /          6 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         39 /
!C     DATA IMACH( 9) / O0007777777777777 /
!C     DATA IMACH(10) /          8 /
!C     DATA IMACH(11) /         13 /
!C     DATA IMACH(12) /        -50 /
!C     DATA IMACH(13) /         76 /
!C     DATA IMACH(14) /         26 /
!C     DATA IMACH(15) /     -32754 /
!C     DATA IMACH(16) /      32780 /
!C
!C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          7 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         64 /
!C     DATA IMACH( 6) /          8 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         63 /
!C     DATA IMACH( 9) / 9223372036854775807 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         47 /
!C     DATA IMACH(12) /      -4095 /
!C     DATA IMACH(13) /       4094 /
!C     DATA IMACH(14) /         94 /
!C     DATA IMACH(15) /      -4095 /
!C     DATA IMACH(16) /       4094 /
!C
!C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          7 /
!C     DATA IMACH( 4) /    6LOUTPUT/
!C     DATA IMACH( 5) /         60 /
!C     DATA IMACH( 6) /         10 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         48 /
!C     DATA IMACH( 9) / 00007777777777777777B /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         47 /
!C     DATA IMACH(12) /       -929 /
!C     DATA IMACH(13) /       1070 /
!C     DATA IMACH(14) /         94 /
!C     DATA IMACH(15) /       -929 /
!C     DATA IMACH(16) /       1069 /
!C
!C     MACHINE CONSTANTS FOR THE CELERITY C1260
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          6 /
!C     DATA IMACH( 4) /          0 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / Z'7FFFFFFF' /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -126 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1022 /
!C     DATA IMACH(16) /       1023 /
!C
!C     MACHINE CONSTANTS FOR THE CONVEX
!C     USING THE -fn COMPILER OPTION
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          7 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -127 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1023 /
!C     DATA IMACH(16) /       1023 /
!C
!C     MACHINE CONSTANTS FOR THE CONVEX
!C     USING THE -fi COMPILER OPTION
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          7 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -125 /
!C     DATA IMACH(13) /        128 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1021 /
!C     DATA IMACH(16) /       1024 /
!C
!C     MACHINE CONSTANTS FOR THE CONVEX
!C     USING THE -p8 COMPILER OPTION
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          7 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         64 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         63 /
!C     DATA IMACH( 9) / 9223372036854775807 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         53 /
!C     DATA IMACH(12) /      -1023 /
!C     DATA IMACH(13) /       1023 /
!C     DATA IMACH(14) /        113 /
!C     DATA IMACH(15) /     -16383 /
!C     DATA IMACH(16) /      16383 /
!C
!C     MACHINE CONSTANTS FOR THE CONVEX
!C     USING THE -pd8 COMPILER OPTION
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          7 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         64 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         63 /
!C     DATA IMACH( 9) / 9223372036854775807 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         53 /
!C     DATA IMACH(12) /      -1023 /
!C     DATA IMACH(13) /       1023 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1023 /
!C     DATA IMACH(16) /       1023 /
!C
!C     MACHINE CONSTANTS FOR THE CRAY
!C     USING THE 46 BIT INTEGER COMPILER OPTION
!C
!C     DATA IMACH( 1) /        100 /
!C     DATA IMACH( 2) /        101 /
!C     DATA IMACH( 3) /        102 /
!C     DATA IMACH( 4) /        101 /
!C     DATA IMACH( 5) /         64 /
!C     DATA IMACH( 6) /          8 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         46 /
!C     DATA IMACH( 9) / 1777777777777777B /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         47 /
!C     DATA IMACH(12) /      -8189 /
!C     DATA IMACH(13) /       8190 /
!C     DATA IMACH(14) /         94 /
!C     DATA IMACH(15) /      -8099 /
!C     DATA IMACH(16) /       8190 /
!C
!C     MACHINE CONSTANTS FOR THE CRAY
!C     USING THE 64 BIT INTEGER COMPILER OPTION
!C
!C     DATA IMACH( 1) /        100 /
!C     DATA IMACH( 2) /        101 /
!C     DATA IMACH( 3) /        102 /
!C     DATA IMACH( 4) /        101 /
!C     DATA IMACH( 5) /         64 /
!C     DATA IMACH( 6) /          8 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         63 /
!C     DATA IMACH( 9) / 777777777777777777777B /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         47 /
!C     DATA IMACH(12) /      -8189 /
!C     DATA IMACH(13) /       8190 /
!C     DATA IMACH(14) /         94 /
!C     DATA IMACH(15) /      -8099 /
!C     DATA IMACH(16) /       8190 /
!C
!C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!C
!C     DATA IMACH( 1) /         11 /
!C     DATA IMACH( 2) /         12 /
!C     DATA IMACH( 3) /          8 /
!C     DATA IMACH( 4) /         10 /
!C     DATA IMACH( 5) /         16 /
!C     DATA IMACH( 6) /          2 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         15 /
!C     DATA IMACH( 9) /      32767 /
!C     DATA IMACH(10) /         16 /
!C     DATA IMACH(11) /          6 /
!C     DATA IMACH(12) /        -64 /
!C     DATA IMACH(13) /         63 /
!C     DATA IMACH(14) /         14 /
!C     DATA IMACH(15) /        -64 /
!C     DATA IMACH(16) /         63 /
!C
!C     MACHINE CONSTANTS FOR THE DEC ALPHA
!C     USING G_FLOAT
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          5 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -127 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1023 /
!C     DATA IMACH(16) /       1023 /
!C
!C     MACHINE CONSTANTS FOR THE DEC ALPHA
!C     USING IEEE_FLOAT
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          6 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -125 /
!C     DATA IMACH(13) /        128 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1021 /
!C     DATA IMACH(16) /       1024 /
!C
!C     MACHINE CONSTANTS FOR THE DEC RISC
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          6 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -125 /
!C     DATA IMACH(13) /        128 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1021 /
!C     DATA IMACH(16) /       1024 /
!C
!C     MACHINE CONSTANTS FOR THE DEC VAX
!C     USING D_FLOATING
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          5 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -127 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         56 /
!C     DATA IMACH(15) /       -127 /
!C     DATA IMACH(16) /        127 /
!C
!C     MACHINE CONSTANTS FOR THE DEC VAX
!C     USING G_FLOATING
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          5 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -127 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1023 /
!C     DATA IMACH(16) /       1023 /
!C
!C     MACHINE CONSTANTS FOR THE ELXSI 6400
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          6 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         32 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -126 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1022 /
!C     DATA IMACH(16) /       1023 /
!C
!C     MACHINE CONSTANTS FOR THE HARRIS 220
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          0 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         24 /
!C     DATA IMACH( 6) /          3 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         23 /
!C     DATA IMACH( 9) /    8388607 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         23 /
!C     DATA IMACH(12) /       -127 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         38 /
!C     DATA IMACH(15) /       -127 /
!C     DATA IMACH(16) /        127 /
!C
!C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /         43 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         36 /
!C     DATA IMACH( 6) /          6 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         35 /
!C     DATA IMACH( 9) / O377777777777 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         27 /
!C     DATA IMACH(12) /       -127 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         63 /
!C     DATA IMACH(15) /       -127 /
!C     DATA IMACH(16) /        127 /
!C
!C     MACHINE CONSTANTS FOR THE HP 730
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          6 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -125 /
!C     DATA IMACH(13) /        128 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1021 /
!C     DATA IMACH(16) /       1024 /
!C
!C     MACHINE CONSTANTS FOR THE HP 2100
!C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          4 /
!C     DATA IMACH( 4) /          1 /
!C     DATA IMACH( 5) /         16 /
!C     DATA IMACH( 6) /          2 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         15 /
!C     DATA IMACH( 9) /      32767 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         23 /
!C     DATA IMACH(12) /       -128 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         39 /
!C     DATA IMACH(15) /       -128 /
!C     DATA IMACH(16) /        127 /
!C
!C     MACHINE CONSTANTS FOR THE HP 2100
!C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          4 /
!C     DATA IMACH( 4) /          1 /
!C     DATA IMACH( 5) /         16 /
!C     DATA IMACH( 6) /          2 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         15 /
!C     DATA IMACH( 9) /      32767 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         23 /
!C     DATA IMACH(12) /       -128 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         55 /
!C     DATA IMACH(15) /       -128 /
!C     DATA IMACH(16) /        127 /
!C
!C     MACHINE CONSTANTS FOR THE HP 9000
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          6 /
!C     DATA IMACH( 4) /          7 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         32 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -126 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1015 /
!C     DATA IMACH(16) /       1017 /
!C
!C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
!C     THE PERKIN ELMER (INTERDATA) 7/32.
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          7 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) /  Z7FFFFFFF /
!C     DATA IMACH(10) /         16 /
!C     DATA IMACH(11) /          6 /
!C     DATA IMACH(12) /        -64 /
!C     DATA IMACH(13) /         63 /
!C     DATA IMACH(14) /         14 /
!C     DATA IMACH(15) /        -64 /
!C     DATA IMACH(16) /         63 /
!C
!C     MACHINE CONSTANTS FOR THE IBM PC
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          0 /
!C     DATA IMACH( 4) /          0 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -125 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1021 /
!C     DATA IMACH(16) /       1023 /
!C
!C     MACHINE CONSTANTS FOR THE IBM RS 6000
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          6 /
!C     DATA IMACH( 4) /          0 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -125 /
!C     DATA IMACH(13) /        128 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1021 /
!C     DATA IMACH(16) /       1024 /
!C
!C     MACHINE CONSTANTS FOR THE INTEL i860
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          6 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -125 /
!C     DATA IMACH(13) /        128 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1021 /
!C     DATA IMACH(16) /       1024 /
!C
!C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          5 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         36 /
!C     DATA IMACH( 6) /          5 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         35 /
!C     DATA IMACH( 9) / "377777777777 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         27 /
!C     DATA IMACH(12) /       -128 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         54 /
!C     DATA IMACH(15) /       -101 /
!C     DATA IMACH(16) /        127 /
!C
!C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          5 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         36 /
!C     DATA IMACH( 6) /          5 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         35 /
!C     DATA IMACH( 9) / "377777777777 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         27 /
!C     DATA IMACH(12) /       -128 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         62 /
!C     DATA IMACH(15) /       -128 /
!C     DATA IMACH(16) /        127 /
!C
!C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!C     32-BIT INTEGER ARITHMETIC.
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          5 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -127 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         56 /
!C     DATA IMACH(15) /       -127 /
!C     DATA IMACH(16) /        127 /
!C
!C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!C     16-BIT INTEGER ARITHMETIC.
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          5 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         16 /
!C     DATA IMACH( 6) /          2 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         15 /
!C     DATA IMACH( 9) /      32767 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -127 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         56 /
!C     DATA IMACH(15) /       -127 /
!C     DATA IMACH(16) /        127 /
!C
!C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          6 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -125 /
!C     DATA IMACH(13) /        128 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1021 /
!C     DATA IMACH(16) /       1024 /
!C
!C     MACHINE CONSTANTS FOR THE SUN
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          6 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -125 /
!C     DATA IMACH(13) /        128 /
!C     DATA IMACH(14) /         53 /
!C     DATA IMACH(15) /      -1021 /
!C     DATA IMACH(16) /       1024 /
!C
!C     MACHINE CONSTANTS FOR THE SUN
!C     USING THE -r8 COMPILER OPTION
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          6 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         32 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         31 /
!C     DATA IMACH( 9) / 2147483647 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         53 /
!C     DATA IMACH(12) /      -1021 /
!C     DATA IMACH(13) /       1024 /
!C     DATA IMACH(14) /        113 /
!C     DATA IMACH(15) /     -16381 /
!C     DATA IMACH(16) /      16384 /
!C
!C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
!C
!C     DATA IMACH( 1) /          5 /
!C     DATA IMACH( 2) /          6 /
!C     DATA IMACH( 3) /          1 /
!C     DATA IMACH( 4) /          6 /
!C     DATA IMACH( 5) /         36 /
!C     DATA IMACH( 6) /          4 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         35 /
!C     DATA IMACH( 9) / O377777777777 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         27 /
!C     DATA IMACH(12) /       -128 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         60 /
!C     DATA IMACH(15) /      -1024 /
!C     DATA IMACH(16) /       1023 /
!C
!C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
!C
!C     DATA IMACH( 1) /          1 /
!C     DATA IMACH( 2) /          1 /
!C     DATA IMACH( 3) /          0 /
!C     DATA IMACH( 4) /          1 /
!C     DATA IMACH( 5) /         16 /
!C     DATA IMACH( 6) /          2 /
!C     DATA IMACH( 7) /          2 /
!C     DATA IMACH( 8) /         15 /
!C     DATA IMACH( 9) /      32767 /
!C     DATA IMACH(10) /          2 /
!C     DATA IMACH(11) /         24 /
!C     DATA IMACH(12) /       -127 /
!C     DATA IMACH(13) /        127 /
!C     DATA IMACH(14) /         56 /
!C     DATA IMACH(15) /       -127 /
!C     DATA IMACH(16) /        127 /
!C
!C***FIRST EXECUTABLE STATEMENT  I1MACH
       IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10
!C
       I1MACH = IMACH(I)
       RETURN
!C
    10 CONTINUE
     ! WRITE (UNIT = OUTPUT, FMT = 9000)
  9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
!C
!C     CALL FDUMP
!C
       STOP
       END FUNCTION   
!*DECK J4SAVE
       FUNCTION J4SAVE (IWHICH, IVALUE, ISET)
!C***BEGIN PROLOGUE  J4SAVE
!C***SUBSIDIARY
!C***PURPOSE  Save or recall global variables needed by error
!C            handling routines.
!C***LIBRARY   SLATEC (XERROR)
!C***TYPE      INTEGER (J4SAVE-I)
!C***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
!C***AUTHOR  Jones, R. E., (SNLA)
!C***DESCRIPTION
!C
!C     Abstract
!C        J4SAVE saves and recalls several global variables needed
!C        by the library error handling routines.
!C
!C     Description of Parameters
!C      --Input--
!C        IWHICH - Index of item desired.
!C                = 1 Refers to current error number.
!C                = 2 Refers to current error control flag.
!C                = 3 Refers to current unit number to which error
!C                    messages are to be sent.  (0 means use standard.)
!C                = 4 Refers to the maximum number of times any
!C                     message is to be printed (as set by XERMAX).
!C                = 5 Refers to the total number of units to which
!C                     each error message is to be written.
!C                = 6 Refers to the 2nd unit for error messages
!C                = 7 Refers to the 3rd unit for error messages
!C                = 8 Refers to the 4th unit for error messages
!C                = 9 Refers to the 5th unit for error messages
!C        IVALUE - The value to be set for the IWHICH-th parameter,
!C                 if ISET is .TRUE. .
!C        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
!C                 given the value, IVALUE.  If ISET=.FALSE., the
!C                 IWHICH-th parameter will be unchanged, and IVALUE
!C                 is a dummy parameter.
!C      --Output--
!C        The (old) value of the IWHICH-th parameter will be returned
!C        in the function value, J4SAVE.
!C
!C***SEE ALSO  XERMSG
!C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!C                 Error-handling Package, SAND82-0800, Sandia
!C                 Laboratories, 1982.
!C***ROUTINES CALLED  (NONE)
!C***REVISION HISTORY  (YYMMDD)
!C   790801  DATE WRITTEN
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900205  Minor modifications to prologue.  (WRB)
!C   900402  Added TYPE section.  (WRB)
!C   910411  Added KEYWORDS section.  (WRB)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  J4SAVE
       LOGICAL ISET
       INTEGER IPARAM(9)
       SAVE IPARAM
       DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
       DATA IPARAM(5)/1/
       DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
!C***FIRST EXECUTABLE STATEMENT  J4SAVE
       J4SAVE = IPARAM(IWHICH)
       IF (ISET) IPARAM(IWHICH) = IVALUE
       RETURN
       END FUNCTION   
!*DECK POLFIT
       SUBROUTINE POLFIT (N, X, Y, W, MAXDEG, NDEG, EPS, R, IERR, A)
!C***BEGIN PROLOGUE  POLFIT
!C***PURPOSE  Fit discrete data in a least squares sense by polynomials
!C            in one variable.
!C***LIBRARY   SLATEC
!C***CATEGORY  K1A1A2
!C***TYPE      SINGLE PRECISION (POLFIT-S, DPOLFT-D)
!C***KEYWORDS  CURVE FITTING, DATA FITTING, LEAST SQUARES, POLYNOMIAL FIT
!C***AUTHOR  Shampine, L. F., (SNLA)
!C           Davenport, S. M., (SNLA)
!C           Huddleston, R. E., (SNLL)
!C***DESCRIPTION
!C
!C     Abstract
!C
!C     Given a collection of points X(I) and a set of values Y(I) which
!C     correspond to some function or measurement at each of the X(I),
!C     subroutine  POLFIT  computes the weighted least-squares polynomial
!C     fits of all degrees up to some degree either specified by the user
!C     or determined by the routine.  The fits thus obtained are in
!C     orthogonal polynomial form.  Subroutine  PVALUE  may then be
!C     called to evaluate the fitted polynomials and any of their
!C     derivatives at any point.  The subroutine  PCOEF  may be used to
!C     express the polynomial fits as powers of (X-C) for any specified
!C     point C.
!C
!C     The parameters for  POLFIT  are
!C
!C     Input --
!C         N -      the number of data points.  The arrays X, Y and W
!C                  must be dimensioned at least  N  (N .GE. 1).
!C         X -      array of values of the independent variable.  These
!C                  values may appear in any order and need not all be
!C                  distinct.
!C         Y -      array of corresponding function values.
!C         W -      array of positive values to be used as weights.  If
!C                  W(1) is negative,  POLFIT  will set all the weights
!C                  to One, which means unweighted least squares error
!C                  will be minimized.  To minimize relative error, the
!C                  user should set the weights to:  W(I) = One/Y(I)**2,
!C                  I = 1,...,N .
!C         MAXDEG - maximum degree to be allowed for polynomial fit.
!C                  MAXDEG  may be any non-negative integer less than  N.
!C                  Note -- MAXDEG  cannot be equal to  N-1  when a
!C                  statistical test is to be used for degree selection,
!C                  i.e., when input value of  EPS  is negative.
!C         EPS -    specifies the criterion to be used in determining
!C                  the degree of fit to be computed.
!C                  (1)  If  EPS  is input negative,  POLFIT  chooses the
!C                       degree based on a statistical F test of
!C                       significance.  One of three possible
!C                       significance levels will be used:  .01, .05 or
!C                       .10.  If  EPS=-One , the routine will
!C                       automatically select one of these levels based
!C                       on the number of data points and the maximum
!C                       degree to be considered.  If  EPS  is input as
!C                       -.01, -.05, or -.10, a significance level of
!C                       .01, .05, or .10, respectively, will be used.
!C                  (2)  If  EPS  is set to 0.,  POLFIT  computes the
!C                       polynomials of degrees 0 through  MAXDEG .
!C                  (3)  If  EPS  is input positive,  EPS  is the RMS
!C                       error tolerance which must be satisfied by the
!C                       fitted polynomial.  POLFIT  will increase the
!C                       degree of fit until this criterion is met or
!C                       until the maximum degree is reached.
!C
!C     Output --
!C         NDEG -   degree of the highest degree fit computed.
!C         EPS -    RMS error of the polynomial of degree  NDEG .
!C         R -      vector of dimension at least NDEG containing values
!C                  of the fit of degree  NDEG  at each of the  X(I) .
!C                  Except when the statistical test is used, these
!C                  values are more accurate than results from subroutine
!C                  PVALUE  normally are.
!C         IERR -   error flag with the following possible values.
!C             1 -- indicates normal execution, i.e., either
!C                  (1)  the input value of  EPS  was negative, and the
!C                       computed polynomial fit of degree  NDEG
!C                       satisfies the specified F test, or
!C                  (2)  the input value of  EPS  was 0., and the fits of
!C                       all degrees up to  MAXDEG  are complete, or
!C                  (3)  the input value of  EPS  was positive, and the
!C                       polynomial of degree  NDEG  satisfies the RMS
!C                       error requirement.
!C             2 -- invalid input parameter.  At least one of the input
!C                  parameters has an illegal value and must be corrected
!C                  before  POLFIT  can proceed.  Valid input results
!C                  when the following restrictions are observed
!C                       N .GE. 1
!C                       0 .LE. MAXDEG .LE. N-1  for  EPS .GE. 0.
!C                       0 .LE. MAXDEG .LE. N-2  for  EPS .LT. 0.
!C                       W(1)=-One  or  W(I) .GT. 0., I=1,...,N .
!C             3 -- cannot satisfy the RMS error requirement with a
!C                  polynomial of degree no greater than  MAXDEG .  Best
!C                  fit found is of degree  MAXDEG .
!C             4 -- cannot satisfy the test for significance using
!C                  current value of  MAXDEG .  Statistically, the
!C                  best fit found is of order  NORD .  (In this case,
!C                  NDEG will have one of the values:  MAXDEG-2,
!C                  MAXDEG-1, or MAXDEG).  Using a higher value of
!C                  MAXDEG  may result in passing the test.
!C         A -      work and output array having at least 3N+3MAXDEG+3
!C                  locations
!C
!C     Note - POLFIT  calculates all fits of degrees up to and including
!C            NDEG .  Any or all of these fits can be evaluated or
!C            expressed as powers of (X-C) using  PVALUE  and  PCOEF
!C            after just one call to  POLFIT .
!C
!C***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
!C                 Curve fitting by polynomials in one variable, Report
!C                 SLA-74-0270, Sandia Laboratories, June 1974.
!C***ROUTINES CALLED  PVALUE, XERMSG
!C***REVISION HISTORY  (YYMMDD)
!C   740601  DATE WRITTEN
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890531  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C   920527  Corrected erroneous statements in DESCRIPTION.  (WRB)
!C***END PROLOGUE  POLFIT
       DOUBLE PRECISION TEMD1,TEMD2
       DIMENSION X(*), Y(*), W(*), R(*), A(*)
       DIMENSION YP(1)
       DIMENSION CO(4,3)
       SAVE CO
       DATA  CO(1,1), CO(2,1), CO(3,1), CO(4,1), CO(1,2), CO(2,2), &
             CO(3,2), CO(4,2), CO(1,3), CO(2,3), CO(3,3), &
         CO(4,3)/-13.086850D0,-2.4648165D0,-3.3846535D0,-1.2973162D0, &
                 -3.3381146D0,-1.7812271D0,-3.2578406D0,-1.6589279D0, &
                 -1.6282703D0,-1.3152745D0,-3.2640179D0,-1.9829776D0/
!C***FIRST EXECUTABLE STATEMENT  POLFIT
       M = ABS(N)
       IF (M .EQ. 0) GO TO 30
       IF (MAXDEG .LT. 0) GO TO 30
       A(1) = MAXDEG
       MOP1 = MAXDEG + 1
       IF (M .LT. MOP1) GO TO 30
       IF (EPS .LT. Zero  .AND.  M .EQ. MOP1) GO TO 30
       XM = M
       ETST = EPS*EPS*XM
       IF (W(1) .LT. Zero) GO TO 2
       DO 1 I = 1,M
         IF (W(I) .LE. Zero) GO TO 30
  1      CONTINUE
       GO TO 4
  2    DO 3 I = 1,M
  3      W(I) = One
  4    IF (EPS .GE. Zero) GO TO 8
!C
!C DETERMINE SIGNIFICANCE LEVEL INDEX TO BE USED IN STATISTICAL TEST FOR
!C CHOOSING DEGREE OF POLYNOMIAL FIT
!C
       IF (EPS .GT. (-.55D0)) GO TO 5
       IDEGF = M - MAXDEG - 1
       KSIG = 1
       IF (IDEGF .LT. 10) KSIG = 2
       IF (IDEGF .LT. 5) KSIG = 3
       GO TO 8
  5    KSIG = 1
       IF (EPS .LT. (-.03D0)) KSIG = 2
       IF (EPS .LT. (-.07D0)) KSIG = 3
!C
!C INITIALIZE INDEXES AND COEFFICIENTS FOR FITTING
!C
  8    K1 = MAXDEG + 1
       K2 = K1 + MAXDEG
       K3 = K2 + MAXDEG + 2
       K4 = K3 + M
       K5 = K4 + M
       DO 9 I = 2,K4
  9      A(I) = Zero
       W11 = Zero
       IF (N .LT. 0) GO TO 11
!C
!C UNCONSTRAINED CASE
!C
       DO 10 I = 1,M
         K4PI = K4 + I
         A(K4PI) = One
  10     W11 = W11 + W(I)
       GO TO 13
!C
!C CONSTRAINED CASE
!C
  11   DO 12 I = 1,M
         K4PI = K4 + I
  12     W11 = W11 + W(I)*A(K4PI)**2
!C
!C COMPUTE FIT OF DEGREE ZERO
!C
  13   TEMD1 = Zero
       DO 14 I = 1,M
         K4PI = K4 + I
         TEMD1 = TEMD1 + DBLE(W(I))*DBLE(Y(I))*DBLE(A(K4PI))
  14     CONTINUE
       TEMD1 = TEMD1/DBLE(W11)
       A(K2+1) = TEMD1
       SIGJ = Zero
       DO 15 I = 1,M
         K4PI = K4 + I
         K5PI = K5 + I
         TEMD2 = TEMD1*DBLE(A(K4PI))
         R(I) = TEMD2
         A(K5PI) = TEMD2 - DBLE(R(I))
  15     SIGJ = SIGJ + W(I)*((Y(I)-R(I)) - A(K5PI))**2
       J = 0
!C
!C SEE IF POLYNOMIAL OF DEGREE 0 SATISFIES THE DEGREE SELECTION CRITERION
!C
       IF (EPS) 24,26,27
!C
!C INCREMENT DEGREE
!C
  16   J = J + 1
       JP1 = J + 1
       K1PJ = K1 + J
       K2PJ = K2 + J
       SIGJM1 = SIGJ
!C
!C COMPUTE NEW B COEFFICIENT EXCEPT WHEN J = 1
!C
       IF (J .GT. 1) A(K1PJ) = W11/W1
!C
!C COMPUTE NEW A COEFFICIENT
!C
       TEMD1 = Zero
       DO 18 I = 1,M
         K4PI = K4 + I
         TEMD2 = A(K4PI)
         TEMD1 = TEMD1 + DBLE(X(I))*DBLE(W(I))*TEMD2*TEMD2
  18     CONTINUE
       A(JP1) = TEMD1/DBLE(W11)
!C
!C EVALUATE ORTHOGONAL POLYNOMIAL AT DATA POINTS
!C
       W1 = W11
       W11 = Zero
       DO 19 I = 1,M
         K3PI = K3 + I
         K4PI = K4 + I
         TEMP = A(K3PI)
         A(K3PI) = A(K4PI)
         A(K4PI) = (X(I)-A(JP1))*A(K3PI) - A(K1PJ)*TEMP
  19     W11 = W11 + W(I)*A(K4PI)**2
!C
!C GET NEW ORTHOGONAL POLYNOMIAL COEFFICIENT USING PARTIAL DOUBLE
!C PRECISION
!C
       TEMD1 = Zero
       DO 20 I = 1,M
         K4PI = K4 + I
         K5PI = K5 + I
         TEMD2 = DBLE(W(I))*DBLE((Y(I)-R(I))-A(K5PI))*DBLE(A(K4PI))
  20     TEMD1 = TEMD1 + TEMD2
       TEMD1 = TEMD1/DBLE(W11)
       A(K2PJ+1) = TEMD1
!C
!C UPDATE POLYNOMIAL EVALUATIONS AT EACH OF THE DATA POINTS, AND
!C ACCUMULATE SUM OF SQUARES OF ERRORS.  THE POLYNOMIAL EVALUATIONS ARE
!C COMPUTED AND STORED IN EXTENDED PRECISION.  FOR THE I-TH DATA POINT,
!C THE MOST SIGNIFICANT BITS ARE STORED IN  R(I) , AND THE LEAST
!C SIGNIFICANT BITS ARE IN  A(K5PI) .
!C
       SIGJ = Zero
       DO 21 I = 1,M
         K4PI = K4 + I
         K5PI = K5 + I
         TEMD2 = DBLE(R(I)) + DBLE(A(K5PI)) + TEMD1*DBLE(A(K4PI))
         R(I) = TEMD2
         A(K5PI) = TEMD2 - DBLE(R(I))
  21     SIGJ = SIGJ + W(I)*((Y(I)-R(I)) - A(K5PI))**2
!C
!C SEE IF DEGREE SELECTION CRITERION HAS BEEN SATISFIED OR IF DEGREE
!C MAXDEG  HAS BEEN REACHED
!C
       IF (EPS) 23,26,27
!C
!C COMPUTE F STATISTICS  (INPUT EPS .LT. Zero)
!C
    23 IF (SIGJ .EQ. Zero) GO TO 29
       DEGF = M - J - 1
       DEN = (CO(4,KSIG)*DEGF + One)*DEGF
       FCRIT = (((CO(3,KSIG)*DEGF) + CO(2,KSIG))*DEGF + CO(1,KSIG))/DEN
       FCRIT = FCRIT*FCRIT
       F = (SIGJM1 - SIGJ)*DEGF/SIGJ
       IF (F .LT. FCRIT) GO TO 25
!C
!C POLYNOMIAL OF DEGREE J SATISFIES F TEST
!C
  24   SIGPAS = SIGJ
       JPAS = J
       NFAIL = 0
       IF (MAXDEG .EQ. J) GO TO 32
       GO TO 16
!C
!C POLYNOMIAL OF DEGREE J FAILS F TEST.  IF THERE HAVE BEEN THREE
!C SUCCESSIVE FAILURES, A STATISTICALLY BEST DEGREE HAS BEEN FOUND.
!C
  25   NFAIL = NFAIL + 1
       IF (NFAIL .GE. 3) GO TO 29
       IF (MAXDEG .EQ. J) GO TO 32
       GO TO 16
!C
!C RAISE THE DEGREE IF DEGREE  MAXDEG  HAS NOT YET BEEN REACHED  (INPUT
!C EPS = 0.)
!C
  26   IF (MAXDEG .EQ. J) GO TO 28
       GO TO 16
!C
!C SEE IF RMS ERROR CRITERION IS SATISFIED  (INPUT EPS .GT. Zero)
!C
  27   IF (SIGJ .LE. ETST) GO TO 28
       IF (MAXDEG .EQ. J) GO TO 31
       GO TO 16
!C
!C RETURNS
!C
  28   IERR = 1
       NDEG = J
       SIG = SIGJ
       GO TO 33
  29   IERR = 1
       NDEG = JPAS
       SIG = SIGPAS
       GO TO 33
  30   IERR = 2
    !  CALL XERMSG ('SLATEC', 'POLFIT', 'INVALID INPUT PARAMETER.', 2, &
    ! +   1)
       CALL Halt('SLATEC POLFIT INVALID INPUT PARAMETER  2')
       GO TO 37
  31   IERR = 3
       NDEG = MAXDEG
       SIG = SIGJ
       GO TO 33
  32   IERR = 4
       NDEG = JPAS
       SIG = SIGPAS
!C
  33   A(K3) = NDEG
!C
!C WHEN STATISTICAL TEST HAS BEEN USED, EVALUATE THE BEST POLYNOMIAL AT
!C ALL THE DATA POINTS IF  R  DOES NOT ALREADY CONTAIN THESE VALUES
!C
       IF(EPS .GE. Zero  .OR.  NDEG .EQ. MAXDEG) GO TO 36
       NDER = 0
       DO 35 I = 1,M
         CALL PVALUE (NDEG,NDER,X(I),R(I),YP,A)
  35     CONTINUE
  36   EPS = SQRT(SIG/XM)
  37   RETURN
       END SUBROUTINE
!
!---------------------------------------------------------------------
! 
!*DECK PVALUE
       SUBROUTINE PVALUE (L, NDER, X, YFIT, YP, A)
!C***BEGIN PROLOGUE  PVALUE
!C***PURPOSE  Use the coefficients generated by POLFIT to evaluate the
!C            polynomial fit of degree L, along with the first NDER of
!C            its derivatives, at a specified point.
!C***LIBRARY   SLATEC
!C***CATEGORY  K6
!C***TYPE      SINGLE PRECISION (PVALUE-S, DP1VLU-D)
!C***KEYWORDS  CURVE FITTING, LEAST SQUARES, POLYNOMIAL APPROXIMATION
!C***AUTHOR  Shampine, L. F., (SNLA)
!C           Davenport, S. M., (SNLA)
!C***DESCRIPTION
!C
!C     Written by L. F. Shampine and S. M. Davenport.
!C
!C     Abstract
!C
!C     The subroutine  PVALUE  uses the coefficients generated by  POLFIT
!C     to evaluate the polynomial fit of degree  L , along with the first
!C     NDER  of its derivatives, at a specified point.  Computationally
!C     stable recurrence relations are used to perform this task.
!C
!C     The parameters for  PVALUE  are
!C
!C     Input --
!C         L -      the degree of polynomial to be evaluated.  L  may be
!C                  any non-negative integer which is less than or equal
!C                  to  NDEG , the highest degree polynomial provided
!C                  by  POLFIT .
!C         NDER -   the number of derivatives to be evaluated.  NDER
!C                  may be 0 or any positive value.  If NDER is less
!C                  than 0, it will be treated as 0.
!C         X -      the argument at which the polynomial and its
!C                  derivatives are to be evaluated.
!C         A -      work and output array containing values from last
!C                  call to  POLFIT .
!C
!C     Output --
!C         YFIT -   value of the fitting polynomial of degree  L  at  X
!C         YP -     array containing the first through  NDER  derivatives
!C                  of the polynomial of degree  L .  YP  must be
!C                  dimensioned at least  NDER  in the calling program.
!C
!C***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
!C                 Curve fitting by polynomials in one variable, Report
!C                 SLA-74-0270, Sandia Laboratories, June 1974.
!C***ROUTINES CALLED  XERMSG
!C***REVISION HISTORY  (YYMMDD)
!C   740601  DATE WRITTEN
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890531  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  PVALUE
       DIMENSION YP(*),A(*)
       CHARACTER*8 XERN1, XERN2
!C***FIRST EXECUTABLE STATEMENT  PVALUE
       IF (L .LT. 0) GO TO 12
       NDO = MAX(NDER,0)
       NDO = MIN(NDO,L)
       MAXORD = A(1) + 0.5D0
       K1 = MAXORD + 1
       K2 = K1 + MAXORD
       K3 = K2 + MAXORD + 2
       NORD = A(K3) + 0.5D0
       IF (L .GT. NORD) GO TO 11
       K4 = K3 + L + 1
       IF (NDER .LT. 1) GO TO 2
       DO 1 I = 1,NDER
  1      YP(I) = Zero
  2    IF (L .GE. 2) GO TO 4
       IF (L .EQ. 1) GO TO 3
!C
!C L IS 0
!C
       VAL = A(K2+1)
       GO TO 10
!C
!C L IS 1
!C
  3    CC = A(K2+2)
       VAL = A(K2+1) + (X-A(2))*CC
       IF (NDER .GE. 1) YP(1) = CC
       GO TO 10
!C
!C L IS GREATER THAN 1
!C
  4    NDP1 = NDO + 1
       K3P1 = K3 + 1
       K4P1 = K4 + 1
       LP1 = L + 1
       LM1 = L - 1
       ILO = K3 + 3
       IUP = K4 + NDP1
       DO 5 I = ILO,IUP
  5      A(I) = Zero
       DIF = X - A(LP1)
       KC = K2 + LP1
       A(K4P1) = A(KC)
       A(K3P1) = A(KC-1) + DIF*A(K4P1)
       A(K3+2) = A(K4P1)
!C
!C EVALUATE RECURRENCE RELATIONS FOR FUNCTION VALUE AND DERIVATIVES
!C
       DO 9 I = 1,LM1
         IN = L - I
         INP1 = IN + 1
         K1I = K1 + INP1
         IC = K2 + IN
         DIF = X - A(INP1)
         VAL = A(IC) + DIF*A(K3P1) - A(K1I)*A(K4P1)
         IF (NDO .LE. 0) GO TO 8
         DO 6 N = 1,NDO
           K3PN = K3P1 + N
           K4PN = K4P1 + N
  6        YP(N) = DIF*A(K3PN) + N*A(K3PN-1) - A(K1I)*A(K4PN)
!C
!C SAVE VALUES NEEDED FOR NEXT EVALUATION OF RECURRENCE RELATIONS
!C
         DO 7 N = 1,NDO
           K3PN = K3P1 + N
           K4PN = K4P1 + N
           A(K4PN) = A(K3PN)
  7        A(K3PN) = YP(N)
  8      A(K4P1) = A(K3P1)
  9      A(K3P1) = VAL
!C
!C NORMAL RETURN OR ABORT DUE TO ERROR
!C
  10   YFIT = VAL
       RETURN
!C
    11 WRITE (XERN1, '(I8)') L
       WRITE (XERN2, '(I8)') NORD
    !  CALL XERMSG ('SLATEC', 'PVALUE', &
    !    'THE ORDER OF POLYNOMIAL EVALUATION, L = ' // XERN1 // &
    !    ' REQUESTED EXCEEDS THE HIGHEST ORDER FIT, NORD = ' // XERN2 //&
    !    ', COMPUTED BY POLFIT -- EXECUTION TERMINATED.', 8, 2)
       CALL Halt('Order of polinomial requested is higher than N-1')
       RETURN
!C
  ! 12 CALL XERMSG ('SLATEC', 'PVALUE', &
  !     'INVALID INPUT PARAMETER.  ORDER OF POLYNOMIAL EVALUATION ' // &
  !     'REQUESTED IS NEGATIVE -- EXECUTION TERMINATED.', 2, 2)
    12 CALL Halt('Order of polinomial requested is higher than N-1')
       RETURN
       END SUBROUTINE
!*DECK XERCNT
       SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
!C***BEGIN PROLOGUE  XERCNT
!C***SUBSIDIARY
!C***PURPOSE  Allow user control over handling of errors.
!C***LIBRARY   SLATEC (XERROR)
!C***CATEGORY  R3C
!C***TYPE      ALL (XERCNT-A)
!C***KEYWORDS  ERROR, XERROR
!C***AUTHOR  Jones, R. E., (SNLA)
!C***DESCRIPTION
!C
!C     Abstract
!C        Allows user control over handling of individual errors.
!C        Just after each message is recorded, but before it is
!C        processed any further (i.e., before it is printed or
!C        a decision to abort is made), a call is made to XERCNT.
!C        If the user has provided his own version of XERCNT, he
!C        can then override the value of KONTROL used in processing
!C        this message by redefining its value.
!C        KONTRL may be set to any value from -2 to 2.
!C        The meanings for KONTRL are the same as in XSETF, except
!C        that the value of KONTRL changes only for this message.
!C        If KONTRL is set to a value outside the range from -2 to 2,
!C        it will be moved back into that range.
!C
!C     Description of Parameters
!C
!C      --Input--
!C        LIBRAR - the library that the routine is in.
!C        SUBROU - the subroutine that XERMSG is being called from
!C        MESSG  - the first 20 characters of the error message.
!C        NERR   - same as in the call to XERMSG.
!C        LEVEL  - same as in the call to XERMSG.
!C        KONTRL - the current value of the control flag as set
!C                 by a call to XSETF.
!C
!C      --Output--
!C        KONTRL - the new value of KONTRL.  If KONTRL is not
!C                 defined, it will remain at its original value.
!C                 This changed value of control affects only
!C                 the current occurrence of the current message.
!C
!C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!C                 Error-handling Package, SAND82-0800, Sandia
!C                 Laboratories, 1982.
!C***ROUTINES CALLED  (NONE)
!C***REVISION HISTORY  (YYMMDD)
!C   790801  DATE WRITTEN
!C   861211  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900206  Routine changed from user-callable to subsidiary.  (WRB)
!C   900510  Changed calling sequence to include LIBRARY and SUBROUTINE
!C           names, changed routine name from XERCTL to XERCNT.  (RWC)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  XERCNT
       CHARACTER*(*) LIBRAR, SUBROU, MESSG
!C***FIRST EXECUTABLE STATEMENT  XERCNT
       RETURN
       END SUBROUTINE
!*DECK XERHLT
       SUBROUTINE XERHLT (MESSG)
!C***BEGIN PROLOGUE  XERHLT
!C***SUBSIDIARY
!C***PURPOSE  Abort program execution and print error message.
!C***LIBRARY   SLATEC (XERROR)
!C***CATEGORY  R3C
!C***TYPE      ALL (XERHLT-A)
!C***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
!C***AUTHOR  Jones, R. E., (SNLA)
!C***DESCRIPTION
!C
!C     Abstract
!C        ***Note*** machine dependent routine
!C        XERHLT aborts the execution of the program.
!C        The error message causing the abort is given in the calling
!C        sequence, in case one needs it for printing on a dayfile,
!C        for example.
!C
!C     Description of Parameters
!C        MESSG is as in XERMSG.
!C
!C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!C                 Error-handling Package, SAND82-0800, Sandia
!C                 Laboratories, 1982.
!C***ROUTINES CALLED  (NONE)
!C***REVISION HISTORY  (YYMMDD)
!C   790801  DATE WRITTEN
!C   861211  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900206  Routine changed from user-callable to subsidiary.  (WRB)
!C   900510  Changed calling sequence to delete length of character
!C           and changed routine name from XERABT to XERHLT.  (RWC)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  XERHLT
       CHARACTER*(*) MESSG
!C***FIRST EXECUTABLE STATEMENT  XERHLT
       STOP
       END SUBROUTINE
!*DECK XERMSG
       SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
!C***BEGIN PROLOGUE  XERMSG
!C***PURPOSE  Process error messages for SLATEC and other libraries.
!C***LIBRARY   SLATEC (XERROR)
!C***CATEGORY  R3C
!C***TYPE      ALL (XERMSG-A)
!C***KEYWORDS  ERROR MESSAGE, XERROR
!C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
!C***DESCRIPTION
!C
!C   XERMSG processes a diagnostic message in a manner determined by the
!C   value of LEVEL and the current value of the library error control
!C   flag, KONTRL.  See subroutine XSETF for details.
!C
!C    LIBRAR   A character constant (or character variable) with the name
!C             of the library.  This will be 'SLATEC' for the SLATEC
!C             Common Math Library.  The error handling package is
!C             general enough to be used by many libraries
!C             simultaneously, so it is desirable for the routine that
!C             detects and reports an error to identify the library name
!C             as well as the routine name.
!C
!C    SUBROU   A character constant (or character variable) with the name
!C             of the routine that detected the error.  Usually it is the
!C             name of the routine that is calling XERMSG.  There are
!C             some instances where a user callable library routine calls
!C             lower level subsidiary routines where the error is
!C             detected.  In such cases it may be more informative to
!C             supply the name of the routine the user called rather than
!C             the name of the subsidiary routine that detected the
!C             error.
!C
!C    MESSG    A character constant (or character variable) with the text
!C             of the error or warning message.  In the example below,
!C             the message is a character constant that contains a
!C             generic message.
!C
!C                   CALL XERMSG ('SLATEC', 'MMPY',
!C                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
!C                  *3, 1)
!C
!C             It is possible (and is sometimes desirable) to generate a
!C             specific message--e.g., one that contains actual numeric
!C             values.  Specific numeric values can be converted into
!C             character strings using formatted WRITE statements into
!C             character variables.  This is called standard Fortran
!C             internal file I/O and is exemplified in the first three
!C             lines of the following example.  You can also catenate
!C             substrings of characters to construct the error message.
!C             Here is an example showing the use of both writing to
!C             an internal file and catenating character strings.
!C
!C                   CHARACTER*5 CHARN, CHARL
!C                   WRITE (CHARN,10) N
!C                   WRITE (CHARL,10) LDA
!C                10 FORMAT(I5)
!C                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
!C                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
!C                  *   CHARL, 3, 1)
!C
!C             There are two subtleties worth mentioning.  One is that
!C             the // for character catenation is used to construct the
!C             error message so that no single character constant is
!C             continued to the next line.  This avoids confusion as to
!C             whether there are trailing blanks at the end of the line.
!C             The second is that by catenating the parts of the message
!C             as an actual argument rather than encoding the entire
!C             message into one large character variable, we avoid
!C             having to know how long the message will be in order to
!C             declare an adequate length for that large character
!C             variable.  XERMSG calls XERPRN to print the message using
!C             multiple lines if necessary.  If the message is very long,
!C             XERPRN will break it into pieces of 72 characters (as
!C             requested by XERMSG) for printing on multiple lines.
!C             Also, XERMSG asks XERPRN to prefix each line with ' *  '
!C             so that the total line length could be 76 characters.
!C             Note also that XERPRN scans the error message backwards
!C             to ignore trailing blanks.  Another feature is that
!C             the substring '$$' is treated as a new line sentinel
!C             by XERPRN.  If you want to construct a multiline
!C             message without having to count out multiples of 72
!C             characters, just use '$$' as a separator.  '$$'
!C             obviously must occur within 72 characters of the
!C             start of each line to have its intended effect since
!C             XERPRN is asked to wrap around at 72 characters in
!C             addition to looking for '$$'.
!C
!C    NERR     An integer value that is chosen by the library routine's
!C             author.  It must be in the range -99 to 999 (three
!C             printable digits).  Each distinct error should have its
!C             own error number.  These error numbers should be described
!C             in the machine readable documentation for the routine.
!C             The error numbers need be unique only within each routine,
!C             so it is reasonable for each routine to start enumerating
!C             errors from 1 and proceeding to the next integer.
!C
!C    LEVEL    An integer value in the range 0 to 2 that indicates the
!C             level (severity) of the error.  Their meanings are
!C
!C            -1  A warning message.  This is used if it is not clear
!C                that there really is an error, but the user's attention
!C                may be needed.  An attempt is made to only print this
!C                message once.
!C
!C             0  A warning message.  This is used if it is not clear
!C                that there really is an error, but the user's attention
!C                may be needed.
!C
!C             1  A recoverable error.  This is used even if the error is
!C                so serious that the routine cannot return any useful
!C                answer.  If the user has told the error package to
!C                return after recoverable errors, then XERMSG will
!C                return to the Library routine which can then return to
!C                the user's routine.  The user may also permit the error
!C                package to terminate the program upon encountering a
!C                recoverable error.
!C
!C             2  A fatal error.  XERMSG will not return to its caller
!C                after it receives a fatal error.  This level should
!C                hardly ever be used; it is much better to allow the
!C                user a chance to recover.  An example of one of the few
!C                cases in which it is permissible to declare a level 2
!C                error is a reverse communication Library routine that
!C                is likely to be called repeatedly until it integrates
!C                across some interval.  If there is a serious error in
!C                the input such that another step cannot be taken and
!C                the Library routine is called again without the input
!C                error having been corrected by the caller, the Library
!C                routine will probably be called forever with improper
!C                input.  In this case, it is reasonable to declare the
!C                error to be fatal.
!C
!C    Each of the arguments to XERMSG is input; none will be modified by
!C    XERMSG.  A routine may make multiple calls to XERMSG with warning
!C    level messages; however, after a call to XERMSG with a recoverable
!C    error, the routine should return to the user.  Do not try to call
!C    XERMSG with a second recoverable error after the first recoverable
!C    error because the error package saves the error number.  The user
!C    can retrieve this error number by calling another entry point in
!C    the error handling package and then clear the error number when
!C    recovering from the error.  Calling XERMSG in succession causes the
!C    old error number to be overwritten by the latest error number.
!C    This is considered harmless for error numbers associated with
!C    warning messages but must not be done for error numbers of serious
!C    errors.  After a call to XERMSG with a recoverable error, the user
!C    must be given a chance to call NUMXER or XERCLR to retrieve or
!C    clear the error number.
!C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!C                 Error-handling Package, SAND82-0800, Sandia
!C                 Laboratories, 1982.
!C***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
!C***REVISION HISTORY  (YYMMDD)
!C   880101  DATE WRITTEN
!C   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
!C           THERE ARE TWO BASIC CHANGES.
!C           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
!C               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
!C               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
!C               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
!C               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
!C               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
!C               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
!C               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
!C           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
!C               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
!C               OF LOWER CASE.
!C   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
!C           THE PRINCIPAL CHANGES ARE
!C           1.  CLARIFY COMMENTS IN THE PROLOGUES
!C           2.  RENAME XRPRNT TO XERPRN
!C           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
!C               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
!C               CHARACTER FOR NEW RECORDS.
!C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
!C           CLEAN UP THE CODING.
!C   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
!C           PREFIX.
!C   891013  REVISED TO CORRECT COMMENTS.
!C   891214  Prologue converted to Version 4.0 format.  (WRB)
!C   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
!C           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
!C           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
!C           XERCTL to XERCNT.  (RWC)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  XERMSG
       CHARACTER*(*) LIBRAR, SUBROU, MESSG
       CHARACTER*8 XLIBR, XSUBR
       CHARACTER*72  TEMP
       CHARACTER*20  LFIRST
!C***FIRST EXECUTABLE STATEMENT  XERMSG
       LKNTRL = J4SAVE (2, 0, .FALSE.)
       MAXMES = J4SAVE (4, 0, .FALSE.)
!C
!C       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
!C       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
!C          SHOULD BE PRINTED.
!C
!C       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
!C          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
!C          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
!C
       IF (NERR.LT.-9999999D0.OR.NERR.GT.99999999D0.OR.NERR.EQ.0 .OR. &
         LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
          CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' // &
            'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '// &
            'JOB ABORT DUE TO FATAL ERROR.', 72)
          CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
          CALL XERHLT (' ***XERMSG -- INVALID INPUT')
          RETURN
       ENDIF
!C
!C       RECORD THE MESSAGE.
!C
       I = J4SAVE (1, NERR, .TRUE.)
       CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
!C
!C       HANDLE PRINT-ONCE WARNING MESSAGES.
!C
       IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN
!C
!C       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
!C
       XLIBR  = LIBRAR
       XSUBR  = SUBROU
       LFIRST = MESSG
       LERR   = NERR
       LLEVEL = LEVEL
       CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
!C
       LKNTRL = MAX(-2, MIN(2,LKNTRL))
       MKNTRL = ABS(LKNTRL)
!C
!C       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
!C       ZERO AND THE ERROR IS NOT FATAL.
!C
       IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30
       IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30
       IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30
       IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
!C
!C       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
!C       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
!C       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
!C       IS NOT ZERO.
!C
       IF (LKNTRL .NE. 0) THEN
          TEMP(1:21) = 'MESSAGE FROM ROUTINE '
          I = MIN(LEN(SUBROU), 16)
          TEMP(22:21+I) = SUBROU(1:I)
          TEMP(22+I:33+I) = ' IN LIBRARY '
          LTEMP = 33 + I
          I = MIN(LEN(LIBRAR), 16)
          TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
          TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
          LTEMP = LTEMP + I + 1
          CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
       ENDIF
!C
!C       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
!C       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
!C       FROM EACH OF THE FOLLOWING THREE OPTIONS.
!C       1.  LEVEL OF THE MESSAGE
!C              'INFORMATIVE MESSAGE'
!C              'POTENTIALLY RECOVERABLE ERROR'
!C              'FATAL ERROR'
!C       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
!C              'PROG CONTINUES'
!C              'PROG ABORTED'
!C       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
!C           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
!C           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
!C              'TRACEBACK REQUESTED'
!C              'TRACEBACK NOT REQUESTED'
!C       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
!C       EXCEED 74 CHARACTERS.
!C       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
!C
       IF (LKNTRL .GT. 0) THEN
!C
!C       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
!C
          IF (LEVEL .LE. 0) THEN
             TEMP(1:20) = 'INFORMATIVE MESSAGE,'
             LTEMP = 20
          ELSEIF (LEVEL .EQ. 1) THEN
             TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
             LTEMP = 30
          ELSE
             TEMP(1:12) = 'FATAL ERROR,'
             LTEMP = 12
          ENDIF
!C
!C       THEN WHETHER THE PROGRAM WILL CONTINUE.
!C
          IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR. &
              (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN
             TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
             LTEMP = LTEMP + 14
          ELSE
             TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
             LTEMP = LTEMP + 16
          ENDIF
!C
!C       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
!C
          IF (LKNTRL .GT. 0) THEN
             TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
             LTEMP = LTEMP + 20
          ELSE
             TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
             LTEMP = LTEMP + 24
          ENDIF
          CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
       ENDIF
!C
!C       NOW SEND OUT THE MESSAGE.
!C
       CALL XERPRN (' *  ', -1, MESSG, 72)
!C
!C       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
!C          TRACEBACK.
!C
       IF (LKNTRL .GT. 0) THEN
          WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
          DO 10 I=16,22
             IF (TEMP(I:I) .NE. ' ') GO TO 20
    10    CONTINUE
!C
    20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
          CALL FDUMP
       ENDIF
!C
!C       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
!C
       IF (LKNTRL .NE. 0) THEN
          CALL XERPRN (' *  ', -1, ' ', 72)
          CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
          CALL XERPRN ('    ',  0, ' ', 72)
       ENDIF
!C
!C       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
!C       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
!C
    30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN
!C
!C       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
!C       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
!C       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
!C
       IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
          IF (LEVEL .EQ. 1) THEN
             CALL XERPRN &
                (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
          ELSE
             CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
          ENDIF
          CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
          CALL XERHLT (' ')
       ELSE
          CALL XERHLT (MESSG)
       ENDIF
       RETURN
       END SUBROUTINE
!*DECK XERPRN
       SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)
!C***BEGIN PROLOGUE  XERPRN
!C***SUBSIDIARY
!C***PURPOSE  Print error messages processed by XERMSG.
!C***LIBRARY   SLATEC (XERROR)
!C***CATEGORY  R3C
!C***TYPE      ALL (XERPRN-A)
!C***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
!C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
!C***DESCRIPTION
!C
!C This routine sends one or more lines to each of the (up to five)
!C logical units to which error messages are to be sent.  This routine
!C is called several times by XERMSG, sometimes with a single line to
!C print and sometimes with a (potentially very long) message that may
!C wrap around into multiple lines.
!C
!C PREFIX  Input argument of type CHARACTER.  This argument contains
!C         characters to be put at the beginning of each line before
!C         the body of the message.  No more than 16 characters of
!C         PREFIX will be used.
!C
!C NPREF   Input argument of type INTEGER.  This argument is the number
!C         of characters to use from PREFIX.  If it is negative, the
!C         intrinsic function LEN is used to determine its length.  If
!C         it is zero, PREFIX is not used.  If it exceeds 16 or if
!C         LEN(PREFIX) exceeds 16, only the first 16 characters will be
!C         used.  If NPREF is positive and the length of PREFIX is less
!C         than NPREF, a copy of PREFIX extended with blanks to length
!C         NPREF will be used.
!C
!C MESSG   Input argument of type CHARACTER.  This is the text of a
!C         message to be printed.  If it is a long message, it will be
!C         broken into pieces for printing on multiple lines.  Each line
!C         will start with the appropriate prefix and be followed by a
!C         piece of the message.  NWRAP is the number of characters per
!C         piece; that is, after each NWRAP characters, we break and
!C         start a new line.  In addition the characters '$$' embedded
!C         in MESSG are a sentinel for a new line.  The counting of
!C         characters up to NWRAP starts over for each new line.  The
!C         value of NWRAP typically used by XERMSG is 72 since many
!C         older error messages in the SLATEC Library are laid out to
!C         rely on wrap-around every 72 characters.
!C
!C NWRAP   Input argument of type INTEGER.  This gives the maximum size
!C         piece into which to break MESSG for printing on multiple
!C         lines.  An embedded '$$' ends a line, and the count restarts
!C         at the following character.  If a line break does not occur
!C         on a blank (it would split a word) that word is moved to the
!C         next line.  Values of NWRAP less than 16 will be treated as
!C         16.  Values of NWRAP greater than 132 will be treated as 132.
!C         The actual line length will be NPREF + NWRAP after NPREF has
!C         been adjusted to fall between 0 and 16 and NWRAP has been
!C         adjusted to fall between 16 and 132.
!C
!C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!C                 Error-handling Package, SAND82-0800, Sandia
!C                 Laboratories, 1982.
!C***ROUTINES CALLED  I1MACH, XGETUA
!C***REVISION HISTORY  (YYMMDD)
!C   880621  DATE WRITTEN
!C   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
!C           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
!C           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
!C           SLASH CHARACTER IN FORMAT STATEMENTS.
!C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
!C           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
!C           LINES TO BE PRINTED.
!C   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
!C           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
!C   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
!C   891214  Prologue converted to Version 4.0 format.  (WRB)
!C   900510  Added code to break messages between words.  (RWC)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  XERPRN
       CHARACTER*(*) PREFIX, MESSG
       INTEGER NPREF, NWRAP
       CHARACTER*148 CBUFF
       INTEGER IU(5), NUNIT
       CHARACTER*2 NEWLIN
       PARAMETER (NEWLIN = '$$')
!C***FIRST EXECUTABLE STATEMENT  XERPRN
       CALL XGETUA(IU,NUNIT)
!C
!C       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
!C       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
!C       ERROR MESSAGE UNIT.
!C
       N = I1MACH(4)
       DO 10 I=1,NUNIT
          IF (IU(I) .EQ. 0) IU(I) = N
    10 CONTINUE
!C
!C       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
!C       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
!C       THE REST OF THIS ROUTINE.
!C
       IF ( NPREF .LT. 0 ) THEN
          LPREF = LEN(PREFIX)
       ELSE
          LPREF = NPREF
       ENDIF
       LPREF = MIN(16, LPREF)
       IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX
!C
!C       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
!C       TIME FROM MESSG TO PRINT ON ONE LINE.
!C
       LWRAP = MAX(16, MIN(132, NWRAP))
!C
!C       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
!C
       LENMSG = LEN(MESSG)
       N = LENMSG
       DO 20 I=1,N
          IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30
          LENMSG = LENMSG - 1
    20 CONTINUE
    30 CONTINUE
!C
!C       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
!C
       IF (LENMSG .EQ. 0) THEN
          CBUFF(LPREF+1:LPREF+1) = ' '
          DO 40 I=1,NUNIT
             WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
    40    CONTINUE
          RETURN
       ENDIF
!C
!C       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
!C       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
!C       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
!C       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
!C
!C       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
!C       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
!C       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
!C       OF THE SECOND ARGUMENT.
!C
!C       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
!C       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
!C       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
!C       POSITION NEXTC.
!C
!C       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
!C                       REMAINDER OF THE CHARACTER STRING.  LPIECE
!C                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
!C                       WHICHEVER IS LESS.
!C
!C       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
!C                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
!C                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
!C                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
!C                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
!C                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
!C                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
!C                       SHOULD BE INCREMENTED BY 2.
!C
!C       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
!C
!C       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
!C                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
!C                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
!C                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
!C                       AT THE END OF A LINE.
!C
       NEXTC = 1
    50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
       IF (LPIECE .EQ. 0) THEN
!C
!C       THERE WAS NO NEW LINE SENTINEL FOUND.
!C
          IDELTA = 0
          LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
          IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
             DO 52 I=LPIECE+1,2,-1
                IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                   LPIECE = I-1
                   IDELTA = 1
                   GOTO 54
                ENDIF
    52       CONTINUE
          ENDIF
    54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
          NEXTC = NEXTC + LPIECE + IDELTA
       ELSEIF (LPIECE .EQ. 1) THEN
!C
!C       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
!C       DON'T PRINT A BLANK LINE.
!C
          NEXTC = NEXTC + 2
          GO TO 50
       ELSEIF (LPIECE .GT. LWRAP+1) THEN
!C
!C       LPIECE SHOULD BE SET DOWN TO LWRAP.
!C
          IDELTA = 0
          LPIECE = LWRAP
          DO 56 I=LPIECE+1,2,-1
             IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                LPIECE = I-1
                IDELTA = 1
                GOTO 58
             ENDIF
    56    CONTINUE
    58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
          NEXTC = NEXTC + LPIECE + IDELTA
       ELSE
!C
!C       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
!C       WE SHOULD DECREMENT LPIECE BY ONE.
!C
          LPIECE = LPIECE - 1
          CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
          NEXTC  = NEXTC + LPIECE + 2
       ENDIF
!C
!C       PRINT
!C
       DO 60 I=1,NUNIT
          WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
    60 CONTINUE
!C
       IF (NEXTC .LE. LENMSG) GO TO 50
       RETURN
       END SUBROUTINE
!*DECK XERSVE
       SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, &
         ICOUNT)
!C***BEGIN PROLOGUE  XERSVE
!C***SUBSIDIARY
!C***PURPOSE  Record that an error has occurred.
!C***LIBRARY   SLATEC (XERROR)
!C***CATEGORY  R3
!C***TYPE      ALL (XERSVE-A)
!C***KEYWORDS  ERROR, XERROR
!C***AUTHOR  Jones, R. E., (SNLA)
!C***DESCRIPTION
!C
!C *Usage:
!C
!C        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
!C        CHARACTER * (len) LIBRAR, SUBROU, MESSG
!C
!C        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
!C
!C *Arguments:
!C
!C        LIBRAR :IN    is the library that the message is from.
!C        SUBROU :IN    is the subroutine that the message is from.
!C        MESSG  :IN    is the message to be saved.
!C        KFLAG  :IN    indicates the action to be performed.
!C                      when KFLAG > 0, the message in MESSG is saved.
!C                      when KFLAG=0 the tables will be dumped and
!C                      cleared.
!C                      when KFLAG < 0, the tables will be dumped and
!C                      not cleared.
!C        NERR   :IN    is the error number.
!C        LEVEL  :IN    is the error severity.
!C        ICOUNT :OUT   the number of times this message has been seen,
!C                      or zero if the table has overflowed and does not
!C                      contain this message specifically.  When KFLAG=0,
!C                      ICOUNT will not be altered.
!C
!C *Description:
!C
!C   Record that this error occurred and possibly dump and clear the
!C   tables.
!C
!C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!C                 Error-handling Package, SAND82-0800, Sandia
!C                 Laboratories, 1982.
!C***ROUTINES CALLED  I1MACH, XGETUA
!C***REVISION HISTORY  (YYMMDD)
!C   800319  DATE WRITTEN
!C   861211  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900413  Routine modified to remove reference to KFLAG.  (WRB)
!C   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
!C           sequence, use IF-THEN-ELSE, make number of saved entries
!C           easily changeable, changed routine name from XERSAV to
!C           XERSVE.  (RWC)
!C   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  XERSVE
       PARAMETER (LENTAB=10)
       INTEGER LUN(5)
       CHARACTER*(*) LIBRAR, SUBROU, MESSG
       CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
       CHARACTER*20 MESTAB(LENTAB), MES
       DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
       SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
       DATA KOUNTX/0/, NMSG/0/
!C***FIRST EXECUTABLE STATEMENT  XERSVE
!C
       IF (KFLAG.LE.0) THEN
!C
!C        Dump the table.
!C
          IF (NMSG.EQ.0) RETURN
!C
!C        Print to each unit.
!C
          CALL XGETUA (LUN, NUNIT)
          DO 20 KUNIT = 1,NUNIT
             IUNIT = LUN(KUNIT)
             IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
!C
!C           Print the table header.
!C
             WRITE (IUNIT,9000)
!C
!C           Print body of table.
!C
             DO 10 I = 1,NMSG
                WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I), &
                   NERTAB(I),LEVTAB(I),KOUNT(I)
    10       CONTINUE
!C
!C           Print number of other errors.
!C
             IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX
             WRITE (IUNIT,9030)
    20    CONTINUE
!C
!C        Clear the error tables.
!C
          IF (KFLAG.EQ.0) THEN
             NMSG = 0
             KOUNTX = 0
          ENDIF
       ELSE
!C
!C        PROCESS A MESSAGE...
!C        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
!C        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
!C
          LIB = LIBRAR
          SUB = SUBROU
          MES = MESSG
          DO 30 I = 1,NMSG
             IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND. &
                MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND. &
                LEVEL.EQ.LEVTAB(I)) THEN
                   KOUNT(I) = KOUNT(I) + 1
                   ICOUNT = KOUNT(I)
                   RETURN
             ENDIF
    30    CONTINUE
!C
          IF (NMSG.LT.LENTAB) THEN
!C
!C           Empty slot found for new message.
!C
             NMSG = NMSG + 1
             LIBTAB(I) = LIB
             SUBTAB(I) = SUB
             MESTAB(I) = MES
             NERTAB(I) = NERR
             LEVTAB(I) = LEVEL
             KOUNT (I) = 1
             ICOUNT    = 1
          ELSE
!C
!C           Table is full.
!C
             KOUNTX = KOUNTX+1
             ICOUNT = 0
          ENDIF
       ENDIF
       RETURN
!C
!C     Formats.
!C
  9000 FORMAT ('0          ERROR MESSAGE SUMMARY' / &
          ' LIBRARY    SUBROUTINE MESSAGE START             NERR', &
          '     LEVEL     COUNT')
  9010 FORMAT (1X,A,3X,A,3X,A,3I10)
  9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
  9030 FORMAT (1X)
       END SUBROUTINE
!*DECK XGETUA
       SUBROUTINE XGETUA (IUNITA, N)
!C***BEGIN PROLOGUE  XGETUA
!C***PURPOSE  Return unit number(s) to which error messages are being
!C            sent.
!C***LIBRARY   SLATEC (XERROR)
!C***CATEGORY  R3C
!C***TYPE      ALL (XGETUA-A)
!C***KEYWORDS  ERROR, XERROR
!C***AUTHOR  Jones, R. E., (SNLA)
!C***DESCRIPTION
!C
!C     Abstract
!C        XGETUA may be called to determine the unit number or numbers
!C        to which error messages are being sent.
!C        These unit numbers may have been set by a call to XSETUN,
!C        or a call to XSETUA, or may be a default value.
!C
!C     Description of Parameters
!C      --Output--
!C        IUNIT - an array of one to five unit numbers, depending
!C                on the value of N.  A value of zero refers to the
!C                default unit, as defined by the I1MACH machine
!C                constant routine.  Only IUNIT(1),...,IUNIT(N) are
!C                defined by XGETUA.  The values of IUNIT(N+1),...,
!C                IUNIT(5) are not defined (for N .LT. 5) or altered
!C                in any way by XGETUA.
!C        N     - the number of units to which copies of the
!C                error messages are being sent.  N will be in the
!C                range from 1 to 5.
!C
!C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!C                 Error-handling Package, SAND82-0800, Sandia
!C                 Laboratories, 1982.
!C***ROUTINES CALLED  J4SAVE
!C***REVISION HISTORY  (YYMMDD)
!C   790801  DATE WRITTEN
!C   861211  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  XGETUA
       DIMENSION IUNITA(5)
!C***FIRST EXECUTABLE STATEMENT  XGETUA
       N = J4SAVE(5,0,.FALSE.)
       DO 30 I=1,N
          INDEX = I+4
          IF (I.EQ.1) INDEX = 3
          IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
    30 CONTINUE
       RETURN
       END SUBROUTINE
!
!--------------------------------------------------------------------
!
!*DECK PCOEF
      SUBROUTINE PCOEF (L, C, TC, A)
!C***BEGIN PROLOGUE  PCOEF
!C***PURPOSE  Convert the POLFIT coefficients to Taylor series form.
!C***LIBRARY   SLATEC
!C***CATEGORY  K1A1A2
!C***TYPE      SINGLE PRECISION (PCOEF-S, DPCOEF-D)
!C***KEYWORDS  CURVE FITTING, DATA FITTING, LEAST SQUARES, POLYNOMIAL FIT
!C***AUTHOR  Shampine, L. F., (SNLA)
!C           Davenport, S. M., (SNLA)
!C***DESCRIPTION
!C
!C     Written BY L. F. Shampine and S. M. Davenport.
!C
!C     Abstract
!C
!C     POLFIT  computes the least squares polynomial fit of degree  L  as
!C     a sum of orthogonal polynomials.  PCOEF  changes this fit to its
!C     Taylor expansion about any point  C , i.e. writes the polynomial
!C     as a sum of powers of (X-C).  Taking  C=0.  gives the polynomial
!C     in powers of X, but a suitable non-zero  C  often leads to
!C     polynomials which are better scaled and more accurately evaluated.
!C
!C     The parameters for  PCOEF  are
!C
!C     INPUT --
!C         L -      Indicates the degree of polynomial to be changed to
!C                  its Taylor expansion.  To obtain the Taylor
!C                  coefficients in reverse order, input  L  as the
!C                  negative of the degree desired.  The absolute value
!C                  of L  must be less than or equal to NDEG, the highest
!C                  degree polynomial fitted by  POLFIT .
!C         C -      The point about which the Taylor expansion is to be
!C                  made.
!C         A -      Work and output array containing values from last
!C                  call to  POLFIT .
!C
!C     OUTPUT --
!C         TC -     Vector containing the first LL+1 Taylor coefficients
!C                  where LL=ABS(L).  If  L.GT.0 , the coefficients are
!C                  in the usual Taylor series order, i.e.
!C                    P(X) = TC(1) + TC(2)*(X-C) + ... + TC(N+1)*(X-C)**N
!C                  If L .LT. 0, the coefficients are in reverse order,
!C                  i.e.
!C                    P(X) = TC(1)*(X-C)**N + ... + TC(N)*(X-C) + TC(N+1)
!C
!C***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
!C                 Curve fitting by polynomials in one variable, Report
!C                 SLA-74-0270, Sandia Laboratories, June 1974.
!C***ROUTINES CALLED  PVALUE
!C***REVISION HISTORY  (YYMMDD)
!C   740601  DATE WRITTEN
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890531  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  PCOEF
!C
      DIMENSION A(*), TC(*)
!C***FIRST EXECUTABLE STATEMENT  PCOEF
      LL = ABS(L)
      LLP1 = LL + 1
      CALL PVALUE (LL,LL,C,TC(1),TC(2),A)
      IF (LL .LT. 2) GO TO 2
      FAC = One    
      DO 1 I = 3,LLP1
        FAC = FAC*(I-1)
 1      TC(I) = TC(I)/FAC
 2    IF (L .GE. 0) GO TO 4
      NR = LLP1/2
      LLP2 = LL + 2
      DO 3 I = 1,NR
        SAVE = TC(I)
        NEW = LLP2 - I
        TC(I) = TC(NEW)
 3      TC(NEW) = SAVE
 4    RETURN
      END SUBROUTINE
!!
!!--------------------------------------------------------------------
!!
!*DECK D1MACH
      DOUBLE PRECISION FUNCTION D1MACH (I)
!C***BEGIN PROLOGUE  D1MACH
!C***PURPOSE  Return floating point machine dependent constants.
!C***LIBRARY   SLATEC
!C***CATEGORY  R1
!C***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D)
!C***KEYWORDS  MACHINE CONSTANTS
!C***AUTHOR  Fox, P. A., (Bell Labs)
!C           Hall, A. D., (Bell Labs)
!C           Schryer, N. L., (Bell Labs)
!C***DESCRIPTION
!C
!C   D1MACH can be used to obtain machine-dependent parameters for the
!C   local machine environment.  It is a function subprogram with one
!C   (input) argument, and can be referenced as follows:
!C
!C        D = D1MACH(I)
!C
!C   where I=1,...,5.  The (output) value of D above is determined by
!C   the (input) value of I.  The results for various values of I are
!C   discussed below.
!C
!C   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
!C   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!C   D1MACH( 3) = B**(-T), the smallest relative spacing.
!C   D1MACH( 4) = B**(1-T), the largest relative spacing.
!C   D1MACH( 5) = LOG10(B)
!C
!C   Assume double precision numbers are represented in the T-digit,
!C   base-B form
!C
!C              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!C
!C   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
!C   EMIN .LE. E .LE. EMAX.
!C
!C   The values of B, T, EMIN and EMAX are provided in I1MACH as
!C   follows:
!C   I1MACH(10) = B, the base.
!C   I1MACH(14) = T, the number of base-B digits.
!C   I1MACH(15) = EMIN, the smallest exponent E.
!C   I1MACH(16) = EMAX, the largest exponent E.
!C
!C   To alter this function for a particular environment, the desired
!C   set of DATA statements should be activated by removing the C from
!C   column 1.  Also, the values of D1MACH(1) - D1MACH(4) should be
!C   checked for consistency with the local operating system.
!C
!C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!C                 a portable library, ACM Transactions on Mathematical
!C                 Software 4, 2 (June 1978), pp. 177-188.
!C***ROUTINES CALLED  XERMSG
!C***REVISION HISTORY  (YYMMDD)
!C   750101  DATE WRITTEN
!C   890213  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!C   900618  Added DEC RISC constants.  (WRB)
!C   900723  Added IBM RS 6000 constants.  (WRB)
!C   900911  Added SUN 386i constants.  (WRB)
!C   910710  Added HP 730 constants.  (SMR)
!C   911114  Added Convex IEEE constants.  (WRB)
!C   920121  Added SUN -r8 compiler option constants.  (WRB)
!C   920229  Added Touchstone Delta i860 constants.  (WRB)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C   920625  Added CONVEX -p8 and -pd8 compiler option constants.
!C           (BKS, WRB)
!C   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
!C***END PROLOGUE  D1MACH
!C
      INTEGER SMALL(4)
      INTEGER LARGE(4)
      INTEGER RIGHT(4)
      INTEGER DIVER(4)
      INTEGER LOG10(4)
!C
      DOUBLE PRECISION DMACH(5)
      SAVE DMACH
!C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
!C
!C     MACHINE CONSTANTS FOR THE AMIGA
!C     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
!C
!C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
!C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
!C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
!C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
!C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
!C
!C     MACHINE CONSTANTS FOR THE AMIGA
!C     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
!C
!C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
!C     DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' /
!C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
!C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
!C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
!C
!C     MACHINE CONSTANTS FOR THE APOLLO
!C
!C     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 /
!C     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF /
!C     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 /
!C     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 /
!C     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF /
!C
!C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
!C
!C     DATA SMALL(1) / ZC00800000 /
!C     DATA SMALL(2) / Z000000000 /
!C     DATA LARGE(1) / ZDFFFFFFFF /
!C     DATA LARGE(2) / ZFFFFFFFFF /
!C     DATA RIGHT(1) / ZCC5800000 /
!C     DATA RIGHT(2) / Z000000000 /
!C     DATA DIVER(1) / ZCC6800000 /
!C     DATA DIVER(2) / Z000000000 /
!C     DATA LOG10(1) / ZD00E730E7 /
!C     DATA LOG10(2) / ZC77800DC0 /
!C
!C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
!C
!C     DATA SMALL(1) / O1771000000000000 /
!C     DATA SMALL(2) / O0000000000000000 /
!C     DATA LARGE(1) / O0777777777777777 /
!C     DATA LARGE(2) / O0007777777777777 /
!C     DATA RIGHT(1) / O1461000000000000 /
!C     DATA RIGHT(2) / O0000000000000000 /
!C     DATA DIVER(1) / O1451000000000000 /
!C     DATA DIVER(2) / O0000000000000000 /
!C     DATA LOG10(1) / O1157163034761674 /
!C     DATA LOG10(2) / O0006677466732724 /
!C
!C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
!C
!C     DATA SMALL(1) / O1771000000000000 /
!C     DATA SMALL(2) / O7770000000000000 /
!C     DATA LARGE(1) / O0777777777777777 /
!C     DATA LARGE(2) / O7777777777777777 /
!C     DATA RIGHT(1) / O1461000000000000 /
!C     DATA RIGHT(2) / O0000000000000000 /
!C     DATA DIVER(1) / O1451000000000000 /
!C     DATA DIVER(2) / O0000000000000000 /
!C     DATA LOG10(1) / O1157163034761674 /
!C     DATA LOG10(2) / O0006677466732724 /
!C
!C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
!C
!C     DATA SMALL(1) / Z"3001800000000000" /
!C     DATA SMALL(2) / Z"3001000000000000" /
!C     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" /
!C     DATA LARGE(2) / Z"4FFE000000000000" /
!C     DATA RIGHT(1) / Z"3FD2800000000000" /
!C     DATA RIGHT(2) / Z"3FD2000000000000" /
!C     DATA DIVER(1) / Z"3FD3800000000000" /
!C     DATA DIVER(2) / Z"3FD3000000000000" /
!C     DATA LOG10(1) / Z"3FFF9A209A84FBCF" /
!C     DATA LOG10(2) / Z"3FFFF7988F8959AC" /
!C
!C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
!C
!C     DATA SMALL(1) / 00564000000000000000B /
!C     DATA SMALL(2) / 00000000000000000000B /
!C     DATA LARGE(1) / 37757777777777777777B /
!C     DATA LARGE(2) / 37157777777777777777B /
!C     DATA RIGHT(1) / 15624000000000000000B /
!C     DATA RIGHT(2) / 00000000000000000000B /
!C     DATA DIVER(1) / 15634000000000000000B /
!C     DATA DIVER(2) / 00000000000000000000B /
!C     DATA LOG10(1) / 17164642023241175717B /
!C     DATA LOG10(2) / 16367571421742254654B /
!C
!C     MACHINE CONSTANTS FOR THE CELERITY C1260
!C
!C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
!C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
!C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
!C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
!C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
!C
!C     MACHINE CONSTANTS FOR THE CONVEX
!C     USING THE -fn OR -pd8 COMPILER OPTION
!C
!C     DATA DMACH(1) / Z'0010000000000000' /
!C     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFF' /
!C     DATA DMACH(3) / Z'3CC0000000000000' /
!C     DATA DMACH(4) / Z'3CD0000000000000' /
!C     DATA DMACH(5) / Z'3FF34413509F79FF' /
!C
!C     MACHINE CONSTANTS FOR THE CONVEX
!C     USING THE -fi COMPILER OPTION
!C
!C     DATA DMACH(1) / Z'0010000000000000' /
!C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!C     DATA DMACH(3) / Z'3CA0000000000000' /
!C     DATA DMACH(4) / Z'3CB0000000000000' /
!C     DATA DMACH(5) / Z'3FD34413509F79FF' /
!C
!C     MACHINE CONSTANTS FOR THE CONVEX
!C     USING THE -p8 COMPILER OPTION
!C
!C     DATA DMACH(1) / Z'00010000000000000000000000000000' /
!C     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
!C     DATA DMACH(3) / Z'3F900000000000000000000000000000' /
!C     DATA DMACH(4) / Z'3F910000000000000000000000000000' /
!C     DATA DMACH(5) / Z'3FFF34413509F79FEF311F12B35816F9' /
!C
!C     MACHINE CONSTANTS FOR THE CRAY
!C
!C     DATA SMALL(1) / 201354000000000000000B /
!C     DATA SMALL(2) / 000000000000000000000B /
!C     DATA LARGE(1) / 577767777777777777777B /
!C     DATA LARGE(2) / 000007777777777777774B /
!C     DATA RIGHT(1) / 376434000000000000000B /
!C     DATA RIGHT(2) / 000000000000000000000B /
!C     DATA DIVER(1) / 376444000000000000000B /
!C     DATA DIVER(2) / 000000000000000000000B /
!C     DATA LOG10(1) / 377774642023241175717B /
!C     DATA LOG10(2) / 000007571421742254654B /
!C
!C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
!C     STATIC DMACH(5)
!C
!C     DATA SMALL /    20K, 3*0 /
!C     DATA LARGE / 77777K, 3*177777K /
!C     DATA RIGHT / 31420K, 3*0 /
!C     DATA DIVER / 32020K, 3*0 /
!C     DATA LOG10 / 40423K, 42023K, 50237K, 74776K /
!C
!C     MACHINE CONSTANTS FOR THE DEC ALPHA
!C     USING G_FLOAT
!C
!C     DATA DMACH(1) / '0000000000000010'X /
!C     DATA DMACH(2) / 'FFFFFFFFFFFF7FFF'X /
!C     DATA DMACH(3) / '0000000000003CC0'X /
!C     DATA DMACH(4) / '0000000000003CD0'X /
!C     DATA DMACH(5) / '79FF509F44133FF3'X /
!C
!C     MACHINE CONSTANTS FOR THE DEC ALPHA
!C     USING IEEE_FORMAT
!C
!C     DATA DMACH(1) / '0010000000000000'X /
!C     DATA DMACH(2) / '7FEFFFFFFFFFFFFF'X /
!C     DATA DMACH(3) / '3CA0000000000000'X /
!C     DATA DMACH(4) / '3CB0000000000000'X /
!C     DATA DMACH(5) / '3FD34413509F79FF'X /
!C
!C     MACHINE CONSTANTS FOR THE DEC RISC
!C
!C     DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000'/
!C     DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF'/
!C     DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000'/
!C     DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000'/
!C     DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413'/
!C
!C     MACHINE CONSTANTS FOR THE DEC VAX
!C     USING D_FLOATING
!C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
!C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
!C
!C     DATA SMALL(1), SMALL(2) /        128,           0 /
!C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
!C     DATA RIGHT(1), RIGHT(2) /       9344,           0 /
!C     DATA DIVER(1), DIVER(2) /       9472,           0 /
!C     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
!C
!C     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
!C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
!C     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
!C     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
!C     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
!C
!C     MACHINE CONSTANTS FOR THE DEC VAX
!C     USING G_FLOATING
!C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
!C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
!C
!C     DATA SMALL(1), SMALL(2) /         16,           0 /
!C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
!C     DATA RIGHT(1), RIGHT(2) /      15552,           0 /
!C     DATA DIVER(1), DIVER(2) /      15568,           0 /
!C     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
!C
!C     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
!C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
!C     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
!C     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
!C     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
!C
!C     MACHINE CONSTANTS FOR THE ELXSI 6400
!C     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION)
!C
!C     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
!C     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X /
!C     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X /
!C     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X /
!C     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X /
!C
!C     MACHINE CONSTANTS FOR THE HARRIS 220
!C
!C     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
!C     DATA LARGE(1), LARGE(2) / '37777777, '37777577 /
!C     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 /
!C     DATA DIVER(1), DIVER(2) / '20000000, '00000334 /
!C     DATA LOG10(1), LOG10(2) / '23210115, '10237777 /
!C
!C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
!C
!C     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 /
!C     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 /
!C     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 /
!C     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 /
!C     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 /
!C
!C     MACHINE CONSTANTS FOR THE HP 730
!C
!C     DATA DMACH(1) / Z'0010000000000000' /
!C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!C     DATA DMACH(3) / Z'3CA0000000000000' /
!C     DATA DMACH(4) / Z'3CB0000000000000' /
!C     DATA DMACH(5) / Z'3FD34413509F79FF' /
!C
!C     MACHINE CONSTANTS FOR THE HP 2100
!C     THREE WORD DOUBLE PRECISION OPTION WITH FTN4
!C
!C     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /
!C     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /
!C     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /
!C     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /
!C     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /
!C
!C     MACHINE CONSTANTS FOR THE HP 2100
!C     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4
!C
!C     DATA SMALL(1), SMALL(2) /  40000B,       0 /
!C     DATA SMALL(3), SMALL(4) /       0,       1 /
!C     DATA LARGE(1), LARGE(2) /  77777B, 177777B /
!C     DATA LARGE(3), LARGE(4) / 177777B, 177776B /
!C     DATA RIGHT(1), RIGHT(2) /  40000B,       0 /
!C     DATA RIGHT(3), RIGHT(4) /       0,    225B /
!C     DATA DIVER(1), DIVER(2) /  40000B,       0 /
!C     DATA DIVER(3), DIVER(4) /       0,    227B /
!C     DATA LOG10(1), LOG10(2) /  46420B,  46502B /
!C     DATA LOG10(3), LOG10(4) /  76747B, 176377B /
!C
!C     MACHINE CONSTANTS FOR THE HP 9000
!C
!C     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B /
!C     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B /
!C     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B /
!C     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B /
!C     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B /
!C
!C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
!C     THE PERKIN ELMER (INTERDATA) 7/32.
!C
!C     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 /
!C     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
!C     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 /
!C     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 /
!C     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF /
!C
!C     MACHINE CONSTANTS FOR THE IBM PC
!C     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION
!C     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087.
!C
!C     DATA SMALL(1) / 2.23D-308  /
!C     DATA LARGE(1) / 1.79D+308  /
!C     DATA RIGHT(1) / 1.11D-16   /
!C     DATA DIVER(1) / 2.22D-16   /
!C     DATA LOG10(1) / 0.301029995663981195D0 /
!C
!C     MACHINE CONSTANTS FOR THE IBM RS 6000
!C
!C     DATA DMACH(1) / Z'0010000000000000' /
!C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!C     DATA DMACH(3) / Z'3CA0000000000000' /
!C     DATA DMACH(4) / Z'3CB0000000000000' /
!C     DATA DMACH(5) / Z'3FD34413509F79FF' /
!C
!C     MACHINE CONSTANTS FOR THE INTEL i860
!C
!C     DATA DMACH(1) / Z'0010000000000000' /
!C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!C     DATA DMACH(3) / Z'3CA0000000000000' /
!C     DATA DMACH(4) / Z'3CB0000000000000' /
!C     DATA DMACH(5) / Z'3FD34413509F79FF' /
!C
!C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
!C
!C     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 /
!C     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 /
!C     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 /
!C     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 /
!C     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 /
!C
!C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
!C
!C     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 /
!C     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 /
!C     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 /
!C     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 /
!C     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 /
!C
!C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
!C
!C     DATA SMALL(1), SMALL(2) /    8388608,           0 /
!C     DATA LARGE(1), LARGE(2) / 2147483647,          -1 /
!C     DATA RIGHT(1), RIGHT(2) /  612368384,           0 /
!C     DATA DIVER(1), DIVER(2) /  620756992,           0 /
!C     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 /
!C
!C     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 /
!C     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 /
!C     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 /
!C     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 /
!C     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 /
!C
!C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
!C
!C     DATA SMALL(1), SMALL(2) /    128,      0 /
!C     DATA SMALL(3), SMALL(4) /      0,      0 /
!C     DATA LARGE(1), LARGE(2) /  32767,     -1 /
!C     DATA LARGE(3), LARGE(4) /     -1,     -1 /
!C     DATA RIGHT(1), RIGHT(2) /   9344,      0 /
!C     DATA RIGHT(3), RIGHT(4) /      0,      0 /
!C     DATA DIVER(1), DIVER(2) /   9472,      0 /
!C     DATA DIVER(3), DIVER(4) /      0,      0 /
!C     DATA LOG10(1), LOG10(2) /  16282,   8346 /
!C     DATA LOG10(3), LOG10(4) / -31493, -12296 /
!C
!C     DATA SMALL(1), SMALL(2) / O000200, O000000 /
!C     DATA SMALL(3), SMALL(4) / O000000, O000000 /
!C     DATA LARGE(1), LARGE(2) / O077777, O177777 /
!C     DATA LARGE(3), LARGE(4) / O177777, O177777 /
!C     DATA RIGHT(1), RIGHT(2) / O022200, O000000 /
!C     DATA RIGHT(3), RIGHT(4) / O000000, O000000 /
!C     DATA DIVER(1), DIVER(2) / O022400, O000000 /
!C     DATA DIVER(3), DIVER(4) / O000000, O000000 /
!C     DATA LOG10(1), LOG10(2) / O037632, O020232 /
!C     DATA LOG10(3), LOG10(4) / O102373, O147770 /
!C
!C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
!C
!C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
!C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
!C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
!C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
!C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
!C
!C     MACHINE CONSTANTS FOR THE SUN
!C
!C     DATA DMACH(1) / Z'0010000000000000' /
!C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!C     DATA DMACH(3) / Z'3CA0000000000000' /
!C     DATA DMACH(4) / Z'3CB0000000000000' /
!C     DATA DMACH(5) / Z'3FD34413509F79FF' /
!C
!C     MACHINE CONSTANTS FOR THE SUN
!C     USING THE -r8 COMPILER OPTION
!C
!C     DATA DMACH(1) / Z'00010000000000000000000000000000' /
!C     DATA DMACH(2) / Z'7FFEFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
!C     DATA DMACH(3) / Z'3F8E0000000000000000000000000000' /
!C     DATA DMACH(4) / Z'3F8F0000000000000000000000000000' /
!C     DATA DMACH(5) / Z'3FFD34413509F79FEF311F12B35816F9' /
!C
!C     MACHINE CONSTANTS FOR THE SUN 386i
!C
!C     DATA SMALL(1), SMALL(2) / Z'FFFFFFFD', Z'000FFFFF' /
!C     DATA LARGE(1), LARGE(2) / Z'FFFFFFB0', Z'7FEFFFFF' /
!C     DATA RIGHT(1), RIGHT(2) / Z'000000B0', Z'3CA00000' /
!C     DATA DIVER(1), DIVER(2) / Z'FFFFFFCB', Z'3CAFFFFF'
!C     DATA LOG10(1), LOG10(2) / Z'509F79E9', Z'3FD34413' /
!C
!C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
!C
!C     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 /
!C     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 /
!C     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 /
!C     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 /
!C     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 /
!C
!C***FIRST EXECUTABLE STATEMENT  D1MACH
      IF (I .LT. 1 .OR. I .GT. 5) CALL XERMSG ('SLATEC', 'D1MACH', &
          'I OUT OF BOUNDS', 1, 2)
!C
      D1MACH = DMACH(I)
      RETURN
!C
      END
!*DECK DCKDER
!
!----------------------------------------------------------------------
!
      SUBROUTINE DCKDER (M, N, X, FVEC, FJAC, LDFJAC, XP, FVECP, MODE,&
                         ERR)
!C***BEGIN PROLOGUE  DCKDER
!C***PURPOSE  Check the gradients of M nonlinear functions in N
!C            variables, evaluated at a point X, for consistency
!C            with the functions themselves.
!C***LIBRARY   SLATEC
!C***CATEGORY  F3, G4C
!C***TYPE      DOUBLE PRECISION (CHKDER-S, DCKDER-D)
!C***KEYWORDS  GRADIENTS, JACOBIAN, MINPACK, NONLINEAR
!C***AUTHOR  Hiebert, K. L. (SNLA)
!C***DESCRIPTION
!C
!C   This subroutine is a companion routine to DNSQ and DNSQE. It may
!C   be used to check the coding of the Jacobian calculation.
!C
!C     SUBROUTINE DCKDER
!C
!C     This subroutine checks the gradients of M nonlinear functions
!C     in N variables, evaluated at a point X, for consistency with
!C     the functions themselves. The user must call DCKDER twice,
!C     first with MODE = 1 and then with MODE = 2.
!C
!C     MODE = 1. On input, X must contain the point of evaluation.
!C               On output, XP is set to a neighboring point.
!C
!C     MODE = 2. On input, FVEC must contain the functions and the
!C                         rows of FJAC must contain the gradients
!C                         of the respective functions each evaluated
!C                         at X, and FVECP must contain the functions
!C                         evaluated at XP.
!C               On output, ERR contains measures of correctness of
!C                          the respective gradients.
!C
!C     The subroutine does not perform reliably if cancellation or
!C     rounding errors cause a severe loss of significance in the
!C     evaluation of a function. Therefore, none of the components
!C     of X should be unusually small (in particular, zero) or any
!C     other value which may cause loss of significance.
!C
!C     The SUBROUTINE statement is
!C
!C       SUBROUTINE DCKDER(M,N,X,FVEC,FJAC,LDFJAC,XP,FVECP,MODE,ERR)
!C
!C     where
!C
!C       M is a positive integer input variable set to the number
!C         of functions.
!C
!C       N is a positive integer input variable set to the number
!C         of variables.
!C
!C       X is an input array of length N.
!C
!C       FVEC is an array of length M. On input when MODE = 2,
!C         FVEC must contain the functions evaluated at X.
!C
!C       FJAC is an M by N array. On input when MODE = 2,
!C         the rows of FJAC must contain the gradients of
!C         the respective functions evaluated at X.
!C
!C       LDFJAC is a positive integer input parameter not less than M
!C         which specifies the leading dimension of the array FJAC.
!C
!C       XP is an array of length N. On output when MODE = 1,
!C         XP is set to a neighboring point of X.
!C
!C       FVECP is an array of length M. On input when MODE = 2,
!C         FVECP must contain the functions evaluated at XP.
!C
!C       MODE is an integer input variable set to 1 on the first call
!C         and 2 on the second. Other values of MODE are equivalent
!C         to MODE = 1.
!C
!C       ERR is an array of length M. On output when MODE = 2,
!C         ERR contains measures of correctness of the respective
!C         gradients. If there is no severe loss of significance,
!C         then if ERR(I) is 1.0 the I-th gradient is correct,
!C         while if ERR(I) is 0.0 the I-th gradient is incorrect.
!C         For values of ERR between 0.0 and 1.0, the categorization
!C         is less certain. In general, a value of ERR(I) greater
!C         than 0.5 indicates that the I-th gradient is probably
!C         correct, while a value of ERR(I) less than 0.5 indicates
!C         that the I-th gradient is probably incorrect.
!C
!C***REFERENCES  M. J. D. Powell, A hybrid method for nonlinear equa-
!C                 tions. In Numerical Methods for Nonlinear Algebraic
!C                 Equations, P. Rabinowitz, Editor.  Gordon and Breach,
!C                 1988.
!C***ROUTINES CALLED  D1MACH
!C***REVISION HISTORY  (YYMMDD)
!C   800301  DATE WRITTEN
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890831  Modified array declarations.  (WRB)
!C   890831  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900326  Removed duplicate information from DESCRIPTION section.
!C           (WRB)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  DCKDER
      INTEGER I, J, LDFJAC, M, MODE, N
      DOUBLE PRECISION D1MACH, EPS, EPSF, EPSLOG, EPSMCH, ERR(*), &
           FACTOR, FJAC(LDFJAC,*), FVEC(*), FVECP(*), ONE, TEMP, X(*), &
           XP(*), ZERO
      SAVE FACTOR, ONE, ZERO
      DATA FACTOR,ONE,ZERO /1.0D2,1.0D0,0.0D0/
!C
!C     EPSMCH IS THE MACHINE PRECISION.
!C
!C***FIRST EXECUTABLE STATEMENT  DCKDER
      EPSMCH = D1MACH(4)
!C
      EPS = SQRT(EPSMCH)
!C
      IF (MODE .EQ. 2) GO TO 20
!C
!C        MODE = 1.
!C
         DO 10 J = 1, N
            TEMP = EPS*ABS(X(J))
            IF (TEMP .EQ. ZERO) TEMP = EPS
            XP(J) = X(J) + TEMP
   10       CONTINUE
         GO TO 70
   20 CONTINUE
!C
!C        MODE = 2.
!C
         EPSF = FACTOR*EPSMCH
         EPSLOG = LOG10(EPS)
         DO 30 I = 1, M
            ERR(I) = ZERO
   30       CONTINUE
         DO 50 J = 1, N
            TEMP = ABS(X(J))
            IF (TEMP .EQ. ZERO) TEMP = ONE
            DO 40 I = 1, M
               ERR(I) = ERR(I) + TEMP*FJAC(I,J)
   40          CONTINUE
   50       CONTINUE
         DO 60 I = 1, M
            TEMP = ONE
            IF (FVEC(I) .NE. ZERO .AND. FVECP(I) .NE. ZERO &
                .AND. ABS(FVECP(I)-FVEC(I)) .GE. EPSF*ABS(FVEC(I))) &
               TEMP = EPS*ABS((FVECP(I)-FVEC(I))/EPS-ERR(I)) &
                      /(ABS(FVEC(I)) + ABS(FVECP(I)))
            ERR(I) = ONE
            IF (TEMP .GT. EPSMCH .AND. TEMP .LT. EPS) &
               ERR(I) = (LOG10(TEMP) - EPSLOG)/EPSLOG
            IF (TEMP .GE. EPS) ERR(I) = ZERO
   60       CONTINUE
   70 CONTINUE
!C
      RETURN
!C
!C     LAST CARD OF SUBROUTINE DCKDER.
!C
      END SUBROUTINE DCKDER
!
!----------------------------------------------------------------------
!
!*DECK DENORM
      FUNCTION DENORM (N, X)
!C***BEGIN PROLOGUE  DENORM
!C***SUBSIDIARY
!C***PURPOSE  Subsidiary to DNSQ and DNSQE
!C***LIBRARY   SLATEC
!C***TYPE      DOUBLE PRECISION (ENORM-S, DENORM-D)
!C***AUTHOR  (UNKNOWN)
!C***DESCRIPTION
!C
!C     Given an N-vector X, this function calculates the
!C     Euclidean norm of X.
!C
!C     The Euclidean norm is computed by accumulating the sum of
!C     squares in three different sums. The sums of squares for the
!C     small and large components are scaled so that no overflows
!C     occur. Non-destructive underflows are permitted. Underflows
!C     and overflows do not occur in the computation of the unscaled
!C     sum of squares for the intermediate components.
!C     The definitions of small, intermediate and large components
!C     depend on two constants, RDWARF and RGIANT. The main
!C     restrictions on these constants are that RDWARF**2 not
!C     underflow and RGIANT**2 not overflow. The constants
!C     given here are suitable for every known computer.
!C
!C     The function statement is
!C
!C       DOUBLE PRECISION FUNCTION DENORM(N,X)
!C
!C     where
!C
!C       N is a positive integer input variable.
!C
!C       X is an input array of length N.
!C
!C***SEE ALSO  DNSQ, DNSQE
!C***ROUTINES CALLED  (NONE)
!C***REVISION HISTORY  (YYMMDD)
!C   800301  DATE WRITTEN
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890831  Modified array declarations.  (WRB)
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900326  Removed duplicate information from DESCRIPTION section.
!C           (WRB)
!C   900328  Added TYPE section.  (WRB)
!C***END PROLOGUE  DENORM
      INTEGER I, N
      REAL(DOUBLE) :: DENORM
      DOUBLE PRECISION AGIANT, FLOATN, ONE, RDWARF, RGIANT, S1, S2, S3,&
           X(*), X1MAX, X3MAX, XABS, ZERO
      SAVE ONE, ZERO, RDWARF, RGIANT
      DATA ONE,ZERO,RDWARF,RGIANT /1.0D0,0.0D0,3.834D-20,1.304D19/
!C***FIRST EXECUTABLE STATEMENT  DENORM
      S1 = ZERO
      S2 = ZERO
      S3 = ZERO
      X1MAX = ZERO
      X3MAX = ZERO
      FLOATN = N
      AGIANT = RGIANT/FLOATN
      DO 90 I = 1, N
         XABS = ABS(X(I))
         IF (XABS .GT. RDWARF .AND. XABS .LT. AGIANT) GO TO 70
            IF (XABS .LE. RDWARF) GO TO 30
!C
!C              SUM FOR LARGE COMPONENTS.
!C
               IF (XABS .LE. X1MAX) GO TO 10
                  S1 = ONE + S1*(X1MAX/XABS)**2
                  X1MAX = XABS
                  GO TO 20
   10          CONTINUE
                  S1 = S1 + (XABS/X1MAX)**2
   20          CONTINUE
               GO TO 60
   30       CONTINUE
!C
!C              SUM FOR SMALL COMPONENTS.
!C
               IF (XABS .LE. X3MAX) GO TO 40
                  S3 = ONE + S3*(X3MAX/XABS)**2
                  X3MAX = XABS
                  GO TO 50
   40          CONTINUE
                  IF (XABS .NE. ZERO) S3 = S3 + (XABS/X3MAX)**2
   50          CONTINUE
   60       CONTINUE
            GO TO 80
   70    CONTINUE
!C
!C           SUM FOR INTERMEDIATE COMPONENTS.
!C
            S2 = S2 + XABS**2
   80    CONTINUE
   90    CONTINUE
!C
!C     CALCULATION OF NORM.
!C
      IF (S1 .EQ. ZERO) GO TO 100
         DENORM = X1MAX*SQRT(S1+(S2/X1MAX)/X1MAX)
         GO TO 130
  100 CONTINUE
         IF (S2 .EQ. ZERO) GO TO 110
            IF (S2 .GE. X3MAX) &
               DENORM = SQRT(S2*(ONE+(X3MAX/S2)*(X3MAX*S3)))
            IF (S2 .LT. X3MAX) &
               DENORM = SQRT(X3MAX*((S2/X3MAX)+(X3MAX*S3)))
            GO TO 120
  110    CONTINUE
            DENORM = X3MAX*SQRT(S3)
  120    CONTINUE
  130 CONTINUE
      RETURN
!C
!C     LAST CARD OF FUNCTION DENORM.
!C
      END FUNCTION DENORM
!*DECK DFDJC3
      SUBROUTINE DFDJC3 (FCN, M, N, X, FVEC, FJAC, LDFJAC, IFLAG, EPSFCN, WA)
!C***BEGIN PROLOGUE  DFDJC3
!C***SUBSIDIARY
!C***PURPOSE  Subsidiary to DNLS1 and DNLS1E
!C***LIBRARY   SLATEC
!C***TYPE      DOUBLE PRECISION (FDJAC3-S, DFDJC3-D)
!C***AUTHOR  (UNKNOWN)
!C***DESCRIPTION
!C
!C  **** Double Precision version of FDJAC3 ****
!C
!C     This subroutine computes a forward-difference approximation
!C     to the M by N Jacobian matrix associated with a specified
!C     problem of M functions in N variables.
!C
!C     The subroutine statement is
!C
!C       SUBROUTINE DFDJC3(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA)
!C
!C     where
!C
!C       FCN is the name of the user-supplied subroutine which
!C         calculates the functions. FCN must be declared
!C         in an external statement in the user calling
!C         program, and should be written as follows.
!C
!C         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
!C         INTEGER LDFJAC,M,N,IFLAG
!C         DOUBLE PRECISION X(N),FVEC(M),FJAC(LDFJAC,N)
!C         ----------
!C         When IFLAG.EQ.1 calculate the functions at X and
!C         return this vector in FVEC.
!C         ----------
!C         RETURN
!C         END
!C
!C         The value of IFLAG should not be changed by FCN unless
!C         the user wants to terminate execution of DFDJC3.
!C         In this case set IFLAG to a negative integer.
!C
!C       M is a positive integer input variable set to the number
!C         of functions.
!C
!C       N is a positive integer input variable set to the number
!C         of variables. N must not exceed M.
!C
!C       X is an input array of length N.
!C
!C       FVEC is an input array of length M which must contain the
!C         functions evaluated at X.
!C
!C       FJAC is an output M by N array which contains the
!C         approximation to the Jacobian matrix evaluated at X.
!C
!C       LDFJAC is a positive integer input variable not less than M
!C         which specifies the leading dimension of the array FJAC.
!C
!C       IFLAG is an integer variable which can be used to terminate
!C         THE EXECUTION OF DFDJC3. See description of FCN.
!C
!C       EPSFCN is an input variable used in determining a suitable
!C         step length for the forward-difference approximation. This
!C         approximation assumes that the relative errors in the
!C         functions are of the order of EPSFCN. If EPSFCN is less
!C         than the machine precision, it is assumed that the relative
!C         errors in the functions are of the order of the machine
!C         precision.
!C
!C       WA is a work array of length M.
!C
!C***SEE ALSO  DNLS1, DNLS1E
!C***ROUTINES CALLED  D1MACH
!C***REVISION HISTORY  (YYMMDD)
!C   800301  DATE WRITTEN
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890831  Modified array declarations.  (WRB)
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900326  Removed duplicate information from DESCRIPTION section.
!C           (WRB)
!C   900328  Added TYPE section.  (WRB)
!C***END PROLOGUE  DFDJC3
      INTEGER M,N,LDFJAC,IFLAG
      DOUBLE PRECISION EPSFCN
      DOUBLE PRECISION X(*),FVEC(*),FJAC(LDFJAC,*),WA(*)
      INTEGER I,J
      DOUBLE PRECISION EPS,EPSMCH,H,TEMP,ZERO
      DOUBLE PRECISION D1MACH
      EXTERNAL FCN
      SAVE ZERO
      DATA ZERO /0.0D0/
!C***FIRST EXECUTABLE STATEMENT  DFDJC3
      EPSMCH = D1MACH(4)
!C
      EPS = SQRT(MAX(EPSFCN,EPSMCH))
!C      SET IFLAG=1 TO INDICATE THAT FUNCTION VALUES
!C           ARE TO BE RETURNED BY FCN.
      IFLAG = 1
      DO 20 J = 1, N
         TEMP = X(J)
         H = EPS*ABS(TEMP)
         IF (H .EQ. ZERO) H = EPS
         X(J) = TEMP + H
         CALL FCN(IFLAG,M,N,X,WA,FJAC,LDFJAC)
         IF (IFLAG .LT. 0) GO TO 30
         X(J) = TEMP
         DO 10 I = 1, M
            FJAC(I,J) = (WA(I) - FVEC(I))/H
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
      RETURN
!C
!C     LAST CARD OF SUBROUTINE DFDJC3.
!C
      END SUBROUTINE DFDJC3
!*DECK DMPAR
      SUBROUTINE DMPAR (N, R, LDR, IPVT, DIAG, QTB, DELTA, PAR, X, SIGMA, WA1, WA2)
!C***BEGIN PROLOGUE  DMPAR
!C***SUBSIDIARY
!C***PURPOSE  Subsidiary to DNLS1 and DNLS1E
!C***LIBRARY   SLATEC
!C***TYPE      DOUBLE PRECISION (LMPAR-S, DMPAR-D)
!C***AUTHOR  (UNKNOWN)
!C***DESCRIPTION
!C
!C   **** Double Precision version of LMPAR ****
!C
!C     Given an M by N matrix A, an N by N nonsingular DIAGONAL
!C     matrix D, an M-vector B, and a positive number DELTA,
!C     the problem is to determine a value for the parameter
!C     PAR such that if X solves the system
!C
!C           A*X = B ,     SQRT(PAR)*D*X = 0 ,
!C
!C     in the least squares sense, and DXNORM is the Euclidean
!C     norm of D*X, then either PAR is zero and
!C
!C           (DXNORM-DELTA) .LE. 0.1*DELTA ,
!C
!C     or PAR is positive and
!C
!C           ABS(DXNORM-DELTA) .LE. 0.1*DELTA .
!C
!C     This subroutine completes the solution of the problem
!C     if it is provided with the necessary information from the
!C     QR factorization, with column pivoting, of A. That is, if
!C     A*P = Q*R, where P is a permutation matrix, Q has orthogonal
!C     columns, and R is an upper triangular matrix with diagonal
!C     elements of nonincreasing magnitude, then DMPAR expects
!C     the full upper triangle of R, the permutation matrix P,
!C     and the first N components of (Q TRANSPOSE)*B. On output
!C     DMPAR also provides an upper triangular matrix S such that
!C
!C            T   T                   T
!C           P *(A *A + PAR*D*D)*P = S *S .
!C
!C     S is employed within DMPAR and may be of separate interest.
!C
!C     Only a few iterations are generally needed for convergence
!C     of the algorithm. If, however, the limit of 10 iterations
!C     is reached, then the output PAR will contain the best
!C     value obtained so far.
!C
!C     The subroutine statement is
!C
!C       SUBROUTINE DMPAR(N,R,LDR,IPVT,DIAG,QTB,DELTA,PAR,X,SIGMA,
!C                        WA1,WA2)
!C
!C     where
!C
!C       N is a positive integer input variable set to the order of R.
!C
!C       R is an N by N array. On input the full upper triangle
!C         must contain the full upper triangle of the matrix R.
!C         On output the full upper triangle is unaltered, and the
!C         strict lower triangle contains the strict upper triangle
!C         (transposed) of the upper triangular matrix S.
!C
!C       LDR is a positive integer input variable not less than N
!C         which specifies the leading dimension of the array R.
!C
!C       IPVT is an integer input array of length N which defines the
!C         permutation matrix P such that A*P = Q*R. Column J of P
!C         is column IPVT(J) of the identity matrix.
!C
!C       DIAG is an input array of length N which must contain the
!C         diagonal elements of the matrix D.
!C
!C       QTB is an input array of length N which must contain the first
!C         N elements of the vector (Q TRANSPOSE)*B.
!C
!C       DELTA is a positive input variable which specifies an upper
!C         bound on the Euclidean norm of D*X.
!C
!C       PAR is a nonnegative variable. On input PAR contains an
!C         initial estimate of the Levenberg-Marquardt parameter.
!C         On output PAR contains the final estimate.
!C
!C       X is an output array of length N which contains the least
!C         squares solution of the system A*X = B, SQRT(PAR)*D*X = 0,
!C         for the output PAR.
!C
!C       SIGMA is an output array of length N which contains the
!C         diagonal elements of the upper triangular matrix S.
!C
!C       WA1 and WA2 are work arrays of length N.
!C
!C***SEE ALSO  DNLS1, DNLS1E
!C***ROUTINES CALLED  D1MACH, DENORM, DQRSLV
!C***REVISION HISTORY  (YYMMDD)
!C   800301  DATE WRITTEN
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890831  Modified array declarations.  (WRB)
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900326  Removed duplicate information from DESCRIPTION section.
!C           (WRB)
!C   900328  Added TYPE section.  (WRB)
!C***END PROLOGUE  DMPAR
      INTEGER N,LDR
      INTEGER IPVT(*)
      DOUBLE PRECISION DELTA,PAR
      DOUBLE PRECISION R(LDR,*),DIAG(*),QTB(*),X(*),SIGMA(*),WA1(*), &
       WA2(*)
      INTEGER I,ITER,J,JM1,JP1,K,L,NSING
      DOUBLE PRECISION DXNORM,DWARF,FP,GNORM,PARC,PARL,PARU,P1,P001, &
       SUM,TEMP,ZERO
      DOUBLE PRECISION D1MACH
      SAVE P1, P001, ZERO
      DATA P1,P001,ZERO /1.0D-1,1.0D-3,0.0D0/
!C***FIRST EXECUTABLE STATEMENT  DMPAR
      DWARF = D1MACH(1)
!C
!C     COMPUTE AND STORE IN X THE GAUSS-NEWTON DIRECTION. IF THE
!C     JACOBIAN IS RANK-DEFICIENT, OBTAIN A LEAST SQUARES SOLUTION.
!C
      NSING = N
      DO 10 J = 1, N
         WA1(J) = QTB(J)
         IF (R(J,J) .EQ. ZERO .AND. NSING .EQ. N) NSING = J - 1
         IF (NSING .LT. N) WA1(J) = ZERO
   10    CONTINUE
      IF (NSING .LT. 1) GO TO 50
      DO 40 K = 1, NSING
         J = NSING - K + 1
         WA1(J) = WA1(J)/R(J,J)
         TEMP = WA1(J)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 30
         DO 20 I = 1, JM1
            WA1(I) = WA1(I) - R(I,J)*TEMP
   20       CONTINUE
   30    CONTINUE
   40    CONTINUE
   50 CONTINUE
      DO 60 J = 1, N
         L = IPVT(J)
         X(L) = WA1(J)
   60    CONTINUE
!C
!C     INITIALIZE THE ITERATION COUNTER.
!C     EVALUATE THE FUNCTION AT THE ORIGIN, AND TEST
!C     FOR ACCEPTANCE OF THE GAUSS-NEWTON DIRECTION.
!C
      ITER = 0
      DO 70 J = 1, N
         WA2(J) = DIAG(J)*X(J)
   70    CONTINUE
      DXNORM = DENORM(N,WA2)
      FP = DXNORM - DELTA
      IF (FP .LE. P1*DELTA) GO TO 220
!C
!C     IF THE JACOBIAN IS NOT RANK DEFICIENT, THE NEWTON
!C     STEP PROVIDES A LOWER BOUND, PARL, FOR THE ZERO OF
!C     THE FUNCTION. OTHERWISE SET THIS BOUND TO ZERO.
!C
      PARL = ZERO
      IF (NSING .LT. N) GO TO 120
      DO 80 J = 1, N
         L = IPVT(J)
         WA1(J) = DIAG(L)*(WA2(L)/DXNORM)
   80    CONTINUE
      DO 110 J = 1, N
         SUM = ZERO
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 100
         DO 90 I = 1, JM1
            SUM = SUM + R(I,J)*WA1(I)
   90       CONTINUE
  100    CONTINUE
         WA1(J) = (WA1(J) - SUM)/R(J,J)
  110    CONTINUE
      TEMP = DENORM(N,WA1)
      PARL = ((FP/DELTA)/TEMP)/TEMP
  120 CONTINUE
!C
!C     CALCULATE AN UPPER BOUND, PARU, FOR THE ZERO OF THE FUNCTION.
!C
      DO 140 J = 1, N
         SUM = ZERO
         DO 130 I = 1, J
            SUM = SUM + R(I,J)*QTB(I)
  130       CONTINUE
         L = IPVT(J)
         WA1(J) = SUM/DIAG(L)
  140    CONTINUE
      GNORM = DENORM(N,WA1)
      PARU = GNORM/DELTA
      IF (PARU .EQ. ZERO) PARU = DWARF/MIN(DELTA,P1)
!C
!C     IF THE INPUT PAR LIES OUTSIDE OF THE INTERVAL (PARL,PARU),
!C     SET PAR TO THE CLOSER ENDPOINT.
!C
      PAR = MAX(PAR,PARL)
      PAR = MIN(PAR,PARU)
      IF (PAR .EQ. ZERO) PAR = GNORM/DXNORM
!C
!C     BEGINNING OF AN ITERATION.
!C
  150 CONTINUE
         ITER = ITER + 1
!C
!C        EVALUATE THE FUNCTION AT THE CURRENT VALUE OF PAR.
!C
         IF (PAR .EQ. ZERO) PAR = MAX(DWARF,P001*PARU)
         TEMP = SQRT(PAR)
         DO 160 J = 1, N
            WA1(J) = TEMP*DIAG(J)
  160       CONTINUE
         CALL DQRSLV(N,R,LDR,IPVT,WA1,QTB,X,SIGMA,WA2)
         DO 170 J = 1, N
            WA2(J) = DIAG(J)*X(J)
  170       CONTINUE
         DXNORM = DENORM(N,WA2)
         TEMP = FP
         FP = DXNORM - DELTA
!C
!C        IF THE FUNCTION IS SMALL ENOUGH, ACCEPT THE CURRENT VALUE
!C        OF PAR. ALSO TEST FOR THE EXCEPTIONAL CASES WHERE PARL
!C        IS ZERO OR THE NUMBER OF ITERATIONS HAS REACHED 10.
!C
         IF (ABS(FP) .LE. P1*DELTA &
             .OR. PARL .EQ. ZERO .AND. FP .LE. TEMP &
                  .AND. TEMP .LT. ZERO .OR. ITER .EQ. 10) GO TO 220
!C
!C        COMPUTE THE NEWTON CORRECTION.
!C
         DO 180 J = 1, N
            L = IPVT(J)
            WA1(J) = DIAG(L)*(WA2(L)/DXNORM)
  180       CONTINUE
         DO 210 J = 1, N
            WA1(J) = WA1(J)/SIGMA(J)
            TEMP = WA1(J)
            JP1 = J + 1
            IF (N .LT. JP1) GO TO 200
            DO 190 I = JP1, N
               WA1(I) = WA1(I) - R(I,J)*TEMP
  190          CONTINUE
  200       CONTINUE
  210       CONTINUE
         TEMP = DENORM(N,WA1)
         PARC = ((FP/DELTA)/TEMP)/TEMP
!C
!C        DEPENDING ON THE SIGN OF THE FUNCTION, UPDATE PARL OR PARU.
!C
         IF (FP .GT. ZERO) PARL = MAX(PARL,PAR)
         IF (FP .LT. ZERO) PARU = MIN(PARU,PAR)
!C
!C        COMPUTE AN IMPROVED ESTIMATE FOR PAR.
!C
         PAR = MAX(PARL,PAR+PARC)
!C
!C        END OF AN ITERATION.
!C
         GO TO 150
  220 CONTINUE
!C
!C     TERMINATION.
!C
      IF (ITER .EQ. 0) PAR = ZERO
      RETURN
!C
!C     LAST CARD OF SUBROUTINE DMPAR.
!C
      END SUBROUTINE DMPAR
!*DECK DNLS1E
      SUBROUTINE DNLS1E (FCN, IOPT, M, N, X, FVEC, TOL,  &
                         NPRINT, INFO, IW, WA, LWA, ERR)
!C***BEGIN PROLOGUE  DNLS1E
!C***PURPOSE  An easy-to-use code which minimizes the sum of the squares
!C            of M nonlinear functions in N variables by a modification
!C            of the Levenberg-Marquardt algorithm.
!C***LIBRARY   SLATEC
!C***CATEGORY  K1B1A1, K1B1A2
!C***TYPE      DOUBLE PRECISION (SNLS1E-S, DNLS1E-D)
!C***KEYWORDS  EASY-TO-USE, LEVENBERG-MARQUARDT, NONLINEAR DATA FITTING,
!C             NONLINEAR LEAST SQUARES
!C***AUTHOR  Hiebert, K. L., (SNLA)
!C***DESCRIPTION
!C
!C 1. Purpose.
!C
!C       The purpose of DNLS1E is to minimize the sum of the squares of M
!C       nonlinear functions in N variables by a modification of the
!C       Levenberg-Marquardt algorithm.  This is done by using the more
!C       general least-squares solver DNLS1.  The user must provide a
!C       subroutine which calculates the functions.  The user has the
!C       option of how the Jacobian will be supplied.  The user can
!C       supply the full Jacobian, or the rows of the Jacobian (to avoid
!C       storing the full Jacobian), or let the code approximate the
!C       Jacobian by forward-differencing.  This code is the combination
!C       of the MINPACK codes (Argonne) LMDER1, LMDIF1, and LMSTR1.
!C
!C
!C 2. Subroutine and Type Statements.
!C
!C       SUBROUTINE DNLS1E(FCN,IOPT,M,N,X,FVEC,TOL,NPRINT,
!C      *                  INFO,IW,WA,LWA)
!C       INTEGER IOPT,M,N,NPRINT,INFO,LWAC,IW(N)
!C       DOUBLE PRECISION TOL,X(N),FVEC(M),WA(LWA)
!C       EXTERNAL FCN
!C
!C
!C  C 3. Parameters. ALL TYPE REAL parameters are DOUBLE PRECISION
!C
!C       Parameters designated as input parameters must be specified on
!C       entry to DNLS1E and are not changed on exit, while parameters
!C       designated as output parameters need not be specified on entry
!C       and are set to appropriate values on exit from DNLS1E.
!C
!C      FCN is the name of the user-supplied subroutine which calculates
!C         the functions.  If the user wants to supply the Jacobian
!C         (IOPT=2 or 3), then FCN must be written to calculate the
!C         Jacobian, as well as the functions.  See the explanation
!C         of the IOPT argument below.
!C         If the user wants the iterates printed (NPRINT positive), then
!C         FCN must do the printing.  See the explanation of NPRINT
!C         below.  FCN must be declared in an EXTERNAL statement in the
!C         calling program and should be written as follows.
!C
!C
!C         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
!C         INTEGER IFLAG,LDFJAC,M,N
!C         DOUBLE PRECISION X(N),FVEC(M)
!C         ----------
!C         FJAC and LDFJAC may be ignored       , if IOPT=1.
!C         DOUBLE PRECISION FJAC(LDFJAC,N)      , if IOPT=2.
!C         DOUBLE PRECISION FJAC(N)             , if IOPT=3.
!C         ----------
!C           If IFLAG=0, the values in X and FVEC are available
!C           for printing.  See the explanation of NPRINT below.
!C           IFLAG will never be zero unless NPRINT is positive.
!C           The values of X and FVEC must not be changed.
!C         RETURN
!C         ----------
!C           If IFLAG=1, calculate the functions at X and return
!C           this vector in FVEC.
!C         RETURN
!C         ----------
!C           If IFLAG=2, calculate the full Jacobian at X and return
!C           this matrix in FJAC.  Note that IFLAG will never be 2 unless
!C           IOPT=2.  FVEC contains the function values at X and must
!C           not be altered.  FJAC(I,J) must be set to the derivative
!C           of FVEC(I) with respect to X(J).
!C         RETURN
!C         ----------
!C           If IFLAG=3, calculate the LDFJAC-th row of the Jacobian
!C           and return this vector in FJAC.  Note that IFLAG will
!C           never be 3 unless IOPT=3.  FVEC contains the function
!C           values at X and must not be altered.  FJAC(J) must be
!C           set to the derivative of FVEC(LDFJAC) with respect to X(J).
!C         RETURN
!C         ----------
!C         END
!C
!C
!C         The value of IFLAG should not be changed by FCN unless the
!C         user wants to terminate execution of DNLS1E.  In this case,
!C         set IFLAG to a negative integer.
!C
!C
!C       IOPT is an input variable which specifies how the Jacobian will
!C         be calculated.  If IOPT=2 or 3, then the user must supply the
!C         Jacobian, as well as the function values, through the
!C         subroutine FCN.  If IOPT=2, the user supplies the full
!C         Jacobian with one call to FCN.  If IOPT=3, the user supplies
!C         one row of the Jacobian with each call.  (In this manner,
!C         storage can be saved because the full Jacobian is not stored.)
!C         If IOPT=1, the code will approximate the Jacobian by forward
!C         differencing.
!C
!C       M is a positive integer input variable set to the number of
!C         functions.
!C
!C       N is a positive integer input variable set to the number of
!C         variables.  N must not exceed M.
!C
!C       X is an array of length N.  On input, X must contain an initial
!C         estimate of the solution vector.  On output, X contains the
!C         final estimate of the solution vector.
!C
!C       FVEC is an output array of length M which contains the functions
!C         evaluated at the output X.
!C
!C       TOL is a non-negative input variable.  Termination occurs when
!C         the algorithm estimates either that the relative error in the
!C         sum of squares is at most TOL or that the relative error
!C         between X and the solution is at most TOL.  Section 4 contains
!C         more details about TOL.
!C
!C       NPRINT is an integer input variable that enables controlled
!C         printing of iterates if it is positive.  In this case, FCN is
!C         called with IFLAG = 0 at the beginning of the first iteration
!C         and every NPRINT iterations thereafter and immediately prior
!C         to return, with X and FVEC available for printing. Appropriate
!C         print statements must be added to FCN (see example) and
!C         FVEC should not be altered.  If NPRINT is not positive, no
!C         special calls of FCN with IFLAG = 0 are made.
!C
!C       INFO is an integer output variable.  If the user has terminated
!C        execution, INFO is set to the (negative) value of IFLAG.  See
!C        description of FCN and JAC. Otherwise, INFO is set as follows.
!C
!C         INFO = 0  improper input parameters.
!C
!C         INFO = 1  algorithm estimates that the relative error in the
!C                   sum of squares is at most TOL.
!C
!C         INFO = 2  algorithm estimates that the relative error between
!C                   X and the solution is at most TOL.
!C
!C         INFO = 3  conditions for INFO = 1 and INFO = 2 both hold.
!C
!C         INFO = 4  FVEC is orthogonal to the columns of the Jacobian to
!C                   machine precision.
!C
!C         INFO = 5  number of calls to FCN has reached 100*(N+1)
!C                   for IOPT=2 or 3 or 200*(N+1) for IOPT=1.
!C
!C         INFO = 6  TOL is too small.  No further reduction in the sum
!C                   of squares is possible.
!C
!C         INFO = 7  TOL is too small.  No further improvement in the
!C                   approximate solution X is possible.
!C
!C         Sections 4 and 5 contain more details about INFO.
!C
!C       IW is an INTEGER work array of length N.
!C
!C       WA is a work array of length LWA.
!C
!C       LWA is a positive integer input variable not less than
!C         N*(M+5)+M for IOPT=1 and 2 or N*(N+5)+M for IOPT=3.
!C
!C
!C  C 4. Successful Completion.
!C
!C       The accuracy of DNLS1E is controlled by the convergence parame-
!C       ter TOL.  This parameter is used in tests which make three types
!C       of comparisons between the approximation X and a solution XSOL.
!C       DNLS1E terminates when any of the tests is satisfied.  If TOL is
!C       less than the machine precision (as defined by the function
!C       R1MACH(4)), then DNLS1E only attempts to satisfy the test
!C       defined by the machine precision.  Further progress is not usu-
!C       ally possible.  Unless high precision solutions are required,
!C       the recommended value for TOL is the square root of the machine
!C       precision.
!C
!C       The tests assume that the functions are reasonably well behaved,
!C       and, if the Jacobian is supplied by the user, that the functions
!C       and the Jacobian are coded consistently.  If these conditions
!C       are not satisfied, then DNLS1E may incorrectly indicate conver-
!C       gence.  If the Jacobian is coded correctly or IOPT=1,
!C       then the validity of the answer can be checked, for example, by
!C       rerunning DNLS1E with tighter tolerances.
!C
!C       First Convergence Test.  If ENORM(Z) denotes the Euclidean norm
!C         of a vector Z, then this test attempts to guarantee that
!C
!C               ENORM(FVEC) .LE. (1+TOL)*ENORM(FVECS),
!C
!C         where FVECS denotes the functions evaluated at XSOL.  If this
!C         condition is satisfied with TOL = 10**(-K), then the final
!C         residual norm ENORM(FVEC) has K significant decimal digits and
!C         INFO is set to 1 (or to 3 if the second test is also satis-
!C         fied).
!C
!C       Second Convergence Test.  If D is a diagonal matrix (implicitly
!C         generated by DNLS1E) whose entries contain scale factors for
!C         the variables, then this test attempts to guarantee that
!C
!C               ENORM(D*(X-XSOL)) .LE.  TOL*ENORM(D*XSOL).
!C
!C         If this condition is satisfied with TOL = 10**(-K), then the
!C         larger components of D*X have K significant decimal digits and
!C         INFO is set to 2 (or to 3 if the first test is also satis-
!C         fied).  There is a danger that the smaller components of D*X
!C         may have large relative errors, but the choice of D is such
!C         that the accuracy of the components of X is usually related to
!C         their sensitivity.
!C
!C       Third Convergence Test.  This test is satisfied when FVEC is
!C         orthogonal to the columns of the Jacobian to machine preci-
!C         sion.  There is no clear relationship between this test and
!C         the accuracy of DNLS1E, and furthermore, the test is equally
!C         well satisfied at other critical points, namely maximizers and
!C         saddle points.  Therefore, termination caused by this test
!C         (INFO = 4) should be examined carefully.
!C
!C
!C  C 5. Unsuccessful Completion.
!C
!C       Unsuccessful termination of DNLS1E can be due to improper input
!C       parameters, arithmetic interrupts, or an excessive number of
!C       function evaluations.
!C
!C       Improper Input Parameters.  INFO is set to 0 if IOPT .LT. 1
!C         or IOPT .GT. 3, or N .LE. 0, or M .LT. N, or TOL .LT. 0.E0,
!C         or for IOPT=1 or 2 LWA .LT. N*(M+5)+M, or for IOPT=3
!C         LWA .LT. N*(N+5)+M.
!C
!C       Arithmetic Interrupts.  If these interrupts occur in the FCN
!C         subroutine during an early stage of the computation, they may
!C         be caused by an unacceptable choice of X by DNLS1E.  In this
!C         case, it may be possible to remedy the situation by not evalu-
!C         ating the functions here, but instead setting the components
!C         of FVEC to numbers that exceed those in the initial FVEC.
!C
!C       Excessive Number of Function Evaluations.  If the number of
!C         calls to FCN reaches 100*(N+1) for IOPT=2 or 3 or 200*(N+1)
!C         for IOPT=1, then this indicates that the routine is converging
!C         very slowly as measured by the progress of FVEC, and INFO is
!C         set to 5.  In this case, it may be helpful to restart DNLS1E,
!C         thereby forcing it to disregard old (and possibly harmful)
!C         information.
!C
!C
!C  C 6. Characteristics of the Algorithm.
!C
!C       DNLS1E is a modification of the Levenberg-Marquardt algorithm.
!C       Two of its main characteristics involve the proper use of
!C       implicitly scaled variables and an optimal choice for the cor-
!C       rection.  The use of implicitly scaled variables achieves scale
!C       invariance of DNLS1E and limits the size of the correction in
!C       any direction where the functions are changing rapidly.  The
!C       optimal choice of the correction guarantees (under reasonable
!C       conditions) global convergence from starting points far from the
!C       solution and a fast rate of convergence for problems with small
!C       residuals.
!C
!C       Timing.  The time required by DNLS1E to solve a given problem
!C         depends on M and N, the behavior of the functions, the accu-
!C         racy requested, and the starting point.  The number of arith-
!C         metic operations needed by DNLS1E is about N**3 to process
!C         each evaluation of the functions (call to FCN) and to process
!C         each evaluation of the Jacobian DNLS1E takes M*N**2 for IOPT=2
!C         (one call to JAC), M*N**2 for IOPT=1 (N calls to FCN) and
!C         1.5*M*N**2 for IOPT=3 (M calls to FCN).  Unless FCN
!C         can be evaluated quickly, the timing of DNLS1E will be
!C         strongly influenced by the time spent in FCN.
!C
!C       Storage.  DNLS1E requires (M*N + 2*M + 6*N) for IOPT=1 or 2 and
!C         (N**2 + 2*M + 6*N) for IOPT=3 single precision storage
!C         locations and N integer storage locations, in addition to
!C         the storage required by the program.  There are no internally
!C         declared storage arrays.
!C
!C  C *Long Description:
!C
!C  C 7. Example.
!C
!C       The problem is to determine the values of X(1), X(2), and X(3)
!C       which provide the best fit (in the least squares sense) of
!C
!C             X(1) + U(I)/(V(I)*X(2) + W(I)*X(3)),  I = 1, 15
!C
!C       to the data
!C
!C             Y = (0.14,0.18,0.22,0.25,0.29,0.32,0.35,0.39,
!C                  0.37,0.58,0.73,0.96,1.34,2.10,4.39),
!C
!C       where U(I) = I, V(I) = 16 - I, and W(I) = MIN(U(I),V(I)).  The
!C       I-th component of FVEC is thus defined by
!C
!C             Y(I) - (X(1) + U(I)/(V(I)*X(2) + W(I)*X(3))).
!C
!C       **********
!C
!C       PROGRAM TEST
!C  C C
!C  C !C     Driver for DNLS1E example.
!C  C C
!C       INTEGER I,IOPT,M,N,NPRINT,JNFO,LWA,NWRITE
!C       INTEGER IW(3)
!C       DOUBLE PRECISION TOL,FNORM,X(3),FVEC(15),WA(75)
!C       DOUBLE PRECISION DENORM,D1MACH
!C       EXTERNAL FCN
!C       DATA NWRITE /6/
!C  C C
!C       IOPT = 1
!C       M = 15
!C       N = 3
!C  C C
!C  C !C     The following starting values provide a rough fit.
!C  C C
!C       X(1) = 1.E0
!C       X(2) = 1.E0
!C       X(3) = 1.E0
!C  C C
!C       LWA = 75
!C       NPRINT = 0
!C  C C
!C  C !C     Set TOL to the square root of the machine precision.
!C  C !C     Unless high precision solutions are required,
!C  C !C     this is the recommended setting.
!C  C C
!C       TOL = SQRT(R1MACH(4))
!C  C C
!C       CALL DNLS1E(FCN,IOPT,M,N,X,FVEC,TOL,NPRINT,
!C      *            INFO,IW,WA,LWA)
!C       FNORM = ENORM(M,FVEC)
!C       WRITE (NWRITE,1000) FNORM,INFO,(X(J),J=1,N)
!C       STOP
!C  C  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
!C      *        5X,' EXIT
!C      *        5X,' FINAL APPROXIMATE SOLUTION' // 5X,3E15.7)
!C       END
!C       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,DUM,IDUM)
!C  C !C     This is the form of the FCN routine if IOPT=1,
!C  C !C     that is, if the user does not calculate the Jacobian.
!C       INTEGER I,M,N,IFLAG
!C       DOUBLE PRECISION X(N),FVEC(M),Y(15)
!C       DOUBLE PRECISION TMP1,TMP2,TMP3,TMP4
!C       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
!C      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
!C      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
!C      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
!C  C C
!C       IF (IFLAG .NE. 0) GO TO 5
!C  C C
!C  C !C     Insert print statements here when NPRINT is positive.
!C  C C
!C       RETURN
!C     5 CONTINUE
!C       DO 10 I = 1, M
!C          TMP1 = I
!C          TMP2 = 16 - I
!C          TMP3 = TMP1
!C          IF (I .GT. 8) TMP3 = TMP2
!C          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
!C    10    CONTINUE
!C       RETURN
!C       END
!C
!C
!C       Results obtained with different compilers or machines
!C       may be slightly different.
!C
!C       FINAL L2 NORM OF THE RESIDUALS  0.9063596E-01
!C
!C       EXIT PARAMETER                         1
!C
!C       FINAL APPROXIMATE SOLUTION
!C
!C        0.8241058E-01  0.1133037E+01  0.2343695E+01
!C
!C
!C       For IOPT=2, FCN would be modified as follows to also
!C       calculate the full Jacobian when IFLAG=2.
!C
!C       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
!C  C C
!C  C !C     This is the form of the FCN routine if IOPT=2,
!C  C !C     that is, if the user calculates the full Jacobian.
!C  C C
!C       INTEGER I,LDFJAC,M,N,IFLAG
!C       DOUBLE PRECISION X(N),FVEC(M),FJAC(LDFJAC,N),Y(15)
!C       DOUBLE PRECISION TMP1,TMP2,TMP3,TMP4
!C       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
!C      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
!C      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
!C      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
!C  C C
!C       IF (IFLAG .NE. 0) GO TO 5
!C  C C
!C  C !C     Insert print statements here when NPRINT is positive.
!C  C C
!C       RETURN
!C     5 CONTINUE
!C       IF(IFLAG.NE.1) GO TO 20
!C       DO 10 I = 1, M
!C          TMP1 = I
!C          TMP2 = 16 - I
!C          TMP3 = TMP1
!C          IF (I .GT. 8) TMP3 = TMP2
!C          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
!C    10    CONTINUE
!C       RETURN
!C  C C
!C  C !C     Below, calculate the full Jacobian.
!C  C C
!C    20    CONTINUE
!C  C C
!C       DO 30 I = 1, M
!C          TMP1 = I
!C          TMP2 = 16 - I
!C          TMP3 = TMP1
!C          IF (I .GT. 8) TMP3 = TMP2
!C          TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
!C          FJAC(I,1) = -1.E0
!C          FJAC(I,2) = TMP1*TMP2/TMP4
!C          FJAC(I,3) = TMP1*TMP3/TMP4
!C    30    CONTINUE
!C       RETURN
!C       END
!C
!C
!C       For IOPT = 3, FJAC would be dimensioned as FJAC(3,3),
!C         LDFJAC would be set to 3, and FCN would be written as
!C         follows to calculate a row of the Jacobian when IFLAG=3.
!C
!C       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
!C  C !C     This is the form of the FCN routine if IOPT=3,
!C  C !C     that is, if the user calculates the Jacobian row by row.
!C       INTEGER I,M,N,IFLAG
!C       DOUBLE PRECISION X(N),FVEC(M),FJAC(N),Y(15)
!C       DOUBLE PRECISION TMP1,TMP2,TMP3,TMP4
!C       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
!C      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
!C      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
!C      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
!C  C C
!C       IF (IFLAG .NE. 0) GO TO 5
!C  C C
!C  C !C     Insert print statements here when NPRINT is positive.
!C  C C
!C       RETURN
!C     5 CONTINUE
!C       IF( IFLAG.NE.1) GO TO 20
!C       DO 10 I = 1, M
!C          TMP1 = I
!C          TMP2 = 16 - I
!C          TMP3 = TMP1
!C          IF (I .GT. 8) TMP3 = TMP2
!C          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
!C    10    CONTINUE
!C       RETURN
!C  C C
!C  C !C     Below, calculate the LDFJAC-th row of the Jacobian.
!C  C C
!C    20 CONTINUE
!C
!C       I = LDFJAC
!C          TMP1 = I
!C          TMP2 = 16 - I
!C          TMP3 = TMP1
!C          IF (I .GT. 8) TMP3 = TMP2
!C          TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
!C          FJAC(1) = -1.E0
!C          FJAC(2) = TMP1*TMP2/TMP4
!C          FJAC(3) = TMP1*TMP3/TMP4
!C       RETURN
!C       END
!C
!C***REFERENCES  Jorge J. More, The Levenberg-Marquardt algorithm:
!C                 implementation and theory.  In Numerical Analysis
!C                 Proceedings (Dundee, June 28 - July 1, 1977, G. A.
!C                 Watson, Editor), Lecture Notes in Mathematics 630,
!C                 Springer-Verlag, 1978.
!C***ROUTINES CALLED  DNLS1, XERMSG
!C***REVISION HISTORY  (YYMMDD)
!C   800301  DATE WRITTEN
!C   890831  Modified array declarations.  (WRB)
!C   891006  Cosmetic changes to prologue.  (WRB)
!C   891006  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  DNLS1E
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER M,N,NPRINT,INFO,LWA,IOPT
      INTEGER INDEX,IW(*)
      DOUBLE PRECISION TOL
      DOUBLE PRECISION X(*),FVEC(*),WA(*),ERR(M)
      EXTERNAL FCN
      INTEGER MAXFEV,MODE,NFEV,NJEV
      DOUBLE PRECISION FACTOR,FTOL,GTOL,XTOL,ZERO,EPSFCN
      SAVE FACTOR, ZERO
      DATA FACTOR,ZERO /1.0D2,0.0D0/
!C***FIRST EXECUTABLE STATEMENT  DNLS1E
      INFO = 0
!C
!C     CHECK THE INPUT PARAMETERS FOR ERRORS.
!C
      IF (IOPT .LT. 1 .OR. IOPT .GT. 3 .OR. &
          N .LE. 0 .OR. M .LT. N .OR. TOL .LT. ZERO &
          .OR. LWA .LT. N*(N+5) + M) GO TO 10
      IF (IOPT .LT. 3 .AND. LWA .LT. N*(M+5) + M) GO TO 10
!C
!C     CALL DNLS1.
!C
      MAXFEV = 100*(N + 1)
      IF (IOPT .EQ. 1) MAXFEV = 2*MAXFEV
      FTOL = TOL
      XTOL = TOL
      GTOL = ZERO
      EPSFCN = ZERO
      MODE = 1
      INDEX = 5*N+M
      CALL DNLS1(FCN,IOPT,M,N,X,FVEC,WA(INDEX+1),M,FTOL,XTOL,GTOL, &
                 MAXFEV,EPSFCN,WA(1),MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,&
                 IW,WA(N+1),WA(2*N+1),WA(3*N+1),WA(4*N+1),WA(5*N+1),ERR)
      IF (INFO .EQ. 8) INFO = 4
   10 CONTINUE
      IF (INFO .EQ. 0) CALL XERMSG ('SLATEC', 'DNLS1E', 'INVALID INPUT PARAMETER.', 2, 1)
      RETURN
!C
!C     LAST CARD OF SUBROUTINE DNLS1E.
!C
      END SUBROUTINE DNLS1E
!*DECK DNLS1
      SUBROUTINE DNLS1 (FCN, IOPT, M, N, X, FVEC, FJAC, LDFJAC, FTOL, &
       XTOL, GTOL, MAXFEV, EPSFCN, DIAG, MODE, FACTOR, NPRINT, INFO, &
       NFEV, NJEV, IPVT, QTF, WA1, WA2, WA3, WA4, ERR)
!C***BEGIN PROLOGUE  DNLS1
!C***PURPOSE  Minimize the sum of the squares of M nonlinear functions
!C            in N variables by a modification of the Levenberg-Marquardt
!C            algorithm.
!C***LIBRARY   SLATEC
!C***CATEGORY  K1B1A1, K1B1A2
!C***TYPE      DOUBLE PRECISION (SNLS1-S, DNLS1-D)
!C***KEYWORDS  LEVENBERG-MARQUARDT, NONLINEAR DATA FITTING,
!C             NONLINEAR LEAST SQUARES
!C***AUTHOR  Hiebert, K. L., (SNLA)
!C***DESCRIPTION
!C
!C  C 1. Purpose.
!C
!C       The purpose of DNLS1 is to minimize the sum of the squares of M
!C       nonlinear functions in N variables by a modification of the
!C       Levenberg-Marquardt algorithm.  The user must provide a subrou-
!C       tine which calculates the functions.  The user has the option
!C       of how the Jacobian will be supplied.  The user can supply the
!C       full Jacobian, or the rows of the Jacobian (to avoid storing
!C       the full Jacobian), or let the code approximate the Jacobian by
!C       forward-differencing.   This code is the combination of the
!C       MINPACK codes (Argonne) LMDER, LMDIF, and LMSTR.
!C
!C
!C  C 2. Subroutine and Type Statements.
!C
!C       SUBROUTINE DNLS1(FCN,IOPT,M,N,X,FVEC,FJAC,LDFJAC,FTOL,XTOL,
!C      *                 GTOL,MAXFEV,EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO
!C      *                 ,NFEV,NJEV,IPVT,QTF,WA1,WA2,WA3,WA4)
!C       INTEGER IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV
!C       INTEGER IPVT(N)
!C       DOUBLE PRECISION FTOL,XTOL,GTOL,EPSFCN,FACTOR
!C       DOUBLE PRECISION X(N),FVEC(M),FJAC(LDFJAC,N),DIAG(N),QTF(N),
!C      *     WA1(N),WA2(N),WA3(N),WA4(M)
!C
!C
!C  C 3. Parameters.
!C
!C       Parameters designated as input parameters must be specified on
!C       entry to DNLS1 and are not changed on exit, while parameters
!C       designated as output parameters need not be specified on entry
!C       and are set to appropriate values on exit from DNLS1.
!C
!C      FCN is the name of the user-supplied subroutine which calculate
!C         the functions.  If the user wants to supply the Jacobian
!C         (IOPT=2 or 3), then FCN must be written to calculate the
!C         Jacobian, as well as the functions.  See the explanation
!C         of the IOPT argument below.
!C         If the user wants the iterates printed (NPRINT positive), then
!C         FCN must do the printing.  See the explanation of NPRINT
!C         below.  FCN must be declared in an EXTERNAL statement in the
!C         calling program and should be written as follows.
!C
!C
!C         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
!C         INTEGER IFLAG,LDFJAC,M,N
!C         DOUBLE PRECISION X(N),FVEC(M)
!C         ----------
!C         FJAC and LDFJAC may be ignored       , if IOPT=1.
!C         DOUBLE PRECISION FJAC(LDFJAC,N)      , if IOPT=2.
!C         DOUBLE PRECISION FJAC(N)             , if IOPT=3.
!C         ----------
!C           If IFLAG=0, the values in X and FVEC are available
!C           for printing.  See the explanation of NPRINT below.
!C           IFLAG will never be zero unless NPRINT is positive.
!C           The values of X and FVEC must not be changed.
!C         RETURN
!C         ----------
!C           If IFLAG=1, calculate the functions at X and return
!C           this vector in FVEC.
!C         RETURN
!C         ----------
!C           If IFLAG=2, calculate the full Jacobian at X and return
!C           this matrix in FJAC.  Note that IFLAG will never be 2 unless
!C           IOPT=2.  FVEC contains the function values at X and must
!C           not be altered.  FJAC(I,J) must be set to the derivative
!C           of FVEC(I) with respect to X(J).
!C         RETURN
!C         ----------
!C           If IFLAG=3, calculate the LDFJAC-th row of the Jacobian
!C           and return this vector in FJAC.  Note that IFLAG will
!C           never be 3 unless IOPT=3.  FVEC contains the function
!C           values at X and must not be altered.  FJAC(J) must be
!C           set to the derivative of FVEC(LDFJAC) with respect to X(J).
!C         RETURN
!C         ----------
!C         END
!C
!C
!C         The value of IFLAG should not be changed by FCN unless the
!C         user wants to terminate execution of DNLS1.  In this case, set
!C         IFLAG to a negative integer.
!C
!C
!C       IOPT is an input variable which specifies how the Jacobian will
!C         be calculated.  If IOPT=2 or 3, then the user must supply the
!C         Jacobian, as well as the function values, through the
!C         subroutine FCN.  If IOPT=2, the user supplies the full
!C         Jacobian with one call to FCN.  If IOPT=3, the user supplies
!C         one row of the Jacobian with each call.  (In this manner,
!C         storage can be saved because the full Jacobian is not stored.)
!C         If IOPT=1, the code will approximate the Jacobian by forward
!C         differencing.
!C
!C       M is a positive integer input variable set to the number of
!C         functions.
!C
!C       N is a positive integer input variable set to the number of
!C         variables.  N must not exceed M.
!C
!C       X is an array of length N.  On input, X must contain an initial
!C         estimate of the solution vector.  On output, X contains the
!C         final estimate of the solution vector.
!C
!C       FVEC is an output array of length M which contains the functions
!C         evaluated at the output X.
!C
!C       FJAC is an output array.  For IOPT=1 and 2, FJAC is an M by N
!C         array.  For IOPT=3, FJAC is an N by N array.  The upper N by N
!C         submatrix of FJAC contains an upper triangular matrix R with
!C         diagonal elements of nonincreasing magnitude such that
!C
!C                T     T           T
!C               P *(JAC *JAC)*P = R *R,
!C
!C         where P is a permutation matrix and JAC is the final calcu-
!C         lated Jacobian.  Column J of P is column IPVT(J) (see below)
!C         of the identity matrix.  The lower part of FJAC contains
!C         information generated during the computation of R.
!C
!C       LDFJAC is a positive integer input variable which specifies
!C         the leading dimension of the array FJAC.  For IOPT=1 and 2,
!C         LDFJAC must not be less than M.  For IOPT=3, LDFJAC must not
!C         be less than N.
!C
!C       FTOL is a non-negative input variable.  Termination occurs when
!C         both the actual and predicted relative reductions in the sum
!C         of squares are at most FTOL.  Therefore, FTOL measures the
!C         relative error desired in the sum of squares.  Section 4 con-
!C         tains more details about FTOL.
!C
!C       XTOL is a non-negative input variable.  Termination occurs when
!C         the relative error between two consecutive iterates is at most
!C         XTOL.  Therefore, XTOL measures the relative error desired in
!C         the approximate solution.  Section 4 contains more details
!C         about XTOL.
!C
!C       GTOL is a non-negative input variable.  Termination occurs when
!C         the cosine of the angle between FVEC and any column of the
!C         Jacobian is at most GTOL in absolute value.  Therefore, GTOL
!C         measures the orthogonality desired between the function vector
!C         and the columns of the Jacobian.  Section 4 contains more
!C         details about GTOL.
!C
!C       MAXFEV is a positive integer input variable.  Termination occurs
!C         when the number of calls to FCN to evaluate the functions
!C         has reached MAXFEV.
!C
!C       EPSFCN is an input variable used in determining a suitable step
!C         for the forward-difference approximation.  This approximation
!C         assumes that the relative errors in the functions are of the
!C         order of EPSFCN.  If EPSFCN is less than the machine preci-
!C         sion, it is assumed that the relative errors in the functions
!C         are of the order of the machine precision.  If IOPT=2 or 3,
!C         then EPSFCN can be ignored (treat it as a dummy argument).
!C
!C       DIAG is an array of length N.  If MODE = 1 (see below), DIAG is
!C         internally set.  If MODE = 2, DIAG must contain positive
!C         entries that serve as implicit (multiplicative) scale factors
!C         for the variables.
!C
!C       MODE is an integer input variable.  If MODE = 1, the variables
!C         will be scaled internally.  If MODE = 2, the scaling is speci-
!C         fied by the input DIAG.  Other values of MODE are equivalent
!C         to MODE = 1.
!C
!C       FACTOR is a positive input variable used in determining the ini-
!C         tial step bound.  This bound is set to the product of FACTOR
!C         and the Euclidean norm of DIAG*X if nonzero, or else to FACTOR
!C         itself.  In most cases FACTOR should lie in the interval
!C         (.1,100.).  100. is a generally recommended value.
!C
!C       NPRINT is an integer input variable that enables controlled
!C         printing of iterates if it is positive.  In this case, FCN is
!C         called with IFLAG = 0 at the beginning of the first iteration
!C         and every NPRINT iterations thereafter and immediately prior
!C         to return, with X and FVEC available for printing. Appropriate
!C         print statements must be added to FCN (see example) and
!C         FVEC should not be altered.  If NPRINT is not positive, no
!C         special calls to FCN with IFLAG = 0 are made.
!C
!C       INFO is an integer output variable.  If the user has terminated
!C        execution, INFO is set to the (negative) value of IFLAG.  See
!C        description of FCN and JAC. Otherwise, INFO is set as follows
!C
!C         INFO = 0  improper input parameters.
!C
!C         INFO = 1  both actual and predicted relative reductions in the
!C                   sum of squares are at most FTOL.
!C
!C         INFO = 2  relative error between two consecutive iterates is
!C                   at most XTOL.
!C
!C         INFO = 3  conditions for INFO = 1 and INFO = 2 both hold.
!C
!C         INFO = 4  the cosine of the angle between FVEC and any column
!C                   of the Jacobian is at most GTOL in absolute value.
!C
!C         INFO = 5  number of calls to FCN for function evaluation
!C                   has reached MAXFEV.
!C
!C         INFO = 6  FTOL is too small.  No further reduction in the sum
!C                   of squares is possible.
!C
!C         INFO = 7  XTOL is too small.  No further improvement in the
!C                   approximate solution X is possible.
!C
!C         INFO = 8  GTOL is too small.  FVEC is orthogonal to the
!C                   columns of the Jacobian to machine precision.
!C
!C         Sections 4 and 5 contain more details about INFO.
!C
!C       NFEV is an integer output variable set to the number of calls to
!C         FCN for function evaluation.
!C
!C       NJEV is an integer output variable set to the number of
!C         evaluations of the full Jacobian.  If IOPT=2, only one call to
!C         FCN is required for each evaluation of the full Jacobian.
!C         If IOPT=3, the M calls to FCN are required.
!C         If IOPT=1, then NJEV is set to zero.
!C
!C       IPVT is an integer output array of length N.  IPVT defines a
!C         permutation matrix P such that JAC*P = Q*R, where JAC is the
!C         final calculated Jacobian, Q is orthogonal (not stored), and R
!C         is upper triangular with diagonal elements of nonincreasing
!C         magnitude.  Column J of P is column IPVT(J) of the identity
!C         matrix.
!C
!C       QTF is an output array of length N which contains the first N
!C         elements of the vector (Q transpose)*FVEC.
!C
!C       WA1, WA2, and WA3 are work arrays of length N.
!C
!C       WA4 is a work array of length M.
!C
!C       ERR(M) errors array
!C
!C
!C  C 4. Successful Completion.
!C
!C       The accuracy of DNLS1 is controlled by the convergence parame-
!C       ters FTOL, XTOL, and GTOL.  These parameters are used in tests
!C       which make three types of comparisons between the approximation
!C       X and a solution XSOL.  DNLS1 terminates when any of the tests
!C       is satisfied.  If any of the convergence parameters is less than
!C       the machine precision (as defined by the function R1MACH(4)),
!C       then DNLS1 only attempts to satisfy the test defined by the
!C       machine precision.  Further progress is not usually possible.
!C
!C       The tests assume that the functions are reasonably well behaved,
!C       and, if the Jacobian is supplied by the user, that the functions
!C       and the Jacobian are coded consistently.  If these conditions
!C       are not satisfied, then DNLS1 may incorrectly indicate conver-
!C       gence.  If the Jacobian is coded correctly or IOPT=1,
!C       then the validity of the answer can be checked, for example, by
!C       rerunning DNLS1 with tighter tolerances.
!C
!C       First Convergence Test.  If ENORM(Z) denotes the Euclidean norm
!C         of a vector Z, then this test attempts to guarantee that
!C
!C               ENORM(FVEC) .LE. (1+FTOL)*ENORM(FVECS),
!C
!C         where FVECS denotes the functions evaluated at XSOL.  If this
!C         condition is satisfied with FTOL = 10**(-K), then the final
!C         residual norm ENORM(FVEC) has K significant decimal digits and
!C         INFO is set to 1 (or to 3 if the second test is also satis-
!C         fied).  Unless high precision solutions are required, the
!C         recommended value for FTOL is the square root of the machine
!C         precision.
!C
!C       Second Convergence Test.  If D is the diagonal matrix whose
!C         entries are defined by the array DIAG, then this test attempts
!C         to guarantee that
!C
!C               ENORM(D*(X-XSOL)) .LE. XTOL*ENORM(D*XSOL).
!C
!C         If this condition is satisfied with XTOL = 10**(-K), then the
!C         larger components of D*X have K significant decimal digits and
!C         INFO is set to 2 (or to 3 if the first test is also satis-
!C         fied).  There is a danger that the smaller components of D*X
!C         may have large relative errors, but if MODE = 1, then the
!C         accuracy of the components of X is usually related to their
!C         sensitivity.  Unless high precision solutions are required,
!C         the recommended value for XTOL is the square root of the
!C         machine precision.
!C
!C       Third Convergence Test.  This test is satisfied when the cosine
!C         of the angle between FVEC and any column of the Jacobian at X
!C         is at most GTOL in absolute value.  There is no clear rela-
!C         tionship between this test and the accuracy of DNLS1, and
!C         furthermore, the test is equally well satisfied at other crit-
!C         ical points, namely maximizers and saddle points.  Therefore,
!C         termination caused by this test (INFO = 4) should be examined
!C         carefully.  The recommended value for GTOL is zero.
!C
!C
!C  C 5. Unsuccessful Completion.
!C
!C       Unsuccessful termination of DNLS1 can be due to improper input
!C       parameters, arithmetic interrupts, or an excessive number of
!C       function evaluations.
!C
!C       Improper Input Parameters.  INFO is set to 0 if IOPT .LT. 1
!C         or IOPT .GT. 3, or N .LE. 0, or M .LT. N, or for IOPT=1 or 2
!C         LDFJAC .LT. M, or for IOPT=3 LDFJAC .LT. N, or FTOL .LT. 0.E0,
!C         or XTOL .LT. 0.E0, or GTOL .LT. 0.E0, or MAXFEV .LE. 0, or
!C         FACTOR .LE. 0.E0.
!C
!C       Arithmetic Interrupts.  If these interrupts occur in the FCN
!C         subroutine during an early stage of the computation, they may
!C         be caused by an unacceptable choice of X by DNLS1.  In this
!C         case, it may be possible to remedy the situation by rerunning
!C         DNLS1 with a smaller value of FACTOR.
!C
!C       Excessive Number of Function Evaluations.  A reasonable value
!C         for MAXFEV is 100*(N+1) for IOPT=2 or 3 and 200*(N+1) for
!C         IOPT=1.  If the number of calls to FCN reaches MAXFEV, then
!C         this indicates that the routine is converging very slowly
!C         as measured by the progress of FVEC, and INFO is set to 5.
!C         In this case, it may be helpful to restart DNLS1 with MODE
!C         set to 1.
!C
!C
!C  C 6. Characteristics of the Algorithm.
!C
!C       DNLS1 is a modification of the Levenberg-Marquardt algorithm.
!C       Two of its main characteristics involve the proper use of
!C       implicitly scaled variables (if MODE = 1) and an optimal choice
!C       for the correction.  The use of implicitly scaled variables
!C       achieves scale invariance of DNLS1 and limits the size of the
!C       correction in any direction where the functions are changing
!C       rapidly.  The optimal choice of the correction guarantees (under
!C       reasonable conditions) global convergence from starting points
!C       far from the solution and a fast rate of convergence for
!C       problems with small residuals.
!C
!C       Timing.  The time required by DNLS1 to solve a given problem
!C         depends on M and N, the behavior of the functions, the accu-
!C         racy requested, and the starting point.  The number of arith-
!C         metic operations needed by DNLS1 is about N**3 to process each
!C         evaluation of the functions (call to FCN) and to process each
!C         evaluation of the Jacobian it takes M*N**2 for IOPT=2 (one
!C         call to FCN), M*N**2 for IOPT=1 (N calls to FCN) and
!C         1.5*M*N**2 for IOPT=3 (M calls to FCN).  Unless FCN
!C         can be evaluated quickly, the timing of DNLS1 will be
!C         strongly influenced by the time spent in FCN.
!C
!C       Storage.  DNLS1 requires (M*N + 2*M + 6*N) for IOPT=1 or 2 and
!C         (N**2 + 2*M + 6*N) for IOPT=3 single precision storage
!C         locations and N integer storage locations, in addition to
!C         the storage required by the program.  There are no internally
!C         declared storage arrays.
!C
!C  C *Long Description:
!C
!C  C 7. Example.
!C
!C       The problem is to determine the values of X(1), X(2), and X(3)
!C       which provide the best fit (in the least squares sense) of
!C
!C             X(1) + U(I)/(V(I)*X(2) + W(I)*X(3)),  I = 1, 15
!C
!C       to the data
!C
!C             Y = (0.14,0.18,0.22,0.25,0.29,0.32,0.35,0.39,
!C                  0.37,0.58,0.73,0.96,1.34,2.10,4.39),
!C
!C       where U(I) = I, V(I) = 16 - I, and W(I) = MIN(U(I),V(I)).  The
!C       I-th component of FVEC is thus defined by
!C
!C             Y(I) - (X(1) + U(I)/(V(I)*X(2) + W(I)*X(3))).
!C
!C       **********
!C
!C       PROGRAM TEST
!C  C C
!C  C !C     Driver for DNLS1 example.
!C  C C
!C       INTEGER J,IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV,
!C      *        NWRITE
!C       INTEGER IPVT(3)
!C       DOUBLE PRECISION FTOL,XTOL,GTOL,FACTOR,FNORM,EPSFCN
!C       DOUBLE PRECISION X(3),FVEC(15),FJAC(15,3),DIAG(3),QTF(3),
!C      *     WA1(3),WA2(3),WA3(3),WA4(15)
!C       DOUBLE PRECISION DENORM,D1MACH
!C       EXTERNAL FCN
!C       DATA NWRITE /6/
!C  C C
!C       IOPT = 1
!C       M = 15
!C       N = 3
!C  C C
!C  C !C     The following starting values provide a rough fit.
!C  C C
!C       X(1) = 1.E0
!C       X(2) = 1.E0
!C       X(3) = 1.E0
!C  C C
!C       LDFJAC = 15
!C  C C
!C  C !C     Set FTOL and XTOL to the square root of the machine precision
!C  C !C     and GTOL to zero.  Unless high precision solutions are
!C  C !C     required, these are the recommended settings.
!C  C C
!C       FTOL = SQRT(R1MACH(4))
!C       XTOL = SQRT(R1MACH(4))
!C       GTOL = 0.E0
!C  C C
!C       MAXFEV = 400
!C       EPSFCN = 0.0
!C       MODE = 1
!C       FACTOR = 1.E2
!C       NPRINT = 0
!C  C C
!C       CALL DNLS1(FCN,IOPT,M,N,X,FVEC,FJAC,LDFJAC,FTOL,XTOL,
!C      *           GTOL,MAXFEV,EPSFCN,DIAG,MODE,FACTOR,NPRINT,
!C      *           INFO,NFEV,NJEV,IPVT,QTF,WA1,WA2,WA3,WA4)
!C       FNORM = ENORM(M,FVEC)
!C       WRITE (NWRITE,1000) FNORM,NFEV,NJEV,INFO,(X(J),J=1,N)
!C       STOP
!C  C  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
!C      *        5X,' NUMBER OF FUNCTION EVALUATIONS',I10 //
!C      *        5X,' NUMBER OF JACOBIAN EVALUATIONS',I10 //
!C      *        5X,' EXIT PARAMETER',16X,I10 //
!C      *        5X,' FINAL APPROXIMATE SOLUTION' // 5X,3E15.7)
!C       END
!C       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,DUM,IDUM)
!C  C !C     This is the form of the FCN routine if IOPT=1,
!C  C !C     that is, if the user does not calculate the Jacobian.
!C       INTEGER I,M,N,IFLAG
!C       DOUBLE PRECISION X(N),FVEC(M),Y(15)
!C       DOUBLE PRECISION TMP1,TMP2,TMP3,TMP4
!C       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
!C      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
!C      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
!C      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
!C  C C
!C       IF (IFLAG .NE. 0) GO TO 5
!C  C C
!C  C !C     Insert print statements here when NPRINT is positive.
!C  C C
!C       RETURN
!C     5 CONTINUE
!C       DO 10 I = 1, M
!C          TMP1 = I
!C          TMP2 = 16 - I
!C          TMP3 = TMP1
!C          IF (I .GT. 8) TMP3 = TMP2
!C          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
!C    10    CONTINUE
!C       RETURN
!C       END
!C
!C
!C       Results obtained with different compilers or machines
!C       may be slightly different.
!C
!C       FINAL L2 NORM OF THE RESIDUALS  0.9063596E-01
!C
!C       NUMBER OF FUNCTION EVALUATIONS        25
!C
!C       NUMBER OF JACOBIAN EVALUATIONS         0
!C
!C       EXIT PARAMETER                         1
!C
!C       FINAL APPROXIMATE SOLUTION
!C
!C        0.8241058E-01  0.1133037E+01  0.2343695E+01
!C
!C
!C       For IOPT=2, FCN would be modified as follows to also
!C       calculate the full Jacobian when IFLAG=2.
!C
!C       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
!C  C C
!C  C !C     This is the form of the FCN routine if IOPT=2,
!C  C !C     that is, if the user calculates the full Jacobian.
!C  C C
!C       INTEGER I,LDFJAC,M,N,IFLAG
!C       DOUBLE PRECISION X(N),FVEC(M),FJAC(LDFJAC,N),Y(15)
!C       DOUBLE PRECISION TMP1,TMP2,TMP3,TMP4
!C       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
!C      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
!C      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
!C      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
!C  C C
!C       IF (IFLAG .NE. 0) GO TO 5
!C  C C
!C  C !C     Insert print statements here when NPRINT is positive.
!C  C C
!C       RETURN
!C     5 CONTINUE
!C       IF(IFLAG.NE.1) GO TO 20
!C       DO 10 I = 1, M
!C          TMP1 = I
!C          TMP2 = 16 - I
!C          TMP3 = TMP1
!C          IF (I .GT. 8) TMP3 = TMP2
!C          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
!C    10    CONTINUE
!C       RETURN
!C  C C
!C  C !C     Below, calculate the full Jacobian.
!C  C C
!C    20    CONTINUE
!C  C C
!C       DO 30 I = 1, M
!C          TMP1 = I
!C          TMP2 = 16 - I
!C          TMP3 = TMP1
!C          IF (I .GT. 8) TMP3 = TMP2
!C          TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
!C          FJAC(I,1) = -1.E0
!C          FJAC(I,2) = TMP1*TMP2/TMP4
!C          FJAC(I,3) = TMP1*TMP3/TMP4
!C    30    CONTINUE
!C       RETURN
!C       END
!C
!C
!C       For IOPT = 3, FJAC would be dimensioned as FJAC(3,3),
!C         LDFJAC would be set to 3, and FCN would be written as
!C         follows to calculate a row of the Jacobian when IFLAG=3.
!C
!C       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
!C  C !C     This is the form of the FCN routine if IOPT=3,
!C  C !C     that is, if the user calculates the Jacobian row by row.
!C       INTEGER I,M,N,IFLAG
!C       DOUBLE PRECISION X(N),FVEC(M),FJAC(N),Y(15)
!C       DOUBLE PRECISION TMP1,TMP2,TMP3,TMP4
!C       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
!C      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
!C      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
!C      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
!C  C C
!C       IF (IFLAG .NE. 0) GO TO 5
!C  C C
!C  C !C     Insert print statements here when NPRINT is positive.
!C  C C
!C       RETURN
!C     5 CONTINUE
!C       IF( IFLAG.NE.1) GO TO 20
!C       DO 10 I = 1, M
!C          TMP1 = I
!C          TMP2 = 16 - I
!C          TMP3 = TMP1
!C          IF (I .GT. 8) TMP3 = TMP2
!C          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
!C    10    CONTINUE
!C       RETURN
!C  C C
!C  C !C     Below, calculate the LDFJAC-th row of the Jacobian.
!C  C C
!C    20 CONTINUE
!C
!C       I = LDFJAC
!C          TMP1 = I
!C          TMP2 = 16 - I
!C          TMP3 = TMP1
!C          IF (I .GT. 8) TMP3 = TMP2
!C          TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
!C          FJAC(1) = -1.E0
!C          FJAC(2) = TMP1*TMP2/TMP4
!C          FJAC(3) = TMP1*TMP3/TMP4
!C       RETURN
!C       END
!C
!C***REFERENCES  Jorge J. More, The Levenberg-Marquardt algorithm:
!C                 implementation and theory.  In Numerical Analysis
!C                 Proceedings (Dundee, June 28 - July 1, 1977, G. A.
!C                 Watson, Editor), Lecture Notes in Mathematics 630,
!C                 Springer-Verlag, 1978.
!C***ROUTINES CALLED  D1MACH, DCKDER, DENORM, DFDJC3, DMPAR, DQRFAC,
!C                    DWUPDT, XERMSG
!C***REVISION HISTORY  (YYMMDD)
!C   800301  DATE WRITTEN
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890831  Modified array declarations.  (WRB)
!C   891006  Cosmetic changes to prologue.  (WRB)
!C   891006  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!C   920205  Corrected XERN1 declaration.  (WRB)
!C   920501  Reformatted the REFERENCES section.  (WRB)
!C***END PROLOGUE  DNLS1
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV
      INTEGER IJUNK,NROW,IPVT(*)
      DOUBLE PRECISION FTOL,XTOL,GTOL,FACTOR,EPSFCN
      DOUBLE PRECISION X(*),FVEC(*),FJAC(LDFJAC,*),DIAG(*),QTF(*), &
           WA1(*),WA2(*),WA3(*),WA4(*)
      LOGICAL SING
      EXTERNAL FCN
      INTEGER I,IFLAG,ITER,J,L,MODECH
      DOUBLE PRECISION ACTRED,DELTA,DIRDER,EPSMCH,FNORM,FNORM1,GNORM, &
           ONE,PAR,PNORM,PRERED,P1,P5,P25,P75,P0001,RATIO,SUM,TEMP, &
           TEMP1,TEMP2,XNORM,ZERO
      DOUBLE PRECISION D1MACH,CHKLIM  !!!!! ERR modified from ERR to ERR(M)
      DOUBLE PRECISION  ERR(M)
      CHARACTER*8 XERN1
      CHARACTER*16 XERN3
      SAVE CHKLIM, ONE, P1, P5, P25, P75, P0001, ZERO
!C
      DATA CHKLIM/.1D0/
      DATA ONE,P1,P5,P25,P75,P0001,ZERO &
           /1.0D0,1.0D-1,5.0D-1,2.5D-1,7.5D-1,1.0D-4,0.0D0/
!C***FIRST EXECUTABLE STATEMENT  DNLS1
      EPSMCH = D1MACH(4)
!C
      INFO = 0
      IFLAG = 0
      NFEV = 0
      NJEV = 0
!C
!C     CHECK THE INPUT PARAMETERS FOR ERRORS.
!C
      IF (IOPT .LT. 1 .OR. IOPT .GT. 3 .OR. N .LE. 0 .OR. &
          M .LT. N .OR. LDFJAC .LT. N .OR. FTOL .LT. ZERO &
          .OR. XTOL .LT. ZERO .OR. GTOL .LT. ZERO &
          .OR. MAXFEV .LE. 0 .OR. FACTOR .LE. ZERO) GO TO 300
      IF (IOPT .LT. 3 .AND. LDFJAC .LT. M) GO TO 300
      IF (MODE .NE. 2) GO TO 20
      DO 10 J = 1, N
         IF (DIAG(J) .LE. ZERO) GO TO 300
   10    CONTINUE
   20 CONTINUE
!C
!C     EVALUATE THE FUNCTION AT THE STARTING POINT
!C     AND CALCULATE ITS NORM.
!C
      IFLAG = 1
      IJUNK = 1
      CALL FCN(IFLAG,M,N,X,FVEC,FJAC,IJUNK)
      NFEV = 1
      IF (IFLAG .LT. 0) GO TO 300
      FNORM = DENORM(M,FVEC)
!C
!C     INITIALIZE LEVENBERG-MARQUARDT PARAMETER AND ITERATION COUNTER.
!C
      PAR = ZERO
      ITER = 1
!C
!C     BEGINNING OF THE OUTER LOOP.
!C
   30 CONTINUE
!C
!C        IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.
!C
         IF (NPRINT .LE. 0) GO TO 40
         IFLAG = 0
         IF (MOD(ITER-1,NPRINT) .EQ. 0) &
            CALL FCN(IFLAG,M,N,X,FVEC,FJAC,IJUNK)
         IF (IFLAG .LT. 0) GO TO 300
   40    CONTINUE
!C
!C        CALCULATE THE JACOBIAN MATRIX.
!C
      IF (IOPT .EQ. 3) GO TO 475
!C
!C     STORE THE FULL JACOBIAN USING M*N STORAGE
!C
      IF (IOPT .EQ. 1) GO TO 410
!C
!C     THE USER SUPPLIES THE JACOBIAN
!C
         IFLAG = 2
         CALL FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
         NJEV = NJEV + 1
!C
!C             ON THE FIRST ITERATION, CHECK THE USER SUPPLIED JACOBIAN
!C
         IF (ITER .LE. 1) THEN
            IF (IFLAG .LT. 0) GO TO 300
!C
!C           GET THE INCREMENTED X-VALUES INTO WA1(*).
!C
            MODECH = 1
            CALL DCKDER(M,N,X,FVEC,FJAC,LDFJAC,WA1,WA4,MODECH,ERR)
!C
!C           EVALUATE FUNCTION AT INCREMENTED VALUE AND PUT IN WA4(*).
!C
            IFLAG = 1
            CALL FCN(IFLAG,M,N,WA1,WA4,FJAC,LDFJAC)
            NFEV = NFEV + 1
            IF(IFLAG .LT. 0) GO TO 300
            DO 350 I = 1, M
               MODECH = 2
               CALL DCKDER(1,N,X,FVEC(I),FJAC(I,1),LDFJAC,WA1, &
                    WA4(I),MODECH,ERR)
               IF (ERR(1) .LT. CHKLIM) THEN !!!! ERR modified from ERR to ERR(1)
                  WRITE (XERN1, '(I8)') I
                  WRITE (XERN3, '(1PE15.6)') ERR
                  CALL XERMSG ('SLATEC', 'DNLS1', 'DERIVATIVE OF ' // &
                     'FUNCTION ' // XERN1 // ' MAY BE WRONG, ERR = ' //&
                     XERN3 // ' TOO CLOSE TO 0.', 7, 0)
               ENDIF
  350       CONTINUE
         ENDIF
!C
         GO TO 420
!C
!C     THE CODE APPROXIMATES THE JACOBIAN
!C
410      IFLAG = 1
         CALL DFDJC3(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA4)
         NFEV = NFEV + N
  420    IF (IFLAG .LT. 0) GO TO 300
!C
!C        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
!C
         CALL DQRFAC(M,N,FJAC,LDFJAC,.TRUE.,IPVT,N,WA1,WA2,WA3)
!C
!C        FORM (Q TRANSPOSE)*FVEC AND STORE THE FIRST N COMPONENTS IN
!C        QTF.
!C
         DO 430 I = 1, M
            WA4(I) = FVEC(I)
  430         CONTINUE
         DO 470 J = 1, N
            IF (FJAC(J,J) .EQ. ZERO) GO TO 460
            SUM = ZERO
            DO 440 I = J, M
               SUM = SUM + FJAC(I,J)*WA4(I)
  440          CONTINUE
            TEMP = -SUM/FJAC(J,J)
            DO 450 I = J, M
               WA4(I) = WA4(I) + FJAC(I,J)*TEMP
  450          CONTINUE
  460       CONTINUE
            FJAC(J,J) = WA1(J)
            QTF(J) = WA4(J)
  470       CONTINUE
         GO TO 560
!C
!C        ACCUMULATE THE JACOBIAN BY ROWS IN ORDER TO SAVE STORAGE.
!C        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN MATRIX
!C        CALCULATED ONE ROW AT A TIME, WHILE SIMULTANEOUSLY
!C        FORMING (Q TRANSPOSE)*FVEC AND STORING THE FIRST
!C        N COMPONENTS IN QTF.
!C
  475    DO 490 J = 1, N
            QTF(J) = ZERO
            DO 480 I = 1, N
               FJAC(I,J) = ZERO
  480          CONTINUE
  490        CONTINUE
         DO 500 I = 1, M
            NROW = I
            IFLAG = 3
            CALL FCN(IFLAG,M,N,X,FVEC,WA3,NROW)
            IF (IFLAG .LT. 0) GO TO 300
!C
!C            ON THE FIRST ITERATION, CHECK THE USER SUPPLIED JACOBIAN.
!C
            IF(ITER .GT. 1) GO TO 498
!C
!C            GET THE INCREMENTED X-VALUES INTO WA1(*).
!C
            MODECH = 1
            CALL DCKDER(M,N,X,FVEC,FJAC,LDFJAC,WA1,WA4,MODECH,ERR)
!C
!C            EVALUATE AT INCREMENTED VALUES, IF NOT ALREADY EVALUATED.
!C
            IF(I .NE. 1) GO TO 495
!C
!C            EVALUATE FUNCTION AT INCREMENTED VALUE AND PUT INTO WA4(*).
!C
            IFLAG = 1
            CALL FCN(IFLAG,M,N,WA1,WA4,FJAC,NROW)
            NFEV = NFEV + 1
            IF(IFLAG .LT. 0) GO TO 300
495         CONTINUE
            MODECH = 2
            CALL DCKDER(1,N,X,FVEC(I),WA3,1,WA1,WA4(I),MODECH,ERR)
            IF (ERR(1) .LT. CHKLIM) THEN !!!! ERR modified from ERR to ERR(1)
               WRITE (XERN1, '(I8)') I
               WRITE (XERN3, '(1PE15.6)') ERR
               CALL XERMSG ('SLATEC', 'DNLS1', 'DERIVATIVE OF FUNCTION '&
                  // XERN1 // ' MAY BE WRONG, ERR = ' // XERN3 // &
                  ' TOO CLOSE TO 0.', 7, 0)
            ENDIF
498         CONTINUE
!C
            TEMP = FVEC(I)
            CALL DWUPDT(N,FJAC,LDFJAC,WA3,QTF,TEMP,WA1,WA2)
  500       CONTINUE
         NJEV = NJEV + 1
!C
!C        IF THE JACOBIAN IS RANK DEFICIENT, CALL DQRFAC TO
!C        REORDER ITS COLUMNS AND UPDATE THE COMPONENTS OF QTF.
!C
         SING = .FALSE.
         DO 510 J = 1, N
            IF (FJAC(J,J) .EQ. ZERO) SING = .TRUE.
            IPVT(J) = J
            WA2(J) = DENORM(J,FJAC(1,J))
  510       CONTINUE
         IF (.NOT.SING) GO TO 560
         CALL DQRFAC(N,N,FJAC,LDFJAC,.TRUE.,IPVT,N,WA1,WA2,WA3)
         DO 550 J = 1, N
            IF (FJAC(J,J) .EQ. ZERO) GO TO 540
            SUM = ZERO
            DO 520 I = J, N
               SUM = SUM + FJAC(I,J)*QTF(I)
  520         CONTINUE
            TEMP = -SUM/FJAC(J,J)
            DO 530 I = J, N
               QTF(I) = QTF(I) + FJAC(I,J)*TEMP
  530          CONTINUE
  540       CONTINUE
            FJAC(J,J) = WA1(J)
  550       CONTINUE
  560    CONTINUE
!C
!C        ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
!C        TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
!C
         IF (ITER .NE. 1) GO TO 80
         IF (MODE .EQ. 2) GO TO 60
         DO 50 J = 1, N
            DIAG(J) = WA2(J)
            IF (WA2(J) .EQ. ZERO) DIAG(J) = ONE
   50       CONTINUE
   60    CONTINUE
!C
!C        ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED X
!C        AND INITIALIZE THE STEP BOUND DELTA.
!C
         DO 70 J = 1, N
            WA3(J) = DIAG(J)*X(J)
   70       CONTINUE
         XNORM = DENORM(N,WA3)
         DELTA = FACTOR*XNORM
         IF (DELTA .EQ. ZERO) DELTA = FACTOR
   80    CONTINUE
!C
!C        COMPUTE THE NORM OF THE SCALED GRADIENT.
!C
         GNORM = ZERO
         IF (FNORM .EQ. ZERO) GO TO 170
         DO 160 J = 1, N
            L = IPVT(J)
            IF (WA2(L) .EQ. ZERO) GO TO 150
            SUM = ZERO
            DO 140 I = 1, J
               SUM = SUM + FJAC(I,J)*(QTF(I)/FNORM)
  140          CONTINUE
            GNORM = MAX(GNORM,ABS(SUM/WA2(L)))
  150       CONTINUE
  160       CONTINUE
  170    CONTINUE
!C
!C        TEST FOR CONVERGENCE OF THE GRADIENT NORM.
!C
         IF (GNORM .LE. GTOL) INFO = 4
         IF (INFO .NE. 0) GO TO 300
!C
!C        RESCALE IF NECESSARY.
!C
         IF (MODE .EQ. 2) GO TO 190
         DO 180 J = 1, N
            DIAG(J) = MAX(DIAG(J),WA2(J))
  180       CONTINUE
  190    CONTINUE
!C
!C        BEGINNING OF THE INNER LOOP.
!C
  200    CONTINUE
!C
!C           DETERMINE THE LEVENBERG-MARQUARDT PARAMETER.
!C
            CALL DMPAR(N,FJAC,LDFJAC,IPVT,DIAG,QTF,DELTA,PAR,WA1,WA2, &
                       WA3,WA4)
!C
!C           STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
!C
            DO 210 J = 1, N
               WA1(J) = -WA1(J)
               WA2(J) = X(J) + WA1(J)
               WA3(J) = DIAG(J)*WA1(J)
  210          CONTINUE
            PNORM = DENORM(N,WA3)
!C
!C           ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
!C
            IF (ITER .EQ. 1) DELTA = MIN(DELTA,PNORM)
!C
!C           EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
!C
            IFLAG = 1
            CALL FCN(IFLAG,M,N,WA2,WA4,FJAC,IJUNK)
            NFEV = NFEV + 1
            IF (IFLAG .LT. 0) GO TO 300
            FNORM1 = DENORM(M,WA4)
!C
!C           COMPUTE THE SCALED ACTUAL REDUCTION.
!C
            ACTRED = -ONE
            IF (P1*FNORM1 .LT. FNORM) ACTRED = ONE - (FNORM1/FNORM)**2
!C
!C           COMPUTE THE SCALED PREDICTED REDUCTION AND
!C           THE SCALED DIRECTIONAL DERIVATIVE.
!C
            DO 230 J = 1, N
               WA3(J) = ZERO
               L = IPVT(J)
               TEMP = WA1(L)
               DO 220 I = 1, J
                  WA3(I) = WA3(I) + FJAC(I,J)*TEMP
  220             CONTINUE
  230          CONTINUE
            TEMP1 = DENORM(N,WA3)/FNORM
            TEMP2 = (SQRT(PAR)*PNORM)/FNORM
            PRERED = TEMP1**2 + TEMP2**2/P5
            DIRDER = -(TEMP1**2 + TEMP2**2)
!C
!C           COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED
!C           REDUCTION.
!C
            RATIO = ZERO
            IF (PRERED .NE. ZERO) RATIO = ACTRED/PRERED
!C
!C           UPDATE THE STEP BOUND.
!C
            IF (RATIO .GT. P25) GO TO 240
               IF (ACTRED .GE. ZERO) TEMP = P5
               IF (ACTRED .LT. ZERO) &
                  TEMP = P5*DIRDER/(DIRDER + P5*ACTRED)
               IF (P1*FNORM1 .GE. FNORM .OR. TEMP .LT. P1) TEMP = P1
               DELTA = TEMP*MIN(DELTA,PNORM/P1)
               PAR = PAR/TEMP
               GO TO 260
  240       CONTINUE
               IF (PAR .NE. ZERO .AND. RATIO .LT. P75) GO TO 250
               DELTA = PNORM/P5
               PAR = P5*PAR
  250          CONTINUE
  260       CONTINUE
!C
!C           TEST FOR SUCCESSFUL ITERATION.
!C
            IF (RATIO .LT. P0001) GO TO 290
!C
!C           SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
!C
            DO 270 J = 1, N
               X(J) = WA2(J)
               WA2(J) = DIAG(J)*X(J)
  270          CONTINUE
            DO 280 I = 1, M
               FVEC(I) = WA4(I)
  280          CONTINUE
            XNORM = DENORM(N,WA2)
            FNORM = FNORM1
            ITER = ITER + 1
  290       CONTINUE
!C
!C           TESTS FOR CONVERGENCE.
!C
            IF (ABS(ACTRED) .LE. FTOL .AND. PRERED .LE. FTOL &
                .AND. P5*RATIO .LE. ONE) INFO = 1
            IF (DELTA .LE. XTOL*XNORM) INFO = 2
            IF (ABS(ACTRED) .LE. FTOL .AND. PRERED .LE. FTOL &
                .AND. P5*RATIO .LE. ONE .AND. INFO .EQ. 2) INFO = 3
            IF (INFO .NE. 0) GO TO 300
!C
!C           TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
!C
            IF (NFEV .GE. MAXFEV) INFO = 5
            IF (ABS(ACTRED) .LE. EPSMCH .AND. PRERED .LE. EPSMCH &
                .AND. P5*RATIO .LE. ONE) INFO = 6
            IF (DELTA .LE. EPSMCH*XNORM) INFO = 7
            IF (GNORM .LE. EPSMCH) INFO = 8
            IF (INFO .NE. 0) GO TO 300
!C
!C           END OF THE INNER LOOP. REPEAT IF ITERATION UNSUCCESSFUL.
!C
            IF (RATIO .LT. P0001) GO TO 200
!C
!C        END OF THE OUTER LOOP.
!C
         GO TO 30
  300 CONTINUE
!C
!C     TERMINATION, EITHER NORMAL OR USER IMPOSED.
!C
      IF (IFLAG .LT. 0) INFO = IFLAG
      IFLAG = 0
      IF (NPRINT .GT. 0) CALL FCN(IFLAG,M,N,X,FVEC,FJAC,IJUNK)
      IF (INFO .LT. 0) CALL XERMSG ('SLATEC', 'DNLS1', &
          'EXECUTION TERMINATED BECAUSE USER SET IFLAG NEGATIVE.', 1, 1)
      IF (INFO .EQ. 0) CALL XERMSG ('SLATEC', 'DNLS1', &
          'INVALID INPUT PARAMETER.', 2, 1)
      IF (INFO .EQ. 4) CALL XERMSG ('SLATEC', 'DNLS1', &
          'THIRD CONVERGENCE CONDITION, CHECK RESULTS BEFORE ACCEPTING.', &
          1, 1)
      IF (INFO .EQ. 5) CALL XERMSG ('SLATEC', 'DNLS1', &
          'TOO MANY FUNCTION EVALUATIONS.', 9, 1)
      IF (INFO .GE. 6) CALL XERMSG ('SLATEC', 'DNLS1', &
          'TOLERANCES TOO SMALL, NO FURTHER IMPROVEMENT POSSIBLE.', 3, 1)
      RETURN
!C
!C     LAST CARD OF SUBROUTINE DNLS1.
!C
      END SUBROUTINE DNLS1
!*DECK DQRFAC
      SUBROUTINE DQRFAC (M, N, A, LDA, PIVOT, IPVT, LIPVT, SIGMA, &
                         ACNORM, WA)
!C***BEGIN PROLOGUE  DQRFAC
!C***SUBSIDIARY
!C***PURPOSE  Subsidiary to DNLS1, DNLS1E, DNSQ and DNSQE
!C***LIBRARY   SLATEC
!C***TYPE      DOUBLE PRECISION (QRFAC-S, DQRFAC-D)
!C***AUTHOR  (UNKNOWN)
!C***DESCRIPTION
!C
!C   **** Double Precision version of QRFAC ****
!C
!C     This subroutine uses Householder transformations with column
!C     pivoting (optional) to compute a QR factorization of the
!C     M by N matrix A. That is, DQRFAC determines an orthogonal
!C     matrix Q, a permutation matrix P, and an upper trapezoidal
!C     matrix R with diagonal elements of nonincreasing magnitude,
!C     such that A*P = Q*R. The Householder transformation for
!C     column K, K = 1,2,...,MIN(M,N), is of the form
!C
!C                           T
!C           I - (1/U(K))*U*U
!C
!C     where U has zeros in the first K-1 positions. The form of
!C     this transformation and the method of pivoting first
!C     appeared in the corresponding LINPACK subroutine.
!C
!C     The subroutine statement is
!C
!C       SUBROUTINE DQRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,SIGMA,ACNORM,WA)
!C
!C     where
!C
!C       M is a positive integer input variable set to the number
!C         of rows of A.
!C
!C       N is a positive integer input variable set to the number
!C         of columns of A.
!C
!C       A is an M by N array. On input A contains the matrix for
!C         which the QR factorization is to be computed. On output
!C         the strict upper trapezoidal part of A contains the strict
!C         upper trapezoidal part of R, and the lower trapezoidal
!C         part of A contains a factored form of Q (the non-trivial
!C         elements of the U vectors described above).
!C
!C       LDA is a positive integer input variable not less than M
!C         which specifies the leading dimension of the array A.
!C
!C       PIVOT is a logical input variable. If pivot is set .TRUE.,
!C         then column pivoting is enforced. If pivot is set .FALSE.,
!C         then no column pivoting is done.
!C
!C       IPVT is an integer output array of length LIPVT. IPVT
!C         defines the permutation matrix P such that A*P = Q*R.
!C         Column J of P is column IPVT(J) of the identity matrix.
!C         If pivot is .FALSE., IPVT is not referenced.
!C
!C       LIPVT is a positive integer input variable. If PIVOT is
!C             .FALSE., then LIPVT may be as small as 1. If PIVOT is
!C             .TRUE., then LIPVT must be at least N.
!C
!C       SIGMA is an output array of length N which contains the
!C         diagonal elements of R.
!C
!C       ACNORM is an output array of length N which contains the
!C         norms of the corresponding columns of the input matrix A.
!C         If this information is not needed, then ACNORM can coincide
!C         with SIGMA.
!C
!C       WA is a work array of length N. If pivot is .FALSE., then WA
!C         can coincide with SIGMA.
!C
!C***SEE ALSO  DNLS1, DNLS1E, DNSQ, DNSQE
!C***ROUTINES CALLED  D1MACH, DENORM
!C***REVISION HISTORY  (YYMMDD)
!C   800301  DATE WRITTEN
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890831  Modified array declarations.  (WRB)
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900326  Removed duplicate information from DESCRIPTION section.
!C           (WRB)
!C   900328  Added TYPE section.  (WRB)
!C***END PROLOGUE  DQRFAC
      INTEGER M,N,LDA,LIPVT
      INTEGER IPVT(*)
      LOGICAL PIVOT
      SAVE ONE, P05, ZERO
      DOUBLE PRECISION A(LDA,*),SIGMA(*),ACNORM(*),WA(*)
      INTEGER I,J,JP1,K,KMAX,MINMN
      DOUBLE PRECISION AJNORM,EPSMCH,ONE,P05,SUM,TEMP,ZERO
      DOUBLE PRECISION D1MACH
      DATA ONE,P05,ZERO /1.0D0,5.0D-2,0.0D0/
!C***FIRST EXECUTABLE STATEMENT  DQRFAC
      EPSMCH = D1MACH(4)
!C
!C     COMPUTE THE INITIAL COLUMN NORMS AND INITIALIZE SEVERAL ARRAYS.
!C
      DO 10 J = 1, N
         ACNORM(J) = DENORM(M,A(1,J))
         SIGMA(J) = ACNORM(J)
         WA(J) = SIGMA(J)
         IF (PIVOT) IPVT(J) = J
   10    CONTINUE
!C
!C     REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS.
!C
      MINMN = MIN(M,N)
      DO 110 J = 1, MINMN
         IF (.NOT.PIVOT) GO TO 40
!C
!C        BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION.
!C
         KMAX = J
         DO 20 K = J, N
            IF (SIGMA(K) .GT. SIGMA(KMAX)) KMAX = K
   20       CONTINUE
         IF (KMAX .EQ. J) GO TO 40
         DO 30 I = 1, M
            TEMP = A(I,J)
            A(I,J) = A(I,KMAX)
            A(I,KMAX) = TEMP
   30       CONTINUE
         SIGMA(KMAX) = SIGMA(J)
         WA(KMAX) = WA(J)
         K = IPVT(J)
         IPVT(J) = IPVT(KMAX)
         IPVT(KMAX) = K
   40    CONTINUE
!C
!C        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE
!C        J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR.
!C
         AJNORM = DENORM(M-J+1,A(J,J))
         IF (AJNORM .EQ. ZERO) GO TO 100
         IF (A(J,J) .LT. ZERO) AJNORM = -AJNORM
         DO 50 I = J, M
            A(I,J) = A(I,J)/AJNORM
   50       CONTINUE
         A(J,J) = A(J,J) + ONE
!C
!C        APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS
!C        AND UPDATE THE NORMS.
!C
         JP1 = J + 1
         IF (N .LT. JP1) GO TO 100
         DO 90 K = JP1, N
            SUM = ZERO
            DO 60 I = J, M
               SUM = SUM + A(I,J)*A(I,K)
   60          CONTINUE
            TEMP = SUM/A(J,J)
            DO 70 I = J, M
               A(I,K) = A(I,K) - TEMP*A(I,J)
   70          CONTINUE
            IF (.NOT.PIVOT .OR. SIGMA(K) .EQ. ZERO) GO TO 80
            TEMP = A(J,K)/SIGMA(K)
            SIGMA(K) = SIGMA(K)*SQRT(MAX(ZERO,ONE-TEMP**2))
            IF (P05*(SIGMA(K)/WA(K))**2 .GT. EPSMCH) GO TO 80
            SIGMA(K) = DENORM(M-J,A(JP1,K))
            WA(K) = SIGMA(K)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
         SIGMA(J) = -AJNORM
  110    CONTINUE
      RETURN
!C
!C     LAST CARD OF SUBROUTINE DQRFAC.
!C
      END SUBROUTINE DQRFAC
!*DECK DQRSLV
      SUBROUTINE DQRSLV (N, R, LDR, IPVT, DIAG, QTB, X, SIGMA, WA)
!C***BEGIN PROLOGUE  DQRSLV
!C***SUBSIDIARY
!C***PURPOSE  Subsidiary to DNLS1 and DNLS1E
!C***LIBRARY   SLATEC
!C***TYPE      DOUBLE PRECISION (QRSOLV-S, DQRSLV-D)
!C***AUTHOR  (UNKNOWN)
!C***DESCRIPTION
!C
!C  C  **** Double Precision version of QRSOLV ****
!C
!C     Given an M by N matrix A, an N by N diagonal matrix D,
!C     and an M-vector B, the problem is to determine an X which
!C     solves the system
!C
!C           A*X = B ,     D*X = 0 ,
!C
!C     in the least squares sense.
!C
!C     This subroutine completes the solution of the problem
!C     if it is provided with the necessary information from the
!C     QR factorization, with column pivoting, of A. That is, if
!C     A*P = Q*R, where P is a permutation matrix, Q has orthogonal
!C     columns, and R is an upper triangular matrix with diagonal
!C     elements of nonincreasing magnitude, then DQRSLV expects
!C     the full upper triangle of R, the permutation matrix P,
!C     and the first N components of (Q TRANSPOSE)*B. The system
!C     A*X = B, D*X = 0, is then equivalent to
!C
!C                  T       T
!C           R*Z = Q *B ,  P *D*P*Z = 0 ,
!C
!C     where X = P*Z. If this system does not have full rank,
!C     then a least squares solution is obtained. On output DQRSLV
!C     also provides an upper triangular matrix S such that
!C
!C            T   T               T
!C           P *(A *A + D*D)*P = S *S .
!C
!C     S is computed within DQRSLV and may be of separate interest.
!C
!C     The subroutine statement is
!C
!C       SUBROUTINE DQRSLV(N,R,LDR,IPVT,DIAG,QTB,X,SIGMA,WA)
!C
!C     where
!C
!C       N is a positive integer input variable set to the order of R.
!C
!C       R is an N by N array. On input the full upper triangle
!C         must contain the full upper triangle of the matrix R.
!C         On output the full upper triangle is unaltered, and the
!C         strict lower triangle contains the strict upper triangle
!C         (transposed) of the upper triangular matrix S.
!C
!C       LDR is a positive integer input variable not less than N
!C         which specifies the leading dimension of the array R.
!C
!C       IPVT is an integer input array of length N which defines the
!C         permutation matrix P such that A*P = Q*R. Column J of P
!C         is column IPVT(J) of the identity matrix.
!C
!C       DIAG is an input array of length N which must contain the
!C         diagonal elements of the matrix D.
!C
!C       QTB is an input array of length N which must contain the first
!C         N elements of the vector (Q TRANSPOSE)*B.
!C
!C       X is an output array of length N which contains the least
!C         squares solution of the system A*X = B, D*X = 0.
!C
!C       SIGMA is an output array of length N which contains the
!C         diagonal elements of the upper triangular matrix S.
!C
!C       WA is a work array of length N.
!C
!C***SEE ALSO  DNLS1, DNLS1E
!C***ROUTINES CALLED  (NONE)
!C***REVISION HISTORY  (YYMMDD)
!C   800301  DATE WRITTEN
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890831  Modified array declarations.  (WRB)
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900326  Removed duplicate information from DESCRIPTION section.
!C           (WRB)
!C   900328  Added TYPE section.  (WRB)
!C***END PROLOGUE  DQRSLV
      INTEGER N,LDR
      INTEGER IPVT(*)
      DOUBLE PRECISION R(LDR,*),DIAG(*),QTB(*),X(*),SIGMA(*),WA(*)
      INTEGER I,J,JP1,K,KP1,L,NSING
      DOUBLE PRECISION COS,COTAN,P5,P25,QTBPJ,SIN,SUM,TAN,TEMP,ZERO
      SAVE P5, P25, ZERO
      DATA P5,P25,ZERO /5.0D-1,2.5D-1,0.0D0/
!C***FIRST EXECUTABLE STATEMENT  DQRSLV
      DO 20 J = 1, N
         DO 10 I = J, N
            R(I,J) = R(J,I)
   10       CONTINUE
         X(J) = R(J,J)
         WA(J) = QTB(J)
   20    CONTINUE
!C
!C     ELIMINATE THE DIAGONAL MATRIX D USING A GIVENS ROTATION.
!C
      DO 100 J = 1, N
!C
!C        PREPARE THE ROW OF D TO BE ELIMINATED, LOCATING THE
!C        DIAGONAL ELEMENT USING P FROM THE QR FACTORIZATION.
!C
         L = IPVT(J)
         IF (DIAG(L) .EQ. ZERO) GO TO 90
         DO 30 K = J, N
            SIGMA(K) = ZERO
   30       CONTINUE
         SIGMA(J) = DIAG(L)
!C
!C        THE TRANSFORMATIONS TO ELIMINATE THE ROW OF D
!C        MODIFY ONLY A SINGLE ELEMENT OF (Q TRANSPOSE)*B
!C        BEYOND THE FIRST N, WHICH IS INITIALLY ZERO.
!C
         QTBPJ = ZERO
         DO 80 K = J, N
!C
!C           DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
!C           APPROPRIATE ELEMENT IN THE CURRENT ROW OF D.
!C
            IF (SIGMA(K) .EQ. ZERO) GO TO 70
            IF (ABS(R(K,K)) .GE. ABS(SIGMA(K))) GO TO 40
               COTAN = R(K,K)/SIGMA(K)
               SIN = P5/SQRT(P25+P25*COTAN**2)
               COS = SIN*COTAN
               GO TO 50
   40       CONTINUE
               TAN = SIGMA(K)/R(K,K)
               COS = P5/SQRT(P25+P25*TAN**2)
               SIN = COS*TAN
   50       CONTINUE
!C
!C           COMPUTE THE MODIFIED DIAGONAL ELEMENT OF R AND
!C           THE MODIFIED ELEMENT OF ((Q TRANSPOSE)*B,0).
!C
            R(K,K) = COS*R(K,K) + SIN*SIGMA(K)
            TEMP = COS*WA(K) + SIN*QTBPJ
            QTBPJ = -SIN*WA(K) + COS*QTBPJ
            WA(K) = TEMP
!C
!C           ACCUMULATE THE TRANSFORMATION IN THE ROW OF S.
!C
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 70
            DO 60 I = KP1, N
               TEMP = COS*R(I,K) + SIN*SIGMA(I)
               SIGMA(I) = -SIN*R(I,K) + COS*SIGMA(I)
               R(I,K) = TEMP
   60          CONTINUE
   70       CONTINUE
   80       CONTINUE
   90    CONTINUE
!C
!C        STORE THE DIAGONAL ELEMENT OF S AND RESTORE
!C        THE CORRESPONDING DIAGONAL ELEMENT OF R.
!C
         SIGMA(J) = R(J,J)
         R(J,J) = X(J)
  100    CONTINUE
!C
!C     SOLVE THE TRIANGULAR SYSTEM FOR Z. IF THE SYSTEM IS
!C     SINGULAR, THEN OBTAIN A LEAST SQUARES SOLUTION.
!C
      NSING = N
      DO 110 J = 1, N
         IF (SIGMA(J) .EQ. ZERO .AND. NSING .EQ. N) NSING = J - 1
         IF (NSING .LT. N) WA(J) = ZERO
  110    CONTINUE
      IF (NSING .LT. 1) GO TO 150
      DO 140 K = 1, NSING
         J = NSING - K + 1
         SUM = ZERO
         JP1 = J + 1
         IF (NSING .LT. JP1) GO TO 130
         DO 120 I = JP1, NSING
            SUM = SUM + R(I,J)*WA(I)
  120       CONTINUE
  130    CONTINUE
         WA(J) = (WA(J) - SUM)/SIGMA(J)
  140    CONTINUE
  150 CONTINUE
!C
!C     PERMUTE THE COMPONENTS OF Z BACK TO COMPONENTS OF X.
!C
      DO 160 J = 1, N
         L = IPVT(J)
         X(L) = WA(J)
  160    CONTINUE
      RETURN
!C
!C     LAST CARD OF SUBROUTINE DQRSLV.
!C
      END SUBROUTINE DQRSLV
!*DECK DWUPDT
      SUBROUTINE DWUPDT (N, R, LDR, W, B, ALPHA, COS, SIN)
!C***BEGIN PROLOGUE  DWUPDT
!C***SUBSIDIARY
!C***PURPOSE  Subsidiary to DNLS1 and DNLS1E
!C***LIBRARY   SLATEC
!C***TYPE      DOUBLE PRECISION (RWUPDT-S, DWUPDT-D)
!C***AUTHOR  (UNKNOWN)
!C***DESCRIPTION
!C
!C     Given an N by N upper triangular matrix R, this subroutine
!C     computes the QR decomposition of the matrix formed when a row
!C     is added to R. If the row is specified by the vector W, then
!C     DWUPDT determines an orthogonal matrix Q such that when the
!C     N+1 by N matrix composed of R augmented by W is premultiplied
!C     by (Q TRANSPOSE), the resulting matrix is upper trapezoidal.
!C     The orthogonal matrix Q is the product of N transformations
!C
!C           G(1)*G(2)* ... *G(N)
!C
!C     where G(I) is a Givens rotation in the (I,N+1) plane which
!C     eliminates elements in the I-th plane. DWUPDT also
!C     computes the product (Q TRANSPOSE)*C where C is the
!C     (N+1)-vector (b,alpha). Q itself is not accumulated, rather
!C     the information to recover the G rotations is supplied.
!C
!C     The subroutine statement is
!C
!C       SUBROUTINE DWUPDT(N,R,LDR,W,B,ALPHA,COS,SIN)
!C
!C     where
!C
!C       N is a positive integer input variable set to the order of R.
!C
!C       R is an N by N array. On input the upper triangular part of
!C         R must contain the matrix to be updated. On output R
!C         contains the updated triangular matrix.
!C
!C       LDR is a positive integer input variable not less than N
!C         which specifies the leading dimension of the array R.
!C
!C       W is an input array of length N which must contain the row
!C         vector to be added to R.
!C
!C       B is an array of length N. On input B must contain the
!C         first N elements of the vector C. On output B contains
!C         the first N elements of the vector (Q TRANSPOSE)*C.
!C
!C       ALPHA is a variable. On input ALPHA must contain the
!C         (N+1)-st element of the vector C. On output ALPHA contains
!C         the (N+1)-st element of the vector (Q TRANSPOSE)*C.
!C
!C       COS is an output array of length N which contains the
!C         cosines of the transforming Givens rotations.
!C
!C       SIN is an output array of length N which contains the
!C         sines of the transforming Givens rotations.
!C
!C     **********
!C
!C***SEE ALSO  DNLS1, DNLS1E
!C***ROUTINES CALLED  (NONE)
!C***REVISION HISTORY  (YYMMDD)
!C   800301  DATE WRITTEN
!C   890531  Changed all specific intrinsics to generic.  (WRB)
!C   890831  Modified array declarations.  (WRB)
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900326  Removed duplicate information from DESCRIPTION section.
!C           (WRB)
!C   900328  Added TYPE section.  (WRB)
!C***END PROLOGUE  DWUPDT
      INTEGER N,LDR
      DOUBLE PRECISION ALPHA
      DOUBLE PRECISION R(LDR,*),W(*),B(*),COS(*),SIN(*)
      INTEGER I,J,JM1
      DOUBLE PRECISION COTAN,ONE,P5,P25,ROWJ,TAN,TEMP,ZERO
      SAVE ONE, P5, P25, ZERO
      DATA ONE,P5,P25,ZERO /1.0D0,5.0D-1,2.5D-1,0.0D0/
!C***FIRST EXECUTABLE STATEMENT  DWUPDT
      DO 60 J = 1, N
         ROWJ = W(J)
         JM1 = J - 1
!C
!C        APPLY THE PREVIOUS TRANSFORMATIONS TO
!C        R(I,J), I=1,2,...,J-1, AND TO W(J).
!C
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            TEMP = COS(I)*R(I,J) + SIN(I)*ROWJ
            ROWJ = -SIN(I)*R(I,J) + COS(I)*ROWJ
            R(I,J) = TEMP
   10       CONTINUE
   20    CONTINUE
!C
!C        DETERMINE A GIVENS ROTATION WHICH ELIMINATES W(J).
!C
         COS(J) = ONE
         SIN(J) = ZERO
         IF (ROWJ .EQ. ZERO) GO TO 50
         IF (ABS(R(J,J)) .GE. ABS(ROWJ)) GO TO 30
            COTAN = R(J,J)/ROWJ
            SIN(J) = P5/SQRT(P25+P25*COTAN**2)
            COS(J) = SIN(J)*COTAN
            GO TO 40
   30    CONTINUE
            TAN = ROWJ/R(J,J)
            COS(J) = P5/SQRT(P25+P25*TAN**2)
            SIN(J) = COS(J)*TAN
   40    CONTINUE
!C
!C        APPLY THE CURRENT TRANSFORMATION TO R(J,J), B(J), AND ALPHA.
!C
         R(J,J) = COS(J)*R(J,J) + SIN(J)*ROWJ
         TEMP = COS(J)*B(J) + SIN(J)*ALPHA
         ALPHA = -SIN(J)*B(J) + COS(J)*ALPHA
         B(J) = TEMP
   50    CONTINUE
   60    CONTINUE
      RETURN
!C
!C     LAST CARD OF SUBROUTINE DWUPDT.
!C
      END SUBROUTINE DWUPDT
!!
!!--------------------------------------------------------------------
!!
!*DECK DCOV
      SUBROUTINE DCOVPolinom (M, N, X, FVEC, R, LDR, INFO, WA1, WA2,&
         WA3, WA4, Design, IPIV, InvIPIV, CondNum_O)
!C***BEGIN PROLOGUE  DCOV
!C***PURPOSE  Calculate the covariance matrix for a nonlinear data
!C            fitting problem.  It is intended to be used after a
!C            successful return from either DNLS1 or DNLS1E.
!C***LIBRARY   SLATEC
!C***CATEGORY  K1B1
!C***TYPE      DOUBLE PRECISION (SCOV-S, DCOV-D)
!C***KEYWORDS  COVARIANCE MATRIX, NONLINEAR DATA FITTING,
!C             NONLINEAR LEAST SQUARES
!C***AUTHOR  Hiebert, K. L., (SNLA)
!C***DESCRIPTION
!C
!C  1. Purpose.
!C
!C     DCOV calculates the covariance matrix for a nonlinear data
!C     fitting problem.  It is intended to be used after a
!C     successful return from either DNLS1 or DNLS1E. DCOV
!C     and DNLS1 (and DNLS1E) have compatible parameters.  The
!C     required external subroutine, FCN, is the same
!C     for all three codes, DCOV, DNLS1, and DNLS1E.
!C
!C  2. Subroutine and Type Statements.
!C
!C     SUBROUTINE DCOV(FCN,IOPT,M,N,X,FVEC,R,LDR,INFO,
!C                     WA1,WA2,WA3,WA4)
!C     INTEGER IOPT,M,N,LDR,INFO
!C     DOUBLE PRECISION X(N),FVEC(M),R(LDR,N),WA1(N),WA2(N),WA3(N),WA4(M)
!C     EXTERNAL FCN
!C
!C  3. Parameters. All TYPE REAL parameters are DOUBLE PRECISION
!C
!C      FCN is the name of the user-supplied subroutine which calculates
!C         the functions.  If the user wants to supply the Jacobian
!C         (IOPT=2 or 3), then FCN must be written to calculate the
!C         Jacobian, as well as the functions.  See the explanation
!C         of the IOPT argument below.
!C         If the user wants the iterates printed in DNLS1 or DNLS1E,
!C         then FCN must do the printing.  See the explanation of NPRINT
!C         in DNLS1 or DNLS1E.  FCN must be declared in an EXTERNAL
!C         statement in the calling program and should be written as
!C         follows.
!C
!C         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
!C         INTEGER IFLAG,LDFJAC,M,N
!C         DOUBLE PRECISION X(N),FVEC(M)
!C         ----------
!C         FJAC and LDFJAC may be ignored       , if IOPT=1.
!C         DOUBLE PRECISION FJAC(LDFJAC,N)      , if IOPT=2.
!C         DOUBLE PRECISION FJAC(N)             , if IOPT=3.
!C         ----------
!C           If IFLAG=0, the values in X and FVEC are available
!C           for printing in DNLS1 or DNLS1E.
!C           IFLAG will never be zero when FCN is called by DCOV.
!C           The values of X and FVEC must not be changed.
!C         RETURN
!C         ----------
!C           If IFLAG=1, calculate the functions at X and return
!C           this vector in FVEC.
!C         RETURN
!C         ----------
!C           If IFLAG=2, calculate the full Jacobian at X and return
!C           this matrix in FJAC.  Note that IFLAG will never be 2 unless
!C           IOPT=2.  FVEC contains the function values at X and must
!C           not be altered.  FJAC(I,J) must be set to the derivative
!C           of FVEC(I) with respect to X(J).
!C         RETURN
!C         ----------
!C           If IFLAG=3, calculate the LDFJAC-th row of the Jacobian
!C           and return this vector in FJAC.  Note that IFLAG will
!C           never be 3 unless IOPT=3.  FJAC(J) must be set to
!C           the derivative of FVEC(LDFJAC) with respect to X(J).
!C         RETURN
!C         ----------
!C         END
!C
!C
!C         The value of IFLAG should not be changed by FCN unless the
!C         user wants to terminate execution of DCOV.  In this case, set
!C         IFLAG to a negative integer.
!C
!C
!C       IOPT is an input variable which specifies how the Jacobian will
!C         be calculated.  If IOPT=2 or 3, then the user must supply the
!C         Jacobian, as well as the function values, through the
!C         subroutine FCN.  If IOPT=2, the user supplies the full
!C         Jacobian with one call to FCN.  If IOPT=3, the user supplies
!C         one row of the Jacobian with each call.  (In this manner,
!C         storage can be saved because the full Jacobian is not stored.)
!C         If IOPT=1, the code will approximate the Jacobian by forward
!C         differencing.
!C
!C       M is a positive integer input variable set to the number of
!C         functions.
!C
!C       N is a positive integer input variable set to the number of
!C         variables.  N must not exceed M.
!C
!C       X is an array of length N.  On input X must contain the value
!C         at which the covariance matrix is to be evaluated.  This is
!C         usually the value for X returned from a successful run of
!C         DNLS1 (or DNLS1E).  The value of X will not be changed.
!C
!C    FVEC is an output array of length M which contains the functions
!C         evaluated at X.
!C
!C       R is an output array.  For IOPT=1 and 2, R is an M by N array.
!C         For IOPT=3, R is an N by N array.  On output, if INFO=1,
!C         the upper N by N submatrix of R contains the covariance
!C         matrix evaluated at X.
!C
!C     LDR is a positive integer input variable which specifies
!C         the leading dimension of the array R.  For IOPT=1 and 2,
!C         LDR must not be less than M.  For IOPT=3, LDR must not
!C         be less than N.
!C
!C    INFO is an integer output variable.  If the user has terminated
!C         execution, INFO is set to the (negative) value of IFLAG.  See
!C         description of FCN. Otherwise, INFO is set as follows.
!C
!C         INFO = 0 Improper input parameters (M.LE.0 or N.LE.0).
!C
!C         INFO = 1 Successful return.  The covariance matrix has been
!C                  calculated and stored in the upper N by N
!C                  submatrix of R.
!C
!C         INFO = 2 The Jacobian matrix is singular for the input value
!C                  of X.  The covariance matrix cannot be calculated.
!C                  The upper N by N submatrix of R contains the QR
!C                  factorization of the Jacobian (probably not of
!C                  interest to the user).
!C
!C WA1,WA2 are work arrays of length N.
!C and WA3
!C
!C     WA4 is a work array of length M.
!C
!C***REFERENCES  (NONE)
!C***ROUTINES CALLED  DENORM, DFDJC3, DQRFAC, DWUPDT, XERMSG
!C***REVISION HISTORY  (YYMMDD)
!C   810522  DATE WRITTEN
!C   890831  Modified array declarations.  (WRB)
!C   891006  Cosmetic changes to prologue.  (WRB)
!C   891006  REVISION DATE from Version 3.2
!C   891214  Prologue converted to Version 4.0 format.  (BAB)
!C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!C   900510  Fixed an error message.  (RWC)
!C***END PROLOGUE  DCOV
!C
!C     REVISED 850601-1100
!C     REVISED YYMMDD HHMM
!C
      INTEGER I,IFLAG,INFO,IOPT,J,K,KP1,LDR,M,N,NM1,NROW
      INTEGER,DIMENSION(M) :: IPIV,InvIPIV
      REAL(DOUBLE)   :: X(*),FVEC(*),WA1(*)
      REAL(DOUBLE)   :: WA2(*),WA3(*), WA4(*)
      EXTERNAL FCN
      DOUBLE PRECISION ONE,SIGMA,TEMP,ZERO
      LOGICAL SING
      REAL(DOUBLE)                       :: ABSR,MaxCond,CondNum
      REAL(DOUBLE),OPTIONAL              :: CondNum_O
      REAL(DOUBLE)                       :: Design(:,:),R(:,:)
      SAVE ZERO, ONE
      DATA ZERO/0.D0/,ONE/1.D0/
!C***FIRST EXECUTABLE STATEMENT  DCOV
      !
      CondNum=1.D-7
      IF(PRESENT(CondNum_O)) CondNum=CondNum_O
      !
      SING=.FALSE.
      IFLAG=0
      IF (M.LE.0 .OR. N.LE.0) GO TO 300
!C
!C     CALCULATE SIGMA = (SUM OF THE SQUARED RESIDUALS) / (M-N)
      IFLAG=1
    ! FVEC is known from input
    ! CALL FCN(IFLAG,M,N,X,FVEC,R,LDR,VectX_O)
      IF (IFLAG.LT.0) GO TO 300
      TEMP=DENORM(M,FVEC)
      SIGMA=ONE
      IF (M.NE.N) SIGMA=TEMP*TEMP/DBLE(M-N)
!C
!C     CALCULATE THE JACOBIAN
    ! IF (IOPT.EQ.3) GO TO 200
!C
!C     STORE THE FULL JACOBIAN USING M*N STORAGE
    ! IF (IOPT.EQ.1) GO TO 100
!C
!C     USER SUPPLIES THE JACOBIAN
    ! IFLAG=2
    ! CALL FCN(IFLAG,M,N,X,FVEC,R,LDR,VectX_O)
    ! GO TO 110
!C
!C     CODE APPROXIMATES THE JACOBIAN
    ! 100   CALL DFDJC3(FCN,M,N,X,FVEC,R,LDR,IFLAG,ZERO,WA4)
    ! 110   IF (IFLAG.LT.0) GO TO 300
!!!!!!!! Fill in Jacobian /Design/ from input
      R=Design  
!C
!C     COMPUTE THE QR DECOMPOSITION
    ! CALL DQRFAC(M,N,R,LDR,.FALSE.,IPIV,1,WA1,WA1,WA1)
      CALL DQRFAC(M,N,R,LDR,.TRUE.,IPIV,N,WA1,WA2,WA3)
      DO 120 I=1,N
120   R(I,I)=WA1(I)
!      GO TO 225
!!C
!!C     COMPUTE THE QR FACTORIZATION OF THE JACOBIAN MATRIX CALCULATED ONE
!!C     ROW AT A TIME AND STORED IN THE UPPER TRIANGLE OF R.
!!C     ( (Q TRANSPOSE)*FVEC IS ALSO CALCULATED BUT NOT USED.)
!200   CONTINUE
!      DO 210 J=1,N
!      WA2(J)=ZERO
!      DO 205 I=1,N
!      R(I,J)=ZERO
!205   CONTINUE
!210   CONTINUE
!      IFLAG=3
!      DO 220 I=1,M
!      NROW = I
!      CALL FCN(IFLAG,M,N,X,FVEC,WA1,NROW,Vect_O)
!      IF (IFLAG.LT.0) GO TO 300
!      TEMP=FVEC(I)
!      CALL DWUPDT(N,R,LDR,WA1,WA2,TEMP,WA3,WA4)
!220   CONTINUE
!!C
!!C     CHECK IF R IS SINGULAR.
225   CONTINUE
      !
      MaxCond=Zero
      DO I=1,N
        ABSR=ABS(R(I,I))
        IF(ABSR<MaxCond) MaxCond=ABSR
        DO J=I+1,N
          R(J,I)=Zero
        ENDDO
      ENDDO
      MaxCond=MaxCond*CondNum
      !
!       DO 230 I=1,N
!         IF (ABS(R(I,I))<MaxCond) THEN
!           R(I,I)=Zero
!         ELSE
!           R(I,I)=One/R(I,I)
!         ENDIF
!230   CONTINUE
!      IF (SING) GO TO 300
!C
!C     R IS UPPER TRIANGULAR.  CALCULATE (R TRANSPOSE) INVERSE AND STORE
!C     IN THE UPPER TRIANGLE OF R.
      IF (N.EQ.1) GO TO 275
      NM1=N-1
      DO 270 K=1,NM1
!C
!C     INITIALIZE THE RIGHT-HAND SIDE (WA1(*)) AS THE K-TH COLUMN OF THE
!C     IDENTITY MATRIX.
      WA1(1:N)=ZERO
      WA1(K)=ONE
!C
      R(K,K)=WA1(K)/R(K,K)
      KP1=K+1
      DO 260 I=KP1,N
!C
!C     SUBTRACT R(K,I-1)*R(I-1,*) FROM THE RIGHT-HAND SIDE, WA1(*).
      DO 250 J=I,N
      WA1(J)=WA1(J)-R(K,I-1)*R(I-1,J)
250   CONTINUE
      R(K,I)=WA1(I)/R(I,I)
260   CONTINUE
270   CONTINUE !!!! end of K loop
275   R(N,N)=ONE/R(N,N)
!C
!C     CALCULATE R-INVERSE * (R TRANSPOSE) INVERSE AND STORE IN THE UPPER
!C     TRIANGLE OF R.
      DO 290 I=1,N
      DO 290 J=I,N
      TEMP=ZERO
      DO 280 K=J,N
      TEMP=TEMP+R(I,K)*R(J,K)
280   CONTINUE
      R(I,J)=TEMP*SIGMA
290   CONTINUE
      INFO=1
!
! Permute R back to original system
!
      DO I=1,N
        InvIPIV(IPIV(I))=I
      ENDDO
      DO I=1,N
        DO J=I,N
          Design(InvIPIV(I),InvIPIV(J))=R(I,J)
          Design(InvIPIV(J),InvIPIV(I))=R(I,J)
        ENDDO
      ENDDO
      R=Design
!C
300   CONTINUE
      IF (M.LE.0 .OR. N.LE.0) INFO=0
      IF (IFLAG.LT.0) INFO=IFLAG
      IF (SING) INFO=2
      IF (INFO .LT. 0) CALL XERMSG ('SLATEC', 'DCOV', &
         'EXECUTION TERMINATED BECAUSE USER SET IFLAG NEGATIVE.', 1, 1)
      IF (INFO .EQ. 0) CALL XERMSG ('SLATEC', 'DCOV', &
         'INVALID INPUT PARAMETER.', 2, 1)
      IF (INFO .EQ. 2) CALL XERMSG ('SLATEC', 'DCOV', &
         'SINGULAR JACOBIAN MATRIX, COVARIANCE MATRIX CANNOT BE ' // &
         'CALCULATED.', 1, 1)
      RETURN
      END SUBROUTINE DCOVPolinom
!
!--------------------------------------------------------------------
!
    SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC,VectX_O)
!C
!!C     This is the form of the FCN routine if IOPT=2,
!!C     that is, if the user calculates the full Jacobian.
!C
    INTEGER I,LDFJAC,M,N,IFLAG
    DOUBLE PRECISION X(N),FVEC(M),FJAC(LDFJAC,N)
    REAL(DOUBLE),DIMENSION(:),OPTIONAL :: VectX_O
!C
    END SUBROUTINE FCN
!!
!!--------------------------------------------------------------------
!!
END MODULE SLATEC
