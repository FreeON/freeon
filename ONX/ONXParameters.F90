MODULE ONXParameters
   USE GlobalScalars
   USE GlobalCharacters
   IMPLICIT NONE
   INTEGER, PARAMETER      :: MXBUF1   = 83000
   INTEGER, PARAMETER      :: MXBUF2   = 6000 
   INTEGER, PARAMETER      :: MXINTS   = 500000
   INTEGER, PARAMETER      :: MXBATCH  = 500 
   INTEGER, PARAMETER      :: MXCNT    = 10
   INTEGER, PARAMETER      :: MXL      = 2
   INTEGER, PARAMETER      :: MXL2     = 2*MXL
   INTEGER, PARAMETER      :: MXL4     = 4*MXL
   INTEGER, PARAMETER      :: MXB      = 10
   INTEGER, PARAMETER      :: MInfo    = 3
   INTEGER, PARAMETER      :: LenS     = MXL2*(MXL2+3)/2
   INTEGER, PARAMETER      :: Len0     = (MXL+1)*(MXL+2)*(MXL+3)/6
   INTEGER, PARAMETER      :: Len1     = (MXL2+1)*(MXL2+2)*(MXL2+3)/6
   INTEGER, PARAMETER      :: Len2     = MXL2*(MXL2+1)*(MXL2+2)/6
   INTEGER, PARAMETER      :: Len3     = MXL2*(MXL2+1)/2
   INTEGER, PARAMETER      :: Len4     = Len2+Len3+MXL2+1
   INTEGER, PARAMETER      :: Lmax(10) = (/0,1,1,2,2,2,3,3,3,3/)
   INTEGER, PARAMETER      :: Leng(15) = (/1,4,3,10,9,6,20,19,16,10,35,24,31,25,15/)
   REAL(DOUBLE), PARAMETER :: re0      = 0.0D0
   REAL(DOUBLE), PARAMETER :: re1      = 1.0D0
   REAL(DOUBLE), PARAMETER :: re2      = 2.0D0
   REAL(DOUBLE), PARAMETER :: rpi      = 3.1415926535898D0       ! Pi
   REAL(DOUBLE), PARAMETER :: PSch     = 4.9738746919472D0
   REAL(DOUBLE), PARAMETER :: Prev     = 34.986836655254D0       ! 2 Pi^(5/2)
   REAL(DOUBLE), PARAMETER :: Prev1    = 5.914967172796D0        ! Sqrt[2 Pi^(5/2)]
   REAL(DOUBLE), PARAMETER :: Prev2    = Prev1*Prev1

   INTEGER, PARAMETER      :: FilePrm  = 1
   INTEGER, PARAMETER      :: FileDis  = 2
   INTEGER, PARAMETER      :: FileVec  = 3
   INTEGER, PARAMETER      :: FileCos  = 4
   INTEGER, PARAMETER      :: FileSin  = 5
   INTEGER, PARAMETER      :: FileiSA  = 6
   INTEGER, PARAMETER      :: FileBfn  = 7
   INTEGER, PARAMETER      :: FileSrt  = 8
   INTEGER, PARAMETER      :: FWrite   = 0
   INTEGER, PARAMETER      :: FRead    = 1
   INTEGER, PARAMETER      :: FErase   = 2
   INTEGER, PARAMETER      :: TypeI    = 0
   INTEGER, PARAMETER      :: TypeR    = 1

   REAL(DOUBLE), SAVE      :: DisRange = 0.0D0
   REAL(DOUBLE), SAVE      :: DenRange = 0.0D0
   REAL(DOUBLE), SAVE      :: ONXRange = 0.0D0
   REAL(DOUBLE), SAVE      :: MatRange = 0.0D0
   INTEGER, SAVE           :: NRows=0,NCols=0,NElem=0
   INTEGER, SAVE           :: MaxN2=0

   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MONDO_SCRATCH

   TYPE BUFL
     INTEGER        :: PrmL
     INTEGER        :: DisL
     INTEGER        :: VecL
     INTEGER        :: CosL
     INTEGER        :: SinL
     INTEGER        :: IsaL
     INTEGER        :: BfnL
     INTEGER        :: SrtL
   END TYPE

END MODULE ONXParameters
