MODULE ONXParameters
   USE GlobalScalars
   USE GlobalCharacters
   IMPLICIT NONE
   REAL(DOUBLE), PARAMETER :: re0      = 0.0D0
   REAL(DOUBLE), PARAMETER :: re1      = 1.0D0
   REAL(DOUBLE), PARAMETER :: re2      = 2.0D0
   REAL(DOUBLE), PARAMETER :: rpi      = 3.1415926535898D0       ! Pi
   REAL(DOUBLE), PARAMETER :: PSch     = 4.9738746919472D0
   REAL(DOUBLE), PARAMETER :: Prev     = 34.986836655254D0       ! 2 Pi^(5/2)
   REAL(DOUBLE), PARAMETER :: Prev1    = 5.914967172796D0        ! Sqrt[2 Pi^(5/2)]
   REAL(DOUBLE), PARAMETER :: Prev2    = Prev1*Prev1
   REAL(DOUBLE), PARAMETER :: Small    = 1.0D-10       
   INTEGER, PARAMETER      :: IDmn(0:4)= (/1,4,17,24,36/)

   REAL(DOUBLE), SAVE      :: DisRange = 0.0D0
   REAL(DOUBLE), SAVE      :: DenRange = 0.0D0
   REAL(DOUBLE), SAVE      :: ONXRange = 0.0D0
   REAL(DOUBLE), SAVE      :: MatRange = 0.0D0

   REAL(DOUBLE), SAVE      :: xNERIs   = 0.0D0
   
   INTEGER, SAVE           :: NRows=0,NCols=0,NElem=0
   INTEGER, SAVE           :: MaxN2=0

   LOGICAL, SAVE           :: Gradient=.FALSE.

END MODULE ONXParameters
