MODULE ONXMemory
   USE GlobalScalars
   USE GlobalCharacters
   IMPLICIT NONE
!------------------ Error codes that we can recover from
   INTEGER, SAVE              :: ErrorCode = -1
   INTEGER, PARAMETER         :: eInit     = -1
   INTEGER, PARAMETER         :: eAOK      =  0   ! Everything is O.K.
   INTEGER, PARAMETER         :: eMAXI     =  1   ! 2-e buffers full 
   INTEGER, PARAMETER         :: eMAXD     =  2   
   INTEGER, PARAMETER         :: eMAXK     =  3   
   INTEGER, PARAMETER         :: eMAXT     =  4
   INTEGER, PARAMETER         :: eMAXDis   =  5
   INTEGER, PARAMETER         :: eMAXPrm   =  6
   INTEGER, PARAMETER         :: eMAXSL    =  7
   INTEGER, PARAMETER         :: eMaxInts  =  8
END MODULE ONXMemory

