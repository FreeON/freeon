MODULE ONXMemory
   USE GlobalScalars
   USE GlobalCharacters
   IMPLICIT NONE
!------------------ Error Codes
   INTEGER, SAVE              :: ErrorCode = -1
   INTEGER, PARAMETER         :: eInit  = -1
   INTEGER, PARAMETER         :: eMXINT =  1
   INTEGER, PARAMETER         :: eMAXD  =  2
   INTEGER, PARAMETER         :: eMAXC  =  3
   INTEGER, PARAMETER         :: eMAXT  =  4
!------------------ Buffer Space Dimensions
   INTEGER, SAVE              :: MXINT,MAXD,MAXC,MAXT,MAXP
    

END MODULE ONXMemory

