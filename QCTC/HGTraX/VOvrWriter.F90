MODULE AuxAux
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE Parse
  USE Indexing
!  USE MemMan
  IMPLICIT NONE
CONTAINS
  SUBROUTINE PREAMBLE(U,Ell)
    INTEGER :: U,Ell
    CHARACTER(LEN=72) :: String,String1,String2,String3
    CHARACTER(LEN=*), PARAMETER :: FMT='(A72)'
    CHARACTER(LEN=6) :: Space = "      "
!
    String=Space//"SUBROUTINE OvrIGen"//TRIM(IntToChar(Ell)) &
                //"(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)"
    WRITE(77,FMT)String

    String=Space//"IMPLICIT DOUBLE PRECISION(G,O,S,F)"
    WRITE (77,FMT)String
!
!    String=Space//"IMPLICIT INTEGER(I-N)"
!    WRITE (77,FMT)String
!
    String=Space//"INTEGER Nc,Ns,Mesh,I,J,K"
    WRITE (77,FMT)String
!
    String=Space//"REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)"
    WRITE (77,FMT)String

    String=Space//"REAL*8  ZetaP,ZetaQ(*),AuxR(0:"//TRIM(IntToChar(Ell))//",*) "
    WRITE (77,FMT)String

    String=Space//"REAL*8 ET,TwoT,T,RTE,RPE,Upq,Pi"
    WRITE (77,FMT)String

    String=Space//"REAL*8 SqrtT,SqrtPi,One,Two"
    WRITE (77,FMT)String

    String=Space//"PARAMETER(Pi=3.1415926535897932385D0)"
    WRITE (77,FMT)String

    String=Space//"PARAMETER(One=1D0)"
    WRITE (77,FMT)String
    String=Space//"PARAMETER(Two=2D0)"
    WRITE (77,FMT)String


  END SUBROUTINE PREAMBLE
!
END MODULE AuxAux


PROGRAM VAuxWriter
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE Indexing
!  USE MemMan
  USE Parse
  USE AuxAux
  IMPLICIT NONE
  CHARACTER(LEN=72) :: FileName
  CHARACTER(LEN=2)  :: ChL
  INTEGER           :: BigEll,Ell,L
  REAL(DOUBLE)      :: Flt
  CHARACTER(LEN=72) :: String,String1,String2,String3
  CHARACTER(LEN=*), PARAMETER :: FMT='(A72)'
  CHARACTER(LEN=6) :: Space2 = "      "

  FileName="Gamma.F"
  OPEN(UNIT=77,FILE=FileName, &
       ACCESS='SEQUENTIAL',FORM='FORMATTED', &
       ERR=11,STATUS='OLD',POSITION='APPEND')         

  BigEll=3+3+3+4 

  DO Ell=0,12 !BigEll
     ChL=IntToChar(Ell)
     CALL Preamble(77,Ell)     
!     String=Space2//'RETURN'
!     WRITE(77,FMT)String
     String=Space2//'DO 100 I=1,Nc'
     WRITE(77,FMT)String
     String=Space2//'   PQx(I)=-(P(1)-Q(1,I))'
     WRITE(77,FMT)String
     String=Space2//'   PQy(I)=-(P(2)-Q(2,I))'
     WRITE(77,FMT)String
     String=Space2//'   PQz(I)=-(P(3)-Q(3,I))'
     WRITE(77,FMT)String
     String=Space2//'   RTE=ZetaP*ZetaQ(I)'
     WRITE(77,FMT)String
     String=Space2//'   RPE=ZetaP+ZetaQ(I)'
     WRITE(77,FMT)String
     String=Space2//'   Omega=RTE/RPE'
     WRITE(77,FMT)String
     String=Space2//'   Upq=(Pi/RPE)**15D-1'
     WRITE(77,FMT)String
     String=Space2//'   T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))'
     WRITE(77,FMT)String
     String=Space2//'   ET=EXP(-T)'
     WRITE(77,FMT)String
     String=Space2//'   TwoT=2.0D0*T'
     WRITE(77,FMT)String
     String=Space2//'   o1=Upq'
     WRITE(77,FMT)String
     IF(Ell.NE.0)THEN
        String=Space2//'   o2=-2D0*Omega'
        WRITE(77,FMT)String
     ENDIF
     DO L=0,Ell
        String=Space2//"   AuxR("//TRIM(IntToChar(L))//',I)=o1*ET'
        WRITE(77,FMT)String
        IF(L.LE.Ell-1)THEN
           String=Space2//"   o1=o2*o1"
           WRITE(77,FMT)String
        ENDIF
     ENDDO
     WRITE(77,*)'100  CONTINUE'
     !
     WRITE(77,*)'     RETURN'
     WRITE(77,*)'     END'
  ENDDO
  
  STOP
  11 STOP "FUCKED"

END PROGRAM VAuxWriter

