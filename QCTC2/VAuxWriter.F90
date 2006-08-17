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
    String=Space//"SUBROUTINE AuxIGen"//TRIM(IntToChar(Ell)) &
                //"(Nc,P,Q,ZetaP,ZetaQ,PQx,PQy,PQz,AuxR)"
    WRITE(77,FMT)String

    String=Space//"IMPLICIT DOUBLE PRECISION(G,O,S,F)"
    WRITE (77,FMT)String
!
!    String=Space//"IMPLICIT INTEGER(I-N)"
!    WRITE (77,FMT)String
!
    IF(Ell==0)THEN
       String=Space//"INTEGER Nc,Ns,Mesh,I,J,K,IF0"
    ELSE
       String=Space//"INTEGER Nc,Ns,Mesh,I,J,K,IF0,"//"IF"//TRIM(IntToChar(Ell))
    ENDIF
    WRITE (77,FMT)String
!
    String=Space//"REAL*8  P(3),Q(3,*),PQx(*),PQy(*),PQz(*)"
    WRITE (77,FMT)String

    String=Space//"REAL*8  ZetaP,ZetaQ(*),AuxR(0:"//TRIM(IntToChar(Ell))//",*) "
    WRITE (77,FMT)String
!
    String=Space//"REAL*8 ET,TwoT,T,RTE,RPE,Upq,TwoPi5x2"
    WRITE (77,FMT)String
!
    String=Space//"REAL*8 SqrtT,SqrtPi,One,Two"
    WRITE (77,FMT)String
!
    String=Space//"PARAMETER(SqrtPi=1.7724538509055160273D0)"
    WRITE (77,FMT)String
    String=Space//"PARAMETER(TwoPi5x2=3.4986836655249725693D1)"
    WRITE (77,FMT)String
!
    String=Space//"PARAMETER(One=1D0)"
    WRITE (77,FMT)String
    String=Space//"PARAMETER(Two=2D0)"
    WRITE (77,FMT)String
!
    String=Space//'INCLUDE "GammaGrid_77.Inc"'
    WRITE (77,FMT)String
!
    String=Space//'INCLUDE "GammaDimensions_77.Inc"'
    WRITE (77,FMT)String
!
    String=Space//'REAL*8 F0_0(0:Gamma_Mesh)'
    WRITE (77,FMT)String
    String=Space//'REAL*8 F0_1(0:Gamma_Mesh)'
    WRITE (77,FMT)String
    String=Space//'REAL*8 F0_2(0:Gamma_Mesh)'
    WRITE (77,FMT)String
    String=Space//'REAL*8 F0_3(0:Gamma_Mesh)'
    WRITE (77,FMT)String
    String=Space//'REAL*8 F0_4(0:Gamma_Mesh)'
    WRITE (77,FMT)String
    String=Space//'REAL*8 F0_5(0:Gamma_Mesh)'
    WRITE (77,FMT)String
    String=Space//'REAL*8 F0_6(0:Gamma_Mesh)'
    WRITE (77,FMT)String
!
    IF(Ell.NE.0)THEN
    String=Space//'REAL*8 F'//TRIM(IntToChar(Ell))//'_0(0:Gamma_Mesh)'
    WRITE (77,FMT)String
    String=Space//'REAL*8 F'//TRIM(IntToChar(Ell))//'_1(0:Gamma_Mesh)'
    WRITE (77,FMT)String
    String=Space//'REAL*8 F'//TRIM(IntToChar(Ell))//'_2(0:Gamma_Mesh)'
    WRITE (77,FMT)String
    String=Space//'REAL*8 F'//TRIM(IntToChar(Ell))//'_3(0:Gamma_Mesh)'
    WRITE (77,FMT)String
    String=Space//'REAL*8 F'//TRIM(IntToChar(Ell))//'_4(0:Gamma_Mesh)'
    WRITE (77,FMT)String
    String=Space//'REAL*8 F'//TRIM(IntToChar(Ell))//'_5(0:Gamma_Mesh)'
    WRITE (77,FMT)String
    String=Space//'REAL*8 F'//TRIM(IntToChar(Ell))//'_6(0:Gamma_Mesh)'
    WRITE (77,FMT)String
    ENDIF
!
    String=Space//'INCLUDE "F0_77.Inc"'
    WRITE (77,FMT)String
    IF(Ell.NE.0)THEN
       String=Space//'INCLUDE "F'//TRIM(IntToChar(Ell))//'_77.Inc"'
       WRITE (77,FMT)String
    ENDIF
!
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
       ERR=11,STATUS='NEW')         

  BigEll=3+3+3+4 

  DO Ell=0,12 !BigEll
     ChL=IntToChar(Ell)
     CALL Preamble(77,Ell)     
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
     String=Space2//'   Upq=TwoPi5x2/(RTE*SQRT(RPE))'
     WRITE(77,FMT)String
     String=Space2//'   T=Omega*(PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I))'
     WRITE(77,FMT)String
     String=Space2//'   IF(T.LT.1D0)THEN'
     WRITE(77,FMT)String
     String=Space2//'      J=AINT(T*Gamma_Grid)'
     WRITE(77,FMT)String
     String=Space2//'      G'//TRIM(ChL)//'=F'//TRIM(ChL)//'_0(J)+T*(F' &
                  //TRIM(ChL)//'_1(J)+T*(F'//TRIM(ChL)//'_2(J)+T*(F'//TRIM(ChL)//'_3(J)'
     WRITE(77,FMT)String
     String='     >        +T*(F'//TRIM(ChL)//'_4(J)+T*(F'//TRIM(ChL)//'_5(J)+T*F'//TRIM(ChL)//'_6(J))))))'
     WRITE(77,FMT)String
     IF(Ell.NE.0)THEN
        String=SPACE2//'      ET=DEXP(-T)'
        WRITE(77,FMT)String
        String=Space2//'      TwoT=2.0D0*T'
        WRITE(77,FMT)String
     ENDIF
     DO L=Ell-1,0,-1        
        IF (L.EQ.0) THEN
           STRING = 'G'//TRIM(IntToChar(L))//'=TwoT*G'//TRIM(IntToChar(L+1))//'+ET'
        ELSE
           Flt=1.0D0/(2.0D0*DBLE(L)+1.0D0)
           STRING = 'G'//TRIM(IntToChar(L))//'='//TRIM(DblToChar(Flt)) & 
                //'*(TwoT*G'//TRIM(IntToChar(L+1))//'+ET)'
        END IF
        String=Space2//"      "//String
        WRITE(77,FMT)String
     END DO


     String=Space2//'   ELSEIF(T.LT.Gamma_Switch)THEN'
     WRITE(77,FMT)String
     String=Space2//'      J=AINT(T*Gamma_Grid)'
     WRITE(77,FMT)String
     String=Space2//'      G0=F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)'
     WRITE(77,FMT)String
     String='     >        +T*(F0_4(J)+T*(F0_5(J)+T*F0_6(J))))))'
     WRITE(77,FMT)String
     IF(Ell.NE.0)THEN
        String=SPACE2//'      ET=DEXP(-T)'
        WRITE(77,FMT)String
        String=Space2//'      OneOvr2T=5D-1/T'
        WRITE(77,FMT)String
     ENDIF
     DO L=1,Ell
        Flt=(2D0*DBLE(L)-1D0)
        STRING = 'G'//TRIM(IntToChar(L))//'=OneOvr2T*('//DblToChar(Flt)// '*G'//TRIM(IntToChar(L-1))//'-ET)'
        String=Space2//"      "//String
        WRITE(77,FMT)String
     END DO
     String=Space2//'   ELSE'
     WRITE(77,FMT)String
!
     String=Space2//'      SqrtT=DSQRT(T)'
     WRITE(77,FMT)String
     String=Space2//'      OneOvT=One/T'
     WRITE(77,FMT)String
!
     String=Space2//'      G0=SqrtPi/(Two*SqrtT)'
     WRITE(77,FMT)String
 
     DO L=0,Ell
        IF(L.GT.0)THEN
           String=Space2//'      G'//TRIM(IntToChar(L))//'=G' &
                //TRIM(IntToChar(L-1))//'*'&
                //TRIM(DblToChar(DBLE(L)-0.5D0))//'*OneOvT'
           WRITE(77,FMT)String
        ENDIF
     ENDDO

     String=Space2//'   ENDIF'
     WRITE(77,FMT)String

     String=Space2//"   "//'o1=Upq'
     WRITE(77,FMT)String

     IF(Ell.NE.0)THEN
        String=Space2//"   "//'o2=-2D0*Omega'
        WRITE(77,FMT)String
     ENDIF

     DO L=0,Ell
        String=Space2//"   AuxR("//TRIM(IntToChar(L))//',I)=o1*G'//TRIM(IntToChar(L))
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

