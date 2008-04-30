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
    String=Space//"SUBROUTINE MD3TRR"//TRIM(IntToChar(Ell)) &
                //"(Nc,PQx,PQy,PQz,AuxR,R)"
    WRITE(77,FMT)String

    String=Space//"IMPLICIT INTEGER(I)"
    WRITE (77,FMT)String
    
!    String=Space//"IMPLICIT NONE"
!    WRITE (77,FMT)String
!
    String=Space//"INTEGER Nc,Mc,Mc1"
    WRITE (77,FMT)String
!
    String=Space//"REAL*8 AuxR(0:"//TRIM(IntToChar(Ell))//",*),R("//TRIM(IntToChar(LHGTF(Ell)))//",*),PQx(*),PQy(*),PQz(*)"

    WRITE (77,FMT)String

!!    WRITE(77,*)'     RETURN'
  END SUBROUTINE PREAMBLE
!
  SUBROUTINE RReq8(IUnroll,Lx_1,XYZ,Lx_2,R8,Lx_3)

    REAL(DOUBLE) :: R8
    INTEGER  :: IUnroll,Lx_1,Lx_2,Lx_3,I
    CHARACTER(LEN=1) :: XYZ
    CHARACTER(LEN=3) :: ICh
    CHARACTER(LEN=5) :: Ch1,Ch2,Ch3
    CHARACTER(LEN=*), PARAMETER :: FMT='(A72)'
    CHARACTER(LEN=72) :: String

    Ch1=IntToChar(Lx_1)
    Ch2=IntToChar(Lx_2)
    Ch3=IntToChar(Lx_3)
    ! 
    DO I=0,IUnroll-1       
       ICh="I_"//IntToChar(I)
       IF(I.GT.0)THEN
          STRING='         '//ICh//'=I_0+'//TRIM(IntToChar(I))
          WRITE(77,FMT)STRING
       ENDIF
       IF(R8==0D0)THEN
!          STRING='         R('//ICh//','//TRIM(Ch1)//')=PQ'//XYZ//'('//ICh//')*R('//ICh//','//TRIM(Ch2)//')'
          STRING='         R('//TRIM(Ch1)//','//TRIM(ICh)//')=PQ'//XYZ//'('//ICh//')*R('//TRIM(Ch2)//','//TRIM(ICh)//')'
       ELSEIF(R8==1D0)THEN
          STRING='         R('//TRIM(Ch1)//','//TRIM(ICh)//')=PQ'//XYZ//'('//ICh//')*R('//TRIM(Ch2)//','//TRIM(ICh) &
               //')+R('//TRIM(Ch3)//','//TRIM(ICh)//')'
       ELSE
!          STRING='         R('//ICh//','//TRIM(Ch1)//')=PQ'//XYZ//'('//ICh//')*R('//ICh//','//TRIM(Ch2) &
!               //')+'//TRIM(DblToMedmChar(R8))//'*R('//ICh//','//TRIM(Ch3)//')'
          STRING='         R('//TRIM(Ch1)//','//TRIM(ICh)//')=PQ'//XYZ//'('//ICh//')*R('//TRIM(Ch2)//','//TRIM(ICh) &
               //')+'//TRIM(DblToShrtChar(R8))//'*R('//TRIM(Ch3)//','//TRIM(ICh)//')'
       ENDIF
       WRITE(77,FMT)STRING
    ENDDO
    !
  END SUBROUTINE RReq8


  SUBROUTINE RAux8(IUnroll,J)

    INTEGER :: IUnroll,I,J
    CHARACTER(LEN=1) :: XYZ
    CHARACTER(LEN=3) :: ICh
    CHARACTER(LEN=5) :: Ch1,Ch2,Ch3
    CHARACTER(LEN=*), PARAMETER :: FMT='(A72)'
    CHARACTER(LEN=72) :: String
    ! 
    DO I=0,IUnroll-1
       ICh="I_"//IntToChar(I)
       STRING='         '//ICh//'=I_0+'//TRIM(IntToChar(I))
       IF(I.GT.0) &
       WRITE(77,FMT)STRING
       STRING='         R(1,'//ICh//')=AuxR('//TRIM(IntToChar(J))//','//ICh//')'
       WRITE(77,FMT)STRING
    ENDDO
    !
  END SUBROUTINE RAux8

END MODULE AuxAux


PROGRAM VRWriter
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
!  USE MemMan
  USE Indexing
  USE Parse
  USE AuxAux
  IMPLICIT NONE
  CHARACTER(LEN=72) :: FileName
  CHARACTER(LEN=2)  :: ChL
  INTEGER           :: BigEll,Ell
  REAL(DOUBLE)      :: R8
  CHARACTER(LEN=72) :: String,String1,String2,String3
  CHARACTER(LEN=*), PARAMETER :: FMT='(A72)'
  CHARACTER(LEN=6) :: Space2 = "      "
  INTEGER :: J,LX_1,LX_2,LX_3,IUnroll,M,N,L

  IUnroll=8

  FileName="MD3TRR.F"
  OPEN(UNIT=77,FILE=FileName, &
       ACCESS='SEQUENTIAL',FORM='FORMATTED', &
       ERR=11,STATUS='NEW')         

  BigEll=3+3+3+4 

  DO Ell=0,BigEll
     IF(Ell.LE.4)THEN
        IUnroll=4
     ELSE
        IUnroll=4
     ENDIF
     ChL=IntToChar(Ell)
     CALL Preamble(77,Ell)          
      !
      IUnroll=1
      String=Space2//'DO 200 I_0=1,Nc'
      WRITE(77,FMT)String
      DO J=Ell,0,-1
         DO L=Ell-J,1,-1
            Lx_1=LMNdex(L,0,0)
            Lx_2=LMNdex(L-1,0,0)
            Lx_3=LMNdex(L-2,0,0)
            R8=DBLE(L-1)
            CALL RReq8(IUnroll,Lx_1,'x',Lx_2,R8,Lx_3)
            DO M=Ell-J-L,1,-1
               Lx_1=LMNdex(M,L,0)
               Lx_2=LMNdex(M-1,L,0)
               Lx_3=LMNdex(M-2,L,0)
               R8=DBLE(M-1)
               CALL RReq8(IUnroll,Lx_1,'x',Lx_2,R8,Lx_3)
               DO N=Ell-J-L-M,1,-1
                  Lx_1=LMNdex(N,M,L)
                  Lx_2=LMNdex(N-1,M,L)
                  Lx_3=LMNdex(N-2,M,L)
                  R8=DBLE(N-1)
                  CALL RReq8(IUnroll,Lx_1,'x',Lx_2,R8,Lx_3)
               ENDDO
               Lx_1=LMNdex(0,M,L)
               Lx_2=LMNdex(0,M-1,L)
               Lx_3=LMNdex(0,M-2,L)
               R8=DBLE(M-1)
               CALL RReq8(IUnroll,Lx_1,'y',Lx_2,R8,Lx_3)
               Lx_1=LMNdex(M,0,L)
               Lx_2=LMNdex(M-1,0,L)
               Lx_3=LMNdex(M-2,0,L)
               R8=DBLE(M-1)
               CALL RReq8(IUnroll,Lx_1,'x',Lx_2,R8,Lx_3)
            ENDDO
            Lx_1=LMNdex(0,L,0)
            Lx_2=LMNdex(0,L-1,0)
            Lx_3=LMNdex(0,L-2,0)
            R8=DBLE(L-1)
            CALL RReq8(IUnroll,Lx_1,'y',Lx_2,R8,Lx_3)
            Lx_1=LMNdex(0,0,L)
            Lx_2=LMNdex(0,0,L-1)
            Lx_3=LMNdex(0,0,L-2)
            R8=DBLE(L-1)
            CALL RReq8(IUnroll,Lx_1,'z',Lx_2,R8,Lx_3)
         ENDDO
         CALL RAux8(IUnroll,J)
      ENDDO
      String='200   CONTINUE'
      WRITE(77,FMT)String
      WRITE(77,*)'     RETURN'
      WRITE(77,*)'     END'
   ENDDO


  STOP
  11 STOP "FUCKED"





END PROGRAM VRWriter
