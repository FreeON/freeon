MODULE VKVK
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
!  USE MemMan
  USE Indexing
  USE Parse
  INTEGER,PARAMETER :: HGEll4=3*HGEll-2
  INTEGER, DIMENSION(-1:HGEll4,-1:HGEll4,-1:HGEll4) :: LMNx

#ifdef POINTERS_IN_DERIVED_TYPES
     TYPE(INT_VECT),DIMENSION(:,:),POINTER   :: PLMNx,QLMNx,PQLMNx
#else
     TYPE(INT_VECT),DIMENSION(:,:),ALLOCATABLE :: PLMNx,QLMNx,PQLMNx
#endif


  CONTAINS

  SUBROUTINE HGTFSetUp
    USE Order
    INTEGER L,M,N,EllP,EllQ,LP,MP,NP,LQ,MQ,NQ,CountP,CountQ,CountPQ,PLen,QLen,PQLen,Ix,Iy,Iz,J
    TYPE(INT_VECT) :: Key,Ord,Imp
    INTEGER(INT8),EXTERNAL:: Interleave,HilbertKey
    character(len=72) :: filename

!!$    !------------------------------------------------------------------------------------          
!!$    WRITE(*,*)' BFEll = ',BFEll
!!$    WRITE(*,*)' HGEll = ',HGEll    
!!$    WRITE(*,*)' HGEll4= ',HGEll4
!!$    WRITE(*,*)' HGLen = ',HGLen,HGLen4,HGLen**2

    ! Precompute HGTF addressing 
    LMNx(-1,:,:)=1
    LMNx(:,-1,:)=1
    LMNx(:,:,-1)=1
    DO L=0,HGEll4
       DO M=0,HGEll4-L
          DO N=0,HGEll4-L-M
             LMNx(L,M,N)=LMNDex(L,M,N)
          ENDDO
       ENDDO
    ENDDO

    ALLOCATE(PLMNx(0:HGEll+1,0:HGEll+1))
    ALLOCATE(QLMNx(0:HGEll+1,0:HGEll+1))
    ALLOCATE(PQLMNx(0:HGEll+1,0:HGEll+1))

    ! Precompute all the 
    DO EllP=0,HGEll+1
       DO EllQ=0,HGEll+1
          PLen=LHGTF(EllP)
          QLen=LHGTF(EllQ)
          PQLen=PLen*QLen
          !
          CALL New(PLMNx(EllP,EllQ),PQLen)
          CALL New(QLMNx(EllP,EllQ),PQLen)
          CALL New(PQLMNx(EllP,EllQ),PQLen)
          !
          CountPQ=0                
          DO LP=0,EllP
             DO MP=0,EllP-LP
                DO NP=0,EllP-LP-MP 
                   DO LQ=0,EllQ
                      DO MQ=0,EllQ-LQ
                         DO NQ=0,EllQ-LQ-MQ
                            CountPQ=CountPQ+1
                            QLMNx(EllP,EllQ)%I(CountPQ)=LMNx(LQ,MQ,NQ)
                            PLMNx(EllP,EllQ)%I(CountPQ)=LMNx(LP,MP,NP)
                            PQLMNx(EllP,EllQ)%I(CountPQ)=LMNx(LP+LQ,MP+MQ,NP+NQ)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

          N=PQLen

          CALL New(Key,PQLen)
          CALL New(Ord,PQLen)
          CALL New(Imp,PQLen)

          DO J=1,N
             IF(EllP>EllQ)THEN
                Ix=0 
                Iy=PLMNx(EllP,EllQ)%I(J)
             ELSE
                Ix=QLMNx(EllP,EllQ)%I(J)
                Iy=0
             ENDIF
             Iz=PQLMNx(EllP,EllQ)%I(J)
!             Key%I(J)=Iz
             Key%I(J)=Interleave(21,Ix,Iy,Iz)
             !             Key%I(J)=HilbertKey(21,Key%I(J))
             Ord%I(J)=J
          END DO

          CALL QSORTI(Ord%I(1),N,Key%I(1))
          !
          DO J=1,N
             IMP%I(J)=QLMNx(EllP,EllQ)%I(Ord%I(J))
          ENDDO
          DO J=1,N
             QLMNx(EllP,EllQ)%I(J)=IMP%I(J)
          ENDDO
          !
          DO J=1,N
             IMP%I(J)=PLMNx(EllP,EllQ)%I(Ord%I(J))
          ENDDO
          DO J=1,N
             PLMNx(EllP,EllQ)%I(J)=IMP%I(J)
          ENDDO
          !
          DO J=1,N
             IMP%I(J)=PQLMNx(EllP,EllQ)%I(Ord%I(J))
          ENDDO
          DO J=1,N
             PQLMNx(EllP,EllQ)%I(J)=IMP%I(J)
          ENDDO
!!$
!!$
!!$          ketcontract 327
!!$          md3trr4     180
!!$          total      1086

!!$         filename='Spc'//TRIM(IntToChar(EllP))//'c'//TRIM(IntToChar(EllQ))//'.dat'
!!$
!!$         OPEN(UNIT=77,FILE=FileName, &
!!$              ACCESS='SEQUENTIAL',FORM='FORMATTED', &
!!$              STATUS='NEW')         
!!$
!!$          DO J=1,N
!!$             WRITE(77,22)J,KEY%I(Ord%I(J)),PLMNx(EllP,EllQ)%I(J), &
!!$                        QLMNx(EllP,EllQ)%I(J),PQLMNx(EllP,EllQ)%I(J)
!!$
!!$             WRITE(*,22)J,KEY%I(Ord%I(J)),PLMNx(EllP,EllQ)%I(J), &
!!$                        QLMNx(EllP,EllQ)%I(J),PQLMNx(EllP,EllQ)%I(J)
!!$             22 FORMAT(I5," ",I10,4(I4,' '))
!!$          ENDDO
!!$          close(77)

          CALL Delete(Key)
          CALL Delete(Ord)
          CALL Delete(IMP)

       ENDDO
    ENDDO
    !
  END SUBROUTINE HGTFSetUp

  SUBROUTINE PREAMBLE(U,EllP,EllQ)
    INTEGER :: U,EllP,EllQ,LenQ
    CHARACTER(LEN=72) :: String,String1,String2,String3
    CHARACTER(LEN=*), PARAMETER :: FMT='(A72)'
    CHARACTER(LEN=6) :: CLen,Space = "      "
!
    String=Space//"SUBROUTINE KTrax_"//TRIM(IntToChar(EllP)) &
                                //"_"//TRIM(IntToChar(EllQ)) &
                //"(Nc,R,QCo,Ket)"
    WRITE(77,FMT)String

    String=Space//"IMPLICIT INTEGER(I)"
    WRITE (77,FMT)String
    
!    String=Space//"IMPLICIT NONE"
!    WRITE (77,FMT)String
!
    String=Space//"INTEGER Nc,Mc"
    WRITE (77,FMT)String
!
    String=Space//"REAL*8 R("//TRIM(IntToChar(LHGTF(EllP+ELLQ))) &
         //",*),QCo("//TRIM(IntToChar(LHGTF(ELLQ)))//",*),Ket(*)"
    WRITE (77,FMT)String

!!    WRITE(77,*)'     RETURN'

  END SUBROUTINE PREAMBLE

  SUBROUTINE KRoll(IUnroll,P,Q,PQ)

    INTEGER :: IUnroll,P,Q,PQ,I,Elll
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
       Ch1=IntToChar(P)
       Ch2=IntToChar(Q)
       Ch3=IntToChar(PQ)

       STRING='         Ket('//TRIM(Ch1)//')=Ket('//TRIM(Ch1)//')+R('// &
            TRIM(Ch3)//','//TRIM(ICh)//')*QCo('//TRIM(Ch2)//','//TRIM(ICh)//')'
       WRITE(77,FMT)STRING
    ENDDO
    !
  END SUBROUTINE KRoll


END MODULE VKVK


PROGRAM VKWriter
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
!  USE MemMan
  USE Indexing
  USE Parse
  USE VKVK
  IMPLICIT NONE
  CHARACTER(LEN=72) :: FileName
  CHARACTER(LEN=2)  :: ChL
  INTEGER           :: BigEll,Ell,EllP,EllQ,IP,IQ,IPQ,J,PQLen,Elll
  REAL(DOUBLE)      :: R8
  CHARACTER(LEN=72) :: String,String1,String2,String3
  CHARACTER(LEN=*), PARAMETER :: FMT='(A72)'
  CHARACTER(LEN=6) :: Space2 = "      "
  INTEGER :: I,LX_1,LX_2,LX_3,IUnroll,M,N,L,MPQ

  CALL HGTFSetUp

  FileName="KTrax.F"
  OPEN(UNIT=77,FILE=FileName, &
       ACCESS='SEQUENTIAL',FORM='FORMATTED', &
       ERR=11,STATUS='NEW')         

  DO Elll=0,8

  DO EllP=0,Elll
     DO EllQ=0,Elll
        IF(EllP+EllQ.EQ.Elll.AND.EllP.LT.HGEll+1.AND.EllQ.LT.HGEll)THEN
           WRITE(*,*)' Elll = ',Elll,'HGEll = ',HGEll,' P = ',EllP,' Q = ',EllQ
        IUnroll=4

        CALL Preamble(77,EllP,EllQ)          
!        String=Space2//"Mc=MOD(Nc,"//TRIM(IntToChar(IUnRoll))//")"
!        WRITE(77,FMT)STRING

        PQLen=LHGTF(EllP)*LHGTF(EllQ)

        String=Space2//'DO 100 I_0=1,Nc'
        WRITE(77,FMT)String

        DO J=1,PQLen
           IP=PLMNx(EllP,EllQ)%I(J)
           IQ=QLMNx(EllP,EllQ)%I(J)
           IPQ=PQLMNx(EllP,EllQ)%I(J)
           CALL KRoll(1,IP,IQ,IPQ)
        ENDDO
        String='100   CONTINUE'
        WRITE(77,FMT)String
!!$
!!$        String=Space2//'DO 200 I_0=Mc+1,Nc,'//TRIM(IntToChar(IUnroll))
!!$        WRITE(77,FMT)String
!!$
!!$        IUnroll=4
!!$
!!$        DO J=1,PQLen
!!$           IP=PLMNx(EllP,EllQ)%I(J)
!!$           IQ=QLMNx(EllP,EllQ)%I(J)
!!$           IPQ=PQLMNx(EllP,EllQ)%I(J)
!!$           CALL KRoll(IUnroll,IP,IQ,IPQ)
!!$        ENDDO
!!$
!!$        String='200   CONTINUE'
!!$        WRITE(77,FMT)String

        WRITE(77,*)'     RETURN'
        WRITE(77,*)'     END'
     ENDIF
      !
     ENDDO
  ENDDO
ENDDO

  STOP
  11 STOP "FUCKED"

END PROGRAM VKWriter
