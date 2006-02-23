MODULE AuxAux
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE Parse
  USE Indexing
!  USE MemMan
  IMPLICIT NONE
  INTEGER,PARAMETER :: HGEll4=3*HGEll-2
  INTEGER, DIMENSION(-1:HGEll4,-1:HGEll4,-1:HGEll4) :: LMNx
#ifdef POINTERS_IN_DERIVED_TYPES
  TYPE(INT_VECT),DIMENSION(:,:),POINTER   :: PLMNx,QLMNx,PQLMNx
#else
  TYPE(INT_VECT),DIMENSION(:,:),ALLOCATABLE :: PLMNx,QLMNx,PQLMNx
#endif

CONTAINS

  SUBROUTINE HGTRAX(Tag,EllP,EllQ)
    INTEGER :: L,Ell,EllP,EllQ
    REAL(DOUBLE) :: Flt
    CHARACTER(LEN=72) :: String,String1,String2,String3,Tag
    CHARACTER(LEN=*), PARAMETER :: FMT='(A72)'
    CHARACTER(LEN=2) :: ChL
    !
    Ell=EllP+EllQ
    ChL=IntToChar(Ell)
      
    !
    String="      SUBROUTINE HGTraX"//TRIM(IntToChar(EllP)) &
         //TRIM(IntToChar(EllQ))//"(PQx,PQy,PQz,U,O,Co,Ket)"
    WRITE(77,FMT)String
    String="      IMPLICIT NONE"
    WRITE (77,FMT)String
    String="      REAL*8 PQx,PQy,PQz,U,O,o1,o2,OneOvT"
    WRITE (77,FMT)String
    String="      REAL*8 G(0:"//TRIM(IntToChar(EllP+EllQ))//")"
    WRITE (77,FMT)String
    String="      REAL*8 AuxR(0:"//TRIM(IntToChar(EllP+EllQ))//")"
    WRITE (77,FMT)String
    String="      REAL*8 R("//TRIM(IntToChar(LHGTF(EllP+EllQ)))//")"
    WRITE (77,FMT)String
    String="      REAL*8 Co("//TRIM(IntToChar(LHGTF(EllQ)))//")"
    WRITE (77,FMT)String
    String="      REAL*8 Ket("//TRIM(IntToChar(LHGTF(EllP)))//")"
    WRITE (77,FMT)String
    String="      INTEGER Mesh,I,J,K"
    WRITE (77,FMT)String
    String="      REAL*8 ET,TwoT,T,T2,T3,T4"
    WRITE (77,FMT)String
    String="      REAL*8 SqrtT,SqrtPi,One,Two"
    WRITE (77,FMT)String
    String="      REAL*8 Switch,Grid,F0Asymp,F1Asymp,F2Asymp,F3Asymp"
    WRITE (77,FMT)String
    String="      REAL*8 F4Asymp,F5Asymp,F6Asymp,F7Asymp,F8Asymp"
    WRITE (77,FMT)String
    String="      REAL*8 F9Asymp,F10Asymp,F11Asymp,F12Asymp"
    WRITE (77,FMT)String
    String="      PARAMETER(SqrtPi=1.7724538509055160273D0)"
    WRITE (77,FMT)String
    String="      PARAMETER(One=1D0)"
    WRITE (77,FMT)String
    String="      PARAMETER(Two=2D0)"
    WRITE (77,FMT)String
    String='      INCLUDE "Gamma_Asymptotics.Inc"'
    WRITE (77,FMT)String
    String='      INCLUDE "Mesh.Inc"'
    WRITE (77,FMT)String
    String='      REAL*8 F'//TRIM(IntToChar(Ell))//'_0(0:Mesh)'
    WRITE (77,FMT)String
    String='      REAL*8 F'//TRIM(IntToChar(Ell))//'_1(0:Mesh)'
    WRITE (77,FMT)String
    String='      REAL*8 F'//TRIM(IntToChar(Ell))//'_2(0:Mesh)'
    WRITE (77,FMT)String
    String='      REAL*8 F'//TRIM(IntToChar(Ell))//'_3(0:Mesh)'
    WRITE (77,FMT)String
    String='      INCLUDE "Gamma_'//TRIM(IntToChar(Ell))//'.Inc"'
    WRITE (77,FMT)String
    String='      T=O*(PQx*PQx+PQy*PQy+PQz*PQz)'
    WRITE(77,FMT)String
    String='      IF(T.LT.Switch)THEN'
    WRITE(77,FMT)String
    String='         T2=T*T'
    WRITE(77,FMT)String
    String='         T3=T*T2'
    WRITE(77,FMT)String
    String='         J=AINT(T*Grid)'
    WRITE(77,FMT)String
    String='         G('//TRIM(ChL)//')=(F'//TRIM(ChL)//'_0(J)' &
         //'+T*F'//TRIM(ChL)//'_1(J)'//'+T2*F' &
         //TRIM(ChL)//'_2(J)'//'+T3*F'//TRIM(ChL)//'_3(J))'
    WRITE(77,FMT)String
    IF(Ell.NE.0)THEN
    String='         ET=DEXP(-T)'
       WRITE(77,FMT)String
    String='         TwoT=2.0D0*T'
       WRITE(77,FMT)String
    ENDIF
    DO L=Ell-1,0,-1        
       IF (L.EQ.0) THEN
    STRING='         G('//TRIM(IntToChar(L))//')=TwoT*G('//TRIM(IntToChar(L+1))//')+ET'
       ELSE
          Flt=1.0D0/(2.0D0*DBLE(L)+1.0D0)
    STRING='         G('//TRIM(IntToChar(L))//')='//TRIM(DblToChar(Flt)) & 
               //'*(TwoT*G('//TRIM(IntToChar(L+1))//')+ET)'
       END IF
       WRITE(77,FMT)String
    END DO

    String='         o1=U'
    WRITE(77,FMT)String
    IF(Ell.NE.0)THEN
    String='         o2=-2D0*O'
       WRITE(77,FMT)String
    ENDIF
    DO L=0,Ell
    String="         AuxR("//TRIM(IntToChar(L))//')=o1*G('//TRIM(IntToChar(L))//')'
       WRITE(77,FMT)String
       IF(L.LE.Ell-1)THEN
    String="         o1=o2*o1"
          WRITE(77,FMT)String
       ENDIF
    ENDDO
    String='      ELSE'
    WRITE(77,FMT)String
    String='         SqrtT=DSQRT(T)'
    WRITE(77,FMT)String
    String='         OneOvT=One/T'
    WRITE(77,FMT)String
    String='         G(0)=SqrtPi/(Two*SqrtT)'
    WRITE(77,FMT)String     
    String="         o1=U"
    WRITE(77,FMT)String
    IF(Ell.NE.0)THEN
    String="         o2=-2D0*O"
       WRITE(77,FMT)String
    ENDIF

    DO L=0,Ell
       IF(L.GT.0)THEN
    String='         G('//TRIM(IntToChar(L))//')=G(' &
               //TRIM(IntToChar(L-1))//')*'&
               //TRIM(DblToChar(DBLE(L)-0.5D0))//'*OneOvT'
          WRITE(77,FMT)String
       ENDIF
    String="         AuxR("//TRIM(IntToChar(L))//')=o1*G('//TRIM(IntToChar(L))//')'
       WRITE(77,FMT)String
       IF(L.LE.Ell-1)THEN
    String="         o1=o2*o1"
          WRITE(77,FMT)String
       ENDIF
    ENDDO
    String='      ENDIF'
    WRITE(77,FMT)String
    String="      CALL MD3TRR__"//TRIM(TAG)//'(PQx,PQy,PQz,AuxR,R)'
    WRITE(77,FMT)String
    String="      CALL KTraX__"//TRIM(TAG)//'(R,Co,Ket)'
    WRITE(77,FMT)String
    String="      RETURN"
    WRITE(77,FMT)String
    String="      END"
    WRITE(77,FMT)String
    WRITE(77,*)' '
  END SUBROUTINE HGTRAX

  SUBROUTINE MD3TRR(Tag,Ell)
    CHARACTER(LEN=72) :: Tag
    CHARACTER(LEN=2)  :: ChL
    INTEGER           :: BigEll,Ell
    REAL(DOUBLE)      :: R8
    CHARACTER(LEN=72) :: String,String1,String2,String3
    CHARACTER(LEN=*), PARAMETER :: FMT='(A72)'
    CHARACTER(LEN=6) :: Space2 = "      "
    INTEGER :: J,LX_1,LX_2,LX_3,IUnroll,M,N,L    
    String="      SUBROUTINE MD3TRR__"//TRIM(Tag)//"(PQx,PQy,PQz,AuxR,R)"
    WRITE(77,FMT)String
    String="      IMPLICIT NONE"
    WRITE (77,FMT)String
    String="      REAL*8 AuxR(0:"//TRIM(IntToChar(Ell))//")"
    WRITE (77,FMT)String
    String="      REAL*8 R("//TRIM(IntToChar(LHGTF(Ell)))//")"
    WRITE (77,FMT)String
    String="      REAL*8 PQx,PQy,PQz"
    WRITE (77,FMT)String
    DO J=Ell,0,-1
       DO L=Ell-J,1,-1
          Lx_1=LMNdex(L,0,0)
          Lx_2=LMNdex(L-1,0,0)
          Lx_3=LMNdex(L-2,0,0)
          R8=DBLE(L-1)
          CALL RReq8(Lx_1,'x',Lx_2,R8,Lx_3)
          DO M=Ell-J-L,1,-1
             Lx_1=LMNdex(M,L,0)
             Lx_2=LMNdex(M-1,L,0)
             Lx_3=LMNdex(M-2,L,0)
             R8=DBLE(M-1)
             CALL RReq8(Lx_1,'x',Lx_2,R8,Lx_3)
             DO N=Ell-J-L-M,1,-1
                Lx_1=LMNdex(N,M,L)
                Lx_2=LMNdex(N-1,M,L)
                Lx_3=LMNdex(N-2,M,L)
                R8=DBLE(N-1)
                CALL RReq8(Lx_1,'x',Lx_2,R8,Lx_3)
             ENDDO
             Lx_1=LMNdex(0,M,L)
             Lx_2=LMNdex(0,M-1,L)
             Lx_3=LMNdex(0,M-2,L)
             R8=DBLE(M-1)
             CALL RReq8(Lx_1,'y',Lx_2,R8,Lx_3)
             Lx_1=LMNdex(M,0,L)
             Lx_2=LMNdex(M-1,0,L)
             Lx_3=LMNdex(M-2,0,L)
             R8=DBLE(M-1)
             CALL RReq8(Lx_1,'x',Lx_2,R8,Lx_3)
          ENDDO
          Lx_1=LMNdex(0,L,0)
          Lx_2=LMNdex(0,L-1,0)
          Lx_3=LMNdex(0,L-2,0)
          R8=DBLE(L-1)
          CALL RReq8(Lx_1,'y',Lx_2,R8,Lx_3)
          Lx_1=LMNdex(0,0,L)
          Lx_2=LMNdex(0,0,L-1)
          Lx_3=LMNdex(0,0,L-2)
          R8=DBLE(L-1)
          CALL RReq8(Lx_1,'z',Lx_2,R8,Lx_3)
       ENDDO
       CALL RAux8(J)
    ENDDO
    WRITE(77,*)'     RETURN'
    WRITE(77,*)'     END'
  END SUBROUTINE MD3TRR

  SUBROUTINE RReq8(Lx_1,XYZ,Lx_2,R8,Lx_3)
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
    IF(R8==0D0)THEN
       STRING='      R('//TRIM(Ch1)//')=PQ'//XYZ//'*R('//TRIM(Ch2)//')'
    ELSEIF(R8==1D0)THEN
       STRING='      R('//TRIM(Ch1)//')=PQ'//XYZ//'*R('//TRIM(Ch2)//')+R('//TRIM(Ch3)//')'
    ELSE
       STRING='      R('//TRIM(Ch1)//')=PQ'//XYZ//'*R('//TRIM(Ch2)//')+' &
                          //TRIM(DblToShrtChar(R8))//'*R('//TRIM(Ch3)//')'
    ENDIF
    WRITE(77,FMT)STRING
  END SUBROUTINE RReq8

  SUBROUTINE RAux8(J)

    INTEGER :: IUnroll,I,J
    CHARACTER(LEN=1) :: XYZ
    CHARACTER(LEN=3) :: ICh
    CHARACTER(LEN=5) :: Ch1,Ch2,Ch3
    CHARACTER(LEN=*), PARAMETER :: FMT='(A72)'
    CHARACTER(LEN=72) :: String
    ! 
    STRING='      R(1)=AuxR('//TRIM(IntToChar(J))//')'
    WRITE(77,FMT)STRING
    !
  END SUBROUTINE RAux8


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


  SUBROUTINE KTRAX(Tag,EllP,EllQ)
    INTEGER :: EllP,EllQ,LenQ,IP,IQ,IPQ,J,PQLen
    CHARACTER(LEN=72) :: Tag,String,String1,String2,String3
    CHARACTER(LEN=*), PARAMETER :: FMT='(A72)'
    CHARACTER(LEN=6) :: Ch1,CH2,CH3
    String="      SUBROUTINE KTrax__"//TRIM(Tag)//"(R,Co,Ket)"
    WRITE(77,FMT)String
    String="      IMPLICIT NONE"
    WRITE (77,FMT)String
    String="      REAL*8 R("//TRIM(IntToChar(LHGTF(EllP+ELLQ)))//')'
    WRITE (77,FMT)String
    String="      REAL*8 Co("//TRIM(IntToChar(LHGTF(ELLQ)))//')'
    WRITE (77,FMT)String
    String="      REAL*8 Ket("//TRIM(IntToChar(LHGTF(ELLP)))//')'
    WRITE (77,FMT)String

    PQLen=LHGTF(EllP)*LHGTF(EllQ)
    DO J=1,PQLen
       IP=PLMNx(EllP,EllQ)%I(J)
       IQ=QLMNx(EllP,EllQ)%I(J)
       IPQ=PQLMNx(EllP,EllQ)%I(J)
       Ch1=IntToChar(IP)
       Ch2=IntToChar(IQ)
       Ch3=IntToChar(IPQ)
       STRING='      Ket('//TRIM(Ch1)//')=Ket('//TRIM(Ch1)//')+R('//TRIM(Ch3)//')*Co('//TRIM(Ch2)//')'
       WRITE(77,FMT)STRING
    ENDDO
    WRITE(77,*)'     RETURN'
    WRITE(77,*)'     END'
  END SUBROUTINE KTRAX


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
  CHARACTER(LEN=72) :: Tag,FileName
  INTEGER           :: Ell,EllP,EllQ
  REAL(DOUBLE)      :: R8
  CHARACTER(LEN=72) :: String,String1,String2,String3
  CHARACTER(LEN=*), PARAMETER :: FMT='(A72)'
  INTEGER :: J,M,N,L    

  CALL HGTFSetUp
  DO EllP=0,HGEll
     FileName="HGTraX"//TRIM(IntToChar(EllP))//".F"
     OPEN(UNIT=77,FILE=FileName, &
          ACCESS='SEQUENTIAL',FORM='FORMATTED', &
          ERR=11,STATUS='NEW')         
     DO EllQ=0,HGEll
        Tag=TRIM(IntToChar(EllP))//'_'//TRIM(IntToChar(EllQ))
        CALL MD3TRR(Tag,EllP+EllQ)
        CALL KTRAX(Tag,EllP,EllQ)
     ENDDO
     DO EllQ=0,HGEll
        Tag=TRIM(IntToChar(EllP))//'_'//TRIM(IntToChar(EllQ))
        CALL HGTRAX(Tag,EllP,EllQ)
     ENDDO
     CLOSE(77)
  ENDDO

  FileName="HGTraX.Inc"
  OPEN(UNIT=77,FILE=FileName, &
       ACCESS='SEQUENTIAL',FORM='FORMATTED', &
       ERR=11,STATUS='NEW')         

  STRING='      SELECT CASE(QC%Prim%Ell)'
  WRITE(77,FMT)STRING
  DO EllP=0,HGEll
     STRING='      CASE('//TRIM(IntToChar(EllP))//')'
     WRITE(77,FMT)STRING     
     STRING='      DO N=1,NNear'
     WRITE(77,FMT)STRING     
        STRING='      Q=>Near(N)%P'
     WRITE(77,FMT)STRING     
        STRING='      DO EllQ=0,Q%HERM%Ell'


     WRITE(77,FMT)STRING     
        STRING='         DO I=Q%HERM%Nq(EllQ)'
     WRITE(77,FMT)STRING     
     


        STRING='         ENDDO'
     WRITE(77,FMT)STRING     
        STRING='      ENDDO'
     WRITE(77,FMT)STRING     
        STRING='   ENDDO'
     WRITE(77,FMT)STRING     


     DO EllQ=0,4




  STOP
11 STOP ' Fucked' 


END PROGRAM VRWriter
