PROGRAM VectorRhoMD
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE ProcessControl
  USE InOut
  USE Parse
  USE Indexing
  IMPLICIT NONE
  INTEGER                              :: I,J,K,L,M,N,Ell,LMN
  CHARACTER(LEN=300), DIMENSION(20,400) :: String
  INTEGER,            DIMENSION(  400) :: LMND
  CHARACTER(LEN=1)                     :: ChEll,ChL,ChL1,ChL2,ChDL1
  CHARACTER(LEN=3)                     :: ChLMN,ChLMNLen
  REAL(DOUBLE) :: R
  CHARACTER(Len=1),DIMENSION(3)        :: Carts=(/'X','Y','Z'/)
  CHARACTER(LEN=6),DIMENSION(14)       :: Floats=(/"      ", &
                                                   " 2.D0*", &
                                                   " 3.D0*", &
                                                   " 4.D0*", &
                                                   " 5.D0*", &
                                                   " 6.D0*", &
                                                   " 7.D0*", &
                                                   " 8.D0*", &
                                                   " 9.D0*", &
                                                   "10.D0*", &
                                                   "11.D0*", &
                                                   "12.D0*", &
                                                   "13.D0*", &
                                                   "14.D0*"/)


!==================================================================================
!
!==================================================================================
  WRITE(*,*)' HiCu HiCu HiCu HiCu HiCu HiCu HiCu HiCu HiCu HiCu '
  WRITE(*,*)' '
  WRITE(*,*)' Generating explicit source for density deposition'
  WRITE(*,*)' '
  String(1,1)='      MaxEll = '//TRIM(IntToChar(MaxEll))  &
            //'      MaxEll = '//TRIM(IntToChar(MaxEll))  &
            //'      MaxEll = '//TRIM(IntToChar(MaxEll))
  WRITE(*,*)String(1,1)
  WRITE(*,*)' '
  WRITE(*,*)' HiCu HiCu HiCu HiCu HiCu HiCu HiCu HiCu HiCu HiCu '
!==================================================================================
!
!==================================================================================
  CALL OpenASCII('ExplicitBraElements.Inc',Out,.TRUE.)
  WRITE(Out,*)(Blnk,I=1,9),'Z=Prim%Z'
  WRITE(Out,*)(Blnk,I=1,9),'TwoZ=Two*Z'
  WRITE(Out,*)(Blnk,I=1,9),'SELECT CASE(Prim%Ell)'
  DO Ell=0,MaxEll
     ChEll   =IntToChar(Ell)
     ChLMNLen=IntToChar(LHGTF(Ell))
     String(1,1)=Squish('CASE('//IntToChar(Ell)//')')
     WRITE(Out,*)(Blnk,I=1,9),TRIM(String(1,1))
     WRITE(Out,*)(Blnk,I=1,12),'DO I=1,NGrid '
     WRITE(Out,*)(Blnk,I=1,15),'RPx=Cube%Grid(1,I)-Prim%P(1)'
     WRITE(Out,*)(Blnk,I=1,15),'RPy=Cube%Grid(2,I)-Prim%P(2)'
     WRITE(Out,*)(Blnk,I=1,15),'RPz=Cube%Grid(3,I)-Prim%P(3)'
     WRITE(Out,*)(Blnk,I=1,15),'RP2=RPx**2+RPy**2+RPz**2'
     WRITE(Out,*)(Blnk,I=1,15),'Dist=Z*RP2'
     WRITE(Out,*)(Blnk,I=1,15),'IF(.TRUE.)THEN'
!     WRITE(Out,*)(Blnk,I=1,15),'IF(Dist<PenetratDistanceThreshold)THEN'
     WRITE(Out,*)(Blnk,I=1,18),'Xpt=EXPInv(Dist)'
     WRITE(Out,*)(Blnk,I=1,18),'LambdaX(0)=Xpt'
     WRITE(Out,*)(Blnk,I=1,18),'LambdaX(1)=TwoZ*RPx*Xpt'
     IF(Ell>=1)THEN
        DO L=2,Ell+1
           String(1,1)='LambdaX('//IntToChar(L)                 &
                     //')=TwoZ*(RPx*LambdaX('//IntToChar(L-1)   &
                     //')-'//DblToShrtChar(DBLE(L-1))           &
                     //'*LambdaX('//IntToChar(L-2)//'))'
           WRITE(Out,*)(Blnk,I=1,18),TRIM(Squish(String(1,1)))
        ENDDO
     ENDIF
     WRITE(Out,*)(Blnk,I=1,18),'LambdaY(0)=One'
     WRITE(Out,*)(Blnk,I=1,18),'LambdaY(1)=TwoZ*RPy'
     IF(Ell>=1)THEN
        DO L=2,Ell+1
           String(1,1)='LambdaY('//IntToChar(L)                 &
                     //')=TwoZ*(RPy*LambdaY('//IntToChar(L-1)   &
                     //')-'//DblToShrtChar(DBLE(L-1))           &
                     //'*LambdaY('//IntToChar(L-2)//'))'
           WRITE(Out,*)(Blnk,I=1,18),TRIM(Squish(String(1,1)))
        ENDDO
     ENDIF
     WRITE(Out,*)(Blnk,I=1,18),'LambdaZ(0)=One'
     WRITE(Out,*)(Blnk,I=1,18),'LambdaZ(1)=TwoZ*RPz'
     IF(Ell>=1)THEN
        DO L=2,Ell+1
           String(1,1)='LambdaZ('//IntToChar(L)                 &
                     //')=TwoZ*(RPz*LambdaZ('//IntToChar(L-1)   &
                     //')-'//DblToShrtChar(DBLE(L-1))               &
                     //'*LambdaZ('//IntToChar(L-2)//'))'
           WRITE(Out,*)(Blnk,I=1,18),TRIM(Squish(String(1,1)))
        ENDDO
     ENDIF
     WRITE(Out,*)(Blnk,I=1,18),'Wght=Cube%Wght(I)'
     WRITE(Out,*)(Blnk,I=1,18),'dEdRho=Cube%Vals(I,1)'
     WRITE(Out,*)(Blnk,I=1,18),'dEdAbsGradRho2=Cube%Vals(I,2)'
     WRITE(Out,*)(Blnk,I=1,18),'GradRhoX=Cube%Vals(I,3)'
     WRITE(Out,*)(Blnk,I=1,18),'GradRhoY=Cube%Vals(I,4)'
     WRITE(Out,*)(Blnk,I=1,18),'GradRhoZ=Cube%Vals(I,5)'
!
     K=0
     DO L=0,Ell
        DO M=0,Ell-L
           DO N=0,Ell-M-L
              LMN=LMNDex(L,M,N)
              K=K+1
              LMND(LMN)=K
              String(1,K)='PrimDist=LambdaX('//IntToChar(L)  &
                               //')*LambdaY('//IntToChar(M)  &
                               //')*LambdaZ('//IntToChar(N)//')'
              String(2,K)='GradPrimDistX=-LambdaX('//IntToChar(L+1)  &
                                     //')*LambdaY('//IntToChar(M)  &
                                     //')*LambdaZ('//IntToChar(N)//')'
              String(3,K)='GradPrimDistY=-LambdaX('//IntToChar(L)  &
                                     //')*LambdaY('//IntToChar(M+1)  &
                                     //')*LambdaZ('//IntToChar(N)//')'
              String(4,K)='GradPrimDistZ=-LambdaX('//IntToChar(L)  &
                                     //')*LambdaY('//IntToChar(M)  &
                                     //')*LambdaZ('//IntToChar(N+1)//')'
              String(5,K)='GradBraRhoDot=GradRhoX*GradPrimDistX &'
              String(6,K)='+GradRhoY*GradPrimDistY+GradRhoZ*GradPrimDistZ'
              String(7,K)='Prim%Ket('//IntToChar(LMN)//')=Prim%Ket('//IntToChar(LMN)//') & '
              String(8,K)='+Wght*(dEdRho*PrimDist+dEdAbsGradRho2*GradBraRhoDot)'
        WRITE(Out,*)(Blnk,I=1,18),TRIM(Squish(String(1,K)))
        WRITE(Out,*)(Blnk,I=1,18),TRIM(Squish(String(2,K)))
        WRITE(Out,*)(Blnk,I=1,18),TRIM(Squish(String(3,K)))
        WRITE(Out,*)(Blnk,I=1,18),TRIM(Squish(String(4,K)))
        WRITE(Out,*)(Blnk,I=1,18),TRIM(Squish(String(5,K)))
        WRITE(Out,*)(Blnk,I=1,18),TRIM(Squish(String(6,K)))
        WRITE(Out,*)(Blnk,I=1,18),TRIM(Squish(String(7,K)))
        WRITE(Out,*)(Blnk,I=1,18),TRIM(Squish(String(8,K)))

             ENDDO
         ENDDO
     ENDDO
     WRITE(Out,*)(Blnk,I=1,15),'ENDIF'
     WRITE(Out,*)(Blnk,I=1,12),'ENDDO'
  ENDDO
  String(1,1)=Squish('CASE('//IntToChar(MaxEll+1)//':)')
  WRITE(Out,*)(Blnk,I=1,9),TRIM(String(1,1))
  WRITE(Out,*)(Blnk,I=1,12),'CALL Halt("Increase MaxEll and remake")'
  WRITE(Out,*)(Blnk,I=1,9),'END SELECT'
  CLOSE(Out)
!=================================================================================
!
!==================================================================================
  CALL OpenASCII('ExplicitLeafContribution.Inc',Out,.TRUE.)
!  WRITE(Out,*)(Blnk,I=1,9),'RhoV=Zero'
  WRITE(Out,*)(Blnk,I=1,6),'RQx=Grid(1,IGrid)-Node%Qx'
  WRITE(Out,*)(Blnk,I=1,6),'RQy=Grid(2,IGrid)-Node%Qy'
  WRITE(Out,*)(Blnk,I=1,6),'RQz=Grid(3,IGrid)-Node%Qz'
  WRITE(Out,*)(Blnk,I=1,6),'RQ2=RQx*RQx+RQy*RQy+RQz*RQz'
  WRITE(Out,*)(Blnk,I=1,6),'X=Node%Zeta*RQ2'
  WRITE(Out,*)(Blnk,I=1,6),'IF(X>Exp_Switch)RETURN'
  WRITE(Out,*)(Blnk,I=1,6),'J=AINT(X*Exp_Grid)'
  WRITE(Out,*)(Blnk,I=1,6),'Xpt=Exp_0(J)+X*(Exp_1(J)+X*(Exp_2(J)+X*(Exp_3(J)+X*Exp_4(J))))'
  WRITE(Out,*)(Blnk,I=1,6),'TwoZ=Two*Node%Zeta'
  WRITE(Out,*)(Blnk,I=1,9),'LambdaX(0)=Xpt'
  WRITE(Out,*)(Blnk,I=1,9),'LambdaX(1)=TwoZ*RQx*Xpt'
  WRITE(Out,*)(Blnk,I=1,9),'LambdaY(0)=One'
  WRITE(Out,*)(Blnk,I=1,9),'LambdaY(1)=TwoZ*RQy'
  WRITE(Out,*)(Blnk,I=1,9),'LambdaZ(0)=One'
  WRITE(Out,*)(Blnk,I=1,9),'LambdaZ(1)=TwoZ*RQz'
  WRITE(Out,*)(Blnk,I=1,6),'SELECT CASE(Node%Ell)'
  DO Ell=0,MaxEll
     ChEll   =IntToChar(Ell)
     ChLMNLen=IntToChar(LHGTF(Ell))
     String(1,1)=Squish('CASE('//IntToChar(Ell)//')')
     WRITE(Out,*)(Blnk,I=1,6),TRIM(String(1,1))
     IF(Ell>=1)THEN
        DO K=1,3
           DO L=2,Ell+1
              String(1,1)='Lambda'//Carts(K)//'('//IntToChar(L)           &
                        //')=TwoZ*(RQ'//Carts(K)//'*Lambda'//Carts(K)     &
                        //'('//IntToChar(L-1)//')-'//Floats(L-1)          &
                        //'Lambda'//Carts(K)//'('//IntToChar(L-2)//'))'
               WRITE(Out,*)(Blnk,I=1,9),TRIM(Squish(String(1,1)))
           ENDDO
        ENDDO
     ENDIF
     K=0
     DO L=0,Ell
        DO M=0,Ell-L
           DO N=0,Ell-M-L
              LMN=LMNDex(L,M,N)
              K=K+1
              LMND(LMN)=K
#ifdef MULTIPLE_DIST
              String(1,LMN)='Co=Node%Co(KC+'//IntToChar(LMN)//')'
#else
              String(1,LMN)='Co=Node%Co('//IntToChar(LMN)//')'
#endif
              String(2,LMN)='RhoV(1,IGrid)=RhoV(1,IGrid)+LambdaX('//IntToChar(L) &
                          //')*LambdaY('//IntToChar(M)//')*LambdaZ('//IntToChar(N)//')*Co'
              String(3,LMN)='RhoV(2,IGrid)=RhoV(2,IGrid)-LambdaX('//IntToChar(L+1) &
                          //')*LambdaY('//IntToChar(M)//')*LambdaZ('//IntToChar(N)//')*Co'
              String(4,LMN)='RhoV(3,IGrid)=RhoV(3,IGrid)-LambdaX('//IntToChar(L) &
                          //')*LambdaY('//IntToChar(M+1)//')*LambdaZ('//IntToChar(N)//')*Co'
              String(5,LMN)='RhoV(4,IGrid)=RhoV(4,IGrid)-LambdaX('//IntToChar(L) &
                          //')*LambdaY('//IntToChar(M)//')*LambdaZ('//IntToChar(N+1)//')*Co'
            ENDDO
         ENDDO
     ENDDO
     DO J=1,LHGTF(Ell)
        K=LMND(J)
        WRITE(Out,*)(Blnk,I=1,9),TRIM(Squish(String(1,J)))
        WRITE(Out,*)(Blnk,I=1,9),TRIM(Squish(String(2,J)))
        WRITE(Out,*)(Blnk,I=1,9),TRIM(Squish(String(3,J)))
        WRITE(Out,*)(Blnk,I=1,9),TRIM(Squish(String(4,J)))
        WRITE(Out,*)(Blnk,I=1,9),TRIM(Squish(String(5,J)))
     ENDDO
#ifdef MULTIPLE_DIST
     WRITE(Out,*)(Blnk,I=1,9),'ENDDO'
#endif     
  ENDDO
#ifdef MULTIPLE_DIST
#else
!  WRITE(Out,*)(Blnk,I=1,6),TRIM(Squish('CASE('//IntToChar(MaxEll+1)//':)'))
!  WRITE(Out,*)(Blnk,I=1,9),'CALL MondoHalt("Increase MaxEll and remake")'
  WRITE(Out,*)(Blnk,I=1,6),'END SELECT'
#endif

  CLOSE(Out)
!==================================================================================
!
!==================================================================================
  CALL OpenASCII('ExplicitLeafPopulation.Inc',Out,.TRUE.)
!  WRITE(Out,*)(Blnk,I=1,9),'Population=Zero'
#ifdef MULTIPLE_DIST
#else
   WRITE(Out,*)(Blnk,I=1,12),'Z=Node%Zeta'
   WRITE(Out,*)(Blnk,I=1,12),'LQx=Box%BndBox(1,1)-Node%Qx'
   WRITE(Out,*)(Blnk,I=1,12),'LQy=Box%BndBox(2,1)-Node%Qy'
   WRITE(Out,*)(Blnk,I=1,12),'LQz=Box%BndBox(3,1)-Node%Qz'
   WRITE(Out,*)(Blnk,I=1,12),'UQx=Box%BndBox(1,2)-Node%Qx'
   WRITE(Out,*)(Blnk,I=1,12),'UQy=Box%BndBox(2,2)-Node%Qy'
   WRITE(Out,*)(Blnk,I=1,12),'UQz=Box%BndBox(3,2)-Node%Qz'
   WRITE(Out,*)(Blnk,I=1,12),'LQ2=LQx*LQx+LQy*LQy+LQz*LQz'
   WRITE(Out,*)(Blnk,I=1,12),'UQ2=UQx*UQx+UQy*UQy+UQz*UQz'
   WRITE(Out,*)(Blnk,I=1,12),'TwoZ=Two*Z'
   WRITE(Out,*)(Blnk,I=1,12),'SqZ=SQRT(Z)'
   WRITE(Out,*)(Blnk,I=1,12),'CoFact=SqrtPi/(Two*SqZ)'
   WRITE(Out,*)(Blnk,I=1,12),'LLambdaX(0)=-CoFact*ERF(SqZ*LQx)'
   WRITE(Out,*)(Blnk,I=1,12),'LLambdaY(0)=-CoFact*ERF(SqZ*LQy)'
   WRITE(Out,*)(Blnk,I=1,12),'LLambdaZ(0)=-CoFact*ERF(SqZ*LQz)'
   WRITE(Out,*)(Blnk,I=1,12),'ULambdaX(0)=-CoFact*ERF(SqZ*UQx)'
   WRITE(Out,*)(Blnk,I=1,12),'ULambdaY(0)=-CoFact*ERF(SqZ*UQy)'
   WRITE(Out,*)(Blnk,I=1,12),'ULambdaZ(0)=-CoFact*ERF(SqZ*UQz)'
   WRITE(Out,*)(Blnk,I=1,6),'SELECT CASE(Node%Ell)'
#endif
  DO Ell=0,MaxEll
     ChEll   =IntToChar(Ell)
     ChLMNLen=IntToChar(LHGTF(Ell))
#ifdef MULTIPLE_DIST
     WRITE(Out,*)(Blnk,I=1,9),TRIM(Squish('JQ=Node%Qdex('//IntToChar(Ell)//')'))
     WRITE(Out,*)(Blnk,I=1,9),TRIM(Squish('JC=Node%Cdex('//IntToChar(Ell)//')'))
     WRITE(Out,*)(Blnk,I=1,9),'DO',Blnk,Squish('IQ=1,Node%NEll('//IntToChar(Ell)//')')
     WRITE(Out,*)(Blnk,I=1,12),'KQ=JQ+IQ'
     WRITE(Out,*)(Blnk,I=1,12),Squish('KC=JC+(IQ-1)*'//IntToChar(LHGTF(Ell)))
     WRITE(Out,*)(Blnk,I=1,12),'Z=Node%Zeta(KQ)'
     WRITE(Out,*)(Blnk,I=1,12),'LQx=Box%BndBox(1,1)-Node%Qx(KQ)'
     WRITE(Out,*)(Blnk,I=1,12),'LQy=Box%BndBox(2,1)-Node%Qy(KQ)'
     WRITE(Out,*)(Blnk,I=1,12),'LQz=Box%BndBox(3,1)-Node%Qz(KQ)'
     WRITE(Out,*)(Blnk,I=1,12),'UQx=Box%BndBox(1,2)-Node%Qx(KQ)'
     WRITE(Out,*)(Blnk,I=1,12),'UQy=Box%BndBox(2,2)-Node%Qy(KQ)'
     WRITE(Out,*)(Blnk,I=1,12),'UQz=Box%BndBox(3,2)-Node%Qz(KQ)'
     WRITE(Out,*)(Blnk,I=1,12),'LQ2=LQx*LQx+LQy*LQy+LQz*LQz'
     WRITE(Out,*)(Blnk,I=1,12),'UQ2=UQx*UQx+UQy*UQy+UQz*UQz'
     WRITE(Out,*)(Blnk,I=1,12),'TwoZ=Two*Z'
     WRITE(Out,*)(Blnk,I=1,12),'SqZ=SQRT(Z)'
     WRITE(Out,*)(Blnk,I=1,12),'CoFact=SqrtPi/(Two*SqZ)'
     WRITE(Out,*)(Blnk,I=1,12),'LLambdaX(0)=-CoFact*ERF(SqZ*LQx)'
     WRITE(Out,*)(Blnk,I=1,12),'LLambdaY(0)=-CoFact*ERF(SqZ*LQy)'
     WRITE(Out,*)(Blnk,I=1,12),'LLambdaZ(0)=-CoFact*ERF(SqZ*LQz)'
     WRITE(Out,*)(Blnk,I=1,12),'ULambdaX(0)=-CoFact*ERF(SqZ*UQx)'
     WRITE(Out,*)(Blnk,I=1,12),'ULambdaY(0)=-CoFact*ERF(SqZ*UQy)'
     WRITE(Out,*)(Blnk,I=1,12),'ULambdaZ(0)=-CoFact*ERF(SqZ*UQz)'

#else
     String(1,1)=Squish('CASE('//IntToChar(Ell)//')')
     WRITE(Out,*)(Blnk,I=1,6),TRIM(String(1,1))
#endif
     IF(Ell/=0)THEN
!        WRITE(Out,*)(Blnk,I=1,12),'LXptX=EXP(-Z*LQx*LQx)'
!        WRITE(Out,*)(Blnk,I=1,12),'UXptX=EXP(-Z*UQx*UQx)'
        WRITE(Out,*)(Blnk,I=1,12),'LXptX=EXPInv(Z*LQx*LQx)'
        WRITE(Out,*)(Blnk,I=1,12),'UXptX=EXPInv(Z*UQx*UQx)'
        WRITE(Out,*)(Blnk,I=1,12),'LLambdaX(1)=LXptX'
        WRITE(Out,*)(Blnk,I=1,12),'ULambdaX(1)=UXptX'
        WRITE(Out,*)(Blnk,I=1,12),'LLambdaX(2)=TwoZ*LQx*LXptX'
        WRITE(Out,*)(Blnk,I=1,12),'ULambdaX(2)=TwoZ*UQx*UXptX'
        DO L=3,Ell
           String(1,1)='LLambdaX('//IntToChar(L)                 &
                     //')=TwoZ*(LQx*LLambdaX('//IntToChar(L-1)   &
                     //')-'//DblToShrtChar(DBLE(L-2))           &
                     //'*LLambdaX('//IntToChar(L-2)//'))'
           String(2,1)='ULambdaX('//IntToChar(L)                 &
                     //')=TwoZ*(UQx*ULambdaX('//IntToChar(L-1)   &
                     //')-'//DblToShrtChar(DBLE(L-2))           &
                     //'*ULambdaX('//IntToChar(L-2)//'))'
           WRITE(Out,*)(Blnk,I=1,12),TRIM(Squish(String(1,1)))
           WRITE(Out,*)(Blnk,I=1,12),TRIM(Squish(String(2,1)))
        ENDDO

        WRITE(Out,*)(Blnk,I=1,12),'LXptY=EXPInv(Z*LQy*LQy)'
        WRITE(Out,*)(Blnk,I=1,12),'UXptY=EXPInv(Z*UQy*UQy)'
!        WRITE(Out,*)(Blnk,I=1,12),'LXptY=DEXP(-Z*LQy*LQy)'
!        WRITE(Out,*)(Blnk,I=1,12),'UXptY=DEXP(-Z*UQy*UQy)'
        WRITE(Out,*)(Blnk,I=1,12),'LLambdaY(1)=LXptY'
        WRITE(Out,*)(Blnk,I=1,12),'ULambdaY(1)=UXptY'
        WRITE(Out,*)(Blnk,I=1,12),'LLambdaY(2)=TwoZ*LQy*LXptY'
        WRITE(Out,*)(Blnk,I=1,12),'ULambdaY(2)=TwoZ*UQy*UXptY'
        DO L=3,Ell
           String(1,1)='LLambdaY('//IntToChar(L)                 &
                     //')=TwoZ*(LQy*LLambdaY('//IntToChar(L-1)   &
                     //')-'//DblToShrtChar(DBLE(L-2))            &
                     //'*LLambdaY('//IntToChar(L-2)//'))'
           String(2,1)='ULambdaY('//IntToChar(L)                 &
                     //')=TwoZ*(UQy*ULambdaY('//IntToChar(L-1)   &
                     //')-'//DblToShrtChar(DBLE(L-2))            &
                     //'*ULambdaY('//IntToChar(L-2)//'))'
           WRITE(Out,*)(Blnk,I=1,12),TRIM(Squish(String(1,1)))
           WRITE(Out,*)(Blnk,I=1,12),TRIM(Squish(String(2,1)))
        ENDDO
!        WRITE(Out,*)(Blnk,I=1,12),'LXptZ=DEXP(-Z*LQz*LQz)'
!        WRITE(Out,*)(Blnk,I=1,12),'UXptZ=DEXP(-Z*UQz*UQz)'
        WRITE(Out,*)(Blnk,I=1,12),'LXptZ=EXPInv(Z*LQz*LQz)'
        WRITE(Out,*)(Blnk,I=1,12),'UXptZ=EXPInv(Z*UQz*UQz)'
        WRITE(Out,*)(Blnk,I=1,12),'LLambdaZ(1)=LXptZ'
        WRITE(Out,*)(Blnk,I=1,12),'ULambdaZ(1)=UXptZ'
        WRITE(Out,*)(Blnk,I=1,12),'LLambdaZ(2)=TwoZ*LQz*LXptZ'
        WRITE(Out,*)(Blnk,I=1,12),'ULambdaZ(2)=TwoZ*UQz*UXptZ'
        DO L=3,Ell
           String(1,1)='LLambdaZ('//IntToChar(L)                 &
                     //')=TwoZ*(LQz*LLambdaZ('//IntToChar(L-1)   &
                     //')-'//DblToShrtChar(DBLE(L-2))            &
                     //'*LLambdaZ('//IntToChar(L-2)//'))'
           String(2,1)='ULambdaZ('//IntToChar(L)                 &
                     //')=TwoZ*(UQZ*ULambdaZ('//IntToChar(L-1)   &
                     //')-'//DblToShrtChar(DBLE(L-2))            &
                     //'*ULambdaZ('//IntToChar(L-2)//'))'
           WRITE(Out,*)(Blnk,I=1,12),TRIM(Squish(String(1,1)))
           WRITE(Out,*)(Blnk,I=1,12),TRIM(Squish(String(2,1)))
        ENDDO
     ENDIF
     DO L=0,Ell
        String(1,1)='LambdaX('//IntToChar(L)//')=LLambdaX('//IntToChar(L)//')-ULambdaX('//IntToChar(L)//')'
        String(2,1)='LambdaY('//IntToChar(L)//')=LLambdaY('//IntToChar(L)//')-ULambdaY('//IntToChar(L)//')'
        String(3,1)='LambdaZ('//IntToChar(L)//')=LLambdaZ('//IntToChar(L)//')-ULambdaZ('//IntToChar(L)//')'
        WRITE(Out,*)(Blnk,I=1,12),TRIM(Squish(String(1,1)))
        WRITE(Out,*)(Blnk,I=1,12),TRIM(Squish(String(2,1)))
        WRITE(Out,*)(Blnk,I=1,12),TRIM(Squish(String(3,1)))
     ENDDO
     K=0
     DO L=0,Ell
        DO M=0,Ell-L
           DO N=0,Ell-M-L
              LMN=LMNDex(L,M,N)
              K=K+1
              LMND(LMN)=K

#ifdef MULTIPLE_DIST
              String(1,LMN)='Pop=Pop+LambdaX('//IntToChar(L) &
                          //')*LambdaY('//IntToChar(M)//')*LambdaZ('//IntToChar(N)//')*Node%Co(KC+'//IntToChar(LMN)//')'

#else
              String(1,LMN)='Pop=Pop+LambdaX('//IntToChar(L) &
                          //')*LambdaY('//IntToChar(M)//')*LambdaZ('//IntToChar(N)//')*Node%Co('//IntToChar(LMN)//')'

#endif
            ENDDO
         ENDDO
     ENDDO
     DO J=1,LHGTF(Ell)
        K=LMND(J)
        WRITE(Out,*)(Blnk,I=1,12),TRIM(Squish(String(1,J)))
     ENDDO
#ifdef MULTIPLE_DIST
     WRITE(Out,*)(Blnk,I=1,9),'ENDDO'
#endif
   ENDDO
#ifdef MULTIPLE_DIST
#else
   WRITE(Out,*)(Blnk,I=1,6),'END SELECT'
#endif
  CLOSE(Out)


END PROGRAM
  
