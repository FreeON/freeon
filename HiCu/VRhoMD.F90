MODULE WriteRoutines
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE ProcessControl
  USE InOut
  USE Parse
  CONTAINS
     SUBROUTINE WriteErf(Unit,NBlnks,Res,Arg,Mult)
        INTEGER                 :: Unit,NBlnks
        CHARACTER(LEN=NBlnks)   :: Tab1
        CHARACTER(LEN=NBlnks+3) :: Tab2
        CHARACTER(LEN=*)        :: Res,Arg,Mult
        CHARACTER(LEN=DEFAULT_CHR_LEN) :: BigStr
        Tab1='';DO I=1,NBlnks; Tab1=Tab1//Blnk; ENDDO
        Tab2=Tab1//Blnk//Blnk//Blnk
        WRITE(Unit,*)Tab1,'W=',Arg
        WRITE(Unit,*)Tab1,'Sgn=One; IF(W<0.0D0)Sgn=-One'
        WRITE(Unit,*)Tab1,'X=Sgn*W'
        WRITE(Unit,*)Tab1,'IF(X>Erf_Switch)THEN'
        WRITE(Unit,*)Tab2,TRIM(Res),'=',TRIM(Mult),'*Sgn'
        WRITE(Unit,*)Tab1,'ELSE'
        WRITE(Unit,*)Tab2,'J=AINT(X*Erf_Grid)'
        BigStr=TRIM(Res)//'='//TRIM(Mult)//'*Sgn*(Erf_0(J)+X*(Erf_1(J)+X*(Erf_2(J)+X*Erf_3(J))))' 
        BigStr=Tab2//BigStr
        WRITE(Unit,*)TRIM(BigStr)
        WRITE(Unit,*)Tab1,'ENDIF'
     END SUBROUTINE
     SUBROUTINE WriteExpInv(Unit,NBlnks,Res,Arg)
        INTEGER                 :: Unit,NBlnks
        CHARACTER(LEN=NBlnks)   :: Tab1
        CHARACTER(LEN=NBlnks+3) :: Tab2
        CHARACTER(LEN=*)        :: Res,Arg
        CHARACTER(LEN=DEFAULT_CHR_LEN) :: BigStr
!        WRITE(Unit,*)TRIM(RES),'=EXPInv(',TRIM(Arg),')'
!        RETURN
        Tab1='';DO I=1,NBlnks; Tab1=Tab1//Blnk; ENDDO
        Tab2=Tab1//Blnk//Blnk//Blnk
        WRITE(Unit,*)Tab1,'X='//Arg
        WRITE(Unit,*)Tab1,'IF(X>Exp_Switch)THEN'
        WRITE(Unit,*)Tab2//Res//'=Zero'
        WRITE(Unit,*)Tab1,'ELSE'
        WRITE(Unit,*)Tab2,'J=AINT(X*Exp_Grid)'
        WRITE(Unit,*)Tab2//Res//'=Exp_0(J)+X*(Exp_1(J)+X*(Exp_2(J)+X*(Exp_3(J)+X*Exp_4(J))))'
        WRITE(Unit,*)Tab1,'ENDIF'
     END SUBROUTINE 
END MODULE
PROGRAM VectorRhoMD
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE ProcessControl
  USE InOut
  USE Parse
  USE Indexing
  USE WriteRoutines
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
  WRITE(Out,*)(Blnk,I=1,6),'Z=Prim%Z'
  WRITE(Out,*)(Blnk,I=1,6),'TwoZ=Two*Z'
  WRITE(Out,*)(Blnk,I=1,6),'SELECT CASE(Prim%Ell)'
  DO Ell=0,MaxEll
     ChEll   =IntToChar(Ell)
     ChLMNLen=IntToChar(LHGTF(Ell))
     String(1,1)=Squish('CASE('//IntToChar(Ell)//')')
     WRITE(Out,*)(Blnk,I=1,6),TRIM(String(1,1))
     WRITE(Out,*)(Blnk,I=1,9),'DO I=1,NGrid '
     WRITE(Out,*)(Blnk,I=1,12),'RPx=Cube%Grid(I,1)-Prim%P(1)'
     WRITE(Out,*)(Blnk,I=1,12),'RPy=Cube%Grid(I,2)-Prim%P(2)'
     WRITE(Out,*)(Blnk,I=1,12),'RPz=Cube%Grid(I,3)-Prim%P(3)'
     WRITE(Out,*)(Blnk,I=1,12),'RP2=RPx**2+RPy**2+RPz**2'
     WRITE(Out,*)(Blnk,I=1,12),'X=Z*RP2'
     WRITE(Out,*)(Blnk,I=1,12),'IF(X<Exp_Switch)THEN'
     WRITE(Out,*)(Blnk,I=1,15),'J=AINT(X*Exp_Grid)'
     WRITE(Out,*)(Blnk,I=1,15),'Xpt=Exp_0(J)+X*(Exp_1(J)+X*(Exp_2(J)+X*(Exp_3(J)+X*Exp_4(J))))'
     WRITE(Out,*)(Blnk,I=1,15),'LambdaX(0)=Xpt'
     WRITE(Out,*)(Blnk,I=1,15),'LambdaY(0)=One'
     WRITE(Out,*)(Blnk,I=1,15),'LambdaZ(0)=One'
     WRITE(Out,*)(Blnk,I=1,15),'LambdaX(1)=TwoZ*RPx*Xpt'
     WRITE(Out,*)(Blnk,I=1,15),'LambdaY(1)=TwoZ*RPy'
     WRITE(Out,*)(Blnk,I=1,15),'LambdaZ(1)=TwoZ*RPz'
     IF(Ell>=1)THEN
        DO K=1,3
           DO L=2,Ell+1
              String(1,1)='Lambda'//Carts(K)//'('//IntToChar(L)           &
                        //')=TwoZ*(RP'//Carts(K)                          &
                        //'*Lambda'//Carts(K)                             &
                        //'('//IntToChar(L-1)//')-'//Floats(L-1)          &
                        //'Lambda'//Carts(K)//'('//IntToChar(L-2)//'))'
              WRITE(Out,*)(Blnk,I=1,15),TRIM(Squish(String(1,1)))
           ENDDO
        ENDDO
     ENDIF
     WRITE(Out,*)(Blnk,I=1,15),'Wght=Cube%Wght(I)'
     WRITE(Out,*)(Blnk,I=1,15),'dEdRho=Cube%Vals(I,1)'
     WRITE(Out,*)(Blnk,I=1,15),'dEdAbsGradRho2=Cube%Vals(I,2)'
     WRITE(Out,*)(Blnk,I=1,15),'GradRhoX=Cube%Vals(I,3)'
     WRITE(Out,*)(Blnk,I=1,15),'GradRhoY=Cube%Vals(I,4)'
     WRITE(Out,*)(Blnk,I=1,15),'GradRhoZ=Cube%Vals(I,5)'
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
        WRITE(Out,*)(Blnk,I=1,15),TRIM(Squish(String(1,K)))
        WRITE(Out,*)(Blnk,I=1,15),TRIM(Squish(String(2,K)))
        WRITE(Out,*)(Blnk,I=1,15),TRIM(Squish(String(3,K)))
        WRITE(Out,*)(Blnk,I=1,15),TRIM(Squish(String(4,K)))
        WRITE(Out,*)(Blnk,I=1,15),TRIM(Squish(String(5,K)))
        WRITE(Out,*)(Blnk,I=1,15),TRIM(Squish(String(6,K)))
        WRITE(Out,*)(Blnk,I=1,15),TRIM(Squish(String(7,K)))
        WRITE(Out,*)(Blnk,I=1,15),TRIM(Squish(String(8,K)))

             ENDDO
         ENDDO
     ENDDO
     WRITE(Out,*)(Blnk,I=1,12),'ENDIF'
     WRITE(Out,*)(Blnk,I=1,9),'ENDDO'
  ENDDO
  String(1,1)=Squish('CASE('//IntToChar(MaxEll+1)//':)')
  WRITE(Out,*)(Blnk,I=1,6),TRIM(String(1,1))
  WRITE(Out,*)(Blnk,I=1,9),'CALL Halt("Increase MaxEll in MondoMods/GlobalScalars.F90 and remake")'
  WRITE(Out,*)(Blnk,I=1,6),'END SELECT'
  CLOSE(Out)
!=================================================================================
!
!==================================================================================
  CALL OpenASCII('ExplicitLeafContribution.Inc',Out,.TRUE.)
  WRITE(Out,*)(Blnk,I=1,6),'SELECT CASE(Node%Ell)'
  DO Ell=0,MaxEll
     ChEll   =IntToChar(Ell)
     ChLMNLen=IntToChar(LHGTF(Ell))
     String(1,1)=Squish('CASE('//IntToChar(Ell)//')')
     WRITE(Out,*)(Blnk,I=1,6),TRIM(String(1,1))
     WRITE(Out,*)(Blnk,I=1,9),'DO IGrid=1,NGrid'
     WRITE(Out,*)(Blnk,I=1,12),'RQx=Grid(IGrid,1)-Node%Qx'
     WRITE(Out,*)(Blnk,I=1,12),'RQy=Grid(IGrid,2)-Node%Qy'
     WRITE(Out,*)(Blnk,I=1,12),'RQz=Grid(IGrid,3)-Node%Qz'
     WRITE(Out,*)(Blnk,I=1,12),'RQ2=RQx*RQx+RQy*RQy+RQz*RQz'
     WRITE(Out,*)(Blnk,I=1,12),'X=Node%Zeta*RQ2'
     WRITE(Out,*)(Blnk,I=1,12),'IF(X<Exp_Switch)THEN'
     WRITE(Out,*)(Blnk,I=1,15),'J=AINT(X*Exp_Grid)'
     WRITE(Out,*)(Blnk,I=1,15),'Xpt=Exp_0(J)+X*(Exp_1(J)+X*(Exp_2(J)+X*(Exp_3(J)+X*Exp_4(J))))'
     WRITE(Out,*)(Blnk,I=1,15),'TwoZ=Two*Node%Zeta'
     WRITE(Out,*)(Blnk,I=1,15),'LambdaX(0)=Xpt'
     WRITE(Out,*)(Blnk,I=1,15),'LambdaX(1)=TwoZ*RQx*Xpt'
     WRITE(Out,*)(Blnk,I=1,15),'LambdaY(0)=One'
     WRITE(Out,*)(Blnk,I=1,15),'LambdaY(1)=TwoZ*RQy'
     WRITE(Out,*)(Blnk,I=1,15),'LambdaZ(0)=One'
     WRITE(Out,*)(Blnk,I=1,15),'LambdaZ(1)=TwoZ*RQz'
     IF(Ell>=1)THEN
        DO K=1,3
           DO L=2,Ell+1
              String(1,1)='Lambda'//Carts(K)//'('//IntToChar(L)           &
                        //')=TwoZ*(RQ'//Carts(K)//'*Lambda'//Carts(K)     &
                        //'('//IntToChar(L-1)//')-'//Floats(L-1)          &
                        //'Lambda'//Carts(K)//'('//IntToChar(L-2)//'))'
               WRITE(Out,*)(Blnk,I=1,15),TRIM(Squish(String(1,1)))
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
              String(1,LMN)='Co=Node%Co('//IntToChar(LMN)//')'
              String(2,LMN)='RhoV(IGrid,1)=RhoV(IGrid,1)+LambdaX('//IntToChar(L) &
                          //')*LambdaY('//IntToChar(M)//')*LambdaZ('//IntToChar(N)//')*Co'
              String(3,LMN)='RhoV(IGrid,2)=RhoV(IGrid,2)-LambdaX('//IntToChar(L+1) &
                          //')*LambdaY('//IntToChar(M)//')*LambdaZ('//IntToChar(N)//')*Co'
              String(4,LMN)='RhoV(IGrid,3)=RhoV(IGrid,3)-LambdaX('//IntToChar(L) &
                          //')*LambdaY('//IntToChar(M+1)//')*LambdaZ('//IntToChar(N)//')*Co'
              String(5,LMN)='RhoV(IGrid,4)=RhoV(IGrid,4)-LambdaX('//IntToChar(L) &
                          //')*LambdaY('//IntToChar(M)//')*LambdaZ('//IntToChar(N+1)//')*Co'
            ENDDO
         ENDDO
     ENDDO
     DO J=1,LHGTF(Ell)
        K=LMND(J)
        WRITE(Out,*)(Blnk,I=1,15),TRIM(Squish(String(1,J)))
        WRITE(Out,*)(Blnk,I=1,15),TRIM(Squish(String(2,J)))
        WRITE(Out,*)(Blnk,I=1,15),TRIM(Squish(String(3,J)))
        WRITE(Out,*)(Blnk,I=1,15),TRIM(Squish(String(4,J)))
        WRITE(Out,*)(Blnk,I=1,15),TRIM(Squish(String(5,J)))
     ENDDO
     WRITE(OUT,*)(Blnk,I=1,12),'ENDIF'
     WRITE(OUT,*)(Blnk,I=1,9),'ENDDO'
  ENDDO
  WRITE(Out,*)(Blnk,I=1,6),'END SELECT'

  CLOSE(Out)
!==================================================================================
!
!==================================================================================
  CALL OpenASCII('ExplicitLeafPopulation.Inc',Out,.TRUE.)
   WRITE(Out,*)(Blnk,I=1,6),'Z=Node%Zeta'
   WRITE(Out,*)(Blnk,I=1,6),'SqZ=SQRT(Z)'
   WRITE(Out,*)(Blnk,I=1,6),'TwoZ=Two*Z'
   WRITE(Out,*)(Blnk,I=1,6),'CoFact=SqrtPi/(Two*SqZ)'
!
   WRITE(Out,*)(Blnk,I=1,6),'LQx=Box%BndBox(1,1)-Node%Qx'
   WRITE(Out,*)(Blnk,I=1,6),'UQx=Box%BndBox(1,2)-Node%Qx'
   CALL WriteErf(Out,6,'LLambdaX(0)','SqZ*LQx','-CoFact')
   CALL WriteErf(Out,6,'ULambdaX(0)','SqZ*UQx','-CoFact')
   WRITE(Out,*)(Blnk,I=1,6),'LQy=Box%BndBox(2,1)-Node%Qy'
   WRITE(Out,*)(Blnk,I=1,6),'UQy=Box%BndBox(2,2)-Node%Qy'
   CALL WriteErf(Out,6,'LLambdaY(0)','SqZ*LQy','-CoFact')
   CALL WriteErf(Out,6,'ULambdaY(0)','SqZ*UQy','-CoFact')
   WRITE(Out,*)(Blnk,I=1,6),'LQz=Box%BndBox(3,1)-Node%Qz'
   WRITE(Out,*)(Blnk,I=1,6),'UQz=Box%BndBox(3,2)-Node%Qz'
   CALL WriteErf(Out,6,'LLambdaZ(0)','SqZ*LQz','-CoFact')
   CALL WriteErf(Out,6,'ULambdaZ(0)','SqZ*UQz','-CoFact')
   WRITE(Out,*)(Blnk,I=1,6),'LQ2=LQx*LQx+LQy*LQy+LQz*LQz'
   WRITE(Out,*)(Blnk,I=1,6),'UQ2=UQx*UQx+UQy*UQy+UQz*UQz'
   WRITE(Out,*)(Blnk,I=1,6),'SELECT CASE(Node%Ell)'
   DO Ell=0,MaxEll
      ChEll   =IntToChar(Ell)
      ChLMNLen=IntToChar(LHGTF(Ell))
      String(1,1)=Squish('CASE('//IntToChar(Ell)//')')
      WRITE(Out,*)(Blnk,I=1,6),TRIM(String(1,1))
      IF(Ell/=0)THEN
        DO K=1,3
           CALL WriteExpInv(Out,9,'LXpt'//Carts(K),'Z*(LQ'//Carts(K)//'**2)')
           CALL WriteExpInv(Out,9,'UXpt'//Carts(K),'Z*(UQ'//Carts(K)//'**2)')
           WRITE(Out,*)(Blnk,I=1,9),'LLambda'//Carts(K)//'(1)=LXpt'//Carts(K)
           WRITE(Out,*)(Blnk,I=1,9),'ULambda'//Carts(K)//'(1)=UXpt'//Carts(K)
           WRITE(Out,*)(Blnk,I=1,9),'LLambda'//Carts(K)//'(2)=TwoZ*LQ'//Carts(K)//'*LXpt'//Carts(K)
           WRITE(Out,*)(Blnk,I=1,9),'ULambda'//Carts(K)//'(2)=TwoZ*UQ'//Carts(K)//'*UXpt'//Carts(K)
        ENDDO
        DO L=3,Ell
           DO K=1,3
              String(1,1)='LLambda'//Carts(K)//'('//IntToChar(L)    &
                        //')=TwoZ*(LQ'//Carts(K)//'*LLambda'//Carts(K)//'('//IntToChar(L-1)   &
                        //')-'//Floats(L-2)//'LLambda'//Carts(K)//'('//IntToChar(L-2)//'))'
              String(2,1)='ULambda'//Carts(K)//'('//IntToChar(L)    &
                        //')=TwoZ*(UQ'//Carts(K)//'*ULambda'//Carts(K)//'('//IntToChar(L-1)   &
                        //')-'//Floats(L-2)//'ULambda'//Carts(K)//'('//IntToChar(L-2)//'))'
               WRITE(Out,*)(Blnk,I=1,9),TRIM(Squish(String(1,1)))
               WRITE(Out,*)(Blnk,I=1,9),TRIM(Squish(String(2,1)))
           ENDDO
          ENDDO
        ENDIF
     DO L=0,Ell
        String(1,1)='LambdaX('//IntToChar(L)//')=LLambdaX('//IntToChar(L)//')-ULambdaX('//IntToChar(L)//')'
        String(2,1)='LambdaY('//IntToChar(L)//')=LLambdaY('//IntToChar(L)//')-ULambdaY('//IntToChar(L)//')'
        String(3,1)='LambdaZ('//IntToChar(L)//')=LLambdaZ('//IntToChar(L)//')-ULambdaZ('//IntToChar(L)//')'
        WRITE(Out,*)(Blnk,I=1,9),TRIM(Squish(String(1,1)))
        WRITE(Out,*)(Blnk,I=1,9),TRIM(Squish(String(2,1)))
        WRITE(Out,*)(Blnk,I=1,9),TRIM(Squish(String(3,1)))
     ENDDO
     K=0
     DO L=0,Ell
        DO M=0,Ell-L
           DO N=0,Ell-M-L
              LMN=LMNDex(L,M,N)
              K=K+1
              LMND(LMN)=K
              String(1,LMN)='Pop=Pop+LambdaX('//IntToChar(L) &
                          //')*LambdaY('//IntToChar(M)//')*LambdaZ('//IntToChar(N) &
                          //')*Node%Co('//IntToChar(LMN)//')'

            ENDDO
         ENDDO
     ENDDO
     DO J=1,LHGTF(Ell)
        K=LMND(J)
        WRITE(Out,*)(Blnk,I=1,9),TRIM(Squish(String(1,J)))
     ENDDO
   ENDDO
   WRITE(Out,*)(Blnk,I=1,6),'END SELECT'
  CLOSE(Out)


END PROGRAM
  
