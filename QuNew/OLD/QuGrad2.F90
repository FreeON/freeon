PROGRAM QuGrad
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE SetXYZ
  USE LinAlg
  USE Functionals
  USE MatFunk
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
  TYPE(ARGMT)                :: Args
  CHARACTER(LEN=6),PARAMETER :: Prog='QuGrad'
  TYPE(BSet)                 :: BS
  TYPE(CRDS)                 :: GM
  TYPE(DBL_VECT)             :: G,X,Co
  INTEGER                    :: N3,I,J,K,IGeom,I1,I2,I3,I4,L,M,NOxygens
  REAL(DOUBLE)               :: A,F,R,E,dRdK,dRdJ,dRdK2,Cross, &
                                DeOO,BOO,ReOO,DeOHb,DeOHnb,BOHb,BOHnb,ReOHb,ReOHnb, &
                                DeHHb,DeHHnb,BHHb,BHHnb,ReHHb,ReHHnb
  REAL(DOUBLE),DIMENSION(6,6)   :: De,B,Re
  TYPE(DBL_RNK2)              :: H
!----------------------------------------------------------------------------

  CALL StartUp(Args,Prog)
  IGeom=Args%I%I(3)
!  Override default blocking...
  N3=3*NAtoms
  NBasF=N3
  PrintFlags%Mat=DEBUG_MATRICES
  BSiz%I=3
  OffS%I(1)=1
  DO I=2,NAtoms
     OffS%I(I)=OffS%I(I-1)+3
  ENDDO
!  Allocations
  CALL New(G,N3)
  CALL New(X,N3)
  CALL New(H,(/N3,N3/))
  CALL Get(GM,Tag_O=CurGeom)

  NOxygens=2
  DeOO=-0.01d0
  DeOHb=0.16d0
  DeOHnb=0.08d0
  DeHHb=0.05d0
  DeHHnb=-0.005d0
  BOO=1.0d0
  BOHb=2.0d0
  BOHnb=0.25d0
  BHHb=0.5d0
  BHHnb=1.0d0
  ReOO=0.0d0*AngstromsToAU
  ReOHb=1.0d0*AngstromsToAU
  ReOHnb=3.0d0*AngstromsToAU
  ReHHb=1.6d0*AngstromsToAU
  ReHHnb=0.0d0*AngstromsToAU

  K=0
  DO I=1,NAtoms
     DO J=1,3
        K=K+1
        X%D(K)=GM%Carts%D(J,I)
     ENDDO
     DO J=I,NAtoms
        IF (I.EQ.J) THEN
           De(I,J)=Zero
           B(I,J)=Zero
           Re(I,J)=Zero
        ELSEIF (I.LE.NOxygens) THEN
           IF (J.LE.NOxygens) THEN
              De(I,J)=DeOO
              B(I,J)=BOO
              Re(I,J)=ReOO
           ELSEIF (J.EQ.2*I+NOxygens-1.OR.J.EQ.2*I+NOxygens)THEN
              De(I,J)=DeOHb
              B(I,J)=BOHb
              Re(I,J)=ReOHb
           ELSE
              De(I,J)=DeOHnb
              B(I,J)=BOHnb
              Re(I,J)=ReOHnb
           ENDIF
        ELSE
           IF (J.EQ.I+1) THEN
              IF (MOD(NAtoms-I,2).EQ.0) THEN
                 De(I,J)=DeHHnb
                 B(I,J)=BHHnb
                 Re(I,J)=ReHHnb
              ELSE
                 De(I,J)=DeHHb
                 B(I,J)=BHHb
                 Re(I,J)=ReHHb
              ENDIF
           ELSE
              De(I,J)=DeHHnb
              B(I,J)=BHHnb
              Re(I,J)=ReHHnb
           ENDIF
        ENDIF
        De(J,I)=De(I,J)
        B(J,I)=B(I,J)
        Re(J,I)=Re(I,J)
     ENDDO        
  ENDDO

  E=Zero
  G%D=Zero
  H%D=Zero
  DO K=1,NAtoms
     DO J=1,NAtoms
        IF (K.NE.J) then
           I1=3*(K-1)
           I2=3*(J-1)
           
           R=SQRT((X%D(I1+1)-X%D(I2+1))**2+(X%D(I1+2)-X%D(I2+2))**2+(X%D(I1+3)-X%D(I2+3))**2)

!     Harmonic Oscillator
!               E=E+0.5d0*fc(k,j)*(r-re(k,j))**2
!     Morse Pot.         
           E=E+Half*De(K,J)*(One - EXP(-B(K,J)*(R-Re(K,J))))**2
           DO L=1,3
              I1=3*(K-1)+L
              I2=3*(J-1)+L
              dRdK=(X%D(I1)-X%D(I2))/R
!     Harmonic Oscillator
!     Ep(i1)=Ep(i1)+2.0d0*fc(k,j)*(r-re(k,j))*
!     +                 drdk
!     Morse Pot.
              G%D(I1)=G%D(I1)+Two*De(K,J)*B(K,J)*(One-EXP(-B(K,J)*(R-Re(K,J))))* &
                   EXP(-B(K,J)*(R-Re(K,J)))*dRdK
!     Compute the Hessian
              DO M=1,3
!     Morse Pot.
                 I3=3*(J-1)+M
                 I4=3*(K-1)+M
                 dRdJ=(X%D(I3)-X%D(I4))/R
                 dRdK2=(X%D(I4)-X%D(I3))/R
                 Cross=0.0d0
                 IF (L.eq.M) THEN
                    Cross=-Two*De(K,J)*B(K,J)*(One/R)*EXP(-B(K,J)*(R-Re(K,J)))* &
                         ((One-EXP(-B(K,J)*(R-Re(K,J)))))
                 ENDIF

                 H%D(I1,I3)=H%D(I1,I3)+Two*De(K,J)*B(K,J)*dRdK*dRdJ*EXP(-B(K,J)*(R-Re(K,J)))*&
                      (Two*B(K,J)*EXP(-B(K,J)*(R-Re(K,J)))-(One-EXP(-B(K,J)*(R-Re(K,J))))/ &
                      R - B(K,J)) + Cross

                 H%D(I1,I4)=H%D(I1,I4)+Two*De(K,J)*B(K,J)*dRdK*dRdK2*EXP(-B(K,J)*(R-Re(K,J)))*&
                      (Two*B(K,J)*EXP(-B(K,J)*(R-Re(K,J)))-(One-EXP(-B(K,J)*(R-Re(K,J))))/ &
                      R - B(K,J)) - Cross
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO

  CALL PPrint(X,'X['//TRIM(IntToChar(IGeom))//']',Unit_O=6)
!
  WRITE(*,*) 'TOTAL ENERGY', E
  CALL PPrint(G,'G['//TRIM(IntToChar(IGeom))//']',Unit_O=6)
!  CALL PPrint(H,'H['//TRIM(IntToChar(IGeom))//']',Unit_O=6)
  CALL Put(G,'GradE',Tag_O=CurGeom)
  CALL Put(H,'MorseH',Tag_O=CurGeom)
  
  CALL Shutdown(Prog)
END PROGRAM QuGrad





