!    Local QCTC thresholds for tuning MAC and PAC
!    Author: Matt Challacombe
!------------------------------------------------------------------------------
MODULE QCTCThresholds
   USE Derivedtypes
   USE GlobalScalars   
   USE GlobalObjects
   USE ProcessControl
   USE Indexing
   USE Parse
   USE InOut
   USE Macros
   USE Thresholding
   USE BoundingBox
   USE PoleNodeType
   IMPLICIT NONE
   REAL(DOUBLE)            :: TauPAC
   REAL(DOUBLE)            :: TauMAC
#ifdef NewPAC
   REAL(DOUBLE)                    :: GFactor
   REAL(DOUBLE),DIMENSION(0:SPEll) :: AACoef
#endif
   CONTAINS
      SUBROUTINE SetLocalThresholds(Tau)
         REAL(DOUBLE) :: Tau
!        Penetration Acceptability Criterion (PAC) threshold
         TauPAC=Tau*0.1D0
!        Multipole Acceptability Criterion (MAC) threshold
         TauMAC=Tau
      END SUBROUTINE SetLocalThresholds
#ifdef NewPAC
!======================================================================================
!
!======================================================================================
      SUBROUTINE SetLocalCoefs(LocalG)
        REAL(DOUBLE)  :: LocalG,X1,X2,X3,F1,F2,F3
        INTEGER       :: L,I
!       Set Up GFactor         
        IF(LocalG .GE. One) THEN
           CALL Halt('GFactor must be less then One')
        ENDIF
        GFactor  = LocalG
!       Solve for AACoef
        AACoef(0)= 1.D0
        DO L=1,BFEll
           X1 = SQRT(0.499D0*DBLE(L)/(1.D0-GFactor))
           X2 = 2.D0*X1
           F1 = Herm(X1,L+1)+2.D0*GFactor*X1*Herm(X1,L)
           F2 = Herm(X2,L+1)+2.D0*GFactor*X2*Herm(X2,L)
           IF(SIGN(1.D0,F1)==SIGN(1.D0,F2)) THEN
              CALL Halt('Failed to Find Root in  SetLocalCoefs')
           ENDIF
           DO I=1,100
              X3 = 0.5D0*(X1+X2)
              F3 = Herm(X3,L+1)+2.D0*GFactor*X3*Herm(X3,L)
              IF(SIGN(1.D0,F3)==SIGN(1.D0,F1)) THEN
                 X1 = X3
              ELSE
                 X2 = X3
              ENDIF
           ENDDO
           AACoef(L) = ABS(Herm(X3,L))*EXP(-(One-GFactor)*X3*X3)
!           WRITE(*,*) 'L=',L,' X3 = ',X3,' F3 = ',F3
!           WRITE(*,*) 'ACoef = ',AACoef(L)
        ENDDO
!
      END SUBROUTINE SetLocalCoefs
!======================================================================================
!
!======================================================================================
      FUNCTION MaxCoef(Ell,Zeta,HGCo) 
!
        INTEGER                       :: Ell,L,M,N,LMN,IJK
        REAL(DOUBLE)                  :: MaxCoef,BFac,Zeta
        REAL(DOUBLE), DIMENSION(1:)   :: HGCo
!
        MaxCoef = Zero
        IF(Ell==0) THEN
           MaxCoef = ABS(HGCo(1))
           RETURN
        END IF
        DO L=0,Ell
           DO M=0,Ell-L
              DO N=0,Ell-L-M
                 LMN = LMNDex(L,M,N)
                 IJK = L+M+N
                 IF(IJK==0) THEN
                    BFac = One
                 ELSE
                    BFac = (Zeta)**(Half*DBLE(IJK))
                 ENDIF
                 MaxCoef = MaxCoef+ABS(HGCo(LMN))*BFac*AACoef(L)*AACoef(M)*AACoef(N)
              ENDDO
           ENDDO
        ENDDO
!
      END FUNCTION MaxCoef
!======================================================================================
!
!======================================================================================
      FUNCTION Erfcc(x)
        REAL(DOUBLE)                :: Erfcc,x,y
        REAL(DOUBLE),PARAMETER      :: A1= -1.1283791670955125739D0
        REAL(DOUBLE),PARAMETER      :: A2= -0.6366197723675813430D0
!
        y     = A1*x+A2*x*x
        Erfcc = EXP(y)
!
      END FUNCTION Erfcc
!======================================================================================
!
!======================================================================================
      FUNCTION Herm(x,L)
        INTEGER                     :: L,I,J
        REAL(DOUBLE)                :: Herm,x
        REAL(DOUBLE),DIMENSION(0:9) :: CC
!
        SELECT CASE(L)
        CASE(0)
           CC(0) = 1.D0
        CASE(1)
           CC(0) =  0.D0
           CC(1) = -2.D0
        CASE(2)
           CC(0) = -2.D0
           CC(1) =  0.D0
           CC(2) =  4.D0
        CASE(3)
           CC(0) =  0.D0
           CC(1) = 12.D0
           CC(2) =  0.D0
           CC(3) = -8.D0
        CASE(4)
           CC(0) =  12.D0
           CC(1) =   0.D0
           CC(2) = -48.D0
           CC(3) =   0.D0
           CC(4) =  16.D0
        CASE(5)
           CC(0) =    0.D0
           CC(1) = -120.D0
           CC(2) =    0.D0
           CC(3) =  160.D0
           CC(4) =    0.D0
           CC(5) =  -32.D0
        CASE(6)
           CC(0) = -120.D0
           CC(1) =    0.D0
           CC(2) =  720.D0
           CC(3) =    0.D0
           CC(4) = -480.D0
           CC(5) =    0.D0
           CC(6) =   64.D0
        CASE(7)
           CC(0) =     0.D0
           CC(1) =  1680.D0
           CC(2) =     0.D0
           CC(3) = -3360.D0
           CC(4) =     0.D0
           CC(5) =  1344.D0
           CC(6) =     0.D0
           CC(7) =  -128.D0
        CASE(8)
           CC(0) =   1680.D0
           CC(1) =      0.D0
           CC(2) = -13440.D0
           CC(3) =      0.D0
           CC(4) =  13440.D0
           CC(5) =      0.D0
           CC(6) =  -3584.D0
           CC(7) =      0.D0
           CC(8) =    256.D0
        CASE(9)
           CC(0) =      0.D0
           CC(1) = -30240.D0 
           CC(2) =      0.D0
           CC(3) =  80640.D0
           CC(4) =      0.D0
           CC(5) = -48384.D0
           CC(6) =      0.D0
           CC(7) =   9216.D0
           CC(8) =      0.D0
           CC(9) =   -512.D0
         CASE(10:)
            CALL Halt('Hermite Function not defined in SetLocalCoefs')
        END SELECT
        Herm = CC(L)
        DO I=0,L-1
           J = L-1-I
           Herm = Herm*x+CC(J)
        ENDDO
      END FUNCTION Herm
!
#endif
END MODULE
