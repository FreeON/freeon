!==============================================================================
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe, 1997 
!
!--  Matt Challacombe and  C. J. Tymczak
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!    MATTS MULTIPOLE MANIPULATION LIBRARY IN F90
!==============================================================================
MODULE MMML
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  USE Moments
  IMPLICIT NONE
!---------------------------------------------------------------------------
! Global variables
!---------------------------------------------------------------------------
!
  TYPE(DBL_VECT)        :: Factorial
  TYPE(DBL_VECT)        :: FactOlm0
  TYPE(DBL_VECT)        :: FactOlm1
  TYPE(DBL_RNK2)        :: FactOlm2
  TYPE(DBL_VECT)        :: FactMlm0
  TYPE(DBL_VECT)        :: FactMlm1
  TYPE(DBL_RNK2)        :: FactMlm2
!
  TYPE(DBL_VECT)        :: Cosine
  TYPE(DBL_VECT)        :: Sine
  TYPE(DBL_VECT)        :: RToTh
  TYPE(DBL_RNK2)        :: LegendreP
!
CONTAINS
!---------------------------------------------------------------------------
! Set Up the Factorials for the Multipoles
!---------------------------------------------------------------------------
  SUBROUTINE FactorialSetUpF90(MaxL)
    INTEGER                           :: MaxL
    INTEGER                           :: L,M
    REAL(DOUBLE)                      :: Sgn,DblFact,TwoTimes
!
    CALL New(Factorial,2*MaxL,0)
!
    CALL New(FactOlm0 ,MaxL  ,0)
    CALL New(FactOlm1 ,MaxL  ,0)
    CALL New(FactMlm0 ,MaxL  ,0)
    CALL New(FactMlm1 ,MaxL  ,0)
!
    CALL New(FactOlm2 ,(/MaxL,MaxL/),(/0,0/))
    CALL New(FactMlm2 ,(/MaxL,MaxL/),(/0,0/))
!
    CALL New(Cosine   ,MaxL,0)
    CALL New(Sine     ,MaxL,0)
    CALL New(RToTh    ,MaxL,0)
    CALL New(LegendreP,(/MaxL,MaxL/),(/0,0/))
!      
    Factorial%D(0)=One
    DO L=1,2*MaxL
       Factorial%D(L)=Factorial%D(L-1)*DBLE(L)         
    ENDDO
!
    DO L=0,MaxL
       FactOlm1%D(L)=DBLE(2*L-1)
       FactMlm1%D(L)=DBLE(2*L-1)
    ENDDO
!
    DO L=2,MaxL
       DO M=0,L-2
          FactOlm2%D(M,L)=1.0D0/DBLE((L+M)*(L-M))
       ENDDO
    ENDDO
!
    DO L=0,MaxL
       DO M=0,MaxL
          FactMlm2%D(M,L)=DBLE((L+M-1)*(L-M-1))
       ENDDO
    ENDDO
!
    Sgn      = One
    DblFact  = One
    TwoTimes = One
    DO L=0,MaxL
       FactOlm0%D(L)=Sgn*DblFact/Factorial%D(2*L)
       FactMlm0%D(L)=Sgn*DblFact
       DblFact=DblFact*TwoTimes
       TwoTimes=TwoTimes+Two
       Sgn=-Sgn
    ENDDO
!
  END SUBROUTINE FactorialSetUpF90
!---------------------------------------------------------------------------------------
! The IrRegular Multipoles
!---------------------------------------------------------------------------------------
  SUBROUTINE IrRegularF90(Ell,LenEll,x,y,z,Cpq,Spq)
    INTEGER                          :: Ell,LenEll
    INTEGER                          :: L,M,LM,L1,L2
    REAL(DOUBLE)                     :: x,y,z
    REAL(DOUBLE)                     :: x2,y2,R,OneOverR,Rxy,OneOvRxy,Rl,CoTan, &
                                        TwoC,Sq,rs,Fact,o
    REAL(DOUBLE),DIMENSION(0:LenEll) :: Cpq,Spq
!
    x2=x*x
    y2=y*y      
    R=SQRT(x2+y2+z*z)
    OneOverR=One/R
!
    Cosine%D(0)=One
    Sine%D(  0)=Zero
!
    Rxy=DSQRT(x2+y2)
    IF(Rxy .NE. Zero)THEN
       OneOvRxy=One/DSQRT(x2+y2)       
       Cosine%D(1)=y*OneOvRxy
       Sine%D(  1)=x*OneOvRxy
    ELSE
       Cosine%D(1)=0.70710678118654752D0
       Sine%D(  1)=0.70710678118654752D0         
    ENDIF
    Rl=OneOverR
    CoTan=z*Rl
!
    TwoC=Two*Cosine%D(1)
    DO L=2,Ell
       L1=L-1
       L2=L-2
       Cosine%D(L)=TwoC*Cosine%D(L1)-Cosine%D(L2)
       Sine%D(  L)=TwoC*Sine%D(  L1)-Sine%D(  L2)
    ENDDO
!             
    DO L=0,Ell
       RToTh%D(L)=Rl
       Rl=Rl*OneOverR
    ENDDO
!
    Sq=SQRT(One-CoTan*CoTan)
    rs=One
    DO L=0,Ell
       LegendreP%D(L,L)=FactMlm0%D(L)*rs
       rs=rs*sq
    ENDDO
!
    DO L=0,Ell-1
       L1=L+1          
       LegendreP%D(L,L1)=CoTan*DBLE(2*L+1)*LegendreP%D(L,L)
    ENDDO
!
    DO L=2,Ell         
       Fact=CoTan*DBLE(2*L-1)
       DO M=0,L-2
          LegendreP%D(M,L)=LegendreP%D(M,L-1)*Fact-LegendreP%D(M,L-2)*FactMlm2%D(M,L)
       ENDDO
    ENDDO
!
    DO L=0,Ell
       o=RToTh%D(l)
       DO M=0,L
          LM = LTD(L)+M
          Cpq(LM)=o*LegendreP%D(M,L)*Cosine%D(M)
          Spq(LM)=o*LegendreP%D(M,L)*Sine%D(  M)
       ENDDO
    ENDDO
!
  END SUBROUTINE IrRegularF90
!---------------------------------------------------------------------------------------
! The Regular Multipoles
!---------------------------------------------------------------------------------------
  SUBROUTINE RegularF90(Ell,LenEll,x,y,z,Cpq,Spq)
    INTEGER                          :: Ell,LenEll
    INTEGER                          :: L,M,LM,L1,L2
    REAL(DOUBLE)                     :: x,y,z
    REAL(DOUBLE)                     :: x2,y2,R,OneOverR,Rxy,OneOvRxy,Rl,CoTan, &
                                        TwoC,Sq,rs,Fact,o
    REAL(DOUBLE),DIMENSION(0:LenEll) :: Cpq,Spq
!
    x2=x*x
    y2=y*y      
    R=SQRT(x2+y2+z*z)
!
    IF(R.EQ.Zero)THEN
       Cpq(0)=One
       Spq(0)=Zero
       DO LM=1,LenEll
          Cpq(LM)=Zero
          Spq(LM)=Zero
       ENDDO
       RETURN
    ENDIF
! 
    CoTan=z/R
    Cosine%D(0)=One
    Sine%D(  0)=Zero
!
    Rxy=SQRT(x2+y2)
    IF(Rxy.NE.Zero)THEN
       OneOvRxy=One/SQRT(x2+y2)       
       Cosine%D(1)=y*OneOvRxy
       Sine%D(  1)=x*OneOvRxy
    ELSE
       Cosine%D(1)=0.70710678118654752D0
       Sine%D(  1)=0.70710678118654752D0         
    ENDIF
!
    TwoC=Two*Cosine%D(1)
    DO L=2,Ell
       L1=L-1
       L2=L-2
       Cosine%D(L)=TwoC*Cosine%D(L1)-Cosine%D(L2)
       Sine%D(  L)=TwoC*Sine%D(  L1)-Sine%D(  L2)
    ENDDO
!
    Rl=One
    RS=One
    Sq=SQRT(One-CoTan*CoTan)
    DO L=0,Ell
       LegendreP%D(L,L)=FactOlm0%D(L)*RS
       RToTh%D(L)=Rl
       Rl=Rl*R
       RS=RS*Sq
    ENDDO
!
    DO L=0,Ell-1
       LegendreP%D(L,L+1)=CoTan*LegendreP%D(L,L)
    ENDDO
!
    DO L=2,Ell         
       Fact=CoTan*DBLE(2*L-1)
       DO M=0,L-2
          LegendreP%D(M,L)=(LegendreP%D(M,L-1)*Fact- LegendreP%D(M,L-2))*FactOlm2%D(M,L)
       ENDDO
    ENDDO
!
    DO L=0,Ell
       o=RToTh%D(L)
       DO M=0,L
          LM=LTD(L)+M
          Cpq(LM)=o*LegendreP%D(M,L)*Cosine%D(M)
          Spq(LM)=o*LegendreP%D(M,L)*Sine%D(  M)
       ENDDO
    ENDDO
!
  END SUBROUTINE RegularF90
!-----------------------------------------------------------------
! Translate the IrRegular Momments M_lm[Q-P] and Contract with Cq
!-----------------------------------------------------------------
  SUBROUTINE CTraXF90(LP,LQ,LPQ,LenP,LenQ,LenPQ,PQx,PQy,PQz,Cp,Sp,Cpq,Spq,Cq,Sq,TOF)
    INTEGER                          :: LP,LQ,LPQ,LenP,LenQ,LenPQ
    REAL(DOUBLE)                     :: PQx,PQy,PQz
    REAL(DOUBLE),DIMENSION(0:LenP)   :: Cp,Sp
    REAL(DOUBLE),DIMENSION(0:LenQ)   :: Cq,Sq
    REAL(DOUBLE),DIMENSION(0:LenPQ)  :: Cpq,Spq
!
    INTEGER                          :: L,LL,M,K,KK,N,LDX,LLKK,LKDX,KDX,LLKM
    REAL(DOUBLE)                     :: Tmp,TmpC,TmpS,Sgn,SgnL,SgnLM
    LOGICAL                          :: TOF
!
    IF(TOF) THEN
       CALL IrRegularF90(LPQ,LenPQ,PQx,PQy,PQz,Cpq,Spq)
    ENDIF
!
!     ZERO
!
      IF(LP == 0)THEN
         Tmp=Zero
         DO L=0,LQ
            LL=LSP(L)
            Cp(0)=Cp(0)+Cpq(LL)*Cq(ll)+Spq(LL)*Sq(LL)            
            DO M=1,L
              LDX=LL+M
              Tmp=Tmp+Cpq(LDX)*Cq(LDX)+Spq(LDX)*Sq(LDX)
           ENDDO
        ENDDO
        Cp(0)=Cp(0)+Two*Tmp
        RETURN
      ENDIF
!
!     ONE
!
      Sgn=One
      DO L=0,LP
         LL=LSP(L)
         TmpC=Zero
         DO K=0,LQ
            LL  =LSP(K)
            LLKK=LSP(L+K)
            TmpC=TmpC+Cpq(LLKK)*Cq(KK)
         ENDDO
         Cp(LL)=Cp(LL)+Sgn*TmpC
         Sgn=-Sgn
      ENDDO
!
!     TWO
!
      Sgn=-Two
      DO L=1,LP
         LL=LSP(L)
         DO K=0,LQ
            KK=LSP(K)
            LLKK=LSP(L+K)
            DO M=1,L                   
               LDX=LL+M
               LLKM=LLKK+M
               Cp(LDX)=Cp(LDX)+Sgn*Cpq(LLKM)*Cq(KK)
               Sp(LDX)=Sp(LDX)+Sgn*Spq(LLKM)*Cq(KK)
            ENDDO
            
         ENDDO
         Sgn=-Sgn
      ENDDO 
!
!     THREE
!
      Sgn=Two
      DO L=0,LP
         LL=LSP(L)
         TmpC=Zero
         DO K=0,LQ
            KK=LSP(K)
            LLKK=LSP(L+K)
            DO N=1,K
               KDX=KK+N
               LKDX=LLKK+N
               TmpC=TmpC+Cpq(LKDX)*Cq(KDX)+Spq(LKDX)*Sq(KDX)
            ENDDO
         ENDDO
         Cp(LL)=Cp(LL)+Sgn*TmpC
         Sgn=-Sgn
      ENDDO 
!
!     FOUR
!
      Sgn=-Two
      DO L=1,LP
         LL=LSP(L)
         DO K=1,LQ
            KK=LSP(K)
            LLKK=LSP(L+K)
            DO M=1,L                   
               LDX=LL+M
               LLKM=LLKK+M
               TmpC =Zero
               TmpS =Zero
               DO N=1,K
                  KDX=KK+N
                  LKDX=LLKM+N
                  TmpC=TmpC+Spq(LKDX)*Sq(KDX)+Cpq(LKDX)*Cq(KDX)
                  TmpS=TmpS-Cpq(LKDX)*Sq(KDX)+Spq(LKDX)*Cq(KDX)
               ENDDO
               Cp(LDX)=Cp(LDX)+Sgn*TmpC
               Sp(LDX)=Sp(LDX)+Sgn*TmpS
            ENDDO
         ENDDO
         Sgn=-Sgn
      ENDDO 
!
!     FIVE
!
      SgnL=One
      DO L=1,LP
         LL=LSP(L)
         DO K=1,LQ
            KK=LSP(K)
            LLKK=LSP(L+K)
            SgnLM=SgnL
            DO M=1,L      
               LDX=LL+M
               LLKM=LLKK+M
               Sgn=Two*SgnL   
               DO N=1,MIN(M,K)
                  KDX=KK+N
                  LKDX=LLKM-N
                  Cp(LDX)=Cp(LDX)+Sgn*Cpq(LKDX)*Cq(KDX)
                  Sp(LDX)=Sp(LDX)+Sgn*Cpq(LKDX)*Sq(KDX)
                  Sgn=-Sgn
               ENDDO
               LLKM=LLKK-M
               TmpC=Zero
               TmpS=Zero
               DO N=M+1,K
                  KDX=KK+N
                  LKDX=LLKM+N
                  TmpC=TmpC+Cpq(LKDX)*Cq(KDX)
                  TmpS=TmpS+Cpq(LKDX)*Sq(KDX)
               ENDDO
               Cp(LDX)=Cp(LDX)+SgnLM*Two*TmpC
               Sp(LDX)=Sp(LDX)+SgnLM*Two*TmpS
               SgnLM=-SgnLM
            ENDDO
         ENDDO
         SgnL=-SgnL
      ENDDO 
!!$C
!!$C     SIX
!!$C
!!$      SgnL=1.0D0
!!$      DO L=1,LP
!!$         LL=LSP(L)
!!$         DO K=1,LQ
!!$            KK=LSP(K)
!!$            LLKK=LSP(L+K)
!!$            SgnLM=SgnL
!!$            DO m=1,L      
!!$               Ldx=LL+m
!!$               LLKm=LLKK+m
!!$               Sgn=-Two*SgnL               
!!$               DO n=1,MIN(m-1,K)
!!$                  Kdx=KK+n
!!$                  LKdx=LLKm-n            
!!$                  Cp(Ldx)=Cp(Ldx)+Sgn*Spq(LKdx)*Sq(Kdx)
!!$                  Sp(Ldx)=Sp(Ldx)-Sgn*Spq(LKdx)*Cq(Kdx)
!!$                  Sgn=-Sgn
!!$               ENDDO
!!$               LLKm=LLKK-m
!!$               TmpC =Zero
!!$               TmpS =Zero
!!$               DO n=m+1,K
!!$                  Kdx=KK+n
!!$                  LKdx=LLKm+n
!!$                  TmpC=TmpC+Spq(LKdx)*Sq(Kdx)
!!$                  TmpS=TmpS+Spq(LKdx)*Cq(Kdx)
!!$               ENDDO
!!$               Cp(Ldx)=Cp(Ldx)+SgnLM*Two*TmpC
!!$               Sp(Ldx)=Sp(Ldx)-SgnLM*Two*TmpS          
!!$               SgnLM=-SgnLM
!!$            ENDDO
!!$         ENDDO
!!$         SgnL=-SgnL
!!$      ENDDO 
!
  END SUBROUTINE CTraXF90
!---------------------------------------------------------------
! Translate the Regular Momments O_lm[r-P]
!---------------------------------------------------------------
  SUBROUTINE XLateF90(LP,LQ,LPQ,LenP,LenQ,LenPQ,PQx,PQy,PQz,Cp,Sp,Cpq,Spq,Cq,Sq)
    INTEGER                          :: LP,LQ,LPQ,LenP,LenQ,LenPQ
    REAL(DOUBLE)                     :: PQx,PQy,PQz
    REAL(DOUBLE),DIMENSION(0:LenP)   :: Cp,Sp
    REAL(DOUBLE),DIMENSION(0:LenQ)   :: Cq,Sq
    REAL(DOUBLE),DIMENSION(0:LenPQ)  :: Cpq,Spq
!!$C
!!$      ID(L)=L*(L+1)/2
!!$C
!!$      CALL Regular(MaxL,LPQ,LenPQ,PQx,PQy,PQz,Cpq,Spq,
!!$     $          LegendreP,Cosine,Sine,RToTh,
!!$     $          FactOlm0,FactOlm1,FactOlm2) 
!!$C
!!$      DO l=0,LP
!!$         ll=LSP(L)
!!$         DO m=0,l
!!$            ldx=ll+m
!!$            DO k=0,MIN(l,LQ)
!!$               lk=l-k
!!$               kk=ID(k)
!!$               lklk=ID(lk)
!!$               nStart=MAX(-k, k-l+m)
!!$               nStop =MIN( k,-k+l+m)
!!$               DO n=nStart,nStop
!!$                  nabs=ABS(n)
!!$                  mnabs=ABS(m-n)
!!$                  kdx=kk+nabs                          
!!$                  lkdx=lklk+mnabs
!!$                  IF(m-n.LT.0)THEN
!!$                     cmn=(-1.0D0)**mnabs
!!$                     smn=(-1.0D0)**(mnabs+1)
!!$                  ELSE
!!$                     cmn=1.0D0
!!$                     smn=1.0D0
!!$                  ENDIF
!!$                  IF(n.LT.0)THEN
!!$                     cn=(-1.0D0)**nabs
!!$                     sn=(-1.0D0)**(nabs+1)
!!$                  ELSE
!!$                     cn=1.0D0
!!$                     sn=1.0D0
!!$                  ENDIF
!!$                  Cp(ldx)=Cp(ldx)+cn*Cq(kdx)*cmn*Cpq(lkdx)
!!$     $                           -sn*Sq(kdx)*smn*Spq(lkdx)        
!!$                  Sp(ldx)=Sp(ldx)+cn*Cq(kdx)*smn*Spq(lkdx)
!!$     $                           +sn*Sq(kdx)*cmn*Cpq(lkdx)
!!$               ENDDO
!!$            ENDDO
!!$         ENDDO
!!$      ENDDO
!
  END SUBROUTINE XLateF90
!---------------------------------------------------------------------------------------
! Given a Gaussian basis Function , calculate the Multipole 
! Momments in Spherical Basis
!---------------------------------------------------------------------------------------
  SUBROUTINE HGTFToSP(LQ,Zeta,RC,RS,Coefs)
    INTEGER                          :: LQ,LenQ
    REAL(DOUBLE)                     :: Zeta,PiZeta
    REAL(DOUBLE),DIMENSION(:)        :: Coefs
    REAL(DOUBLE),DIMENSION(0:)       :: RC,RS
!
    INCLUDE 'HGTFToSP.inc' 
!
  END SUBROUTINE HGTFToSP
!---------------------------------------------------------------------------------------
! Given a Gaussian basis Function , calculate the Multipole 
! Moments in Cartesion Basis
!---------------------------------------------------------------------------------------
  SUBROUTINE HGTFToCarts(LQ,LP,Zeta,RX,RY,RZ,QTemp,Coefs) 
    INTEGER                           :: LQ,LP,LCode
    REAL(DOUBLE)                      :: Zeta,PiZeta,RX,RY,RZ
    REAL(DOUBLE),DIMENSION(:)         :: Coefs
    REAL(DOUBLE),DIMENSION(0:)        :: QTemp
!
    INCLUDE 'HGTFToCarts.inc' 
!
  END SUBROUTINE HGTFToCarts
!
END MODULE MMML


