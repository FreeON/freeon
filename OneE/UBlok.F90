MODULE ECPBlock
  USE Parse
  USE Indexing
  USE AtomPairs
  USE DerivedTypes
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  IMPLICIT NONE
  ! Global temporary arrays for explicit fortran intermediates
  REAL(DOUBLE),DIMENSION(500) :: W,V
  ! Global solid angle.  Soliiidddd.
  REAL(DOUBLE),PARAMETER    :: FourPi=12.56637061435917295385D0
  ! Global array of binomial coefficients
  REAL(DOUBLE),DIMENSION(0:6,0:6) :: Binomial=RESHAPE((/     &
       1, 0, 0, 0, 0, 0, 0,    &
       1, 1, 0, 0, 0, 0, 0,    &
       1, 2, 1, 0, 0, 0, 0,    &
       1, 3, 3, 1, 0, 0, 0,    &
       1, 4, 6, 4, 1, 0, 0,    &
       1, 5,10,10, 5, 1, 0,    &
       1, 6,15,20,15, 6, 1 /),(/7,7/))
  !---------------------------------------------------------------------
CONTAINS  !
  FUNCTION UBlock(BS,GM,Pair) RESULT(UVck)
    TYPE(BSET)                              :: BS
    TYPE(CRDS)                              :: GM
    TYPE(AtomPair)                          :: Pair
    TYPE(PrimPair)                          :: Prim
    REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB) :: UVck
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB) :: ECP
    INTEGER :: NBFA,NBFB,KA,KB,KK,KC,CFA,CFB,PFA, &
         PFB,IndexA,IndexB,StartLA,StartLB, &
         StopLA,StopLB,MaxLA,MaxLB,IA,IB,   &
         LMNA,LMNB,LA,LB,MA,MB,NA,NB,AtC,   &
         PFC,ECPGL,ProjL,EllA,EllB,Lambda,LenAB
    REAL(DOUBLE) :: Ax,Ay,Az,Bx,By,Bz,AB2,Cx,Cy,Cz,   &
         ACx,ACy,ACz,BCx,BCy,BCz,AC2,BC2,AC,BC,&
         ZetaA,ZetaB,ZetaAB,ZetaC,Alpha,       &
         ZzAC2,ZzBC2,Px,Py,Pz,Gauss,CA,CB,CC,  & 
         Kx,Ky,Kz,KAx,KAy,KAz,KBx,KBy,KBz,     &
         Kappa,KappaA,KappaB, tmp
    ! Static work arrays; dimensions defined in MondoMods/GlobalScalars
    REAL(DOUBLE),DIMENSION(0:HGEll,1:HGLen) :: Omega1
    REAL(DOUBLE),DIMENSION(0:HGEll,0:HGEll+ECPEll) :: Keew1
    REAL(DOUBLE),DIMENSION(0:PrjEll+BFEll,-PrjEll:PrjEll,1:BFLen) :: OmegaA,OmegaB
    REAL(DOUBLE),DIMENSION(0:PrjEll+BFEll,0:PrjEll+BFEll,0:PrjEll+ECPEll) :: Keew2
    TYPE(DBL_RNK2) :: d
    !--------------------------------------------------------------------------------------------------
    KA   = Pair%KA
    KB   = Pair%KB
    NBFA = Pair%NA
    NBFB = Pair%NB
    
!    CALL New(d,(/NBFA,NBFB/))

    Ax = Pair%A(1)
    Ay = Pair%A(2)
    Az = Pair%A(3)
    Bx = Pair%B(1)
    By = Pair%B(2)
    Bz = Pair%B(3)
    AB2=Pair%AB2
    ECP=Zero
    ! Go over ECP centers
    DO AtC=1,NAtoms
       KC=GM%AtTyp%I(AtC)
       IF(BS%NTyp1PF%I(KC)>0)THEN
          Cx=GM%Carts%D(1,AtC) 
          Cy=GM%Carts%D(2,AtC) 
          Cz=GM%Carts%D(3,AtC) 
          ACx=Ax-Cx
          ACy=Ay-Cy
          ACz=Az-Cz
          BCx=Bx-Cx
          BCy=By-Cy
          BCz=Bz-Cz
          AC2=(ACx**2+ACy**2+ACz**2)
          BC2=(BCx**2+BCy**2+BCz**2)
          AC=SQRT(AC2)
          BC=SQRT(BC2)
          ! Some testing of overlaps
          IF(.TRUE.)THEN
             ! Go over basis functions for this atom-atom blok
             DO CFA=1,BS%NCFnc%I(KA)           
                IndexA=CFBlokDex(BS,CFA,KA)
                StartLA=BS%LStrt%I(CFA,KA)        
                StopLA =BS%LStop%I(CFA,KA)
                MaxLA=BS%ASymm%I(2,CFA,KA)
                DO CFB=1,BS%NCFnc%I(KB)        
                   IndexB  = CFBlokDex(BS,CFB,KB)
                   StartLB = BS%LStrt%I(CFB,KB)
                   StopLB  = BS%LStop%I(CFB,KB)
                   MaxLB   = BS%ASymm%I(2,CFB,KB)       
                   LenAB   = LHGTF(MaxLA+MaxLB)
                   DO PFA=1,BS%NPFnc%I(CFA,KA) 
                      DO PFB=1,BS%NPFnc%I(CFB,KB)        
                         ZetaA=BS%Expnt%D(PFA,CFA,KA)
                         ZetaB=BS%Expnt%D(PFB,CFB,KB)
                         ZetaAB=ZetaA+ZetaB
                         ZzAC2=ZetaA*AC2
                         ZzBC2=ZetaB*BC2
                         Px=(zetaA*Ax+zetaB*Bx)/ZetaAB
                         Py=(zetaA*Ay+zetaB*By)/ZetaAB
                         Pz=(zetaA*Az+zetaB*Bz)/ZetaAB
                         Kx=Px-Cx
                         Ky=Py-Cy
                         Kz=Pz-Cz
                         Kappa=Two*ZetaAB*SQRT(Kx*Kx+Ky*Ky+Kz*Kz)
                         Gauss=EXP(-ZzAC2-ZzBC2)
                         ! Type one angular integrals, including 4 Pi and with mu terms summed
                         CALL AngularOne(MaxLA+MaxLB,Kx,Ky,Kz,Omega1) 
!                         WRITE(*,*)' Omega = ',Omega1(0:MaxLA+MaxLB,1:LenAB)
                         ! Go over primitives in the unprojected, type one integrals
                         DO PFC=1,BS%NTyp1PF%I(KC)
                            ! Here is the type 1 primitive contraction coefficient,
                            CC=BS%Typ1CCo%D(PFC,KC)
                            ! the exponent, 
                            ZetaC=BS%Typ1Exp%D(PFC,KC)
                            ! and the radial "ell" value of the primitive Gaussian ECP
                            ECPGL=BS%Typ1Ell%I(PFC,KC)
!                            WRITE(*,33)'a',Ax,Ay,Az
!                            WRITE(*,33)'b',Bx,By,Bz
!                            WRITE(*,33)'c',Cx,Cy,Cz
!                            WRITE(*,44)'A',ZetaA
!                            WRITE(*,44)'B',ZetaB
!                            WRITE(*,44)'C',ZetaC
!                            WRITE(*,*)' EllC = ',ECPGL,';'
!33                          format(1x,A1,' = {',F20.10,', ',F20.10,', ',F20.10,'};')
!44                          format(1x,'Zeta',A1,' = ',F20.10,';')
                            ! Compute type one radial integrals
                            Alpha=ZetaA+ZetaB+ZetaC
                            CALL RadialOne(MaxLA+MaxLB,MaxLA+MaxLB+ECPGL,Alpha,Kappa,ZzAC2+ZzBC2,Keew1)
!                            WRITE(*,*)' Keew1 = ',Keew1(0:MaxLA+MaxLB,MaxLA+MaxLB+ECPGL)
                            ! Go over basis function symmetries
                            IA=IndexA
                            DO LMNA=StartLA,StopLA
                               LA=BS%LxDex%I(LMNA)
                               MA=BS%LyDex%I(LMNA) 
                               NA=BS%LzDex%I(LMNA)
                               EllA=LA+MA+NA
                               CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                               IA=IA+1
                               IB=IndexB
                               DO LMNB=StartLB,StopLB
                                  LB=BS%LxDex%I(LMNB)
                                  MB=BS%LyDex%I(LMNB) 
                                  NB=BS%LzDex%I(LMNB)
                                  EllB=LB+MB+NB
                                  CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                                  IB=IB+1

!                                  WRITE(*,*)'EllA = {',LA,', ',MA,', ',NA,'};'
!                                  WRITE(*,*)'EllB = {',LB,', ',MB,', ',NB,'};'

                                  ECP(IA,IB)=ECP(IA,IB) &
                                       +CA*CB*CC*ECPOneSum(LA,MA,NA,LB,MB,NB,ECPGL,         &
                                        ACx,ACy,ACz,BCx,BCy,BCz,Omega1,Keew1) 

!                                  TMP=ECPOneSum(LA,MA,NA,LB,MB,NB,ECPGL, &
!                                       ACx,ACy,ACz,BCx,BCy,BCz,Omega1,Keew1) 
!                                  WRITE(*,*)' Type1 = ',Tmp

                               ENDDO
                            ENDDO
                         ENDDO
!                         STOP
                         ! Loop on ECP projector symmetries 
                         DO ProjL=0,BS%ProjEll%I(KC)
                            ! Go over primitives in the projected, type two integrals
                            DO PFC=1,BS%NTyp2PF%I(ProjL,KC)                     
                               ! Here is the type 2 primitive contraction coefficient,
                               CC=BS%Typ2CCo%D(PFC,ProjL,KC)
                               ! the corresponding exponent, 
                               ZetaC=BS%Typ2Exp%D(PFC,ProjL,KC)     
                               ! and the radial "ell" value of the ECP primitive 
                               ECPGL=BS%Typ2Ell%I(PFC,ProjL,KC)
                               ! Compute the type two radial integrals                            
                               Alpha=ZetaA+ZetaB+ZetaC
                               KappaA=Two*ZetaA*AC
                               KappaB=Two*ZetaB*BC

!                            WRITE(*,33)'a',Ax,Ay,Az
!                            WRITE(*,33)'b',Bx,By,Bz
!                            WRITE(*,33)'c',Cx,Cy,Cz
!                            WRITE(*,44)'A',ZetaA
!                            WRITE(*,44)'B',ZetaB
!                            WRITE(*,44)'C',ZetaC
!                            WRITE(*,*)' EllC = ',ECPGL,';'
!33                          format(1x,A1,' = {',F20.10,', ',F20.10,', ',F20.10,'};')
!44                          format(1x,'Zeta',A1,' = ',F20.10,';')


                               CALL RadialTwo(ProjL+MaxLA,ProjL+MaxLB,MaxLA+MaxLB+ECPGL, &
                                    Alpha,ZzAC2,KappaA,ZzBC2,KappaB,Keew2) 
                               ! Compute the type two angular integrals
                               CALL AngularTwo(ProjL,MaxLA,KAx,KAy,KAz,OmegaA)
                               CALL AngularTwo(ProjL,MaxLB,KBx,KBy,KBz,OmegaB)
                               ! Go over basis function symmetries
                               IA=IndexA
                               DO LMNA=StartLA,StopLA
                                  LA=BS%LxDex%I(LMNA)
                                  MA=BS%LyDex%I(LMNA) 
                                  NA=BS%LzDex%I(LMNA)
                                  EllA=LA+MA+NA
                                  CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                                  IA=IA+1
                                  IB=IndexB
                                  DO LMNB=StartLB,StopLB
                                     LB=BS%LxDex%I(LMNB)
                                     MB=BS%LyDex%I(LMNB) 
                                     NB=BS%LzDex%I(LMNB)
                                     EllB=LB+MB+NB
                                     CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                                     IB=IB+1
                                     ECP(IA,IB)=ECP(IA,IB)+CA*CB*CC &
                                         *ECPTwoSum(LA,MA,NA,LB,MB,NB,ProjL,ECPGL,   &
                                          ACx,ACy,ACz,BCx,BCy,BCz,OmegaA,OmegaB,Keew2) 
!                                     WRITE(*,*)' CA = ',CA
!                                     WRITE(*,*)' CB = ',CB
!                                     WRITE(*,*)' CC = ',CC
!                                     WRITE(*,*)' ECPTwoSum = ',ECPTwoSum(LA,MA,NA,LB,MB,NB,ProjL,ECPGL, &
!                                           ACx,ACy,ACz,BCx,BCy,BCz,OmegaA,OmegaB,Keew2)
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDIF
    ENDDO
    
!    D%d=ECP
!    CALL PPrint(d,"D",Unit_O=6)
!    CALL Delete(D)
!   STOP
    UVck = BlockToVect(Pair%NA,Pair%NB,ECP)
  END FUNCTION UBlock

  FUNCTION ECPOneSum(LA,MA,NA,LB,MB,NB,ECPGL,ACx,ACy,ACz,BCx,BCy,BCz,O,Q) RESULT(V1)
    INTEGER :: LA,MA,NA,LB,MB,NB,ECPGL,EllAB,  &
         ProjL,Emm,Lambda,LambdaA,LambdaB,AlphabetL, &
         a,b,c,d,e,f,ABCDEF, lll,vvv
    REAL(DOUBLE) :: ACx,ACy,ACz,BCx,BCy,BCz,V1, &
         BinA,BinB,BinC,BinD,BinE,BinF
    REAL(DOUBLE),DIMENSION(0:,1:) :: O
    REAL(DOUBLE),DIMENSION(0:,0:) :: Q

!    WRITE(*,*)' ======================================================================'
!    DO LLL=0,LA+MA+NA+MB+NB
!       DO vvv=0,LA+MA+NA+MB+NB+ECPGL       
!          WRITE(*,*)LLL,VVV,' Q2 = ',Q(LLL,VVV)
!       ENDDO
!    ENDDO

    ACx=ACx+1D-100
    ACy=ACy+1D-100
    ACz=ACz+1D-100
    BCx=BCx+1D-100
    BCy=BCy+1D-100
    BCz=BCz+1D-100
    !
    V1=Zero
    EllAB=LA+MA+NA+LB+MB+NB
    DO a=0,LA
       BinA=ACx**(LA-a)*Binomial(LA,a)
       DO b=0,MA
          BinB=BinA*ACy**(MA-b)*Binomial(MA,b)
          DO c=0,NA
             BinC=BinB*ACz**(NA-c)*Binomial(NA,c)
             DO d=0,LB
                BinD=BinC*BCx**(LB-d)*Binomial(LB,d)
                DO e=0,MB
                   BinE=BinD*BCy**(MB-e)*Binomial(MB,e)
                   DO f=0,NB
                      BinF=BinE*BCz**(NB-f)*Binomial(NB,f)
                      ! All the funny indexing ...
                      AlphabetL=a+b+c+d+e+f+ECPGL
                      ABCDEF=LMNDex(a+d,b+e,c+f)
                      !->
                      DO Lambda=0,EllAB
                         V1=V1+BinF*Q(Lambda,AlphabetL)*O(Lambda,ABCDEF)
!IF(LA==LB)THEN
!                         WRITE(*,*)Lambda,AlphabetL,ABCDEF
!                         WRITE(*,*)' BinF = ',BinF 
!                         WRITE(*,*)' O = ',O(Lambda,ABCDEF)
!                         WRITE(*,*)' Q = ',Q(Lambda,AlphabetL)
!ENDIF
                      ENDDO
                      !<-
                   ENDDO
                ENDDO
             ENDDO 
          ENDDO
       ENDDO
    ENDDO
    V1=V1*FourPi
  END FUNCTION ECPOneSum

  FUNCTION ECPTwoSum(LA,MA,NA,LB,MB,NB,ProjL,ECPGL,ACx,ACy,ACz,BCx,BCy,BCz,OA,OB,Q) RESULT(V2)
    INTEGER :: LA,MA,NA,LB,MB,NB,ECPGL,EllAB,  &
         ProjL,Emm,LambdaA,LambdaB,AlphabetL, &
         a,b,c,d,e,f,LMNabc,LMNdef
    REAL(DOUBLE) :: ACx,ACy,ACz,BCx,BCy,BCz,V2, &
         BinA,BinB,BinC,BinD,BinE,BinF

    REAL(DOUBLE),DIMENSION(0:,-ProjL:,1:) :: OA,OB
    REAL(DOUBLE),DIMENSION(0:,0:,0:) :: Q
!
    ACx=ACx+1D-100
    ACy=ACy+1D-100
    ACz=ACz+1D-100
    BCx=BCx+1D-100
    BCy=BCy+1D-100
    BCz=BCz+1D-100
    !
    V2=Zero
    EllAB=LA+MA+NA+LB+MB+NB
    DO a=0,LA
       BinA=ACx**(LA-a)*Binomial(LA,a)
       DO b=0,MA
          BinB=BinA*ACy**(MA-b)*Binomial(MA,b)
          DO c=0,NA
             BinC=BinB*ACz**(NA-c)*Binomial(NA,c)
             LMNabc=LMNDex(a,b,c)
             DO d=0,LB
                BinD=BinC*BCx**(LB-d)*Binomial(LB,d)
                DO e=0,MB
                   BinE=BinD*BCy**(MB-e)*Binomial(MB,e)
                   DO f=0,NB
                      LMNdef=LMNDex(d,e,f)
                      AlphabetL=a+b+c+d+e+f+ECPGL
                      BinF=BinE*BCz**(NB-f)*Binomial(NB,f)
                      ! Contract the type two (projected) contribution
!WRITE(*,*)'============================================================================='
!WRITE(*,*)' ProjL = ',ProjL
                      DO LambdaA=0,ProjL+a+b+c
                         DO LambdaB=0,ProjL+d+e+f
                            DO Emm=-ProjL,ProjL
                               V2=V2+BinF*OA(Emm,LambdaA,LMNabc)*OB(Emm,LambdaB,LMNdef)*Q(LambdaA,LambdaB,AlphabetL) 

!WRITE(*,*)'A ',Emm,LambdaA,LMNabc,OA(Emm,LambdaA,LMNabc)
!WRITE(*,*)'B ',Emm,LambdaB,LMNdef,OB(Emm,LambdaB,LMNdef)
!WRITE(*,*)'R ',LambdaA,LambdaB,AlphabetL,Q(LambdaA,LambdaB,AlphabetL) 

                            ENDDO
                         ENDDO
                      ENDDO
                      !<-
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    V2=V2*FourPi
  END FUNCTION ECPTwoSum

  SUBROUTINE AngularOne(Ell,X,Y,Z,O) 
    INTEGER                       :: Ell,  l,lmn
    REAL(DOUBLE)                  :: K,X,Y,Z,Kx,Ky,Kz
    REAL(DOUBLE),DIMENSION(0:,1:) :: O
    K=SQRT(X**2+Y**2+Z**2)
    IF(K/=Zero)THEN
       Kx=X/K
       Ky=Y/K
       Kz=Z/K
    ELSE
       Kx=Zero
       Ky=Zero
       Kz=One
    ENDIF
    SELECT CASE(Ell)
    INCLUDE "Omega1.Inc"
    CASE DEFAULT
       CALL Halt(' No explicit code for case Ell = '//TRIM(IntToChar(Ell))//' in AngularOne.')
    END SELECT

!    DO L=0,Ell
!       DO LMN=1,LHGTF(Ell)
!          WRITE(*,*)L,LMN,' Omega = ',O(L,LMN)
!       ENDDO
!    ENDDO
  END SUBROUTINE AngularOne

  SUBROUTINE AngularTwo(ProjL,Ell,X,Y,Z,O)
    INTEGER                                    :: ProjL,Ell,Index
    REAL(DOUBLE)                               :: X,Y,Z,Kx,Ky,Kz,K
    REAL(DOUBLE), DIMENSION(0:,-ProjL:,1:)     :: O

!    integer :: L,P,LMN

    K=SQRT(X**2+Y**2+Z**2)
    IF(K/=Zero)THEN
       Kx=X/K
       Ky=Y/K
       Kz=Z/K
    ELSE
       Kx=Zero
       Ky=Zero
       Kz=Zero
    ENDIF
    Index=100*Ell+ProjL
    SELECT CASE(Index)
    INCLUDE "Omega2.Inc"
    CASE DEFAULT
       CALL Halt(' No explicit code for case Ell = '//TRIM(IntToChar(Ell))//' in AngularTwo.')
    END SELECT
!    DO L=0,Ell
!       DO P=-ProjL,ProjL
!          DO LMN=1,LHGTF(Ell)
!              WRITE(*,*)L,P,LMN,' Omega = ',O(L,P,LMN)
!           ENDDO
!        ENDDO
!     ENDDO
  END SUBROUTINE AngularTwo

  SUBROUTINE RadialOne(Ell,N,Alpha,K,ZAC2pZBC2,Q)
    INTEGER,PARAMETER                      :: NPts=10000
    INTEGER                                :: Ell,Ehn,N,I,Lambda
    REAL(DOUBLE)                           :: K,Alpha,ZAC2pZBC2,R,Delta,Left,Rhgt,Cntr,DistToZero,AlphaInv
    REAL(DOUBLE),DIMENSION(1:NPts,0:BFEll) :: M
    REAL(DOUBLE),DIMENSION(1:NPts)         :: KV,EX,RV
    REAL(DOUBLE),DIMENSION(0:,0:)          :: Q
    !
    IF(K==Zero)THEN
       Q=Zero
       AlphaInv=One/Alpha
       ! Only values with Lambda=0 survive
       Q(0,0)=Half*SQRT(Pi*AlphaInv)
       Q(0,1)=Half*AlphaInv
       DO Ehn=2,N-1
          Q(0,Ehn+1)=Half*AlphaInv*DBLE(Ehn)*Q(0,Ehn-1)
       ENDDO
    ELSE
       DistToZero=-LOG(1D-10)
       Left=MAX(Zero,(K-SQRT(K**2+4D0*Alpha*DistToZero))/(Two*Alpha))
       Rhgt=(K+SQRT(K**2+4D0*Alpha*DistToZero))/(Two*Alpha)
       Delta=(Rhgt-Left)/(NPts)
       DO I=1,NPts
          R=Left+I*Delta
          RV(I)=R
          KV(I)=K*R
          EX(I)=EXP(-Alpha*R*R)
       ENDDO
       DO Lambda=0,Ell
          CALL BesselM(Lambda,NPts,ZAC2pZBC2,KV,M(1,Lambda))
       ENDDO
       DO Ehn=0,N
          DO Lambda=0,Ell
             Q(Lambda,Ehn)=Zero
             DO I=1,NPts
                Q(Lambda,Ehn)=Q(Lambda,Ehn)+M(I,Lambda)*EX(I)*RV(I)**(Ehn)
                !             WRITE(*,*)I,RV(I),M(I,Lambda)*EX(I)*RV(I)**(Two+Ehn)
             ENDDO
          ENDDO
       ENDDO
       Q=Delta*Q
    ENDIF
!    write(*,*)'=============================================================='
!    write(*,*)' Ell = ',Ell,' N = ',N
!    DO Ehn=0,N
!       DO Lambda=0,Ell
!          WRITE(*,*)Lambda,Ehn,Q(Lambda,Ehn)
!       ENDDO
!    ENDDO
  END SUBROUTINE RadialOne

  SUBROUTINE RadialTwo(EllA,EllB,N,Alpha,ZAC2,Ka,ZBC2,Kb,Q)
    INTEGER,PARAMETER                         :: NPts=10000
    INTEGER                                   :: EllA,EllB,N,LambdaA,LambdaB,I,Ehn
    REAL(DOUBLE),DIMENSION(0:,0:,0:)          :: Q
    REAL(DOUBLE),DIMENSION(1:NPts,0:BFEll)    :: MA,MB
    REAL(DOUBLE),DIMENSION(1:NPts)            :: XA,XB,EX,RV
    REAL(DOUBLE)                              :: ZAC2,ZBC2,R,Delta,Ka,Kb,Alpha,Left,Rhgt,Cntr,KaKb,DistToZero
    !
    DistToZero=-LOG(1D-10)
    KaKb=(Ka+Kb)
    Left=MAX(Zero,(KaKb-SQRT(KaKb**2+4D0*Alpha*DistToZero))/(Two*Alpha))
    Rhgt=(KaKb+SQRT(KaKb**2+4D0*Alpha*DistToZero))/(Two*Alpha)
    Delta=(Rhgt-Left)/(NPts)
    DO I=1,NPts
       R=Left+I*Delta
       RV(I)=R
       XA(I)=KA*R
       XB(I)=KB*R
       EX(I)=EXP(-Alpha*R*R)
    ENDDO
    DO LambdaA=0,EllA
       CALL BesselM(LambdaA,NPts,ZAC2,XA,MA(1,LambdaA))
    ENDDO
    DO LambdaB=0,EllB
       CALL BesselM(LambdaB,NPts,ZBC2,XB,MB(1,LambdaB))
    ENDDO
    DO Ehn=0,N
       DO LambdaA=0,EllA
          DO LambdaB=0,EllB
             Q(LambdaA,LambdaB,Ehn)=Zero
             DO I=1,NPts
                Q(LambdaA,LambdaB,Ehn)=Q(LambdaA,LambdaB,Ehn)+MA(I,LambdaA)*MB(I,LambdaB)*EX(I)*RV(I)**(Ehn)
!                WRITE(*,*)RV(I),MA(I,LambdaA)*MB(I,LambdaB)*EX(I)*RV(I)**(Ehn)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    Q=Delta*Q
!    write(*,*)'=============================================================='
!    write(*,*)' Ell = ',EllA,' EllB = ',EllB,' N = ',N
!    DO Ehn=0,N
!       DO LambdaA=0,EllA
!          DO LambdaB=0,EllB
!             WRITE(*,*)Ehn,LambdaA,LambdaB,' Q = ',Q(LambdaA,LambdaB,Ehn)
!          ENDDO
!       ENDDO
!    ENDDO
  END SUBROUTINE RadialTwo

  SUBROUTINE BesselM(L,NPts,G2,X,M)
    INTEGER                        :: L,NPts,I
    REAL(DOUBLE),DIMENSION(1:NPts) :: M,X
    REAL(DOUBLE)                   :: G2,EPlus,EMnus,X2,ThreeX2,FifteenX3
    SELECT CASE(L)         
    CASE(0)       
       DO I=1,NPts
          IF(X(I)==Zero)THEN
             M(I)=One
          ELSE
             EPlus=Exp( X(I)-G2)
             EMnus=Exp(-X(I)-G2)
             M(I)=Half*(EPlus-EMnus)/X(I)
          ENDIF
       ENDDO
    CASE(1)
       DO I=1,NPts
          IF(X(I)==Zero)THEN
             M(I)=Zero
          ELSE
             EPlus=Exp( X(I)-G2)
             EMnus=Exp(-X(I)-G2)
             M(I)=Half*(EPlus*(X(I)-One)+EMnus*(X(I)+One))/(X(I)*X(I))
          ENDIF
       ENDDO
    CASE(2)
       DO I=1,NPts
          IF(X(I)==Zero)THEN
             M(I)=Zero
          ELSE
             X2=X(I)*X(I)
             EPlus=Exp( X(I)-G2)
             EMnus=Exp(-X(I)-G2)
             ThreeX2=Three+X2
             M(I)=Half*(EPlus*(ThreeX2-Three*X(I))+EMnus*(ThreeX2+Three*X(I)))/(X(I)*X2)
          ENDIF
       ENDDO
    CASE(3)
       DO I=1,NPts 
          IF(X(I)==Zero)THEN
             M(I)=Zero
          ELSE
             X2=X(I)*X(I)
             EPlus=Exp( X(I)-G2)
             EMnus=Exp(-X(I)-G2)
             FifteenX3=(1.5D1+X2)*X(I)
             M(I)=Half*(EPlus*(FifteenX3-6.D0*X2-1.5D1)+EMnus*(FifteenX3+6.D0*X2+1.5D1))/(X2*X2)
          ENDIF
       ENDDO
    CASE DEFAULT
       CALL Halt(' No code for case L = '//TRIM(IntToChar(L))//' in BesselM.')
    END SELECT
!    WRITE(*,*)' ======================================================'
!    WRITE(*,*)' Ell = ',l,' Ell = ',l,' Ell = ',l,' Ell = ',l
!    WRITE(*,*)' X = ',X
!    WRITE(*,*)' M = ',M
  END SUBROUTINE BesselM

END MODULE ECPBlock
