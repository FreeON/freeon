MODULE ECPBlock
  USE Parse
  USE Indexing
  USE AtomPairs
  USE DerivedTypes
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  IMPLICIT NONE
  ! Global temporary scatch array
  REAL(DOUBLE),DIMENSION(150) :: W
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
         PFC,ECPGL,ProjL,EllA,EllB,Lambda
    REAL(DOUBLE) :: Ax,Ay,Az,Bx,By,Bz,AB2,Cx,Cy,Cz,   &
         ACx,ACy,ACz,BCx,BCy,BCz,AC2,BC2,     &
         ZetaA,ZetaB,ZetaAB,ZetaC,Alpha,      &
         ZzAC2,ZzBC2,Px,Py,Pz,Gauss,CA,CB,CC, & 
         Kx,Ky,Kz,KAx,KAy,KAz,KBx,KBy,KBz,    &
         Kappa,KappaA,KappaB
    ! Static work arrays; dimensions defined in MondoMods/GlobalScalars
    REAL(DOUBLE),DIMENSION(0:HGEll) :: Omega1
    REAL(DOUBLE),DIMENSION(0:HGEll,0:HGEll+ECPEll) :: Keew1
    REAL(DOUBLE),DIMENSION(-PrjEll:PrjEll,0:PrjEll+BFEll,1:BFLen) :: OmegaA,OmegaB
    REAL(DOUBLE),DIMENSION(0:PrjEll+BFEll,0:PrjEll+BFEll,0:PrjEll+ECPEll) :: Keew2
    !--------------------------------------------------------------------------------------------------
    KA   = Pair%KA
    KB   = Pair%KB
    NBFA = Pair%NA
    NBFB = Pair%NB
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
                         ! Type one angular integrals, including 4 Pi and with mu terms
                         ! summed, removing the dependency on Cartesian angular indecies
                         CALL AngularOne(MaxLA+MaxLB,Kx,Ky,Kz,Omega1) 
                         ! Go over primitives in the unprojected, type one integrals
                         DO PFC=1,BS%NTyp1PF%I(KC)
                            ! Here is the type 1 primitive contraction coefficient,
                            CC=BS%Typ1CCo%D(PFC,KC)
                            ! the exponent, 
                            ZetaC=BS%Typ1Exp%D(PFC,KC)
                            ! and the radial "ell" value of the primitive Gaussian ECP
                            ECPGL=BS%Typ1Ell%I(PFC,KC)
                            ! Compute type one radial integrals
                            Alpha=ZetaA+ZetaB+ZetaC
                            CALL RadialOne(MaxLA+MaxLB+ECPGL,Kappa,Alpha,Keew1)
                            ! Go over basis function symmetries
                            DO LMNA=StartLA,StopLA
                               LA=BS%LxDex%I(LMNA)
                               MA=BS%LyDex%I(LMNA) 
                               NA=BS%LzDex%I(LMNA)
                               EllA=LA+MA+NA
                               CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                               DO LMNB=StartLB,StopLB
                                  LB=BS%LxDex%I(LMNB)
                                  MB=BS%LyDex%I(LMNB) 
                                  NB=BS%LzDex%I(LMNB)
                                  EllB=LB+MB+NB
                                  CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                                  ECP(LMNA,LMNB)=ECP(LMNA,LMNB)+CA*CB*CC*Gauss      &
                                       *ECPOneSum(LA,MA,NA,LB,MB,NB,ECPGL, &
                                       ACx,ACy,ACz,BCx,BCy,BCz,Omega1,Keew1) 
                               ENDDO
                            ENDDO
                         ENDDO
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
                               KappaA=Zero !! TEMP
                               KappaB=Zero !! TEMP
                               CALL RadialTwo(ProjL+MaxLA,ProjL+MaxLB,MaxLA+MaxLB+ECPGL, &
                                    Alpha,ZzAC2,KappaA,ZzBC2,KappaB,Keew2) 
                               ! Compute the type two angular integrals
                               CALL AngularTwo(ProjL,MaxLA,KAx,KAy,KAz,OmegaA)
                               CALL AngularTwo(ProjL,MaxLB,KBx,KBy,KBz,OmegaB)
                               ! Go over basis function symmetries
                               DO LMNA=StartLA,StopLA
                                  LA=BS%LxDex%I(LMNA)
                                  MA=BS%LyDex%I(LMNA) 
                                  NA=BS%LzDex%I(LMNA)
                                  EllA=LA+MA+NA
                                  CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                                  DO LMNB=StartLB,StopLB
                                     LB=BS%LxDex%I(LMNB)
                                     MB=BS%LyDex%I(LMNB) 
                                     NB=BS%LzDex%I(LMNB)
                                     EllB=LB+MB+NB
                                     CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                                     ECP(LMNA,LMNB)=ECP(LMNA,LMNB)+CA*CB*CC &
                                          *ECPTwoSum(LA,MA,NA,LB,MB,NB,ProjL,ECPGL, &
                                          ACx,ACy,ACz,BCx,BCy,BCz,OmegaA,OmegaB,Keew2)
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
    UVck = BlockToVect(Pair%NA,Pair%NB,ECP)
  END FUNCTION UBlock

  FUNCTION ECPOneSum(LA,MA,NA,LB,MB,NB,ECPGL,ACx,ACy,ACz,BCx,BCy,BCz,O,Q) RESULT(V1)
    INTEGER :: LA,MA,NA,LB,MB,NB,ECPGL,EllAB,  &
         ProjL,Emm,Lambda,LambdaA,LambdaB,AlphabetL, &
         a,b,c,d,e,f
    REAL(DOUBLE) :: ACx,ACy,ACz,BCx,BCy,BCz,V1, &
         BinA,BinB,BinC,BinD,BinE,BinF
    REAL(DOUBLE),DIMENSION(0:HGEll) :: O
    REAL(DOUBLE),DIMENSION(0:HGEll,0:HGEll+ECPEll) :: Q
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
                      AlphabetL=a+b+c+d+e+f+ECPGL
                      !->
                      DO Lambda=0,EllAB
                         V1=V1+BinF*Q(Lambda,AlphabetL)*O(Lambda)
                      ENDDO
                      !<-
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  END FUNCTION ECPOneSum

  FUNCTION ECPTwoSum(LA,MA,NA,LB,MB,NB,ProjL,ECPGL,ACx,ACy,ACz,BCx,BCy,BCz,OA,OB,Q) RESULT(V2)
    INTEGER :: LA,MA,NA,LB,MB,NB,ECPGL,EllAB,  &
         ProjL,Emm,LambdaA,LambdaB,AlphabetL, &
         a,b,c,d,e,f,LMNabc,LMNdef
    REAL(DOUBLE) :: ACx,ACy,ACz,BCx,BCy,BCz,V2, &
         BinA,BinB,BinC,BinD,BinE,BinF

    REAL(DOUBLE),DIMENSION(-PrjEll:PrjEll,0:PrjEll+BFEll,1:BFLen) :: OA,OB
    REAL(DOUBLE),DIMENSION(0:PrjEll+BFEll,0:PrjEll+BFEll,0:PrjEll+ECPEll) :: Q
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
                      DO LambdaA=0,ProjL+a+b+c
                         DO LambdaB=0,ProjL+d+e+f
                            DO Emm=-ProjL,ProjL
                               V2=V2+OA(Emm,LambdaA,LMNabc)*OB(Emm,LambdaB,LMNdef)*Q(LambdaA,LambdaB,AlphabetL) 
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
  END FUNCTION ECPTwoSum


  SUBROUTINE AngularOne(Ell,Kx,Ky,Kz,OmegaL) 
    INTEGER                     :: Ell
    REAL(DOUBLE)                :: Kx,Ky,Kz
    REAL(DOUBLE),DIMENSION(:)   :: OmegaL
    OmegaL=Zero
  END SUBROUTINE AngularOne

  SUBROUTINE RadialOne(Ell,K,Alpha,Q)
    INTEGER                     :: Ell
    REAL(DOUBLE)                :: K,Alpha
    REAL(DOUBLE),DIMENSION(:,:) :: Q
    Q=Zero
  END SUBROUTINE RadialOne

  SUBROUTINE AngularTwo(ProjL,Ell,Kx,Ky,Kz,Omega)
    INTEGER                                    :: ProjL,Ell
    REAL(DOUBLE)                               :: Kx,Ky,Kz
    REAL(DOUBLE), DIMENSION(:,:,:) :: Omega
    Omega=Zero
  END SUBROUTINE AngularTwo

  SUBROUTINE RadialTwo(EllAL,EllBL,EllABC,Alpha,ZAC2,KA,ZBC2,KB,Q)
    INTEGER                                          :: EllAL,EllBL,EllABC
    REAL(DOUBLE)                                     :: Alpha,ZAC2,KA,ZBC2,KB
    REAL(DOUBLE),DIMENSION(:,:,:) :: Q
    Q=Zero    
  END SUBROUTINE RadialTwo

  SUBROUTINE QTwo(EllA,EllB,N,Alpha,ZAC2,Ka,ZBC2,Kb,Q)
    INTEGER,PARAMETER                      :: NPts=100
    INTEGER                                :: EllA,EllB,N,LambdaA,LambdaB,I
    REAL(DOUBLE),DIMENSION(0:EllA,0:EllB)  :: Q
    REAL(DOUBLE),DIMENSION(1:NPts,0:BFEll) :: MA,MB
    REAL(DOUBLE),DIMENSION(1:NPts)         :: XA,XB,EX
    REAL(DOUBLE)                           :: ZAC2,ZBC2,R,Delta,Ka,Kb,Alpha,Left,Rhgt,Cntr,KaKb,DistToZero
    !
    DistToZero=-LOG(1D-8)
    KaKb=(Ka+Kb)
    Left=MAX(Zero,(KaKb-SQRT(KaKb**2+4D0*Alpha*DistToZero))/(Two*Alpha))
    Rhgt=(KaKb+SQRT(KaKb**2+4D0*Alpha*DistToZero))/(Two*Alpha)
    Delta=(Rhgt-Left)/(NPts)
    DO I=1,NPts
       R=Left+I*Delta
       XA(I)=KA*R
       XB(I)=KB*R
       EX(I)=EXP(-Alpha*R*R)*R**(Two+N)
    ENDDO
    DO LambdaA=0,EllA
       CALL BesselM(LambdaA,NPts,ZAC2,XA,MA(1,LambdaA))
    ENDDO
    DO LambdaB=0,EllB
       CALL BesselM(LambdaB,NPts,ZBC2,XB,MB(1,LambdaB))
    ENDDO
    DO LambdaA=0,EllA
       DO LambdaB=0,EllB
          Q(LambdaA,LambdaB)=Zero
          DO I=1,NPts
             Q(LambdaA,LambdaB)=Q(LambdaA,LambdaB)+MA(I,LambdaA)*MB(I,LambdaB)*EX(I)
          ENDDO
       ENDDO
    ENDDO
    Q=Delta*Q
  END SUBROUTINE QTwo

  SUBROUTINE BesselM(L,NPts,G2,X,M)
    INTEGER                        :: L,NPts,I
    REAL(DOUBLE),DIMENSION(1:NPts) :: M,X
    REAL(DOUBLE)                   :: G2,EPlus,EMnus,X2,ThreeX2,FifteenX3
    SELECT CASE(L)         
    CASE(0)
       DO I=1,NPts
          EPlus=Exp( X(I)-G2)
          EMnus=Exp(-X(I)-G2)
          M(I)=Half*(EPlus-EMnus)/X(I)
       ENDDO
    CASE(1)
       DO I=1,NPts

          M(I)=Half*(EPlus*(X(I)-One)+EMnus*(X(I)+One))/(X(I)*X(I))
       ENDDO
    CASE(2)
       DO I=1,NPts
          X2=X(I)*X(I)
          EPlus=Exp( X(I)-G2)
          EMnus=Exp(-X(I)-G2)
          ThreeX2=Three+X2
          M(I)=Half*(EPlus*(ThreeX2-Three*X(I))+EMnus*(ThreeX2+Three*X(I)))/(X(I)*X2)
       ENDDO
    CASE(3)
       DO I=1,NPts 
          X2=X(I)*X(I)
          EPlus=Exp( X(I)-G2)
          EMnus=Exp(-X(I)-G2)
          FifteenX3=(1.5D1+X2)*X(I)
          M(I)=Half*(EPlus*(FifteenX3-6.D0*X2-1.5D1)+EMnus*(FifteenX3+6.D0*X2+1.5D1))/(X2*X2)
       ENDDO
    CASE DEFAULT
       CALL Halt(' No code for case L = '//TRIM(IntToChar(L))//' in BesselM.')
    END SELECT
  END SUBROUTINE BesselM


#ifdef DLJFLSJFDLSJKF
  FUNCTION Omega1(Lambda,Kx,Ky,Kz)
    INTEGER      :: Lambda
    REAL(DOUBLE) :: Omega1,Kx,Ky,Kz
    SELECT CASE(Lambda)
       INCLUDE "Omega1.Inc"
    CASE DEFAULT
       CALL Halt(' No explicit code for case Lambda = '//TRIM(IntToChar(Lambda))//' in Omega1.')
    END SELECT
  END FUNCTION Omega1

  FUNCTION Omega2(l,a,b,c,KAx,KAy,KAz,e,d,f,KBx,KBy,KBz)
    INTEGER      :: l,a,b,c,e,d,f,labcedf
    REAL(DOUBLE) :: Omega2,KAx,KAy,KAz,KBx,KBy,KBz
    labcedf=10000000*l+1000000*a+10000*b+1000*c+100*d+10*e+f
    SELECT CASE(labcedf)
       INCLUDE "Omega2.Inc"
    CASE DEFAULT
       CALL Halt(' No explicit code for case labcdef = '//TRIM(IntToChar(labcedf))//' in Omega2.')
    END SELECT
  END FUNCTION Omega2
#endif

END MODULE ECPBlock
