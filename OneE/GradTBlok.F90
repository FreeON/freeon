
MODULE GradTBlock
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  USE BraKetBloks
  IMPLICIT NONE
  CONTAINS
  FUNCTION GradTBlok(BS,MD,AtA,AtB,DiffAtm,P,Pair) RESULT(TPVck)
!--------------------------------------------

!--------------------------------------------------------------------------------

  TYPE(BSET),                       INTENT(IN)    :: BS
  INTEGER                              :: NBFA,NBFB,KA,KB,      &
                                                     AtD,AtA,AtB,DiffAtm
  TYPE(DBL_RNK4),                   INTENT(INOUT) :: MD
  TYPE(AtomPair)                          :: Pair
  REAL(DOUBLE),DIMENSION(3)                       :: TPVck
  REAL(DOUBLE),DIMENSION(:)                       :: P
  REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)        :: GXP,GYP,GZP
  REAL(DOUBLE), EXTERNAL            :: BlkTrace_2
  REAL(DOUBLE)                         :: Ax,Ay,Az,Bx,By,Bz,AB2 
  REAL(DOUBLE)                                    :: Px,Py,Pz,PAx,PAy,PAz, &
                                                     PBx,PBy,PBz,ZetaA,    &
                                                     ZetaB,EtaAB,EtaIn,    &
                                                     XiAB,ExpAB,Ov,PiE32,  &
                                                     CA,CB,RLA,RLB,        &
                                                     RMA,RMB,RNA,RNB,      &
                                                     MDDerivAx,MDDerivBx,  &
                                                     MDDerivAy,MDDerivBy,  &
                                                     MDDerivAz,MDDerivBz,  &
                                                     OnA,OnB,              &
                                                     MDDerivDx,MDDerivDy,  &
                                                     MDDerivDz
  INTEGER                                         :: CFA,CFB,PFA,PFB,      &
                                                     IndexA,IndexB,        &
                                                     StartLA,StartLB,      &
                                                     StopLA,StopLB,        &
                                                     MaxLA,MaxLB,IA,IB,    &
                                                     LMNA,LMNB,LA,LB,MA,MB,&
                                                     NA,NB,DA
  CHARACTER                            :: OneE
!-------------------------------------------------------------------------------------- 
  GXP=Zero
  GYP=Zero
  GZP=Zero
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
!
  DO CFA=1,BS%NCFnc%I(KA)                     ! Loop over contracted function A 
     MaxLA=BS%ASymm%I(2,CFA,KA)
     DO CFB=1,BS%NCFnc%I(KB)                   ! Loop over contracted function B
        MaxLB=BS%ASymm%I(2,CFB,KB)
        DO PFA=1,BS%NPFnc%I(CFA,KA)            ! Loops over primitives in 
           DO PFB=1,BS%NPFnc%I(CFB,KB)         ! contracted functions A and B
              ZetaA=BS%Expnt%D(PFA,CFA,KA)
              ZetaB=BS%Expnt%D(PFB,CFB,KB)
              EtaAB=ZetaA+ZetaB 
              EtaIn=One/EtaAB
              XiAB =ZetaA*ZetaB*EtaIn
              ExpAB=EXP(-XiAB*AB2)
              PiE32=(Pi/EtaAB)**(1.5D0)
              Ov=ExpAB*PiE32
              Px=(ZetaA*Ax+ZetaB*Bx)*EtaIn
              Py=(ZetaA*Ay+ZetaB*By)*EtaIn
              Pz=(ZetaA*Az+ZetaB*Bz)*EtaIn
              PAx=Px-Ax
              PAy=Py-Ay
              PAz=Pz-Az
              PBx=Px-Bx
              PBy=Py-By
              PBz=Pz-Bz

              CALL MD2TRR(BS%NASym+4,-2,MaxLA+2,MaxLB+2,EtaAB,MD%D, &
                   PAx,PBx,PAy,PBy,PAz,PBz) 
!
!             Get the Derivative of S_AB
!
              OneE='T'
              CALL dOneE(DiffAtm,AtA,AtB,CFA,PFA,KA,CFB,PFB,KB,Ov,BS,MD,OneE, &
                   GXP,GYP,GZP,Pair)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  TPVck(1)=  BlkTrace_2(Pair%NB,Pair%NA,P,GXP(1,1))
  TPVck(2)=  BlkTrace_2(Pair%NB,Pair%NA,P,GYP(1,1))
  TPVck(3)=  BlkTrace_2(Pair%NB,Pair%NA,P,GZP(1,1))
END FUNCTION GradTBlok
END MODULE GradTBlock
