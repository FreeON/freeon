
MODULE GradSBlock
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  IMPLICIT NONE
  CONTAINS
  FUNCTION GradSBlok(BS,MD,AtA,AtB,DiffAtm,W,Pair) RESULT(SPVck)
!--------------------------------------------

!--------------------------------------------------------------------------------

  TYPE(BSET),                       INTENT(IN)    :: BS
  INTEGER                              :: NBFA,NBFB,KA,KB,      &
                                                     AtD,AtA,AtB,DiffAtm
  TYPE(DBL_RNK4),                   INTENT(INOUT) :: MD
  TYPE(AtomPair)                          :: Pair
  REAL(DOUBLE),DIMENSION(3)                       :: SPVck
  REAL(DOUBLE),DIMENSION(:)                       :: W
  REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)        :: SXP,SYP,SZP
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
!-------------------------------------------------------------------------------------- 
  SXP=Zero
  SYP=Zero
  SZP=Zero
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
  IF (DiffAtm==AtA.AND.AtA==AtB) THEN
     OnA=One
     OnB=One
  ELSEIF (DiffAtm==AtA) THEN
     OnA=One
     OnB=Zero
  ELSEIF (DiffAtm==AtB) THEN
     OnB=One
     OnA=Zero
  ENDIF
!
  IndexA=0
  DO CFA=1,BS%NCFnc%I(KA)                     ! Loop over contracted function A 
     IndexB=0
     StartLA=BS%LStrt%I(CFA,KA)        
     StopLA=BS%LStop%I(CFA,KA)
     MaxLA=BS%ASymm%I(2,CFA,KA)
     DO CFB=1,BS%NCFnc%I(KB)                   ! Loop over contracted function B
        StartLB=BS%LStrt%I(CFB,KB)
        StopLB =BS%LStop%I(CFB,KB)
        MaxLB=BS%ASymm%I(2,CFB,KB)

        DO PFA=1,BS%NPFnc%I(CFA,KA)            ! Loops over primitives in 
           DO PFB=1,BS%NPFnc%I(CFB,KB)            ! contracted functions A and B
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

              CALL MD2TRR(BS%NASym+1,-1,MaxLA+1,MaxLB+1,EtaAB,MD%D, &
                   PAx,PBx,PAy,PBy,PAz,PBz) 
              
              IB=IndexB
              DO LMNB=StartLB,StopLB
                 LB=BS%LxDex%I(LMNB)
                 MB=BS%LyDex%I(LMNB)
                 NB=BS%LzDex%I(LMNB)
                 CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                 IB=IB+1
                 IA=IndexA
                 DO LMNA=StartLA,StopLA
                    IA=IA+1
                    LA=BS%LxDex%I(LMNA)
                    MA=BS%LyDex%I(LMNA)
                    NA=BS%LzDex%I(LMNA)
                    CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                    
                    RLA=DBLE(LA)
                    RMA=DBLE(MA)
                    RNA=DBLE(NA)
                    RLB=DBLE(LB)
                    RMB=DBLE(MB)
                    RNB=DBLE(NB)

                    MDDerivAx=Two*ZetaA*MD%D(1,LA+1,LB,0) &
                         -RLA*MD%D(1,LA-1,LB,0)
                    MDDerivBx=Two*ZetaB*MD%D(1,LA,LB+1,0) &
                         -RLB*MD%D(1,LA,LB-1,0)
                    MDDerivAy=Two*ZetaA*MD%D(2,MA+1,MB,0) &
                         -RMA*MD%D(2,MA-1,MB,0)
                    MDDerivBy=Two*ZetaB*MD%D(2,MA,MB+1,0) &
                         -RMB*MD%D(2,MA,MB-1,0)
                    MDDerivAz=Two*ZetaA*MD%D(3,NA+1,NB,0) &
                         -RNA*MD%D(3,NA-1,NB,0)
                    MDDerivBz=Two*ZetaB*MD%D(3,NA,NB+1,0) &
                         -RNB*MD%D(3,NA,NB-1,0)

                    SXP(IA,IB)=SXP(IA,IB) &
                         +(OnA*MDDerivAx+OnB*MDDerivBx) &
                         *CA*CB*Ov*MD%D(2,MA,MB,0)*MD%D(3,NA,NB,0)
                    SYP(IA,IB)=SYP(IA,IB) &
                         +(OnA*MDDerivAy+OnB*MDDerivBy) &
                         *CA*CB*Ov*MD%D(1,LA,LB,0)*MD%D(3,NA,NB,0)
                    SZP(IA,IB)=SZP(IA,IB) &
                         +(OnA*MDDerivAz+OnB*MDDerivBz) &
                         *CA*CB*Ov*MD%D(2,MA,MB,0)*MD%D(1,LA,LB,0)
                 ENDDO
              ENDDO
           ENDDO 
        ENDDO
        IndexB=IndexB+(StopLB-StartLB+1)
     ENDDO
     IndexA=IndexA+(StopLA-StartLA+1)
  ENDDO

  SPVck(1)=  BlkTrace_2(Pair%NB,Pair%NA,W,SXP(1,1))
  SPVck(2)=  BlkTrace_2(Pair%NB,Pair%NA,W,SYP(1,1))
  SPVck(3)=  BlkTrace_2(Pair%NB,Pair%NA,W,SZP(1,1))

END FUNCTION GradSBlok
END MODULE
