MODULE GradTBlock
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  IMPLICIT NONE
CONTAINS
  FUNCTION GradTBlok(BS,MD,AtA,AtB,DiffAtm,P,Pair) RESULT(TPVck)
    TYPE(BSET),                 INTENT(IN)    :: BS
    TYPE(DBL_RNK4),             INTENT(INOUT) :: MD
    TYPE(AtomPair)                            :: Pair
!
    REAL(DOUBLE),DIMENSION(3)                 :: TPVck
    REAL(DOUBLE),DIMENSION(:)                 :: P
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)   :: TBlk,TXP,TYP,TZP
    REAL(DOUBLE), EXTERNAL            :: BlkTrace_2
!
    INTEGER                                   :: NBFA,NBFB,KA,KB,KK,AtA,AtB,DiffAtm
    REAL(DOUBLE)                              :: Ax,Ay,Az,Bx,By,Bz,AB2 
    REAL(DOUBLE)                              :: Px,Py,Pz,PAx,PAy,PAz,PBx,PBy,PBz,ZetaA,ZetaB,  &
                                                 EtaAB,EtaIn,ZAB,XiAB,ExpAB,Ov,PiE32,CA,CB, &
                                                 TKx,TKy,TKz,xTKxA,yTKyA,zTKzA,xTKxB,yTKyB, &
                                                 zTKzB,DTKxA,DTKyA,DTKzA,DTKxB,DTKyB,DTKzB, &
                                                 DLA,DLB,DMA,DMB,DNA,DNB,OnA,OnB

    INTEGER                                   :: CFA,CFB,PFA,PFB,IndexA,IndexB,StartLA,StartLB, &
                                                 StopLA,StopLB,MaxLA,MaxLB,IA,IB,LMNA,LMNB,     &
                                                 LA,LB,MA,MB,NA,NB
!
    TXP=Zero
    TYP=Zero
    TZP=Zero
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

    IndexA=0                  
    DO CFA=1,BS%NCFnc%I(KA)                       ! Loop over contracted function A
       IndexB=0
       StartLA=BS%LStrt%I(CFA,KA)        
       StopLA =BS%LStop%I(CFA,KA)
       MaxLA=BS%ASymm%I(2,CFA,KA)
       DO CFB=1,BS%NCFnc%I(KB)                    ! Loop over contracted function B
          StartLB=BS%LStrt%I(CFB,KB)
          StopLB =BS%LStop%I(CFB,KB)
          MaxLB=BS%ASymm%I(2,CFB,KB)
          
          DO PFA=1,BS%NPFnc%I(CFA,KA)             ! Loops over primitives in 
             DO PFB=1,BS%NPFnc%I(CFB,KB)          ! contracted functions A and B
                ZetaA=BS%Expnt%D(PFA,CFA,KA)
                ZetaB=BS%Expnt%D(PFB,CFB,KB)
                EtaAB=ZetaA+ZetaB 
                EtaIn=One/EtaAB
                ZAB  =ZetaA*ZetaB
                XiAB =ZAB*EtaIn
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

                      DLA=DBLE(LA)
                      DLB=DBLE(LB)
                      DMA=DBLE(MA)
                      DMB=DBLE(MB)
                      DNA=DBLE(NA)
                      DNB=DBLE(NB)

                      TKx=     DLA*DLB*MD%D(1,LA-1,LB-1,0) &
                           -ZetaA*Two*DLB*MD%D(1,LA+1,LB-1,0) &
                           -ZetaB*Two*DLA*MD%D(1,LA-1,LB+1,0) &
                           +     Four*ZetaA*ZetaB*MD%D(1,LA+1,LB+1,0)
                      TKy=     DMA*DMB*MD%D(2,MA-1,MB-1,0) &
                           -ZetaA*Two*DMB*MD%D(2,MA+1,MB-1,0) &
                           -ZetaB*Two*DMA*MD%D(2,MA-1,MB+1,0) &
                           +     Four*ZetaA*ZetaB*MD%D(2,MA+1,MB+1,0)
                      TKz=     DNA*DNB*MD%D(3,NA-1,NB-1,0) &
                           -ZetaA*Two*DNB*MD%D(3,NA+1,NB-1,0) &
                           -ZetaB*Two*DNA*MD%D(3,NA-1,NB+1,0) &
                           +     Four*ZetaA*ZetaB*MD%D(3,NA+1,NB+1,0)

                      xTKxA=     DLA*DLB*MD%D(1,LA,LB-1,0) &
                           -ZetaA*Two*DLB*MD%D(1,LA+2,LB-1,0) &
                           -ZetaB*Two*DLA*MD%D(1,LA,LB+1,0) &
                           +     Four*ZetaA*ZetaB*MD%D(1,LA+2,LB+1,0)
                      yTKyA=     DMA*DMB*MD%D(2,MA,MB-1,0) &
                           -ZetaA*Two*DMB*MD%D(2,MA+2,MB-1,0) &
                           -ZetaB*Two*DMA*MD%D(2,MA,MB+1,0) &
                           +     Four*ZetaA*ZetaB*MD%D(2,MA+2,MB+1,0)
                      zTKzA=     DNA*DNB*MD%D(3,NA,NB-1,0) &
                           -ZetaA*Two*DNB*MD%D(3,NA+2,NB-1,0) &
                           -ZetaB*Two*DNA*MD%D(3,NA,NB+1,0) &
                           +     Four*ZetaA*ZetaB*MD%D(3,NA+2,NB+1,0)

                      xTKxB=     DLA*DLB*MD%D(1,LA-1,LB,0) &
                           -ZetaA*Two*DLB*MD%D(1,LA+1,LB,0) &
                           -ZetaB*Two*DLA*MD%D(1,LA-1,LB+2,0) &
                           +     Four*ZetaA*ZetaB*MD%D(1,LA+1,LB+2,0)
                      yTKyB=     DMA*DMB*MD%D(2,MA-1,MB,0) &
                           -ZetaA*Two*DMB*MD%D(2,MA+1,MB,0) &
                           -ZetaB*Two*DMA*MD%D(2,MA-1,MB+2,0) &
                           +     Four*ZetaA*ZetaB*MD%D(2,MA+1,MB+2,0)
                      zTKzB=     DNA*DNB*MD%D(3,NA-1,NB,0) &
                           -ZetaA*Two*DNB*MD%D(3,NA+1,NB,0) &
                           -ZetaB*Two*DNA*MD%D(3,NA-1,NB+2,0) &
                           +     Four*ZetaA*ZetaB*MD%D(3,NA+1,NB+2,0)

                      DTKxA=DLA*DLB*(DLA-One)*MD%D(1,LA-2,LB-1,0) &
                           -ZetaA*Two*DLB*(DLA+One)*MD%D(1,LA,LB-1,0) &
                           -ZetaB*Two*DLA*(DLA-One)*MD%D(1,LA-2,LB+1,0) &
                           +     Four*ZetaA*ZetaB*(DLA+One)*MD%D(1,LA,LB+1,0)

                      DTKyA=DMA*DMB*(DMA-One)*MD%D(2,MA-2,MB-1,0) &
                           -ZetaA*Two*DMB*(DMA+One)*MD%D(2,MA,MB-1,0) &
                           -ZetaB*Two*DMA*(DMA-One)*MD%D(2,MA-2,MB+1,0) &
                           +     Four*ZetaA*ZetaB*(DMA+One)*MD%D(2,MA,MB+1,0)

                      DTKzA=DNA*DNB*(DNA-One)*MD%D(3,NA-2,NB-1,0) &
                           -ZetaA*Two*DNB*(DNA+One)*MD%D(3,NA,NB-1,0) &
                           -ZetaB*Two*DNA*(DNA-One)*MD%D(3,NA-2,NB+1,0) &
                           +     Four*ZetaA*ZetaB*(DNA+One)*MD%D(3,NA,NB+1,0)

                      DTKxB=DLA*DLB*(DLB-One)*MD%D(1,LA-1,LB-2,0) &
                           -ZetaA*Two*DLB*(DLB-One)*MD%D(1,LA+1,LB-2,0) &
                           -ZetaB*Two*DLA*(DLB+One)*MD%D(1,LA-1,LB,0) &
                           +     Four*ZetaA*ZetaB*(DLB+One)*MD%D(1,LA+1,LB,0)

                      DTKyB=DMA*DMB*(DMB-One)*MD%D(2,MA-1,MB-2,0) &
                           -ZetaA*Two*DMB*(DMB-One)*MD%D(2,MA+1,MB-2,0) &
                           -ZetaB*Two*DMA*(DMB+One)*MD%D(2,MA-1,MB,0) &
                           +     Four*ZetaA*ZetaB*(DMB+One)*MD%D(2,MA+1,MB,0)

                      DTKzB=DNA*DNB*(DNB-One)*MD%D(3,NA-1,NB-2,0) &
                           -ZetaA*Two*DNB*(DNB-One)*MD%D(3,NA+1,NB-2,0) &
                           -ZetaB*Two*DNA*(DNB+One)*MD%D(3,NA-1,NB,0) &
                           +     Four*ZetaA*ZetaB*(DNB+One)*MD%D(3,NA+1,NB,0)

                      TXP(IA,IB)=TXP(IA,IB) - Half*CA*CB*Ov* &
                           (OnA*((DTKxA-Two*ZetaA*xTKxA)*     &
                           MD%D(2,MA,MB,0)*MD%D(3,NA,NB,0) + &
                           (DLA*MD%D(1,LA-1,LB,0) - Two*ZetaA* &
                           MD%D(1,LA+1,LB,0))*(TKy*MD%D(3,NA,NB,0)+ &
                           TKz*MD%D(2,MA,MB,0))) + &
                           OnB*((DTKxB-Two*ZetaB*xTKxB)*     &
                           MD%D(2,MA,MB,0)*MD%D(3,NA,NB,0) + &
                           (DLB*MD%D(1,LA,LB-1,0) - Two*ZetaB* &
                           MD%D(1,LA,LB+1,0))*(TKy*MD%D(3,NA,NB,0)+ &
                           TKz*MD%D(2,MA,MB,0))))

                      TYP(IA,IB)=TYP(IA,IB) - Half*CA*CB*Ov* &
                           (OnA*((DTKyA-Two*ZetaA*yTKyA)*     &
                           MD%D(1,LA,LB,0)*MD%D(3,NA,NB,0) + &
                           (DMA*MD%D(2,MA-1,MB,0) - Two*ZetaA* &
                           MD%D(2,MA+1,MB,0))*(TKx*MD%D(3,NA,NB,0)+ &
                           TKz*MD%D(1,LA,LB,0))) + &
                           OnB*((DTKyB-Two*ZetaB*yTKyB)*     &
                           MD%D(1,LA,LB,0)*MD%D(3,NA,NB,0) + &
                           (DMB*MD%D(2,MA,MB-1,0) - Two*ZetaB* &
                           MD%D(2,MA,MB+1,0))*(TKx*MD%D(3,NA,NB,0)+ &
                           TKz*MD%D(1,LA,LB,0))))

                      TZP(IA,IB)=TZP(IA,IB) - Half*CA*CB*Ov* &
                           (OnA*((DTKzA-Two*ZetaA*zTKzA)*     &
                           MD%D(1,LA,LB,0)*MD%D(2,MA,MB,0) + &
                           (DNA*MD%D(3,NA-1,NB,0) - Two*ZetaA* &
                           MD%D(3,NA+1,NB,0))*(TKx*MD%D(2,MA,MB,0)+ &
                           TKy*MD%D(1,LA,LB,0))) + &
                           OnB*((DTKzB-Two*ZetaB*zTKzB)*     &
                           MD%D(1,LA,LB,0)*MD%D(2,MA,MB,0) + &
                           (DNB*MD%D(3,NA,NB-1,0) - Two*ZetaB* &
                           MD%D(3,NA,NB+1,0))*(TKx*MD%D(2,MA,MB,0)+ &
                           TKy*MD%D(1,LA,LB,0))))
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          IndexB=IndexB+(StopLB-StartLB+1)
       ENDDO
       IndexA=IndexA+(StopLA-StartLA+1)
    ENDDO

    TPVck(1)=  BlkTrace_2(Pair%NB,Pair%NA,P,TXP(1,1))
    TPVck(2)=  BlkTrace_2(Pair%NB,Pair%NA,P,TYP(1,1))
    TPVck(3)=  BlkTrace_2(Pair%NB,Pair%NA,P,TZP(1,1))

  END FUNCTION GradTBlok
END MODULE GradTBlock



















