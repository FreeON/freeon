MODULE TBlock
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  IMPLICIT NONE
  CONTAINS
  FUNCTION TBlok(BS,MD,Pair) RESULT(TVck)
    TYPE(BSET)                                :: BS
    TYPE(DBL_RNK4)                            :: MD
    TYPE(AtomPair)                            :: Pair
!
    REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB)   :: TVck
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)   :: TBlk
!
    REAL(DOUBLE),DIMENSION(0:BS%NASym,0:BS%NASym) :: TKx,TKy,TKz
!
    INTEGER                                   :: NBFA,NBFB,KA,KB,KK,AtA,AtB
    REAL(DOUBLE)                              :: Ax,Ay,Az,Bx,By,Bz,AB2 
    REAL(DOUBLE)                              :: Px,Py,Pz,PAx,PAy,PAz,PBx,PBy,PBz,ZetaA,ZetaB,  &
                                                 EtaAB,EtaIn,ZAB,XiAB,ExpAB,Ov,PiE32,CA,CB,     &
                                                 TBlkx,TBlky,TBlkz,DLA,DLB
    INTEGER                                   :: CFA,CFB,PFA,PFB,IndexA,IndexB,StartLA,StartLB, &
                                                 StopLA,StopLB,MaxLA,MaxLB,IA,IB,LMNA,LMNB,     &
                                                 LA,LB,MA,MB,NA,NB
!
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
    TBlk=Zero
!                 
    DO CFA=1,BS%NCFnc%I(KA)                       ! Loop over contracted function A 
       IndexA  = CFBlokDex(BS,CFA,KA)
       StartLA = BS%LStrt%I(CFA,KA)        
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       DO CFB=1,BS%NCFnc%I(KB)                    ! Loop over contracted function B
          IndexB  = CFBlokDex(BS,CFB,KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)   
          DO PFA=1,BS%NPFnc%I(CFA,KA)             ! Loops over primitives in 
             DO PFB=1,BS%NPFnc%I(CFB,KB)          ! contracted functions A and B
                ZetaA=BS%Expnt%D(PFA,CFA,KA)
                ZetaB=BS%Expnt%D(PFB,CFB,KB)
                EtaAB=ZetaA+ZetaB 
                EtaIn=One/EtaAB
                ZAB  =ZetaA*ZetaB
                XiAB =ZAB*EtaIn
!-------------------------------------------------------------
                IF(TestPrimPair(MaxLA+MaxLB,XiAB,AB2)) THEN
                   ExpAB = EXP(-XiAB*AB2)
                   PiE32=(Pi/EtaAB)**(1.5D0)
                   Ov=ExpAB*PiE32
!
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
                   DO LA=0,MaxLA
                      DLA=DBLE(LA)
                      DO LB=0,MaxLB
                         DLB=DBLE(LB)
                         TKx(LA,LB)=     DLA*DLB*MD%D(1,LA-1,LB-1,0) &
                              -ZetaA*Two*DLB*MD%D(1,LA+1,LB-1,0) &
                              -ZetaB*Two*DLA*MD%D(1,LA-1,LB+1,0) &
                              +     Four*ZAB*MD%D(1,LA+1,LB+1,0)
                         TKy(LA,LB)=     DLA*DLB*MD%D(2,LA-1,LB-1,0) &
                              -ZetaA*Two*DLB*MD%D(2,LA+1,LB-1,0) &
                              -ZetaB*Two*DLA*MD%D(2,LA-1,LB+1,0) &
                              +     Four*ZAB*MD%D(2,LA+1,LB+1,0)
                         TKz(LA,LB)=     DLA*DLB*MD%D(3,LA-1,LB-1,0) &
                              -ZetaA*Two*DLB*MD%D(3,LA+1,LB-1,0) &
                              -ZetaB*Two*DLA*MD%D(3,LA-1,LB+1,0) &
                              +     Four*ZAB*MD%D(3,LA+1,LB+1,0)
                      ENDDO
                   ENDDO
!
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
!
                         TBlkx = TKx(LA,LB)     *MD%D(2,MA,MB,0)*MD%D(3,NA,NB,0) 
                         TBlky = MD%D(1,LA,LB,0)*TKy(MA,MB)     *MD%D(3,NA,NB,0) 
                         TBlkz = MD%D(1,LA,LB,0)*MD%D(2,MA,MB,0)*TKz(NA,NB) 
                         TBlk(IA,IB)=TBlk(IA,IB)+Half*CA*CB*Ov*(TBlkx+TBlky+TBlkz)
!
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    TVck = BlockToVect(NBFA,NBFB,TBlk)
  END FUNCTION TBlok
END MODULE TBlock



















