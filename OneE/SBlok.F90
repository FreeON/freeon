MODULE OverlapBlock
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  USE BraBloks
  IMPLICIT NONE
  CONTAINS
  FUNCTION SBlok(BS,Pair) RESULT(SVck)
    TYPE(BSET)                              :: BS
    TYPE(AtomPair)                          :: Pair
    TYPE(PrimPair)                          :: Prim

    REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB) :: SVck

    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB) :: SBlk
    INTEGER                                 :: NBFA,NBFB,KA,KB,KK
    REAL(DOUBLE)                            :: Ax,Ay,Az,Bx,By,Bz,AB2 
    REAL(DOUBLE)                            :: PiE32,Amp
    INTEGER                                 :: CFA,CFB,PFA,PFB,IndexA,IndexB,StartLA,StartLB, &
                                               StopLA,StopLB,MaxLA,MaxLB,IA,IB,LMNA,LMNB,     &
                                               LA,LB,MA,MB,NA,NB
!--------------------------------------------------------------------------------------------------
    KA   = Pair%KA
    KB   = Pair%KB
    NBFA = Pair%NA
    NBFB = Pair%NB
    AB2=Pair%AB2
    SBlk=Zero
!----------------------------------
    Prim%A=Pair%A
    Prim%B=Pair%B
    Prim%AB2=Pair%AB2
    Prim%KA=Pair%KA
    Prim%KB=Pair%KB
!----------------------------------
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
!----------------------------------
       Prim%CFA=CFA
       Prim%CFB=CFB
       Prim%Ell=MaxLA+MaxLB
!----------------------------------
       DO PFA=1,BS%NPFnc%I(CFA,KA) 
       DO PFB=1,BS%NPFnc%I(CFB,KB)        
!----------------------------------
          Prim%ZA=BS%Expnt%D(PFA,CFA,KA)
          Prim%ZB=BS%Expnt%D(PFB,CFB,KB)
          Prim%Zeta=Prim%ZA+Prim%ZB
          Prim%Xi=Prim%ZA*Prim%ZB/Prim%Zeta
          IF(TestPrimPair(Prim%Xi,Prim%AB2))THEN
             Prim%PFA=PFA
             Prim%PFB=PFB
             Amp=SetBraBlok(Prim,BS)
!----------------------------------
             PiE32=(Pi/Prim%Zeta)**(1.5D0) 
             IB=IndexB
             DO LMNB=StartLB,StopLB
                IB=IB+1
                IA=IndexA
                DO LMNA=StartLA,StopLA
                   IA=IA+1
                   SBlk(IA,IB)=SBlk(IA,IB)+HGBra%D(1,IA,IB)*PiE32
                ENDDO
             ENDDO
          ENDIF
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    SVck = BlockToVect(Pair%NA,Pair%NB,SBlk)
  END FUNCTION SBlok
END MODULE OverlapBlock
