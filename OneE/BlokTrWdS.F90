!    BLOKWISE ACCUMULATION OF Tr{W_(A,B).dS^T_(A,B)/dA}
!    Author: Matt Challacombe
!------------------------------------------------------------------------------
MODULE BlokTrWdS
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE LinAlg
  USE AtomPairs
  USE BraBloks
  IMPLICIT NONE
  CONTAINS
     FUNCTION TrWdS(BS,Pair,W) RESULT(Vck)
        TYPE(BSET)                                :: BS
        TYPE(AtomPair)                            :: Pair
        TYPE(PrimPair)                            :: Prim
        REAL(DOUBLE),DIMENSION(3)                 :: Vck
        REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)   :: W
        REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB,3) :: dS
        INTEGER                                   :: NBFA,NBFB,KA,KB,KK
        REAL(DOUBLE)                              :: PiE32,Amp
        INTEGER                                   :: CFA,CFB,PFA,PFB,IndexA,IndexB,StartLA,StartLB, &
                                                      StopLA,StopLB,MaxLA,MaxLB,IA,IB,LMNA,LMNB,    &
                                                      LA,LB,MA,MB,NA,NB,K,KX,KY,KZ
        REAL(DOUBLE), EXTERNAL                    :: BlkTrace_2 
!---------------------------------------------------------------------------------------------------------
    IF(Pair%SameAtom)THEN
       Vck=Zero
       RETURN
    ENDIF 
!--------------------------------------------------------------------------------------------------
    Prim%A=Pair%A
    Prim%B=Pair%B
    Prim%KA=Pair%KA
    Prim%KB=Pair%KB
    Prim%AB2=Pair%AB2
!----------------------------------
    KA=Prim%KA
    KB=Prim%KB
    dS=Zero
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
          IF(TestPrimPair(Prim%Ell,Prim%Xi,Prim%AB2)) THEN
             Prim%PFA=PFA
             Prim%PFB=PFB
             K=K+1
             Amp=SetBraBlok(Prim,BS,Gradients_O=Pair%SameAtom)
!-------------------------------------------------------------------
             PiE32=(Pi/Prim%Zeta)**(1.5D0) 
             IB=IndexB
             DO LMNB=StartLB,StopLB
                IB=IB+1
                IA=IndexA
                DO LMNA=StartLA,StopLA
                   IA=IA+1
                   DO K=1,3
                     dS(IA,IB,K)=dS(IA,IB,K)+dHGBra%D(1,IA,IB,K)*PiE32
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDDO
       ENDDO
    ENDDO
    ENDDO
!    WRITE(*,*)'=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-'
!    PrintFlags%Fmt=DEBUG_DBLSTYLE
!    CALL  Print_DBL_Rank2A(dS(:,:,1),'dSx',Unit_O=6)
!    CALL  Print_DBL_Rank2A(W,'W',Unit_O=6)
!    CALL  Print_DBL_Rank2A(MATMUL(W,dS(:,:,3)),'W.dS',Unit_O=6)
    DO K=1,3
       Vck(K)=BlkTrace_2(Pair%NA,Pair%NB,W,TRANSPOSE(dS(:,:,K)))
    ENDDO
!   WRITE(*,*)' SBlkFrc3 = ',Vck(3)
  END FUNCTION TrWdS
END MODULE BlokTrWdS
