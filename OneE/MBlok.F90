MODULE MBlock
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  IMPLICIT NONE
CONTAINS
  FUNCTION DBlok(BS,MD,Pair,IXYZ,COrig) RESULT(DVck)
    TYPE(BSET)                                :: BS
    TYPE(DBL_RNK4)                            :: MD
    TYPE(AtomPair)                            :: Pair
!
    REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB)   :: DVck
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)   :: DBlk
!
    REAL(DOUBLE),DIMENSION(0:BS%NASym,0:BS%NASym) :: TKx,TKy,TKz
!
    INTEGER                                   :: NBFA,NBFB,KA,KB,KK,AtA,AtB
    REAL(DOUBLE)                              :: Ax,Ay,Az,Bx,By,Bz,AB2 
    REAL(DOUBLE)                              :: Px,Py,Pz,PAx,PAy,PAz,PBx,PBy,PBz,ZetaA,ZetaB,  &
         &                                       EtaAB,EtaIn,ZAB,XiAB,ExpAB,Ov,PiE32,CA,CB,     &
         &                                       TBlkx,TBlky,TBlkz,DLA,DLB
    INTEGER                                   :: CFA,CFB,PFA,PFB,IndexA,IndexB,StartLA,StartLB, &
         &                                       StopLA,StopLB,MaxLA,MaxLB,IA,IB,LMNA,LMNB,     &
         &                                       LA,LB,MA,MB,NA,NB
!vw--->
    INTEGER     , INTENT(IN)   :: IXYZ
    REAL(DOUBLE)               :: POx,POy,POz
    TYPE(DBL_VECT), INTENT(IN) :: COrig
!vw<---
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

    DBlk(:,:)=Zero
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
                !
                CALL MD2TRR(BS%NASym+1,-1,MaxLA,MaxLB,EtaAB,MD%D, &
                     &      PAx,PBx,PAy,PBy,PAz,PBz)
                !vw--->
                POx = Px-COrig%D(1)
                POy = Py-COrig%D(2)
                POz = Pz-COrig%D(3)
                !vw<---
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
                      !vw--->
                      SELECT CASE(IXYZ)
                      CASE(1)
                         DBlk(IA,IB)=DBlk(IA,IB) &
                              &     +CA*CB*Ov*((POx*MD%D(1,LA,LB,0)+MD%D(1,LA,LB,1))*MD%D(2,MA,MB,0)*MD%D(3,NA,NB,0))
                      CASE(2)
                         DBlk(IA,IB)=DBlk(IA,IB) &
                              &     +CA*CB*Ov*((POy*MD%D(2,MA,MB,0)+MD%D(2,MA,MB,1))*MD%D(1,LA,LB,0)*MD%D(3,NA,NB,0))
                      CASE(3)
                         DBlk(IA,IB)=DBlk(IA,IB) &
                              &     +CA*CB*Ov*((POz*MD%D(3,NA,NB,0)+MD%D(3,NA,NB,1))*MD%D(1,LA,LB,0)*MD%D(2,MA,MB,0))
                      END SELECT
                      !vw<---
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    DVck = BlockToVect(NBFA,NBFB,DBlk)
    !
  END FUNCTION DBlok
  !
END MODULE MBlock



















