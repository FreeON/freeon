MODULE RhoBlok
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  USE BraBloks
  USE RhoTools
  IMPLICIT NONE
  CONTAINS
!--------------------------------------------------------------
! Count the Number of Primative basis funtions contained in atom A and B
!--------------------------------------------------------------
  SUBROUTINE PrimCount(BS,Pair,NDist,NCoef)
    TYPE(BSET)                 :: BS
    TYPE(AtomPair)             :: Pair
!
    INTEGER                    :: KA,KB,CFA,CFB,PFA,PFB,MaxLA,MaxLB,NDist,NCoef
    REAL(DOUBLE)               :: AB2,ZetaA,ZetaB,ZetaAB,XiAB,ExpAB
    LOGICAL                    :: AEQB
!      
    KA   = Pair%KA
    KB   = Pair%KB  
    AB2  = Pair%AB2 
    AEQB = Pair%SameAtom
!
    DO CFA=1,BS%NCFnc%I(KA)                       ! Loop over contracted function A 
       DO CFB=1,BS%NCFnc%I(KB)                    ! Loop over contracted function B          
          DO PFA=1,BS%NPFnc%I(CFA,KA)             ! Loops over primitives in A
             DO PFB=1,BS%NPFnc%I(CFB,KB)          ! Loops over primitives in A
                ZetaA = BS%Expnt%D(PFA,CFA,KA)
                ZetaB = BS%Expnt%D(PFB,CFB,KB)
                ZetaAB= ZetaA+ZetaB 
                XiAB  = ZetaA*ZetaB/ZetaAB
                IF(TestPrimPair(XiAB,AB2))THEN
                   MaxLA   = BS%ASymm%I(2,CFA,KA)
                   MaxLB   = BS%ASymm%I(2,CFB,KB)
                   NDist   = NDist + 1
                   NCoef   = NCoef + LHGTF(MaxLA+MaxLB)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
!
  END SUBROUTINE PrimCount
!--------------------------------------------------------------
! Calculate the Primative Distributions Generated from atom A and B
!--------------------------------------------------------------
  SUBROUTINE RhoBlk(BS,Dmat,Pair,NDist,NCoef,Rho)
    TYPE(BSET)                              :: BS
    TYPE(AtomPair)                          :: Pair
    TYPE(PrimPair)                          :: Prim
    TYPE(HGRho_new)                         :: Rho
    INTEGER                                 :: NDist,NCoef
!
    REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB) :: Dmat
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB) :: DD
    INTEGER                                 :: KA,KB,NBFA,NBFB,CFA,CFB,PFA,PFB, &
                                               IndexA,StartLA,StopLA,MaxLA, &
                                               IndexB,StartLB,StopLB,MaxLB, &
                                               IE,OffCo,LMN,LMNA,LMNB, &
                                               IA,IB,EllA,EllB,LAB,MAB,NAB, &
                                               LenKet,AtA,AtB
    REAL(DOUBLE)                            :: ZetaA,ZetaB,ZetaAB,ZetaIn,XiAB,ExpAB, &
                                               AB2,Ax,Ay,Az,Bx,By,Bz, &
                                               Px,Py,Pz,PAx,PAy,PAz,PBx,PBy,PBz, &
                                               Ex,Exy,Exyz,CA,CB,MaxAmp
!

    LOGICAL                                 :: AEQB
!
    KA   = Pair%KA
    KB   = Pair%KB  
    NBFA = Pair%NA
    NBFB = Pair%NB
    AEQB = Pair%SameAtom
!----------------------------------
    Prim%A=Pair%A
    Prim%B=Pair%B
    Prim%AB2=Pair%AB2
    Prim%KA=Pair%KA
    Prim%KB=Pair%KB
!
    DD = VectToBlock(NBFA,NBFB,Dmat)
    DO CFA=1,BS%NCFnc%I(KA)       
       IndexA  = CFBlokDex(BS,CFA,KA)
       StartLA = BS%LStrt%I(CFA,KA)        
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       DO CFB=1,BS%NCFnc%I(KB)                       
          IndexB  = CFBlokDex(BS,CFB,KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          Prim%CFA=CFA
          Prim%CFB=CFB
          Prim%Ell=MaxLA+MaxLB
          DO PFA=1,BS%NPFnc%I(CFA,KA)                ! 
             Prim%PFA=PFA 
             DO PFB=1,BS%NPFnc%I(CFB,KB) 
                Prim%PFB=PFB
                Prim%ZA=BS%Expnt%D(PFA,CFA,KA)
                Prim%ZB=BS%Expnt%D(PFB,CFB,KB)
                Prim%Zeta=Prim%ZA+Prim%ZB
                Prim%Xi=Prim%ZA*Prim%ZB/Prim%Zeta
                IF(TestPrimPair(Prim%Xi,Prim%AB2)) THEN
!--------------------------------------------------------------
!                  Set primitive values
!                  Primitive coefficients in a HG representation
!--------------------------------------------------------------
                   MaxAmp=SetBraBlok(Prim,BS)
                   LenKet = LHGTF(Prim%Ell)
!
!                  Update the Counters
!
                   NDist = NDist+1
                   NCoef = NCoef+LenKet
                   OffCo = NCoef-LenKet+1
!
!                  Store the distribution
!
                   Rho%Ell%I(NDist)   = Prim%Ell
                   Rho%Zeta%D(NDist)  = Prim%Zeta
                   Rho%Qx%D(NDist)    = Prim%P(1)
                   Rho%Qy%D(NDist)    = Prim%P(2)
                   Rho%Qz%D(NDist)    = Prim%P(3)
                   Rho%OffCo%I(NDist) = OffCo
!
!                  Calculate and Store the Coefficients of the Distribution
!
                   Rho%Co%D(OffCo:OffCo+LenKet-1) = Zero 
                   IA=IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      EllA = BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                         EllB = BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)
                         DO LMN=1,LHGTF(EllA+EllB)
                            Rho%Co%D(OffCo+LMN-1)=Rho%Co%D(OffCo+LMN-1)+HGBra%D(LMN,IA,IB)*DD(IA,IB)
                         ENDDO
                      ENDDO
                   ENDDO                  
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE RhoBlk
!--------------------------------------------------------------
! Add nuclear centers to the density
!--------------------------------------------------------------
  SUBROUTINE AddNukes(GM,Rho,IOff_O)
    TYPE(CRDS)                    :: GM
    TYPE(HGRho)                   :: Rho
    INTEGER                       :: IA,Iq,Ir,IOffQ,IOffR
    INTEGER,DIMENSION(2),OPTIONAL :: IOff_O
    REAL(DOUBLE)                  :: DDelta
!---------------------------------------------------
! Initialize the Nuclear centers and Constants
    DDelta = Half*(Rho%Expt%D(Rho%NExpt)/Pi)**(ThreeHalves)
    IF(PRESENT(IOff_O))THEN
       IOffQ=IOff_O(1)
       IOffR=IOff_O(2)
    ELSE
       IOffQ=Rho%OffQ%I(Rho%NExpt)
       IOffR=Rho%OffR%I(Rho%NExpt)
    ENDIF

    DO IA = 1,GM%Natms       
       Iq=IOffQ+IA
       Ir=IOffR+IA
       Rho%Qx%D(Iq)= GM%Carts%D(1,IA)
       Rho%Qy%D(Iq)= GM%Carts%D(2,IA)
       Rho%Qz%D(Iq)= GM%Carts%D(3,IA)
       Rho%Co%D(Ir)=-GM%AtNum%D(IA)*DDelta
    ENDDO
  END SUBROUTINE AddNukes

!
END MODULE RhoBlok
