!    MODULE FOR GENERIC THE McMurchie Davidson APPROACH TO COMPUTATION OF
!    HERMITE GAUSSIAN ERIS VIA RECURENCE RELATIONS
!    Author: Matt Challacombe 
!------------------------------------------------------------------------------
MODULE McMurchie
   USE DerivedTypes
   USE GlobalScalars
   USE MemMan
   IMPLICIT NONE
   CONTAINS 
!-----------------------------------------------------------     
!     McMurchie-Davidson 2-term recurence relation
!
      SUBROUTINE MD2TRR(NASym,MD0,MaxLA,MaxLB,EtaAB,MD, &
                        PAx,PBx,PAy,PBy,PAz,PBz) 
         REAL(DOUBLE), INTENT(IN)  :: EtaAB,PAx,PBx,PAy,PBy,PAz,PBz
         REAL(DOUBLE)              :: RL1,TwoZ
         INTEGER,      INTENT(IN)  :: NASym,MD0,MaxLA,MaxLB
         REAL(DOUBLE), INTENT(OUT) :: MD(3,MD0:NASym,MD0:NASym,MD0:2*NASym)
         INTEGER                   :: LTot,LA,LB,LAB
         LTot=MaxLA+MaxLB
         DO LAB=MD0,LTot
            DO LB=MD0,MaxLB
               DO LA=MD0,MaxLA
                  MD(1,LA,LB,LAB)=Zero
                  MD(2,LA,LB,LAB)=Zero
                  MD(3,LA,LB,LAB)=Zero
               ENDDO
            ENDDO
         ENDDO                                                          
         MD(1,0,0,0)=One
         MD(2,0,0,0)=One
         MD(3,0,0,0)=One
         IF(LTot.EQ.0)RETURN
         TwoZ=Half/EtaAB
         DO LA=1,MaxLA
            MD(1,LA,0,0)=PAx*MD(1,LA-1,0,0)+MD(1,LA-1,0,1)
            MD(2,LA,0,0)=PAy*MD(2,LA-1,0,0)+MD(2,LA-1,0,1)
            MD(3,LA,0,0)=PAz*MD(3,LA-1,0,0)+MD(3,LA-1,0,1)
            DO LAB=1,LA-1           
               RL1=DBLE(LAB+1)                                      
               MD(1,LA,0,LAB)=TwoZ*MD(1,LA-1,0,LAB-1) &
                             + PAx*MD(1,LA-1,0,LAB  ) &           
                             + RL1*MD(1,LA-1,0,LAB+1)                         
               MD(2,LA,0,LAB)=TwoZ*MD(2,LA-1,0,LAB-1) &
                             + PAy*MD(2,LA-1,0,LAB  ) &            
                             + RL1*MD(2,LA-1,0,LAB+1)                         
               MD(3,LA,0,LAB)=TwoZ*MD(3,LA-1,0,LAB-1) &
                             + PAz*MD(3,LA-1,0,LAB  ) &            
                             + RL1*MD(3,LA-1,0,LAB+1)                         
            ENDDO
            MD(1,LA,0,LA)=TwoZ*MD(1,LA-1,0,LA-1)+PAx*MD(1,LA-1,0,LA)               
            MD(2,LA,0,LA)=TwoZ*MD(2,LA-1,0,LA-1)+PAy*MD(2,LA-1,0,LA)               
            MD(3,LA,0,LA)=TwoZ*MD(3,LA-1,0,LA-1)+PAz*MD(3,LA-1,0,LA)               
         ENDDO 
         DO LB=1,MaxLB
            DO LA=0,MaxLA  
               MD(1,LA,LB,0)=PBx*MD(1,LA,LB-1,0)+MD(1,LA,LB-1,1)
               MD(2,LA,LB,0)=PBy*MD(2,LA,LB-1,0)+MD(2,LA,LB-1,1)
               MD(3,LA,LB,0)=PBz*MD(3,LA,LB-1,0)+MD(3,LA,LB-1,1)
               DO LAB=1,LTot-1
                  RL1=DBLE(LAB+1)
                  MD(1,LA,LB,LAB)=TwoZ*MD(1,LA,LB-1,LAB-1) &
                                 + PBx*MD(1,LA,LB-1,LAB  ) &        
                                 + RL1*MD(1,LA,LB-1,LAB+1)
                  MD(2,LA,LB,LAB)=TwoZ*MD(2,LA,LB-1,LAB-1) &
                                 + PBy*MD(2,LA,LB-1,LAB  ) &         
                                 + RL1*MD(2,LA,LB-1,LAB+1)
                  MD(3,LA,LB,LAB)=TwoZ*MD(3,LA,LB-1,LAB-1) & 
                                 + PBz*MD(3,LA,LB-1,LAB  ) &        
                                 + RL1*MD(3,LA,LB-1,LAB+1)
               ENDDO
               MD(1,LA,LB,LTot)=TwoZ*MD(1,LA,LB-1,LAB-1)+PBx*MD(1,LA,LB-1,LAB)
               MD(2,LA,LB,LTot)=TwoZ*MD(2,LA,LB-1,LAB-1)+PBy*MD(2,LA,LB-1,LAB)
               MD(3,LA,LB,LTot)=TwoZ*MD(3,LA,LB-1,LAB-1)+PBz*MD(3,LA,LB-1,LAB)
            ENDDO
          ENDDO

!         DO LAB=0,LTot
!            DO LB=MD0,MaxLB
!               DO LA=MD0,MaxLA
!                  WRITE(*,22)LA,LB,LAB,MD(:,LA,LB,LAB)
!               22 FORMAT(3(I3,","),3(D14.6,","))
!               ENDDO
!            ENDDO
!         ENDDO                                                          
      END SUBROUTINE MD2TRR
!-----------------------------------------------------------     
!     McMurchie-Davidson 3-term recurence relation
!
      SUBROUTINE MD3TRR(MaxL,LTot,R,AuxR,Upq,PQx,PQy,PQz)
         INTEGER,                                INTENT(IN)    :: LTot,MaxL
         REAL(DOUBLE), DIMENSION(0:LTot),        INTENT(IN)    :: AuxR 
         REAL(DOUBLE), DIMENSION(0:MaxL,0:MaxL, &
                                 0:MaxL,0:MaxL), INTENT(INOUT) :: R
         REAL(DOUBLE),                           INTENT(IN)    :: PQx,PQy,PQz,Upq
         INTEGER                                               :: J,J1,L,L1,L2, &
                                                                  M,M1,M2,N,N1,N2
         REAL(DOUBLE)                                          :: REM1,REN1,REL1
         DO J=0,LTot
            R(0,0,0,J)=Upq*AuxR(J)
         ENDDO
         DO J=0,LTot-1
            J1=J+1
            R(0,0,1,J)=R(0,0,0,J1)*PQz
         ENDDO
         DO N=2,LTot
            N1=N-1
            N2=N-2
            REN1=DBLE(N1)
            DO J=0,LTot-N
               J1=J+1
               R(0,0,N,J)=R(0,0,N1,J1)*PQz+R(0,0,N2,J1)*REN1
            ENDDO
         ENDDO
         DO N=0,LTot
            DO J=0,LTot-N-1
               J1=J+1
               R(0,1,N,J)=R(0,0,N,J1)*PQy
            ENDDO
            DO M=2,LTot-N
               M1=M-1
               M2=M-2
               REM1=DBLE(M1)
               DO J=0,LTot-N-M
                  J1=J+1
                  R(0,M,N,J)=R(0,M1,N,J1)*PQy+R(0,M2,N,J1)*REM1
               ENDDO
            ENDDO
         ENDDO
         DO N=0,LTot
            DO M=0,LTot-N
               DO J=0,LTot-N-M-1
                  J1=J+1
                  R(1,M,N,J)=R(0,M,N,J1)*PQx
               ENDDO
               DO L=2,LTot-N-M
                  L1=L-1
                  L2=L-2
                  REL1=DBLE(L1)
                  DO J=0,LTot-N-M-L
                     J1=J+1
                     R(L,M,N,J)=R(L1,M,N,J1)*PQx+R(L2,M,N,J1)*REL1
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      END SUBROUTINE MD3TRR
!-------------------------------------------------------------
!     Compute the auxiliary integrals R_{0,0,0,J}(T)
!
      SUBROUTINE AuxInts(MaxL,LTot,AuxR,Omega,T)
         REAL(DOUBLE),                  INTENT(IN)  :: Omega,T
         INTEGER,                       INTENT(IN)  :: MaxL,LTot 
         REAL(DOUBLE),DIMENSION(0:MaxL),INTENT(OUT) :: AuxR
         REAL(DOUBLE),PARAMETER                     :: Switch=26.0D0
         INTEGER,PARAMETER                          :: LPlus=50
         INTEGER,PARAMETER                          :: L2=12+LPlus
         REAL(DOUBLE),DIMENSION(0:L2)               :: F
         REAL(DOUBLE)                               :: SqrtT,ET,OneOvT,FJ,TwoT, &
                                                       OmegaJ,TwoO
         INTEGER                                    :: J
!---------------------------------------------------------------------------------
!        Compute the incomplete gamma functions F_j(T)
!
         IF(T==Zero)THEN
            OmegaJ=One
            TwoO=-Two*Omega
            DO J=0,LTot
               AuxR(J)=OmegaJ/DBLE(2*J+1)
               OmegaJ=TwoO*OmegaJ
            ENDDO
            RETURN            
         ELSEIF(T.LT.Switch) THEN
!---------------------------------------------------
!           Downward recursion:
!           F_{j}(T)=(2*T*F_{j+1}+Exp[-T])/(2*j+1)
!
            TwoT=Two*T
            ET=EXP(-T)
            FJ=Zero
            DO J=LTot+LPlus,0,-1
               F(J)=FJ
               FJ=(TwoT*F(J)+ET)/(Two*DBLE(J)-One)
            ENDDO
         ELSE
!----------------------------------------------------
!           Multipole approx and upward recursion
!
            SqrtT=SQRT(T)
            OneOvT=One/T
            F(0)=SqrtPi/(Two*SqrtT) 
            DO J=1,LTot
               F(J)=F(J-1)*(DBLE(J)-Half)*OneOvT
            ENDDO
         ENDIF
!------------------------------------------------------
!        Generate the auxiliary integrals 
!        R_{000j}=(-2*omega)^j F_{j}(T)
!
         OmegaJ=One
         TwoO=-Two*Omega
         DO J=0,LTot
            AuxR(J)=OmegaJ*F(J)
            OmegaJ=TwoO*OmegaJ
         ENDDO
      END SUBROUTINE AuxInts




!-------------------------------------------------------------
!     Compute the auxiliary integrals R_{0,0,0,J}(T)
!
      SUBROUTINE ErrInts(MaxL,LTot,ErrR,Omega,T)
         REAL(DOUBLE),                  INTENT(IN)  :: Omega,T
         INTEGER,                       INTENT(IN)  :: MaxL,LTot 
         REAL(DOUBLE),DIMENSION(0:MaxL),INTENT(OUT) :: ErrR
         REAL(DOUBLE),PARAMETER                     :: Switch=35.0D0
         INTEGER,PARAMETER                          :: LPlus=150
         INTEGER,PARAMETER                          :: L2=12+LPlus
         REAL(DOUBLE),DIMENSION(0:MaxL)             :: E
         REAL(DOUBLE),DIMENSION(0:L2)               :: M,F
         REAL(DOUBLE)                               :: SqrtT,ET,OneOvT,FJ,TwoT, &
                                                       OmegaJ,TwoO
         INTEGER                                    :: J
!---------------------------------------------------------------------------------
!        Compute the incomplete gamma functions F_j(T)
!
         IF(T==Zero) &
            CALL Halt('Infinity in ErrInts')
         IF(T>Switch)THEN
            ErrR=Zero
            RETURN
         ENDIF
!---------------------------------------------------
!        Downward recursion:
!        F_{j}(T)=(2*T*F_{j+1}+Exp[-T])/(2*j+1)
         TwoT=Two*T
         ET=EXP(-T)
         FJ=Zero
         DO J=LTot+LPlus,0,-1
            F(J)=FJ
            FJ=(TwoT*F(J)+ET)/(Two*DBLE(J)-One)
         ENDDO
!----------------------------------------------------
!        Multipole approx and upward recursion
         SqrtT=SQRT(T)
         OneOvT=One/T
         M(0)=SqrtPi/(Two*SqrtT) 
         DO J=1,LTot
            M(J)=M(J-1)*(DBLE(J)-Half)*OneOvT
         ENDDO
!------------------------------------------------------
!        Generate the auxiliary error integrals 
!        R_{000j}=(-2*omega)^j [F_{j}(T)-M_{j}(T)]
         OmegaJ=One
         TwoO=-Two*Omega
         DO J=0,LTot
            ErrR(J)=OmegaJ*ABS(M(J)-F(J))
            OmegaJ=TwoO*OmegaJ
         ENDDO
      END SUBROUTINE ErrInts
!----------------------------------------------------------------------------------------------------------------
      SUBROUTINE PrintIntegrals(MaxL,Max2L,Max4L,MaxLMN,Max4LMN,BS,GM,pInt,iInt) 
        TYPE(CRDS) :: GM
        TYPE(BSET) :: BS
        INTEGER :: MaxL, Max2L, Max4L, MaxLMN, Max4LMN, pInt, iInt
        REAL(DOUBLE), DIMENSION(0:Max4L)                         :: AuxR 
        REAL(DOUBLE), DIMENSION(1:3,0:MaxL,0:MaxL,0:Max2L)       :: Eac,Ebd
        REAL(DOUBLE), DIMENSION(0:Max4L,0:Max4L,0:Max4L,0:Max4L) :: MDR
        REAL(DOUBLE), DIMENSION(Max4LMN)                         :: Integrals      
        INTEGER :: AtA,CFA,PFA,StartLA,StopLA,StrideA
        INTEGER :: AtB,CFB,PFB,StartLB,StopLB,StrideB
        INTEGER :: AtC,CFC,PFC,StartLC,StopLC,StrideC
        INTEGER :: AtD,CFD,PFD,StartLD,StopLD,StrideD
        INTEGER :: StrideAB,StrideCD

        INTEGER :: IndexC1,IndexC,IndexA1,IndexA
        INTEGER :: IndexB1,IndexB,IndexD1,IndexD
        INTEGER :: IndexK1,IndexK
        INTEGER :: KC,KA,KB,KD
        INTEGER :: MaxLA,MaxLC,MaxLB,MaxLD,LTot
        INTEGER :: LA,MA,NA,LB,MB,NB,LC,NC,MC,LD,MD,ND
        INTEGER :: I,IK,IA,IB,IC,ID
        INTEGER :: LMNA,LMNB,LMNC,LMND,LAC,MAC,NAC,LBD,MBD,NBD
        REAL(DOUBLE) :: Cx,Cy,Cz,Ax,Ay,Az
        REAL(DOUBLE) :: Bx,By,Bz,Dx,Dy,Dz
        REAL(DOUBLE) :: Px,Py,Pz,Qx,Qy,Qz
        REAL(DOUBLE) :: BDx,BDy,BDz,BD2,ACx,ACy,ACz,AC2
        REAL(DOUBLE) :: CDx,CDy,CDz,CD2,ABx,ABy,ABz,AB2
        REAL(DOUBLE) :: PAx,PAy,PAz,PCx,PCy,PCz
        REAL(DOUBLE) :: QBx,QBy,QBz,QDx,QDy,QDz
        REAL(DOUBLE) :: PQx,PQy,PQz,PQ2
        REAL(DOUBLE) :: ZetaA,ZetaB,ZetaC,ZetaD,EtaAC,RhoBD
        REAL(DOUBLE) :: Zac,Zbd,EtaIn,RhoIn,XiBD,XiAC,ExpAC,ExpBD,W,U
        REAL(DOUBLE) :: CA,CB,CC,CD,CCoAC,CCoBD

        IndexC1=0
        DO 101 AtC=1,NAtoms
           Cx=GM%Carts%D(1,AtC)
           Cy=GM%Carts%D(2,AtC)
           Cz=GM%Carts%D(3,AtC)
           KC=GM%AtTyp%I(AtC)
           IndexD1=0
           DO 102 AtD=1,NAtoms
              Dx=GM%Carts%D(1,AtD)
              Dy=GM%Carts%D(2,AtD)
              Dz=GM%Carts%D(3,AtD)
              KD=GM%AtTyp%I(AtD)
              CDx=Cx-Dx
              CDy=Cy-Dy
              CDz=Cz-Dz
              CD2=CDx*CDx+CDy*CDy+CDz*CDz
              IndexK1=0
              IndexA1=0
              DO 103 AtA=1,NAtoms
                 Ax=GM%Carts%D(1,AtA)
                 Ay=GM%Carts%D(2,AtA)
                 Az=GM%Carts%D(3,AtA)
                 KA=GM%AtTyp%I(AtA)
                 ACx=Ax-Cx
                 ACy=Ay-Cy
                 ACz=Az-Cz
                 AC2=ACx*ACx+ACy*ACy+ACz*ACz
                 IndexB1=0
                 DO 104 AtB=1,NAtoms
                    Bx=GM%Carts%D(1,AtB)
                    By=GM%Carts%D(2,AtB)
                    Bz=GM%Carts%D(3,AtB)
                    KB=GM%AtTyp%I(AtB)
                    ABx=Ax-Bx
                    ABy=Ay-By
                    ABz=Az-Bz
                    AB2=ABx*ABx+ABy*ABy+ABz*ABz
                    BDx=Bx-Dx
                    BDy=By-Dy
                    BDz=Bz-Dz
                    BD2=BDx*BDx+BDy*BDy+BDz*BDz
                    IndexC=IndexC1
                    DO 201 CFC=1,BS%NCFnc%I(KC)
                       StartLC=BS%LStrt%I(CFC,KC)
                       StopLC =BS%LStop%I(CFC,KC)
                       StrideC=StopLC-StartLC+1
                       MaxLC=BS%ASymm%I(2,CFC,KC)
                       IndexD=IndexD1
                       DO 202 CFD=1,BS%NCFnc%I(KD)
                          StartLD=BS%LStrt%I(CFD,KD)
                          StopLD =BS%LStop%I(CFD,KD)
                          StrideD=StopLD-StartLD+1
                          MaxLD=BS%ASymm%I(2,CFD,KD)
                          StrideCD=StrideC*StrideD
                          IndexK=IndexK1
                          IndexA=IndexA1
                          DO 203 CFA=1,BS%NCFnc%I(KA)
                             StartLA=BS%LStrt%I(CFA,KA)
                             StopLA =BS%LStop%I(CFA,KA)
                             StrideA=StopLA-StartLA+1
                             MaxLA=BS%ASymm%I(2,CFA,KA)
                             IndexB=IndexB1
                             DO 204 CFB=1,BS%NCFnc%I(KB)
                                StartLB=BS%LStrt%I(CFB,KB)
                                StopLB =BS%LStop%I(CFB,KB)
                                StrideB=StopLB-StartLB+1
                                MaxLB=BS%ASymm%I(2,CFB,KB)
                                StrideAB=StrideA*StrideB
                                Integrals=0.0D0
                                DO 300 PFA=1,BS%NPFnc%I(CFA,KA)
                                   DO 300 PFC=1,BS%NPFnc%I(CFC,KC)
                                      ZetaA=BS%Expnt%D(PFA,CFA,KA)
                                      ZetaC=BS%Expnt%D(PFC,CFC,KC)
                                      EtaAC=ZetaA+ZetaC
                                      ZAC  =ZetaA*ZetaC
                                      EtaIn=1.0D0/EtaAC
                                      XiAC =ZetaA*ZetaC*EtaIn
                                      ExpAC=DEXP(-XiAC*AC2)
                                      Px=(ZetaA*Ax+ZetaC*Cx)*EtaIn
                                      Py=(ZetaA*Ay+ZetaC*Cy)*EtaIn
                                      Pz=(ZetaA*Az+ZetaC*Cz)*EtaIn
                                      PAx=Px-Ax
                                      PAy=Py-Ay
                                      PAz=Pz-Az
                                      PCx=Px-Cx
                                      PCy=Py-Cy
                                      PCz=Pz-Cz
                                      ! Compute McMurchie Davidson E coefficients for HG Primitives
                                      CALL MD2TRR(BS%NASym,0,MaxLA,MaxLC,EtaAC,Eac,PAx,PCx,PAy,PCy,PAz,PCz)
                                      ! CALL MD2TRRx(NASym,0,MaxLA,MaxLC,EtaAC,PAx,PCx,PAy,PCy,PAz,PCz,EACx,EACy,EACz)
                                      DO 400 PFB=1,BS%NPFnc%I(CFB,KB)
                                         DO 400 PFD=1,BS%NPFnc%I(CFD,KD)
                                            ZetaB=BS%Expnt%D(PFB,CFB,KB)
                                            ZetaD=BS%Expnt%D(PFD,CFD,KD)
                                            RhoBD=ZetaB+ZetaD
                                            ZBD  =ZetaB*ZetaD
                                            RhoIn=1.0D0/RhoBD
                                            XiBD =ZetaB*ZetaD*RhoIn
                                            ExpBD=DEXP(-XiBD*BD2)
                                            Qx=(ZetaB*Bx+ZetaD*Dx)*RhoIn
                                            Qy=(ZetaB*By+ZetaD*Dy)*RhoIn
                                            Qz=(ZetaB*Bz+ZetaD*Dz)*RhoIn
                                            QBx=Qx-Bx
                                            QBy=Qy-By
                                            QBz=Qz-Bz
                                            QDx=Qx-Dx
                                            QDy=Qy-Dy
                                            QDz=Qz-Dz
                                            ! Compute McMurchie Davidson E coefficients for HG Primitives
                                            CALL MD2TRR(BS%NASym,0,MaxLB,MaxLD,RhoBD,Ebd,QBx,QDx,QBy,QDy,QBz,QDz)
                                            ! CALL MD2TRRx(NASym,0,QBx,QDx,QBy,QDy,QBz,QDz,EBDx,EBDy,EBDz)
                                            PQx=Px-Qx
                                            PQy=Py-Qy
                                            PQz=Pz-Qz
                                            PQ2=PQx*PQx+PQy*PQy+PQz*PQz
                                            W=EtaAC*RhoBD/(EtaAC+RhoBD)
                                            U=34.986836655249725693D0/(EtaAC*RhoBD*DSQRT(EtaAC+RhoBD))
                                            LTot=MaxLA+MaxLB+MaxLC+MaxLD
                                            CALL AuxInts(Max4L,LTot,AuxR,W,W*PQ2)
                                            CALL MD3TRR(Max4L,LTot,MDR,AuxR,U,PQx,PQy,PQz)
                                            I=0
                                            IA=IndexA
                                            DO 500 LMNA=StartLA,StopLA
                                               IA=IA+1
                                               LA=BS%LxDex%I(LMNA)
                                               MA=BS%LyDex%I(LMNA)
                                               NA=BS%LzDex%I(LMNA)
                                               CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                                               IB=IndexB
                                               DO 500 LMNB=StartLB,StopLB
                                                  IB=IB+1
                                                  LB=BS%LxDex%I(LMNB)
                                                  MB=BS%LyDex%I(LMNB)
                                                  NB=BS%LzDex%I(LMNB)
                                                  CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                                                  IC=IndexC
                                                  DO 500 LMNC=StartLC,StopLC
                                                     IC=IC+1
                                                     LC=BS%LxDex%I(LMNC)
                                                     MC=BS%LyDex%I(LMNC)
                                                     NC=BS%LzDex%I(LMNC)
                                                     CC=BS%CCoef%D(LMNC,PFC,CFC,KC)
                                                     ID=IndexD
                                                     DO 500 LMND=StartLD,StopLD
                                                        ID=ID+1
                                                        LD=BS%LxDex%I(LMND)
                                                        MD=BS%LyDex%I(LMND)
                                                        ND=BS%LzDex%I(LMND)
                                                        CD=BS%CCoef%D(LMND,PFD,CFD,KD)
                                                        CCoAC=CA*CC*ExpAC
                                                        CCoBD=CB*CD*ExpBD
                                                        I=I+1
                                                        DO LAC=0,LA+LC
                                                           DO MAC=0,MA+MC
                                                              DO NAC=0,NA+NC
                                                                 DO LBD=0,LB+LD
                                                                    DO MBD=0,MB+MD
                                                                       DO NBD=0,NB+ND
                                                                          Integrals(I)=Integrals(I)+CCoAC*CCoBD &
                                                                               *Eac(1,LA,LC,LAC)       &
                                                                               *Eac(2,MA,MC,MAC)       &
                                                                               *Eac(3,NA,NC,NAC)       &
                                                                               *Ebd(1,LB,LD,LBD)       &
                                                                               *Ebd(2,MB,MD,MBD)       &
                                                                               *Ebd(3,NB,ND,NBD)       &
                                                                               *(-1.D0)**(LBD+MBD+NBD)   &
                                                                               *MDR(LAC+LBD,MAC+MBD,NAC+NBD,0)
                                                                       ENDDO
                                                                    ENDDO
                                                                 ENDDO
                                                              ENDDO
                                                           ENDDO
                                                        ENDDO
500                                         CONTINUE
400                                   CONTINUE
300                             CONTINUE
                                I=0
                                IK=IndexK
                                IA=IndexA
                                DO LMNA=StartLA,StopLA
                                   IA=IA+1
                                   IB=IndexB
                                   DO LMNB=StartLB,StopLB
                                      IB=IB+1
                                      IC=IndexC
                                      DO LMNC=StartLC,StopLC
                                         IC=IC+1
                                         ID=IndexD
                                         DO LMND=StartLD,StopLD
                                            ID=ID+1
                                            I=I+1
                                            CALL iPrint(Integrals(I),IA,IC,IB,ID,pInt,iInt)
                                         ENDDO
                                      ENDDO
                                   ENDDO
                                ENDDO
                                IndexK=IndexK+StrideAB
                                IndexB=IndexB+StrideB
204                          CONTINUE
                             IndexA=IndexA+StrideA
203                       CONTINUE
                          IndexD=IndexD+StrideD
202                    CONTINUE
                       IndexC=IndexC+StrideC
201                 CONTINUE
                    IndexK1=IndexK1+BS%BFKnd%I(KA)*BS%BFKnd%I(KB)
                    IndexB1=IndexB1+BS%BFKnd%I(KB)
104              CONTINUE
                 IndexA1=IndexA1+BS%BFKnd%I(KA)
103           CONTINUE
              IndexD1=IndexD1+BS%BFKnd%I(KD)
102        CONTINUE
           IndexC1=IndexC1+BS%BFKnd%I(KC)
101     CONTINUE
END SUBROUTINE PrintIntegrals


SUBROUTINE iPrint(Int,I,J,K,L,Mode,iOut)
  IMPLICIT REAL*8 (a-h,o-z)
  IMPLICIT INTEGER (i-n)
  REAL*8  Int
  INTEGER E
  IF(ABS(Int)<1D-12)RETURN
  IF(Mode.EQ.1)THEN
     WRITE(iOut,101)I,J,K,L,Int
  ELSEIF(Mode.EQ.2)THEN
     WRITE(iOut,201)I,J,K,L,Int
  ELSE
     WRITE(iOut,301)I,J,K,L,FRACTION(Int),EXPONENT(Int)
  ENDIF
101 FORMAT(' (',I3,',',I3,'|',I3,',',I3,') = ',F20.16)
201 FORMAT(' (',I3,',',I3,'|',I3,',',I3,') = ',D22.16)
301 FORMAT(1x,'Int[[',I3,',',I3,',',I3,',',I3,']] = ',F19.16,'*2^(',I3,');')
END SUBROUTINE iPrint

!----------------------------------------------------------------------------------------------------------------
     SUBROUTINE IntegralsC(BS,GM,CS,Integrals) 
        TYPE(CRDS)     :: GM
        TYPE(BSET)     :: BS
        TYPE(CellSet)  :: CS
        TYPE(DBL_RNK4) :: Integrals
        INTEGER :: MaxL, Max2L, Max4L, MaxLMN, Max4LMN
        REAL(DOUBLE), DIMENSION(0:4*BS%NASym)                                        :: AuxR 
        REAL(DOUBLE), DIMENSION(1:3,0:BS%NASym,0:BS%NASym,0:2*BS%NASym)              :: Eac,Ebd
        REAL(DOUBLE), DIMENSION(0:4*BS%NASym,0:4*BS%NASym,0:4*BS%NASym,0:4*BS%NASym) :: MDR
!      
        INTEGER :: AtA,CFA,PFA,StartLA,StopLA,StrideA
        INTEGER :: AtB,CFB,PFB,StartLB,StopLB,StrideB
        INTEGER :: AtC,CFC,PFC,StartLC,StopLC,StrideC
        INTEGER :: AtD,CFD,PFD,StartLD,StopLD,StrideD
        INTEGER :: StrideAB,StrideCD
        INTEGER :: NC1,NC2,NC3,NC4

        INTEGER :: IndexC1,IndexC,IndexA1,IndexA
        INTEGER :: IndexB1,IndexB,IndexD1,IndexD
        INTEGER :: IndexK1,IndexK
        INTEGER :: KC,KA,KB,KD
        INTEGER :: MaxLA,MaxLC,MaxLB,MaxLD,LTot
        INTEGER :: LA,MA,NA,LB,MB,NB,LC,NC,MC,LD,MD,ND
        INTEGER :: I,IK,IA,IB,IC,ID
        INTEGER :: LMNA,LMNB,LMNC,LMND,LAC,MAC,NAC,LBD,MBD,NBD
        REAL(DOUBLE) :: Cx,Cy,Cz,Ax,Ay,Az
        REAL(DOUBLE) :: Bx,By,Bz,Dx,Dy,Dz
        REAL(DOUBLE) :: Px,Py,Pz,Qx,Qy,Qz
        REAL(DOUBLE) :: BDx,BDy,BDz,BD2,ACx,ACy,ACz,AC2
        REAL(DOUBLE) :: CDx,CDy,CDz,CD2,ABx,ABy,ABz,AB2
        REAL(DOUBLE) :: PAx,PAy,PAz,PCx,PCy,PCz
        REAL(DOUBLE) :: QBx,QBy,QBz,QDx,QDy,QDz
        REAL(DOUBLE) :: PQx,PQy,PQz,PQ2
        REAL(DOUBLE) :: ZetaA,ZetaB,ZetaC,ZetaD,EtaAC,RhoBD
        REAL(DOUBLE) :: Zac,Zbd,EtaIn,RhoIn,XiBD,XiAC,ExpAC,ExpBD,W,U
        REAL(DOUBLE) :: CA,CB,CC,CD,CCoAC,CCoBD
!
        MaxL    =   BS%NASym
        Max2L   = 2*BS%NASym
        Max4L   = 4*BS%NASym
        MaxLMN  = BS%LMNLen
        Max4LMN = BS%LMNLen**4
!
        WRITE(*,*) CS%CellCarts%D(1,7), CS%CellCarts%D(1,9)
!
        Integrals%D(:,:,:,:)=Zero
        I = 0
        IndexC1 = 0
        DO 101  AtC=1,NAtoms
        DO 1010 NC1=7,7!1,CS%NCells
           Cx=GM%Carts%D(1,AtC)+CS%CellCarts%D(1,NC1)
           Cy=GM%Carts%D(2,AtC)+CS%CellCarts%D(2,NC1)
           Cz=GM%Carts%D(3,AtC)+CS%CellCarts%D(3,NC1)
           KC=GM%AtTyp%I(AtC)
           IndexD1=0
           DO 102  AtD=1,NAtoms
           DO 1020 NC2 = 9,9!1,CS%NCells
              Dx=GM%Carts%D(1,AtD)+CS%CellCarts%D(1,NC2)
              Dy=GM%Carts%D(2,AtD)+CS%CellCarts%D(2,NC2)
              Dz=GM%Carts%D(3,AtD)+CS%CellCarts%D(3,NC2)
              KD=GM%AtTyp%I(AtD)
              CDx=Cx-Dx
              CDy=Cy-Dy
              CDz=Cz-Dz
              CD2=CDx*CDx+CDy*CDy+CDz*CDz
              IndexK1=0
              IndexA1=0
              DO 103 AtA=1,NAtoms
!              DO 1030 NC3=8,8!1,CS%NCells
                 Ax=GM%Carts%D(1,AtA)!+CS%CellCarts%D(1,NC3)
                 Ay=GM%Carts%D(2,AtA)!+CS%CellCarts%D(2,NC3)
                 Az=GM%Carts%D(3,AtA)!+CS%CellCarts%D(3,NC3)
                 KA=GM%AtTyp%I(AtA)
                 ACx=Ax-Cx
                 ACy=Ay-Cy
                 ACz=Az-Cz
                 AC2=ACx*ACx+ACy*ACy+ACz*ACz
                 IndexB1=0
                 DO 104  AtB = 1,NAtoms
!                 DO 1040 NC4=9,9!1,CS%NCells
                    Bx=GM%Carts%D(1,AtB)!+CS%CellCarts%D(1,NC4)
                    By=GM%Carts%D(2,AtB)!+CS%CellCarts%D(2,NC4)
                    Bz=GM%Carts%D(3,AtB)!+CS%CellCarts%D(3,NC4)
                    KB=GM%AtTyp%I(AtB)
                    ABx=Ax-Bx
                    ABy=Ay-By
                    ABz=Az-Bz
                    AB2=ABx*ABx+ABy*ABy+ABz*ABz
                    BDx=Bx-Dx
                    BDy=By-Dy
                    BDz=Bz-Dz
                    BD2=BDx*BDx+BDy*BDy+BDz*BDz
                    IndexC=IndexC1
                    DO 201 CFC=1,BS%NCFnc%I(KC)
                       StartLC=BS%LStrt%I(CFC,KC)
                       StopLC =BS%LStop%I(CFC,KC)
                       StrideC=StopLC-StartLC+1
                       MaxLC=BS%ASymm%I(2,CFC,KC)
                       IndexD=IndexD1
                       DO 202 CFD=1,BS%NCFnc%I(KD)
                          StartLD=BS%LStrt%I(CFD,KD)
                          StopLD =BS%LStop%I(CFD,KD)
                          StrideD=StopLD-StartLD+1
                          MaxLD=BS%ASymm%I(2,CFD,KD)
                          StrideCD=StrideC*StrideD
                          IndexK=IndexK1
                          IndexA=IndexA1
                          DO 203 CFA=1,BS%NCFnc%I(KA)
                             StartLA=BS%LStrt%I(CFA,KA)
                             StopLA =BS%LStop%I(CFA,KA)
                             StrideA=StopLA-StartLA+1
                             MaxLA=BS%ASymm%I(2,CFA,KA)
                             IndexB=IndexB1
                             DO 204 CFB=1,BS%NCFnc%I(KB)
                                StartLB=BS%LStrt%I(CFB,KB)
                                StopLB =BS%LStop%I(CFB,KB)
                                StrideB=StopLB-StartLB+1
                                MaxLB=BS%ASymm%I(2,CFB,KB)
                                StrideAB=StrideA*StrideB
                                DO 300 PFA=1,BS%NPFnc%I(CFA,KA)
                                   DO 301 PFC=1,BS%NPFnc%I(CFC,KC)
                                      ZetaA=BS%Expnt%D(PFA,CFA,KA)
                                      ZetaC=BS%Expnt%D(PFC,CFC,KC)
                                      EtaAC=ZetaA+ZetaC
                                      ZAC  =ZetaA*ZetaC
                                      EtaIn=1.0D0/EtaAC
                                      XiAC =ZetaA*ZetaC*EtaIn
                                      ExpAC=DEXP(-XiAC*AC2)
                                      Px=(ZetaA*Ax+ZetaC*Cx)*EtaIn
                                      Py=(ZetaA*Ay+ZetaC*Cy)*EtaIn
                                      Pz=(ZetaA*Az+ZetaC*Cz)*EtaIn
                                      PAx=Px-Ax
                                      PAy=Py-Ay
                                      PAz=Pz-Az
                                      PCx=Px-Cx
                                      PCy=Py-Cy
                                      PCz=Pz-Cz
                                      ! Compute McMurchie Davidson E coefficients for HG Primitives
                                      CALL MD2TRR(BS%NASym,0,MaxLA,MaxLC,EtaAC,Eac,PAx,PCx,PAy,PCy,PAz,PCz)
                                      ! CALL MD2TRRx(NASym,0,MaxLA,MaxLC,EtaAC,PAx,PCx,PAy,PCy,PAz,PCz,EACx,EACy,EACz)
                                      DO 400 PFB=1,BS%NPFnc%I(CFB,KB)
                                         DO 401 PFD=1,BS%NPFnc%I(CFD,KD)
                                            ZetaB=BS%Expnt%D(PFB,CFB,KB)
                                            ZetaD=BS%Expnt%D(PFD,CFD,KD)
                                            RhoBD=ZetaB+ZetaD
                                            ZBD  =ZetaB*ZetaD
                                            RhoIn=1.0D0/RhoBD
                                            XiBD =ZetaB*ZetaD*RhoIn
                                            ExpBD=DEXP(-XiBD*BD2)
                                            Qx=(ZetaB*Bx+ZetaD*Dx)*RhoIn
                                            Qy=(ZetaB*By+ZetaD*Dy)*RhoIn
                                            Qz=(ZetaB*Bz+ZetaD*Dz)*RhoIn
                                            QBx=Qx-Bx
                                            QBy=Qy-By
                                            QBz=Qz-Bz
                                            QDx=Qx-Dx
                                            QDy=Qy-Dy
                                            QDz=Qz-Dz
                                            ! Compute McMurchie Davidson E coefficients for HG Primitives
                                            CALL MD2TRR(BS%NASym,0,MaxLB,MaxLD,RhoBD,Ebd,QBx,QDx,QBy,QDy,QBz,QDz)
                                            ! CALL MD2TRRx(NASym,0,QBx,QDx,QBy,QDy,QBz,QDz,EBDx,EBDy,EBDz)
                                            PQx=Px-Qx
                                            PQy=Py-Qy
                                            PQz=Pz-Qz
                                            PQ2=PQx*PQx+PQy*PQy+PQz*PQz
                                            W=EtaAC*RhoBD/(EtaAC+RhoBD)
                                            U=34.986836655249725693D0/(EtaAC*RhoBD*DSQRT(EtaAC+RhoBD))
                                            LTot=MaxLA+MaxLB+MaxLC+MaxLD
                                            CALL AuxInts(Max4L,LTot,AuxR,W,W*PQ2)
                                            CALL MD3TRR(Max4L,LTot,MDR,AuxR,U,PQx,PQy,PQz)
                                            I=0
                                            IA=IndexA
                                            DO 500 LMNA=StartLA,StopLA
                                               IA=IA+1
                                               LA=BS%LxDex%I(LMNA)
                                               MA=BS%LyDex%I(LMNA)
                                               NA=BS%LzDex%I(LMNA)
                                               CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                                               IB=IndexB
                                               DO 501 LMNB=StartLB,StopLB
                                                  IB=IB+1
                                                  LB=BS%LxDex%I(LMNB)
                                                  MB=BS%LyDex%I(LMNB)
                                                  NB=BS%LzDex%I(LMNB)
                                                  CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                                                  IC=IndexC
                                                  DO 502 LMNC=StartLC,StopLC
                                                     IC=IC+1
                                                     LC=BS%LxDex%I(LMNC)
                                                     MC=BS%LyDex%I(LMNC)
                                                     NC=BS%LzDex%I(LMNC)
                                                     CC=BS%CCoef%D(LMNC,PFC,CFC,KC)
                                                     ID=IndexD
                                                     DO 503 LMND=StartLD,StopLD
                                                        ID=ID+1
                                                        LD=BS%LxDex%I(LMND)
                                                        MD=BS%LyDex%I(LMND)
                                                        ND=BS%LzDex%I(LMND)
                                                        CD=BS%CCoef%D(LMND,PFD,CFD,KD)
                                                        CCoAC=CA*CC*ExpAC
                                                        CCoBD=CB*CD*ExpBD
                                                        I=I+1
                                                        DO LAC=0,LA+LC
                                                           DO MAC=0,MA+MC
                                                              DO NAC=0,NA+NC
                                                                 DO LBD=0,LB+LD
                                                                    DO MBD=0,MB+MD
                                                                       DO NBD=0,NB+ND
                                                                          Integrals%D(AtA,AtC,AtB,AtD)=Integrals%D(AtA,AtC,AtB,AtD)+CCoAC*CCoBD &
!                                                                          Integrals%D(IA,IB,IC,ID)=Integrals%D(IA,IB,IC,ID)+CCoAC*CCoBD &
                                                                               *Eac(1,LA,LC,LAC)       &
                                                                               *Eac(2,MA,MC,MAC)       &
                                                                               *Eac(3,NA,NC,NAC)       &
                                                                               *Ebd(1,LB,LD,LBD)       &
                                                                               *Ebd(2,MB,MD,MBD)       &
                                                                               *Ebd(3,NB,ND,NBD)       &
                                                                               *(-1.D0)**(LBD+MBD+NBD)   &
                                                                               *MDR(LAC+LBD,MAC+MBD,NAC+NBD,0)
                                                                       ENDDO
                                                                    ENDDO
                                                                 ENDDO
                                                              ENDDO
                                                           ENDDO
                                                        ENDDO
503                                                  CONTINUE
502                                               CONTINUE
501                                            CONTINUE
500                                         CONTINUE
401                                      CONTINUE
400                                   CONTINUE
301                                CONTINUE
300                             CONTINUE
                                IndexK=IndexK+StrideAB
                                IndexB=IndexB+StrideB
204                          CONTINUE
                             IndexA=IndexA+StrideA
203                       CONTINUE
                          IndexD=IndexD+StrideD
202                    CONTINUE
                       IndexC=IndexC+StrideC
201                 CONTINUE
1040       CONTINUE
                    IndexK1=IndexK1+BS%BFKnd%I(KA)*BS%BFKnd%I(KB)
                    IndexB1=IndexB1+BS%BFKnd%I(KB)
104              CONTINUE
1030       CONTINUE
                 IndexA1=IndexA1+BS%BFKnd%I(KA)
103           CONTINUE
1020       CONTINUE
              IndexD1=IndexD1+BS%BFKnd%I(KD)
102        CONTINUE
1010       CONTINUE
           IndexC1=IndexC1+BS%BFKnd%I(KC)
101     CONTINUE
!
END SUBROUTINE IntegralsC



END MODULE
