!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 1999, The University of California
!
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
         REAL(DOUBLE), INTENT(OUT) :: MD(3,MD0:NASym,MD0:NASym,0:2*NASym)
         INTEGER                   :: LTot,LA,LB,LAB
         LTot=MaxLA+MaxLB
         DO LAB=0,LTot
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
!---------------------------------------------------
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

END MODULE
