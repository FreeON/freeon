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


!!$
!         DO LAB=0,LTot
!            DO LB=MD0,MaxLB
!               DO LA=MD0,MaxLA
!                  WRITE(*,22)LA,LB,LAB,MD(:,LA,LB,LAB)
!               22 FORMAT(3(I3,","),3(D14.6,","))
!               ENDDO
!            ENDDO
!         ENDDO                                                          
      END SUBROUTINE MD2TRR
!
!!%======================================================================================================
!
!    NOTE THAT D[e,a] is NOT the same as f!!  SEE FOR EXAMPLE Helgaker and Taylor, TCA 83, p177 (1992), 
!    Eq (14) VS. Eqs. (19) and (20).  
!
!!$  e[0,0,0]:=Exp[-chi(a-b)^2];
!!$  e[i_,j_,n_]:=0/;(TrueQ[n>i+j]||TrueQ[n<0]);
!!$  e[i_,j_,n_]:=(e[i-1,j,n-1]/(2z)+PA*e[i-1,j,n]+(n+1)*e[i-1,j,n+1])/;TrueQ[j==0];
!!$  e[i_,j_,n_]:=(e[i,j-1,n-1]/(2z)+PB*e[i,j-1,n]+(n+1)*e[i,j-1,n+1])/;TrueQ[j!=0];
!!$
!!$  de[0,0,0]:=-2chi(a-b)Exp[-chi(a-b)^2];
!!$  de[i_,j_,n_]:=0/;(TrueQ[n>i+j]||TrueQ[n<0]);
!!$  de[i_,j_,n_]:=(de[i-1,j,n-1]/(2z)+(PA*de[i-1,j,n]-(zb/z)e[i-1,j,n])+(n+1)*de[i-1,j,n+1])/;TrueQ[j==0];
!!$  de[i_,j_,n_]:=(de[i,j-1,n-1]/(2z)+(PB*de[i,j-1,n]+(za/z)e[i,j-1,n])+(n+1)*de[i,j-1,n+1])/;TrueQ[j!=0];
!!$
!!$  f[i_, j_, n_] := 2 za e[i + 1, j, n] - i*e[i - 1, j, n];
!!$  
!!%======================================================================================================
      SUBROUTINE dMD2TRR(NASym,MD0,MaxLA,MaxLB,Za,Zb,EtaAB,MD,dMD, &
                         PAx,PBx,PAy,PBy,PAz,PBz) 
         REAL(DOUBLE), INTENT(IN)  :: EtaAB,PAx,PBx,PAy,PBy,PAz,PBz
         REAL(DOUBLE)              :: RL1,TwoZ,Za,Zb,ChiAB,ZaZ,ZbZ
         INTEGER,      INTENT(IN)  :: NASym,MD0,MaxLA,MaxLB
         REAL(DOUBLE), INTENT(OUT) :: MD(3,MD0:NASym,MD0:NASym,MD0:2*NASym)
         REAL(DOUBLE), INTENT(OUT) ::dMD(3,MD0:NASym,MD0:NASym,MD0:2*NASym)
         INTEGER                   :: LTot,LA,LB,LAB
         !
         ChiAB=Za*Zb/EtaAB        
         !
         LTot=MaxLA+MaxLB
         DO LAB=MD0,LTot
            DO LB=MD0,MaxLB
               DO LA=MD0,MaxLA
                  MD(1,LA,LB,LAB)=Zero
                  MD(2,LA,LB,LAB)=Zero
                  MD(3,LA,LB,LAB)=Zero
                  dMD(1,LA,LB,LAB)=Zero
                  dMD(2,LA,LB,LAB)=Zero
                  dMD(3,LA,LB,LAB)=Zero
               ENDDO
            ENDDO
         ENDDO          

         !
         MD(1,0,0,0)=One
         MD(2,0,0,0)=One
         MD(3,0,0,0)=One
         !
         dMD(1,0,0,0)=-Two*ChiAB*(PBx-PAx)
         dMD(2,0,0,0)=-Two*ChiAB*(PBy-PAy)
         dMD(3,0,0,0)=-Two*ChiAB*(PBz-PAz)
         !
         IF(LTot.EQ.0)RETURN
         !
         TwoZ=Half/EtaAB
         ZaZ=Za/EtaAB
         ZbZ=Zb/EtaAB
         !
         DO LA=1,MaxLA
            !
            MD(1,LA,0,0)=PAx*MD(1,LA-1,0,0)+MD(1,LA-1,0,1)
            MD(2,LA,0,0)=PAy*MD(2,LA-1,0,0)+MD(2,LA-1,0,1)
            MD(3,LA,0,0)=PAz*MD(3,LA-1,0,0)+MD(3,LA-1,0,1)
            !
            dMD(1,LA,0,0)=PAx*dMD(1,LA-1,0,0)-ZbZ*MD(1,LA-1,0,0)+dMD(1,LA-1,0,1)
            dMD(2,LA,0,0)=PAy*dMD(2,LA-1,0,0)-ZbZ*MD(2,LA-1,0,0)+dMD(2,LA-1,0,1)
            dMD(3,LA,0,0)=PAz*dMD(3,LA-1,0,0)-ZbZ*MD(3,LA-1,0,0)+dMD(3,LA-1,0,1)
            !
            DO LAB=1,LA-1           
               !
               RL1=DBLE(LAB+1)                                      
               !
               MD(1,LA,0,LAB)=TwoZ*MD(1,LA-1,0,LAB-1) &
                             + PAx*MD(1,LA-1,0,LAB  ) &           
                             + RL1*MD(1,LA-1,0,LAB+1)                         
               MD(2,LA,0,LAB)=TwoZ*MD(2,LA-1,0,LAB-1) &
                             + PAy*MD(2,LA-1,0,LAB  ) &            
                             + RL1*MD(2,LA-1,0,LAB+1)                         
               MD(3,LA,0,LAB)=TwoZ*MD(3,LA-1,0,LAB-1) &
                             + PAz*MD(3,LA-1,0,LAB  ) &            
                             + RL1*MD(3,LA-1,0,LAB+1)                         
               !
               dMD(1,LA,0,LAB)=TwoZ*dMD(1,LA-1,0,LAB-1) &
                             +  PAx*dMD(1,LA-1,0,LAB  ) &    
                             -  ZbZ* MD(1,LA-1,0,LAB  ) &           
                             +  RL1*dMD(1,LA-1,0,LAB+1)                         

               dMD(2,LA,0,LAB)=TwoZ*dMD(2,LA-1,0,LAB-1) &
                             +  PAy*dMD(2,LA-1,0,LAB  ) &            
                             -  ZbZ* MD(2,LA-1,0,LAB  ) &            
                             +  RL1*dMD(2,LA-1,0,LAB+1)                         

               dMD(3,LA,0,LAB)=TwoZ*dMD(3,LA-1,0,LAB-1) &
                             +  PAz*dMD(3,LA-1,0,LAB  ) &            
                             -  ZbZ* MD(3,LA-1,0,LAB  ) &            
                             +  RL1*dMD(3,LA-1,0,LAB+1)                         
               !
            ENDDO                                                
            !
            MD(1,LA,0,LA)=TwoZ*MD(1,LA-1,0,LA-1)+PAx*MD(1,LA-1,0,LA)               
            MD(2,LA,0,LA)=TwoZ*MD(2,LA-1,0,LA-1)+PAy*MD(2,LA-1,0,LA)               
            MD(3,LA,0,LA)=TwoZ*MD(3,LA-1,0,LA-1)+PAz*MD(3,LA-1,0,LA)               
            !
            dMD(1,LA,0,LA)=TwoZ*dMD(1,LA-1,0,LA-1)+PAx*dMD(1,LA-1,0,LA)-ZbZ*MD(1,LA-1,0,LA)               
            dMD(2,LA,0,LA)=TwoZ*dMD(2,LA-1,0,LA-1)+PAy*dMD(2,LA-1,0,LA)-ZbZ*MD(2,LA-1,0,LA)               
            dMD(3,LA,0,LA)=TwoZ*dMD(3,LA-1,0,LA-1)+PAz*dMD(3,LA-1,0,LA)-ZbZ*MD(3,LA-1,0,LA)               
            !
         ENDDO 
         !
         DO LB=1,MaxLB
            DO LA=0,MaxLA  
               !
               MD(1,LA,LB,0)=PBx*MD(1,LA,LB-1,0)+MD(1,LA,LB-1,1)
               MD(2,LA,LB,0)=PBy*MD(2,LA,LB-1,0)+MD(2,LA,LB-1,1)
               MD(3,LA,LB,0)=PBz*MD(3,LA,LB-1,0)+MD(3,LA,LB-1,1)
               !
               dMD(1,LA,LB,0)=PBx*dMD(1,LA,LB-1,0)+ZaZ*MD(1,LA,LB-1,0)+dMD(1,LA,LB-1,1)
               dMD(2,LA,LB,0)=PBy*dMD(2,LA,LB-1,0)+ZaZ*MD(2,LA,LB-1,0)+dMD(2,LA,LB-1,1)
               dMD(3,LA,LB,0)=PBz*dMD(3,LA,LB-1,0)+ZaZ*MD(3,LA,LB-1,0)+dMD(3,LA,LB-1,1)
               !
               DO LAB=1,LTot-1
                  !
                  RL1=DBLE(LAB+1)
                  !
                  MD(1,LA,LB,LAB)=TwoZ*MD(1,LA,LB-1,LAB-1) &
                                 + PBx*MD(1,LA,LB-1,LAB  ) &        
                                 + RL1*MD(1,LA,LB-1,LAB+1)
                  MD(2,LA,LB,LAB)=TwoZ*MD(2,LA,LB-1,LAB-1) &
                                 + PBy*MD(2,LA,LB-1,LAB  ) &         
                                 + RL1*MD(2,LA,LB-1,LAB+1)
                  MD(3,LA,LB,LAB)=TwoZ*MD(3,LA,LB-1,LAB-1) & 
                                 + PBz*MD(3,LA,LB-1,LAB  ) &        
                                 + RL1*MD(3,LA,LB-1,LAB+1)
                  !
                  dMD(1,LA,LB,LAB)=TwoZ*dMD(1,LA,LB-1,LAB-1) &
                                  + PBx*dMD(1,LA,LB-1,LAB  ) &        
                                  + ZaZ* MD(1,LA,LB-1,LAB  ) &        
                                  + RL1*dMD(1,LA,LB-1,LAB+1)
                  dMD(2,LA,LB,LAB)=TwoZ*dMD(2,LA,LB-1,LAB-1) &
                                  + PBy*dMD(2,LA,LB-1,LAB  ) &         
                                  + ZaZ* MD(2,LA,LB-1,LAB  ) &         
                                  + RL1*dMD(2,LA,LB-1,LAB+1)
                  dMD(3,LA,LB,LAB)=TwoZ*dMD(3,LA,LB-1,LAB-1) & 
                                  + PBz*dMD(3,LA,LB-1,LAB  ) &        
                                  + ZaZ* MD(3,LA,LB-1,LAB  ) &        
                                  + RL1*dMD(3,LA,LB-1,LAB+1)
                  !
               ENDDO
               !
               MD(1,LA,LB,LTot)=TwoZ*MD(1,LA,LB-1,LAB-1)+PBx*MD(1,LA,LB-1,LAB)
               MD(2,LA,LB,LTot)=TwoZ*MD(2,LA,LB-1,LAB-1)+PBy*MD(2,LA,LB-1,LAB)
               MD(3,LA,LB,LTot)=TwoZ*MD(3,LA,LB-1,LAB-1)+PBz*MD(3,LA,LB-1,LAB)
               !
               dMD(1,LA,LB,LTot)=TwoZ*dMD(1,LA,LB-1,LAB-1)+PBx*dMD(1,LA,LB-1,LAB)+ZaZ*MD(1,LA,LB-1,LAB)
               dMD(2,LA,LB,LTot)=TwoZ*dMD(2,LA,LB-1,LAB-1)+PBy*dMD(2,LA,LB-1,LAB)+ZaZ*MD(2,LA,LB-1,LAB)
               dMD(3,LA,LB,LTot)=TwoZ*dMD(3,LA,LB-1,LAB-1)+PBz*dMD(3,LA,LB-1,LAB)+ZaZ*MD(3,LA,LB-1,LAB)
               !
            ENDDO
          ENDDO
!!$
!!$         DO LAB=0,LTot
!!$            DO LB=MD0,MaxLB
!!$               DO LA=MD0,MaxLA
!!$                  WRITE(*,22)LA,LB,LAB,MD(1,LA,LB,LAB)
!!$               22 FORMAT(3(I3,","),3(D14.6,","))
!!$               ENDDO
!!$            ENDDO
!!$         ENDDO                                                          

      END SUBROUTINE dMD2TRR


      SUBROUTINE MD2TRR2(NASym,MD0,MaxLA,MaxLB,EtaAB,MD, &
                        PAx,PBx,PAy,PBy,PAz,PBz) 
         REAL(DOUBLE), INTENT(IN)  :: EtaAB,PAx,PBx,PAy,PBy,PAz,PBz
         REAL(DOUBLE)              :: RL1,TwoZ
         INTEGER,      INTENT(IN)  :: NASym,MD0,MaxLA,MaxLB
         REAL(DOUBLE), INTENT(OUT) :: MD(3,MD0:NASym,MD0:NASym,MD0:2*NASym)
         INTEGER                   :: LTot,LA,LB,LAB
         LTot=MaxLA+MaxLB

!         WRITE(*,*)' PAX = ',PAX
!         WRITE(*,*)' PBX = ',PBX
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


!          WRITE(*,*)' Ex(1,0,1) = ',MD(1,1,0,1)
!          WRITE(*,*)' Ex(0,1,1) = ',MD(1,0,1,1)
!          WRITE(*,*)' Ex(1,0,0) = ',MD(1,1,0,0)
!          WRITE(*,*)' Ex(0,1,0) = ',MD(1,0,1,0)

!         DO LAB=0,LTot
!            DO LB=MD0,MaxLB
!               DO LA=MD0,MaxLA
!                  WRITE(*,22)LA,LB,LAB,MD(:,LA,LB,LAB)
!               22 FORMAT(3(I3,","),3(D14.6,","))
!               ENDDO
!            ENDDO
!         ENDDO                                                          
      END SUBROUTINE MD2TRR2
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
         INTEGER,PARAMETER                          :: LPlus=300
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

      SUBROUTINE OvrInts(MaxL,LTot,AuxR,Omega,T)
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
            OmegaJ=One
            TwoO=-Two*Omega
            ET=EXP(-T)
            DO J=0,LTot
               AuxR(J)=OmegaJ*ET
               OmegaJ=TwoO*OmegaJ
            ENDDO
            RETURN            
          END SUBROUTINE OvrInts




!-------------------------------------------------------------
!     Compute the auxiliary integrals R_{0,0,0,J}(T)
!
      SUBROUTINE ErrInts(MaxL,LTot,ErrR,Omega,T)
         REAL(DOUBLE),                  INTENT(IN)  :: Omega,T
         INTEGER,                       INTENT(IN)  :: MaxL,LTot 
         REAL(DOUBLE),DIMENSION(0:MaxL),INTENT(OUT) :: ErrR
         REAL(DOUBLE),PARAMETER                     :: Switch=35.0D0
         INTEGER,PARAMETER                          :: LPlus=250
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


END MODULE
