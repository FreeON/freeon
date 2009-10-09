!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
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

SUBROUTINE Integrals2E(BS,GM,TwoE)

  TYPE(CRDS)     :: GM
  TYPE(BSET)     :: BS
  !        TYPE(CellSet)  :: CS
  INTEGER :: MaxL, Max2L, Max4L, MaxLMN, Max4LMN
  REAL(DOUBLE), DIMENSION(0:4*BS%NASym)                                        :: AuxR
  REAL(DOUBLE), DIMENSION(1:3,0:BS%NASym,0:BS%NASym,0:2*BS%NASym)              :: Eab,Ecd
  REAL(DOUBLE), DIMENSION(0:4*BS%NASym,0:4*BS%NASym,0:4*BS%NASym,0:4*BS%NASym) :: MDR
  REAL(DOUBLE), DIMENSION(BS%NBasF,BS%NBasF,BS%NBasF,BS%NBasF)                 :: TwoE

!  INTEGER :: NC1,NC2,NC3,NC4
  !
  INTEGER :: AtA,CFA,PFA,StartLA,StopLA,StrideA
  INTEGER :: AtB,CFB,PFB,StartLB,StopLB,StrideB
  INTEGER :: AtC,CFC,PFC,StartLC,StopLC,StrideC
  INTEGER :: AtD,CFD,PFD,StartLD,StopLD,StrideD
  INTEGER :: StrideAB,StrideCD

  INTEGER :: IndexA1,IndexA,IndexB1,IndexB
  INTEGER :: IndexC1,IndexC,IndexD1,IndexD
  INTEGER :: KC,KA,KB,KD
  INTEGER :: MaxLA,MaxLB,MaxLC,MaxLD,LTot
  INTEGER :: LA,MA,NA,LB,MB,NB,LC,NC,MC,LD,MD,ND
  INTEGER :: I,IK,IA,IB,IC,ID
  INTEGER :: LMNA,LMNB,LMNC,LMND,LAB,MAB,NAB,LCD,MCD,NCD
  REAL(DOUBLE) :: Ax,Ay,Az,Bx,By,Bz
  REAL(DOUBLE) :: Cx,Cy,Cz,Dx,Dy,Dz
  REAL(DOUBLE) :: Px,Py,Pz,Qx,Qy,Qz
  REAL(DOUBLE) :: CDx,CDy,CDz,CD2,ABx,ABy,ABz,AB2
  REAL(DOUBLE) :: PAx,PAy,PAz,PBx,PBy,PBz
  REAL(DOUBLE) :: QCx,QCy,QCz,QDx,QDy,QDz
  REAL(DOUBLE) :: PQx,PQy,PQz,PQ2
  REAL(DOUBLE) :: ZetaA,ZetaB,ZetaC,ZetaD,EtaAB,RhoCD
  REAL(DOUBLE) :: Zab,Zcd,EtaIn,RhoIn,XiCD,XiAB,ExpAB,ExpCD,W,U
  REAL(DOUBLE) :: CA,CB,CC,CD,CCoAB,CCoCD
  !
  MaxL    =   BS%NASym
  Max2L   = 2*BS%NASym
  Max4L   = 4*BS%NASym
  MaxLMN  = BS%LMNLen
  Max4LMN = BS%LMNLen**4
  !
  TwoE(:,:,:,:)=Zero

  IndexA1=0
  DO AtA=1,NAtoms
     Ax=GM%Carts%D(1,AtA)
     Ay=GM%Carts%D(2,AtA)
     Az=GM%Carts%D(3,AtA)
     KA=GM%AtTyp%I(AtA)

     IndexB1=0
     DO AtB = 1,NAtoms
        Bx=GM%Carts%D(1,AtB)
        By=GM%Carts%D(2,AtB)
        Bz=GM%Carts%D(3,AtB)
        KB=GM%AtTyp%I(AtB)

        IndexC1 = 0
        DO AtC=1,NAtoms
           Cx=GM%Carts%D(1,AtC)
           Cy=GM%Carts%D(2,AtC)
           Cz=GM%Carts%D(3,AtC)
           KC=GM%AtTyp%I(AtC)

           IndexD1=0
           DO AtD=1,NAtoms
              Dx=GM%Carts%D(1,AtD)
              Dy=GM%Carts%D(2,AtD)
              Dz=GM%Carts%D(3,AtD)
              KD=GM%AtTyp%I(AtD)

!!$              ACx=Ax-Cx
!!$              ACy=Ay-Cy
!!$              ACz=Az-Cz
!!$              AC2=ACx*ACx+ACy*ACy+ACz*ACz
              ABx=Ax-Bx
              ABy=Ay-By
              ABz=Az-Bz
              AB2=ABx*ABx+ABy*ABy+ABz*ABz
!!$              BDx=Bx-Dx
!!$              BDy=By-Dy
!!$              BDz=Bz-Dz
!!$              BD2=BDx*BDx+BDy*BDy+BDz*BDz

              CDx=Cx-Dx
              CDy=Cy-Dy
              CDz=Cz-Dz
              CD2=CDx*CDx+CDy*CDy+CDz*CDz

              IndexA=IndexA1
              DO CFA=1,BS%NCFnc%I(KA)
                 StartLA=BS%LStrt%I(CFA,KA)
                 StopLA =BS%LStop%I(CFA,KA)
                 StrideA=StopLA-StartLA+1
                 MaxLA=BS%ASymm%I(2,CFA,KA)

                 IndexB=IndexB1
                 DO CFB=1,BS%NCFnc%I(KB)
                    StartLB=BS%LStrt%I(CFB,KB)
                    StopLB =BS%LStop%I(CFB,KB)
                    StrideB=StopLB-StartLB+1
                    MaxLB=BS%ASymm%I(2,CFB,KB)

                    IndexC=IndexC1
                    DO CFC=1,BS%NCFnc%I(KC)
                       StartLC=BS%LStrt%I(CFC,KC)
                       StopLC =BS%LStop%I(CFC,KC)
                       StrideC=StopLC-StartLC+1
                       MaxLC=BS%ASymm%I(2,CFC,KC)

                       IndexD=IndexD1
                       DO CFD=1,BS%NCFnc%I(KD)
                          StartLD=BS%LStrt%I(CFD,KD)
                          StopLD =BS%LStop%I(CFD,KD)
                          StrideD=StopLD-StartLD+1
                          MaxLD=BS%ASymm%I(2,CFD,KD)

                          DO PFA=1,BS%NPFnc%I(CFA,KA)
                             DO PFB=1,BS%NPFnc%I(CFB,KB)
                                ZetaA=BS%Expnt%D(PFA,CFA,KA)
                                ZetaB=BS%Expnt%D(PFB,CFB,KB)
                                EtaAB=ZetaA+ZetaB
                                ZAB  =ZetaA*ZetaB
                                EtaIn=1.0D0/EtaAB
                                XiAB =ZetaA*ZetaB*EtaIn
                                ExpAB=DEXP(-XiAB*AB2)
                                Px=(ZetaA*Ax+ZetaB*Bx)*EtaIn
                                Py=(ZetaA*Ay+ZetaB*By)*EtaIn
                                Pz=(ZetaA*Az+ZetaB*Bz)*EtaIn
                                PAx=Px-Ax
                                PAy=Py-Ay
                                PAz=Pz-Az
                                PBx=Px-Bx
                                PBy=Py-By
                                PBz=Pz-Bz
                                CALL MD2TRR(BS%NASym,0,MaxLA,MaxLB,EtaAB,Eab,PAx,PBx,PAy,PBy,PAz,PBz)
                                DO PFC=1,BS%NPFnc%I(CFC,KC)
                                   DO PFD=1,BS%NPFnc%I(CFD,KD)
                                      ZetaC=BS%Expnt%D(PFC,CFC,KC)
                                      ZetaD=BS%Expnt%D(PFD,CFD,KD)
                                      RhoCD=ZetaC+ZetaD
                                      ZCD  =ZetaC*ZetaD
                                      RhoIn=1.0D0/RhoCD
                                      XiCD =ZetaC*ZetaD*RhoIn
                                      ExpCD=DEXP(-XiCD*CD2)
                                      Qx=(ZetaC*Cx+ZetaD*Dx)*RhoIn
                                      Qy=(ZetaC*Cy+ZetaD*Dy)*RhoIn
                                      Qz=(ZetaC*Cz+ZetaD*Dz)*RhoIn
                                      QCx=Qx-Cx
                                      QCy=Qy-Cy
                                      QCz=Qz-Cz
                                      QDx=Qx-Dx
                                      QDy=Qy-Dy
                                      QDz=Qz-Dz
                                      CALL MD2TRR(BS%NASym,0,MaxLC,MaxLD,RhoCD,Ecd,QCx,QDx,QCy,QDy,QCz,QDz)
                                      !
                                      PQx=Px-Qx
                                      PQy=Py-Qy
                                      PQz=Pz-Qz
                                      PQ2=PQx*PQx+PQy*PQy+PQz*PQz
                                      W=EtaAB*RhoCD/(EtaAB+RhoCD)
                                      U=34.986836655249725693D0/(EtaAB*RhoCD*DSQRT(EtaAB+RhoCD))
                                      LTot=MaxLA+MaxLB+MaxLC+MaxLD
                                      CALL AuxInts(Max4L,LTot,AuxR,W,W*PQ2)
                                      CALL MD3TRR(Max4L,LTot,MDR,AuxR,U,PQx,PQy,PQz)

                                      IA=IndexA
                                      DO LMNA=StartLA,StopLA
                                         IA=IA+1
                                         LA=BS%LxDex%I(LMNA)
                                         MA=BS%LyDex%I(LMNA)
                                         NA=BS%LzDex%I(LMNA)
                                         CA=BS%CCoef%D(LMNA,PFA,CFA,KA)

                                         IB=IndexB
                                         DO LMNB=StartLB,StopLB
                                            IB=IB+1
                                            LB=BS%LxDex%I(LMNB)
                                            MB=BS%LyDex%I(LMNB)
                                            NB=BS%LzDex%I(LMNB)
                                            CB=BS%CCoef%D(LMNB,PFB,CFB,KB)

                                            IC=IndexC
                                            DO LMNC=StartLC,StopLC
                                               IC=IC+1
                                               LC=BS%LxDex%I(LMNC)
                                               MC=BS%LyDex%I(LMNC)
                                               NC=BS%LzDex%I(LMNC)
                                               CC=BS%CCoef%D(LMNC,PFC,CFC,KC)

                                               ID=IndexD
                                               DO LMND=StartLD,StopLD
                                                  ID=ID+1
                                                  LD=BS%LxDex%I(LMND)
                                                  MD=BS%LyDex%I(LMND)
                                                  ND=BS%LzDex%I(LMND)
                                                  CD=BS%CCoef%D(LMND,PFD,CFD,KD)

                                                  CCoAB=CA*CB*ExpAB
                                                  CCoCD=CC*CD*ExpCD

                                                  DO LAB=0,LA+LB
                                                     DO MAB=0,MA+MB
                                                        DO NAB=0,NA+NB

                                                           DO LCD=0,LC+LD
                                                              DO MCD=0,MC+MD
                                                                 DO NCD=0,NC+ND

                                                                   TwoE(IA,IB,IC,ID)=TwoE(IA,IB,IC,ID)+CCoAB*CCoCD &
                                                                         *Eab(1,LA,LB,LAB)       &
                                                                         *Eab(2,MA,MB,MAB)       &
                                                                         *Eab(3,NA,NB,NAB)       &
                                                                         *Ecd(1,LC,LD,LCD)       &
                                                                         *Ecd(2,MC,MD,MCD)       &
                                                                         *Ecd(3,NC,ND,NCD)       &
                                                                         *(-1.D0)**(LCD+MCD+NCD)   &
                                                                         *MDR(LAB+LCD,MAB+MCD,NAB+NCD,0)
                                                                 ENDDO
                                                              ENDDO
                                                           ENDDO
                                                        ENDDO
                                                     ENDDO
                                                  ENDDO

                                               ENDDO
                                            ENDDO
                                         ENDDO
                                      ENDDO
                                   ENDDO
                                ENDDO
                             ENDDO
                          ENDDO
                          IndexD=IndexD+StrideD
                       ENDDO
                       IndexC=IndexC+StrideC
                    ENDDO
                    IndexB=IndexB+StrideB
                 ENDDO
                 IndexA=IndexA+StrideA
              ENDDO
              IndexD1=IndexD1+BS%BFKnd%I(KD)
           ENDDO
           IndexC1=IndexC1+BS%BFKnd%I(KC)
        ENDDO
        IndexB1=IndexB1+BS%BFKnd%I(KB)
     ENDDO
     IndexA1=IndexA1+BS%BFKnd%I(KA)
  ENDDO
  !
END SUBROUTINE Integrals2E




END MODULE
