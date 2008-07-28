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
MODULE RayleighQuotientIteration
  USE InOut
  USE Macros
  USE MemMan
  USE Indexing
  USE AtomPairs
  USE PrettyPrint
  USE GlobalScalars
  USE ControlStructures
  USE OptionKeys
  USE McMurchie
  USE MondoLogger
  IMPLICIT NONE
CONTAINS

  SUBROUTINE RQI(N,M,Nam,O,S,G,B)
    INTEGER            :: N,M,I,J,K,L
    TYPE(FileNames)    :: Nam
    TYPE(State)        :: S
    TYPE(Options)      :: O
    TYPE(CRDS)         :: G
    TYPE(BSET)         :: B
    !
    TYPE(BCSR)                      :: sP,sQ,sF,sZ
    TYPE(DBL_RNK2)                  :: P,Q,F,Z
    REAL(DOUBLE)                    :: Ek,EkOld,Beta,Lambda
    REAL(DOUBLE),DIMENSION(N,N)     :: Xk,Gk,Pk,LXk,LPk
    REAL(DOUBLE),DIMENSION(N,N)     :: XkOld,PkOld,GkOld
    REAL(DOUBLE),DIMENSION(N)       :: Values
    REAL(DOUBLE),DIMENSION(N,N,M)   :: Vectors
    REAL(DOUBLE),DIMENSION(N,N,N,N) :: TwoE


    CALL Get(sP,TrixFile("OrthoD",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME, &
         Stats_O=S%Current%I,OffSet_O=0))
    CALL Get(sF,TrixFile("OrthoF",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME, &
         Stats_O=S%Current%I,OffSet_O=0))
    CALL Get(sZ,TrixFile("X",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I))
    !
    CALL SetEq(P,sP)
    CALL SetEq(F,sF)
    CALL SetEq(Z,sZ)
    !
    CALL SetEq(Q,sP)
    DO I=1,N 
       Q%D(I,I)=One-P%D(I,I)
    ENDDO
    !
    CALL Integrals2E(B,G,TwoE)
    !
!!$    WRITE(*,*)TrixFile("OrthoD",PWD_O=N%M_SCRATCH,Name_O=Nam%SCF_NAME, &
!!$         Stats_O=S%Current%I,OffSet_O=0)
!!$    WRITE(*,*)TrixFile("OrthoF",PWD_O=N%M_SCRATCH,Name_O=Nam%SCF_NAME, &
!!$         Stats_O=S%Current%I,OffSet_O=0)
!!$    WRITE(*,*)TrixFile("X",PWD_O=N%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I)
!!$    DO I=1,N
!!$       DO J=1,N
!!$          DO K=1,N
!!$             DO L=1,N
!!$                CALL iPrint(TwoE(I,J,K,L),I,J,K,L,1,6)
!!$             ENDDO
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
    !
    DO I=1,1 !M
       !
       CALL RPAGuess(I,Xk)
       CALL Anihilate(N,P%D,Q%D,Xk)
       CALL LOn2(N,F%D,P%D,Z%D,TwoE,Xk,LXk)
       CALL Anihilate(N,P%D,Q%D,LXk)
       !
       Beta=Zero
       Ek=ThoulessQ(N,P%D,Xk,LXk) 
       !
       DO K=0,10
          !
          Gk=Two*(LXk-Ek*Xk)
          IF(K>0)THEN
             Beta=Pdot1(N,P%D,Gk-Gkold,Gk)/Pdot1(N,P%D,GkOld,GkOld)    			
          ENDIF
          Pk=Gk+Beta*PkOld  
          !
          CALL LOn2(NBasF,F%D,P%D,Z%D,TwoE,Pk,LPk)
          CALL RQILineSearch(N,P%D,Pk,Xk,LXk,LPk,Lambda)          
          !
          XkOld=Xk
          EkOld=Ek
          GkOld=Gk
          PkOld=Pk	   
          !
          Xk=XkOld+Lambda*Pk
          !
          CALL Anihilate(N,P%D,Q%D,Xk)
          CALL ReNorm(N,P%D,Xk)
          CALL LOn2(N,F%D,P%D,Z%D,TwoE,Xk,LXk)
          CALL Anihilate(N,P%D,Q%D,LXk)
          Ek=ThoulessQ(N,P%D,Xk,LXk) 

          WRITE(*,*)I,Ek
          !
       ENDDO
       !
       Values(I)=Ek 
       Vectors(:,:,I)=Xk
       !
    ENDDO
    !
  END SUBROUTINE RQI

     SUBROUTINE Anihilate(N,P,Q,X)
       INTEGER :: N
       REAL(DOUBLE),DIMENSION(N,N) :: P,Q,X
       X=ProjectPH(N,Q,P,X)
       !!       X=Project(N,Q,P,X)
     END SUBROUTINE Anihilate

     SUBROUTINE LOn2(N,F,P,Z,TwoE,X,LX)
       INTEGER :: N
       REAL(DOUBLE),DIMENSION(N,N)     :: F,P,Z,X,LX,AA,BB,Temp
       REAL(DOUBLE),DIMENSION(N,N,N,N) :: TwoE
       LX=LiouvAO(N,F  ,P  ,Z,TwoE,X )  
!         LiouvAO(N,For,Por,X,DSao,AA) 

     END SUBROUTINE LOn2

     FUNCTION ThoulessQ(N,P,X,LX) RESULT(Ek)
       !! So, where is  the denominator??
       INTEGER :: N
       REAL(DOUBLE) :: Ek
       REAL(DOUBLE),DIMENSION(N,N)     :: F,P,X,LX,Tmp1
       Ek=Pdot1(N,P,X,LX)
     END FUNCTION ThoulessQ


     FUNCTION LiouvDot(N,BB,DSao,temp2)  RESULT(temp1) 
       ! Calculates action of the Coulomb operator in AO space temp1=BB * (ij||kl) 
       IMPLICIT NONE
       INTEGER :: I,J,K,N,one
       REAL (DOUBLE),DIMENSION(N*N):: BB,temp2
       REAL (DOUBLE),DIMENSION(N,N):: temp1
       REAL(DOUBLE),DIMENSION(N*N,N*N)::	DSao                
       REAL(DOUBLE) :: ddot

       one=1 
       K=0    
       DO I=1,N
          DO J=1,N
             K=K+1
             temp2=DSao(:,K)
             !		temp1(J,I)= ddot(N*N,BB,one,temp2,one)     ! This line is 
             temp1(J,I)=DOT_PRODUCT(BB,Temp2)           ! the most CPU consuming step
          ENDDO
       END DO

     END FUNCTION LiouvDot

     !-------------------------------------------------------------------------------      
     FUNCTION LiouvAO(N,For,Por,X,DSao,AA)  RESULT(BB)
       ! Calculates action of the Liouville operator in AO space BB=L AA, (ij||kl) 
       IMPLICIT NONE
       INTEGER :: I,J,M,K,L,N,one
       REAL (DOUBLE),DIMENSION(N,N)::For,Por,AA,BB,X
       REAL(DOUBLE),DIMENSION(N,N,N,N)::	DSao                
       REAL(DOUBLE),DIMENSION(N,N) :: temp1,temp2
       REAL(DOUBLE) :: E,ddot

       ! AA to AO
       one=1 
       BB=MATMUL(TRANSPOSE(X),(MATMUL(AA,X)))      

       !        DO I=1,N
       !           DO J=1,N
       !              temp1(I,J)= 0.0
       !              DO K=1,N
       !                 DO L=1,N
       !                    temp1(I,J)=temp1(I,J)+BB(K,L)*DSao(K,L,I,J)
       !                 END DO
       !              END DO
       !!            temp2=DSao(:,:,I,J)       ! This and the next line is equivalent to Liouvdot 
       !!		temp1(I,J)= ddot(N*N,BB,one,temp2,one)
       !           END DO
       !        END DO   

       ! Extremely inefficient procedure above, trying to replace with DOT_PRODUCT or ddot
       temp1= LiouvDot(N,BB,DSao,temp2)  
       BB=MATMUL(For,AA)-MATMUL(AA,For)
       ! temp back to orthog
       temp2=MATMUL(TRANSPOSE(X),(MATMUL(temp1,X)))
       BB=BB+MATMUL(temp2,Por)-MATMUL(Por,temp2)

     END FUNCTION LiouvAO

     SUBROUTINE LiouvilleAO(N,For,Por,DSao,AA,BB,temp,X)  
       ! Calculates action of the Liouville operator in AO space BB=L AA, (ij||kl) 
       IMPLICIT NONE
       INTEGER :: I,J,M,K,L,N
       REAL (DOUBLE),DIMENSION(N,N)::For,Por,AA,BB,temp,X
       REAL(DOUBLE),DIMENSION(N,N,N,N)::	DSao                
       REAL(DOUBLE) :: E

       ! AA to AO 
       BB=MATMUL(TRANSPOSE(X),(MATMUL(AA,X)))      
       DO I=1,N
          DO J=1,N
             temp(I,J)=0.0  
             DO K=1,N
                DO L=1,N
                   temp(I,J)=temp(I,J)+BB(K,L)*DSao(K,L,I,J)
                END DO
             END DO
          END DO
       END DO
       BB=MATMUL(For,AA)-MATMUL(AA,For)
       ! temp back to orthog
       AA=MATMUL(TRANSPOSE(X),(MATMUL(temp,X)))
       BB=BB+MATMUL(AA,Por)-MATMUL(Por,AA)

     END SUBROUTINE LiouvilleAO

     SUBROUTINE RPAGuess(N,X)
       INTEGER :: N
       REAL(DOUBLE), DIMENSION(N,N) :: X
       X=One

!!$       do i=1,N
!!$          do j=1,N
!!$             Xk(i,j)= 2.0*rranf()-1.0  
!!$	  enddo
!!$       enddo
!!$       temp1 = ProjectPH(N,Qor,Por,Xk)       
!!$       Xk = temp1/sqrt(abs(Pdot1(N,Por,temp1,temp1,tmp1)))
     END SUBROUTINE RPAGuess

     SUBROUTINE ReNorm(N,P,X)
       INTEGER :: N
       REAL(DOUBLE) :: Norm
       REAL(DOUBLE),DIMENSION(N,N) :: P,X
       Norm=sqrt(abs(Pdot1(N,P,X,X)))
       X=X/Norm
     END SUBROUTINE ReNorm

     !************************************************************************      
     FUNCTION Pdot(N,P,AA,BB,CC) RESULT(Tr) 
       ! Calculates RPA scalar product Tr=Tr([AA^+,P],BB)

       IMPLICIT NONE
       INTEGER :: N
       REAL(DOUBLE) :: Tr
       REAL (DOUBLE),DIMENSION(N,N)::P,AA,BB,CC

       CC=MATMUL((MATMUL(TRANSPOSE(AA),P)-MATMUL(P,TRANSPOSE(AA))),BB)
       Tr=0.5*Trace2(CC,N)

     END FUNCTION Pdot

     !************************************************************************      
     FUNCTION Pdot1(N,P,AA,BB) RESULT(Tr) 
       ! Calculates RPA scalar product Tr=Tr([AA,P],BB^+)

       IMPLICIT NONE
       INTEGER :: N
       REAL(DOUBLE) :: Tr
       REAL (DOUBLE),DIMENSION(N,N)::P,AA,BB,CC

       CC=MATMUL((MATMUL(AA,P)-MATMUL(P,AA)),TRANSPOSE(BB))
       Tr=0.5*Trace2(CC,N)

     END FUNCTION Pdot1

     !-------------------------------------------------------------------------------      
     FUNCTION Project(N,P,Q,AA)  RESULT(BB)
       ! BB=P AA Q + Q AA P
       ! calculates  projection to p-h an h-p space using Q and P

       IMPLICIT NONE
       INTEGER :: N
       REAL (DOUBLE),DIMENSION(N,N)::P,Q,AA,BB

       BB=0.25*(MATMUL(MATMUL(P,AA),Q)+MATMUL(MATMUL(Q,AA),P))

     END FUNCTION Project

     !-------------------------------------------------------------------------------      
     FUNCTION ProjectPH(N,P,Q,AA)  RESULT(BB)
       ! BB=Q AA P (X-component, large)   (0  Y)
       ! BB=P AA Q (Y-component, small)   (X  0)
       ! calculates  projection to p-h OR h-p space using Q and P

       IMPLICIT NONE
       INTEGER :: N
       REAL (DOUBLE),DIMENSION(N,N)::P,Q,AA,BB
       BB=0.25*(MATMUL(MATMUL(P,AA),Q))

     END FUNCTION ProjectPH



   Subroutine RQILineSearch(N,P,Pk,Xk,LXk,LPk,Lambda)          
     INTEGER :: N
     REAL(DOUBLE) :: Lambda,Lambda_p,Lambda_m
     REAL (DOUBLE),DIMENSION(N,N)::P,Pk,Xk,LXk,LPk,Tmp1
     REAL(DOUBLE) :: XX,PP,XP,PX,PLP,XLP,XLX,PLX,AA,BB,CC

     XX =Pdot1(N,P,Xk,Xk) ! 1 - normalized by definition
     PP =Pdot1(N,P,Pk,Pk)
     XP =Pdot1(N,P,Xk,Pk)
     PX =Pdot1(N,P,Pk,Xk)
     PLP=Pdot1(N,P,Pk,LPk)
     XLX=Pdot1(N,P,Xk,LXk)
     XLP=Pdot1(N,P,Xk,LPk)
     PLX=Pdot1(N,P,Pk,LXk)
     
     AA=PLP*(PX+XP)-PP*(PLX+XLP)
     BB=2.0*PLP*XX-2.0*XLX*PP
     CC=XX*(PLX+XLP)-XLX*(PX+XP)
     lambda_p=(-BB+SQRT(BB*BB-4.0*AA*CC))/(2.0*AA)
     lambda_m=(-BB-SQRT(BB*BB-4.0*AA*CC))/(2.0*AA)	  
     Lambda=Lambda_P
   END Subroutine RQILineSearch

    FUNCTION Trace2(Matrix,N) RESULT(Tr)
      IMPLICIT NONE
      INTEGER :: I,N
      REAL(DOUBLE) :: Tr
      REAL(DOUBLE),DIMENSION(N,N)::Matrix
      Tr=0D0
      DO I=1,N
         Tr=Tr+Matrix(I,I)
      END DO
    END FUNCTION Trace2
    


   END MODULE RayleighQuotientIteration

