!------------------------------------------------------------------------------
!--  This code is part of the MondoSCF suite of programs for linear scaling 
!    electronic structure theory and ab initio molecular dynamics.
!
!--  Copyright (c) 2001, the Regents of the University of California.  
!    This SOFTWARE has been authored by an employee or employees of the 
!    University of California, operator of the Los Alamos National Laboratory 
!    under Contract No. W-7405-ENG-36 with the U.S. Department of Energy.  
!    The U.S. Government has rights to use, reproduce, and distribute this 
!    SOFTWARE.  The public may copy, distribute, prepare derivative works 
!    and publicly display this SOFTWARE without charge, provided that this 
!    Notice and any statement of authorship are reproduced on all copies.  
!    Neither the Government nor the University makes any warranty, express 
!    or implied, or assumes any liability or responsibility for the use of 
!    this SOFTWARE.  If SOFTWARE is modified to produce derivative works, 
!    such modified SOFTWARE should be clearly marked, so as not to confuse 
!    it with the version available from LANL.  The return of derivative works
!    to the primary author for integration and general release is encouraged. 
!    The first publication realized with the use of MondoSCF shall be
!    considered a joint work.  Publication of the results will appear
!    under the joint authorship of the researchers nominated by their
!    respective institutions. In future publications of work performed
!    with MondoSCF, the use of the software shall be properly acknowledged,
!    e.g. in the form "These calculations have been performed using MondoSCF, 
!    a suite of programs for linear scaling electronic structure theory and
!    ab initio molecular dynamics", and given appropriate citation.  
!------------------------------------------------------------------------------
!    Author: Matt Challacombe
!   GENERAL CODE FOR COMPUTING THE ELECTRON POPULATION (AND POSSIBLY 
!   DERIVATIVES) NODE BY NODE
!------------------------------------------------------------------------------
             Z=Node%Zeta
             TwoZ=Two*Z
             SqZ=SQRT(Z)
             CoFact=SqrtPi/(Two*SqZ)
!            Lower bounds
             LQx=Box%BndBox(1,1)-Node%Qx
             LQy=Box%BndBox(2,1)-Node%Qy
             LQz=Box%BndBox(3,1)-Node%Qz
             LQ2=LQx*LQx+LQy*LQy+LQz*LQz
!            Lower Lambdas
             LLambdaX(0)=-CoFact*ERF(SqZ*LQx)
             LLambdaY(0)=-CoFact*ERF(SqZ*LQy)
             LLambdaZ(0)=-CoFact*ERF(SqZ*LQz)
             LLambdaX(1)=EXP(-Z*LQx*LQx)
             LLambdaY(1)=EXP(-Z*LQy*LQy)
             LLambdaZ(1)=EXP(-Z*LQz*LQz)
!            Upper bounds
             UQx=Box%BndBox(1,2)-Node%Qx
             UQy=Box%BndBox(2,2)-Node%Qy
             UQz=Box%BndBox(3,2)-Node%Qz
             UQ2=UQx*UQx+UQy*UQy+UQz*UQz
!            Upper Lambdas
             ULambdaX(0)=-CoFact*ERF(SqZ*UQx)
             ULambdaY(0)=-CoFact*ERF(SqZ*UQy)
             ULambdaZ(0)=-CoFact*ERF(SqZ*UQz)
             ULambdaX(1)=EXP(-Z*UQx*UQx)             
             ULambdaY(1)=EXP(-Z*UQy*UQy)             
             ULambdaZ(1)=EXP(-Z*UQz*UQz)             
!            Generic Lambdas
             DO L=2,Node%Ell+1
                L1=L-1
                L2=L-2
                RL2=DBLE(L2)
                LLambdaX(L)=TwoZ*(LQx*LLambdaX(L1)-RL2*LLambdaX(L2))
                LLambdaY(L)=TwoZ*(LQy*LLambdaY(L1)-RL2*LLambdaY(L2))
                LLambdaZ(L)=TwoZ*(LQz*LLambdaZ(L1)-RL2*LLambdaZ(L2))
                ULambdaX(L)=TwoZ*(UQx*ULambdaX(L1)-RL2*ULambdaX(L2))
                ULambdaY(L)=TwoZ*(UQy*ULambdaY(L1)-RL2*ULambdaY(L2))
                ULambdaZ(L)=TwoZ*(UQz*ULambdaZ(L1)-RL2*ULambdaZ(L2))
             ENDDO
             DO L=0,Node%Ell+1
                LambdaX(L)=LLambdaX(L)-ULambdaX(L)
                LambdaY(L)=LLambdaY(L)-ULambdaY(L)
                LambdaZ(L)=LLambdaZ(L)-ULambdaZ(L)
             ENDDO
             DO L=0,Node%Ell
                DO M=0,Node%Ell-L
                   DO N=0,Node%Ell-M-L
                      LMN=LMNDex(L,M,N)
                      Pop=Pop+Node%Co(LMN)*LambdaX(L)*LambdaY(M)*LambdaZ(N)
                      dPopX=dPopX-Node%Co(LMN)*LambdaX(L+1)*LambdaY(M)*LambdaZ(N)
                      dPopY=dPopY-Node%Co(LMN)*LambdaX(L)*LambdaY(M+1)*LambdaZ(N)
                      dPopZ=dPopZ-Node%Co(LMN)*LambdaX(L)*LambdaY(M)*LambdaZ(N+1)
                   ENDDO
                ENDDO
             ENDDO 
!===============================================================================
!
