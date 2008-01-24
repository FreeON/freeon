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
!   GENERAL CODE FOR COMPUTING THE ELECTRON DENSITY ON THE CUBE GRID 
!------------------------------------------------------------------------------
         DO IGrid=1,NGrid
            RQx=Grid(IGrid,1)-Node%Qx            
            RQy=Grid(IGrid,2)-Node%Qy
            RQz=Grid(IGrid,3)-Node%Qz
            RQ2=RQx*RQx+RQy*RQy+RQz*RQz
            Xpt=EXP(-Node%Zeta*RQ2)
            IF(Xpt>1.D-16)THEN
               TwoZ=Two*Node%Zeta
               LambdaX(0)=Xpt
               LambdaY(0)=One
               LambdaZ(0)=One
               LambdaX(1)=TwoZ*RQx*Xpt
               Lambday(1)=TwoZ*RQy
               LambdaZ(1)=TwoZ*RQz
               DO L=2,Node%Ell+1
                  L1=L-1
                  L2=L-2
                  RL1=DBLE(L1)
                  LambdaX(L)=TwoZ*(RQx*LambdaX(L1)-RL1*LambdaX(L2))
                  LambdaY(L)=TwoZ*(RQy*LambdaY(L1)-RL1*LambdaY(L2))
                  LambdaZ(L)=TwoZ*(RQz*LambdaZ(L1)-RL1*LambdaZ(L2))
               ENDDO
               DO L=0,Node%Ell
                  DO M=0,Node%Ell-L
                     DO N=0,Node%Ell-M-L
                        DO iSDen=0,NSDen-1
                           LMN=LMNDex(L,M,N)+iSDen*OffSDen
                           I0=IGrid+iSDen*NGrid
                           !
                           RhoV(I0,1)=RhoV(I0,1)+LambdaX(L  )*LambdaY(M  )*LambdaZ(N  )*Node%Co(LMN)
                           RhoV(I0,2)=RhoV(I0,2)-LambdaX(L+1)*LambdaY(M  )*LambdaZ(N  )*Node%Co(LMN)
                           RhoV(I0,3)=RhoV(I0,3)-LambdaX(L  )*LambdaY(M+1)*LambdaZ(N  )*Node%Co(LMN)
                           RhoV(I0,4)=RhoV(I0,4)-LambdaX(L  )*LambdaY(M  )*LambdaZ(N+1)*Node%Co(LMN)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
