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
!    COMPUTE THE COULOMB MATRIX IN O(N Lg N) CPU TIME
!    Authors: Matt Challacombe and C.J. Tymczak
!------------------------------------------------------------------------------
MODULE JGen
  USE DerivedTypes
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  USE Indexing
  USE InOut
  USE Thresholding
  USE Globals
  USE AtomPairs
  USE BraBloks
  USE PoleTree
  USE TreeWalk
  IMPLICIT NONE
  LOGICAL PrintFlag
!----------!
  CONTAINS !
!=============================================================================================
!
!=============================================================================================
    SUBROUTINE MakeJ(J)
      TYPE(BCSR)             :: J
      TYPE(AtomPair)         :: Pair
      INTEGER                :: AtA,AtB
      INTEGER                :: JP,K,NA,NB,NAB,P,Q,R,I,L
#ifdef PERIODIC        
      INTEGER                :: NCA,NCB
      REAL(DOUBLE)           :: Ax,Ay,Az,Bx,By,Bz
      LOGICAL                :: NotMakeBlock
#endif     
!---------------------------------------------- 
!     Initialize the matrix and associated indecies
      P=1; 
      R=1; 
      J%RowPt%I(1)=1
      CALL SetEq(J%MTrix,Zero)
!     Loop over atom pairs
      J%NAtms=NAtoms
      DO AtA=1,NAtoms
         DO AtB=1,NAtoms
            IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
               NAB = Pair%NA*Pair%NB
               IF(AtB<=AtA)THEN
!              Compute only the lower triangle of symmetric J
#ifdef PERIODIC
                  Ax = Pair%A(1)
                  Ay = Pair%A(2)
                  Az = Pair%A(3)
                  Bx = Pair%B(1)
                  By = Pair%B(2)           
                  Bz = Pair%B(3)
                  NotMakeBlock = .TRUE.
                  DO NCA = 1,CS%NCells
                     Pair%A(1) = Ax+CS%CellCarts%D(1,NCA)
                     Pair%A(2) = Ay+CS%CellCarts%D(2,NCA)
                     Pair%A(3) = Az+CS%CellCarts%D(3,NCA)
                     DO NCB = 1,CS%NCells
                        Pair%B(1) = Bx+CS%CellCarts%D(1,NCB)
                        Pair%B(2) = By+CS%CellCarts%D(2,NCB)
                        Pair%B(3) = Bz+CS%CellCarts%D(3,NCB)
                        Pair%AB2  = (Pair%A(1)-Pair%B(1))**2 &
                                  + (Pair%A(2)-Pair%B(2))**2 &
                                  + (Pair%A(3)-Pair%B(3))**2
                        IF(TestAtomPair(Pair)) THEN
                           J%MTrix%D(R:R+NAB-1)=J%MTrix%D(R:R+NAB-1)+Two*JBlock(Pair,PoleRoot)
                           NotMakeBlock = .FALSE.
                        ENDIF
                     ENDDO
                  ENDDO
#else
                  J%MTrix%D(R:R+NAB-1)=Two*JBlock(Pair,PoleRoot)
#endif
               ENDIF
               J%ColPt%I(P)=AtB
               J%BlkPt%I(P)=R
               R=R+NAB
               P=P+1 
               J%RowPt%I(AtA+1)=P        
               IF(R>MaxNon0.OR.P>MaxBlks) &
                  CALL Halt(' BCSR dimensions blown in J ')
#ifdef PERIODIC
               IF(NotMakeBlock) &
                  CALL Halt(' Making a Zero Block in JGen ')
#endif
            ENDIF
         ENDDO
      ENDDO
      J%NBlks=P-1
      J%NNon0=R-1
!     Fill the upper triangle of J
      DO I=1,NAtoms
         DO JP=J%RowPt%I(I),J%RowPt%I(I+1)-1
            L=J%ColPt%I(JP)
            IF(I>L)THEN
               DO K=J%RowPt%I(L),J%RowPt%I(L+1)-1
                  IF(J%ColPt%I(K)==I)THEN
                     Q=J%BlkPt%I(K)                     
                     EXIT
                  ENDIF
               ENDDO               
               P=J%BlkPt%I(JP)
               NA=BS%BFKnd%I(GM%AtTyp%I(I))
               NB=BS%BFKnd%I(GM%AtTyp%I(L))
               NAB=NA*NB
               CALL XPose(NA,NB,J%MTrix%D(P:P+NAB-1),J%MTrix%D(Q:Q+NAB-1))
            ENDIF
         ENDDO
      ENDDO
!
    END SUBROUTINE MakeJ
!=============================================================================================
!
!=============================================================================================
     SUBROUTINE XPose(M,N,A,AT)
       INTEGER                      :: I,J,M,N,IDex,JDex
       REAL(DOUBLE), DIMENSION(M*N) :: A,AT
       DO I=1,N
          DO J=1,M
             IDex=(I-1)*M+J
             JDex=(J-1)*N+I
             AT(JDex)=A(IDex)
          ENDDO
       ENDDO
     END SUBROUTINE XPose
!=======================================================================================
!
!=======================================================================================
     FUNCTION JBlock(Pair,PoleRoot) RESULT(Jvct)
       TYPE(AtomPair)                           :: Pair
       TYPE(PoleNode), POINTER                  :: PoleRoot
!
       REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)  :: JBlk
       REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB)  :: Jvct
       REAL(DOUBLE),DIMENSION(0:SPLen)          :: SPBraC,SPBraS 
       REAL(DOUBLE)                             :: ZetaA,ZetaB,EtaAB,EtaIn,    &
                                                   XiAB,ExpAB,CA,CB,CC,Ov,     &
                                                   PAx,PAy,PAz,PBx,PBy,PBz,    &
                                                   MDx,MDxy,MDxyz,Amp2,MaxAmp,Tau,OmegaMin
       INTEGER                                  :: KA,KB,CFA,CFB,PFA,PFB,      &
                                                   IndexA,IndexB,              &
                                                   StartLA,StartLB,            &
                                                   StopLA,StopLB
       INTEGER                                  :: I,J,MaxLA,MaxLB,IA,IB,    &
                                                   LMNA,LMNB,LA,LB,MA,MB,    &
                                                   NA,NB,LAB,MAB,NAB,LM,LMN,Ell,EllA,EllB
!-------------------------------------------------------------------------------------- 
       JBlk=Zero
       KA=Pair%KA
       KB=Pair%KB
!----------------------------------
       Prim%A=Pair%A
       Prim%B=Pair%B
       Prim%AB2=Pair%AB2
       Prim%KA=Pair%KA
       Prim%KB=Pair%KB
!----------------------------------
       IndexA=0                  
       DO CFA=1,BS%NCFnc%I(KA)    
       DO CFB=1,BS%NCFnc%I(KB) 
          IndexA  = CFBlokDex(BS,CFA,KA)                
          IndexB  = CFBlokDex(BS,CFB,KB)                  
          StartLA = BS%LStrt%I(CFA,KA)        
          StartLB = BS%LStrt%I(CFB,KB)
          StopLA  = BS%LStop%I(CFA,KA)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLA   = BS%ASymm%I(2,CFA,KA)
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
             IF(TestPrimPair(Prim%Xi,Prim%AB2))THEN
                Prim%PFA=PFA 
                Prim%PFB=PFB
                MaxAmp=SetBraBlok(Prim,BS)
!--------------------------------------------------
!               Evaluate this primitives Ket contribution to J_ab
                HGKet=Zero;SPKetC=Zero;SPKetS=Zero
!               Set MAC and PAC parameters
                Tau=Thresholds%TwoE
                DP2=(MaxAmp/Tau)**(Two/DBLE(SPEll+1))
                PoleSwitch=PFunk(Prim%Ell+MaxL,Tau/MaxAmp)
!               Walk the walk

                CALL JWalk(PoleRoot)
!               Contract <Bra|Ket> bloks to compute matrix elements of J
                IA = IndexA
                DO LMNA=StartLA,StopLA
                   IA=IA+1
                   IB=IndexB
                   EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)                         
                   DO LMNB=StartLB,StopLB
                      IB=IB+1
                      EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)                         
                      Ell=EllA+EllB
                      DO LMN=1,LHGTF(Ell)
                         JBlk(IA,IB)=JBlk(IA,IB)+Phase%D(LMN)*HGBra%D(LMN,IA,IB)*HGKet(LMN)
                      ENDDO

                      CALL HGToSP(Prim,HGBra%D(:,IA,IB),SPBraC,SPBraS)
                      DO LM=0,LSP(Ell)
                         JBlk(IA,IB)=JBlk(IA,IB)+SPBraC(LM)*SPKetC(LM) &
                                                +SPBraS(LM)*SPKetS(LM)
                      ENDDO
                  ENDDO
                  ENDDO
             ENDIF !End primitive thresholding
          ENDDO 
          ENDDO
       ENDDO
       ENDDO
       Jvct=BlockToVect(Pair%NA,Pair%NB,Jblk)
     END FUNCTION JBlock
!====================================================================================================
!
!====================================================================================================
END MODULE JGen
