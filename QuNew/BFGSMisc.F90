!------------------------------------------------------------------------------
!--  This source code is part of the MondoSCF suite of programs for 
!    linear scaling electronic structure theory and ab initio molecular 
!    dynamics.
!
!--  Matt Challacombe
!    Los Alamos National Laboratory
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
!    SUPPORTING ROUTINES FOR BFGSHess
!    Author(s):  Matt Challacombe
!------------------------------------------------------------------------------
MODULE BFGSMisc
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalCharacters
   USE InOut
   USE PrettyPrint
   USE MemMan
   USE Parse
   USE Macros
   USE SetXYZ
   USE LinAlg
   USE Functionals
   USE MatFunk
   USE AtomPairs
#ifdef PARALLEL
   USE MondoMPI
#endif
   CONTAINS 
!--------------------------------------------------------------
!     Find min diagonal element of B[i+1]
!--------------------------------------------------------------
      FUNCTION MinDiagEl(B,DXB,DG,D1,D2)
         TYPE(BCSR)                  :: B
         TYPE(DBL_VECT)              :: DG,DXB
         REAL(DOUBLE),DIMENSION(3,3) :: AdBlk
         REAL(DOUBLE)                :: D1,D2,MinDiagEl
         INTEGER                     :: AtA,JP,JG,KP,KG, &
                                        PDXB,QDXB,PB,OA,OB
!-------------------------------------------------------------------
!        Find minimum diagonal element of B
!-------------------------------------------------------------------
         MinDiagEl=BIG_DBL
         DO AtA=1,NAtoms
            OA=OffS%I(AtA)
            DO JP=B%RowPt%I(AtA),B%RowPt%I(AtA+1)-1
               JG=B%ColPt%I(JP)
!              Find diagonal blocks of B
               IF(AtA==JG)THEN
                  OB=OffS%I(JG)            
                  PB=B%BlkPt%I(JP)
                  AdBlk=BFGSBlk(DG%D(OA:),DG%D(OB:),DXB%D(OA:),DXB%D(OB:),D1,D2,B%MTrix%D(PB:))
                  DO I=1,3            
                     MinDiagEl=MIN(MinDiagEl,ABS(AdBlk(I,I)))
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      END FUNCTION MinDiagEl
!-------------------------------------------------------------------
!     Sparse update of the BFGS Hessian 
!-------------------------------------------------------------------
      SUBROUTINE ReUpB(GM,BS,BOld,DG,DXB,D1,D2,BNew)
         TYPE(CRDS)                  :: GM
         TYPE(BSET)                  :: BS 
         TYPE(BCSR)                  :: BOld,BNew
         TYPE(DBL_VECT)              :: DG,DXB
         REAL(DOUBLE),DIMENSION(3,3) :: AdBlkd
         REAL(DOUBLE)                :: D1,D2,MinDiagEl
         INTEGER                     :: AtA,AtB,JP,JG,KP,KG, &
                                        PDXB,QDXB,PB,OA,OB, &
                                        P,Q,R,S,M,N,MN,MN1
         TYPE(INT_VECT)              :: Flag
         REAL(DOUBLE),DIMENSION(3,3) :: AdBlk
         TYPE(AtomPair)            :: Pair
!---------------------------------------------------------------------
         Q=1
         R=1
         CALL New(Flag,NAtoms+1)
         Flag%I=0
         BNew%NAtms=BOld%NAtms
         DO AtA=1,BOld%NAtms
            BNew%RowPt%I(AtA)=R
            OA=OffS%I(AtA)
            M=BSiz%I(AtA)
!           Copy 
            DO JP=BOld%RowPt%I(AtA),BOld%RowPt%I(AtA+1)-1
               JG=BOld%ColPt%I(JP)
               P=BOld%BlkPt%I(JP)
               MN=M*BSiz%I(JG)
               MN1=MN-1
               BNew%MTrix%D(Q:Q+MN1)=BOld%MTrix%D(P:P+MN1)
               BNew%ColPt%I(R)=JG
               BNew%BlkPt%I(R)=Q
               Flag%I(JG)=Q
               Q=Q+MN
               R=R+1
            ENDDO
!           Update
            DO AtB=1,NAtoms
               IF(SetAtomPair(GM,BS,AtA,AtB,Pair))THEN
                  OB=OffS%I(AtB)
                  S=Flag%I(AtB)    
                  N=BSiz%I(AtB)
                  MN=M*N
                  MN1=MN-1
                  IF(S/=0)THEN
                     AdBlk=BFGSBlk(DG%D(OA:),DG%D(OB:),DXB%D(OA:),DXB%D(OB:),D1,D2,BNew%MTrix%D(S:))
                     BNew%MTrix%D(S:S+MN1)=BlockToVect(M,N,AdBlk) 
                  ELSE                     
                     AdBlk=BFGSBlk(DG%D(OA:),DG%D(OB:),DXB%D(OA:),DXB%D(OB:),D1,D2)
!                    Add new block? 
                     IF(SQRT(DBL_Dot(MN,AdBlk,AdBlk))>Thresholds%Trix)THEN
                        BNew%MTrix%D(Q:Q+MN1)=BlockToVect(M,N,AdBlk) 
                        BNew%ColPt%I(R)=AtB
                        BNew%BlkPt%I(R)=Q
                        Flag%I(AtB)=Q
                        Q=Q+MN
                        R=R+1
                     ENDIF          
                  ENDIF
               ENDIF
            ENDDO
            DO KL=BNew%RowPt%I(AtA),R-1
               Flag%I(BNew%ColPt%I(KL))=0
            ENDDO
         ENDDO
         BNew%NBlks=R-1
         BNew%NNon0=Q-1              
         BNew%RowPt%I(BNew%NAtms+1)=R
         CALL Delete(Flag)
      END SUBROUTINE ReUpB
!----------------------------------------------------------------
!     Block BFGS update 
!----------------------------------------------------------------
      FUNCTION BFGSBlk(DGA,DGB,DXBA,DXBB,D1,D2,B_O)
         REAL(DOUBLE)                :: D1,D2
         REAL(DOUBLE),DIMENSION(:)   :: DGA,DGB,DXBA,DXBB
         REAL(DOUBLE),DIMENSION(3,3) :: BFGSBlk
         REAL(DOUBLE),DIMENSION(3,3),OPTIONAL :: B_O
!-----------------------------------------------------------------
         IF(PRESENT(B_O))THEN
            DO I=1,3
               DO J=1,3
                  BFGSBlk(I,J)=B_O(I,J)+DGA(I)*DGB(J)/D1-DXBA(I)*DXBB(J)/D2
              ENDDO
            ENDDO
          ELSE
            DO I=1,3
               DO J=1,3
                  BFGSBlk(I,J)=DGA(I)*DGB(J)/D1-DXBA(I)*DXBB(J)/D2
              ENDDO
            ENDDO
          ENDIF
      END FUNCTION BFGSBlk
END MODULE 

