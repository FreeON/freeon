!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!    COMPUTE A SPARSE HESSIAN USING BFGS UPDATES
!
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
!--------------------------------------------------------------
!        DXB=[/\X].[B]
!
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
!               IF(SetAtomPair(GM,BS,AtA,AtB,Pair))THEN
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
!               ENDIF
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

PROGRAM BFGSHess
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
#ifdef PARALLEL
   USE MondoMPI
#endif
   USE BFGSMisc
   IMPLICIT NONE
   TYPE(ARGMT)                :: Args
   CHARACTER(LEN=8),PARAMETER :: Prog='BFGSHess'
   TYPE(AtomPair)             :: Pair
   TYPE(BSet)                 :: BS
   TYPE(CRDS)                 :: GM0,GM1
   TYPE(BCSR)                 :: B0,B1,sDX,sDXB
   TYPE(DBL_VECT)             :: G0,G1,X0,X1,DX,DG,DXB
   TYPE(DBL_RNK2)             :: X, G, DnsB,BInv
   INTEGER                    :: N3,I,J,K,IGeom
   REAL(DOUBLE)               :: D1,D2,F,AdEl,MinEl,A
   REAL(DOUBLE),DIMENSION(6)  :: AA
!----------------------------------------------------------------------------
!   CALL StartUp(Args,Prog)
!
!   IGeom=Args%I%I(3)
!
!  Override default blocking...
   NAtoms=2
   N3=3*NAtoms
   NBasF=N3
   MaxAtms=NAtoms+1
   MaxBlks=NAtoms**2+1
   MaxNon0=N3**2+1
   PrintFlags%Mat=DEBUG_MATRICES
   CALL SetDSYEVWork(N3)
!
   CALL New(BSiz,NAtoms)
   CALL New(OffS,NAtoms)
!
   BSiz%I=3
   OffS%I(1)=1
   DO I=2,NAtoms
      OffS%I(I)=OffS%I(I-1)+3
   ENDDO      
! 
   CALL New(G0,N3)
   CALL New(G1,N3)
   CALL New(X0,N3)
   CALL New(X1,N3)
   CALL New(DG,N3)
   CALL New(DX,N3)
   CALL New(DXB,N3)
   CALL New(B0)
   CALL New(B1)
   CALL New(sDX)
   CALL New(sDXB)
   CALL New(X,(/N3,1000/))
   CALL New(G,(/N3,1000/))
   CALL New(DnsB,(/N3,N3/))
   CALL New(BInv,(/N3,N3/))
!
   X%D(:,1)=(/1.D0,2.D0,3.D0,4.D0,5.D0,6.D0/)

   AA=(/1.D6,1.D4,1.D2,1.D0,1.D-2,1.D-4/)

   WRITE(*,*)' 1 '

   DO IGeom=1,30

      F=Zero
      DO J=1,6
         F=F+J*AA(J)*X%D(J,IGeom)**2
         G%D(J,IGeom)=2.0D0*AA(J)*X%D(J,IGeom)
      ENDDO

   IF(IGeom==1)THEN
      CALL SetToI(B1)
      GOTO 10
   ENDIF
!
   X0%D=X%D(:,IGeom-1)
   X1%D=X%D(:,IGeom)
   G0%D=G%D(:,IGeom-1)
   G1%D=G%D(:,IGeom)
!
   DX%D=X1%D-X0%D
   DG%D=G1%D-G0%D
!
   CALL PPrint(DX,'DX['//TRIM(IntToChar(IGeom))//']',Unit_O=6)   
   CALL PPrint(DG,'DG['//TRIM(IntToChar(IGeom))//']',Unit_O=6)
!
!   CALL Get(B,'Hessian')
!
   D1=DOT_PRODUCT(DX%D,DG%D)             !  D1=[/\G^T].[/\X]

   CALL Set_BCSR_EQ_VECT(sDX,DX)         !  sDX=SparseM[DX]


   CALL Multiply(B0,sDX,sDXB)            !  sDXB=[/\X^T].[B_I]

   CALL Set_VECT_EQ_BCSR(DXB,sDXB)       !  DXB=DenseV[sDXB]

   D2=DOT_PRODUCT(DX%D,DXB%D)             !  D2=[/\X^T].[B_I].[/\X]

   MinEl=MinDiagEl(B0,DXB,DG,D1,D2)      !

   Thresholds%Trix=1.D-10 !*MinEl        !

   CALL ReUpB(GM1,BS,B0,DG,DXB,D1,D2,B1) ! [B_I+1]=BFGS([B_I],[/\X],[/\G])
!
10 CONTINUE

   CALL SetEq(B0,B1)

   CALL PPrint(B0,'B['//TRIM(IntToChar(IGeom+1))//']',Unit_O=6)

   CALL SetEq(DnsB,B1)

   CALL FunkOnSqMat(N3,Inverse,DnsB%D,BInv%D)

   X%D(:,IGeom+1)=X%D(:,IGeom)-MATMUL(BInv%D,G%D(:,IGeom))

   CALL Filter(B0,B1)
!
   ENDDO
!
   CALL Delete(G0)
   CALL Delete(G1)
   CALL Delete(X0)
   CALL Delete(X1)
   CALL Delete(DG)
   CALL Delete(DX)
   CALL Delete(DXB)
   CALL Delete(B0)
   CALL Delete(B1)
!   CALL Shutdown(Prog)
!
END PROGRAM 
