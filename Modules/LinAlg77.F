C------------------------------------------------------------------------------
C    This code is part of the MondoSCF suite of programs for linear scaling
C    electronic structure theory and ab initio molecular dynamics.
C
C    Copyright (2004). The Regents of the University of California. This
C    material was produced under U.S. Government contract W-7405-ENG-36
C    for Los Alamos National Laboratory, which is operated by the University
C    of California for the U.S. Department of Energy. The U.S. Government has
C    rights to use, reproduce, and distribute this software.  NEITHER THE
C    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
C    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
C
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by the
C    Free Software Foundation; either version 2 of the License, or (at your
C    option) any later version. Accordingly, this program is distributed in
C    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
C    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
C    PURPOSE. See the GNU General Public License at www.gnu.org for details.
C
C    While you may do as you like with this software, the GNU license requires
C    that you clearly mark derivative software.  In addition, you are encouraged
C    to return derivative works to the MondoSCF group for review, and possible
C    disemination in future releases.
C------------------------------------------------------------------------------
C=========================================================================
C     Generic F77 style (D,B)CSR symbolic matrix multiply: C=C+Beta*A.B
C=========================================================================
C23456789012345678901234567890123456789012345678901234567890123456789012
      INTEGER FUNCTION SymbolikMM_GENERIC77(NAtoms,
     >                                      ANAtms,ANBlks,ANAtoms,
     >                                      BNAtms,BNBlks,BNAtoms,
     >                                      DNAtms,DNBlks,MxBlks,MxNon0,
     >                                      CNAtms,CNBlks,CNNon0,UpDate,
     >                                      AOffSt,ARowPt,AColPt,
     >                                      BRowPt,BColPt,
     >                                      CRowPt,CColPt,CBlkPt,
     >                                      DRowPt,DColPt,DBlkPt,
     >                                      BSiz,Flag)
         IMPLICIT NONE
         INTEGER NAtoms,ANAtms,ANBlks,ANAtoms,BNAtms,BNBlks,BNAtoms,
     >           DNAtms,DNBlks,MxBlks,MxNon0,CNAtms,CNBlks,CNNon0,
     >           AOffSt,JP,KP,S,T,IL,IG,JG,KG,MA,MB,MM,
     >           IStrtA,IStopA,IStrtB,IStopB,IStrtD,IStopD
         INTEGER ARowPt(ANAtoms+1),AColPt(ANBlks),
     >           BRowPt(BNAtoms+1),BColPt(BNBlks),
     >           CRowPt(CNAtms+1),CColPt(MxBlks),CBlkPt(MxBlks),
     >           DRowPt(DNAtms+1),DColPt(DNBlks),DBlkPt(DNBlks),
     >           BSiz(NAtoms),Flag(NAtoms)
         LOGICAL UpDate
C-----------------------------------------------------------------------------------
         T=1
         IF(UpDate)THEN
            S=CNNon0+1
            CALL INT_VECT_EQ_INT_VECT(CNAtms+1,DRowPt,CRowPt)
            CALL INT_VECT_EQ_INT_VECT(CNBlks  ,DColPt,CColPt)
            CALL INT_VECT_EQ_INT_VECT(CNBlks  ,DBlkPt,CBlkPt)
         ELSE
            S=1
            CNAtms=ANAtms
         ENDIF
         DO IL=1,ANAtms
            IG=IL+AOffSt
            CRowPt(IL)=T
            MA=BSiz(IG)
            IStrtA=ARowPt(IG)
            IStopA=ARowPt(IG+1)-1
            IF(UpDate)THEN
               IStrtD=DRowPt(IL)
               IStopD=DRowPt(IL+1)-1
               DO JP=IStrtD,IStopD
                  JG=DColPt(JP)
                  IF(Flag(JG).EQ.0)THEN
                     Flag(JG)=1
                     CColPt(T)=JG
                     CBlkPt(T)=DBlkPt(JP)
                     T=T+1
                  ENDIF
               ENDDO
            ENDIF
            DO JP=IStrtA,IStopA
               JG=AColPt(JP)
               IStrtB=BRowPt(JG)
               IStopB=BRowPt(JG+1)-1
               IF(IStrtB.NE.0.AND.IStopB.NE.0)THEN
                  MB=BSiz(JG)
                  MM=MA*MB
                  DO KP=IStrtB,IStopB
                     KG=BColPt(KP)
                     IF(Flag(KG).EQ.0)THEN
                        CColPt(T)=KG
                        CBlkPt(T)=S
                        Flag(KG) =1
                        T=T+1
                        S=S+MA*BSiz(KG)
                        IF(T.GT.MxBlks.OR.S.GT.MxNon0)THEN
                           SymbolikMM_Generic77=-1
                           RETURN
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            DO KP=CRowPt(IL),T-1
               Flag(CColPt(KP))=0
            ENDDO
         ENDDO
         CNNon0=S-1
         CNBlks=T-1
         CRowPt(CNAtms+1)=T
         SymbolikMM_Generic77=0
         RETURN
      END !FUNCTION SymbolikMM_Generic
!======================================================================
!     Generic F77 style (D,B)CSR numeric matrix multiply: C=C+A.B
!======================================================================
      SUBROUTINE NumerikMM_GENERIC77(ANAtms,AOffSt,
     >                               ARowPt,AColPt,ABlkPt,AMTrix,
     >                               BRowPt,BColPt,BBlkPt,BMTrix,
     >                               CRowPt,CColPt,CBlkPt,CMTrix,
     >                               BSiz,Flag,FlOp)

         IMPLICIT NONE
         REAL*8  AMTrix(*),BMTrix(*),CMTrix(*)
         REAL*8  FlOp,Op
         INTEGER ANAtms,AOffSt,
     >           K,JP,KP,P,Q,R,IL,IG,JG,KG,
     >           MA,MB,NB,MNA,IStrtA,IStopA,
     >           IStrtB,IStopB,IStrtC,IStopC
         INTEGER ARowPt(*),AColPt(*),ABlkPt(*),BRowPt(*),BColPt(*),
     >           BBlkPt(*),BSiz(*),CRowPt(*),CColPt(*),CBlkPt(*),Flag(*)
!------------------------------------------------------------------------------
         Op=0.0D0
         DO IL=1,ANAtms
            IG=IL+AOffSt
            MA=BSiz(IG)
            IStrtA=ARowPt(IG)
            IStopA=ARowPt(IG+1)-1
            IStrtC=CRowPt(IL)
            IStopC=CRowPt(IL+1)-1
            DO JP=IStrtC,IStopC
               Flag(CColPt(JP))=CBlkPt(JP)
            ENDDO
            DO JP=IStrtA,IStopA
               JG=AColPt(JP)
               IStrtB=BRowPt(JG)
               IStopB=BRowPt(JG+1)-1
               IF(IStrtB.NE.0.AND.IStopB.NE.0)THEN
                  MB=BSiz(JG)
                  MNA=MA*MB
                  P=ABlkPt(JP)
                  DO KP=IStrtB,IStopB
                     KG=BColPt(KP)
                     Q=BBlkPt(KP)
                     R=Flag(KG)
                     NB=BSiz(KG)
                     Op=Op+DBLE(MNA*NB)
                     CALL DGEMM_NN(MA,MB,NB,1.0D0,AMTrix(P),
     >                             BMTrix(Q),CMTrix(R))
                  ENDDO
               ENDIF
            ENDDO
            DO K=IStrtC,IStopC
               Flag(CColPt(K))=0
            ENDDO
         ENDDO
         FlOp=FlOp+2.0D0*Op
      END !SUBROUTINE NumerikMM_GENERIC77
#ifdef PARALLEL
!===============================================================================
!     Post non-blocking ISend and IRecv for DBCSR data exchange
!===============================================================================
      INTEGER FUNCTION PostRecv77(NAtoms,MyId,From,Tag,COMM,ANAtms,
     >                       ANBlks,BNSMat,BNBlks,UNAtms,UNBlks,UNNon0,
     >                       ARowPt,AColPt,BGRwPt,BGClPt,
     >                       URowPt,UColPt,UBlkPt,UMTrix,
     >                       Flag,BSiz,FromBeg,FromEnd,
     >                       MyOff,FromOff)
         IMPLICIT NONE
         INTEGER I,S,T,IG,IL,KG,KP,JL,JP,JG,JH,
     >           NAtoms,MyId,From,Tag,COMM,
     >           ANAtms,ANBlks,BNSMat,BNBlks,UNAtms,UNBlks,UNNon0,
     >           FromBeg,FromEnd,MyOff,FromOff,IErr
         REAL*8 UMTrix(UNNon0)
         INTEGER ARowPt(ANAtms+1),AColPt(ANBlks),
     >           BGRwPt(NAtoms+1),BGClPt(BNBlks),
     >           URowPt(UNAtms+1),UColPt(UNBlks),UBlkPt(UNBlks),
     >           Flag(NAtoms),BSiz(NAtoms)
         INTEGER BIG_INT
         PARAMETER(BIG_INT=100000)
         INCLUDE 'mpif.h'
C
         DO I=1,NAtoms !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Flag(I)=0
         ENDDO
!--------------------------------------------------------------------------
!        Find overlaping row-col blocks and compute pointers
!
         T=1
         S=1
         CALL INT_VECT_EQ_INT_SCLR(UNAtms,URowPt,BIG_INT)
         DO JH=FromBeg,FromEnd
            DO IL=1,ANAtms
               IG=IL+MyOff
               DO JP=ARowPt(IL),ARowPt(IL+1)-1
                  JG=AColPt(JP)
                  IF(JG.EQ.JH)THEN
                     JL=JG-FromOff
                     IF(Flag(JL).EQ.0)THEN
                        Flag(JL)=1
                        URowPt(JL)=T
                        DO KP=BGRwPt(JG),BGRwPt(JG+1)-1
                           KG=BGClPt(KP)
                           UColPt(T)=KG
                           UBlkPt(T)=S
                           T=T+1
                           S=S+BSiz(JG)*BSiz(KG)*BNSMat!<<<<< add spin here?
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         CALL INT_VECT_EQ_INT_SCLR(NAtoms,Flag,0)
C!-------------------------------------------------
         URowPt(UNAtms+1)=T
         DO I=UNAtms+1,2,-1
            IF(URowPt(I-1).EQ.BIG_INT)THEN
               URowPt(I-1)=URowPt(I)
            ENDIF
         ENDDO
C!-----------------------------------------------------------------------
         CALL MPI_IRECV(UMTrix,UNNon0,MPI_DOUBLE_PRECISION,From,Tag,
     >                  COMM,PostRecv77,IErr)
C
         RETURN
      END ! FUNCTION PostRecv77
C=========================================================================
C
C=========================================================================
      INTEGER FUNCTION PostSend77(NAtoms,ANBlks,BNSMat,BNAtms,BNBlks,
     >                            BNNon0,V,To,Tag,MyId,VType,COMM,
     >                            AGrwPt,AGClPt,BRowPt,BColPt,
     >                            BBlkPt,BMTrix,VBlks,VDisp,Flag,
     >                            BSiz,MyOff,MyBeg,MyEnd,ToBeg,ToEnd)
         IMPLICIT NONE
         INTEGER NAtoms,ANBlks,BNSMat,BNAtms,BNBlks,BNNon0,V,To,
     >           Tag,MyId,VType,COMM,MyOff,MyBeg,MyEnd,ToBeg,ToEnd
         REAL*8  BMTrix(BNNon0)
         INTEGER IErr,NAtms,I,K,MN,S,T,TT,JP,KP,IG,JG,JL,JH,KG
         INTEGER BSiz(NAtoms),Flag(NAtoms),VBlks(V),VDisp(V),
     >           AGrwPt(NAtoms+1),AGClPt(ANBlks),BRowPt(BNAtms+1),
     >           BColPt(BNBlks),BBlkPt(BNBlks)
C        INCLUDE 'MONDO_MPI_INCLUDE.Inc'
         INCLUDE 'mpif.h'
C
         DO I=1,NAtoms !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Flag(I)=0
         ENDDO


C!-----------------------------------------------------------------------
C!        Find overlaping row-col blocks in BMTrix and mark for SEND
C!
         T=0
         S=1
         NAtms=MyEnd-MyBeg+1
         DO JH=MyBeg,MyEnd
            DO IG=ToBeg,ToEnd
               DO JP=AGRwPt(IG),AGRwPt(IG+1)-1
                  JG=AGClPt(JP)
                  IF(JG.EQ.JH)THEN
                     JL=JG-MyOff
                     IF(Flag(JL).EQ.0)THEN
                        Flag(JL)=1
                        DO KP=BRowPt(JL),BRowPt(JL+1)-1
                           KG=BColPt(KP)
                           MN=BSiz(JG)*BSiz(KG)*BNSMat  !<<<<< add spin here?
C!                          Check for contigous blocks ..
C                          IF(T.NE.0.AND.VDisp(T)+MN.EQ.BBlkPt(KP)-1)THEN
C                              VBlks(T)=VBlks(T)+MN
C                           ELSE
                              T=T+1
                              VDisp(T)=BBlkPt(KP)-1
                              VBlks(T)=MN
C                           ENDIF
                           S=S+MN
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
C
         CALL INT_VECT_EQ_INT_SCLR(NAtoms,Flag,0)
         CALL MPI_TYPE_INDEXED(T,VBlks,VDisp,
     >        MPI_DOUBLE_PRECISION,VType,IErr)
         CALL MPI_TYPE_COMMIT(VType,IErr)
#ifdef NONBLOCKING
C        Post a non-blocking send
         CALL MPI_ISEND(BMTrix,1,VType,To,Tag,COMM,PostSend77,IErr)
#else
C        Post a blocking send
         CALL MPI_SEND(BMTrix,1,VType,To,Tag,COMM,IErr)
         PostSend77=MPI_REQUEST_NULL
#endif
         RETURN
      END !FUNCTION PostSend77
C==========================================================
C
C==========================================================

      SUBROUTINE WaitAll77(N,Reqs,Stats)
C        INCLUDE 'MONDO_MPI_INCLUDE.Inc'
         INCLUDE 'mpif.h'
         INTEGER I,J,N,M,IErr
         INTEGER Reqs(N),Stats(MPI_STATUS_SIZE,N)
C------------------------------------------------------------
         DO I=1,N
            IF(Reqs(I).NE.MPI_REQUEST_NULL)GOTO 10
         ENDDO
         RETURN
 10      CONTINUE
         WRITE(*,*)' WAITING ....'
         CALL MPI_WAITALL(N,Reqs,Stats,IErr)
         WRITE(*,*)' DONE WAITING IN WAITALL '
         RETURN
      END













#endif
C==========================================================
C
C==========================================================
      REAL*8 FUNCTION BlkTrace_1(M,A)
         IMPLICIT NONE
         INTEGER I,M
         REAL*8  A(M,M)
         BlkTrace_1=0.0D0
         DO I=1,M
            BlkTrace_1=BlkTrace_1+A(I,I)
         ENDDO
      RETURN
      END
C==========================================================
C
C==========================================================
      REAL*8 FUNCTION BlkTrace_2(M,N,A,B)
         IMPLICIT NONE
         INTEGER I,K,M,N
         REAL*8  A(M,N),B(N,M)
         BlkTrace_2=0.0D0
         DO I=1,M
            DO K=1,N
               BlkTrace_2=BlkTrace_2+A(I,K)*B(K,I)
            ENDDO
         ENDDO
      RETURN
      END
C=================================================================================
C     Scaling of a REAL*8 vector by a constant: A=b*A
C=================================================================================
      SUBROUTINE DBL_Scale(N,A,B)
         IMPLICIT NONE
         INTEGER  I,I1,I2,I3,I4,I5,I6,M,N
         REAL*8   A(N),B
!---------------------------------------------------------------------------
         M=MOD(N,7)
         DO I=1,M
            A(I)=A(I)*B
         ENDDO
         IF(N.LT.7)RETURN
         DO I=M+1,N,7
            I1=I+1
            I2=I+2
            I3=I+3
            I4=I+4
            I5=I+5
            I6=I+6
            A(I )=A(I )*B
            A(I1)=A(I1)*B
            A(I2)=A(I2)*B
            A(I3)=A(I3)*B
            A(I4)=A(I4)*B
            A(I5)=A(I5)*B
            A(I6)=A(I6)*B
         ENDDO
      END !SUBROUTINE DBL_Scale
C==========================================================
C     Dot product of two REAL*8 vectors: Dot=(A,B)
C==========================================================
      REAL*8 FUNCTION DBL_Dot(N,A,B)
      IMPLICIT NONE
      INTEGER N
      REAL*8  A(N),B(N)
      INTEGER M,I,I1,I2,I3,I4,I5,I6
C----------------------------------------------
      M=MOD(N,7)
      DBL_Dot=0.0D0
      DO I=1,M
         DBL_Dot=DBL_Dot+A(I)*B(I)
      ENDDO
      IF(N.LT.7)RETURN
      DO I=M+1,N,7
         I1=I+1
         I2=I+2
         I3=I+3
         I4=I+4
         I5=I+5
         I6=I+6
         DBL_Dot=DBL_Dot+A(I )*B(I )
         DBL_Dot=DBL_Dot+A(I1)*B(I1)
         DBL_Dot=DBL_Dot+A(I2)*B(I2)
         DBL_Dot=DBL_Dot+A(I3)*B(I3)
         DBL_Dot=DBL_Dot+A(I4)*B(I4)
         DBL_Dot=DBL_Dot+A(I5)*B(I5)
         DBL_Dot=DBL_Dot+A(I6)*B(I6)
      ENDDO
      RETURN
      END
!====================================================================
!     Add a scalar to the diagonal of a symmetric block
!====================================================================
      SUBROUTINE AddToDiag(N,A,B)
         IMPLICIT NONE
         INTEGER I,N
         REAL*8 A(N,N),B
         DO I=1,N
            A(I,I)=A(I,I)+B
         ENDDO
         RETURN
      END !SUBROUTINE AddToDiag
