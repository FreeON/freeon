!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 1999, The University of California
!
!    INDEXING FUNCTIONS
!
MODULE Indexing
   USE DerivedTypes
   IMPLICIT NONE
   CONTAINS

      FUNCTION LHGTF(L)
         INTEGER, INTENT(IN) :: L
         INTEGER :: LHGTF
         LHGTF=(L+1)*(L+2)*(L+3)/6                                    
      END FUNCTION LHGTF

      FUNCTION LSP(L)
        INTEGER, INTENT(IN) :: L
        INTEGER :: LSP
        LSP=L*(L+3)/2   
      END FUNCTION LSP

      FUNCTION LTD(L)
        INTEGER, INTENT(IN) :: L
        INTEGER :: LTD
        LTD=L*(L+1)/2   
      END FUNCTION LTD

      FUNCTION IStride(I,J)
         INTEGER, INTENT(IN) :: I,J
         INTEGER :: IStride
         IStride=LHGTF(J)-LHGTF(I-1)
      END FUNCTION IStride

      FUNCTION LBegin(L)
         INTEGER, INTENT(IN) :: L
         INTEGER :: LBegin
         LBegin=(L*(L+1)*(L+2))/6+1         
      END FUNCTION LBegin

      FUNCTION LEnd(L)
         INTEGER, INTENT(IN) :: L
         INTEGER :: LEnd
         LEnd=LBegin(L)+L*(L+1)/2+L
      END FUNCTION LEnd

      FUNCTION LMNDex(L,M,N)
         INTEGER, INTENT(IN) :: L,M,N
         INTEGER :: LMNDex
         LMNDex=LBegin(L+M+N)+N*(2*(L+M+N)-N+3)/2+M
      END FUNCTION LMNDex

      FUNCTION MatIndx(I,J,N)
         INTEGER, INTENT(IN) :: I,J,N
         INTEGER :: MatIndx
         MatIndx=I+(J-1)*N
      END FUNCTION MatIndx

      FUNCTION CFBlokDex(BS,Istate,Nstate)
        TYPE(BSET)    :: BS
        INTEGER       :: I,CFBlokDex,Istate,Nstate
        CFBlokDex = 0
        DO I=1,Istate-1
           CFBlokDex = CFBlokDex+(BS%LStop%I(I,Nstate)-BS%LStrt%I(I,Nstate)+1)
        ENDDO
      END FUNCTION CFBlokDex

      SUBROUTINE BSetIndx(BS)
         TYPE(BSET) :: BS
         INTEGER    :: K,L,M,N,NC,LMN,MinL,MaxL
         BS%LMNLen=(BS%NASym+1)*(BS%NASym+2)*(BS%NASym+3)/6
         DO K=1,BS%NKind
            BS%BFKnd%I(K)=0
            DO NC=1,BS%NCFnc%I(K)
               MinL=BS%ASymm%I(1,NC,K)
               MaxL=BS%ASymm%I(2,NC,K)      
               BS%LStrt%I(NC,K)=LBegin(MinL)
               BS%LStop%I(NC,K)=LEnd(MaxL)
               BS%BFKnd%I(K)=BS%BFKnd%I(K) &
                            +BS%LStop%I(NC,K)-BS%LStrt%I(NC,K)+1
            ENDDO
         ENDDO
         DO L=0,BS%NASym
            DO M=0,BS%NASym-L
               DO N=0,BS%NASym-L-M
                  LMN=LMNDex(L,M,N)
                  BS%LxDex%I(LMN)=L
                  BS%LyDex%I(LMN)=M
                  BS%LzDex%I(LMN)=N
               ENDDO
            ENDDO
         ENDDO
    END SUBROUTINE BSetIndx
!
END MODULE
