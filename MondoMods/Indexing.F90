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

    SUBROUTINE BlockBuild(G,B,BS,OS)
      TYPE(CRDS)      :: G
      TYPE(BSET)      :: B         
      TYPE(INT_VECT)  :: BS,OS
      INTEGER         :: NA,NK,NC,Stride
      !-------------------------------------------------------------------------------------!
      B%NBasF=0
      ! Off set starts at 1
      OS%I(1)=1      
      DO NA=1,G%NAtms
         ! Block size is total number of basis functions per atom 
         BS%I(NA)=0
         NK=G%AtTyp%I(NA)
         ! Go over contracted functions
         DO NC=1,B%NCFnc%I(NK)
            ! Add in size of each contraction
            Stride=B%LStop%I(NC,NK)-B%LStrt%I(NC,NK)+1
            BS%I(NA)=BS%I(NA)+Stride
            ! Oh yeah, accumulate basis function counter too...
            B%NBasF=B%NBasF+Stride
         ENDDO
         ! Off set counter from block sizes
         IF(NA.GE.2)OS%I(NA)=OS%I(NA-1)+BS%I(NA-1)         
      ENDDO
    END SUBROUTINE BlockBuild

!
END MODULE
