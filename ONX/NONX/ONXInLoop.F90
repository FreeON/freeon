MODULE ONXInLoop
!H=================================================================================
!H MODULE ONXInLoop
!H This MODULE contains:
!H  PUBLIC:
!H  o SUB Scatter
!H  o SUB Digest
!H  o SUB Contract
!H
!H  PRIVATE:
!H
!H Comments:
!H
!H---------------------------------------------------------------------------------
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE ONXParameters
  USE ONXMemory
  USE ONXGet , ONLY: GetAdrB
  USE ONXPut , ONLY: PutSubBlk
#ifdef PARALLEL
  USE FastMatrices
#endif
  !
  IMPLICIT NONE
  PRIVATE
  !
!--------------------------------------------------------------------------------- 
! PUBLIC DECLARATIONS
!--------------------------------------------------------------------------------- 
  PUBLIC :: Scatter
  PUBLIC :: Digest
  PUBLIC :: Contract
  !
CONTAINS
  !
  SUBROUTINE Scatter(N,NA,NB,IndexA,SB,SubInd,DB,KB,K)
!H--------------------------------------------------------------------------------- 
!H SUBROUTINE Scatter(N,NA,NB,IndexA,SB,SubInd,DB,KB,K)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------------------
#ifdef PARALLEL
    TYPE(FASTMAT ), POINTER       :: K
#else
    TYPE(BCSR    ), INTENT(INOUT) :: K
#endif
    TYPE(DSL     ), INTENT(IN   ) :: SB
    TYPE(INT_RNK2), INTENT(IN   ) :: SubInd
    TYPE(DBuf    ), INTENT(IN   ) :: DB
    INTEGER       , INTENT(IN   ) :: N,NA,NB,IndexA
    REAL(DOUBLE)  , INTENT(IN   ) :: KB(N,NA,NB)
    !-------------------------------------------------------------------------------
#ifdef PARALLEL
    TYPE(FASTMAT ), POINTER       :: C
    TYPE(SRST    ), POINTER       :: P
#endif
    INTEGER                       :: I,I0,Ind,IndexB,Ioff
    INTEGER                       :: AtA,NBFA,RS
    INTEGER                       :: AtB,NBFB,CS
    !-------------------------------------------------------------------------------
    !
    AtA  = SubInd%I(1,IndexA)
    NBFA = SubInd%I(2,IndexA)
    RS   = SubInd%I(3,IndexA)
#ifdef PARALLEL
    C => FindFastMatRow_1(K,AtA)    ! Should be improved.
#endif
    !
    DO I = 1,N
       I0     = SB%SLDis%I(I)-3             !some magic number!
       IndexB = INT(ABS(DB%DisBuf%D(I0)))
       AtB    = SubInd%I(1,IndexB)
       NBFB   = SubInd%I(2,IndexB)
       CS     = SubInd%I(3,IndexB)
       !
#ifdef PARALLEL
       SRSTCount = C%Nodes
       P => InsertSRSTNode(C%RowRoot,AtB)
       IF(.NOT.ASSOCIATED(P%MTrix)) THEN
          ALLOCATE(P%MTrix(NBFA,NBFB),STAT=MemStatus)
          CALL IncMem(MemStatus,0,NBFA*NBFB,'AddFASTMATBlok')
          P%MTrix = Zero
       ENDIF
       C%Nodes = SRSTCount
       CALL PutSubBlk(I,N,NBFA,NBFB,NA,NB,RS,CS,P%MTrix,KB)
#else
       CALL GetAdrB(AtA,AtB,Ind,K,0)
       Ioff = K%BlkPt%I(Ind)
       IF(Ind>0) CALL PutSubBlk(I,N,NBFA,NBFB,NA,NB,RS,CS, &
            &                   K%MTrix%D(Ioff),KB)
#endif
    ENDDO
    !
  END SUBROUTINE Scatter
  !
  !
  SUBROUTINE Digest(N,NA,NB,NC,ND,L1,L2,L3,L4,IntSwitch,K,W,D)
!H--------------------------------------------------------------------------------- 
!H SUBROUTINE Digest(N,NA,NB,NC,ND,L1,L2,L3,L4,IntSwitch,K,W,D)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------------------
    INTEGER,INTENT(IN)        :: N,IntSwitch
    INTEGER,INTENT(IN)        :: NA,NB,NC,ND
    INTEGER,INTENT(IN)        :: L1,L2,L3,L4
    REAL(DOUBLE),INTENT(OUT)  :: K(N,NA,NB)
    REAL(DOUBLE),INTENT(IN)   :: W(N,L2,L1,L4,L3)
    REAL(DOUBLE),INTENT(IN)   :: D(NC,ND)
    !-------------------------------------------------------------------------------
    INTEGER                   :: I,IA,IB,IC,ID
    REAL(DOUBLE)              :: Dcd
    !-------------------------------------------------------------------------------
    K=0.0D0

    SELECT CASE (IntSwitch)
    CASE (11)
       DO IB=1,L3
          DO ID=1,L4
             DO IA=1,L1
                DO IC=1,L2
                   Dcd=D(IC,ID)
                   DO I=1,N
                      K(I,IA,IB)=K(I,IA,IB)-W(I,IC,IA,ID,IB)*Dcd
                   END DO
                END DO
             END DO
          END DO
       END DO
    CASE (01)
       DO IB=1,L3
          DO ID=1,L4
             DO IC=1,L1
                Dcd=D(IC,ID)
                DO IA=1,L2
                   DO I=1,N
                      K(I,IA,IB)=K(I,IA,IB)-W(I,IA,IC,ID,IB)*Dcd
                   END DO
                END DO
             END DO
          END DO
       END DO
    CASE (10)
       DO ID=1,L3
          DO IB=1,L4
             DO IA=1,L1
                DO IC=1,L2
                   Dcd=D(IC,ID)
                   DO I=1,N
                      K(I,IA,IB)=K(I,IA,IB)-W(I,IC,IA,IB,ID)*Dcd
                   END DO
                END DO
             END DO
          END DO
       END DO
    CASE (00)
       DO ID=1,L3
          DO IB=1,L4
             DO IC=1,L1
                Dcd=D(IC,ID)
                DO IA=1,L2
                   DO I=1,N
                      K(I,IA,IB)=K(I,IA,IB)-W(I,IA,IC,IB,ID)*Dcd
                   END DO
                END DO
             END DO
          END DO
       END DO
    CASE DEFAULT
       WRITE(*,*) "IntSwitch=",IntSwitch
       CALL Halt(' Illegal IntSwitch in ONX:Digest')
    END SELECT
    !
  END SUBROUTINE Digest
  !
  !
  SUBROUTINE Contract(N,KBra,KKet,NVRR,LngDrv,CDrv,CB,CK,C,U)
!H--------------------------------------------------------------------------------- 
!H SUBROUTINE Contract(N,KBra,KKet,NVRR,LngDrv,CDrv,CB,CK,C,U)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------------------
    INTEGER, INTENT(IN)        :: N,KBra,KKet,NVRR,LngDrv
    INTEGER, INTENT(IN)        :: CDrv(4,LngDrv)
    REAL(DOUBLE), INTENT(IN)   :: CB(KBra,3)
    REAL(DOUBLE), INTENT(IN)   :: CK(N,KKet,3)
    REAL(DOUBLE), INTENT(OUT)  :: C(N,NVRR)
    REAL(DOUBLE), INTENT(IN)   :: U(N,KBra,KKet,NVRR)
    !-------------------------------------------------------------------------------
    INTEGER                    :: I,J,K,L,iC,iP,iQ,iU
    !-------------------------------------------------------------------------------
    IF (KBra==1.AND.KKet==1) THEN
       DO L=1,LngDrv
          iC  = CDrv(1,L)
          iP  = CDrv(2,L)
          iQ  = CDrv(3,L)
          iU  = CDrv(4,L)
          IF (iP==0.AND.iQ==0) THEN
             DO I=1,N
                C(I,iC) = U(I,1,1,iU)
             END DO
          ELSEIF(iP.GT.0.AND.iQ.EQ.0) THEN
             DO I=1,N
                C(I,iC) = U(I,1,1,iU)*CB(1,iP)
             END DO
          ELSEIF(iP.EQ.0.AND.iQ.GT.0) THEN
             DO I=1,N
                C(I,iC) = U(I,1,1,iU)*CK(I,1,iQ)
             END DO
          ELSE
             DO I=1,N
                C(I,iC) = U(I,1,1,iU)*CB(1,iP)*CK(I,1,iQ)
             END DO
          END IF
       END DO ! L
    ELSE
       DO L=1,LngDrv
          iC  = CDrv(1,L)
          iP  = CDrv(2,L)
          iQ  = CDrv(3,L)
          iU  = CDrv(4,L)
    
          DO I=1,N
             C(I,iC)=0.0D0
          END DO
          IF(iP.EQ.0.AND.iQ.EQ.0) THEN
             DO K=1,KKet
                DO J=1,KBra
                   DO I=1,N
                      C(I,iC) = C(I,iC) + U(I,J,K,iU)
                   END DO
                END DO
             END DO
          ELSEIF(iP.GT.0.AND.iQ.EQ.0) THEN
             DO K=1,KKet
                DO J=1,KBra
                   DO I=1,N
                      C(I,iC) = C(I,iC) + U(I,J,K,iU)*CB(J,iP)
                   END DO
                END DO
             END DO
          ELSEIF(iP.EQ.0.AND.iQ.GT.0) THEN
             DO K=1,KKet
                DO J=1,KBra
                   DO I=1,N
                      C(I,iC) = C(I,iC) + U(I,J,K,iU)*CK(I,K,iQ)
                   END DO
                END DO
             END DO
          ELSE
             DO K=1,KKet
                DO J=1,KBra
                   DO I=1,N
                      C(I,iC) = C(I,iC) + U(I,J,K,iU)*CB(J,iP)*CK(I,K,iQ)
                   END DO
                END DO
             END DO
          END IF
       END DO
    END IF
    !
  END SUBROUTINE Contract
  !
END MODULE ONXInLoop
