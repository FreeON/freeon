MODULE BlokOpsDns
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE PrettyPrint
  USE DGEMMS
!-------------------------------------------------
   INTERFACE XPose
      MODULE PROCEDURE XPose_BCSR
   END INTERFACE
!-------------------------------------------------------------------------
   INTEGER,SAVE                             :: BIGBLOK,LWORK
   REAL(DOUBLE), ALLOCATABLE, DIMENSION(:)  :: WORK,VALS
   REAL(DOUBLE), ALLOCATABLE, DIMENSION(:,:):: VECT,MAT1,MAT2
   CONTAINS
!-------------------------------------------------------------------------

   SUBROUTINE PPPrint(M,N,A,S)
      INTEGER :: I,J,K,M,N
      REAL(DOUBLE),DIMENSION(M,N) :: A
      CHARACTER(LEN=*)            :: S
      WRITE(*,*)' --------------------------------------'
      WRITE(*,*)' '
      WRITE(*,*)' ',TRIM(S),' = '
      DO I=1,M
         WRITE(*,33)(A(I,J),J=1,N)
      ENDDO
      WRITE(*,*)' '
   33 FORMAT(20(2x,D12.5))
  END SUBROUTINE PPPrint

   SUBROUTINE PPPrintT(M,N,A,S)
      INTEGER :: I,J,K,M,N
      REAL(DOUBLE),DIMENSION(M,N) :: A
      REAL(DOUBLE),DIMENSION(N,M) :: AT
      CHARACTER(LEN=*)            :: S
      WRITE(*,*)' --------------------------------------'
      WRITE(*,*)' '
      WRITE(*,*)' ',TRIM(S),' = '
      DO I=1,M
         DO J=1,N
            AT(J,I)=A(I,J)
         ENDDO
      ENDDO
      DO I=1,N
         WRITE(*,33)(AT(I,J),J=1,M)
      ENDDO
      WRITE(*,*)' '
   33 FORMAT(20(2x,D12.5))
  END SUBROUTINE PPPrintT

   SUBROUTINE SetDSYEVWork(MaxBlokSize)
      INTEGER, INTENT(IN) :: MaxBlokSize
      BIGBLOK=MaxBlokSize
      WRITE(*,*)' BigBlok = ',BigBlok
      LWORK=MAX(1,3*BIGBLOK)
      ALLOCATE(WORK(LWORK))
      ALLOCATE(VALS(BIGBLOK))
      ALLOCATE(VECT(BIGBLOK,BIGBLOK)) 
      ALLOCATE(MAT1(BIGBLOK,BIGBLOK)); MAT1=Zero
      ALLOCATE(MAT2(BIGBLOK,BIGBLOK)) 
   END SUBROUTINE SetDSYEVWork
!
   SUBROUTINE UnSetDSYEVWork()
      DEALLOCATE(WORK)
      DEALLOCATE(VALS)
      DEALLOCATE(VECT) 
      DEALLOCATE(MAT1) 
      DEALLOCATE(MAT2) 
   END SUBROUTINE UnSetDSYEVWork
!
   SUBROUTINE SetToI(A,B)
      TYPE(BCSR) :: A
      INTEGER    :: I,N,N2,Q,R       
      REAL(DOUBLE),DIMENSION(:,:),OPTIONAL :: B
      Q=1
      R=1
      A%NAtms=0
      A%RowPt%I(1)=1
      IF(PRESENT(B))THEN
         DO I=1,NAtoms
            A%NAtms=A%NAtms+1  
            N2=BSiz%I(I)**2
            A%ColPt%I(Q)=I
            A%BlkPt%I(Q)=R
            A%MTrix%D(R:R+N2-1)=B(1:N2,I)
            Q=Q+1 
            R=R+N2     
            A%RowPt%I(A%NAtms+1)=Q
         ENDDO
      ELSE
         DO I=1,NAtoms
            A%NAtms=A%NAtms+1  
            N=BSiz%I(I)
            A%ColPt%I(Q)=I
            A%BlkPt%I(Q)=R
            CALL DiagI(N,A%MTrix%D(R:))
            Q=Q+1 
            R=R+N**2
            A%RowPt%I(A%NAtms+1)=Q
         ENDDO
      ENDIF
      A%NBlks=Q-1 
      A%NNon0=R-1
   END SUBROUTINE SetToI

   SUBROUTINE DiagI(N,A)
      INTEGER                      :: I,N
      REAL(DOUBLE), DIMENSION(N,N) :: A
      A=Zero
      DO I=1,N
         A(I,I)=One
      ENDDO
   END SUBROUTINE DiagI
!
   SUBROUTINE FunkOnSqMat(N,Funk,A,FOnA)
      INTEGER,                    INTENT(IN)    :: N
      REAL(DOUBLE),DIMENSION(N,N),INTENT(IN)    :: A
      REAL(DOUBLE),DIMENSION(N,N),INTENT(OUT)   :: FOnA
      INTEGER                                   :: I,INFO      
      REAL(DOUBLE), PARAMETER                        :: ZeroTol=1.D-8
      REAL(DOUBLE), EXTERNAL                    :: Funk
      INTERFACE DSYEV
         SUBROUTINE DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
            USE GlobalScalars
            CHARACTER(LEN=1), INTENT(IN)    :: JOBZ, UPLO
            INTEGER,          INTENT(IN)    :: LDA,  LWORK, N
            INTEGER,          INTENT(OUT)   :: INFO
            REAL(DOUBLE),     INTENT(INOUT) :: A(LDA,*)
            REAL(DOUBLE),     INTENT(OUT)   :: W(*)
            REAL(DOUBLE),     INTENT(OUT)   :: WORK(*)
         END SUBROUTINE DSYEV
      END INTERFACE
      VECT(1:N,1:N)=A(1:N,1:N)
      CALL DSYEV('V','U',N,VECT,BIGBLOK,VALS,WORK,LWORK,INFO)
      IF(INFO/=SUCCEED) &
         CALL Halt('DSYEV hosed in FunkOnSqMat. INFO='//TRIM(IntToChar(INFO)))  
      MAT1=Zero
      DO I=1,N
         IF(VALS(I)<1.D-7)THEN
            MAT1(I,I)=Zero
         ELSE
            MAT1(I,I)=Funk(VALS(I))
         ENDIF
      ENDDO     
      MAT2=Zero;
      FOnA=Zero;
      CALL DGEMM__NN(N,N,N,One,VECT(1:N,1:N),MAT1(1:N,1:N),MAT2(1:N,1:N))
      CALL DGEMM__NT(N,N,N,One,MAT2(1:N,1:N),VECT(1:N,1:N),FOnA(1:N,1:N))
   END SUBROUTINE FunkOnSqMat

   FUNCTION Inverse(A) RESULT(InvA)
      REAL(DOUBLE), INTENT(IN)  :: A
      REAL(DOUBLE) :: InvA
      InvA=One/A
   END FUNCTION Inverse

   FUNCTION SqRoot(A) RESULT(SqRtA)
      REAL(DOUBLE), INTENT(IN)  :: A
      REAL(DOUBLE) :: SqRtA
      SqRtA=SQRT(A)
   END FUNCTION SqRoot

   SUBROUTINE XPose_BCSR(A,B)
      TYPE(BCSR)           :: A
      TYPE(BCSR), OPTIONAL :: B
      TYPE(INT_VECT)       :: RowPt,ColPt,BlkPt
      LOGICAL :: SymbolicOnly
      IF(PRESENT(B))THEN
         CALL XPose_GENERIC(A%RowPt%I,A%ColPt%I,A%BlkPt%I, &
                            B%RowPt%I,B%ColPt%I,B%BlkPt%I,A%MTrix%D,B%MTrix%D)
         B%NAtms=A%NAtms
         B%NBlks=A%NBlks
         B%NNon0=A%NNon0
      ELSE
         CALL New(RowPt,NAtoms+1)
         CALL New(ColPt,A%NBlks)
         CALL New(BlkPt,A%NBlks)
         CALL XPose_GENERIC(A%RowPt%I,A%ColPt%I,A%BlkPt%I, &
                            RowPt%I,ColPt%I,BlkPt%I)
         A%RowPt%I(1:NAtoms+1)=RowPt%I(1:NAtoms+1)
         A%ColPt%I(1:A%NBlks)=ColPt%I(1:A%NBlks)      
         A%BlkPt%I(1:A%NBlks)=BlkPt%I(1:A%NBlks)      
         CALL Delete(RowPt)
         CALL Delete(ColPt)
         CALL Delete(BlkPt)
      ENDIF

   END SUBROUTINE XPose_BCSR

   SUBROUTINE XPose_GENERIC(ARowPt,AColPt,ABlkPt,BRowPt,BColPt,BBlkPt,AMTrix,BMTrix)
      REAL(DOUBLE), DIMENSION(:), OPTIONAL :: AMTrix,BMTrix
      INTEGER,      DIMENSION(:)           :: ARowPt,AColPt,ABlkPt, &
                                              BRowPt,BColPt,BBlkPt
      INTEGER                              :: I,J,K,MI,NJ,MN,P,Q,NEXT
      LOGICAL                              :: SymbolicOnly
      
!
      IF(PRESENT(AMTrix).AND.PRESENT(BMTrix))THEN
         SymbolicOnly=.FALSE.
      ELSE
         SymbolicOnly=.TRUE.
      ENDIF
!
      DO I=1,NAtoms
         BRowPt(I)=0
      ENDDO
      DO I=1,NAtoms
         IF((ARowPt(I).NE.0).AND.(ARowPt(I+1)-1.NE.0))THEN
            DO K=ARowPt(I),ARowPt(I+1)-1
               J=AColPt(K)+1
               BRowPt(J)=BRowPt(J)+1
            ENDDO
         ENDIF
      ENDDO
      BRowPt(1)=1
      DO I=1,NAtoms
         BRowPt(I+1)=BRowPt(I)+BRowPt(I+1)
      ENDDO
      Q=1
      IF(SymbolicOnly)THEN
         DO I=1,NAtoms
            MI=BSiz%I(I)
            IF((ARowPt(I).NE.0).AND.(ARowPt(I+1)-1.NE.0))THEN
               DO K=ARowPt(I),ARowPt(I+1)-1
                  J=AColPt(K)
                  P=ABlkPt(K)
                  NJ=BSiz%I(J)
                  NEXT=BRowPt(J)
                  MN=MI*NJ
                  BBlkPt(NEXT)=P
                  BColPt(NEXT)=I
                  BRowPt(J)=NEXT+1
                  Q=Q+MN
               ENDDO
            ENDIF
         ENDDO
      ELSE
         DO I=1,NAtoms
            MI=BSiz%I(I)
            IF((ARowPt(I).NE.0).AND.(ARowPt(I+1)-1.NE.0))THEN
               DO K=ARowPt(I),ARowPt(I+1)-1
                  J=AColPt(K)
                  P=ABlkPt(K)
                  NJ=BSiz%I(J)
                  NEXT=BRowPt(J)
                  MN=MI*NJ
                  CALL XPoseSqMat(MI,NJ,AMTrix(P:),BMTrix(Q:))
                  BBlkPt(NEXT)=Q
                  BColPt(NEXT)=I
                  BRowPt(J)=NEXT+1
                  Q=Q+MN
               ENDDO
            ENDIF
         ENDDO
      ENDIF
      DO I=NAtoms,1,-1
         BRowPt(I+1)=BRowPt(I)
      ENDDO
      BRowPt(1)=1
   END SUBROUTINE XPose_GENERIC

   SUBROUTINE XPoseSqMat(M,N,A,AT)
      INTEGER,                     INTENT(IN)  :: M,N
      REAL(DOUBLE), DIMENSION(M,N),INTENT(IN)  :: A
      REAL(DOUBLE), DIMENSION(N,M),INTENT(OUT) :: AT
      AT=TRANSPOSE(A)
   END SUBROUTINE XPoseSqMat        

END MODULE BlokOpsDns
