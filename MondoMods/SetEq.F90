MODULE SetXYZ
   USE DerivedTypes
   USE GlobalScalars   
   USE GlobalObjects
   USE MemMan
#ifdef PARALLEL
   USE MondoMPI
#endif
   IMPLICIT NONE
   INTERFACE SetEq
      MODULE PROCEDURE Set_TIME_EQ_TIME,          &
                       Set_BCSR_EQ_BCSR,          &
                       Set_RNK2_EQ_BCSR,          &
                       Set_BCSR_EQ_RNK2,          &
#ifdef PARALLEL
                       Set_DBCSR_EQ_BCSR,         &
                       Set_BCSR_EQ_DBCSR,         &
                       Set_RNK2_EQ_DBCSR,         &
                       Set_DBCSR_EQ_DBCSR,         &
#endif
                       Set_INT_VECT_EQ_INT_SCLR,  &  
                       Set_DBL_VECT_EQ_DBL_SCLR,  &
                       Set_INT_VECT_EQ_INT_VECT,  &
                       Set_DBL_VECT_EQ_DBL_VECT
   END INTERFACE
   EXTERNAL bcsr_to_dens
   CONTAINS
!
   SUBROUTINE SetToI(A,B)
      TYPE(BCSR)              :: A
      TYPE(DBL_RNK2),OPTIONAL :: B
      INTEGER                 :: I,N,N2,Q,R       
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
            A%MTrix%D(R:R+N2-1)=B%D(1:N2,I)
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



!======================================================================
!     Copy a Timer
!======================================================================
      SUBROUTINE Set_TIME_EQ_TIME(T2,T1)
         TYPE(TIME), INTENT(IN)  :: T1        
         TYPE(TIME), INTENT(OUT) :: T2        
         T2%CPUS=T1%CPUS
         T2%Wall=T1%Wall
         T2%FLOP=T1%FLOP
      END SUBROUTINE Set_TIME_EQ_TIME
!======================================================================
!     Copy a BCSR matrix
!======================================================================
      SUBROUTINE Set_BCSR_EQ_BCSR(B,A)
         TYPE(BCSR), INTENT(INOUT) :: A        
         TYPE(BCSR), INTENT(INOUT) :: B        
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif 
            IF(AllocQ(B%Alloc) .AND. &
               (B%NAtms<A%NAtms.OR.B%NBlks<A%NBlks.OR.B%NNon0<A%NNon0) )THEN
               CALL Delete(B)
               CALL New(B)
               B%NAtms=A%NAtms; B%NBlks=A%NBlks; B%NNon0=A%NNon0
            ELSE
               IF(.NOT.AllocQ(B%Alloc))CALL New(B)
               B%NAtms=A%NAtms; B%NBlks=A%NBlks; B%NNon0=A%NNon0
            ENDIF
            B%RowPt%I(1:A%NAtms+1)=A%RowPt%I(1:A%NAtms+1)
            B%ColPt%I(1:A%NBlks)  =A%ColPt%I(1:A%NBlks)
            B%BlkPt%I(1:A%NBlks)  =A%BlkPt%I(1:A%NBlks)
            B%MTrix%D(1:A%NNon0)  =A%MTrix%D(1:A%NNon0)
#ifdef PARALLEL
         ENDIF
#endif
      END SUBROUTINE Set_BCSR_EQ_BCSR
#ifdef PARALLEL
!======================================================================
!     Copy a DBCSR matrix
!======================================================================
      SUBROUTINE Set_DBCSR_EQ_DBCSR(B,A)
         TYPE(DBCSR), INTENT(INOUT) :: A        
         TYPE(DBCSR), INTENT(INOUT) :: B        
         LOGICAL                    :: LimitsQ
         IF(AllocQ(B%Alloc))THEN
            LimitsQ=.NOT.                         &
                  (B%NAtms<=SIZE(A%RowPt%I)).AND. &
                  (B%NBlks<=SIZE(A%ColPt%I)).AND. &
                  (B%NBlks<=SIZE(A%BlkPt%I)).AND. &
                  (B%NNon0<=SIZE(A%MTrix%D))
            IF(LimitsQ)THEN
               CALL Delete(B)
               CALL New(B,(/A%NAtms,A%NBlks,A%NNon0/))
            ENDIF
         ELSE
            CALL New(B,(/A%NAtms,A%NBlks,A%NNon0/))
         ENDIF
!        Local
         B%NAtms=A%NAtms
         B%NBlks=A%NBlks
         B%NNon0=A%NNon0
         B%Node=A%Node        
         CALL SetEq(B%RowPt,A%RowPt,A%NAtms+1)
         CALL SetEq(B%ColPt,A%ColPt,A%NBlks)
         CALL SetEq(B%BlkPt,A%BlkPt,A%NBlks)
         CALL SetEq(B%MTrix,A%MTrix,A%NNon0)
!        Global
         B%GUpDate=A%GUpDate
         CALL SetEq(B%GRwPt,A%GRwPt,NAtoms+1)
         CALL SetEq(B%GClPt,A%GClPt) ! ??
      END SUBROUTINE Set_DBCSR_EQ_DBCSR
#endif
!======================================================================
!     Transform a BCSR matrix into a dense matrix (Rank 2 array) 
!======================================================================
      SUBROUTINE Set_RNK2_EQ_BCSR(B,A)
         TYPE(BCSR),    INTENT(IN)    :: A        
         TYPE(DBL_RNK2),INTENT(INOUT) :: B        
         INTERFACE 
            SUBROUTINE BCSR_TO_DENS(NBasF,NAtoms,NBlks,NNon0,BSiz,OffS, &
                                    A,MTrix,RowPt,ColPt,BlkPt)
               INTEGER                            ,INTENT(IN)  :: NBasF,NAtoms, &
                                                                  NBlks,NNon0
               INTEGER, PARAMETER                              :: DOUBLE=KIND(0.D0)
               REAL(DOUBLE),DIMENSION(NBasF,NBasF),INTENT(OUT) :: A
               REAL(DOUBLE),DIMENSION(NNon0),      INTENT(IN)  :: MTrix
               INTEGER     ,DIMENSION(NAtoms+1),   INTENT(IN)  :: RowPt  
               INTEGER     ,DIMENSION(NBlks),      INTENT(IN)  :: ColPt,BlkPt
               INTEGER     ,DIMENSION(NAtoms),     INTENT(IN)  :: BSiz,OffS
            END SUBROUTINE BCSR_TO_DENS
         END INTERFACE
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            IF(.NOT.AllocQ(B%Alloc))THEN
               CALL New(B,(/NBasF,NBasF/))
            ELSEIF(SIZE(B%D,1)/=NBasF.OR.SIZE(B%D,2)/=NBasF)THEN
               CALL Delete(B)
               CALL New(B,(/NBasF,NBasF/))
            ENDIF
            CALL BCSR_To_DENS(NBasF,NAtoms,A%NBlks,A%NNon0,BSiz%I,OffS%I, &
                              B%D,A%MTrix%D,A%RowPt%I,A%ColPt%I,A%BlkPt%I)
#ifdef PARALLEL
         ENDIF
#endif
      END SUBROUTINE Set_RNK2_EQ_BCSR
!======================================================================
!     Transform a dense matrix (Rank 2 array) into a BCSR matrix 
!======================================================================
      SUBROUTINE Set_BCSR_EQ_RNK2(B,A)
         TYPE(DBL_RNK2), INTENT(INOUT) :: A        
         TYPE(BCSR),     INTENT(INOUT) :: B        
         INTEGER                       :: I,J,P,Q,OI,OJ,MA,NA
         IF(AllocQ(B%Alloc))  &
         CALL Delete(B)
         CALL New(B)
         P=1
         Q=1
         OI=0
         B%RowPt%I(1)=1 
         DO I=1,NAtoms
            OJ=0
            MA=BSiz%I(I) 
            DO J=1,NAtoms
               NA=BSiz%I(J)
               CALL BlockToBlock(NBasF,MA,NA,OI,OJ,A%D,B%MTrix%D(Q:))
               B%BlkPt%I(P)=Q
               B%ColPt%I(P)=J               
               Q=Q+MA*NA
               P=P+1
               B%RowPt%I(I+1)=P
               OJ=OJ+NA
            ENDDO
            OI=OI+MA
         ENDDO
         B%NAtms=NAtoms
         B%NBlks=P-1
         B%NNon0=Q-1
      END SUBROUTINE Set_BCSR_EQ_RNK2

      SUBROUTINE BlockToBlock(NBasF,M,N,OI,OJ,A,B)
         IMPLICIT NONE
         INTEGER,                             INTENT(IN)    :: NBasF,M,N,OI,OJ
         REAL(DOUBLE),DIMENSION(NBasF,NBasF), INTENT(IN)    :: A
         REAL(DOUBLE),DIMENSION(M,N),         INTENT(INOUT) :: B
         INTEGER                                            :: I,J

         DO J=1,N
            DO I=1,M
               B(I,J)=A(OI+I,OJ+J)
            ENDDO
         ENDDO
      END SUBROUTINE BlockToBlock
!
#ifdef PARALLEL
!===================================================================
!     Load the Global Row Pointer from the local Row Pointer
!===================================================================
      SUBROUTINE Load(A,RowPt,Id_O)
         TYPE(DBCSR),     INTENT(IN)    :: A
         TYPE(INT_VECT),  INTENT(INOUT) :: RowPt
         INTEGER,OPTIONAL,INTENT(IN)    :: Id_O
         INTEGER                        :: I,K,Id,I1,I2
!--------------------------------------------------------------------
         Id=A%Node; IF(PRESENT(Id_O))Id=Id_O
         I1=OffSt%I(Id)+1
         I2=I1+A%NAtms
         IF(I2>SIZE(RowPt%I))CALL Halt(' Indexing hosed in Load ')
         K=1
         DO I=I1,I2
            RowPt%I(I)=A%RowPt%I(K)
            K=K+1
         ENDDO
      END SUBROUTINE Load
!===================================================================
!     Zero out the Global Row Pointer
!===================================================================
      SUBROUTINE Clear(A,RowPt,Id_O)
         TYPE(DBCSR),     INTENT(IN)    :: A 
         TYPE(INT_VECT),  INTENT(INOUT) :: RowPt
         INTEGER,OPTIONAL,INTENT(IN)    :: Id_O
         INTEGER                        :: I,Id,I1,I2
!--------------------------------------------------------------------
         Id=A%Node; IF(PRESENT(Id_O))Id=Id_O
         I1=OffSt%I(Id)+1
         I2=I1+A%NAtms
         DO I=I1,I2
            RowPt%I(I)=0
         ENDDO
      END SUBROUTINE Clear
!======================================================================
!     Transform a DBCSR matrix into locally a dense matrices, 
!     one for each node.  Usefull only for printing.  
!======================================================================
      SUBROUTINE Set_RNK2_EQ_DBCSR(B,A,Id_O)
         TYPE(DBCSR),    INTENT(INOUT)       :: A        
         TYPE(DBL_RNK2),INTENT(INOUT)     :: B        
         INTEGER,OPTIONAL                 :: Id_O
         INTEGER                          :: K,L,JP,P,II,JJ,MA,MA1,NA, &
                                             MOff,KOff,IStrt,IStop     
         INTEGER                          :: IL,IG,JG,Id
         IF(B%Alloc/=ALLOCATED_TRUE.AND.B%Alloc/=ALLOCATED_TEMP)THEN
            CALL New(B,(/NBasF,NBasF/))
         ELSEIF(SIZE(B%D,1)/=NBasF.OR.SIZE(B%D,2)/=NBasF)THEN
            CALL Delete(B)
            CALL New(B,(/NBasF,NBasF/))
         ENDIF          
         Id=A%Node; IF(PRESENT(Id_O))Id=Id_O
         B%D(1:NBasF,1:NBasF)=Zero
         DO IL=1,A%NAtms
            IG=OffSt%I(Id)+IL
            MA=BSiz%I(IG)
            II=OffS%I(IG)
            MA1=MA-1
            IStrt=A%RowPt%I(IL)
            IStop=A%RowPt%I(IL+1)-1
            IF(IStrt/=0.AND.IStop/=0)THEN
            DO JP=IStrt,IStop
               JG=A%ColPt%I(JP)
               P =A%BlkPt%I(JP)
               NA=BSiz%I(JG)
               JJ=OffS%I(JG)
               DO K=0,NA-1
                  KOff=JJ+K
                  MOff=P+K*MA
                  DO L=0,MA1
                     B%D(II+L,KOff)=A%MTrix%D(MOff+L)                 
                  ENDDO
               ENDDO
            ENDDO
            ENDIF
         ENDDO
      END SUBROUTINE Set_RNK2_EQ_DBCSR
!======================================================================
!     Gather a distributed BCSR matrix to a serial BCSR matrix on ROOT
!======================================================================
      SUBROUTINE Set_BCSR_EQ_DBCSR(B,A)
         TYPE(DBCSR), INTENT(INOUT) :: A
         TYPE(BCSR),  INTENT(INOUT) :: B
         INTEGER                    :: I,K,NAtms
         TYPE(INT_VECT)             :: MA,MB,MN,NA,NB,NN 
!-----------------------------------------------------------------------
!        Allocate indecies
!
         IF(MyId==ROOT)THEN
            CALL New(MA,NPrc,M_O=0)
            CALL New(MB,NPrc,M_O=0)
            CALL New(MN,NPrc,M_O=0)
            CALL New(NA,NPrc,M_O=0)
            CALL New(NB,NPrc,M_O=0)
            CALL New(NN,NPrc,M_O=0)
         ENDIF
!----------------------------------------------------------
!        Allocate the serial matrix
!
         IF(AllocQ(B%Alloc))CALL Delete(B)
         B%NAtms=NAtoms
         B%NBlks=Reduce(A%NBlks)
         B%NNon0=Reduce(A%NNon0)
!-1
!MINUS     
         CALL New(B,(/B%NAtms,B%NBlks,B%NNon0/))
!------------------------------------
!        Number of atoms per node
!
         NAtms=End%I(MyId)-Beg%I(MyId)+1
         IF(MyID==NPrc-1)NAtms=NAtms+1
!------------------------------------
!        Gather the indecies to root
!
         CALL Gather(NAtms  ,MA)
         CALL Gather(A%NBlks,MB)
         CALL Gather(A%NNon0,MN)
         IF(MyID==ROOT)THEN
!---------------------------------------------------
!           Calculate index offsets (displacements)
!
            NA%I(0)=  NAtms
            NB%I(0)=A%NBlks
            NN%I(0)=A%NNon0
            DO I=1,NPrc-1
               NA%I(I)=MA%I(I)+NA%I(I-1)
               NB%I(I)=MB%I(I)+NB%I(I-1)
               NN%I(I)=MN%I(I)+NN%I(I-1)
            ENDDO
            DO I=NPrc,1,-1
               NA%I(I)=NA%I(I-1)
               NB%I(I)=NB%I(I-1)
               NN%I(I)=NN%I(I-1)
            ENDDO
            NA%I(0)=0
            NB%I(0)=0
            NN%I(0)=0
!---------------------------------------------
!           Copy local portions of the matirx
!
            CALL SetEq(B%RowPt,A%RowPt,  NAtms)
            CALL SetEq(B%ColPt,A%ColPt,A%NBlks)
            CALL SetEq(B%BlkPt,A%BlkPt,A%NBlks)
            CALL SetEq(B%MTrix,A%MTrix,A%NNon0)           
         ENDIF
!----------------------------------
!        Gather the matrix to ROOT
!        
         CALL Gather(A%RowPt,B%RowPt,  NAtms,MA,NA)
         CALL Gather(A%ColPt,B%ColPt,A%NBlks,MB,NB)
         CALL Gather(A%BlkPt,B%BlkPt,A%NBlks,MB,NB)
         CALL Gather(A%MTrix,B%MTrix,A%NNon0,MN,NN)
!--------------------------------------------------
!        Add offsets to achieve correct indexing
!        
         IF(MyId==ROOT)THEN
            DO I=1,NPrc-1
               DO K=NA%I(I)+1,NA%I(I+1)
                  B%RowPt%I(K)=B%RowPt%I(K)+NB%I(I)
               ENDDO
               DO K=NB%I(I)+1,NB%I(I+1)
                  B%BlkPt%I(K)=B%BlkPt%I(K)+NN%I(I)
               ENDDO
            ENDDO
         ENDIF
!------------------
!        Tidy up
!
         IF(MyId==ROOT)THEN
            CALL Delete(MA)
            CALL Delete(MB)
            CALL Delete(MN)
            CALL Delete(NA)
            CALL Delete(NB)
            CALL Delete(NN)
         ENDIF
      END SUBROUTINE Set_BCSR_EQ_DBCSR
!============================================================================
!     Scatter a serial BCSR matrix from ROOT to a distributed BCSR matrix
!============================================================================
      SUBROUTINE Set_DBCSR_EQ_BCSR(B,A)
         TYPE(BCSR),  INTENT(INOUT)  :: A
         TYPE(DBCSR), INTENT(INOUT)  :: B
         INTEGER                     :: I,Id,J,JG,M,MN,MN1,P,  &
                                        NAtms,NBlks,NNon0
         LOGICAL                     :: ReAllocate
!-----------------------------------------------------------------------
!        Allocate if required
!
         IF(.NOT.AllocQ(B%Alloc))CALL New(B)
!------------------------------------------------
!        Distribute to each processor
!
         DO Id=NPrc-1,0,-1
            IF(MyId==ROOT)THEN
               B%NAtms=0
               B%NBlks=1
               B%NNon0=1
               B%RowPt%I(1)=1
               DO I=Beg%I(Id),End%I(Id)
                  M=BSiz%I(I)
                  B%NAtms=B%NAtms+1
                  DO J=A%RowPt%I(I),A%RowPt%I(I+1)-1
                     JG=A%ColPt%I(J)
                     B%ColPt%I(B%NBlks)=JG
                     B%BlkPt%I(B%NBlks)=B%NNon0
                     B%NBlks=B%NBlks+1
                     MN=M*BSiz%I(JG);MN1=MN-1
                     P=A%BlkPt%I(J)         
                     B%MTrix%D(B%NNon0:B%NNon0+MN1)=A%MTrix%D(P:P+MN1) 
                     B%NNon0=B%NNon0+MN
                  ENDDO
                  B%RowPt%I(B%NAtms+1)=B%NBlks
               ENDDO
               B%NBlks=B%NBlks-1
               B%NNon0=B%NNon0-1
! MINUS 
               IF(Id/=ROOT)THEN
                  CALL Send(B%NAtms,Id,1)               
                  CALL Send(B%NBlks,Id,2)               
                  CALL Send(B%NNon0,Id,3)
                  CALL Send(B%RowPt,B%NAtms+1,Id,4)               
                  CALL Send(B%ColPt,B%NBlks,Id,5)               
                  CALL Send(B%BlkPt,B%NBlks,Id,6)               
                  CALL Send(B%MTrix,B%NNon0,Id,7)
               ENDIF               
            ELSEIF(MyId==Id)THEN
               CALL Recv(B%NAtms,ROOT,1)               
               CALL Recv(B%NBlks,ROOT,2)               
               CALL Recv(B%NNon0,ROOT,3)
               ReAllocate=(SIZE(B%RowPt%I)<B%NAtms+1).OR. &
                          (SIZE(B%ColPt%I)<B%NBlks)  .OR. &
                          (SIZE(B%BlkPt%I)<B%NBlks)  .OR. &
                          (SIZE(B%MTrix%D)<B%NNon0)
               IF(ReAllocate)THEN     
                  NAtms=B%NAtms; NBlks=B%NBlks; NNon0=B%NNon0
                  CALL Delete(B)     
                  CALL New(B,(/NAtms,NBlks,NNon0/))
               ENDIF
               CALL Recv(B%RowPt,B%NAtms+1,ROOT,4)               
               CALL Recv(B%ColPt,B%NBlks,ROOT,5)               
               CALL Recv(B%BlkPt,B%NBlks,ROOT,6)               
               CALL Recv(B%MTrix,B%NNon0,ROOT,7)               
            ENDIF
         ENDDO
!-------------------------------------------------------
         CALL BCast(A%NBlks)
         IF(MyId==ROOT)THEN
           CALL SetEq(B%GRwPt,A%RowPt,NAtoms+1)
           CALL SetEq(B%GClPt,A%ColPt,A%NBlks)
         ENDIF
         CALL BCast(B%GRwPt,NAtoms+1)
         CALL BCast(B%GClPt,A%NBlks)
         B%Node=MyId
         B%GUpDate=STATUS_TRUE
!
      END SUBROUTINE Set_DBCSR_EQ_BCSR
!===================================================================
!     Load the Global Row Pointer from the local Row Pointer
!===================================================================
      SUBROUTINE Load_DBCSR_GRwPt(A,Id_O)
         TYPE(DBCSR),     INTENT(INOUT) :: A
         INTEGER,OPTIONAL,INTENT(IN)    :: Id_O          
         INTEGER                        :: I,K,Id,I1,I2
!--------------------------------------------------------------------
         IF(PRESENT(Id_O))THEN
            Id=Id_O
         ELSE
            Id=MyId 
         ENDIF
         K=1
         I1=OffSt%I(Id)+1
         I2=I1+A%NAtms
         IF(I2>NAtoms+1)CALL Halt(' Indexing hosed in Load_DBCSR_GRwPt ')
         DO I=I1,I2
            A%GRwPt%I(I)=A%RowPt%I(K)
            K=K+1
         ENDDO
         CALL AlignNodes()
      END SUBROUTINE Load_DBCSR_GRwPt
!===================================================================
!     Zero out the Global Row Pointer
!===================================================================
      SUBROUTINE Clear_DBCSR_GRwPt(A,Id_O)
         TYPE(DBCSR),     INTENT(INOUT) :: A
         INTEGER,OPTIONAL,INTENT(IN)    :: Id_O          
         INTEGER                        :: I,Id,I1,I2
!--------------------------------------------------------------------
         IF(PRESENT(Id_O))THEN
            Id=Id_O
         ELSE
            Id=MyId 
         ENDIF
         I1=OffSt%I(Id)+1
         I2=I1+A%NAtms
         DO I=I1,I2
            A%GRwPt%I(I)=0
         ENDDO
      END SUBROUTINE Clear_DBCSR_GRwPt
#endif
!===============================================================
!     Copy a double vector to a double vector: Y(1:N)=X(1:N)
!===============================================================
      SUBROUTINE Set_DBL_VECT_EQ_DBL_VECT(Y,X,N_O)
         TYPE(DBL_VECT), INTENT(IN)  :: X        
         TYPE(DBL_VECT), INTENT(INOUT) :: Y        
         INTEGER,OPTIONAL,INTENT(IN) :: N_O
         INTEGER                     :: N         
         INTERFACE 
            SUBROUTINE DBL_VECT_EQ_DBL_VECT(N,Y,X)
               INTEGER, INTENT(IN)                   :: N
               INTEGER, PARAMETER                    :: DOUBLE=KIND(0.D0)
               REAL(DOUBLE),DIMENSION(N),INTENT(IN)  :: X
               REAL(DOUBLE),DIMENSION(N),INTENT(OUT) :: Y
            END SUBROUTINE DBL_VECT_EQ_DBL_VECT
         END INTERFACE            
         N=SIZE(Y%D); IF(PRESENT(N_O))N=N_O
         IF(SIZE(X%D)<N.OR.SIZE(Y%D)<N) & 
            CALL Halt(' Dimensions off in Set_DBL_VECT_DBL_VECT')
         CALL DBL_VECT_EQ_DBL_VECT(N,Y%D,X%D)
      END SUBROUTINE Set_DBL_VECT_EQ_DBL_VECT
!===============================================================
!     Copy an integer vector to an integer vector: Y(1:N)=X(1:N)
!===============================================================
      SUBROUTINE Set_INT_VECT_EQ_INT_VECT(Y,X,N_O)
         TYPE(INT_VECT), INTENT(IN)    :: X        
         TYPE(INT_VECT), INTENT(INOUT) :: Y        
         INTEGER,OPTIONAL,INTENT(IN)   :: N_O
         INTEGER                       :: N         
         INTERFACE 
            SUBROUTINE INT_VECT_EQ_INT_VECT(N,Y,X)
               INTEGER, INTENT(IN)              :: N
               INTEGER,DIMENSION(N),INTENT(IN)  :: X
               INTEGER,DIMENSION(N),INTENT(OUT) :: Y
            END SUBROUTINE INT_VECT_EQ_INT_VECT
         END INTERFACE            
         N=SIZE(Y%I); IF(PRESENT(N_O))N=N_O
         IF(SIZE(X%I)<N.OR.SIZE(Y%I)<N) & 
            CALL Halt(' Dimensions off in Set_INT_VECT_INT_VECT')
         CALL INT_VECT_EQ_INT_VECT(N,Y%I,X%I)
      END SUBROUTINE Set_INT_VECT_EQ_INT_VECT
!===============================================================
!     Set an integer scalar to an integer vector: Y(1:N)=X
!===============================================================
      SUBROUTINE Set_INT_VECT_EQ_INT_SCLR(Y,X,N_O)
         INTEGER,        INTENT(IN)    :: X        
         TYPE(INT_VECT), INTENT(INOUT) :: Y        
         INTEGER,OPTIONAL,INTENT(IN)   :: N_O
         INTEGER                       :: N         
         INTERFACE 
            SUBROUTINE INT_VECT_EQ_INT_SCLR(N,Y,X)
               INTEGER,INTENT(IN)               :: N
               INTEGER,INTENT(IN)               :: X
               INTEGER,DIMENSION(N),INTENT(OUT) :: Y
            END SUBROUTINE INT_VECT_EQ_INT_SCLR
         END INTERFACE            
         N=SIZE(Y%I); IF(PRESENT(N_O))N=N_O
         IF(SIZE(Y%I)<N) & 
            CALL Halt(' Dimensions off in Set_INT_SCLR_INT_VECT')
         CALL INT_VECT_EQ_INT_SCLR(N,Y%I,X)
      END SUBROUTINE Set_INT_VECT_EQ_INT_SCLR
!===============================================================
!     Set a double scalar to a double vector: Y(1:N)=X
!===============================================================
      SUBROUTINE Set_DBL_VECT_EQ_DBL_SCLR(Y,X,N_O)
         REAL(DOUBLE),   INTENT(IN)    :: X        
         TYPE(DBL_VECT), INTENT(INOUT) :: Y        
         INTEGER,OPTIONAL,INTENT(IN)   :: N_O
         INTEGER                       :: N         
         INTERFACE 
            SUBROUTINE DBL_VECT_EQ_DBL_SCLR(N,Y,X)
               INTEGER,INTENT(IN)                    :: N
               INTEGER, PARAMETER                    :: DOUBLE=KIND(0.D0)
               REAL(DOUBLE),INTENT(IN)               :: X
               REAL(DOUBLE),DIMENSION(N),INTENT(OUT) :: Y
            END SUBROUTINE DBL_VECT_EQ_DBL_SCLR
         END INTERFACE            
         N=SIZE(Y%D); IF(PRESENT(N_O))N=N_O
         IF(SIZE(Y%D)<N) & 
            CALL Halt(' Dimensions off in Set_DBL_SCLR_DBL_VECT')
         CALL DBL_VECT_EQ_DBL_SCLR(N,Y%D,X)
      END SUBROUTINE Set_DBL_VECT_EQ_DBL_SCLR
END MODULE
