MODULE MemMan
   USE DerivedTypes
   USE GlobalScalars   
   USE GlobalCharacters
   USE GlobalObjects
   USE ProcessControl    
   IMPLICIT NONE
   INTERFACE New
      MODULE PROCEDURE New_INT_VECT, New_INT_RNK2, &
                       New_INT_RNK3, New_INT_RNK4, &
                       New_DBL_VECT, New_DBL_RNK2, &
                       New_DBL_RNK3, New_DBL_RNK4, &
                       New_DBL_RNK6, &
                       New_CHR_VECT, New_CRDS,     &
#ifdef PARALLEL 
                       New_DBCSR,    New_MPI_INDX, &
#endif    
                       New_INTC,     New_BMATR, &
                       New_BCSR,     New_BSET,     &
                       New_ARGMT,    New_HGRho,    &
                       New_DBuf,     New_IBuf,     &
                       New_IDrv,     New_DSL,      &
                       New_GradD,    New_CMPoles
   END INTERFACE
   INTERFACE Delete
      MODULE PROCEDURE Delete_INT_VECT, Delete_INT_RNK2, &
                       Delete_INT_RNK3, Delete_INT_RNK4, &
                       Delete_DBL_VECT, Delete_DBL_RNK2, &
                       Delete_DBL_RNK3, Delete_DBL_RNK4, &
                       Delete_DBL_RNK6, &
                       Delete_CHR_VECT, Delete_CRDS,     &
#ifdef PARALLEL 
                       Delete_DBCSR,    Delete_MPI_INDX, &
#endif    
                       Delete_INTC,     Delete_BMATR, &
                       Delete_BCSR,     Delete_BSET,     &
                       Delete_ARGMT,    Delete_HGRho,    &
                       Delete_DBuf,     Delete_IBuf,     &
                       Delete_IDrv,     Delete_DSL,      &
                       Delete_GradD,    Delete_CMPoles
   END INTERFACE
!
   INTERFACE SetToBig
      MODULE PROCEDURE SetBig_INT_VECT, SetBig_DBL_VECT
   END INTERFACE
!
!-------------------------------------------------  
!  Allocation keys
!
   INTEGER                :: MemStatus
   INTEGER, PARAMETER     :: ALLOCATED_TRUE =STATUS_TRUE
   INTEGER, PARAMETER     :: ALLOCATED_SKIP_POINTERS_ON =243241
   INTEGER, PARAMETER     :: ALLOCATED_SKIP_POINTERS_OFF=956521
   INTEGER, PARAMETER     :: ALLOCATED_TEMP =19347
   INTEGER, PARAMETER     :: ALLOCATED_FALSE=STATUS_FALSE
!------------------------------------------------------------------------

   CONTAINS 

      SUBROUTINE SetBig_INT_VECT(A)
         TYPE(INT_VECT) :: A
         A%I=BIG_INT
      END SUBROUTINE SetBig_INT_VECT

      SUBROUTINE SetBig_DBL_VECT(A)
         TYPE(DBL_VECT) :: A
         A%D=BIG_DBL
      END SUBROUTINE SetBig_DBL_VECT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE New_INT_VECT(A,N,M_O)
         TYPE(INT_VECT),  INTENT(INOUT) :: A
         INTEGER                        :: Ints
         INTEGER                        :: M,N
         INTEGER,OPTIONAL,INTENT(IN)    :: M_O
         CALL AllocChk(A%Alloc)
         M=1; IF(PRESENT(M_O))M=M_O
         Ints=N-M+1
         ALLOCATE(A%I(M:N),STAT=MemStatus)
         CALL IncMem(MemStatus,Ints,0)
         A%Alloc=ALLOCATED_TRUE
         CALL SetToBig(A)
      END SUBROUTINE New_INT_VECT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE New_INT_RNK2(A,N,M_O)
         TYPE(INT_RNK2),  INTENT(INOUT)  :: A
         INTEGER                         :: Ints
         INTEGER,DIMENSION(2)            :: M,N
         INTEGER,OPTIONAL, &
                 DIMENSION(2),INTENT(IN) :: M_O
         CALL AllocChk(A%Alloc)
         M=1; IF(PRESENT(M_O))M=M_O
         ALLOCATE(A%I(M(1):N(1), M(2):N(2)),STAT=MemStatus)
         Ints=(N(1)-M(1)+1)*(N(2)-M(2)+1)
         CALL IncMem(MemStatus,Ints,0)
         A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE New_INT_RNK2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE New_INT_RNK3(A,N,M_O)
         TYPE(INT_RNK3),  INTENT(INOUT)  :: A
         INTEGER                         :: Ints
         INTEGER,DIMENSION(3)            :: M,N
         INTEGER,OPTIONAL, &
                 DIMENSION(3),INTENT(IN) :: M_O
         CALL AllocChk(A%Alloc)
         M=1; IF(PRESENT(M_O))M=M_O
         ALLOCATE(A%I(M(1):N(1),M(2):N(2),M(3):N(3)),STAT=MemStatus)
         Ints=(N(1)-M(1)+1)*(N(2)-M(2)+1) &
             *(N(3)-M(3)+1)
         CALL IncMem(MemStatus,Ints,0)
         A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE New_INT_RNK3
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE New_INT_RNK4(A,N,M_O)
         TYPE(INT_RNK4),  INTENT(INOUT)  :: A
         INTEGER                         :: Ints
         INTEGER,DIMENSION(4)            :: M,N
         INTEGER,OPTIONAL, &
                 DIMENSION(4),INTENT(IN) :: M_O
         CALL AllocChk(A%Alloc)
         M=1; IF(PRESENT(M_O))M=M_O
         ALLOCATE(A%I(M(1):N(1),M(2):N(2), &
                      M(3):N(3),M(4):N(4)),STAT=MemStatus)
         Ints=(N(1)-M(1)+1)*(N(2)-M(2)+1) &
             *(N(3)-M(3)+1)*(N(4)-M(4)+1)
         CALL IncMem(MemStatus,Ints,0)
         A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE New_INT_RNK4
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE New_DBL_VECT(A,N,M_O)
         TYPE(DBL_VECT),  INTENT(INOUT) :: A
         INTEGER                        :: Dbls
         INTEGER                        :: M,N
         INTEGER,OPTIONAL,INTENT(IN)    :: M_O
         CALL AllocChk(A%Alloc)
         M=1; IF(PRESENT(M_O))M=M_O
         ALLOCATE(A%D(M:N),STAT=MemStatus)
         Dbls=N-M+1
         CALL IncMem(MemStatus,0,Dbls)
         A%Alloc=ALLOCATED_TRUE
         CALL SetToBig(A)
      END SUBROUTINE New_DBL_VECT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE New_DBL_RNK2(A,N,M_O)
         TYPE(DBL_RNK2),  INTENT(INOUT)  :: A
         INTEGER                         :: Dbls
         INTEGER,DIMENSION(2)            :: M,N
         INTEGER,OPTIONAL, &
                 DIMENSION(2),INTENT(IN) :: M_O
         CALL AllocChk(A%Alloc)
         M=1; IF(PRESENT(M_O))M=M_O
         ALLOCATE(A%D(M(1):N(1), M(2):N(2)),STAT=MemStatus)
         Dbls=(N(1)-M(1)+1)*(N(2)-M(2)+1)
         CALL IncMem(MemStatus,0,Dbls)
         A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE New_DBL_RNK2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE New_DBL_RNK3(A,N,M_O)
         TYPE(DBL_RNK3),  INTENT(INOUT)  :: A
         INTEGER                         :: Dbls
         INTEGER,DIMENSION(3)            :: M,N
         INTEGER,OPTIONAL, &
                 DIMENSION(3),INTENT(IN) :: M_O
         CALL AllocChk(A%Alloc)
         M=1; IF(PRESENT(M_O))M=M_O
         ALLOCATE(A%D(M(1):N(1),M(2):N(2),M(3):N(3)),STAT=MemStatus)
         Dbls=(N(1)-M(1)+1)*(N(2)-M(2)+1) &
             *(N(3)-M(3)+1)
         CALL IncMem(MemStatus,0,Dbls)
         A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE New_DBL_RNK3
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE New_DBL_RNK4(A,N,M_O)
         TYPE(DBL_RNK4),  INTENT(INOUT)  :: A
         INTEGER                         :: Dbls
         INTEGER,DIMENSION(4)            :: M,N
         INTEGER,OPTIONAL, &
                 DIMENSION(4),INTENT(IN) :: M_O
         CALL AllocChk(A%Alloc)
         M=1; IF(PRESENT(M_O))M=M_O
         ALLOCATE(A%D(M(1):N(1),M(2):N(2), &
                      M(3):N(3),M(4):N(4)),STAT=MemStatus)
         Dbls=(N(1)-M(1)+1)*(N(2)-M(2)+1) &
             *(N(3)-M(3)+1)*(N(4)-M(4)+1)
         CALL IncMem(MemStatus,0,Dbls)
         A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE New_DBL_RNK4

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE New_DBL_RNK6(A,N,M_O)
         TYPE(DBL_RNK6),  INTENT(INOUT)  :: A
         INTEGER                         :: Dbls
         INTEGER,DIMENSION(6)            :: M,N
         INTEGER,OPTIONAL, &
                 DIMENSION(6),INTENT(IN) :: M_O
         CALL AllocChk(A%Alloc)
         M=1; IF(PRESENT(M_O))M=M_O
         ALLOCATE(A%D(M(1):N(1),M(2):N(2), &
                      M(3):N(3),M(4):N(4), &
                      M(5):N(5),M(6):N(6)),STAT=MemStatus)
         Dbls=(N(1)-M(1)+1)*(N(2)-M(2)+1) &
             *(N(3)-M(3)+1)*(N(4)-M(4)+1) &
             *(N(5)-M(5)+1)*(N(6)-M(6)+1)
         CALL IncMem(MemStatus,0,Dbls)
         A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE New_DBL_RNK6
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE New_CHR_VECT(A,N,M_O)
         TYPE(CHR_VECT),  INTENT(OUT) :: A
         INTEGER,         INTENT(IN)  :: N
         INTEGER,OPTIONAL,INTENT(IN)  :: M_O
         INTEGER                      :: M
         CALL AllocChk(A%Alloc)
         M=1; IF(PRESENT(M_O))M=M_O
         ALLOCATE(A%C(M:N),STAT=MemStatus)
         CALL IncMem(MemStatus,0,0)
         A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE New_CHR_VECT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
      SUBROUTINE New_INTC(A,N,M_O)
         TYPE(INTC),      INTENT(OUT) :: A
         INTEGER,         INTENT(IN)  :: N
         INTEGER,OPTIONAL,INTENT(IN)  :: M_O
         INTEGER                      :: M
         CALL AllocChk(A%Alloc)
         M=1; IF(PRESENT(M_O))M=M_O
         ALLOCATE(A%Def(M:N),STAT=MemStatus)
         ALLOCATE(A%Atoms(M:N,1:4),STAT=MemStatus)
         ALLOCATE(A%Value(M:N),STAT=MemStatus)
         ALLOCATE(A%Constraint(M:N),STAT=MemStatus)
         ALLOCATE(A%Active(M:N),STAT=MemStatus)
         CALL IncMem(MemStatus,0,0)
         A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE New_INTC
!     
!------------------------------------------------------
!     
      SUBROUTINE New_BMATR(A,N,M_O)
         TYPE(BMATR),     INTENT(OUT) :: A
         INTEGER,         INTENT(IN)  :: N
         INTEGER,OPTIONAL,INTENT(IN)  :: M_O
         INTEGER                      :: M
         CALL AllocChk(A%Alloc)
         M=1; IF(PRESENT(M_O))M=M_O
         ALLOCATE(A%IB(M:N,1:12),STAT=MemStatus)
         ALLOCATE(A%B(M:N,1:12),STAT=MemStatus)
         CALL IncMem(MemStatus,0,0)
         A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE New_BMATR
!     
!-----------------------------------------------------
!     
      SUBROUTINE New_CRDS(A)
         TYPE(CRDS),INTENT(INOUT)       :: A
         CALL AllocChk(A%Alloc)
         CALL New(A%BndBox,(/3,2/))
         CALL New(A%AtTyp,A%NAtms)
         CALL New(A%AtNum,A%NAtms)
         CALL New(A%AtNam,A%NAtms)
         CALL New(A%AtMMTyp,A%NAtms)
         CALL New(A%AtMss,A%NAtms)
         CALL New(A%Carts,(/3,A%NAtms/))
         CALL New(A%Vects,(/3,A%NAtms/))
#ifdef PERIODIC
         CALL New(A%BoxCarts,(/3,A%NAtms/))
         CALL New(A%AbBoxCarts,(/3,A%NAtms/))
         CALL New(A%BoxVects,(/3,A%NAtms/))
         CALL New(A%AbCarts,(/3,A%NAtms/))
#endif
         A%Alloc=ALLOCATED_TRUE
         A%ETotal=Zero
      END SUBROUTINE New_CRDS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
      SUBROUTINE New_BCSR(A,N_O,OnAll_O)
         TYPE(BCSR),INTENT(INOUT)             :: A
         INTEGER,OPTIONAL,DIMENSION(3)        :: N_O
         LOGICAL,OPTIONAL                     :: OnAll_O
         LOGICAL                              :: OnAll
         IF(PRESENT(OnAll_O))THEN
            OnAll=OnAll_O
         ELSE
            OnAll=.FALSE.
         ENDIF
         CALL AllocChk(A%Alloc)
#ifdef PARALLEL
         IF(MyId==ROOT.OR.OnAll)THEN
#endif
            IF(PRESENT(N_O))THEN
               A%NAtms=N_O(1)
               A%NBlks=N_O(2) 
               A%NNon0=N_O(3)
            ELSE
               A%NAtms=NAtoms! MaxAtms
               A%NBlks=MaxBlks
               A%NNon0=MaxNon0
            ENDIF
            CALL New(A%RowPt,A%NAtms+1)
            CALL New(A%ColPt,A%NBlks  )
            CALL New(A%BlkPt,A%NBlks  )
            CALL New(A%MTrix,A%NNon0  )
            A%Alloc=ALLOCATED_TRUE
#ifdef PARALLEL
         ENDIF
#endif
      END SUBROUTINE New_BCSR
#ifdef PARALLEL
!----------------------------------------------------------------------------
!     
!
      SUBROUTINE New_DBCSR(A,N_O,NoGlobals_O,Node_O)
         TYPE(DBCSR),INTENT(INOUT)      :: A
         INTEGER,OPTIONAL,DIMENSION(3)  :: N_O
         INTEGER,OPTIONAL               :: Node_O
         LOGICAL,OPTIONAL               :: NoGlobals_O
!------------------------------------------------------
!        Check for a previous allocation
         CALL AllocChk(A%Alloc)
!        Set scalars
         IF(PRESENT(Node_O))THEN
            A%Node=Node_O
         ELSE
            A%Node=MyId
         ENDIF
         A%Alloc=ALLOCATED_TRUE
!        Allocate local variables
         IF(PRESENT(N_O))THEN
            CALL New(A%RowPt,N_O(1)+1)
            CALL New(A%ColPt,N_O(2))
            CALL New(A%BlkPt,N_O(2))
            CALL New(A%MTrix,N_O(3))
         ELSE
            CALL New(A%RowPt,MaxAtmsNode)
            CALL New(A%ColPt,MaxBlksNode)
            CALL New(A%BlkPt,MaxBlksNode)
            CALL New(A%MTrix,MaxNon0Node)
         ENDIF
!        Allocate global variables
         IF(PRESENT(NoGlobals_O))THEN 
            IF(NoGlobals_O)THEN
               RETURN
            ELSE
               CALL Halt(' Logic error in New_DBCSR ')
            ENDIF
         ELSE         
            CALL New(A%GRwPt,NAtoms+1 )
            CALL New(A%GClPt,MaxBlks  )
            A%GRwPt%I(1:NAtoms+1)=BIG_INT
         ENDIF
      END SUBROUTINE New_DBCSR
!----------------------------------------------------------------------------
!
      SUBROUTINE New_MPI_INDX(A)
         TYPE(MPI_INDX),INTENT(INOUT) :: A
         CALL AllocChk(A%Alloc)
         A%Alloc=ALLOCATED_TRUE
         A%Type=0
         CALL New(A%Blks,A%NBlks)
         CALL New(A%Disp,A%NBlks)
      END SUBROUTINE New_MPI_INDX
#endif
!----------------------------------------------------------------------------
!     Allocate a basis set
!
      SUBROUTINE New_BSET(A)
         TYPE(BSET),INTENT(INOUT)       :: A
         CALL AllocChk(A%Alloc)
         CALL New(A%Kinds,A%NKind)
         CALL New(A%NCFnc,A%NKind)
         CALL New(A%BFKnd,A%NAtms)
         CALL New(A%LxDex,A%LMNLen)
         CALL New(A%LyDex,A%LMNLen)
         CALL New(A%LzDex,A%LMNLen)
         CALL New(A%LStrt,(/A%NCtrt,A%NKind/))
         CALL New(A%LStop,(/A%NCtrt,A%NKind/))
         CALL New(A%NPFnc,(/A%NCtrt,A%NKind/))
         CALL New(A%ASymm,(/2,A%NCtrt,A%NKind/))
         CALL New(A%Expnt,(/A%NPrim,A%NCtrt,A%NKind/))
         CALL New(A%CCoef,(/A%LMNLen,A%NPrim,A%NCtrt,A%NKind/))
         A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE New_BSET
!----------------------------------------------------------------------------
!     Allocate ONX distribution buffers
!
      SUBROUTINE New_DBuf(A)
         TYPE(DBuf),INTENT(INOUT)       :: A
         CALL AllocChk(A%Alloc)
         CALL New(A%TCode,A%MAXT)
         CALL New(A%CCode,A%MAXK)
         CALL New(A%TCPop,(/A%MAXT,A%MAXK/))
         CALL New(A%DisPtr,(/3,A%NShells,A%MAXT,A%MAXK/))
         CALL New(A%DisBuf,A%MAXDis)
         CALL New(A%PrmBuf,A%MAXPrm)
         CALL New(A%TBufC,(/A%MAXC,A%MAXD/))
         CALL New(A%TBufP,(/A%MAXP,A%NPrim*A%NPrim+A%MInfo,A%MAXD/))
         A%DisPtr%I=0
         A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE New_DBuf
!----------------------------------------------------------------------------
!     Allocate ONX integral buffers
!
      SUBROUTINE New_IBuf(A)
         TYPE(IBuf),INTENT(INOUT)       :: A
         CALL AllocChk(A%Alloc)
         CALL New(A%W1,A%MAXI)
         CALL New(A%W2,A%MAXI)
         CALL New(A%CB,(/A%NPrim*A%NPrim,3/))
         CALL New(A%CK,(/A%MaxInts,A%NPrim*A%NPrim,3/))
         CALL New(A%WR,(/A%MaxInts*A%NPrim**4,12/))
         CALL New(A%WZ,(/A%MaxInts*A%NPrim**4,5/))
         A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE New_IBuf
!----------------------------------------------------------------------------
!     Allocate ONX integral drivers
!
      SUBROUTINE New_IDrv(A)
         TYPE(IDrv),INTENT(INOUT)       :: A
         INTEGER                        :: LngTmp
         CALL AllocChk(A%Alloc)
         LngTmp=A%LngLoc/3
         CALL New(A%VLOC,A%LngVRR)
         CALL New(A%CDrv,A%LngCC)
         CALL New(A%SLOC,(/3,LngTmp/))
         A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE New_IDrv
!----------------------------------------------------------------------------
!     Allocate ONX ML buffer space
!
      SUBROUTINE New_DSL(A)
         TYPE(DSL),INTENT(INOUT)       :: A
         CALL AllocChk(A%Alloc)
         CALL New(A%SLDis,A%MAXSL)
         CALL New(A%SLPrm,A%MAXSL)
         CALL New(A%SLKey,A%MAXSL)
         A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE New_DSL
!----------------------------------------------------------------------------
!     Allocate ONX gradient driver space
!
      SUBROUTINE New_GradD(A)
        TYPE(GradD),INTENT(INOUT)       :: A
        CALL AllocChk(A%Alloc)
        CALL New(A%GDrv1,(/4,2250/))
        CALL New(A%GDrv2,(/4,10/))
        CALL New(A%GDrv3,(/4,10/))
        CALL New(A%GDrv4,(/4,2500/))
        CALL New(A%GDrv5,(/6,10/))
        A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE New_GradD
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE New_ARGMT(A,N_O)
         TYPE(ARGMT)                    :: A
         INTEGER,OPTIONAL,DIMENSION(2)  :: N_O
         CALL AllocChk(A%Alloc)
         IF(PRESENT(N_O))THEN
            A%NC=N_O(1)
            A%NI=N_O(2)
         ENDIF
         CALL New(A%C,A%NC)
         CALL New(A%I,A%NI)
         A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE New_ARGMT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE Delete_INT_VECT(A)
         TYPE(INT_VECT),INTENT(INOUT) :: A
         INTEGER                      :: Ints
         Ints=SIZE(A%I)
         DEALLOCATE(A%I,STAT=MemStatus)
         CALL DecMem(MemStatus,Ints,0)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_INT_VECT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE Delete_INT_RNK2(A)
         TYPE(INT_RNK2),INTENT(INOUT) :: A
         INTEGER                      :: Ints
         Ints=SIZE(A%I)
         DEALLOCATE(A%I,STAT=MemStatus)
         CALL DecMem(MemStatus,Ints,0)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_INT_RNK2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE Delete_INT_RNK3(A)
         TYPE(INT_RNK3),INTENT(INOUT) :: A
         INTEGER                      :: Ints
         Ints=SIZE(A%I)
         DEALLOCATE(A%I,STAT=MemStatus)
         CALL DecMem(MemStatus,Ints,0)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_INT_RNK3
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE Delete_INT_RNK4(A)
         TYPE(INT_RNK4),INTENT(INOUT) :: A
         INTEGER                      :: Ints
         Ints=SIZE(A%I)
         DEALLOCATE(A%I,STAT=MemStatus)
         CALL DecMem(MemStatus,Ints,0)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_INT_RNK4
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE Delete_DBL_VECT(A)
         TYPE(DBL_VECT),INTENT(INOUT) :: A
         INTEGER                      :: Dbls
         Dbls=SIZE(A%D)
         DEALLOCATE(A%D,STAT=MemStatus)
         CALL DecMem(MemStatus,0,Dbls)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_DBL_VECT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE Delete_DBL_RNK2(A)
         TYPE(DBL_RNK2),INTENT(INOUT) :: A
         INTEGER                      :: Dbls
         Dbls=SIZE(A%D)
         DEALLOCATE(A%D,STAT=MemStatus)
         CALL DecMem(MemStatus,0,Dbls)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_DBL_RNK2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE Delete_DBL_RNK3(A)
         TYPE(DBL_RNK3),INTENT(INOUT) :: A
         INTEGER                      :: Dbls
         Dbls=SIZE(A%D)
         DEALLOCATE(A%D,STAT=MemStatus)
         CALL DecMem(MemStatus,0,Dbls)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_DBL_RNK3
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE Delete_DBL_RNK4(A)
         TYPE(DBL_RNK4),INTENT(INOUT) :: A
         INTEGER                      :: Dbls
         Dbls=SIZE(A%D)
         DEALLOCATE(A%D,STAT=MemStatus)
         CALL DecMem(MemStatus,0,Dbls)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_DBL_RNK4
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE Delete_DBL_RNK6(A)
         TYPE(DBL_RNK6),INTENT(INOUT) :: A
         INTEGER                      :: Dbls
         Dbls=SIZE(A%D)
         DEALLOCATE(A%D,STAT=MemStatus)
         CALL DecMem(MemStatus,0,Dbls)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_DBL_RNK6
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE Delete_CHR_VECT(A)
         TYPE(CHR_VECT) :: A
         INTEGER        :: MemStatus
         DEALLOCATE(A%C,STAT=MemStatus)
         CALL DecMem(MemStatus,0,0)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_CHR_VECT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      SUBROUTINE Delete_INTC(A)
         TYPE(INTC)     :: A
         INTEGER        :: MemStatus
         DEALLOCATE(A%Def,STAT=MemStatus)
         DEALLOCATE(A%Atoms,STAT=MemStatus)
         DEALLOCATE(A%Value,STAT=MemStatus)
         DEALLOCATE(A%Constraint,STAT=MemStatus)
         DEALLOCATE(A%Active,STAT=MemStatus)
         CALL DecMem(MemStatus,0,0)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_INTC
!     
      SUBROUTINE Delete_BMATR(A)
         TYPE(BMATR)     :: A
         INTEGER        :: MemStatus
         DEALLOCATE(A%IB,STAT=MemStatus)
         DEALLOCATE(A%B,STAT=MemStatus)
         CALL DecMem(MemStatus,0,0)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_BMATR
!
!------------------------------------------------------
!
      SUBROUTINE Delete_CRDS(A)
         TYPE(CRDS),INTENT(INOUT)       :: A
         CALL Delete(A%BndBox)
         CALL Delete(A%AtTyp)
         CALL Delete(A%AtNum)
         CALL Delete(A%AtNam)
         CALL Delete(A%AtMMTyp)
         CALL Delete(A%AtMss)
         CALL Delete(A%Carts)
         CALL Delete(A%Vects)
#ifdef PERIODIC
         CALL Delete(A%BoxCarts)
         CALL Delete(A%AbBoxCarts)
         CALL Delete(A%BoxVects)
         CALL Delete(A%AbCarts)
#endif 
         A%NAtms=0
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_CRDS 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE Delete_BCSR(A)
         TYPE(BCSR),INTENT(INOUT)       :: A
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            CALL Delete(A%RowPt)
            CALL Delete(A%ColPt)
            CALL Delete(A%BlkPt)
            CALL Delete(A%MTrix)
            A%NAtms=0
            A%NBlks=0
            A%NNon0=0
            A%Alloc=ALLOCATED_FALSE            
#ifdef PARALLEL
         ENDIF
#endif
      END SUBROUTINE Delete_BCSR
#ifdef PARALLEL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE Delete_DBCSR(A)
         TYPE(DBCSR),INTENT(INOUT) :: A
!        Delete local variables
         IF(AllocQ(A%RowPt%Alloc))CALL Delete(A%RowPt)
         IF(AllocQ(A%ColPt%Alloc))CALL Delete(A%ColPt)
         IF(AllocQ(A%BlkPt%Alloc))CALL Delete(A%BlkPt)
         IF(AllocQ(A%MTrix%Alloc))CALL Delete(A%MTrix)
!        Delete global variables
         IF(AllocQ(A%GRwPt%Alloc))CALL Delete(A%GRwPt)
         IF(AllocQ(A%GClPt%Alloc))CALL Delete(A%GClPt)
         A%NAtms=0
         A%NBlks=0
         A%NNon0=0
         A%Alloc=ALLOCATED_FALSE            
         A%Node=BIG_INT
      END SUBROUTINE Delete_DBCSR

      SUBROUTINE Delete_MPI_INDX(A)
         TYPE(MPI_INDX),INTENT(INOUT) :: A
         INTEGER                      :: IErr
         CALL Delete(A%Blks)
         CALL Delete(A%Disp)
         A%Alloc=ALLOCATED_FALSE            
      END SUBROUTINE Delete_MPI_INDX
#endif
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
      SUBROUTINE Delete_BSET(A)
         TYPE(BSET),INTENT(INOUT)       :: A
         A%NAtms=0; A%NBasF=0; A%NKind=0
         A%NCtrt=0; A%NPrim=0; A%NASym=0
         A%LMNLen=0
         CALL Delete(A%Kinds)
         CALL Delete(A%NCFnc)
         CALL Delete(A%BFKnd)
         CALL Delete(A%LxDex)
         CALL Delete(A%LyDex)
         CALL Delete(A%LzDex)
         CALL Delete(A%LStrt)
         CALL Delete(A%LStop)
         CALL Delete(A%NPFnc)
         CALL Delete(A%ASymm)
         CALL Delete(A%Expnt)
         CALL Delete(A%CCoef)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_BSET
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Delete_DBuf(A)
         TYPE(DBuf),INTENT(INOUT)       :: A
         CALL Delete(A%TCode)
         CALL Delete(A%CCode)
         CALL Delete(A%TCPop)
         CALL Delete(A%DisPtr)
         CALL Delete(A%DisBuf)
         CALL Delete(A%PrmBuf)
         CALL Delete(A%TBufC)
         CALL Delete(A%TBufP)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_DBuf
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Delete_IBuf(A)
         TYPE(IBuf),INTENT(INOUT)       :: A
         CALL Delete(A%W1)
         CALL Delete(A%W2)
         CALL Delete(A%CB)
         CALL Delete(A%CK)
         CALL Delete(A%WR)
         CALL Delete(A%WZ)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_IBuf
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Delete_IDrv(A)
         TYPE(IDrv),INTENT(INOUT)       :: A
         CALL Delete(A%VLOC)
         CALL Delete(A%CDrv)
         CALL Delete(A%SLOC)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_IDrv
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Delete_DSL(A)
         TYPE(DSL),INTENT(INOUT)       :: A
         CALL Delete(A%SLDis)
         CALL Delete(A%SLPrm)
         CALL Delete(A%SLKey)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_DSL
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Delete_GradD(A)
        TYPE(GradD),INTENT(INOUT)       :: A
        CALL Delete(A%GDrv1)
        CALL Delete(A%GDrv2)
        CALL Delete(A%GDrv3)
        CALL Delete(A%GDrv4)
        CALL Delete(A%GDrv5)
        A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_GradD
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
      SUBROUTINE Delete_ARGMT(A)
         TYPE(ARGMT),INTENT(INOUT)      :: A
         A%NI=0; A%NC=0
         CALL Delete(A%I)
         CALL Delete(A%C)
         A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE Delete_ARGMT
!========================================================================================
! Allocate the Multipoles
!========================================================================================
  SUBROUTINE New_CMPoles(A)
    TYPE(CMPoles)                  :: A
!
    CALL AllocChk(A%Alloc)
    CALL New(A%DPole,3)
    CALL New(A%QPole,6)
    A%Alloc=ALLOCATED_TRUE
!
  END SUBROUTINE New_CMPoles
!========================================================================================
! Delete the Multipoles
!========================================================================================
  SUBROUTINE Delete_CMPoles(A)
    TYPE(CMPoles)                  :: A
!
    CALL Delete(A%DPole)
    CALL Delete(A%QPole)
    A%Alloc=ALLOCATED_FALSE
!
  END SUBROUTINE Delete_CMPoles
!========================================================================================
! Allocate the Density
!========================================================================================
  SUBROUTINE New_HGRho(A,N_O)
    TYPE(HGRho)                     :: A
    INTEGER,OPTIONAL,DIMENSION(3)   :: N_O
!
    IF(AllocQ(A%Alloc)) THEN
       IF(PRESENT(N_O)) THEN
          IF(N_O(1) /= A%NExpt) THEN
             A%NExpt = N_O(1)
             CALL Delete(A%NQ)
             CALL Delete(A%Lndx)
             CALL Delete(A%OffQ)
             CALL Delete(A%OffR)
             CALL Delete(A%Expt)
             CALL New(A%NQ  ,A%NExpt)
             CALL New(A%Lndx,A%NExpt)
             CALL New(A%OffQ,A%NExpt)
             CALL New(A%OffR,A%NExpt)
             CALL New(A%Expt,A%NExpt)
          ENDIF
          IF(N_O(2) /= A%NDist) THEN
             A%NDist = N_O(2)
             CALL Delete(A%Qx)
             CALL Delete(A%Qy)
             CALL Delete(A%Qz)
             CALL New(A%Qx   ,A%NDist)
             CALL New(A%Qy   ,A%NDist)
             CALL New(A%Qz   ,A%NDist)
          ENDIF
          IF(N_O(3) /= A%NCoef) THEN
             A%NCoef = N_O(3)
             CALL Delete(A%Co)
             CALL New(A%Co,A%NCoef)
          ENDIF
       ENDIF
    ELSE
       A%Alloc=ALLOCATED_TRUE
       IF(PRESENT(N_O)) THEN
          A%NExpt=N_O(1)
          A%NDist=N_O(2)
          A%NCoef=N_O(3)
       ELSE
          A%NExpt=0
          A%NDist=0
          A%NCoef=0
       ENDIF
       CALL New(A%NQ  ,A%NExpt)
       CALL New(A%Lndx,A%NExpt)
       CALL New(A%OffQ,A%NExpt)
       CALL New(A%OffR,A%NExpt)
       CALL New(A%Expt,A%NExpt)
       CALL New(A%Qx  ,A%NDist)
       CALL New(A%Qy  ,A%NDist)
       CALL New(A%Qz  ,A%NDist)
       CALL New(A%Co  ,A%NCoef)
    ENDIF
!
  END SUBROUTINE New_HGRho
!========================================================================================
! Delete the density
!========================================================================================
  SUBROUTINE Delete_HGRho(A)
    TYPE(HGRho)            :: A
!
    IF(AllocQ(A%Alloc)) THEN
       CALL Delete(A%NQ)
       CALL Delete(A%Lndx)
       CALL Delete(A%OffQ)
       CALL Delete(A%OffR)
       CALL Delete(A%Expt)
       CALL Delete(A%Qx)
       CALL Delete(A%Qy)
       CALL Delete(A%Qz)
       CALL Delete(A%Co)
    ENDIF
!
  END SUBROUTINE Delete_HGRho
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
     FUNCTION AllocQ(Alloc)
         INTEGER :: Alloc
         LOGICAL :: AllocQ
         AllocQ=.FALSE.
         IF(Alloc==ALLOCATED_TRUE)THEN
            AllocQ=.TRUE.
         ENDIF
         RETURN
     END FUNCTION AllocQ
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
     SUBROUTINE AllocChk(Alloc)
         INTEGER,INTENT(IN)              :: Alloc
#ifdef PARALLEL
         CHARACTER(LEN=INTERNAL_INT_LEN) :: ChMyId
#endif
         IF(AllocQ(Alloc))THEN
#ifdef PARALLEL
            WRITE(ChMyId,INTERNAL_INT_FMT)MyId
            ChMyId=ADJUSTL(ChMyId)
            CALL Halt(' On node '//TRIM(ChMyId)  &
                   //', Attempt to allocate memory already allocated.')
#else
            CALL Halt(' Attempt to allocate memory already allocated.')
#endif
         ENDIF
      END SUBROUTINE AllocChk
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
      SUBROUTINE InitMEMS(A_O)
         TYPE(MEMS), INTENT(INOUT), OPTIONAL :: A_O
         IF(PRESENT(A_O))THEN
            A_O%Allocs=0
            A_O%DeAllocs=0
            A_O%MemTab=0
            A_O%MaxMem=0
            A_O%MaxAlloc=0
         ELSE
            MemStats%Allocs=0
            MemStats%DeAllocs=0
            MemStats%MemTab=0
            MemStats%MaxMem=0
            MemStats%MaxAlloc=0
         ENDIF
      END SUBROUTINE InitMEMS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
      SUBROUTINE IncMem(Stats,Ints,Dbls,Proc_O) 
         INTEGER, INTENT(IN)                  :: Stats,Ints,Dbls
         INTEGER                              :: MemInBytes
         CHARACTER(LEN=*), OPTIONAL           :: Proc_O
         CHARACTER(LEN=2*DEFAULT_CHR_LEN)     :: Mssg
         CHARACTER(LEN=INTERNAL_INT_LEN)      :: ChMem,ChTab
         MemInBytes=ToBytes(Ints,Dbls)
         IF(MemStatus/=SUCCEED)THEN
            WRITE(*,*)' Bombed in IncMem, MemStatus = ',Stats
            WRITE(ChMem,INTERNAL_INT_FMT)MemInBytes
            WRITE(ChTab,INTERNAL_DBL_FMT)BToMB(MemStats%MemTab)
            IF(PRESENT(Proc_O))THEN
               Mssg='>>>ALLOCATE error at '//TRIM(Proc_O)//' : '& 
                  //' Attempting to allocate '//ChMem           &
                  //' bytes.'//Rtrn//'    To here, '//ChTab     &
                  //' MB were in use.'
            ELSE
               Mssg='>>>ALLOCATE error: '                       & 
                  //' Attempting to allocate '//ChMem           &
                  //' bytes.'//Rtrn//'    To here, '//ChTab     &
                  //' MB were in use.'
            ENDIF
            CALL Halt(Mssg)
         ENDIF
         MemStats%Allocs=MemStats%Allocs+1         
         MemStats%MemTab=MemStats%MemTab+MemInBytes
         MemStats%MaxMem=MAX(MemStats%MaxMem,MemStats%MemTab)
         MemStats%MaxAlloc=MAX(MemStats%MaxAlloc,MemInBytes)
      END SUBROUTINE IncMem

      FUNCTION BToMB(Bytes)
         INTEGER :: Bytes
         REAL(DOUBLE) :: BToMB
         REAL(DOUBLE),PARAMETER :: Ten24Inv2=1.0D0/(1024.D0)**2
         BToMB=DBLE(Bytes)*Ten24Inv2
      END FUNCTION BToMB
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     
!
      SUBROUTINE DecMem(MemStatus,Ints,Dbls) 
         INTEGER, INTENT(IN)                  :: MemStatus,Ints,Dbls
         INTEGER                              :: MemInBytes
         CHARACTER(LEN=2*DEFAULT_CHR_LEN)     :: Mssg
         CHARACTER(LEN=INTERNAL_INT_LEN)      :: ChMem,ChTab,ChMyId
         MemInBytes=ToBytes(Ints,Dbls)
         IF(MemStatus/=SUCCEED)THEN
            WRITE(ChMem,INTERNAL_INT_FMT)MemInBytes
            WRITE(ChTab,INTERNAL_INT_FMT)MemStats%MemTab
#ifdef PARALLEL
            WRITE(ChMyId,INTERNAL_INT_FMT)MyId
            Mssg='>>>DEALLOCATE error, MyId='//TRIM(ChMyId)//' : '  &
               //' Attempting to deallocate '//TRIM(ChMem)          &
               //' bytes.'//Rtrn//'    To here, '//TRIM(ChTab)      &
               //' bytes were in use.'
#else
            Mssg='>>>DEALLOCATE error:'                        & 
               //' Attempting to deallocate '//TRIM(ChMem)     &
               //' bytes.'//Rtrn//'    To here, '//TRIM(ChTab) &
               //' bytes were in use.'
#endif
            CALL Halt(Mssg)
         ENDIF
         MemStats%MemTab=MemStats%MemTab-MemInBytes
         MemStats%DeAllocs=MemStats%DeAllocs+1         
      END SUBROUTINE DecMem
!----------------------------------------------------------------------------
!
!
      FUNCTION ToBytes(NInt,NDbl)
         INTEGER,INTENT(IN) :: NInt,NDbl
         INTEGER            :: ToBytes
         ToBytes=NInt*BytesPerInt()+NDbl*BytesPerDbl()
      END FUNCTION ToBytes
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
     FUNCTION ToMB(NInt,NDbl)
         INTEGER,INTENT(IN) :: NInt,NDbl
         REAL(DOUBLE)       :: ToMB
         ToMB=NInt*IntToMB+NDbl*DblToMB
      END FUNCTION ToMB
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
      FUNCTION BytesPerInt()
         INTEGER :: BytesPerInt
         IF(KIND(BytesPerInt)==INT4)THEN
            BytesPerInt=4
         ELSEIF(KIND(BytesPerInt)==INT8)THEN
            BytesPerInt=8
         ELSE
            CALL Halt('Unkown byte count in BytesPerInt.')
         ENDIF
      END FUNCTION BytesPerInt
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
      FUNCTION BytesPerDbl()
         INTEGER :: BytesPerDbl
         BytesPerDbl=8
      END FUNCTION BytesPerDbl
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!! ckgan: added NewBCSR and DeleteBCSR because all procs need that.
      SUBROUTINE NewBCSR(A,N_O)
         TYPE(BCSR),INTENT(INOUT)             :: A
         INTEGER,OPTIONAL,DIMENSION(3)        :: N_O
         CALL AllocChk(A%Alloc)
            IF(PRESENT(N_O))THEN
               A%NAtms=N_O(1)
               A%NBlks=N_O(2)
               A%NNon0=N_O(3)
            ELSE
               A%NAtms=NAtoms! MaxAtms
               A%NBlks=MaxBlks
               A%NNon0=MaxNon0
            ENDIF
            CALL New(A%RowPt,A%NAtms+1)
            CALL New(A%ColPt,A%NBlks  )
            CALL New(A%BlkPt,A%NBlks  )
            CALL New(A%MTrix,A%NNon0  )
            A%Alloc=ALLOCATED_TRUE
      END SUBROUTINE NewBCSR
      SUBROUTINE DeleteBCSR(A)
         TYPE(BCSR),INTENT(INOUT)       :: A
            CALL Delete(A%RowPt)
            CALL Delete(A%ColPt)
            CALL Delete(A%BlkPt)
            CALL Delete(A%MTrix)
            A%NAtms=0
            A%NBlks=0
            A%NNon0=0
            A%Alloc=ALLOCATED_FALSE
      END SUBROUTINE DeleteBCSR
END MODULE
