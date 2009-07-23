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
                      New_DBL_RNK6, New_CHR_VECT, &
                      New_LOG_VECT, New_CHR10_VECT,&
                      New_BONDDATA, New_ATOMBONDS,&
                      New_PBCInfo,  New_CRDS,     &
                      New_PBCFit ,                &
#ifdef PARALLEL
                      New_DBCSR,    New_MPI_INDX, &
#endif
                      New_INTC,     New_BMATR,    &
                      New_Chol,                   &
                      New_IntCBox,  New_ANGLEDATA,&
                      New_OUTPDATA, New_Sp1x1,    &
                      New_BCSR,     New_BSET,     &
                      New_ARGMT,    New_HGRho,    &
                      New_DBuf,     New_IBuf,     &
                      New_IDrv,     New_DSL,      &
                      New_GradD,    New_CMPoles,  &
                      New_CellSet
  END INTERFACE

  INTERFACE Initialize
     MODULE PROCEDURE Initialize_INT_VECT, Initialize_INT_RNK2,  &
                      Initialize_INT_RNK3, Initialize_INT_RNK4, &
                      Initialize_DBL_VECT, Initialize_DBL_RNK2, &
                      Initialize_DBL_RNK3, Initialize_DBL_RNK4, &
                      Initialize_DBL_RNK6, Initialize_CHR_VECT, &
                      Initialize_LOG_VECT, Initialize_CHR10_VECT,&
                      Initialize_BONDDATA, Initialize_ATOMBONDS,&
                      Initialize_PBCInfo,  Initialize_CRDS,     &
                      Initialize_PBCFit ,                &
#ifdef PARALLEL
                      Initialize_DBCSR,    Initialize_MPI_INDX, &
#endif
                      Initialize_INTC,     Initialize_BMATR,    &
                      Initialize_Chol,                   &
                      Initialize_IntCBox,  Initialize_ANGLEDATA,&
                      Initialize_OUTPDATA, Initialize_Sp1x1,    &
                      Initialize_BCSR,     Initialize_BSET,     &
                      Initialize_ARGMT,    Initialize_HGRho,    &
                      Initialize_DBuf,     Initialize_IBuf,     &
                      Initialize_IDrv,     Initialize_DSL,      &
                      Initialize_GradD,    Initialize_CMPoles,  &
                      Initialize_CellSet
  END INTERFACE

  INTERFACE Delete
     MODULE PROCEDURE Delete_INT_VECT, Delete_INT_RNK2, &
                      Delete_INT_RNK3, Delete_INT_RNK4, &
                      Delete_DBL_VECT, Delete_DBL_RNK2, &
                      Delete_DBL_RNK3, Delete_DBL_RNK4, &
                      Delete_DBL_RNK6, Delete_CHR_VECT, &
                      Delete_LOG_VECT, Delete_CHR10_VECT, &
                      Delete_BONDDATA, Delete_ATOMBONDS,&
                      Delete_PBCInfo,  Delete_CRDS,     &
                      Delete_PBCFit  ,                   &
#ifdef PARALLEL
                      Delete_DBCSR,    Delete_MPI_INDX, &
#endif
                      Delete_INTC,     Delete_BMATR,    &
                      Delete_Chol,                      &
                      Delete_IntCBox,  Delete_ANGLEDATA,&
                      Delete_OUTPDATA, Delete_Sp1x1,    &
                      Delete_BCSR,     Delete_BSET,     &
                      Delete_ARGMT,    Delete_HGRho,    &
                      Delete_DBuf,     Delete_IBuf,     &
                      Delete_IDrv,     Delete_DSL,      &
                      Delete_GradD,    Delete_CMPoles,  &
                      Delete_CellSet
  END INTERFACE

  INTERFACE SetToBig
     MODULE PROCEDURE SetBig_INT_VECT, SetBig_DBL_VECT
  END INTERFACE

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

  SUBROUTINE Initialize_INT_VECT(A)
    TYPE(INT_VECT), INTENT(INOUT) :: A
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_INT_VECT

  SUBROUTINE New_INT_RNK2(A,N,M_O)
    TYPE(INT_RNK2),  INTENT(INOUT)           :: A
    INTEGER                                  :: Ints
    INTEGER,DIMENSION(2)                     :: M,N
    INTEGER,OPTIONAL,DIMENSION(2),INTENT(IN) :: M_O
    CALL AllocChk(A%Alloc)
    M=1; IF(PRESENT(M_O))M=M_O
    ALLOCATE(A%I(M(1):N(1), M(2):N(2)),STAT=MemStatus)
    Ints=(N(1)-M(1)+1)*(N(2)-M(2)+1)
    CALL IncMem(MemStatus,Ints,0)
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_INT_RNK2

  SUBROUTINE Initialize_INT_RNK2(A)
    TYPE(INT_RNK2), INTENT(INOUT) :: A
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_INT_RNK2

  SUBROUTINE New_INT_RNK3(A,N,M_O)
    TYPE(INT_RNK3),  INTENT(INOUT)           :: A
    INTEGER                                  :: Ints
    INTEGER,DIMENSION(3)                     :: M,N
    INTEGER,OPTIONAL,DIMENSION(3),INTENT(IN) :: M_O
    CALL AllocChk(A%Alloc)
    M=1; IF(PRESENT(M_O))M=M_O
    ALLOCATE(A%I(M(1):N(1),M(2):N(2),M(3):N(3)),STAT=MemStatus)
    Ints=(N(1)-M(1)+1)*(N(2)-M(2)+1) &
         *(N(3)-M(3)+1)
    CALL IncMem(MemStatus,Ints,0)
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_INT_RNK3

  SUBROUTINE Initialize_INT_RNK3(A)
    TYPE(INT_RNK3), INTENT(INOUT) :: A
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_INT_RNK3

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

  SUBROUTINE Initialize_INT_RNK4(A)
    TYPE(INT_RNK4), INTENT(INOUT) :: A
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_INT_RNK4

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

  SUBROUTINE Initialize_DBL_VECT(A)
    TYPE(DBL_VECT), INTENT(INOUT) :: A
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_DBL_VECT

  SUBROUTINE New_DBL_RNK2(A,N,M_O)
    TYPE(DBL_RNK2),  INTENT(INOUT)           :: A
    INTEGER                                  :: Dbls
    INTEGER,DIMENSION(2)                     :: M,N
    INTEGER,OPTIONAL,DIMENSION(2),INTENT(IN) :: M_O

    CALL AllocChk(A%Alloc)
    M=1; IF(PRESENT(M_O))M=M_O
    ALLOCATE(A%D(M(1):N(1), M(2):N(2)),STAT=MemStatus)
    Dbls=(N(1)-M(1)+1)*(N(2)-M(2)+1)
    CALL IncMem(MemStatus,0,Dbls)
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_DBL_RNK2

  SUBROUTINE Initialize_DBL_RNK2(A)
    TYPE(DBL_RNK2), INTENT(INOUT) :: A
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_DBL_RNK2

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

  SUBROUTINE Initialize_DBL_RNK3(A)
    TYPE(DBL_RNK3), INTENT(INOUT) :: A
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_DBL_RNK3

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

  SUBROUTINE Initialize_DBL_RNK4(A)
    TYPE(DBL_RNK4), INTENT(INOUT) :: A
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_DBL_RNK4

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

  SUBROUTINE Initialize_DBL_RNK6(A)
    TYPE(DBL_RNK6), INTENT(INOUT) :: A
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_DBL_RNK6

  SUBROUTINE New_CHR10_VECT(A,N,M_O)
    TYPE(CHR10_VECT),  INTENT(OUT) :: A
    INTEGER,         INTENT(IN)  :: N
    INTEGER,OPTIONAL,INTENT(IN)  :: M_O
    INTEGER                      :: M
    CALL AllocChk(A%Alloc)
    M=1; IF(PRESENT(M_O))M=M_O
    ALLOCATE(A%C(M:N),STAT=MemStatus)
    CALL IncMem(MemStatus,0,0)
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_CHR10_VECT

  SUBROUTINE Initialize_CHR10_VECT(A)
    TYPE(CHR10_VECT), INTENT(INOUT) :: A
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_CHR10_VECT

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

  SUBROUTINE Initialize_CHR_VECT(A)
    TYPE(CHR_VECT), INTENT(INOUT) :: A
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_CHR_VECT

  SUBROUTINE New_LOG_VECT(A,N)
    TYPE(LOG_VECT),  INTENT(OUT) :: A
    INTEGER,         INTENT(IN)  :: N
    CALL AllocChk(A%Alloc)
    ALLOCATE(A%L(1:N),STAT=MemStatus)
    A%L=.FALSE.
    CALL IncMem(MemStatus,0,0)
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_LOG_VECT

  SUBROUTINE Initialize_LOG_VECT(A)
    TYPE(LOG_VECT), INTENT(INOUT) :: A
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_LOG_VECT

  !
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  SUBROUTINE New_INTC(A,N)
    TYPE(INTC),      INTENT(OUT) :: A
    INTEGER,         INTENT(IN)  :: N
    !
    IF(AllocQ(A%Alloc)) CALL Delete(A)
    A%N=N
    IF(N==0)RETURN
    !CALL AllocChk(A%Alloc)
    CALL New(A%Def,N)
    A%Def%C(:)(1:10)='XXXXXXXXXX'
    CALL New(A%Atoms,(/N,4/))
    A%Atoms%I=0
    CALL New(A%Cells,(/N,12/))
    A%Cells%I=0
    CALL New(A%Value,N)
    A%Value%D=Zero
    CALL New(A%Constraint,N)
    A%Constraint%L=.FALSE.
    CALL New(A%ConstrValue,N)
    A%ConstrValue%D=Zero
    CALL New(A%Active,N)
    A%Active%L=.TRUE.
    CALL New(A%PredVal,N)
    A%PredVal%D=Zero
    CALL New(A%PredGrad,N)
    A%PredGrad%D=Zero
    CALL New(A%InvHess,N)
    A%InvHess%D=Zero
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_INTC

  SUBROUTINE Initialize_INTC(A)
    TYPE(INTC), INTENT(INOUT) :: A
    CALL Initialize(A%Def)
    CALL Initialize(A%Atoms)
    CALL Initialize(A%Cells)
    CALL Initialize(A%Value)
    CALL Initialize(A%Constraint)
    CALL Initialize(A%ConstrValue)
    CALL Initialize(A%Active)
    CALL Initialize(A%PredVal)
    CALL Initialize(A%PredGrad)
    CALL Initialize(A%InvHess)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_INTC

  SUBROUTINE New_BMATR(A,N)
    TYPE(BMATR) :: A
    INTEGER     :: N
    CALL AllocChk(A%Alloc)
    CALL New(A%IB,(/N,4/))
    CALL New(A%B,(/N,12/))
    CALL New(A%BLI,N)
    CALL New(A%BL,(/N,9/))
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_BMATR

  SUBROUTINE Initialize_BMATR(A)
    TYPE(BMATR), INTENT(INOUT) :: A
    CALL Initialize(A%IB)
    CALL Initialize(A%B)
    CALL Initialize(A%BLI)
    CALL Initialize(A%BL)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_BMATR

  SUBROUTINE New_ANGLEDATA(A,N)
    TYPE(ANGLEDATA) :: A
    INTEGER         :: N
    A%N=N
    CALL New(A%IJK,(/3,N/))
    CALL New(A%Type,N)
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_ANGLEDATA

  SUBROUTINE Initialize_ANGLEDATA(A)
    TYPE(ANGLEDATA), INTENT(INOUT) :: A
    CALL Initialize(A%IJK)
    CALL Initialize(A%Type)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_ANGLEDATA

  SUBROUTINE Delete_ANGLEDATA(A)
    TYPE(ANGLEDATA) :: A
    CALL Delete(A%IJK)
    CALL Delete(A%Type)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_ANGLEDATA

  SUBROUTINE New_OUTPDATA(A,N)
    TYPE(OUTPDATA) :: A
    INTEGER         :: N
    A%N=N
    CALL New(A%IJKL,(/4,N/))
    CALL New(A%Type,N)
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_OUTPDATA

  SUBROUTINE Initialize_OUTPDATA(A)
    TYPE(OUTPDATA), INTENT(INOUT) :: A
    CALL Initialize(A%IJKL)
    CALL Initialize(A%Type)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_OUTPDATA

  SUBROUTINE Delete_OUTPDATA(A)
    TYPE(OUTPDATA) :: A
    CALL Delete(A%IJKL)
    CALL Delete(A%Type)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_OUTPDATA

  SUBROUTINE New_IntCBox(A,NBox,NX,NY,NZ,NatmsLoc)
    INTEGER        :: NBox,NatmsLoc,NX,NY,NZ
    TYPE(IntCBox)  :: A
    A%N=NBox
    A%NX=NX
    A%NY=NY
    A%NZ=NZ
    A%NatmsLoc=NatmsLoc
    CALL New(A%I,NBox+1)
    CALL New(A%J,NatmsLoc)
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_IntCBox

  SUBROUTINE Initialize_IntCBox(A)
    TYPE(IntCBox), INTENT(INOUT) :: A
    CALL Initialize(A%I)
    CALL Initialize(A%J)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_IntCBox

  SUBROUTINE Delete_IntCBox(A)
    TYPE(IntCBox) :: A
    CALL Delete(A%I)
    CALL Delete(A%J)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_IntCBox

  SUBROUTINE New_ATOMBONDS(A,NatmsLoc,MaxBonds)
    TYPE(ATOMBONDS) :: A
    INTEGER         :: NatmsLoc,MaxBonds
    !
    IF(AllocQ(A%Alloc)) CALL Delete(A)
    IF(NatmsLoc==0) RETURN
    A%N1=NatmsLoc
    A%N2=MaxBonds
    CALL New(A%Count,NatmsLoc)
    CALL New(A%Bonds,(/NatmsLoc,MaxBonds/))
    CALL New(A%Atoms,(/NatmsLoc,MaxBonds/))
    IF(NatmsLoc/=0) THEN
       CALL INT_VECT_EQ_INT_SCLR(NatmsLoc,A%Count%I(1),0)
       IF(MaxBonds/=0) THEN
          CALL INT_VECT_EQ_INT_SCLR(NatmsLoc*MaxBonds,A%Bonds%I(1,1),0)
          CALL INT_VECT_EQ_INT_SCLR(NatmsLoc*MaxBonds,A%Atoms%I(1,1),0)
       ENDIF
    ENDIF
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_ATOMBONDS

  SUBROUTINE Initialize_ATOMBONDS(A)
    TYPE(ATOMBONDS), INTENT(INOUT) :: A
    CALL Initialize(A%Count)
    CALL Initialize(A%Bonds)
    CALL Initialize(A%Atoms)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_ATOMBONDS

  SUBROUTINE Delete_ATOMBONDS(A)
    TYPE(ATOMBONDS) :: A
    !
    IF(A%N1==0) RETURN
    A%N1=0
    A%N2=0
    CALL Delete(A%Count)
    CALL Delete(A%Bonds)
    CALL Delete(A%Atoms)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_ATOMBONDS

  SUBROUTINE New_Sp1x1(A,NRow,NZ,Symb_O)
    TYPE(Sp1x1)       :: A
    INTEGER           :: NRow,NZ
    LOGICAL,OPTIONAL  :: Symb_O
    !
    CALL New(A%IA,NRow+1)
    CALL New(A%JA,NZ)
    IF(PRESENT(Symb_O)) THEN
       IF(.NOT.Symb_O) THEN
          CALL New(A%AN,NZ)
       ENDIF
    ELSE
       CALL New(A%AN,NZ)
    ENDIF
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_Sp1x1

  SUBROUTINE Initialize_Sp1x1(A)
    TYPE(Sp1x1), INTENT(INOUT) :: A
    CALL Initialize(A%IA)
    CALL Initialize(A%JA)
    CALL Initialize(A%AN)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_Sp1x1

  SUBROUTINE Delete_Sp1x1(A)
    TYPE(Sp1x1) :: A
    CALL Delete(A%IA)
    CALL Delete(A%JA)
    IF(AllocQ(A%AN%Alloc)) CALL Delete(A%AN)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_Sp1x1

  SUBROUTINE New_BONDDATA(A,NBond)
    TYPE(BONDDATA) :: A
    INTEGER NBond,NatmsLoc
    IF(AllocQ(A%Alloc)) CALL Delete(A)
    A%N=NBond
    IF(NBond==0) RETURN
    CALL New(A%IJ,(/2,NBond/))
    CALL New(A%Length,NBond)
    CALL New(A%Type,NBond)
    CALL New(A%HBExtraSN,NBond)
    CALL New(A%HBExtraNC,NBond)
    CALL New(A%LonelyAtom,NBond)
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_BONDDATA

  SUBROUTINE Initialize_BONDDATA(A)
    TYPE(BONDDATA), INTENT(INOUT) :: A
    CALL Initialize(A%IJ)
    CALL Initialize(A%Length)
    CALL Initialize(A%Type)
    CALL Initialize(A%HBExtraSN)
    CALL Initialize(A%HBExtraNC)
    CALL Initialize(A%LonelyAtom)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_BONDDATA

  SUBROUTINE Delete_BONDDATA(A)
    TYPE(BONDDATA) :: A
    IF(A%N==0) THEN
       A%Alloc=ALLOCATED_FALSE
       RETURN
    ENDIF
    A%N=0
    CALL Delete(A%IJ)
    CALL Delete(A%Length)
    CALL Delete(A%Type)
    CALL Delete(A%HBExtraSN)
    CALL Delete(A%HBExtraNC)
    CALL Delete(A%LonelyAtom)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_BONDDATA

  SUBROUTINE New_Chol(A,NCart,ChNon0)
    INTEGER NCart,ChNon0
    TYPE(Cholesky) :: A
    CALL New(A%GcScale,NCart)
    CALL New(A%Perm,NCart)
    CALL New(A%IPerm,NCart)
    CALL New(A%ChRowPt,NCart+1)
    CALL New(A%ChColPt,ChNon0)
    CALL New(A%ChDiag,NCart)
    CALL New(A%ChFact,ChNon0)
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_Chol

  SUBROUTINE Initialize_Chol(A)
    TYPE(Cholesky), INTENT(INOUT) :: A
    CALL Initialize(A%GcScale)
    CALL Initialize(A%Perm)
    CALL Initialize(A%IPerm)
    CALL Initialize(A%ChRowPt)
    CALL Initialize(A%ChColPt)
    CALL Initialize(A%ChDiag)
    CALL Initialize(A%ChFact)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_Chol

  SUBROUTINE Delete_Chol(A)
    TYPE(Cholesky) :: A
    CALL Delete(A%GcScale)
    CALL Delete(A%Perm)
    CALL Delete(A%IPerm)
    CALL Delete(A%ChRowPt)
    CALL Delete(A%ChColPt)
    CALL Delete(A%ChDiag)
    CALL Delete(A%ChFact)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_Chol

  SUBROUTINE New_PBCInfo(A)
    TYPE(PBCInfo),INTENT(INOUT)       :: A
    CALL AllocChk(A%Alloc)
    CALL New(A%AutoW     ,3)
    CALL New(A%SuperCell ,3)
    CALL New(A%CellCenter,3)
    CALL New(A%TransVec  ,3)
    CALL New(A%BoxShape  ,(/3,3/))
    CALL New(A%InvBoxSh  ,(/3,3/))
    CALL New(A%LatFrc    ,(/3,3/))
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_PBCInfo

  SUBROUTINE Initialize_PBCInfo(A)
    TYPE(PBCInfo), INTENT(INOUT) :: A
    CALL Initialize(A%AutoW)
    CALL Initialize(A%SuperCell)
    CALL Initialize(A%CellCenter)
    CALL Initialize(A%TransVec)
    CALL Initialize(A%BoxShape)
    CALL Initialize(A%InvBoxSh)
    CALL Initialize(A%LatFrc)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_PBCInfo

  SUBROUTINE New_CRDS(A)
    TYPE(CRDS),INTENT(INOUT)       :: A
    CALL AllocChk(A%Alloc)
    CALL New(A%ETotalPerSCF, 256, 0)
    CALL New(A%BndBox,(/3,2/))
    CALL New(A%PBC)
    ! We are not allocating OvCells and InCells here since we don't know
    ! yet what NCells is going to be.
    CALL Initialize(A%OvCells)
    CALL Initialize(A%InCells)
    CALL New(A%AtNum,A%NAtms)
    CALL New(A%AtTyp,A%NAtms)
    CALL New(A%AtNam,A%NAtms)
    CALL New(A%AtMss,A%NAtms)
    CALL New(A%CConstrain,A%NAtms)
    CALL New(A%DoFreq,A%NAtms)
    CALL New(A%Carts,(/3,A%NAtms/))
    CALL New(A%BoxCarts,(/3,A%NAtms/))
    CALL New(A%Velocity,(/3,A%NAtms/))
    CALL New(A%Gradients,(/3,A%NAtms/))
    CALL New(A%Fext,(/3,A%NAtms/))
    CALL New(A%Displ,(/3,A%NAtms/))
    CALL New(A%PBCDispl,(/3,3/))

    A%ETotalPerSCF%D = 0.0D0
    A%ETotal = Zero
    A%Alloc = ALLOCATED_TRUE
  END SUBROUTINE New_CRDS

  SUBROUTINE Initialize_CRDS(A)
    TYPE(CRDS), INTENT(INOUT) :: A
    CALL Initialize(A%ETotalPerSCF)
    CALL Initialize(A%BndBox)
    CALL Initialize(A%PBC)
    CALL Initialize(A%OvCells)
    CALL Initialize(A%InCells)
    CALL Initialize(A%AtNum)
    CALL Initialize(A%AtTyp)
    CALL Initialize(A%AtNam)
    CALL Initialize(A%AtMss)
    CALL Initialize(A%CConstrain)
    CALL Initialize(A%DoFreq)
    CALL Initialize(A%Carts)
    CALL Initialize(A%BoxCarts)
    CALL Initialize(A%Velocity)
    CALL Initialize(A%Gradients)
    CALL Initialize(A%Fext)
    CALL Initialize(A%Displ)
    CALL Initialize(A%PBCDispl)
    A%NAtms = Zero
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_CRDS

  SUBROUTINE New_BCSR(A,N_O,OnAll_O,NSMat_O)
    TYPE(BCSR),INTENT(INOUT)             :: A
    INTEGER,OPTIONAL,DIMENSION(3)        :: N_O
    LOGICAL,OPTIONAL                     :: OnAll_O
    INTEGER,OPTIONAL                     :: NSMat_O
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
       A%NSMat=1
       IF(PRESENT(NSMat_O))A%NSMat=NSMat_O
       IF(PRESENT(N_O))THEN
          A%NAtms=N_O(1)
          A%NBlks=N_O(2)
          A%NNon0=N_O(3)
          CALL New(A%RowPt,A%NAtms+1)
          CALL New(A%ColPt,A%NBlks  )
          CALL New(A%BlkPt,A%NBlks  )!*A%NSMat)
          CALL New(A%MTrix,A%NNon0  )
       ELSE
          A%NAtms=NAtoms! MaxAtms
          A%NBlks=MaxBlks
          A%NNon0=MaxNon0*A%NSMat
          CALL New(A%RowPt,A%NAtms+1)
          CALL New(A%ColPt,A%NBlks  )
          CALL New(A%BlkPt,A%NBlks  )!*A%NSMat)
          CALL New(A%MTrix,A%NNon0  )
       ENDIF
       A%Alloc=ALLOCATED_TRUE
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE New_BCSR

  SUBROUTINE Initialize_BCSR(A)
    TYPE(BCSR), INTENT(INOUT) :: A
    CALL Initialize(A%RowPt)
    CALL Initialize(A%ColPt)
    CALL Initialize(A%BlkPt)
    CALL Initialize(A%MTrix)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_BCSR

#ifdef PARALLEL

  SUBROUTINE New_DBCSR(A,N_O,NoGlobals_O,Node_O,NSMat_O)
    TYPE(DBCSR),INTENT(INOUT)      :: A
    INTEGER,OPTIONAL,DIMENSION(3)  :: N_O
    INTEGER,OPTIONAL               :: Node_O,NSMat_O
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
    A%NSMat=1
    IF(PRESENT(NSMat_O))A%NSMat=NSMat_O
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
       CALL New(A%MTrix,MaxNon0Node*A%NSMat)
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

  SUBROUTINE Initialize_DBCSR(A)
    TYPE(DBCSR), INTENT(INOUT) :: A
    CALL Intialize(A%RowPt)
    CALL Intialize(A%ColPt)
    CALL Intialize(A%BlkPt)
    CALL Intialize(A%MTrix)
    CALL Intialize(A%GRwPt)
    CALL Intialize(A%GClPt)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_DBCSR

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

  SUBROUTINE Initialize_MPI_INDX(A)
    TYPE(MPI_INDX), INTENT(INOUT) :: A
    CALL Initialize(A%Disp)
    CALL Initialize(A%Blks)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_MPI_INDX

#endif

  !----------------------------------------------------------------------------
  !     Allocate a basis set
  !
  SUBROUTINE New_BSET(A)
    TYPE(BSET),INTENT(INOUT)       :: A
    CALL AllocChk(A%Alloc)
    CALL New(A%Kinds,A%NKind)
    CALL New(A%AtNam,A%NKind)
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
    IF(A%HasECPs)THEN
       CALL New(A%NCoreEl,A%NKind)
       CALL New(A%NTyp1PF,A%NKind)
       CALL New(A%ProjEll,A%NKind)
       CALL New(A%Typ1Ell,(/A%Typ1Fnk,A%NKind/))
       CALL New(A%Typ1Exp,(/A%Typ1Fnk,A%NKind/))
       CALL New(A%Typ1CCo,(/A%Typ1Fnk,A%NKind/))
       CALL New(A%NTyp2PF,(/A%MxProjL,A%NKind/),(/0,1/))
       CALL New(A%Typ2Ell,(/A%Typ2Fnk,A%MxProjL,A%NKind/),(/1,0,1/))
       CALL New(A%Typ2Exp,(/A%Typ2Fnk,A%MxProjL,A%NKind/),(/1,0,1/))
       CALL New(A%Typ2CCo,(/A%Typ2Fnk,A%MxProjL,A%NKind/),(/1,0,1/))
    ENDIF
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_BSET

  SUBROUTINE Initialize_BSET(A)
    TYPE(BSET), INTENT(INOUT) :: A
    CALL Initialize(A%Kinds)
    CALL Initialize(A%AtNam)
    CALL Initialize(A%NCFnc)
    CALL Initialize(A%BFKnd)
    CALL Initialize(A%LxDex)
    CALL Initialize(A%LyDex)
    CALL Initialize(A%LzDex)
    CALL Initialize(A%LStrt)
    CALL Initialize(A%LStop)
    CALL Initialize(A%NPFnc)
    CALL Initialize(A%ASymm)
    CALL Initialize(A%Expnt)
    CALL Initialize(A%CCoef)
    CALL Initialize(A%NCoreEl)
    CALL Initialize(A%NTyp1PF)
    CALL Initialize(A%ProjEll)
    CALL Initialize(A%Typ1Ell)
    CALL Initialize(A%Typ1Exp)
    CALL Initialize(A%Typ1CCo)
    CALL Initialize(A%NTyp2PF)
    CALL Initialize(A%Typ2Ell)
    CALL Initialize(A%Typ2Exp)
    CALL Initialize(A%Typ2CCo)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_BSET

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

  SUBROUTINE Initialize_DBuf(A)
    TYPE(DBuf), INTENT(INOUT) :: A
    CALL Initialize(A%TCode)
    CALL Initialize(A%CCode)
    CALL Initialize(A%TCPop)
    CALL Initialize(A%TBufC)
    CALL Initialize(A%TBufP)
    CALL Initialize(A%DisPtr)
    CALL Initialize(A%DisBuf)
    CALL Initialize(A%PrmBuf)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_DBuf

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

  SUBROUTINE Initialize_IBuf(A)
    TYPE(IBuf), INTENT(INOUT) :: A
    CALL Initialize(A%W1)
    CALL Initialize(A%W2)
    CALL Initialize(A%CB)
    CALL Initialize(A%CK)
    CALL Initialize(A%WR)
    CALL Initialize(A%WZ)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_IBuf

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

  SUBROUTINE Initialize_IDrv(A)
    TYPE(IDrv), INTENT(INOUT) :: A
    CALL Initialize(A%VLOC)
    CALL Initialize(A%CDrv)
    CALL Initialize(A%SLOC)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_IDrv

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

  SUBROUTINE Initialize_DSL(A)
    TYPE(DSL), INTENT(INOUT) :: A
    CALL Initialize(A%SLDis)
    CALL Initialize(A%SLPrm)
    CALL Initialize(A%SLKey)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_DSL

  !----------------------------------------------------------------------------
  !     Allocate ONX gradient driver space
  !
  SUBROUTINE New_GradD(A)
    TYPE(GradD),INTENT(INOUT)       :: A
    CALL AllocChk(A%Alloc)
    CALL New(A%GDrv1,(/4,2250/))
    CALL New(A%GDrv2,(/5,10/))
    CALL New(A%GDrv3,(/5,10/))
    CALL New(A%GDrv4,(/4,2500/))
    CALL New(A%GDrv5,(/6,10/))
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_GradD

  SUBROUTINE Initialize_GradD(A)
    TYPE(GradD), INTENT(INOUT) :: A
    CALL Initialize(A%GDrv1)
    CALL Initialize(A%GDrv2)
    CALL Initialize(A%GDrv3)
    CALL Initialize(A%GDrv4)
    CALL Initialize(A%GDrv5)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_GradD

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

  SUBROUTINE Initialize_ARGMT(A)
    TYPE(ARGMT), INTENT(INOUT) :: A
    CALL Initialize(A%I)
    CALL Initialize(A%C)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_ARGMT

  SUBROUTINE Delete_INT_VECT(A)
    TYPE(INT_VECT),INTENT(INOUT) :: A
    INTEGER                      :: Ints
    Ints=SIZE(A%I)
    DEALLOCATE(A%I,STAT=MemStatus)
    CALL DecMem(MemStatus,Ints,0)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_INT_VECT

  SUBROUTINE Delete_INT_RNK2(A)
    TYPE(INT_RNK2),INTENT(INOUT) :: A
    INTEGER                      :: Ints
    Ints=SIZE(A%I)
    DEALLOCATE(A%I,STAT=MemStatus)
    CALL DecMem(MemStatus,Ints,0)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_INT_RNK2

  SUBROUTINE Delete_INT_RNK3(A)
    TYPE(INT_RNK3),INTENT(INOUT) :: A
    INTEGER                      :: Ints
    Ints=SIZE(A%I)
    DEALLOCATE(A%I,STAT=MemStatus)
    CALL DecMem(MemStatus,Ints,0)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_INT_RNK3

  SUBROUTINE Delete_INT_RNK4(A)
    TYPE(INT_RNK4),INTENT(INOUT) :: A
    INTEGER                      :: Ints
    Ints=SIZE(A%I)
    DEALLOCATE(A%I,STAT=MemStatus)
    CALL DecMem(MemStatus,Ints,0)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_INT_RNK4

  SUBROUTINE Delete_DBL_VECT(A)
    TYPE(DBL_VECT),INTENT(INOUT) :: A
    INTEGER                      :: Dbls
    Dbls=SIZE(A%D)
    DEALLOCATE(A%D,STAT=MemStatus)
    CALL DecMem(MemStatus,0,Dbls)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_DBL_VECT

  SUBROUTINE Delete_DBL_RNK2(A)
    TYPE(DBL_RNK2),INTENT(INOUT) :: A
    INTEGER                      :: Dbls
    Dbls=SIZE(A%D)
    DEALLOCATE(A%D,STAT=MemStatus)
    CALL DecMem(MemStatus,0,Dbls)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_DBL_RNK2

  SUBROUTINE Delete_DBL_RNK3(A)
    TYPE(DBL_RNK3),INTENT(INOUT) :: A
    INTEGER                      :: Dbls
    Dbls=SIZE(A%D)
    DEALLOCATE(A%D,STAT=MemStatus)
    CALL DecMem(MemStatus,0,Dbls)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_DBL_RNK3

  SUBROUTINE Delete_DBL_RNK4(A)
    TYPE(DBL_RNK4),INTENT(INOUT) :: A
    INTEGER                      :: Dbls
    Dbls=SIZE(A%D)
    DEALLOCATE(A%D,STAT=MemStatus)
    CALL DecMem(MemStatus,0,Dbls)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_DBL_RNK4

  SUBROUTINE Delete_DBL_RNK6(A)
    TYPE(DBL_RNK6),INTENT(INOUT) :: A
    INTEGER                      :: Dbls
    Dbls=SIZE(A%D)
    DEALLOCATE(A%D,STAT=MemStatus)
    CALL DecMem(MemStatus,0,Dbls)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_DBL_RNK6

  SUBROUTINE Delete_CHR10_VECT(A)
    TYPE(CHR10_VECT) :: A
    INTEGER        :: MemStatus
    DEALLOCATE(A%C,STAT=MemStatus)
    CALL DecMem(MemStatus,0,0)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_CHR10_VECT

  SUBROUTINE Delete_CHR_VECT(A)
    TYPE(CHR_VECT) :: A
    INTEGER        :: MemStatus
    DEALLOCATE(A%C,STAT=MemStatus)
    CALL DecMem(MemStatus,SIZE(A%C)*DEFAULT_CHR_LEN,0)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_CHR_VECT

  SUBROUTINE Delete_LOG_VECT(A)
    TYPE(LOG_VECT) :: A
    INTEGER        :: MemStatus
    DEALLOCATE(A%L,STAT=MemStatus)
    CALL DecMem(MemStatus,0,0)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_LOG_VECT

  SUBROUTINE Delete_INTC(A)
    TYPE(INTC)     :: A
    IF(AllocQ(A%Alloc)) A%Alloc=ALLOCATED_FALSE
    IF(A%N==0) RETURN
    A%N=0
    CALL Delete(A%Def)
    CALL Delete(A%Atoms)
    CALL Delete(A%Cells)
    CALL Delete(A%Value)
    CALL Delete(A%Constraint)
    CALL Delete(A%ConstrValue)
    CALL Delete(A%Active)
    CALL Delete(A%PredVal)
    CALL Delete(A%PredGrad)
    CALL Delete(A%InvHess)
    CALL DecMem(MemStatus,0,0)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_INTC

  SUBROUTINE Delete_BMATR(A)
    TYPE(BMATR)    :: A
    CALL Delete(A%IB)
    CALL Delete(A%B)
    CALL Delete(A%BL)
    CALL Delete(A%BLI)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_BMATR

  SUBROUTINE Delete_PBCInfo(A)
    TYPE(PBCInfo),INTENT(INOUT)       :: A
    CALL Delete(A%AutoW)
    CALL Delete(A%SuperCell)
    CALL Delete(A%CellCenter)
    CALL Delete(A%TransVec)
    CALL Delete(A%BoxShape)
    CALL Delete(A%InvBoxSh)
    CALL Delete(A%LatFrc)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_PBCInfo

  SUBROUTINE Delete_CRDS(A)
    TYPE(CRDS),INTENT(INOUT)       :: A
    CALL Delete(A%BndBox)
    CALL Delete(A%AtNum)
    CALL Delete(A%AtTyp)
    CALL Delete(A%AtNam)
    CALL Delete(A%DoFreq)
    CALL Delete(A%AtMss)
    CALL Delete(A%CConstrain)
    CALL Delete(A%Carts)
    CALL Delete(A%Velocity)
    CALL Delete(A%Gradients)
    CALL Delete(A%Fext)
    CALL Delete(A%PBC)
    CALL Delete(A%Displ)
    CALL Delete(A%PBCDispl)
    CALL Delete(A%BoxCarts)
    CALL Delete(A%ETotalPerSCF)
    CALL Delete(A%InCells)
    CALL Delete(A%OvCells)
    A%NAtms=0
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_CRDS

  SUBROUTINE Delete_BCSR(A,OnAll_O)
    TYPE(BCSR),INTENT(INOUT) :: A
    LOGICAL,OPTIONAL         :: OnAll_O
    LOGICAL                  :: OnAll
    IF(PRESENT(OnAll_O))THEN
       OnAll=OnAll_O
    ELSE
       OnAll=.FALSE.
    ENDIF
#ifdef PARALLEL
    IF(MyId==ROOT.OR.OnAll)THEN
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

  SUBROUTINE Delete_BSET(A)
    TYPE(BSET),INTENT(INOUT)       :: A
    A%NAtms=0; A%NBasF=0; A%NKind=0
    A%NCtrt=0; A%NPrim=0; A%NASym=0
    A%LMNLen=0
    CALL Delete(A%Kinds)
    CALL Delete(A%AtNam)
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
    IF(A%HasECPs)THEN
       CALL Delete(A%NCoreEl)
       CALL Delete(A%NTyp1PF)
       CALL Delete(A%ProjEll)
       CALL Delete(A%Typ1Ell)
       CALL Delete(A%Typ1Exp)
       CALL Delete(A%Typ1CCo)
       CALL Delete(A%NTyp2PF)
       CALL Delete(A%Typ2Ell)
       CALL Delete(A%Typ2Exp)
       CALL Delete(A%Typ2CCo)
    ENDIF
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_BSET

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

  SUBROUTINE Delete_IDrv(A)
    TYPE(IDrv),INTENT(INOUT)       :: A
    CALL Delete(A%VLOC)
    CALL Delete(A%CDrv)
    CALL Delete(A%SLOC)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_IDrv

  SUBROUTINE Delete_DSL(A)
    TYPE(DSL),INTENT(INOUT)       :: A
    CALL Delete(A%SLDis)
    CALL Delete(A%SLPrm)
    CALL Delete(A%SLKey)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_DSL

  SUBROUTINE Delete_GradD(A)
    TYPE(GradD),INTENT(INOUT)       :: A
    CALL Delete(A%GDrv1)
    CALL Delete(A%GDrv2)
    CALL Delete(A%GDrv3)
    CALL Delete(A%GDrv4)
    CALL Delete(A%GDrv5)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_GradD

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
    CALL AllocChk(A%Alloc)
    CALL New(A%DPole,3)
    CALL New(A%QPole,6)
    A%MPole=0D0
    A%DPole%D=0D0
    A%QPole%D=0D0
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_CMPoles

  SUBROUTINE Initialize_CMPoles(A)
    TYPE(CMPoles), INTENT(INOUT) :: A
    CALL Initialize(A%DPole)
    CALL Initialize(A%QPole)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_CMPoles

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
    INTEGER,OPTIONAL,DIMENSION(4)   :: N_O
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
          IF(N_O(3) /= A%NCoef.OR.N_O(4) /= A%NSDen) THEN
             A%NCoef = N_O(3)
             A%NSDen = N_O(4)
             CALL Delete(A%Co)
             CALL New(A%Co,A%NCoef*A%NSDen) !<<< SPIN
          ENDIF
       ENDIF
    ELSE
       A%Alloc=ALLOCATED_TRUE
       IF(PRESENT(N_O)) THEN
          A%NExpt=N_O(1)
          A%NDist=N_O(2)
          A%NCoef=N_O(3)
          A%NSDen=N_O(4)!<<< SPIN
       ELSE
          A%NExpt=0
          A%NDist=0
          A%NCoef=0
          A%NSDen=1!<<< SPIN
       ENDIF
       CALL New(A%NQ  ,A%NExpt)
       CALL New(A%Lndx,A%NExpt)
       CALL New(A%OffQ,A%NExpt)
       CALL New(A%OffR,A%NExpt)
       CALL New(A%Expt,A%NExpt)
       CALL New(A%Qx  ,A%NDist)
       CALL New(A%Qy  ,A%NDist)
       CALL New(A%Qz  ,A%NDist)
       CALL New(A%Co  ,A%NCoef*A%NSDen)!<<< SPIN
    ENDIF
    !
  END SUBROUTINE New_HGRho

  SUBROUTINE Initialize_HGRho(A)
    TYPE(HGRho), INTENT(INOUT) :: A
    CALL Initialize(A%NQ)
    CALL Initialize(A%Lndx)
    CALL Initialize(A%OffQ)
    CALL Initialize(A%OffR)
    CALL Initialize(A%Expt)
    CALL Initialize(A%Qx)
    CALL Initialize(A%Qy)
    CALL Initialize(A%Qz)
    CALL Initialize(A%Co)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_HGRho

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

  FUNCTION AllocQ(Alloc)
    INTEGER :: Alloc
    LOGICAL :: AllocQ
    AllocQ=.FALSE.
    IF(Alloc==ALLOCATED_TRUE)THEN
       AllocQ=.TRUE.
    ENDIF
    RETURN
  END FUNCTION AllocQ

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
       CALL Halt('Attempt to allocate memory already allocated.')
#endif
    ENDIF
  END SUBROUTINE AllocChk

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

  FUNCTION ToBytes(NInt,NDbl)
    INTEGER,INTENT(IN) :: NInt,NDbl
    INTEGER            :: ToBytes
    ToBytes=NInt*BytesPerInt()+NDbl*BytesPerDbl()
  END FUNCTION ToBytes

  FUNCTION ToMB(NInt,NDbl)
    INTEGER,INTENT(IN) :: NInt,NDbl
    REAL(DOUBLE)       :: ToMB
    ToMB=NInt*IntToMB+NDbl*DblToMB
  END FUNCTION ToMB

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

  FUNCTION BytesPerDbl()
    INTEGER :: BytesPerDbl
    BytesPerDbl=8
  END FUNCTION BytesPerDbl

  !--------------------------------------------------------------------------
  ! Create the CellSet
  !--------------------------------------------------------------------------
  SUBROUTINE New_CellSet(CS, NCells)
    TYPE(CellSet), INTENT(INOUT) :: CS
    INTEGER, INTENT(IN)          :: NCells

    CS%NCells = NCells
    CALL New(CS%CellCarts,(/3,CS%NCells/))
    CS%Alloc=ALLOCATED_TRUE

  END SUBROUTINE New_CellSet

  SUBROUTINE Initialize_CellSet(A)
    TYPE(CellSet), INTENT(INOUT) :: A

    CALL Initialize(A%CellCarts)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_CellSet

  !--------------------------------------------------------------------------
  ! Delete the CellSet
  !--------------------------------------------------------------------------
  SUBROUTINE Delete_CellSet(CS)
    TYPE(CellSet)                  :: CS

    !CALL MondoLog(DEBUG_NONE, "Delete_CellSet", "deleting CellSet")
    IF(AllocQ(CS%Alloc)) THEN
       CS%Alloc  = ALLOCATED_FALSE
       CS%NCells = 0
       CALL Delete(CS%CellCarts)
    ENDIF

  END SUBROUTINE Delete_CellSet

  SUBROUTINE New_PBCFit(A,MaxMem)
    TYPE(PBCFits) :: A
    INTEGER       :: MaxMem
    !
    A%MaxMem=MaxMem
    A%ActMem=0
    CALL New(A%AWeights,A%MaxMem)
    CALL New(A%PBCValues,(/9,A%MaxMem/))
    CALL New(A%PBCGrads,(/9,A%MaxMem/))
    A%Alloc=ALLOCATED_TRUE
  END SUBROUTINE New_PBCFit

  SUBROUTINE Initialize_PBCFit(A)
    TYPE(PBCFits), INTENT(INOUT) :: A
    CALL Initialize(A%AWeights)
    CALL Initialize(A%PBCValues)
    CALL Initialize(A%PBCGrads)
    A%Alloc = ALLOCATED_FALSE
  END SUBROUTINE Initialize_PBCFit

  SUBROUTINE Delete_PBCFit(A)
    TYPE(PBCFits) :: A
    !
    A%MaxMem=0
    A%ActMem=0
    CALL Delete(A%AWeights)
    CALL Delete(A%PBCValues)
    CALL Delete(A%PBCGrads)
    A%Alloc=ALLOCATED_FALSE
  END SUBROUTINE Delete_PBCFit

END MODULE MemMan
