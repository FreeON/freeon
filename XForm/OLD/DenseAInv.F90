PROGRAM AINV
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE BlokOpsDns
  USE DGEMMS
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
  TYPE(BCSR)                                :: S,SZ
  TYPE(DBL_RNK2)                            :: A,Z,ZT,D
  TYPE(BSET)                                :: BS
  TYPE(CRDS)                                :: GM
  TYPE(INT_VECT)                            :: Flag1,Flag2,RowPt,ColPt,BlkPt
  TYPE(DBL_RNK4)                            :: MD
  TYPE(ARGMT)                               :: Args
  REAL(DOUBLE),ALLOCATABLE,DIMENSION(:,:)   :: Blk1,Blk2,Blk3,Blk8
  REAL(DOUBLE),ALLOCATABLE,DIMENSION(:,:,:) :: P
  INTEGER                                   :: I,J,K,L,M,N,DI,DJ,NI,NJ,OI,OJ,NK,OK,DK,MaxBlkSize
  CHARACTER(LEN=7),PARAMETER                :: Prog='DnsAInv'
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: ZIChar,ZJChar,AIChar, AJChar
!---------------------------------------------------------------------------------------- 
! Start up macro
!
  CALL StartUp(Args,Prog)
!----------------------------------------------
! Get basis set and geometry
!
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
!---------------------------------------------- 
! Allocations 
!
!   NBasF=3
!   NAtoms=2
!
  CALL New(A,(/NBasF,NBasF/))
  CALL New(Z,(/NBasF,NBasF/))
  CALL New(ZT,(/NBasF,NBasF/))
  CALL New(D,(/NBasF,NBasF/))

   CALL New(S)
   CALL Get(S,TrixFile('S',Args))
   CALL SetEq(A,S)

   CALL Delete(S)    
!  

!  CALL New(BSiz,NAtoms)
!  CALL New(OffS,NAtoms)
!  OffS%I=(/1,2/)    
!  BSiz%I=(/1,2/)
  MaxBlkSize=0

  DO I=1,NAtoms
     MaxBlkSize=MAX(MaxBlkSize,BSiz%I(I))
  ENDDO

  CALL SetDSYEVWork(MaxBlkSize)

  ALLOCATE(Blk1(MaxBlkSize,MaxBlkSize))
  ALLOCATE(Blk2(MaxBlkSize,MaxBlkSize))
  ALLOCATE(Blk3(MaxBlkSize,MaxBlkSize))
  ALLOCATE(Blk8(MaxBlkSize,MaxBlkSize))
  ALLOCATE(P(MaxBlkSize,MaxBlkSize,NAtoms))
!
! Z=I
  N=NBasF
  DO I=1,N
     Z%D(I,I)=One
  ENDDO
  DO I=1,NAtoms
     NI=BSiz%I(I); OI=OffS%I(I); DI=OI+NI-1
     DO J=1,I-1
        WRITE(*,*)'=========================================================='
        WRITE(*,*)' ICol = ',I,' JCol = ',J
        WRITE(*,*)'=========================================================='
        NJ=BSiz%I(J); OJ=OffS%I(J); DJ=OJ+NJ-1
!       Blk1=P^(j-1)_i=[A^t_j].[Z^(j-1)_i]
        Blk1=Zero

        CALL PPPrintT(N,NJ,A%D(1:N,OJ:DJ),'A^T_J')
        CALL PPPrint(N,NI,Z%D(1:N,OI:DI),'Z_I')

        CALL DGEMM__TN(N,NJ,NI,One,A%D(1:N,OJ:DJ),Z%D(1:N,OI:DI),Blk1(1:NJ,1:NI))

        CALL PPPrint(NJ,NI,Blk1(1:NJ,1:NI),'Blk1bb')
!       Blk2=[P^(j-1)_j]^(-1).[P^(j-1)_i]

        Blk2=Zero
        CALL DGEMM__NN(NJ,NJ,NI,One,P(1:NJ,1:NJ,J),Blk1(1:NJ,1:NI),Blk2(1:NJ,1:NI))
        CALL PPPrint(NJ,NI,Blk2(1:NJ,1:NI),'Blk2')
!       Z^j_i=Z^(j-1)_i-[Z^(j-1)_j].[P^(j-1)_j]^(-1).[P^(j-1)_i]
!        CALL PPPrint(N,NJ,Z%D(1:N,OJ:DJ),'Z(J-1)_J')
!        CALL PPPrint(N,NI,Z%D(1:N,OI:DI),'Z(J-1)_I')
!        CALL DGEMM__NN(N,NJ,NI,-One,Z%D(1:N,OJ:DJ),Blk2(1:NJ,1:NI),Z%D(1:N,OI:DI))
!        CALL PPPrint(N,NI,Z%D(1:N,OI:DI),'Z(J)_I')


        DO K=1,NAtoms
           NK=BSiz%I(K); OK=OffS%I(K); DK=OK+NK-1
           IF(DDot(NK*NJ,Z%D(OK:DK,OJ:DJ),1,Z%D(OK:DK,OJ:DJ),1)/=Zero)THEN
           CALL DGEMM__NN(NK,NJ,NI,-One,Z%D(OK:DK,OJ:DJ),Blk2(1:NJ,1:NI),Z%D(OK:DK,OI:DI))
           ZJChar='Z^(J-1)_('//TRIM(IntToChar(J))//','//TRIM(IntToChar(K))//')'
           ZIChar='Z^(J)_('//TRIM(IntToChar(I))//','//TRIM(IntToChar(K))//')'
           CALL PPPrint(NK,NJ,Z%D(OK:DK,OJ:DJ),ZJChar)
           CALL PPPrint(NK,NI,Z%D(OK:DK,OI:DI),ZIChar)
           ENDIF
         ENDDO

      ENDDO
!     Blk1=P^(i-1)_i=[A^t_i].[Z^(i-1)_i]
#ifdef CORRECT
      Blk1=Zero
      CALL PPPrintT(N,NI,A%D(1:N,OI:DI),'A^T_i')
      CALL PPPrint(N,NI,Z%D(1:N,OI:DI),'Z_i')
      CALL DGEMM__TN(N,NI,NI,One,A%D(1:N,OI:DI),Z%D(1:N,OI:DI),Blk1(1:NI,1:NI))
!     [P^(i-1)_i]^(-1)
      CALL PPPrint(NI,NI,Blk1(1:NI,1:NI),'Blk1')
#else
      Blk1=Zero
      DO K=1,NAtoms
         NK=BSiz%I(K); OK=OffS%I(K); DK=OK+NK-1
         Blk8=Zero
         AIChar='A^T_('//TRIM(IntToChar(I))//','//TRIM(IntToChar(K))//')'
         ZIChar='Z^(I-1)_('//TRIM(IntToChar(I))//','//TRIM(IntToChar(K))//')'
         CALL PPPrintT(NK,NI,A%D(OK:DK,OI:DI),AIChar)
         CALL PPPrint(NK,NI,Z%D(OK:DK,OI:DI),ZIChar)
         CALL DGEMM__TN(NK,NI,NI,One,A%D(OK:DK,OI:DI),Z%D(OK:DK,OI:DI),Blk8(1:NI,1:NI))
         Blk1(1:NI,1:NI)=Blk1(1:NI,1:NI)+Blk8(1:NI,1:NI)
!        [P^(i-1)_i]^(-1)
         CALL PPPrint(NI,NI,Blk8(1:NI,1:NI),TRIM(AIChar)//'.'//TRIM(ZIChar))
      ENDDO
      CALL PPPrint(NI,NI,Blk1(1:NI,1:NI),'Blk1')
#endif
      CALL FunkOnSqMat(NI,Inverse,Blk1(1:NI,1:NI),P(1:NI,1:NI,I))
      CALL PPPrint(NI,NI,P(1:NI,1:NI,I),'InvBlk1')
   ENDDO

  PrintFlags%Mat=DEBUG_MATRICES
!  Z%D=TRANSPOSE(Z%D)
  CALL PPRINT(Z,'ZBefore',Unit_O=6)
!  Z%D=TRANSPOSE(Z%D)
!  STOP 'LDJFDLFJ'

!
   D%D=Zero
   DO I=1,NAtoms
      NI=BSiz%I(I); OI=OffS%I(I); DI=OI+NI-1
      CALL PPPrint(NI,NI,P(1:NI,1:NI,I),' DBlok Before ')    
      CALL FunkOnSqMat(NI,SqRoot,P(1:NI,1:NI,I),D%D(OI:DI,OI:DI))
      CALL PPPrint(NI,NI,D%D(OI:DI,OI:DI),' DBlok After ')
   ENDDO

!   D%D=TRANSPOSE(D%D)
   CALL PPrint(D,'[D]^(-1)',Unit_O=6)   
!   D%D=TRANSPOSE(D%D)
!
   CALL Put(Z,'DnsZ')
   CALL Put(D,'DnsD')

   ZT%D=Z%D
   Z%D=MATMUL(ZT%D,D%D)

   CALL Put(Z,'DnsZD')

   CALL PPrint(Z,'[Z]',Unit_O=6)


   ZT%D=TRANSPOSE(Z%D)
   CALL PPrint(ZT,'[ZT]',Unit_O=6)
   D%D=MATMUL(Z%D,MATMUL(ZT%D,A%D))
   Z%D=Zero
   DO I=1,NBasF
      Z%D(I,I)=-1.0D0
   ENDDO
   D%D=D%D+Z%D
   CALL PPrint(D,'[D0]',Unit_O=6)

!  Consistency check
   DO I=1,N
      IF((D%D(I,I)-One)>1.D-10)THEN
         WRITE(*,*)' D(I,I) = ',D%D(I,I)
         CALL Halt(' Consistency check failed in DenseAInv! ')
      ENDIF
   ENDDO
   CALL Delete(D)
!------------------------------------------------------------
!  Put Z and ZT to disk
!  
   CALL New(S)
   CALL New(SZ)
   CALL SetEq(S,Z)
   CALL Filter(SZ,S)
   CALL Put(SZ,TrixFile('Z',Args))
!-----------------------------------------------------------
! Printing
!
   CALL PChkSum(SZ,'Z',Prog)
   CALL PPrint( SZ,'Z')
   CALL Plot(   SZ,'Z')
!-----------------------------------------------------------
   CALL SetEq(S,ZT)
   CALL Filter(SZ,S)
   CALL Put(SZ,TrixFile('ZT',Args))
!
   CALL ShutDown(Prog)
END PROGRAM
