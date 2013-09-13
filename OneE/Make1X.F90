!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
MODULE Make1X
!H=================================================================================
!H MODULE Make1X
!H This MODULE contains:
!H  PUBLIC:
!H  o SUB
!H
!H  PRIVATE:
!H  o SUB
!H
!H  OPTIONS:
!H
!H
!H  Comments:
!H
!H---------------------------------------------------------------------------------
  !
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE AtomPairs
  USE Int1E
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
  !
  PRIVATE
!---------------------------------------------------------------------------------
! PUBLIC DECLARATIONS
!---------------------------------------------------------------------------------
  PUBLIC  :: MakeS,MakeT,MakeD,MakeV
  !
!---------------------------------------------------------------------------------
! PRIVATE DECLARATIONS
!---------------------------------------------------------------------------------
  PRIVATE :: PreAlloc1E
  !
CONTAINS
  !
  !
  SUBROUTINE PreAlloc1E(GM,BS,NBlks,NNon0)
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BSET)     :: BS
    TYPE(CRDS)     :: GM
    INTEGER        :: NBlks,NNon0
    !-------------------------------------------------------------------
    TYPE(AtomPair) :: Pair
    INTEGER        :: AtA,AtB,NN,R,P
    !-------------------------------------------------------------------
    P=1;R=1;
#ifdef PARALLEL
    DO AtA=Beg%I(MyId),End%I(MyId)
#else
    DO AtA=1,NAtoms
#endif
       DO AtB=1,NAtoms
          IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
             NN = Pair%NA*Pair%NB
             ! May do more agresive test here.
             R=R+NN
             P=P+1
          ENDIF
       ENDDO
    ENDDO
    NBlks=P
    NNon0=R
    !
    WRITE(*,*) 'Default setting'
    WRITE(*,*) 'MaxBlks=',MaxBlks,' MaxNon0=',MaxNon0,' NBasF',NBasF
    WRITE(*,*) 'Optimized setting'
    WRITE(*,*) 'NBlks  =',NBlks  ,' NNon0  =',NNon0
    !
  END SUBROUTINE PreAlloc1E
  !
  !
  SUBROUTINE MakeS(Args,GM,BS)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)                :: S
#else
    TYPE(BCSR)                 :: S
#endif
    INTEGER                    :: NC
    REAL(DOUBLE),DIMENSION(3)  :: B
    TYPE(AtomPair)             :: Pair
    TYPE(BSET)                 :: BS
    TYPE(CRDS)                 :: GM
    TYPE(ARGMT)                :: Args
    INTEGER                    :: P,R,AtA,AtB,NN,NCols,NElem
    CHARACTER(LEN=*),PARAMETER :: Prog='MakeS'
!----------------------------------------------
! Allocations
    CALL PreAlloc1E(GM,BS,NCols,NElem)
    CALL New(S,(/NAtoms+1,NCols,NElem/))
!-----------------------------------------------
! Initialize the matrix and associated indecies

    P=1; R=1; S%RowPt%I(1)=1
    CALL SetEq(S%MTrix,Zero)
!-----------------------------------------------
! Main loops
#ifdef PARALLEL
    S%NAtms=0
    DO AtA=Beg%I(MyId),End%I(MyId)
       S%NAtms=S%NAtms+1
#else
    S%NAtms=NAtoms
    DO AtA=1,NAtoms
#endif
       DO AtB=1,NAtoms
          IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
             NN = Pair%NA*Pair%NB
             B = Pair%B
             DO NC = 1,CS_OUT%NCells
                Pair%B = B+CS_OUT%CellCarts%D(:,NC)
                Pair%AB2  = (Pair%A(1)-Pair%B(1))**2 &
                     + (Pair%A(2)-Pair%B(2))**2 &
                     + (Pair%A(3)-Pair%B(3))**2
                IF(TestAtomPair(Pair)) THEN
                   CALL I1E_Ovlap_AtBlk(BS,Pair,S%MTrix%D(R))
                ENDIF
             ENDDO
             S%ColPt%I(P)=AtB
             S%BlkPt%I(P)=R
             R=R+NN
             P=P+1
#ifdef PARALLEL
             S%RowPt%I(S%NAtms+1)=P
             IF(R>NElem.OR.P>NCols) THEN
                !IF(R>MaxNon0Node.OR.P>MaxBlksNode)THEN
                WRITE(*,*)' MyId = ',MyId,MaxBlksNode,MaxNon0Node
                WRITE(*,*)' MyId = ',MyId,' P = ',P,' R = ',R
                CALL Halt(' DBCSR dimensions blown in MakeS ')
             ENDIF
#else
             S%RowPt%I(AtA+1)=P
             IF(R>NElem.OR.P>NCols) THEN
                !IF(R>MaxNon0.OR.P>MaxBlks) THEN
                WRITE(*,*) 'MakeS: MaxNon0=',MaxNon0
                WRITE(*,*) 'MakeS: R=',R
                WRITE(*,*) 'MakeS: MaxBlks=',MaxBlks
                WRITE(*,*) 'MakeS: P=',P
                CALL Halt(' BCSR dimensions blown in MakeS ')
             ENDIF
#endif
          ENDIF
       ENDDO
    ENDDO
    S%NBlks=P-1
    S%NNon0=R-1
!------------------------------------------------------------
! Put S to disk
    Thresholds%Trix = Thresholds%Trix*1.D-2
    CALL Filter(S)
    write(*,*) 'Filtered: S%NBlks',S%NBlks,'S%NNon0',S%NNon0
    Thresholds%Trix = Thresholds%Trix*1.D2
    IF(SCFActn=='RestartGeomSwitch') THEN
       CALL Put(S,TrixFile('S',Args,Stats_O=(/Current(1),Current(2),Current(3)-1/)))
    ELSE
       CALL Put(S,TrixFile('S',Args))
    ENDIF
!-----------------------------------------------------------
! Printing
    CALL PChkSum(S,'S',Prog)
    CALL PPrint( S,'S')
    CALL Plot(   S,'S')
!------------------------------------------------------------
! Tidy up
    CALL Delete(S)
! didn't count flops, any accumulation is residual
! from matrix routines
    PerfMon%FLOP=Zero
  END SUBROUTINE MakeS
  !
  !
  SUBROUTINE MakeT(Args,GM,BS)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)         :: T
#else
    TYPE(BCSR)          :: T
#endif
    INTEGER                   :: NC
    REAL(DOUBLE),DIMENSION(3) :: B
    TYPE(AtomPair)            :: Pair
!
    TYPE(BSET)          :: BS
    TYPE(CRDS)          :: GM
!
    TYPE(ARGMT)         :: Args
    INTEGER             :: P,R,AtA,AtB,NN,NCols,NElem
    CHARACTER(LEN=*),PARAMETER :: Prog='MakeT'
!----------------------------------------------
! Allocations
!
    CALL PreAlloc1E(GM,BS,NCols,NElem)
    CALL New(T,(/NAtoms+1,NCols,NElem/))
!-----------------------------------------------
! Initialize the matrix and associated indecies
!
    P=1; R=1; T%RowPt%I(1)=1
    CALL SetEq(T%MTrix,Zero)
!-----------------------------------------------
! Main loops
!
#ifdef PARALLEL
    T%NAtms=0
    DO AtA=Beg%I(MyId),End%I(MyId)
       T%NAtms=T%NAtms+1
#else
    T%NAtms=NAtoms
    DO AtA=1,NAtoms
#endif
       DO AtB=1,NAtoms
          IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
             NN = Pair%NA*Pair%NB
             B = Pair%B
             DO NC = 1,CS_OUT%NCells
                Pair%B = B+CS_OUT%CellCarts%D(:,NC)
                Pair%AB2  = (Pair%A(1)-Pair%B(1))**2 &
                     + (Pair%A(2)-Pair%B(2))**2 &
                     + (Pair%A(3)-Pair%B(3))**2
                IF(TestAtomPair(Pair)) THEN
                   CALL I1E_Kinet_AtBlk(BS,Pair,T%MTrix%D(R))
                ENDIF
             ENDDO
             T%ColPt%I(P)=AtB
             T%BlkPt%I(P)=R
             R=R+NN
             P=P+1
#ifdef PARALLEL
             T%RowPt%I(T%NAtms+1)=P
             IF(R>NElem.OR.P>NCols) THEN
                !IF(R>MaxNon0Node.OR.P>MaxBlksNode)THEN
                WRITE(*,*)' MyId = ',MyId,MaxBlksNode,MaxNon0Node
                WRITE(*,*)' MyId = ',MyId,' P = ',P,' R = ',R
                CALL Halt(' DBCSR dimensions blown in MakeT ')
             ENDIF
#else
             T%RowPt%I(AtA+1)=P
             IF(R>NElem.OR.P>NCols) &
                  !IF(R>MaxNon0.OR.P>MaxBlks) &
                  CALL Halt(' BCSR dimensions blown in MakeT ')
#endif
          ENDIF
       ENDDO
    ENDDO
    T%NBlks=P-1
    T%NNon0=R-1
!------------------------------------------------------------
! Put T to disk
!
    CALL Filter(T)
    !write(*,*) 'Filtered: T%NBlks',T%NBlks,'T%NNon0',T%NNon0
    CALL Put(T,TrixFile('T',Args))
!-----------------------------------------------------------
! Printing
!
    CALL PChkSum(T,'T',Prog)
    CALL PPrint( T,'T')
    CALL Plot(   T,'T')
!------------------------------------------------------------
! Tidy up
!
    CALL Delete(T)
! didn't count flops, any accumulation is residual
! from matrix routines
    PerfMon%FLOP=Zero
  END SUBROUTINE MakeT
  !
  !
  SUBROUTINE MakeD(Args,GM,BS)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)                :: M
#else
    TYPE(BCSR)                 :: M
#endif
    INTEGER                    :: NCA,NCB
    REAL(DOUBLE), DIMENSION(3) :: B,A
    TYPE(AtomPair)             :: Pair
!
    TYPE(BSET)                 :: BS
    TYPE(CRDS)                 :: GM
!
    TYPE(ARGMT)                :: Args
    TYPE(DBL_VECT)             :: COrig
    INTEGER                    :: P,R,AtA,AtB,NN,NCols,NElem
    INTEGER                    :: IXYZ
    CHARACTER(LEN=*)              , PARAMETER :: Prog='MakeD'
    CHARACTER(LEN=*), DIMENSION(3), PARAMETER :: Cart=(/'X','Y','Z'/)
!
! Get Multipole origine. TODO
    CALL New(COrig,3)
    CALL SetEQ(COrig,0.0d0)
  !COrig%D(:)=GM%PBC%CellCenter%D(:)
!----------------------------------------------
! Allocations
!
    CALL PreAlloc1E(GM,BS,NCols,NElem)
    CALL New(M,(/NAtoms+1,NCols,NElem/))
!-----------------------------------------------
! Run over cartesian componants
!
    DO IXYZ=1,3
!-----------------------------------------------
! Initialize the matrix and associated indecies
!
       P=1; R=1; M%RowPt%I(1)=1
       CALL SetEq(M%MTrix,Zero)
!-----------------------------------------------
! Main loops
!
#ifdef PARALLEL
       M%NAtms=0
       DO AtA=Beg%I(MyId),End%I(MyId)
          M%NAtms=M%NAtms+1
#else
       M%NAtms=NAtoms
       DO AtA=1,NAtoms
#endif
          DO AtB=1,NAtoms
             IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
                NN = Pair%NA*Pair%NB
!--------------------------------------------------------
                A(:) = Pair%A(:)
                B(:) = Pair%B(:)
                DO NCA = 1,CS_OUT%NCells
                   Pair%A(:) = A(:)+CS_OUT%CellCarts%D(:,NCA)
                   Pair%AB2  = (Pair%A(1)-Pair%B(1))**2 &
                             + (Pair%A(2)-Pair%B(2))**2 &
                             + (Pair%A(3)-Pair%B(3))**2
!--------------------------------------------------------
                   IF(TestAtomPair(Pair)) THEN
                      CALL I1E_EDipl_AtBlk(BS,Pair,IXYZ,COrig%D(1),M%MTrix%D(R))
                   ENDIF
                ENDDO
                M%ColPt%I(P)=AtB
                M%BlkPt%I(P)=R

                R=R+NN
                P=P+1
#ifdef PARALLEL
                M%RowPt%I(M%NAtms+1)=P
                IF(R>NElem.OR.P>NCols) THEN
                !IF(R>MaxNon0Node.OR.P>MaxBlksNode)THEN
                   WRITE(*,*)' MyId = ',MyId,MaxBlksNode,MaxNon0Node
                   WRITE(*,*)' MyId = ',MyId,' P = ',P,' R = ',R
                   CALL Halt(' DBCSR dimensions blown in MakeD ')
                ENDIF
#else
                M%RowPt%I(AtA+1)=P
                IF(R>NElem.OR.P>NCols) &
                     !IF(R>MaxNon0.OR.P>MaxBlks) &
                     CALL Halt(' BCSR dimensions blown in MakeD ')
#endif
             ENDIF
          ENDDO
       ENDDO
       M%NBlks=P-1
       M%NNon0=R-1
!------------------------------------------------------------
! Put D to disk
!
       CALL Filter(M)
       !write(*,*) 'Filtered: M%NBlks',M%NBlks,'M%NNon0',M%NNon0
       CALL Put(M,TrixFile('Dipole'//Cart(IXYZ),Args))
!------------------------------------------------------------
! Printing
!
       CALL PChkSum(M,'Dipole'//Cart(IXYZ),Prog)
       CALL PPrint( M,'Dipole'//Cart(IXYZ))
       CALL Plot(   M,'Dipole'//Cart(IXYZ))
!     IF(Cart(IXYZ)=='X') CALL Print_BCSR(M,'Dipole X',Unit_O=6)
!     IF(Cart(IXYZ)=='Y') CALL Print_BCSR(M,'Dipole Y',Unit_O=6)
!     IF(Cart(IXYZ)=='Z') CALL Print_BCSR(M,'Dipole Z',Unit_O=6)
!------------------------------------------------------------
    ENDDO ! End Loop over Cartesian Componants
!------------------------------------------------------------
! Tidy up
!
    CALL Delete(M )
    CALL Delete(COrig)
!
! didn't count flops, any accumulation is residual
! from matrix routines
    PerfMon%FLOP=Zero
  END SUBROUTINE MakeD
  !
  !
  SUBROUTINE MakeV(Args,GM,BS)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)         :: V,T1,T2
#else
    TYPE(BCSR)          :: V,T1,T2
#endif
    INTEGER                   :: NC
    REAL(DOUBLE),DIMENSION(3) :: B
    TYPE(AtomPair)            :: Pair
!
    TYPE(BSET)          :: BS
    TYPE(CRDS)          :: GM
!
    TYPE(ARGMT)         :: Args
    INTEGER             :: P,R,AtA,AtB,NN,NCols,NElem
    CHARACTER(LEN=*),PARAMETER :: Prog='MakeV'
real(double) :: tt1,tt2
!----------------------------------------------
! Allocations
!
    CALL PreAlloc1E(GM,BS,NCols,NElem)
    CALL New(V,(/NAtoms+1,NCols,NElem/))
!-----------------------------------------------
! Initialize the matrix and associated indecies
!
    P=1; R=1; V%RowPt%I(1)=1
    CALL SetEq(V%MTrix,Zero)
!-----------------------------------------------
! Main loops
!
call cpu_time(tt1)
#ifdef PARALLEL
    V%NAtms=0
    DO AtA=Beg%I(MyId),End%I(MyId)
       V%NAtms=V%NAtms+1
#else
    V%NAtms=NAtoms
    DO AtA=1,NAtoms
#endif
       DO AtB=1,NAtoms
          IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
             NN = Pair%NA*Pair%NB
             !B = Pair%B
             !DO NC = 1,CS_OUT%NCells
             !   Pair%B = B+CS_OUT%CellCarts%D(:,NC)
             !   Pair%AB2  = (Pair%A(1)-Pair%B(1))**2 &
             !        + (Pair%A(2)-Pair%B(2))**2 &
             !        + (Pair%A(3)-Pair%B(3))**2
             IF(CS_OUT%NCells.NE.1) CALL Halt(' MakeV doesn''t work with PBC ')
             !IF(TestAtomPair(Pair)) THEN
             CALL I1E_NucAtt_AtBlk(GM,BS,Pair,V%MTrix%D(R))
             !call I1E_MssVel_AtBlk(BS,Pair,V%MTrix%D(R))
             !ENDIF
             !ENDDO
             V%ColPt%I(P)=AtB
             V%BlkPt%I(P)=R
             R=R+NN
             P=P+1
#ifdef PARALLEL
             V%RowPt%I(V%NAtms+1)=P
             IF(R>NElem.OR.P>NCols) THEN
                !IF(R>MaxNon0Node.OR.P>MaxBlksNode)THEN
                WRITE(*,*)' MyId = ',MyId,MaxBlksNode,MaxNon0Node
                WRITE(*,*)' MyId = ',MyId,' P = ',P,' R = ',R
                CALL Halt(' DBCSR dimensions blown in MakeV ')
             ENDIF
#else
             V%RowPt%I(AtA+1)=P
             IF(R>NElem.OR.P>NCols) &
                  !IF(R>MaxNon0.OR.P>MaxBlks) &
                  CALL Halt(' BCSR dimensions blown in MakeV ')
#endif
          ENDIF
       ENDDO
    ENDDO
    V%NBlks=P-1
    V%NNon0=R-1

call cpu_time(tt2)
write(*,*) 'time=',tt2-tt1
!------------------------------------------------------------
! Put T to disk
!
    CALL Filter(V)
    !write(*,*) 'Filtered: T%NBlks',T%NBlks,'T%NNon0',T%NNon0
    CALL Put(V,TrixFile('V',Args))
!-----------------------------------------------------------
! Printing
!
    CALL PChkSum(V,'V',Prog)
    CALL PPrint( V,'V')
    CALL Plot(   V,'V')
!------------------------------------------------------------
! Tidy up
!
    !call print_bcsr(V,'XXXX',Unit_o=6)
    !CALL Get(T1,TrixFile('T',Args))
    !CALL Add(V,T1,T2)
    !call print_bcsr(T2,'T+V',Unit_o=6)
    !call delete(T1)
    !call delete(T2)

    CALL Delete(V)
! didn't count flops, any accumulation is residual
! from matrix routines
    PerfMon%FLOP=Zero
  END SUBROUTINE MakeV
  !
END MODULE Make1X
