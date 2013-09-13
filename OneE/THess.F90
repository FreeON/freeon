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
PROGRAM THess
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
  USE BlokTrPd2T
#ifdef PARALLEL
  USE MondoMPI
  TYPE(DBL_RNK2)              :: THssTmp
  INTEGER                     :: IErr
#endif
  TYPE(BCSR)                  :: P
  INTEGER                     :: NC,I,J
  REAL(DOUBLE)                :: B(3),F_nlm(21)
  TYPE(AtomPair)              :: Pair
  TYPE(BSET)                  :: BS
  TYPE(CRDS)                  :: GM
  TYPE(ARGMT)                 :: Args
  INTEGER                     :: Q,R,AtA,AtB,JP,MB,MA,NB,MN1,A1,A2
  TYPE(HGRho)                 :: Rho
  TYPE(DBL_RNK2)              :: THss
  REAL(DOUBLE)                :: TFrcChk
  CHARACTER(LEN=*),PARAMETER  :: Prog='THess'
!-------------------------------------------------------------------------------------
! Start up macro

  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
!----------------------------------------------
! Get basis set and geometry

  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
#ifdef PARALLEL
  CALL New(P,OnAll_O=.TRUE.)
#endif
  CALL Get(P,TrixFile('D',Args,1),BCast_O=.TRUE.)
  CALL New(THss,(/3*NAtoms,3*NAtoms/))
  CALL DBL_VECT_EQ_DBL_SCLR(9*NAtoms**2,THss%D(1,1),0.0D0)
!--------------------------------------------------------------------------------
! THess=2*Tr{P.dT} (Extra 2 to account for symmetry of T in the trace)
!--------------------------------------------------------------------------------
#ifdef PARALLEL
  DO AtA=Beg%I(MyId),End%I(MyId)
#else
  DO AtA=1,NAtoms
#endif
     A1=3*(AtA-1)+1
     MA=BSiz%I(AtA)
     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
        AtB=P%ColPt%I(JP)
        A2=3*(AtB-1)+1
        IF(SetAtomPair(GM,BS,AtA,AtB,Pair))THEN
           Q=P%BlkPt%I(JP)
           NB=BSiz%I(AtB)
           MN1=MA*NB-1
           B=Pair%B
           DO NC=1,CS_OUT%NCells
              Pair%B=B+CS_OUT%CellCarts%D(:,NC)
              Pair%AB2=(Pair%A(1)-Pair%B(1))**2  &
                      +(Pair%A(2)-Pair%B(2))**2  &
                      +(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                 IF(.NOT.Pair%SameAtom) THEN
                    F_nlm(1:21) = TrPd2T(BS,Pair,P%MTrix%D(Q:Q+MN1))
                    !AB terms
                    THss%D(A1  ,A2  )=THss%D(A1  ,A2  )+4D0*F_nlm(1) !AxBx
                    THss%D(A1+1,A2  )=THss%D(A1+1,A2  )+4D0*F_nlm(2) !AxBy
                    THss%D(A1+2,A2  )=THss%D(A1+2,A2  )+4D0*F_nlm(3) !AxBz
                    THss%D(A1  ,A2+1)=THss%D(A1  ,A2+1)+4D0*F_nlm(4) !AyBx
                    THss%D(A1+1,A2+1)=THss%D(A1+1,A2+1)+4D0*F_nlm(5) !AyBy
                    THss%D(A1+2,A2+1)=THss%D(A1+2,A2+1)+4D0*F_nlm(6) !AyBz
                    THss%D(A1  ,A2+2)=THss%D(A1  ,A2+2)+4D0*F_nlm(7) !AzBx
                    THss%D(A1+1,A2+2)=THss%D(A1+1,A2+2)+4D0*F_nlm(8) !AzBy
                    THss%D(A1+2,A2+2)=THss%D(A1+2,A2+2)+4D0*F_nlm(9) !AzBz
                    !AA terms
                    THss%D(A1  ,A1  )=THss%D(A1  ,A1  )+2D0*F_nlm(10) !AxAx
                    THss%D(A1  ,A1+1)=THss%D(A1  ,A1+1)+2D0*F_nlm(11) !AxAy
                    THss%D(A1  ,A1+2)=THss%D(A1  ,A1+2)+2D0*F_nlm(12) !AxAz
                    THss%D(A1+1,A1  )=THss%D(A1+1,A1  )+2D0*F_nlm(11) !AyAx
                    THss%D(A1+1,A1+1)=THss%D(A1+1,A1+1)+2D0*F_nlm(13) !AyAy
                    THss%D(A1+1,A1+2)=THss%D(A1+1,A1+2)+2D0*F_nlm(14) !AyAz
                    THss%D(A1+2,A1  )=THss%D(A1+2,A1  )+2D0*F_nlm(12) !AzAx
                    THss%D(A1+2,A1+1)=THss%D(A1+2,A1+1)+2D0*F_nlm(14) !AzAy
                    THss%D(A1+2,A1+2)=THss%D(A1+2,A1+2)+2D0*F_nlm(15) !AzAz
                    !BB terms
                    THss%D(A2  ,A2  )=THss%D(A2  ,A2  )+2D0*F_nlm(16) !BxBx
                    THss%D(A2  ,A2+1)=THss%D(A2  ,A2+1)+2D0*F_nlm(17) !BxBy
                    THss%D(A2  ,A2+2)=THss%D(A2  ,A2+2)+2D0*F_nlm(18) !BxBz
                    THss%D(A2+1,A2  )=THss%D(A2+1,A2  )+2D0*F_nlm(17) !ByBx
                    THss%D(A2+1,A2+1)=THss%D(A2+1,A2+1)+2D0*F_nlm(19) !ByBy
                    THss%D(A2+1,A2+2)=THss%D(A2+1,A2+2)+2D0*F_nlm(20) !ByBz
                    THss%D(A2+2,A2  )=THss%D(A2+2,A2  )+2D0*F_nlm(18) !BzBx
                    THss%D(A2+2,A2+1)=THss%D(A2+2,A2+1)+2D0*F_nlm(20) !BzBy
                    THss%D(A2+2,A2+2)=THss%D(A2+2,A2+2)+2D0*F_nlm(21) !BzBz
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO
!--------------------------------------------------------------------------------
#ifdef PARALLEL
  CALL New(THssTmp,(/3*NAtoms,3*NAtoms/))
  CALL MPI_Reduce(THss%D(1,1),THssTmp%D(1,1),9*NAtoms**2,MPI_DOUBLE_PRECISION,MPI_SUM,0,MONDO_COMM,IErr)
  IF(MyID == 0) THEN
     CALL DCOPY(9*NAtoms**2,THssTmp%D(1,1),1,THss%D(1,1),1)
  ENDIF
  CALL Delete(THssTmp)
#endif
! Do some checksumming, resumming and IO
  !CALL PChkSum(THss,'d2T/dR2',Proc_O=Prog)   !<-- should write a check sum for dbl_rnk2
  CALL Print_DBL_RNK2(THss,'THessian',Unit_O=6)
  !Should save the hessian somewhere.
  !
!------------------------------------------------------------------------------
! Tidy up
  CALL Delete(THss)
  CALL Delete(P)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL ShutDown(Prog)
END PROGRAM THess
