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
! SPARSE-BLOCKED O(N) AINV WITH DISTANCE THRESHOLDING
! Author: Matt Challacombe
!-----------------------------------------------------------------
MODULE AInv
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE MemMan
  USE LinAlg
  USE MatFunk
  IMPLICIT NONE
CONTAINS
  SUBROUTINE BlockedAInv(A,TrixThresh,GM,DstncThresh,Z,Zt,Perf)
    TYPE(BCSR)          :: A,Z,Zt,DiagD
#ifdef FIND_CONDA
    TYPE(DBL_RNK2)      :: B,C
#endif
    TYPE(BSET)          :: BS
    TYPE(CRDS)          :: GM
    TYPE(TIME)          :: Perf
    TYPE(INT_VECT)      :: AiFlg,ZiFlg,ColPt,BlkPt
    TYPE(ARGMT)         :: Args
    INTEGER             :: I,J,Q,R,IDex,JDex,ZDex,ZBlk,NIJ, &
                           n,ni,msiz,strtai,stopai,strtaj,stopaj, &
                           strtzi,stopzi,nj,strtzj,stopzj,jcol,k,kdex, &
                           aiblk,ajblk,zjblk,m,ziblk,icol,zrowpt,zcolpt, &
                           zblkpt,NewBloks,EndBloks,IRow,JRow,ZBlksPreFilter,ZBlksPostFilter
    TYPE(DBL_VECT)      :: Blk1,Blk2
    TYPE(DBL_RNK2)      :: P,DA
    REAL(DOUBLE)        :: Op,Mx0,B2Norm
    REAL(DOUBLE)        :: TrixThresh,DstncThresh,IRowX,IRowY,IRowZ
    TYPE(AtomPair)      :: Pair
    REAL(DOUBLE), &
         EXTERNAL       :: DDOT
    LOGICAL             :: TEST_AINV
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
    CHARACTER(LEN=8),&
         PARAMETER  :: Prog='BlokAInv'
    INTEGER :: II
!------------------------------------------------------------------------------------------------
#ifdef FIND_CONDA
    ! Useful if AInv is behaving badly.  If A is singular, you
    ! will have problems.
    CALL New(B,(/NBasF,NBasF/))
    CALL New(C,(/NBasF,NBasF/))
    CALL SetEq(B,A)
    CALL SetDSYEVWork(NBasF)
    CALL FunkOnSqMat(NBasF,Inverse,B%D,C%D,PrintCond_O=.TRUE.)
    CALL Delete(B)
    CALL Delete(C)
    CALL UnSetDSYEVWork()
#endif
    IF(.NOT.AllocQ(Z%Alloc)) &
    CALL New(Z)
#ifdef USE_METIS
    CALL MetisReorder(A)
#endif
    ! Set global workspace for FunkOnSqMat
    CALL SetDSYEVWork(MaxBlkSize)
    ! Allocate intermediate blocks
    IF(.NOT.AllocQ(Blk1%Alloc)) CALL New(Blk1,MaxBlkSize*MaxBlkSize)
    IF(.NOT.AllocQ(Blk2%Alloc)) CALL New(Blk2,MaxBlkSize*MaxBlkSize)
    ! Allocate diagonal "pivot" blocks
    IF(.NOT.AllocQ(P%Alloc)) CALL New(P,(/MaxBlkSize*MaxBlkSize,NAtoms/))
    ! Allocate coloumn flags
    IF(.NOT.AllocQ(AiFlg%Alloc)) CALL New(AiFlg,NAtoms); AiFlg%I=0
    IF(.NOT.AllocQ(ZiFlg%Alloc)) CALL New(ZiFlg,NAtoms); ZiFlg%I=0
    ! Allocate temporaries for symbolic Z_I
    IF(.NOT.AllocQ(ColPt%Alloc)) CALL New(ColPt,MaxBlks)
    IF(.NOT.AllocQ(BlkPt%Alloc)) CALL New(BlkPt,MaxBlks)
    ! Start with the identity; Z=I
    CALL SetToI(Z)
    ! Index for new Z bloks
    ZBlk=Z%NNon0+1
    ! Main outer loop down rows
    DO IRow=1,Natoms !!!!! DANGER !!!!! mixing Natoms and GM%Natms
      IF(IRow>GM%Natms) THEN
          II=Z%RowPt%I(GM%Natms+1)
          Z%RowPt%I(IRow+1)=II
        GO TO 200 !!!! Quick fix, temporary, good for specific cases
      ENDIF
#ifdef SPATIAL_THRESHOLDING
       ! Set IRow coordinates for distance based screening
       IRowX=GM%Carts%D(1,IRow)
       IRowY=GM%Carts%D(2,IRow)
       IRowZ=GM%Carts%D(3,IRow)
#endif
       NI=BSiz%I(IRow)
       StrtAI=A%RowPt%I(IRow);StopAI=A%RowPt%I(IRow+1)-1
       StrtZI=Z%RowPt%I(IRow);StopZI=Z%RowPt%I(IRow+1)-1
       DO I=StrtAI,StopAI
          AiFlg%I(A%ColPt%I(I))=I
       ENDDO
       ! Store symbolic structure of Z_I
       ZDex=StrtZI
       DO I=StrtZI,StopZI
          IDex=Z%ColPt%I(I)
          ZiFlg%I(IDex)=ZDex
          ColPt%I(ZDex)=IDex
          BlkPt%I(ZDex)=Z%BlkPt%I(I)
          ZDex=ZDex+1
       ENDDO
       ! Inner loop down rows ...
       DO JRow=1,IRow-1
#ifdef SPATIAL_THRESHOLDING
          IF(((IRowX-GM%Carts%D(1,JRow))**2+ &
               (IRowY-GM%Carts%D(2,JRow))**2+ &
               (IRowZ-GM%Carts%D(3,JRow))**2)<DstncThresh)THEN
#endif
             NJ=BSiz%I(JRow)
             StrtAJ=A%RowPt%I(JRow);StopAJ=A%RowPt%I(JRow+1)-1
             StrtZJ=Z%RowPt%I(JRow);StopZJ=Z%RowPt%I(JRow+1)-1
             ! Blk1=P^(j-1)_i=[A^t_j].[Z^(j-1)_i];
             ! Going down rows over N: (NxNJ)^T.(NxNI)
             NIJ=NI*NJ
             Blk1%D(1:NIJ)=Zero
             DO J=StrtAJ,StopAJ
                JDex=A%ColPt%I(J)
#ifdef SPATIAL_THRESHOLDING
                IF(((IRowX-GM%Carts%D(1,JDex))**2+ &
                     (IRowY-GM%Carts%D(2,JDex))**2+ &
                     (IRowZ-GM%Carts%D(3,JDex))**2)<DstncThresh)THEN
#endif
                   IDex=ZiFlg%I(JDex)
                   IF(IDex/=0)THEN
                      ZiBlk=BlkPt%I(IDex)
                      AjBlk=A%BlkPt%I(J)
                      M=BSiz%I(JDex)
                      CALL DGEMM_NN(NJ,M,NI,One,A%MTrix%D(AjBlk),Z%MTrix%D(ZiBlk),Blk1%D)
                      Perf%FLOP=Perf%FLOP+DBLE(NIJ*M)
                   ENDIF
#ifdef SPATIAL_THRESHOLDING
                ENDIF
#endif
             ENDDO
             ! Blk2=[P^(j-1)_j]^(-1).[P^(j-1)_i]
             CALL DGEMM_NNc(NJ,NJ,NI,One,Zero,P%D(:,JRow),Blk1%D,Blk2%D)
             Perf%FLOP=Perf%FLOP+DBLE(NIJ*NJ)
             !  Check the magintude of Blk2.  Update Z_I only if Blk2 is "large" enough.
             B2Norm=SQRT(DDOT(NI*NJ,Blk2%D,1,Blk2%D,1))
             Perf%FLOP=Perf%FLOP+DBLE(NIJ)
             IF(B2Norm>TrixThresh*1.D-1)THEN
                ! Z^j_i=Z^(j-1)_i-[Z^(j-1)_j].{[P^(j-1)_j]^(-1).[P^(j-1)_i]}
                ! Update going down rows:(NxNI)=(NxNI)+(NxNJ).(NJxNI)
                DO JDex=StrtZJ,StopZJ
                   JCol=ColPt%I(JDex)
#ifdef SPATIAL_THRESHOLDING
                   IF(((IRowX-GM%Carts%D(1,JCol))**2+ &
                       (IRowY-GM%Carts%D(2,JCol))**2+ &
                       (IRowZ-GM%Carts%D(3,JCol))**2)<DstncThresh)THEN
#endif
                      ZjBlk=BlkPt%I(JDex)
                      IDex =ZiFlg%I(JCol)
                      M=BSiz%I(JCol)
                      Perf%FLOP=Perf%FLOP+DBLE(NIJ*M)
                      IF(IDex/=0)THEN
                         ZiBlk=BlkPt%I(IDex)
                         CALL DGEMM_NNc(M,NJ,NI,-One,One,Z%MTrix%D(ZjBlk),Blk2%D,Z%MTrix%D(ZiBlk))
                      ELSE
                         ZiBlk=ZBlk
                         CALL DGEMM_NNc(M,NJ,NI,-One,Zero,Z%MTrix%D(ZjBlk),Blk2%D,Z%MTrix%D(ZiBlk))
                         ZiFlg%I(JCol)=ZDex
                         ColPt%I(ZDex)=JCol
                         BlkPt%I(ZDex)=ZiBlk
                         ZDex=ZDex+1
                         ZBlk=ZBlk+M*NI
                      ENDIF
#ifdef SPATIAL_THRESHOLDING
                   ENDIF
#endif
                ENDDO
#ifdef SPATIAL_THRESHOLDING
             ENDIF
#endif
          ENDIF
       ENDDO ! end inner loop over JRow
200    CONTINUE
       ! Reup symbolic structure of Z
       IF(IRow>1)THEN
          NewBloks=ZDex-StopZI-1
          EndBloks=Z%NBlks-StopZJ-1
          Z%RowPt%I(IRow+1:NAtoms+1)=Z%RowPt%I(IRow+1:NAtoms+1)+NewBloks
          ColPt%I(ZDex:ZDex+EndBloks-1)=Z%ColPt%I(StopZI+1:Z%NBlks)
          BlkPt%I(ZDex:ZDex+EndBloks-1)=Z%BlkPt%I(StopZI+1:Z%NBlks)
          Z%NBlks=Z%NBlks+NewBloks
          Z%NNon0=ZBlk
          Z%ColPt%I(StrtZI:Z%NBlks)=ColPt%I(StrtZI:Z%NBlks)
          Z%BlkPt%I(StrtZI:Z%NBlks)=BlkPt%I(StrtZI:Z%NBlks)
       ENDIF
       !  Blk1=P^(i-1)_i=[A^t_i].[Z^(i-1)_i]; Going down rows:(NxNI)^T.(NxNI)
       Blk1%D=Zero
       StrtZI=Z%RowPt%I(IRow);StopZI=Z%RowPt%I(IRow+1)-1
       DO IDex=StrtZI,StopZI
          ICol=Z%ColPt%I(IDex)
          KDex=AiFlg%I(ICol)
          IF(KDex/=0)THEN
             ZiBlk=Z%BlkPt%I(IDex)
             AiBlk=A%BlkPt%I(KDex)
             M=BSiz%I(ICol)
             CALL DGEMM_NN(NI,M,NI,One,A%MTrix%D(AiBlk),Z%MTrix%D(ZiBlk),Blk1%D)
             Perf%FLOP=Perf%FLOP+DBLE(NI*M*NI)
          ENDIF
       ENDDO
       ! P(I)=[P^(i-1)_i]^(-1)
       CALL FunkOnSqMat(NI,Inverse,Blk1%D,P%D(:,IRow))
       ! Estimated performance; 2 DGEMMS+1 DSYEV
       Perf%FLOP=Perf%FLOP+DBLE((2+6)*NI**3)
       ! Reset flags for column flags
       DO I=StrtZI,ZDex-1
          ZiFlg%I(ColPt%I(I))=0
       ENDDO
       DO J=StrtAI,StopAI
          AiFlg%I(A%ColPt%I(J))=0
       ENDDO
    ENDDO ! end main loop over IRow
    !------------------------------------------------------------------------------------------
    ! Finishing touches on Z
    Z%NAtms=NAtoms
    Z%NBlks=ZDex-1
    Z%NNon0=ZBlk-1
    ! Free some memory
    CALL Delete(ColPt)
    CALL Delete(BlkPt)
    ! Compute dimensions of DiagD and allocate it
    DiagD%NBlks=NAtoms
    DiagD%NNon0=0
    DO I=1,NAtoms
       DiagD%NNon0=DiagD%NNon0+BSiz%I(I)**2
    ENDDO
    IF(.NOT.AllocQ(DiagD%Alloc)) CALL New(DiagD,(/NAtoms,DiagD%NBlks,DiagD%NNon0/))
    ! DiagD=P^(-1/2) in BCSR format
    DO I=1,NAtoms
       N=BSiz%I(I)
       CALL FunkOnSqMat(N,SqRoot,P%D(:,I),Blk1%D)
       !    Estimated performance; 2 DGEMMS+1 DSYEV
       Perf%FLOP=Perf%FLOP+DBLE((2+6)*NI**3)
       P%D(1:N*N,I)=Blk1%D(1:N*N)
    ENDDO
    CALL SetToI(DiagD,P)
    ! Free some more memory
    CALL Delete(P)
    CALL Delete(Blk1)
    CALL Delete(Blk2)
    CALL UnSetDSYEVWork()
    ! This is workspace for Z^t
    IF(.NOT.AllocQ(Zt%Alloc)) &
    CALL New(Zt)
    ! Symbolic transpose only, bloks in place
    CALL XPose(Z)
    ZBlksPreFilter=Z%NBlks
    ! Final Z=P^(-1/2).Z
    CALL Multiply(Z,DiagD,Zt)
    !  CALL PlotDecay(Zt,GM,'Z')
    CALL Filter(Z,Zt)
    ZBlksPostFilter=Z%NBlks
    ! Full transpose
    CALL XPose(Z,Zt)

    ! Account for multiplies AND adds in DGEMMs
    Perf%FLOP=Perf%Flop*Two
    CALL Elapsed_TIME(Perf,'Accum')


    CALL Delete(DiagD)

  END SUBROUTINE BlockedAInv
END MODULE AInv


