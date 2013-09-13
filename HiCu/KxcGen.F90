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
!    COMPUTE THE EXCHANGE CORRELATION MATRIX $K_{xc}$ IN O(N)
!    USING HIERARCHICAL CUBATURE INVOLVING K-D BINARY TREES
!    Author: Matt Challacombe
!-------------------------------------------------------------------------------
MODULE KxcGen
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalObjects
  USE ProcessControl
  USE Indexing
  USE McMurchie
  USE InOut
  USE Thresholding
  USE HiCuThresholds
  USE AtomPairs
  USE BraBloks
  USE CubeTree
  USE TreeWalk
#ifdef PARALLEL
  USE FastMatrices
#endif
  IMPLICIT NONE
!-------------------------------------------------------------------------------!
  CONTAINS !
!===============================================================================

!===============================================================================
  SUBROUTINE MakeKxc(Kxc,CubeRoot)
#ifdef PARALLEL
    TYPE(FastMat),POINTER  :: Kxc
#else
    TYPE(BCSR)             :: Kxc
#endif
    TYPE(CubeNode),POINTER :: CubeRoot

    TYPE(DBL_RNK2)         :: Temp
    TYPE(AtomPair)         :: Pair
    INTEGER                :: AtA,AtB
    INTEGER                :: JP,K,NA,NB,NAB,P,Q,R,I,J,iSMat
    INTEGER                   :: NCA,NCB
    REAL(DOUBLE)              :: Fac
    REAL(DOUBLE),DIMENSION(3) :: A,B
!-------------------------------------------------------------------------------
!   Initialize the matrix and associated indecies
#ifdef PARALLEL
#else
    P=1;
    R=1;
    Kxc%RowPt%I(1)=1
    CALL SetEq(Kxc%MTrix,Zero)
    Kxc%NAtms=NAtoms
#endif
!   Loop over atom pairs
    DO AtA=1,NAtoms
      DO AtB=1,NAtoms
        IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
          NAB = Pair%NA*Pair%NB*Kxc%NSMat
          IF(AtB<=AtA)THEN
#ifdef PARALLEL
             Fac=1D0
             IF(AtA.EQ.AtB)Fac=0.5D0
#endif
!         Compute only the lower triangle of symmetric Kxc
            A = Pair%A
            B = Pair%B
            DO NCA = 1,CS_OUT%NCells
              Pair%A = A+CS_OUT%CellCarts%D(:,NCA)
              DO NCB = 1,CS_OUT%NCells
                Pair%B = B+CS_OUT%CellCarts%D(:,NCB)
                Pair%AB2  = (Pair%A(1)-Pair%B(1))**2 &
                            + (Pair%A(2)-Pair%B(2))**2 &
                            + (Pair%A(3)-Pair%B(3))**2
                IF(TestAtomPair(Pair,CubeRoot%Box)) THEN
#ifdef PARALLEL
                  CALL AddFASTMATBlok(Kxc,AtA,AtB,Pair%NA,Pair%NB*Kxc%NSMat,Fac*KxcBlock(Pair,CubeRoot,Kxc%NSMat))
#else
                  Kxc%MTrix%D(R:R+NAB-1)=Kxc%MTrix%D(R:R+NAB-1)+KxcBlock(Pair,CubeRoot,Kxc%NSMat)
#endif
                ENDIF
              ENDDO
            ENDDO
!vw#ifndef PARALLEL
          ENDIF
!vw#endif

!! the logic here strictly follows MakeS
#ifdef PARALLEL
#else
          Kxc%ColPt%I(P)=AtB
          Kxc%BlkPt%I(P)=R
          R=R+NAB
          P=P+1
!! deleted DBCSR part
          Kxc%RowPt%I(AtA+1)=P
          IF(R>MaxNon0*Kxc%NSMat.OR.P>MaxBlks) THEN
             CALL Halt(' BCSR dimensions blown in Kxc ')
          ENDIF
#endif
        ENDIF
      ENDDO
    ENDDO

#ifdef PARALLEL
#else
    Kxc%NBlks=P-1
    Kxc%NNon0=R-1
#endif


#ifdef PARALLEL
#else
!   Fill the upper triangle of Kxc
    DO I=1,NAtoms
      DO JP=Kxc%RowPt%I(I),Kxc%RowPt%I(I+1)-1
        J=Kxc%ColPt%I(JP)
        IF(I>J)THEN
           DO K=Kxc%RowPt%I(J),Kxc%RowPt%I(J+1)-1
              IF(Kxc%ColPt%I(K)==I)THEN
                 Q=Kxc%BlkPt%I(K)
                 EXIT
              ENDIF
           ENDDO
           P=Kxc%BlkPt%I(JP)
           NA=BS%BFKnd%I(GM%AtTyp%I(I))
           NB=BS%BFKnd%I(GM%AtTyp%I(J))
           NAB=NA*NB
           DO iSMat=1,Kxc%NSMat
              CALL XPose(NA,NB,Kxc%MTrix%D(P:P+NAB-1),Kxc%MTrix%D(Q:Q+NAB-1))
              P=P+NAB
              Q=Q+NAB
           ENDDO
        ENDIF
      ENDDO
    ENDDO
#endif
  END SUBROUTINE MakeKxc
!===============================================================================

!===============================================================================
#ifdef PARALLEL
     FUNCTION KxcBlock(Pair,CubeRoot,NSMat) RESULT(Kblk)
#else
     FUNCTION KxcBlock(Pair,CubeRoot,NSMat) RESULT(Kvct)
#endif
       TYPE(AtomPair)                           :: Pair
       TYPE(CubeNode), POINTER                  :: CubeRoot
       INTEGER                                  :: NSMat

       REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB*NSMat)  :: Kblk
#ifdef PARALLEL
#else
       REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB*NSMat)  :: Kvct
#endif

       REAL(DOUBLE)                             :: ZetaA,ZetaB,EtaAB,EtaIn,    &
                                                   XiAB,ExpAB,CA,CB,CC,Ov,     &
                                                   PAx,PAy,PAz,PBx,PBy,PBz,    &
                                                   MDx,MDxy,MDxyz,Amp2,MaxAmp
       INTEGER                                  :: KA,KB,CFA,CFB,PFA,PFB,      &
                                                   IndexA,IndexB,              &
                                                   StartLA,StartLB,            &
                                                   StopLA,StopLB,LenHG
       INTEGER                                  :: I,J,MaxLA,MaxLB,IA,IB,    &
                                                   LMNA,LMNB,LA,LB,MA,MB,    &
                                                   NA,NB,LAB,MAB,NAB,LMN,EllA,EllB
real(double):: Pextent_Old
!-------------------------------------------------------------------------------
       KBlk=Zero
       KA=Pair%KA
       KB=Pair%KB
!-------------------------------------------------------------------------------
       Prim%A=Pair%A
       Prim%B=Pair%B
       Prim%AB2=Pair%AB2
       Prim%KA=Pair%KA
       Prim%KB=Pair%KB
!-------------------------------------------------------------------------------
       IndexA=0
       DO CFA=1,BS%NCFnc%I(KA)
       DO CFB=1,BS%NCFnc%I(KB)
          IndexA  = CFBlokDex(BS,CFA,KA)
          IndexB  = CFBlokDex(BS,CFB,KB)
          StartLA = BS%LStrt%I(CFA,KA)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLA  = BS%LStop%I(CFA,KA)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLA   = BS%ASymm%I(2,CFA,KA)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
!-------------------------------------------------------------------------------
          Prim%CFA=CFA
          Prim%CFB=CFB
          Prim%Ell=MaxLA+MaxLB
!-------------------------------------------------------------------------------
          DO PFA=1,BS%NPFnc%I(CFA,KA)
          DO PFB=1,BS%NPFnc%I(CFB,KB)
!-------------------------------------------------------------------------------
             Prim%ZA=BS%Expnt%D(PFA,CFA,KA)
             Prim%ZB=BS%Expnt%D(PFB,CFB,KB)
             Prim%Zeta=Prim%ZA+Prim%ZB
             Prim%Xi=Prim%ZA*Prim%ZB/Prim%Zeta
             IF(TestPrimPair(Prim%Xi,Prim%AB2))THEN
                Prim%PFA=PFA
                Prim%PFB=PFB
!               Set primitive bra blok, find its extent
                PExtent=SetBraBlok(Prim,BS,Tau_O=TauRho,ExtraEll_O=1)
                PBox%BndBox(:,1)=Prim%P
                PBox%BndBox(:,2)=Prim%P
                PBox=ExpandBox(PBox,PExtent)
!               Quick check to see if primitive touches the grid
                IF(PExtent>Zero.AND.(.NOT.BoxOutSideBox(PBox,CubeRoot%Box)))THEN
!                  Recompute extent now on case by case basis
!                  Perhaps using a better extent
                   PExtent=Zero
                   IA = IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                         EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)
                         LenHG=LHGTF(EllA+EllB)
!                        Using Extent1 here halves cpu time of Kxc build
!                        Numbers for loose seem somewhat less robust though
!                        Need alot more testing (leave for now unitll forces are 100%).
                         PExtent=MAX(PExtent,Extent(EllA+EllB,Prim%Zeta, &
                                     HGBra%D(1:LenHG,IA,IB),TauRho,ExtraEll_O=1))
                      ENDDO
                   ENDDO
!                  Set BBox for this primitive
                   PBox%BndBox(:,1)=Prim%P
                   PBox%BndBox(:,2)=Prim%P
                   PBox=ExpandBox(PBox,PExtent)
!                  Walk the walk
                   Ket=Zero
                   CALL KxcWalk(CubeRoot)
!                  Contract <Bra|Ket> bloks to compute matrix elements of Kxc
                   IA = IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                         EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)
                         DO LMN=1,LHGTF(EllA+EllB)

                            IF(NSMat.EQ.1)THEN
                               Kblk(IA,IB)=Kblk(IA,IB)+HGBra%D(LMN,IA,IB)*Ket(LMN)
                            ELSE
                               Kblk(IA,IB)=Kblk(IA,IB)+HGBra%D(LMN,IA,IB)*Ket(LMN)
                               Kblk(IA,IB+Pair%NB)=Kblk(IA,IB+Pair%NB)+HGBra%D(LMN,IA,IB)*Ket(LMN+LHGTF(Prim%Ell))
                            ENDIF

                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
!-------------------------------------------------------------------------------
             ENDIF
          ENDDO
          ENDDO
       ENDDO
       ENDDO
#ifdef PARALLEL
#else
       Kvct=BlockToVect(Pair%NA,Pair%NB*NSMat,Kblk)
#endif
     END FUNCTION KxcBlock
!===============================================================================

!===============================================================================
END MODULE KxcGen
