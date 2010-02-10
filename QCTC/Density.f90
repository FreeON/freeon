!------------------------------------------------------------------------------
!    This code is part of the FreeON suite of programs for linear scaling
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
!    to return derivative works to the FreeON group for review, and possible
!    dissemination in future releases.
!------------------------------------------------------------------------------
MODULE Density
  USE Derivedtypes
  USE GlobalScalars
  USE GlobalObjects
  USE ProcessControl
  USE QCTCThresholds
  USE Order
  USE Clock
  USE BraBloks
  USE PBC
  IMPLICIT NONE
  !
  TYPE HGLL
     INTEGER                   :: Ell      ! The highest angular symmetry in this cluster
     INTEGER                   :: LNum     ! The naive id # of this distribution in the linked list
     REAL(DOUBLE)              :: Zeta     ! The true zeta for the distributions
     REAL(DOUBLE)              :: Schw     ! The value to use in the Schartz inequality: ( | )^1/2
     REAL(DOUBLE),DIMENSION(3) :: Cent     ! This is the center of the distribution
     TYPE(DBL_VECT)            :: Coef     ! Coefficients of the HGTF density
     TYPE(HGLL),POINTER        :: Next     ! Next link in the chain
  END TYPE HGLL
  !
  TYPE HGLLP
     TYPE(HGLL),POINTER :: L
  END TYPE HGLLP
  !
  TYPE(HGRho)       :: Rho
  TYPE(CMPoles)     :: RhoPoles

  REAL(DOUBLE) :: Density_Time_Start,Density_Time
  !
CONTAINS
  !
  SUBROUTINE MakeRhoList(GM,BS,P,NLink,RhoHead,Prog,NoWrap_O)
    !
    TYPE(CRDS)                      :: GM
    TYPE(BSET)                      :: BS
    TYPE(BCSR)                      :: P
    TYPE(HGLL), POINTER             :: RhoHead,HGLink
    TYPE(AtomPair)                  :: Pair
    TYPE(DBL_VECT)                  :: Psv
    REAL(DOUBLE),DIMENSION(3)       :: B
    INTEGER                         :: AtA,AtB,NC,PBeg,PEnd,BlkP,RowP,NN,SpinM,NLink,NNaive
    LOGICAL, OPTIONAL               :: NoWrap_O
    LOGICAL                         :: NoWrap,TestTru
    CHARACTER(LEN=*)                :: Prog
    CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Mssg
    !
    IF(PRESENT(NoWrap_O))THEN
       NoWrap=NoWrap_O
    ELSE
       NoWrap=.FALSE.
    ENDIF

    CALL NewBraBlok(BS)
    !
    NLink=0
    HGLink=>RhoHead
    !
    SELECT CASE(P%NSMat)
    CASE(1)
       SpinM=1
    CASE(2)
       SpinM=3
    CASE(4)
       SpinM=4
    CASE DEFAULT
       CALL Halt('MakeRho: NSMat not valid!')
    END SELECT
    !
    CALL New(Psv,SpinM*MaxBlkSize**2)
    !
    DO AtA=1,NAtoms
       Pbeg = P%RowPt%I(AtA)
       Pend = P%RowPt%I(AtA+1)-1
       DO RowP=Pbeg,Pend
          AtB=P%ColPt%I(RowP)
          IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
!          TestTru=SetAtomPair(GM,BS,AtA,AtB,Pair)
!          IF(.TRUE.)THEN
             ! Sum this part of the DM
             BlkP=P%BlkPt%I(RowP)
             NN=BSiz%I(AtA)*BSiz%I(AtB)
             SELECT CASE(P%NSMat)
             CASE(1)
                !We need to copy the matrix!
                CALL DCOPY(NN,P%MTrix%D(BlkP),1,Psv%D(1),1)
             CASE(2)
                !We copy the first matrix!
                CALL DCOPY(NN,P%MTrix%D(BlkP),1,Psv%D(1),1)
                !(Pa+Pb)/2->Ptot
                CALL DAXPY(NN,1D0,P%MTrix%D(BlkP+NN),1,Psv%D(1),1)
                ! Scale the total density matrix.
                CALL DSCAL(NN,0.5D0,Psv%D(1),1)
                !Pa
                CALL DCOPY(NN,P%MTrix%D(BlkP),1,Psv%D(NN+1),1)
                !Pb
                CALL DCOPY(NN,P%MTrix%D(BlkP+NN),1,Psv%D(2*NN+1),1)
             CASE(4)
                !We copy the first matrix!
                CALL DCOPY(NN,P%MTrix%D(BlkP),1,Psv%D(1),1)
                !(Pa+Pb)/2->Ptot
                CALL DAXPY(NN,1D0,P%MTrix%D(BlkP+3*NN),1,Psv%D(1),1)
                ! Scale the total density matrix.
                CALL DSCAL(NN,0.5D0,Psv%D(1),1)
             CASE DEFAULT
                CALL Halt(' MakeRho: P%NSMat doesnt have an expected value! ')
             END SELECT
             B=Pair%B
             DO NC=1,GM%OvCells%NCells
                Pair%B=B+GM%OvCells%CellCarts%D(:,NC)
                Pair%AB2=(Pair%A(1)-Pair%B(1))**2 &
                       + (Pair%A(2)-Pair%B(2))**2 &
                       + (Pair%A(3)-Pair%B(3))**2
                IF(TestAtomPair(Pair))THEN
                   CALL RhoPop(GM,BS,GM%InCells,NoWrap,SpinM,Psv%D(1),Pair,HGLink,NLink)
                ENDIF
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    !
    CALL Delete(Psv)
    ! Now eliminate redundancies in the density    !
    NNaive=NLink

    CALL RhoEcon(RhoHead,NLink)

    CALL MondoLog(DEBUG_MAXIMUM,Prog,'RhoEcon = '//TRIM(FltToShrtChar(DBLE(NNaive)/DBLE(NLink))) &
                                     //', Gaussians = '//TRIM(IntToChar(NLink)))

    ! Next, calculate the Schwartz inequality and prune small links ???

    ! Free bra block memory, since we may need it again
    CALL DeleteBraBlok()

  END SUBROUTINE MakeRhoList
  !
  SUBROUTINE Collate(GM,RhoHead,R,Prog,P,NLink)
    TYPE(CRDS)                      :: GM
    TYPE(HGLL),POINTER              :: RhoHead,HGLink    !
    TYPE(CMPoles)                   :: P
    TYPE(HGRho)                     :: R
    INTEGER                         :: NLink
    TYPE(HGLLP),DIMENSION(NLink)    :: LinkArray
    REAL(DOUBLE),DIMENSION(NLink)   :: RealArray
    INTEGER,DIMENSION(NLink)        :: ListArray
    INTEGER                         :: N,NExpt,NCoef,Status
    INTEGER                         :: I,J,K,IAdd,JAdd,NQ,NR,EllZ,LenZ,LenCo
    INTEGER,PARAMETER               :: MaxExp=2000
    INTEGER,DIMENSION(MaxExp)       :: NZ,NL,El,IZ
    REAL(DOUBLE), DIMENSION(MaxExp) :: Ze
    REAL(DOUBLE)                    :: Zeta,ZZ,ChTot,ChNu,ChEl,RX,RY,RZ
    CHARACTER(LEN=*)                :: Prog
    CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Mssg

    REAL(DOUBLE),DIMENSION(3) :: TMP1
   !
    HGLink=>RhoHead%Next
    !
    NExpt=1
    NL(1)=1
    NZ(1)=1
    El(1)=HGLink%Ell
    Ze(1)=HGLink%Zeta
    !
    HGLink=>HGLink%Next

!!    WRITE(*,*)' Fucked in collate (SEE BELOW COMMENTED OUT IF!!!!!!!!!!!! '

    DO WHILE(ASSOCIATED(HGLink))
       DO I=1,NExpt
          IF(ABS(Ze(I)-HGLink%Zeta)<1D-12)THEN
!          IF(.FALSE.)THEN
             NZ(I)=NZ(I)+1
             NL(I)=NL(I)+1
             El(I)=MAX(El(I),HGLink%Ell)
             GOTO 101
          ENDIF
       ENDDO
       NExpt=NExpt+1
       IF(NExpt>MaxExp)CALL Halt(' To many exponents in IntegrateNCollate ')
       NZ(NExpt)=1
       El(NExpt)=HGLink%Ell
       Ze(NExpt)=HGLink%Zeta
101    CONTINUE
       HGLink=>HGLink%Next
    ENDDO
    !
    IF(SUM(NZ(1:NExpt)).NE.NLink) &
       CALL Halt(' Counting error in IntegrateNCollate ')
    !
    NCoef=0
    DO I=1,NExpt
       IZ(I)=I
       NCoef=NCoef+NZ(I)*LHGTF(El(I))
    ENDDO
    !
    CALL New(R,(/NExpt,NLink,NCoef,1/))
    !
    CALL DblIntSort77(NExpt,Ze,IZ,2)
    !
    I=0
    HGLink=>RhoHead%Next
    DO WHILE(ASSOCIATED(HGLink))
       I=I+1
       ListArray(I)=I
       LinkArray(I)%L=>HGLink
       RealArray(I)=HGLink%Zeta
       HGLink=>HGLink%Next
    ENDDO
    !
    CALL DblIntSort77(NLink,RealArray,ListArray,2)
    DO I=1,NLink
       HGLink=>LinkArray(ListArray(I))%L
!       WRITE(33,23)I,HGLink%Ell,HGLink%Zeta
23     FORMAT(I10,(I3,', '),1(F14.5,', '))
    ENDDO
    !
    NQ=0
    NR=0
    DO I=1,NExpt
       R%Expt%D(I)=Ze(I)
       R%NQ%I(I)=NZ(IZ(I))
       R%Lndx%I(I)=El(IZ(I))
       R%OffQ%I(I)=NQ
       R%OffR%I(I)=NR
       EllZ=El(IZ(I))
       LenZ=LHGTF(EllZ)
       DO J=1,NZ(IZ(I))
          HGLink=>LinkArray(ListArray(Nq+1))%L
!!          WRITE(*,22)I,J,El(IZ(I)),HGLink%Ell,Ze(I),HGLink%Zeta
!!22        FORMAT(4(I3,', '),2(F14.5,', '))
          R%Qx%D(Nq+1)=HGLink%Cent(1)
          R%Qy%D(Nq+1)=HGLink%Cent(2)
          R%Qz%D(Nq+1)=HGLink%Cent(3)
          IF(HGLink%Ell==EllZ)THEN
             R%Co%D(NR+1:NR+LenZ)=HGLink%Coef%D
          ELSE
             R%Co%D(NR+1:NR+LenZ)=Zero
             LenCo=LHGTF(HGLink%Ell)
             R%Co%D(NR+1:NR+LenCo)=HGLink%Coef%D(1:LenCo)
          ENDIF
          NQ=NQ+1
          NR=NR+LenZ
       ENDDO
    ENDDO
    !
    ChEl=0D0
    ChNu=0D0
    CALL New(P)

    DO I=1,R%NExpt
       EllZ=R%Lndx%I(I)
       LenZ=LHGTF(EllZ)
       Zeta=R%Expt%D(I)
       ZZ=(Pi/R%Expt%D(I))**1.5D0
       DO J=1,R%NQ%I(I)
          IAdd=R%OffQ%I(I)+J
          JAdd=R%OffR%I(I)+(J-1)*LenZ+1
          RX=Rho%Qx%D(IAdd)
          RY=Rho%Qy%D(IAdd)
          RZ=Rho%Qz%D(IAdd)
          ! Note factors of two to account for closed shell book-keeping
          P%MPole=P%MPole+Two*R%Co%D(JAdd)*ZZ
          P%DPole%D=P%DPole%D+Two*DPole(EllZ,ZZ,RX,RY,RZ,Rho%Co%D(JAdd:JAdd+LenZ-1))
          P%QPole%D=P%QPole%D+Two*QPole(EllZ,Zeta,ZZ,RX,RY,RZ,Rho%Co%D(JAdd:JAdd+LenZ-1))
       ENDDO
    ENDDO

!!$    P%MPole=Zero
!!$    DO I=1,R%NExpt-1
!!$       EllZ=R%Lndx%I(I)
!!$       LenZ=LHGTF(EllZ)
!!$       Zeta=R%Expt%D(I)
!!$       ZZ=(Pi/R%Expt%D(I))**1.5D0
!!$       DO J=1,R%NQ%I(I)
!!$          IAdd=R%OffQ%I(I)+J
!!$          JAdd=R%OffR%I(I)+(J-1)*LenZ+1
!!$           ! Note factors of two to account for closed shell book-keeping
!!$          P%MPole=P%MPole+Two*R%Co%D(JAdd)*ZZ
!!$       ENDDO
!!$    ENDDO
!!$    WRITE(*,*)' Qel  = ',P%MPole
!!$    P%MPole=Zero
!!$    DO I=R%NExpt,R%NExpt
!!$       EllZ=R%Lndx%I(I)
!!$       LenZ=LHGTF(EllZ)
!!$       Zeta=R%Expt%D(I)
!!$       ZZ=(Pi/R%Expt%D(I))**1.5D0
!!$       DO J=1,R%NQ%I(I)
!!$          IAdd=R%OffQ%I(I)+J
!!$          JAdd=R%OffR%I(I)+(J-1)*LenZ+1
!!$           ! Note factors of two to account for closed shell book-keeping
!!$          P%MPole=P%MPole+Two*R%Co%D(JAdd)*ZZ
!!$       ENDDO
!!$    ENDDO
!!$    WRITE(*,*)' Qnuc  = ',P%MPole
!!$    STOP


    Mssg='<q> = '//TRIM(DblToShrtChar(P%MPole))        &
         //', <r> = ('//TRIM(DblToShrtChar(P%DPole%D(1))) &
         //', '//TRIM(DblToShrtChar(P%DPole%D(2)))        &
         //', '//TRIM(DblToShrtChar(P%DPole%D(3)))        &
         //'), <r^2> = '//TRIM(DblToShrtChar(             &
         P%QPole%D(1)+P%QPole%D(2)+P%QPole%D(3)))


    CALL MondoLog(DEBUG_MAXIMUM, Prog, Mssg, "DensityBuild")


  END SUBROUTINE Collate
  !
  SUBROUTINE DeleteHGLL(RhoHead)
    TYPE(HGLL),POINTER  :: RhoHead,HGLink,RMLink
    HGLink=>RhoHead%Next
    DEALLOCATE(RhoHead)
    NULLIFY(RhoHead)
    DO WHILE(ASSOCIATED(HGLink))
       RMLink=>HGLink
       CALL Delete(RMLink%Coef)
       DEALLOCATE(RMLink)
       HGLink=>HGLink%Next
    ENDDO
  END SUBROUTINE DeleteHGLL
  !
  SUBROUTINE AddNukes(GM,RhoHead,NoWrap_O)
    TYPE(CRDS)          :: GM
    TYPE(HGLL),POINTER  :: RhoHead,HGLink    !
    REAL(DOUBLE)        :: GSpike
    INTEGER             :: N,Status,NLink
    LOGICAL, OPTIONAL         :: NoWrap_O
    LOGICAL                   :: NoWrap
    TYPE(PrimPair)                          :: Prim
    !
    IF(PRESENT(NoWrap_O))THEN
       NoWrap=NoWrap_O
    ELSE
       NoWrap=.FALSE.
    ENDIF

    ! Find the last link
    HGLink=>RhoHead%Next
    DO WHILE(ASSOCIATED(HGLink%Next))
       HGLink=>HGLink%Next
    ENDDO
    NLink=HGLink%LNum
    GSpike=(NuclearExpnt/Pi)**(ThreeHalves)
    ! Add nuclear charges to the density
    DO N=1,GM%NAtms
       ALLOCATE(HGLink%Next,STAT=Status)
       IF(Status/=SUCCEED)CALL Halt(' Link allocation failed in QCTC2/RhoBuild.F90:AddNukes ')
       HGLink=>HGLink%Next
       NLink=NLink+1
       HGLink%LNum=NLink
       HGLink%Ell=0
       HGLink%Zeta=NuclearExpnt

       Prim%P=GM%Carts%D(:,N)
       CALL PWrap(GM,Prim,.NOT.NoWrap)
       HGLink%Cent=Prim%Pw

       CALL New(HGLink%Coef,1)
       ! Note factor of 1/2
       HGLink%Coef%D(1)=-Half*GM%AtNum%D(N)*GSpike
    ENDDO
    NULLIFY(HGLink%Next)
  END SUBROUTINE AddNukes

  SUBROUTINE RhoEcon(RhoHead,NLink)
    TYPE(HGLL),POINTER                :: RhoHead,HGLink,HGClone
    INTEGER                           :: NLink,ILink
    TYPE(HGLLP),DIMENSION(NLink)      :: LinkArray1
    REAL(DOUBLE),DIMENSION(NLink)     :: RealArray1,RealArray2
    INTEGER,DIMENSION(NLink)          :: ListArray1,ListArray2
    INTEGER                           :: I,J,K,L,JLst,KLst,KLstM1,KLstP1
    INTEGER                           :: EllJ,EllK,LenJ,LenK,ISort
    REAL(DOUBLE)                      :: MinX,MaxX,MinY,MaxY,MinZ,MaxZ
    REAL(DOUBLE)                      :: ZK,ZJ,ZZJK,RJK2,DX,Bhatt
    REAL(DOUBLE),DIMENSION(3)         :: RJ,RK

    IF(NLink < 1) THEN
      CALL MondoLog(DEBUG_NONE, "RhoEcon", "called with NLink = "//TRIM(IntToChar(NLink)))
      CALL Halt("illegal value")
    ENDIF

    ILink=0
    HGLink=>RhoHead%Next

    MaxX=-1D20
    MinX=1D20
    MaxY=-1D20
    MinY=1D20
    MaxZ=-1D20
    MinZ=1D20

    DO WHILE(ASSOCIATED(HGLink))
       ILink=ILink+1
       ListArray1(ILink)=ILink
       LinkArray1(ILink)%L=>HGLink
       MaxX=MAX(MaxX,HGLink%Cent(1))
       MinX=MIN(MinX,HGLink%Cent(1))
       MaxY=MAX(MaxY,HGLink%Cent(2))
       MinY=MIN(MinY,HGLink%Cent(2))
       MaxZ=MAX(MaxZ,HGLink%Cent(3))
       MinZ=MIN(MinZ,HGLink%Cent(3))
       HGLink=>HGLink%Next
    ENDDO
    !
    IF(MaxX-MinX>=MAX(MaxY-MinY,MaxZ-MinZ))THEN
       ISort=1
    ELSEIF(MaxY-MinY>=MAX(MaxX-MinX,MaxZ-MinZ))THEN
       ISort=2
    ELSE
       ISort=3
    ENDIF
    !
    DO I=1,NLink
       RealArray1(I)=LinkArray1(I)%L%Cent(ISort)
    ENDDO
    !
    CALL DblIntSort77(NLink,RealArray1,ListArray1,2)
    !
    ! Look for redundancies in the (approximate) Bhattacharyya measure
    !
    J=1
    DO WHILE(J.LE.NLink)

       IF(.NOT.ASSOCIATED(LinkArray1(ListArray1(J))%L))THEN
          J=J+1
          CYCLE
       ENDIF

       ZJ=LinkArray1(ListArray1(J))%L%Zeta
       RJ=LinkArray1(ListArray1(J))%L%Cent

       L=0
       DO K=J+1,NLink
          IF(.NOT.ASSOCIATED(LinkArray1(ListArray1(K))%L))CYCLE
          DX=ABS(RJ(ISort)-LinkArray1(ListArray1(K))%L%Cent(ISort))
          IF(DX.GT.1D-6)EXIT
          ZK=LinkArray1(ListArray1(K))%L%Zeta
          RK=LinkArray1(ListArray1(K))%L%Cent
          ZZJK=ZK*ZJ/(ZK+ZJ)
          RJK2=DOT_PRODUCT(RJ-RK,RJ-RK)
          Bhatt=ABS(1D0-EXP(-0.25D0*RJK2*ZZJK)*SQRT(2D0*SQRT(ZK*ZJ)/(ZK+ZJ)))
          L=L+1
          RealArray1(L)=Bhatt
          ListArray2(L)=ListArray1(K)
       ENDDO

       IF(L>0)THEN
          CALL DblIntSort77(L,RealArray1,ListArray2,2)

          DO I=1,L
             IF(RealArray1(I).GT.1D-12)EXIT
             !
             JLst=ListArray1(J)
             KLst=ListArray2(I)

             IF(KLst==1.OR.KLst==NLink)CYCLE

             EllJ=LinkArray1(JLst)%L%Ell
             EllK=LinkArray1(KLst)%L%Ell
             LenJ=LHGTF(EllJ)
             LenK=LHGTF(EllK)
             !
             ! Add density K to J link
             !
             IF(EllJ>=EllK)THEN
                LinkArray1(JLst)%L%Coef%D(1:LenK)=LinkArray1(JLst)%L%Coef%D(1:LenK) &
                                                 +LinkArray1(KLst)%L%Coef%D(1:LenK)
             ELSE
                LinkArray1(KLst)%L%Coef%D(1:LenJ)=LinkArray1(JLst)%L%Coef%D(1:LenJ) &
                                                 +LinkArray1(KLst)%L%Coef%D(1:LenJ)
                CALL Delete(LinkArray1(JLst)%L%Coef)
                CALL New(LinkArray1(JLst)%L%Coef,LenK)
                LinkArray1(JLst)%L%Ell=EllK
                EllJ=EllK
                LinkArray1(JLst)%L%Coef%D(1:LenK)=LinkArray1(KLst)%L%Coef%D(1:LenK)
             ENDIF
             !
             ! Remove link K from the list
             !
             CALL Delete(LinkArray1(KLst)%L%Coef)
             DEALLOCATE(LinkArray1(KLst)%L)
             NULLIFY(LinkArray1(KLst)%L)
             !
             KLstM1=KLst-1
             IF(KLstM1<1)CYCLE
             DO WHILE(.NOT.ASSOCIATED(LinkArray1(KLstM1)%L))
                KLstM1=KLstM1-1
                IF(KLstM1<1)THEN
                   KLstM1=1
                   EXIT
                ENDIF
             ENDDO
             !
             KLstP1=KLst+1
             IF(KLstP1>NLink)CYCLE
             DO WHILE(.NOT.ASSOCIATED(LinkArray1(KLstP1)%L))
                KLstP1=KLstP1+1
                IF(KLstP1>NLink)THEN
                   KLstP1=NLink
                   EXIT
                ENDIF
             ENDDO

             IF(KLstM1==KLstP1)GOTO 101
             !
             LinkArray1(KLstM1)%L%Next=>LinkArray1(KLstP1)%L
             !
          ENDDO
       ENDIF
       J=J+1
    ENDDO
    !
    NLink=0
    HGLink=>RhoHead%Next
    DO WHILE(ASSOCIATED(HGLink))
       NLink=NLink+1
       HGLink=>HGLink%Next
    ENDDO
    RETURN
    !
101 CONTINUE
    WRITE(*,*)' KLst = ',KLst
    WRITE(*,*)' KLst-1 = ',KLstM1
    WRITE(*,*)' KLst+1 = ',KLstP1
    WRITE(*,*)' Associated K-1 ',ASSOCIATED(LinkArray1(KLstM1)%L)
    WRITE(*,*)' Associated K+1 ',ASSOCIATED(LinkArray1(KLstP1)%L)
    WRITE(*,*)' L-1  = ',LinkArray1(KLstM1)%L%LNum
    WRITE(*,*)' L+1  = ',LinkArray1(KLstP1)%L%LNum
    STOP ' Bad logic in RhoEcon '

  END SUBROUTINE RhoEcon

  SUBROUTINE RhoPop(GM,BS,CS,NoWrap,SpinM,Psv,Pair,HGLink,NLink)
    !
    TYPE(CRDS)                              :: GM
    TYPE(BSET)                              :: BS
    TYPE(CellSet)                           :: CS
    INTEGER                                 :: SpinM,NLink
    TYPE(AtomPair)                          :: Pair
    TYPE(PrimPair)                          :: Prim
    TYPE(HGLL),POINTER                      :: HGLink,HGClone
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB,SpinM) :: Psv !<<<SPIN
    !
    REAL(DOUBLE),DIMENSION(3)               :: Q,P

    INTEGER                                 :: KA,KB,NBFA,NBFB,CFA,CFB,PFA,PFB, &
                                               IndexA,StartLA,StopLA,MaxLA, &
                                               IndexB,StartLB,StopLB,MaxLB, &
                                               IE,OffCo,LMN,LMNA,LMNB, &
                                               IA,IB,EllA,EllB,LAB,MAB,NAB, &
                                               LenKet,AtA,AtB,iSMat,NC,STATUS
    !
    REAL(DOUBLE)                            :: ZetaA,ZetaB,ZetaAB,ZetaIn,XiAB,ExpAB, &
                                               AB2,Ax,Ay,Az,Bx,By,Bz, &
                                               Px,Py,Pz,PAx,PAy,PAz,PBx,PBy,PBz, &
                                               Ex,Exy,Exyz,CA,CB,MaxAmp           , zz,MaxCo
    !
    LOGICAL                                 :: AEQB,NoWrap

    REAL(DOUBLE),DIMENSION(3) :: R,DIPOLE
    REAL(DOUBLE),DIMENSION(HGLen) :: TmpCo

    !
    KA   = Pair%KA
    KB   = Pair%KB
    NBFA = Pair%NA
    NBFB = Pair%NB
    AEQB = Pair%SameAtom
    !
    Prim%A=Pair%A
    Prim%B=Pair%B

    Prim%AB2=Pair%AB2
    Prim%KA=Pair%KA
    Prim%KB=Pair%KB
    !
    DO CFA=1,BS%NCFnc%I(KA)
       IndexA  = CFBlokDex(BS,CFA,KA)
       StartLA = BS%LStrt%I(CFA,KA)
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       DO CFB=1,BS%NCFnc%I(KB)
          IndexB  = CFBlokDex(BS,CFB,KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          Prim%CFA=CFA
          Prim%CFB=CFB
          Prim%Ell=MaxLA+MaxLB
          DO PFA=1,BS%NPFnc%I(CFA,KA)                !
             Prim%PFA=PFA
             DO PFB=1,BS%NPFnc%I(CFB,KB)
                Prim%PFB=PFB
                Prim%ZA=BS%Expnt%D(PFA,CFA,KA)
                Prim%ZB=BS%Expnt%D(PFB,CFB,KB)
                Prim%Zeta=Prim%ZA+Prim%ZB
                Prim%Xi=Prim%ZA*Prim%ZB/Prim%Zeta
                !
                IF(TestPrimPair(Prim%Xi,Prim%AB2)) THEN
                   !
                   Prim%P=(Prim%ZA*Prim%A+Prim%ZB*Prim%B)/Prim%Zeta
                   !
                   CALL PWrap(GM,Prim,.NOT.NoWrap)
                   !
                   CALL SetBraBlocks(Prim,BS,CompPrim_O=.FALSE.)
                   !
                   LenKet = LHGTF(Prim%Ell)
                   !
                   TmpCo(1:LenKet)=Zero
                   IA=IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      EllA = BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                         EllB = BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)
                         DO iSMat=1,1 !SpinM
                            DO LMN=1,LHGTF(EllA+EllB)
                               TmpCo(LMN)=TmpCo(LMN)+HGBra%D(LMN,IA,IB)*Psv(IA,IB,iSMat)
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                   !
                   ! Fullup cutoff for zero coefficients
                   MaxCo=-BIG_DBL
                   DO LMN=1,LenKet
                      MaxCo=MAX(MaxCo,ABS(TmpCo(LMN)))
                   ENDDO
                   ! See similar criteria in MondoMods/Thresholding.f90 in
                   ! Extent/PFunk:  Basically, we are looking for HGTFs that
                   ! are zero relative to the penetration error.  The factor
                   ! of 1D18 is used instead of 1D20, since economization
                   ! can change things somewhat.  This should probably all get
                   ! tightened up at some point.
                   IF((TauPAC/MaxCo)*(Prim%Zeta/Pi)**(1.5D0)<1D18)THEN
                      ALLOCATE(HGLink%Next,STAT=Status)
                      IF(Status/=SUCCEED)CALL Halt(' Node ALLOCATE failed in RhoPop ')
                      HGLink=>HGLink%Next
                      NLink=NLink+1
                      HGLink%LNum=NLink
                      HGLink%Ell=Prim%Ell
                      HGLink%Zeta=Prim%Zeta
                      HGLink%Cent=Prim%Pw
                      NULLIFY(HGLink%Next)
                      CALL New(HGLink%Coef,LenKet)
                      HGLink%Coef%D(1:LenKet)=TmpCo(1:LenKet)
                   ENDIF
                ENDIF
                !
             ENDDO
          ENDDO
       ENDDO
    ENDDO


  END SUBROUTINE RhoPop

  FUNCTION DPole(LQ,PiExpt,RX,RY,RZ,Coef)
    INTEGER                   :: LQ
    REAL(DOUBLE)              :: RX,RY,RZ,Expt,PiExpt
    REAL(DOUBLE),DIMENSION(3) :: DPole
    REAL(DOUBLE),DIMENSION(:) :: Coef
!
    SELECT CASE(LQ)
    CASE (0)
       DPole(1) = -PiExpt*Coef(1)*RX
       DPole(2) = -PiExpt*Coef(1)*RY
       DPole(3) = -PiExpt*Coef(1)*RZ
    CASE(1:)
       DPole(1) = -PiExpt*(Coef(2)+Coef(1)*RX)
       DPole(2) = -PiExpt*(Coef(3)+Coef(1)*RY)
       DPole(3) = -PiExpt*(Coef(4)+Coef(1)*RZ)
    END SELECT
    ! !
  END FUNCTION DPole
  !========================================================================================
  ! From HG Coefs, calculate the Quadrupole
  !========================================================================================
  FUNCTION QPole(LQ,Expt,PiExpt,RX,RY,RZ,Coef)
    INTEGER                   :: LQ
    REAL(DOUBLE)              :: RX,RY,RZ,PiExpt,Expt,RX2,RY2,RZ2
    REAL(DOUBLE),DIMENSION(6) :: QPole
    REAL(DOUBLE),DIMENSION(:) :: Coef
!
    RX2 = Half/Expt + RX*RX
    RY2 = Half/Expt + RY*RY
    RZ2 = Half/Expt + RZ*RZ
!!$

!!$    WRITE(*,*)' half'
!!$    !
!!$    RX2=+RX*RX
!!$    RY2=+RY*RY
!!$    RZ2=+RZ*RZ
    !
    SELECT CASE(LQ)
    CASE(0)
       QPole(1) = PiExpt*Coef(1)*RX2
       QPole(2) = PiExpt*Coef(1)*RY2
       QPole(3) = PiExpt*Coef(1)*RZ2
       QPole(4) = PiExpt*Coef(1)*RX*RY
       QPole(5) = PiExpt*Coef(1)*RX*RZ
       QPole(6) = PiExpt*Coef(1)*RY*RZ
    CASE(1)
       QPole(1) = PiExpt*(Coef(1)*RX2 + Two*Coef(2)*RX)
       QPole(2) = PiExpt*(Coef(1)*RY2 + Two*Coef(3)*RY)
       QPole(3) = PiExpt*(Coef(1)*RZ2 + Two*Coef(4)*RZ)
       QPole(4) = PiExpt*(Coef(1)*RX*RY + Coef(2)*RY + Coef(3)*RX)
       QPole(5) = PiExpt*(Coef(1)*RX*RZ + Coef(2)*RZ + Coef(4)*RX)
       QPole(6) = PiExpt*(Coef(1)*RY*RZ + Coef(3)*RZ + Coef(4)*RY)
    CASE(2:)
       QPole(1) = PiExpt*(Coef(1)*RX2 + Two*Coef(2)*RX + Two*Coef(5))
       QPole(2) = PiExpt*(Coef(1)*RY2 + Two*Coef(3)*RY + Two*Coef(7))
       QPole(3) = PiExpt*(Coef(1)*RZ2 + Two*Coef(4)*RZ + Two*Coef(10))
       QPole(4) = PiExpt*(Coef(1)*RX*RY + Coef(2)*RY + Coef(3)*RX + Coef(6))
       QPole(5) = PiExpt*(Coef(1)*RX*RZ + Coef(2)*RZ + Coef(4)*RX + Coef(8))
       QPole(6) = PiExpt*(Coef(1)*RY*RZ + Coef(3)*RZ + Coef(4)*RY + Coef(9))
    END SELECT

  END FUNCTION QPole

  !
END MODULE Density
!  0.6904445367316159
!  0.6904445367316714
