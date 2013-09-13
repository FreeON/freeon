MODULE QCTCIndexing
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalObjects
  USE ProcessControl
  USE PoleGlobals
  USE MemMan
  USE Indexing
  IMPLICIT NONE
  !---------------------------------------------------------------------
  INTEGER,PARAMETER :: HGEll4=2*HGEll
  INTEGER, DIMENSION(-1:HGEll4,-1:HGEll4,-1:HGEll4) :: LMNx
#ifdef POINTERS_IN_DERIVED_TYPES
  TYPE(INT_VECT),DIMENSION(:,:),POINTER   :: PLMNx,QLMNx,PQLMNx
  INTEGER, DIMENSION(:,:), POINTER        :: XLLen
  TYPE(INT_RNK2),DIMENSION(:,:),POINTER   :: XLIdx
  TYPE(DBL_RNK2),DIMENSION(:,:),POINTER   :: XLSgn
  INTEGER, DIMENSION(:,:,:), POINTER      :: CTLen
#else
  TYPE(INT_VECT),DIMENSION(:,:),ALLOCATABLE   :: PLMNx,QLMNx,PQLMNx
  INTEGER, DIMENSION(:,:), ALLOCATABLE        :: XLLen
  TYPE(INT_RNK2),DIMENSION(:,:),ALLOCATABLE   :: XLIdx
  TYPE(DBL_RNK2),DIMENSION(:,:),ALLOCATABLE   :: XLSgn
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE        :: CTLen
#endif
CONTAINS

  SUBROUTINE TensorIndexingSetUp()
    CALL XLSetUp()
    CALL HGSetUp()
  END SUBROUTINE TensorIndexingSetUp
  !===================================================================================================
  !
  !===================================================================================================

SUBROUTINE CTSetUp
    INTEGER                         :: Count
    INTEGER                         :: L1,L2,L3,M1,M2,M3,LDX1,LDX2,LDX3,ABSM2,ABSM3
    REAL(DOUBLE)                    :: CN,SN,CMN,SMN
    REAL(DOUBLE)                    ::  TMP,SGN,TMPC,TWO,TMPS,SGNL,SGNLM
    INTEGER  :: L,LP,LQ,LL,M,LDX,K,KK,LLKK,LLKM,N,KDX,LKDX,ID,CCCount,SSCount,SCCount,CSCount
    INTEGER  :: Count0,Count1,Count2,Count3,Count4,Count5,MinA,MinB,MinC
    REAL(DOUBLE) :: TwoC
    ID(L)=L*(L+1)/2


    ALLOCATE(CTLen(0:HGEll,0:MaxPoleEll,0:5))

    DO LP=0,HGEll
       LQ=MaxPoleEll

       COUNT0=0
       COUNT1=0
       COUNT2=0
       COUNT3=0
       COUNT4=0
       COUNT5=0

      IF(LP.EQ.0)THEN
         DO l=0,LQ
            ll=ID(l)
            Count0=Count0+1
!            Cp(0)=Cp(0)+Cpq(ll)*Cq(ll)+Spq(ll)*Sq(ll)
            DO m=1,l
              ldx=ll+m
              Count0=Count0+1
!              Tmp=Tmp+Cpq(ldx)*Cq(ldx)+Spq(ldx)*Sq(ldx)
           ENDDO
        ENDDO
!        Cp(0)=Cp(0)+2.0D0*Tmp
      ENDIF

      Sgn=1.0D0
      DO l=0,LP
         ll=ID(l)
         TmpC=0.0D0
         DO k=0,LQ
            kk=ID(k)
            llkk=ID(l+k)

            Count2=Count2+1
         ENDDO
         Sgn=-Sgn
      ENDDO

      Two=-2.0D0
      DO l=1,LP
         ll=ID(l)
         DO k=0,LQ
            kk=ID(k)
            llkk=ID(l+k)
            DO m=1,l
               ldx=ll+m
               llkm=llkk+m

               Count2=Count2+1
               Count4=Count4+1

            ENDDO
         ENDDO
         Two=-Two
      ENDDO
!
      Two=2.0D0
      DO l=0,LP
         ll=ID(l)
         TmpC=0.0D0
         DO k=0,LQ
            kk=ID(k)
            llkk=ID(l+k)
            DO n=1,k
               kdx=kk+n
               lkdx=llkk+n

               Count2=Count2+1
               Count3=Count3+1

            ENDDO
         ENDDO
         Two=-Two
      ENDDO

      SgnL=1.0D0
      TwoC=-2.0D0
      DO l=1,LP
         ll=ID(l)
         DO k=1,LQ
            kk=ID(k)
            llkk=ID(l+k)
            SgnLM=SgnL
            DO m=1,l
               ldx=ll+m
               llkm=llkk+m
               Two=2.0D0*SgnL
               MinA=MIN(m,k)
               MinB=MIN(m-1,k)
               MinC=K
               DO n=1,MinB
                  kdx=kk+n
                  lkdx=llkm-n

                  Count1=Count1+1
                  Two=-Two
               ENDDO

               DO n=1,k
                  kdx=kk+n
                  lkdx=llkm+n

                  Count1=Count1+1
               ENDDO

               DO n=MinB+1,MinA
                  kdx=kk+n
                  lkdx=llkm-n

                  Count2=Count2+1
                  Count5=Count5+1

                  Two=-Two
               ENDDO

               llkm=llkk-m
               DO n=m+1,k
                  kdx=kk+n
                  lkdx=llkm+n

                  Count1=Count1+1

               ENDDO
               SgnLM=-SgnLM
            ENDDO
         ENDDO
         SgnL=-SgnL
         TwoC=-TwoC
      ENDDO

      WRITE(*,22)LP,LQ,Count0,Count1,Count2,Count3,Count4,Count5
22    format(5(I5,", "))

#ifdef LSDJFLSDJFLSJFE

      IF(LP.EQ.0)THEN
         Tmp=0.0D0
         DO l=0,LQ
            ll=ID(l)
            Cp(0)=Cp(0)+Cpq(ll)*Cq(ll)+Spq(ll)*Sq(ll)
            DO m=1,l
              ldx=ll+m
              Tmp=Tmp+Cpq(ldx)*Cq(ldx)+Spq(ldx)*Sq(ldx)
           ENDDO
        ENDDO
        Cp(0)=Cp(0)+2.0D0*Tmp
        RETURN
      ENDIF

      Sgn=1.0D0
      DO l=0,LP
         ll=ID(l)
         TmpC=0.0D0
         DO k=0,LQ
            kk=ID(k)
            llkk=ID(l+k)

            Count2=Count2+1
            CTIdx2(LP,LQ)%I(1,Count2)=ll
            CTIdx2(LP,LQ)%I(2,Count2)=llkk
            CTIdx2(LP,LQ)%I(3,Count2)=kk
            CTFlt2(LP,LQ)%I(Count2)=Sgn

!            TmpC=TmpC+Cpq(llkk)*Cq(kk)
         ENDDO
!         Cp(ll)=Cp(ll)+Sgn*TmpC
         Sgn=-Sgn
      ENDDO

      Two=-2.0D0
      DO l=1,LP
         ll=ID(l)
         DO k=0,LQ
            kk=ID(k)
            llkk=ID(l+k)
            DO m=1,l
               ldx=ll+m
               llkm=llkk+m

               Count2=Count2+1
               CTIdx2(LP,LQ)%I(1,Count2)=ldx
               CTIdx2(LP,LQ)%I(2,Count2)=llkm
               CTIdx2(LP,LQ)%I(3,Count2)=kk
               CTFlt2(LP,LQ)%I(Count2)=Two

!               Cp(ldx)=Cp(ldx)+Two*Cpq(llkm)*Cq(kk)

               Count4=Count4+1
               CTIdx4(LP,LQ)%I(1,Count4)=ldx
               CTIdx4(LP,LQ)%I(2,Count4)=llkm
               CTIdx4(LP,LQ)%I(3,Count4)=kk
               CTFlt4(LP,LQ)%I(Count4)=Two

!               Sp(ldx)=Sp(ldx)+Two*Spq(llkm)*Cq(kk)


            ENDDO
         ENDDO
         Two=-Two
      ENDDO
C
      Two=2.0D0
      DO l=0,LP
         ll=ID(l)
         TmpC=0.0D0
         DO k=0,LQ
            kk=ID(k)
            llkk=ID(l+k)
            DO n=1,k
               kdx=kk+n
               lkdx=llkk+n

               Count2=Count2+1
               CTIdx2(LP,LQ)%I(1,Count2)=ll
               CTIdx2(LP,LQ)%I(2,Count2)=lkdx
               CTIdx2(LP,LQ)%I(3,Count2)=kdx
               CTFlt2(LP,LQ)%I(Count2)=Two
!
!               Cp(ll)=Cp(ll)+Two*Cpq(lkdx)*Cq(kdx)

               Count3=Count3+1
               CTIdx3(LP,LQ)%I(1,Count3)=ll
               CTIdx3(LP,LQ)%I(2,Count3)=lkdx
               CTIdx3(LP,LQ)%I(3,Count3)=kdx
               CTFlt3(LP,LQ)%I(Count3)=Two
!
!               Cp(ll)=Cp(ll)+Two*Spq(lkdx)*Sq(kdx)

            ENDDO
         ENDDO
         Two=-Two
      ENDDO

      SgnL=1.0D0
      TwoC=-2.0D0
      DO l=1,LP
         ll=ID(l)
         DO k=1,LQ
            kk=ID(k)
            llkk=ID(l+k)
            SgnLM=SgnL
            DO m=1,l
               ldx=ll+m
               llkm=llkk+m
               Two=2.0D0*SgnL
               MinA=MIN(m,k)
               MinB=MIN(m-1,k)
               MinC=K
               DO n=1,MinB
                  kdx=kk+n
                  lkdx=llkm-n

                  Count1=Count1+1
                  CTIdx1(LP,LQ)%I(1,Count1)=ldx
                  CTIdx1(LP,LQ)%I(2,Count1)=lkdx
                  CTIdx1(LP,LQ)%I(3,Count1)=kdx
                  CTFlt1(LP,LQ)%I(1,Count1)=Two
                  CTFlt1(LP,LQ)%I(2,Count1)=-Two
                  CTFlt1(LP,LQ)%I(3,Count1)=Two
                  CTFlt1(LP,LQ)%I(4,Count1)=Two
!
!                  Cp(ldx)=Cp(ldx)+Two*Cpq(lkdx)*Cq(kdx)
!                  Cp(ldx)=Cp(ldx)-Two*Spq(lkdx)*Sq(kdx)
!                  Sp(ldx)=Sp(ldx)+Two*Cpq(lkdx)*Sq(kdx)
!                  Sp(ldx)=Sp(ldx)+Two*Spq(lkdx)*Cq(kdx)
!
                  Two=-Two
               ENDDO

               DO n=1,k
                  kdx=kk+n
                  lkdx=llkm+n

                  Count1=Count1+1
                  CTIdx1(LP,LQ)%I(1,Count1)=ldx
                  CTIdx1(LP,LQ)%I(2,Count1)=lkdx
                  CTIdx1(LP,LQ)%I(3,Count1)=kdx
                  CTFlt1(LP,LQ)%I(1,Count1)=TwoC
                  CTFlt1(LP,LQ)%I(2,Count1)=TwoC
                  CTFlt1(LP,LQ)%I(3,Count1)=-TwoC
                  CTFlt1(LP,LQ)%I(4,Count1)=TwoC
!
!                  Cp(ldx)=Cp(ldx)+TwoC*Cpq(lkdx)*Cq(kdx)
!                  Cp(ldx)=Cp(ldx)+TwoC*Spq(lkdx)*Sq(kdx)
!                  Sp(ldx)=Sp(ldx)-TwoC*Cpq(lkdx)*Sq(kdx)
!                  Sp(ldx)=Sp(ldx)+TwoC*Spq(lkdx)*Cq(kdx)
               ENDDO

               DO n=MinB+1,MinA
                  kdx=kk+n
                  lkdx=llkm-n

                  Count2=Count2+1
                  CTIdx2(LP,LQ)%I(1,Count2)=ldx
                  CTIdx2(LP,LQ)%I(2,Count2)=lkdx
                  CTIdx2(LP,LQ)%I(3,Count2)=kdx
                  CTFlt2(LP,LQ)%I(Count2)=Two
!
!                  Cp(ldx)=Cp(ldx)+Two*Cpq(lkdx)*Cq(kdx)

                  Count5=Count5+1
                  CTIdx5(LP,LQ)%I(1,Count5)=ldx
                  CTIdx5(LP,LQ)%I(2,Count5)=lkdx
                  CTIdx5(LP,LQ)%I(3,Count5)=kdx
                  CTFlt5(LP,LQ)%I(Count5)=Two
!
!                  Sp(ldx)=Sp(ldx)+Two*Cpq(lkdx)*Sq(kdx)
                  Two=-Two
               ENDDO

               llkm=llkk-m
               DO n=m+1,k
                  kdx=kk+n
                  lkdx=llkm+n

                  Count1=Count1+1
                  CTIdx1(LP,LQ)%I(1,Count1)=ldx
                  CTIdx1(LP,LQ)%I(2,Count1)=lkdx
                  CTIdx1(LP,LQ)%I(3,Count1)=kdx
                  CTFlt1(LP,LQ)%I(1,Count1)=SgnLM*2.0D0
                  CTFlt1(LP,LQ)%I(2,Count1)=SgnLM*2.0D0
                  CTFlt1(LP,LQ)%I(3,Count1)=SgnLM*2.0D0
                  CTFlt1(LP,LQ)%I(4,Count1)=-SgnLM*2.0D0

!                  Cp(ldx)=Cp(ldx)+SgnLM*2.0D0*Cpq(lkdx)*Cq(kdx)
!                  Cp(ldx)=Cp(ldx)+SgnLM*2.0D0*Spq(lkdx)*Sq(kdx)
!                  Sp(ldx)=Sp(ldx)+SgnLM*2.0D0*Cpq(lkdx)*Sq(kdx)
!                  Sp(ldx)=Sp(ldx)-SgnLM*2.0D0*Spq(lkdx)*Cq(kdx)
               ENDDO
               SgnLM=-SgnLM
            ENDDO
         ENDDO
         SgnL=-SgnL
         TwoC=-TwoC
      ENDDO

#endif

   ENDDO

 END SUBROUTINE CTSetUp

 SUBROUTINE XLSetUp
   INTEGER                         :: LP,LQ,Count
   INTEGER                         :: L1,L2,L3,M1,M2,M3,LDX1,LDX2,LDX3,ABSM2,ABSM3
   REAL(DOUBLE)                    :: CN,SN,CMN,SMN


   ALLOCATE(XLLen(0:MaxPFFFEll,0:MaxPoleEll))
   ALLOCATE(XLSgn(0:MaxPFFFEll,0:MaxPoleEll))
   ALLOCATE(XLIdx(0:MaxPFFFEll,0:MaxPoleEll))

   DO LP=MaxPoleEll,MaxPFFFEll
      DO LQ=0,MaxPoleEll

         IF(LQ.LE.HGEll.OR.LQ.EQ.LP)THEN

            Count=0
            DO L1 = 0,LP
               DO M1 = 0,L1
                  LDX1 = LTD(L1)+M1
                  DO L2 = 0,MIN(L1,LQ)
                     L3    = L1-L2
                     DO M2 = -L2,L2
                        ABSM2 = ABS(M2)
                        LDX2  = LTD(L2)+ABSM2
                        M3    = M1-M2
                        ABSM3 = ABS(M3)
                        LDX3  = LTD(L3)+ABSM3
                        IF(ABSM3 .LE. L3)Count=Count+1
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            !
!            WRITE(*,*)LP,LQ,Count
            XLLen(LP,LQ)=Count
            CALL New(XLSgn(LP,LQ),(/4,Count/))
            CALL New(XLIdx(LP,LQ),(/3,Count/))
            !
            Count=0
            DO L1 = 0,LP
               DO M1 = 0,L1
                  LDX1 = LTD(L1)+M1
                  DO L2 = 0,MIN(L1,LQ)
                     L3    = L1-L2
                     DO M2 = -L2,L2
                        ABSM2 = ABS(M2)
                        LDX2  = LTD(L2)+ABSM2
                        M3    = M1-M2
                        ABSM3 = ABS(M3)
                        LDX3  = LTD(L3)+ABSM3
                        IF(ABSM3 .LE. L3) THEN
                           Count=Count+1
                           IF(M2 .LT. 0) THEN
                              CN = (-One)**(ABSM2)
                              SN = -CN
                           ELSE
                              CN = One
                              SN = One
                           ENDIF
                           IF(M3 .LT. 0) THEN
                              CMN = (-One)**(ABSM3)
                              SMN = -CMN
                           ELSE
                              CMN = One
                              SMN = One
                           ENDIF
                           XLSgn(LP,LQ)%D(1,Count)=CN*CMN
                           XLSgn(LP,LQ)%D(2,Count)=-SN*SMN
                           XLSgn(LP,LQ)%D(3,Count)=CN*SMN
                           XLSgn(LP,LQ)%D(4,Count)=SN*CMN
                           XLIdx(LP,LQ)%I(1,Count)=LDX1
                           XLIdx(LP,LQ)%I(2,Count)=LDX2
                           XLIdx(LP,LQ)%I(3,Count)=LDX3
                           !                       Cp(LDX1) = Cp(LDX1)+CN*CMN*Cq(LDX2)*Cpq(LDX3) &
                           !                                          -SN*SMN*Sq(LDX2)*Spq(LDX3)
                           !                       Sp(LDX1) = Sp(LDX1)+CN*SMN*Cq(LDX2)*Spq(LDX3) &
                           !                                          +SN*CMN*Sq(LDX2)*Cpq(LDX3)

                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO
   ENDDO
 END SUBROUTINE XLSetUp
  !===================================================================================================
  !
  !===================================================================================================
  SUBROUTINE HGSetUp
    USE Order
    INTEGER L,M,N,EllP,EllQ,LP,MP,NP,LQ,MQ,NQ,CountP,CountQ,CountPQ,PLen,QLen,PQLen,Ix,Iy,Iz,J
    TYPE(INT_VECT) :: Key,Ord,Imp
    INTEGER(INT8),EXTERNAL:: Interleave,HilbertKey
    character(len=72) :: filename

    ! Precompute HGTF addressing
    LMNx(-1,:,:)=1
    LMNx(:,-1,:)=1
    LMNx(:,:,-1)=1
    DO L=0,HGEll4
       DO M=0,HGEll4-L
          DO N=0,HGEll4-L-M
             LMNx(L,M,N)=LMNDex(L,M,N)
          ENDDO
       ENDDO
    ENDDO

    ALLOCATE(PLMNx(0:HGEll+1,0:HGEll))
    ALLOCATE(QLMNx(0:HGEll+1,0:HGEll))
    ALLOCATE(PQLMNx(0:HGEll+1,0:HGEll))

    ! Precompute all the
    DO EllP=0,HGEll
       DO EllQ=0,HGEll
          PLen=LHGTF(EllP)
          QLen=LHGTF(EllQ)
          PQLen=PLen*QLen
          !
          CALL New(PLMNx(EllP,EllQ),PQLen)
          CALL New(QLMNx(EllP,EllQ),PQLen)
          CALL New(PQLMNx(EllP,EllQ),PQLen)
          !
          CountPQ=0
          DO LP=0,EllP
             DO MP=0,EllP-LP
                DO NP=0,EllP-LP-MP
                   DO LQ=0,EllQ
                      DO MQ=0,EllQ-LQ
                         DO NQ=0,EllQ-LQ-MQ
                            CountPQ=CountPQ+1
                            QLMNx(EllP,EllQ)%I(CountPQ)=LMNx(LQ,MQ,NQ)
                            PLMNx(EllP,EllQ)%I(CountPQ)=LMNx(LP,MP,NP)
                            PQLMNx(EllP,EllQ)%I(CountPQ)=LMNx(LP+LQ,MP+MQ,NP+NQ)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

          N=PQLen

          CALL New(Key,PQLen)
          CALL New(Ord,PQLen)
          CALL New(Imp,PQLen)


          DO J=1,N
             IF(EllP>EllQ)THEN
                Ix=0
                Iy=PLMNx(EllP,EllQ)%I(J)
             ELSE
                Ix=QLMNx(EllP,EllQ)%I(J)
                Iy=0
             ENDIF
             Iz=PQLMNx(EllP,EllQ)%I(J)
             !             Key%I(J)=Iz
             Key%I(J)=Interleave(21,Ix,Iy,Iz)
             !             Key%I(J)=HilbertKey(21,Key%I(J))
             Ord%I(J)=J
          END DO

          CALL QSORTI(Ord%I(1),N,Key%I(1))
          !
          DO J=1,N
             IMP%I(J)=QLMNx(EllP,EllQ)%I(Ord%I(J))
          ENDDO
          DO J=1,N
             QLMNx(EllP,EllQ)%I(J)=IMP%I(J)
          ENDDO
          !
          DO J=1,N
             IMP%I(J)=PLMNx(EllP,EllQ)%I(Ord%I(J))
          ENDDO
          DO J=1,N
             PLMNx(EllP,EllQ)%I(J)=IMP%I(J)
          ENDDO
          !
          DO J=1,N
             IMP%I(J)=PQLMNx(EllP,EllQ)%I(Ord%I(J))
          ENDDO
          DO J=1,N
             PQLMNx(EllP,EllQ)%I(J)=IMP%I(J)
          ENDDO
!!$
!!$
!!$
!!$          ketcontract 327
!!$          md3trr4     180
!!$          total      1086

!!$         filename='Spc'//TRIM(IntToChar(EllP))//'c'//TRIM(IntToChar(EllQ))//'.dat'
!!$
!!$         OPEN(UNIT=77,FILE=FileName, &
!!$              ACCESS='SEQUENTIAL',FORM='FORMATTED', &
!!$              STATUS='NEW')
!!$
!!$          DO J=1,N
!!$             WRITE(77,22)J,KEY%I(Ord%I(J)),PLMNx(EllP,EllQ)%I(J), &
!!$                        QLMNx(EllP,EllQ)%I(J),PQLMNx(EllP,EllQ)%I(J)
!!$
!!$             WRITE(*,22)J,KEY%I(Ord%I(J)),PLMNx(EllP,EllQ)%I(J), &
!!$                        QLMNx(EllP,EllQ)%I(J),PQLMNx(EllP,EllQ)%I(J)
!!$             22 FORMAT(I5," ",I10,4(I4,' '))
!!$          ENDDO
!!$          close(77)

          CALL Delete(Key)
          CALL Delete(Ord)
          CALL Delete(IMP)

       ENDDO
    ENDDO

    !
  END SUBROUTINE HGSetUp


END MODULE QCTCIndexing
