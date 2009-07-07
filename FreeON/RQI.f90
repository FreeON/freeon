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
MODULE RayleighQuotientIteration
  USE InOut
  USE Macros
  USE MemMan
  USE Overlay
  USE Indexing
  USE AtomPairs
  USE PunchHDF
  USE PrettyPrint
  USE GlobalScalars
  USE ControlStructures
  USE OptionKeys
  USE McMurchie
  USE MondoLogger

  IMPLICIT NONE

  TYPE(TIME)          :: TimeTotal,TimeONX,TimeQCTC,TimeBCSR
CONTAINS


  SUBROUTINE TDSCF(C)
    IMPLICIT NONE
    TYPE(Controls)                             :: C
    IF(.NOT. C%POpt%Resp%TD_SCF) RETURN
    CALL SetFrontEndMacros(C%Geos,C%Sets)
    CALL RQI(NBasF,4,C%Nams,C%Opts,C%Stat,C%MPIs,C%Geos, &
             C%Sets%BSets(1,C%Sets%NBSets))
  END SUBROUTINE TDSCF
  !===============================================================================
  !   Subroutine to load the global macro parameters used mostly in the backend
  !   into the front end, so that we can do kluge work. Assumes that we are
  !   at the end of our basis set list, first clone etc. NOT APPROPRIATE for
  !   general work.
  !===============================================================================
  SUBROUTINE SetFrontEndMacros(G,B)
    TYPE(BasisSets)  :: B
    TYPE(Geometries) :: G
    INTEGER          :: II
    !
    MyClone=1
    NBasF=B%BSets(1,B%NBSets)%NBasF
    NAToms=G%Clone(1)%NAtms
    MaxAtms=B%MxAts(B%NBSets)
    MaxBlks=B%MxBlk(B%NBSets)
    MaxNon0=B%MxN0s(B%NBSets)
    CALL New(BSiz,NAtoms)
    CALL New(OffS,NAtoms)
    BSiz%I=B%BSiz(1,B%NBSets)%I
    OffS%I=B%OffS(1,B%NBSets)%I
    MaxBlkSize=0
    DO II=1,G%Clone(1)%NAtms
       MaxBlkSize=MAX(MaxBlkSize,BSiz%I(II))
    ENDDO
  END SUBROUTINE SetFrontEndMacros

  SUBROUTINE RQI(N,M,Nam,O,S,MPI,G,B)
    INTEGER            :: N,M,I,J,K,L,U,V,JTDA,cBAS
    TYPE(FileNames)    :: Nam,RQINams
    TYPE(State)        :: S,RQIStat
    TYPE(Parallel)     :: MPI
    TYPE(Options)      :: O
    TYPE(Geometries)   :: G
    TYPE(BSET)         :: B
    LOGICAL            :: DoTDA
    !
    TYPE(BCSR)                      :: sP,sQ,sF,sZ
    TYPE(DBL_RNK2)                  :: P,Q,F,Z
    REAL(DOUBLE)                    :: Ek,EkOld,dEk,Beta,Lambda,ErrRel,ErrAbs,Shift,dNorm
    INTEGER                         :: XkNon0s,PkNon0s
    REAL(DOUBLE)                    :: XkThreshold,PkThreshold,PMax


!!$    REAL(DOUBLE),DIMENSION(N,N)     :: Xk,Gk,Pk,LXk,LPk
!!$    REAL(DOUBLE),DIMENSION(N,N)     :: XkOld,PkOld,GkOld
!!$    REAL(DOUBLE),DIMENSION(N)       :: Values
!!$    REAL(DOUBLE),DIMENSION(N,N,M)   :: Vectors
!!$    REAL(DOUBLE),DIMENSION(N,N,N,N) :: TwoE,DoubleSlash
    !
    RQIStat=S
    RQINams=Nam
    cBAS=RQIStat%Current%I(2)
    !
    CALL Elapsed_TIME(TimeTotal,Init_O='Init')
    CALL Elapsed_TIME(TimeONX,Init_O='Init')
    CALL Elapsed_TIME(TimeQCTC,Init_O='Init')
    CALL Elapsed_TIME(TimeBCSR,Init_O='Init')
    CALL Elapsed_TIME(TimeTotal,Init_O='Start')

!!$    !
!!$    CALL Get(sP,TrixFile("OrthoD",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0))
!!$    CALL Get(sF,TrixFile("OrthoF",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0))
!!$    CALL Get(sZ,TrixFile("X",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I))
!!$    !
!!$    CALL SetEq(P,sP)
!!$    CALL SetEq(Q,sP)
!!$    P%D=Two*P%D
!!$    Q%D=-Two*Q%D
!!$    DO I=1,N
!!$       Q%D(I,I)=Q%D(I,I)+Two
!!$    ENDDO
!!$    !
!!$    CALL SetEq(F,sF)
!!$    CALL SetEq(Z,sZ)
!!$    !
!!$    CALL Integrals2E(B,G,TwoE)
!!$    DO I=1,N
!!$       DO J=1,N
!!$          DO K=1,N
!!$             DO L=1,N
!!$                DoubleSlash(I,J,K,L)=TwoE(I,J,K,L)-TwoE(I,K,J,L)/2D0
!!$             ENDDO
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
!!$
!!$    goto 111

    DO I=1,1
       !       write(*,*)' ek = ',EK
       !       GOTO 111 ! STOP
       DO JTDA=0,0!1
          !
          IF(JTDA==0)THEN
!             CALL ResetThresholds(cBAS,Nam,O,G,1D2)
             DoTDA=.TRUE.
             XkThreshold=O%Thresholds(cBAS)%Trix
             PkThreshold=O%Thresholds(cBAS)%Trix
             CALL Nihilate0(N,I,Nam,S,MPI)
             CALL LOn2BakEnd(N,I,0,Shift,'Xk',Nam,S,MPI,XkThreshold,XkNon0s)
             CALL Nihilate1(I,Nam,S,MPI,Ek,TDA_O=.TRUE.)
          ELSE
!             CALL ResetThresholds(cBAS,Nam,O,G,1D-2)
             DoTDA=.FALSE.
             XkThreshold=O%Thresholds(cBAS)%Trix
             PkThreshold=O%Thresholds(cBAS)%Trix
             CALL LOn2BakEnd(N,I,0,Shift,'Xk',Nam,S,MPI,XkThreshold,XkNon0s)
             CALL Nihilate1(I,Nam,S,MPI,Ek,TDA_O=.FALSE.)
          ENDIF
          !
          EkOld=BIG_DBL
          DO K=0,200
             ! The non-linear congjugate gradient
             CALL NLCGBakEnd(I,K,Ek,Nam,S,MPI,Beta)
             ! Compute L[Pk]
             CALL LOn2BakEnd(N,I,K,Shift,'Pk',Nam,S,MPI,PkThreshold,PkNon0s,PMax)
             ! Line Search: Min_Lambda{ E[Xk+Lambda*Pk] }
             CALL RQLSBakEnd(I,Nam,S,Lambda)
             ! Anhiliate and renorm Xk
             CALL NihilateXk(I,Nam,S,MPI,Lambda,dNorm,TDA_O=DoTDA)
             ! Compute L[Xk]
             CALL LOn2BakEnd(N,I,K,Shift,'Xk',Nam,S,MPI,XkThreshold,XkNon0s)
             ! Anihilate L[Xk], compute Ek and its relative error
             CALL NihilateLXk(I,Nam,S,MPI,Ek,dEk,TDA_O=DoTDA)
             !
             CALL OpenASCII(Nam%OFile,Out)
             WRITE(*  ,33)I,K,Ek*27.21139613182D0,dEk,ABS(dNorm),TimeONX%Wall,TimeQCTC%Wall,TimeBCSR%Wall, &
                        PMax,100D0*DBLE(XkNon0s)/DBLE(N*N),100D0*DBLE(PkNon0s)/DBLE(N*N)
             WRITE(Out,33)I,K,Ek*27.21139613182D0,dEk,ABS(dNorm),TimeONX%Wall,TimeQCTC%Wall,TimeBCSR%Wall, &
                        PMax,100D0*DBLE(XkNon0s)/DBLE(N*N),100D0*DBLE(PkNon0s)/DBLE(N*N)
             CLOSE(Out,STATUS='KEEP')
             !
             IF(K>3.AND.dNorm<1D-2)THEN
                ! Look for bad behavior
                IF( Ek > EkOld .AND. ABS((Ek-EkOld)/Ek) > O%Thresholds(cBAS)%ETol )THEN

                   ! Sign of variational principle broken, ostensibly due to N-scaling
                   ! approximaitons.  If this happens, we are DONE!
                   WRITE(*,*)' Converged due to variational violation ',Ek,EkOld,  &
                              (Ek-EkOld)/Ek, O%Thresholds(cBAS)%ETol*1D2
                   Ek=EkOld
                   EXIT
                ENDIF
                ! Look for convergence (may be too tight)
                IF(ABS(dEk)  <O%Thresholds(cBAS)%ETol*1D2.AND. &
                   ABS(dNorm)<O%Thresholds(cBAS)%DTol )THEN
                   WRITE(*,*)' Met convergence criteria ',dEk,dNorm,O%Thresholds(cBAS)%ETol, &
                              O%Thresholds(cBAS)%DTol
                   EXIT
                ENDIF
             ENDIF
             EkOld=MIN(EkOld,Ek)
             !
          ENDDO
          CALL Elapsed_TIME(TimeTOTAL,Init_O='Accum')
          IF(DoTDA)THEN
             WRITE(*,44)'TDA:',I,K,Ek*27.21139613182D0,dEk,TimeTotal%Wall
          ELSE
             WRITE(*,44)'RPA:',I,K,Ek*27.21139613182D0,dEk,TimeTotal%Wall
          ENDIF
       ENDDO
    ENDDO
    !
33  FORMAT('St=',I2,', It=',I3,', Ev=',F10.6,', dE=',D8.2,', dN=',D7.2, &
           ', Tk=',D10.4,', Tj=',D10.4,', Tm=',D10.4,', |Gk|=',D8.2,', %Xk=',F6.2,', %Pk=',F6.2)


44  FORMAT(A4,' State=',I2,', Nk=',I3,', Ev=',F9.6,', dE=',D7.2,', WallSec=',D12.4)

!!$
!!$
!!$111 CONTINUE
!!$    DO I=1,1
!!$       !
!!$       DO JTDA=0,1
!!$          IF(JTDA==0)THEN
!!$             DoTDA=.TRUE.
!!$             CALL RPAGuess(N,Xk)
!!$          ELSE
!!$             DoTDA=.FALSE.
!!$             Xk=Vectors(:,:,I)
!!$          ENDIF
!!$
!!$          CALL Anihilate(N,P%D,Q%D,Xk,TDA_O=DoTDA)
!!$          CALL Renorm(N,P%D,Xk)
!!$          CALL LOn2(N,I,Shift,F%D,P%D,Z%D,DoubleSlash,Values,Vectors,Xk,LXk)
!!$          CALL Anihilate(N,P%D,Q%D,LXk,TDA_O=DoTDA)
!!$
!!$          Beta=Zero
!!$          XkOld=Zero
!!$          PkOld=Zero
!!$          Ek=ThoulessQ(N,P%D,Xk,LXk)
!!$          DO K=0,200
!!$             !
!!$             Gk=Two*(LXk-Ek*Xk)
!!$             IF(K>0)Beta=Pdot1(N,P%D,Gk,Gk-Gkold)/Pdot1(N,P%D,GkOld,GkOld)
!!$
!!$             Pk=Gk+Beta*PkOld
!!$             !
!!$             CALL LOn2(N,I,Shift,F%D,P%D,Z%D,DoubleSlash,Values,Vectors,Pk,LPk)
!!$             CALL RQILineSearch(N,P%D,Pk,Xk,LXk,LPk,Lambda)
!!$             !
!!$             EkOld=Ek
!!$             XkOld=Xk
!!$             EkOld=Ek
!!$             GkOld=Gk
!!$             PkOld=Pk
!!$             !
!!$             Xk=XkOld+Lambda*Pk
!!$             !
!!$             CALL Anihilate(N,P%D,Q%D,Xk,TDA_O=DoTDA)
!!$             dNorm=One-sqrt(abs(Pdot1(N,P%D,Xk,Xk)))
!!$             CALL ReNorm(N,P%D,Xk)
!!$             CALL LOn2(N,I,Shift,F%D,P%D,Z%D,DoubleSlash,Values,Vectors,Xk,LXk)
!!$             CALL Anihilate(N,P%D,Q%D,LXk,TDA_O=DoTDA)
!!$             Ek=ThoulessQ(N,P%D,Xk,LXk)
!!$             ErrAbs=Ek-EkOld
!!$             ErrRel=-1D10
!!$             DO U=1,N
!!$                DO V=1,N
!!$                   ErrRel=MAX(ErrRel,ABS(Ek*Xk(U,V)-LXk(U,V))/Ek)
!!$                ENDDO
!!$             ENDDO
!!$             WRITE(*,33)I,K,Ek*27.21139613182D0,Beta,Lambda,ErrRel,dNorm
!!$
!!$             IF(ErrRel<1D-4)EXIT
!!$             !
!!$          ENDDO
!!$          WRITE(*,*)I,K,Ek*27.21139613182D0,ErrRel,ErrAbs
!!$          Values(I)=Ek
!!$          Vectors(:,:,I)=Xk
!!$
!!$       ENDDO
!!$    ENDDO
!!$
  END SUBROUTINE RQI


  SUBROUTINE Nihilate0(N,I,Nam,S,MPI)
    !
    INTEGER               :: I,N,  ii,jj
    REAL(DOUBLE)          :: Norm
    TYPE(FileNames)       :: Nam
    TYPE(State)           :: S
    TYPE(Parallel)        :: MPI
    TYPE(BCSR)            :: sP,sQ,sXk,sI,sT,sT1 ! Nihilate0 delete list
    INTEGER, DIMENSION(3) :: Cur
    CHARACTER(LEN=DCL)    :: XkName,PName,QName
    !
    TYPE(DBL_RNK2)        :: X

    !-----------------------------------------------------------------------------
    CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
    ! Status same as ground state, but now SCF cycle# is RQI state#
    Cur=S%Current%I
    Cur(1)=I
    ! Naming of things
    XkName=TrixFile('OrthoXk',PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=Cur,OffSet_O=0)
    PName= TrixFile("OrthoD", PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QName= TrixFile("OrthoQ", PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    ! Get occupied projector
    CALL Get(sP,PName)
    ! Generate virtual space projector if first time through
    IF(I==1)THEN
       CALL SetEq(sQ,sP)
       CALL Multiply(sQ,-One)
       CALL Add(sQ,One)
       CALL Put(sQ,QName)
    ELSE
       CALL Get(sQ,QName)
    ENDIF
    ! Start guess with particle-hole space ONLY. Starting in this space
    ! GUARANTEES that we will only ever have positive line search solutions.
    ! This is equivalent to the TDA.
    !
!!$    CALL New(X,(/N,N/))
!!$    X%D=Zero
!!$    do ii=1,N
!!$       do jj=MAX(1,ii-5),MIN(N,II+5)
!!$          X%D(ii,jj)=RANDOM_DBL((/-One,One/))
!!$       enddo
!!$    enddo
!!$    CALL SetEq(sXk,X)
!!$    CALL Delete(X)
!!
    CALL SetEq(sXk,sP)
    DO II=1,sXk%NNon0
       sXk%MTrix%D(II)=sXk%MTrix%D(II)*Zero + sXk%MTrix%D(II)*RANDOM_DBL((/-One,One/))
    ENDDO

    CALL Multiply(sQ,sXk,sT1)
    CALL Multiply(sT1,sP,sXk)
    CALL Delete(sQ)
    CALL Delete(sT1)

    ! Normalize the guess transition density.  Note factor
    ! of two has to do with normalization of P and Q.  Here,
    ! these projectors have eigenvalues == 1, not 2.
    CALL XPose(sXk,sT)
    Norm=SQRT(Two*ABS(OneDot(sP,sXk,sT)))
    sXk%MTrix%D=sXk%MTrix%D*(One/Norm)

!    CALL PPrint(sXk,"NORMED",Unit_O=6)

    ! Put guess Xk to disk
    CALL Put(sXk,XkName)

    ! Tidy up
    CALL Delete(sT)
    CALL Delete(sP)
    CALL Delete(sXk)
     ! All done
    CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
  END SUBROUTINE Nihilate0


  SUBROUTINE Nihilate1(I,N,S,MPI,Ek,TDA_O)
    LOGICAL, OPTIONAL     :: TDA_O
    LOGICAL               :: TDA
    INTEGER               :: I
    REAL(DOUBLE)          :: Ek
    TYPE(FileNames)       :: N
    TYPE(State)           :: S
    TYPE(Parallel)        :: MPI
    TYPE(BCSR)            :: sP,sQ,sXk,sLXk,sT,sT1,sT2,sT3 ! Nihilate1 delete list
    INTEGER, DIMENSION(3) :: Cur
    CHARACTER(LEN=DCL)    :: XkName,LXkName,PName,QName
    !-----------------------------------------------------------------------------
    CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
    IF(PRESENT(TDA_O))THEN
       TDA=TDA_O
    ELSE
       TDA=.FALSE.
    ENDIF
    ! Status same as ground state, but now SCF cycle# is RQI state#
    Cur=S%Current%I
    Cur(1)=I
    ! Naming of things
    XkName=TrixFile('OrthoXk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
    LXkName=TrixFile('LXk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
    PName= TrixFile("OrthoD", PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QName= TrixFile("OrthoQ", PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    ! Get occupied projector
    CALL Get(sP,PName)
    CALL Get(sQ,QName)
    ! Kluge due to poor memory management:
    CALL New(ST1)
    CALL New(SLXK)
    ! Get that should have good management (but doesnt)
    CALL Get(sLXk,LXkName)
    ! Anihilate via TDA or TD-SCF symmetry.
    ! Note that factor of (1/4) is not present, due to the
    ! normalization of P to eigenvalues with 1s or 0s.
    IF(TDA)THEN
       ! LXk=Q.LXk.P (only p-h zapped)
       CALL Multiply(sQ,sLXk,sT1)
       CALL Multiply(sT1,sP,sLXk)
       CALL Delete(sT1)
    ELSE
       ! LXk=(P.LXk.Q+Q.LXk.P) (both h-p and p-h zapped)
       CALL Multiply(sP,sLXk,sT1)
       CALL Multiply(sT1,sQ,sT2)
       CALL Multiply(sQ,sLXk,sT1)
       CALL Multiply(sT1,sP,sT3)
       CALL Add(sT2,sT3,sLXk)
       CALL Delete(sT1)
       CALL Delete(sT2)
       CALL Delete(sT3)
    ENDIF
    !
    CALL Delete(sQ)
    !
    CALL Put(sLXk,LXkName)
    !
    CALL Get(sXk,XkName)
    CALL XPose(sLXk,sT)
    !
    ! Again a factor of two for normalization of P
    Ek=Two*OneDot(sP,sXk,sT)
    !
    CALL Delete(sP)
    CALL Delete(sT)
    CALL Delete(sXk)
    CALL Delete(sLXk)
    ! All done
    CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
  END SUBROUTINE Nihilate1

  SUBROUTINE NihilateXk(I,N,S,MPI,Lambda,dNorm,TDA_O)
    LOGICAL, OPTIONAL     :: TDA_O
    LOGICAL               :: TDA
    INTEGER               :: I
    REAL(DOUBLE)          :: Lambda,Norm,dNorm
    TYPE(FileNames)       :: N
    TYPE(State)           :: S
    TYPE(Parallel)        :: MPI
    TYPE(BCSR)            :: sP,sQ,sXk,sPk,sLXk,sT,sT1,sT2,sT3 ! NihilateXk delete list
    INTEGER, DIMENSION(3) :: Cur
    CHARACTER(LEN=DCL)    :: XkName,PkName,LXkName,PName,QName
    !-----------------------------------------------------------------------------
    CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
    IF(PRESENT(TDA_O))THEN
       TDA=TDA_O
    ELSE
       TDA=.FALSE.
    ENDIF
    ! Status same as ground state, but now SCF cycle# is RQI state#
    Cur=S%Current%I
    Cur(1)=I

    ! Naming of things

    XkName=TrixFile('OrthoXk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
    PkName=TrixFile('OrthoPk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
    PName= TrixFile("OrthoD", PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QName= TrixFile("OrthoQ", PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)

    CALL New(sT1)
    CALL New(sPk)
    CALL New(sXk)

    CALL Get(sT1,XkName)
    CALL Get(sPk,PkName)
    sPk%MTrix%D=sPk%MTrix%D*Lambda
    CALL Add(sT1,sPk,sXk)
    !
    CALL Delete(sPk)

    CALL Get(sP,PName)
    CALL Get(sQ,QName)

    ! Anihilate via TDA or TD-SCF symmetry.
    ! Note that factor of (1/4) is not present, due to the
    ! normalization of P to eigenvalues with 1s or 0s.
    IF(TDA)THEN
       ! Xk=Q.Xk.P (only p-h zapped)
       CALL Multiply(sQ,sXk,sT1)
       CALL Multiply(sT1,sP,sXk)
       CALL Delete(sT1)
    ELSE
       ! Xk=(P.Xk.Q+Q.Xk.P) (both h-p and p-h zapped)
       CALL Multiply(sP,sXk,sT1)
       CALL Multiply(sT1,sQ,sT2)
       CALL Multiply(sQ,sXk,sT1)
       CALL Multiply(sT1,sP,sT3)
       CALL Add(sT2,sT3,sXk)
       CALL Delete(sT1)
       CALL Delete(sT2)
       CALL Delete(sT3)
    ENDIF
    !
    CALL Delete(sQ)
    ! Normalize the guess transition density.  Note factor
    ! of two has to do with normalization of P and Q.  Here,
    ! these projectors have eigenvalues == 1, not 2.
    CALL XPose(sXk,sT)
    Norm=SQRT(Two*ABS(OneDot(sP,sXk,sT)))
    dNorm=One-Norm
    sXk%MTrix%D=sXk%MTrix%D*(One/Norm)
    ! Put guess Xk to disk
    CALL Put(sXk,XkName)
    !
    CALL Delete(sT)
    CALL Delete(sP)
    CALL Delete(sXk)
    ! All done
    CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
  END SUBROUTINE NihilateXk
  !
  SUBROUTINE NihilateLXk(I,N,S,MPI,Ek,dEk,TDA_O)
    LOGICAL, OPTIONAL     :: TDA_O
    LOGICAL               :: TDA
    INTEGER               :: I,L
    REAL(DOUBLE)          :: dEk,Ek
    TYPE(FileNames)       :: N
    TYPE(State)           :: S
    TYPE(Parallel)        :: MPI
    TYPE(BCSR)            :: sP,sQ,sXk,sPk,sLXk,sT,sT1,sT2,sT3 ! NihilateLXk delete list
    INTEGER, DIMENSION(3) :: Cur
    CHARACTER(LEN=DCL)    :: XkName,PkName,LXkName,PName,QName
    !-----------------------------------------------------------------------------
    CALL Elapsed_TIME(TimeBCSR,Init_O='Start')

    IF(PRESENT(TDA_O))THEN
       TDA=TDA_O
    ELSE
       TDA=.FALSE.
    ENDIF
    ! Status same as ground state, but now SCF cycle# is RQI state#
    Cur=S%Current%I
    Cur(1)=I

    ! Naming of things

    LXkName=TrixFile('LXk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
    XkName=TrixFile('OrthoXk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
    PkName=TrixFile('OrthoPk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
    PName= TrixFile("OrthoD", PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    QName= TrixFile("OrthoQ", PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    !

    CALL New(sLXk)
    !
    CALL Get(sP,PName)
    CALL Get(sQ,QName)
    CALL Get(sLXk,LXkName)
    ! Anihilate via TDA or TD-SCF symmetry.
    ! Note that factor of (1/4) is not present, due to the
    ! normalization of P to eigenvalues with 1s or 0s.
    IF(TDA)THEN
       ! Xk=Q.Xk.P (only p-h zapped)
       CALL Multiply(sQ,sLXk,sT1)
       CALL Multiply(sT1,sP,sLXk)
       CALL Delete(sT1)
    ELSE
       ! Xk=(P.Xk.Q+Q.Xk.P) (both h-p and p-h zapped)
       CALL Multiply(sP,sLXk,sT1)
       CALL Multiply(sT1,sQ,sT2)
       CALL Multiply(sQ,sLXk,sT1)
       CALL Multiply(sT1,sP,sT3)
       CALL Add(sT2,sT3,sLXk)
       CALL Delete(sT1)
       CALL Delete(sT2)
       CALL Delete(sT3)
    ENDIF
    !
    CALL Delete(sQ)
    ! Put the anihilated LXk to disk
    CALL Put(sLXk,LXkName)
    !
    CALL Get(sXk,XkName)
    CALL XPose(sLXk,sT)
    ! Again a factor of two for normalization of P
    Ek=Two*OneDot(sP,sXk,sT)
    CALL Delete(sP)
    CALL Delete(sT)
    !
    sXk%MTrix%D=-Ek*sXk%MTrix%D
    CALL Add(sXk,sLXk,sT1)
    !
    dEk=-1D10
    DO L=1,sT1%NNon0
       dEk=MAX(dEk,ABS(sT1%MTrix%D(L)))
    ENDDO
    dEk=dEk/Ek
    !
    CALL Delete(sT1)
    CALL Delete(sXk)
    CALL Delete(sLXk)
    ! All done

    CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')

  END SUBROUTINE NihilateLXk

  SUBROUTINE NLCGBakEnd(I,K,Ek,N,S,MPI,Beta)
    !
    INTEGER               :: I,K
    REAL(DOUBLE)          :: Ek
    TYPE(FileNames)       :: N
    TYPE(State)           :: S
    TYPE(Parallel)        :: MPI
    TYPE(BCSR),SAVE       :: sXk,sLXk,sGk,sGkOld,sPkOld,sP,sT,sT1 ! NLCGradient delete list
    INTEGER, DIMENSION(3) :: Cur
    REAL(DOUBLE)          :: Num,Den,Beta
    CHARACTER(LEN=DCL)    :: XkName,LXkName,GkName,PName,GkOldName,PkOldName,PkName

    CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
    !
    Cur=S%Current%I
    Cur(1)=I
    !
    XkName=   TrixFile('OrthoXk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
    LXkName=  TrixFile('LXk',    PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
    GkName=   TrixFile('Gk',     PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
    GkOldName=TrixFile('GkOld',  PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)
    PkName=   TrixFile('OrthoPk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0)

    PkOldName=PkName

    PName=    TrixFile("OrthoD", PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0)
    !
    ! Get Xk and LXk
    CALL Get(sXk,XkName)
    CALL Get(sLXk,LXkName)

    ! Gk=Two*(LXk-Ek*Xk)

    sXk%MTrix%D=-Ek*sXk%MTrix%D

    CALL Add(sLXk,sXk,sGk)

    sGk%MTrix%D=Two*sGk%MTrix%D

!    CALL PPrint(sGk,   ' Gk    ',Unit_O=6)

    ! Gk to disk
    CALL Put(sGk,GkName)
    ! Done with Xk and LXk
    CALL Delete(sXk)
    CALL Delete(sLXk)
    ! Beta=0 for K=0
    IF(K==0)THEN
       CALL Put(sGk,PkName)
       CALL Put(sGk,GkOldName)
       CALL Delete(sGk)
       RETURN
    ENDIF
    ! Get Pk and GkOld
    CALL Get(sP,PName)
    CALL Get(sGkOld,GkOldName)

!    CALL PPrint(sGkOld,' GkOld ',Unit_O=6)

!
    ! Beta=(Gk-Gkold,Gk)_p/(GkOld,GkOld)_p
    sGkOld%MTrix%D=-sGkOld%MTrix%D
    CALL Add(sGk,sGkOld,sT1)
    sGkOld%MTrix%D=-sGkOld%MTrix%D
    CALL XPose(sGk,sT)
    Num=OneDot(sP,sT1,sT)
    CALL XPose(sGkOld,sT)
    Den=OneDot(sP,sGkOld,sT)
    Beta=Num/Den
!    WRITE(*,*)' Beta = ',Beta


    ! Done with P and T
    CALL Delete(sT)
    CALL Delete(sP)
    CALL Delete(sGkOld)
    ! Pk=Gk+Beta*PkOld
    CALL Get(sPkOld,PkOldName)

    sPkOld%MTrix%D=Beta*sPkOld%MTrix%D

    CALL Add(sGk,sPkOld,sT1)
    ! Put GkOld to disk
    CALL Put(sGk,GkOldName)
    ! Put Pk to disk
    CALL Put(sT1,PkName)
    ! Clean up and done
    CALL Delete(sT1)
    CALL Delete(sGk)
    CALL Delete(sPkOld)

    CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
    !
  END SUBROUTINE NLCGBAKEND

  SUBROUTINE LOn2BakEnd(N,I,K,Shift,Trgt,Nam,S,MPI,MatrixThreshold,MatrixNon0s,PMax_O)
    INTEGER                       :: N,I,J,K,MatrixNon0s
    TYPE(FileNames)               :: Nam
    TYPE(State)                   :: S,RQIStat
    TYPE(Parallel)                :: MPI
    REAL(DOUBLE),DIMENSION(N,N)   :: LX,X
    TYPE(BCSR)                    :: sF,sX,sJ,sK,sZ,sP,sJK,sT1,sT2,sT3
    REAL(DOUBLE)                  :: Shift,MatrixThreshold,LocalThreshold,PMax
    CHARACTER(LEN=*)              :: Trgt
    REAL(DOUBLE),OPTIONAL         :: PMax_O
    !----------------------------------------------------------------------------
    CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
    ! Set up local Invokation parameters based on ground state values
    CALL New(RQIStat%Action,2)
    CALL New(RQIStat%Current,3)
    CALL New(RQIStat%Previous,3)
    RQIStat%Current%I=S%Current%I
    RQIStat%Previous%I=S%Previous%I
    ! Action is TD-SCF with secondary parameter the product LX or LP (L[Xk], L[Pk])
    RQIStat%Action%C(1)="TD-SCF"
    RQIStat%Action%C(2)=TRIM(Trgt) !//TRIM(IntToChar(I))
    ! "SCF cycle" is the RQI State number
    RQIStat%Current%I(1)=I
    RQIStat%Previous%I(1)=I
    ! Get the orthogonal transition density matrix Xk (or CG gradient Pk) corresponding to the
    ! resultant of this subroutine, namely L[Xk] or L[Pk] in an orthongal representation
    CALL New(sX) ! Kluge.  Somehow, dimensioning not quite right here:
    CALL Get(sX,TrixFile('Ortho'//TRIM(Trgt),PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=RQIStat%Current%I,OffSet_O=0))
    ! Ground state fockian in an orthogonal representation
    CALL Get(sF,TrixFile("OrthoF",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0))
    ! T2=[Xor,For]
    CALL Multiply(sX,sF,sT2)
    CALL Multiply(sF,sX,sT2,-One)
    ! Done with F
    CALL Delete(sF)
    ! Z is the sparse inverse factor of S
    CALL Get(sZ,TrixFile("X",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I))
    ! Xao = Z^t.Xor.Z
    CALL Multiply(sZ,sX,sT1)
    CALL Multiply(sT1,sZ,sX)

    IF(PRESENT(PMax_O))THEN
!       PMax_O=MAX(sX)
       PMax_O=FNorm(sX)
       LocalThreshold=1D1*MatrixThreshold*PMax_O
       LocalThreshold=MatrixThreshold
    ELSE
       LocalThreshold=MatrixThreshold
    ENDIF
    ! Filter small blocks and return the # of non zero elements
    CALL Filter(sX,Tol_O=LocalThreshold)
    MatrixNon0s=sX%NNon0
    ! This is the AO transition density matrix (or CG gradient)
    CALL Put(sX,TrixFile(Trgt,PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=RQIStat%Current%I,OffSet_O=0))
    ! Done with sX
    CALL Delete(sX)
    CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
    ! Build JK[X] in the AO basis
    CALL Elapsed_TIME(TimeONX,Init_O='Start')
    CALL Invoke('ONX',Nam,RQIStat,MPI)
    CALL Elapsed_TIME(TimeONX,Init_O='Accum')
    CALL Elapsed_TIME(TimeQCTC,Init_O='Start')
    CALL Invoke('QCTC',Nam,RQIStat,MPI)
    CALL Elapsed_TIME(TimeQCTC,Init_O='Accum')
    CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
    ! Pick up J and K
    CALL Get(sJ,TrixFile("J",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=RQIStat%Current%I,OffSet_O=0))
    CALL Get(sK,TrixFile("K",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=RQIStat%Current%I,OffSet_O=0))
    ! JK=Jao[X]+Kao[X]
    CALL Add(sJ,sK,sJK)
    ! Done with J and K
    CALL Delete(sJ)
    CALL Delete(sK)
    ! JK[X]=Zt.JKao[X].Z==JKor
    CALL Multiply(sZ,sJK,sT1)
    CALL Multiply(sT1,sZ,sJK)
    ! Done with Z
    CALL Delete(sZ)
    ! Get some P
    CALL Get(sP,TrixFile("OrthoD",PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0))
    ! T1=[JKor,Por]
    CALL Multiply(sP,sJK,sT1)
    CALL Multiply(sJK,sP,sT1,-One)
    ! Done with P and JK
    CALL Delete(sP)
    CALL Delete(sJK)
    ! L[Xk]=[F,Xk]+[P,JK[X]] (orthogonal)
    CALL Add(sT1,sT2,sT3)
    ! Done with temporaries 1 and 2
    CALL Delete(sT1)
    CALL Delete(sT2)
    ! Put orthogonal L[Xk] or L[Pk] to disk (*.LX or *.LP)
    CALL Put(sT3,TrixFile("L"//TRIM(Trgt),PWD_O=Nam%M_SCRATCH,Name_O=Nam%SCF_NAME,Stats_O=RQIStat%Current%I,OffSet_O=0))
    ! Done with temp #3
    CALL Delete(sT3)
    ! Done with invokation parameters
    CALL Delete(RQIStat%Action)
    CALL Delete(RQIStat%Current)
    CALL Delete(RQIStat%Previous)
    CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
    !----------------------------------------------------------------------------
    IF(I>1)STOP
  END SUBROUTINE LOn2BakEnd

  SUBROUTINE RQLSBakEnd(I,N,S,Lambda)
    INTEGER               :: I
    TYPE(FileNames)       :: N
    TYPE(State)           :: S
    TYPE(Parallel)        :: MPI
    TYPE(BCSR)            :: sP,sXk,sPk,sLXk,sLPk,sT
    INTEGER, DIMENSION(3) :: Cur
    REAL(DOUBLE)          :: Lambda,Lambda_p,Lambda_m
    REAL(DOUBLE)          :: XX,PP,XP,PX,PLP,XLP,XLX,PLX,AA,BB,CC
    !
    CALL Elapsed_TIME(TimeBCSR,Init_O='Start')
    !

    Cur=S%Current%I
    Cur(1)=I
    !
    CALL Get(sP ,TrixFile("OrthoD",PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=S%Current%I,OffSet_O=0))
    CALL Get(sXk,TrixFile('OrthoXk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0))
    CALL Get(sPk,TrixFile('OrthoPk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0))
    !
    CALL XPose(sXk,sT)
    XX =OneDot(sP,sXk,sT)
    PX =OneDot(sP,sPk,sT)
    !
    CALL XPose(sPk,sT)
    PP =OneDot(sP,sPk,sT)
    XP =OneDot(sP,sXk,sT)
    !
    CALL Get(sLXk,TrixFile('LXk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0))
    CALL XPose(sLXk,sT)
    CALL Delete(sLXk)
    XLX=OneDot(sP,sXk,sT)
    PLX=OneDot(sP,sPk,sT)
    !
    CALL Get(sLPk,TrixFile('LPk',PWD_O=N%M_SCRATCH,Name_O=N%SCF_NAME,Stats_O=Cur,OffSet_O=0))
    CALL XPose(sLPk,sT)
    CALL Delete(sLPk)
    PLP=OneDot(sP,sPk,sT)
    XLP=OneDot(sP,sXk,sT)
    !
    CALL Delete(sP)
    CALL Delete(sT)
    CALL Delete(sXk)
    CALL Delete(sPk)
    !
    AA=PLP*(PX+XP)-PP*(PLX+XLP)
    BB=2.0*PLP*XX-2.0*XLX*PP
    CC=XX*(PLX+XLP)-XLX*(PX+XP)
    !
    lambda_p=(-BB+SQRT(BB*BB-4.0*AA*CC))/(2.0*AA)
    lambda_m=(-BB-SQRT(BB*BB-4.0*AA*CC))/(2.0*AA)
    !
    Lambda=Lambda_P
    !
    CALL Elapsed_TIME(TimeBCSR,Init_O='Accum')
    !
  END SUBROUTINE RQLSBakEnd
  !===============================================================================
  ! HERE ARE THE CONVENTIONAL DENSE MATRIX RQI ROUTINES
  !===============================================================================
  !---------------------------------------------------------------------
  ! Calculates TD-SCF scalar product Tr=Tr([A,P],B^+)
  ! Assumes B^t on input
  !---------------------------------------------------------------------
  FUNCTION OneDot(sP,sA,sBt) RESULT(Tr)
    TYPE(BCSR) :: sP,sA,sBt,sT1
    REAL(DOUBLE) :: Tr
    CALL Multiply(sP,sA,sT1)
    CALL Multiply(sA,sP,sT1,-One)
    Tr=Half*Trace(sT1,sBt)
    CALL Delete(sT1)
  END FUNCTION OneDot
  !
  SUBROUTINE Anihilate(N,P,Q,X,TDA_O)
    LOGICAL, OPTIONAL :: TDA_O
    LOGICAL :: TDA
    INTEGER :: N
    REAL(DOUBLE),DIMENSION(N,N) :: P,Q,X
    IF(PRESENT(TDA_O))THEN
       TDA=TDA_O
    ELSE
       TDA=.FALSE.
    ENDIF
    IF(TDA)THEN
       X=ProjectPH(N,P,Q,X)
    ELSE
       X=Project(N,P,Q,X)
    ENDIF
  END SUBROUTINE Anihilate

  SUBROUTINE LOn2(N,M,Shift,F,P,Z,TwoE,Values,Vectors,X,LX)
    INTEGER :: N,M,J
    REAL(DOUBLE),DIMENSION(N,N)     :: F,P,Z,X,LX,AA,BB,Temp,Com
    REAL(DOUBLE),DIMENSION(N,N,N,N) :: TwoE
    REAL(DOUBLE)                  :: Shift,OmegaPls,OmegaMns,WS
    REAL(DOUBLE),DIMENSION(:)     :: Values
    REAL(DOUBLE),DIMENSION(:,:,:) :: Vectors
    LX=LiouvAO(N,F  ,P  ,Z,TwoE,X )
    IF(M==1)RETURN
    WS=Values(M-1)-Values(1)+Shift
    Com=MATMUL(X,P)-MATMUL(P,X)
    DO J=1,M-1
       OmegaPls=Trace2(MATMUL(TRANSPOSE(Vectors(:,:,J)),Com),N)
       OmegaMns=Trace2(MATMUL(Vectors(:,:,J),Com),N)
       LX=LX+WS*(OmegaMns*TRANSPOSE(Vectors(:,:,J))+OmegaPls*Vectors(:,:,J))
    ENDDO
  END SUBROUTINE LOn2
  !
  Subroutine RQILineSearch(N,P,Pk,Xk,LXk,LPk,Lambda)
    INTEGER :: N
    REAL(DOUBLE) :: Lambda,Lambda_p,Lambda_m
    REAL (DOUBLE),DIMENSION(N,N)::P,Pk,Xk,LXk,LPk,Tmp1
    REAL(DOUBLE) :: XX,PP,XP,PX,PLP,XLP,XLX,PLX,AA,BB,CC
    XX =Pdot1(N,P,Xk,Xk) ! 1 - normalized by definition
    PP =Pdot1(N,P,Pk,Pk)
    XP =Pdot1(N,P,Xk,Pk)
    PX =Pdot1(N,P,Pk,Xk)
    PLP=Pdot1(N,P,Pk,LPk)
    XLX=Pdot1(N,P,Xk,LXk)
    XLP=Pdot1(N,P,Xk,LPk)
    PLX=Pdot1(N,P,Pk,LXk)
    AA=PLP*(PX+XP)-PP*(PLX+XLP)
    BB=2.0*PLP*XX-2.0*XLX*PP
    CC=XX*(PLX+XLP)-XLX*(PX+XP)
    lambda_p=(-BB+SQRT(BB*BB-4.0*AA*CC))/(2.0*AA)
    lambda_m=(-BB-SQRT(BB*BB-4.0*AA*CC))/(2.0*AA)
    Lambda=Lambda_P
  END Subroutine RQILineSearch

  FUNCTION ThoulessQ(N,P,X,LX) RESULT(Ek)
    !! So, where is  the denominator??
    INTEGER :: N
    REAL(DOUBLE) :: Ek
    REAL(DOUBLE),DIMENSION(N,N)     :: F,P,X,LX,Tmp1
    Ek=Pdot1(N,P,X,LX)
  END FUNCTION ThoulessQ

  FUNCTION LiouvDot(N,BB,DSao,temp2)  RESULT(temp1)
    ! Calculates action of the Coulomb operator in AO space temp1=BB * (ij||kl)
    IMPLICIT NONE
    INTEGER :: I,J,K,N,one
    REAL (DOUBLE),DIMENSION(N*N):: BB,temp2
    REAL (DOUBLE),DIMENSION(N,N):: temp1
    REAL(DOUBLE),DIMENSION(N*N,N*N)::	DSao
    REAL(DOUBLE) :: ddot

    one=1
    K=0
    DO I=1,N
       DO J=1,N
          K=K+1
          temp2=DSao(:,K)
          !		temp1(J,I)= ddot(N*N,BB,one,temp2,one)     ! This line is
          temp1(J,I)=DOT_PRODUCT(BB,Temp2)           ! the most CPU consuming step
       ENDDO
    END DO

  END FUNCTION LiouvDot

  FUNCTION LiouvAO(N,For,Por,X,DSao,AA)  RESULT(BB)
    ! Calculates action of the Liouville operator in AO space BB=L AA, (ij||kl)
    IMPLICIT NONE
    INTEGER :: I,J,M,K,L,N,one
    REAL (DOUBLE),DIMENSION(N,N)::For,Por,AA,BB,temp1,temp2,X
    REAL(DOUBLE),DIMENSION(N,N,N,N)::	DSao
    REAL(DOUBLE) :: E,ddot

    ! AA to AO
    one=1
    BB=MATMUL(TRANSPOSE(X),(MATMUL(AA,X)))
!    CALL PPrint(BB,'INPUT 2',Unit_O=6)
    DO I=1,N
       DO J=1,N
          temp1(I,J)= 0.0
          DO K=1,N
             DO L=1,N
                temp1(I,J)=temp1(I,J)+BB(K,L)*DSao(K,L,I,J)
             END DO
          END DO
       END DO
    END DO
 !   CALL PPrint(temp1,'AO_JK[X]',Unit_O=6)
    BB=MATMUL(For,AA)-MATMUL(AA,For)
    ! temp back to orthog
    temp2=MATMUL(TRANSPOSE(X),(MATMUL(temp1,X)))
    BB=BB+MATMUL(temp2,Por)-MATMUL(Por,temp2)
  END FUNCTION LiouvAO

  SUBROUTINE RPAGuess(N,X)
    INTEGER :: N,I,J
    REAL(DOUBLE), DIMENSION(N,N) :: X
    X=One
!!$
    do i=1,N
       do j=1,N
          X(i,j)= One/Two**(I+J)
       enddo
    enddo
!!$       temp1 = ProjectPH(N,Qor,Por,Xk)
!!$       Xk = temp1/sqrt(abs(Pdot1(N,Por,temp1,temp1,tmp1)))
  END SUBROUTINE RPAGuess

  SUBROUTINE ReNorm(N,P,X)
    INTEGER :: N
    REAL(DOUBLE) :: Norm
    REAL(DOUBLE),DIMENSION(N,N) :: P,X
    Norm=sqrt(abs(Pdot1(N,P,X,X)))
    X=X/Norm
  END SUBROUTINE ReNorm

  !************************************************************************
  FUNCTION Pdot(N,P,AA,BB,CC) RESULT(Tr)
    ! Calculates RPA scalar product Tr=Tr([AA^+,P],BB)

    IMPLICIT NONE
    INTEGER :: N
    REAL(DOUBLE) :: Tr
    REAL (DOUBLE),DIMENSION(N,N)::P,AA,BB,CC

    CC=MATMUL((MATMUL(TRANSPOSE(AA),P)-MATMUL(P,TRANSPOSE(AA))),BB)
    Tr=0.5*Trace2(CC,N)

  END FUNCTION Pdot

  !************************************************************************
  FUNCTION Pdot1(N,P,AA,BB) RESULT(Tr)
    ! Calculates RPA scalar product Tr=Tr([AA,P],BB^+)

    IMPLICIT NONE
    INTEGER :: N
    REAL(DOUBLE) :: Tr
    REAL (DOUBLE),DIMENSION(N,N)::P,AA,BB,CC

    CC=MATMUL((MATMUL(AA,P)-MATMUL(P,AA)),TRANSPOSE(BB))
    Tr=0.5*Trace2(CC,N)

  END FUNCTION Pdot1
  !-------------------------------------------------------------------------------
  FUNCTION Project(N,P,Q,AA)  RESULT(BB)
    ! BB=P AA Q + Q AA P
    ! calculates  projection to p-h an h-p space using Q and P

    IMPLICIT NONE
    INTEGER :: N
    REAL (DOUBLE),DIMENSION(N,N)::P,Q,AA,BB

    BB=0.25*(MATMUL(MATMUL(P,AA),Q)+MATMUL(MATMUL(Q,AA),P))

  END FUNCTION Project

  !-------------------------------------------------------------------------------
  FUNCTION ProjectPH(N,P,Q,AA)  RESULT(BB)
    ! BB=Q AA P (X-component, large)   (0  Y)
    ! BB=P AA Q (Y-component, small)   (X  0)
    ! calculates  projection to p-h OR h-p space using Q and P

    IMPLICIT NONE
    INTEGER :: N
    REAL (DOUBLE),DIMENSION(N,N)::P,Q,AA,BB
    BB=0.25D0*(MATMUL(MATMUL(Q,AA),P))
  END FUNCTION ProjectPH

  FUNCTION Trace2(Matrix,N) RESULT(Tr)
    IMPLICIT NONE
    INTEGER :: I,N
    REAL(DOUBLE) :: Tr
    REAL(DOUBLE),DIMENSION(N,N)::Matrix
    Tr=0D0
    DO I=1,N
       Tr=Tr+Matrix(I,I)
    END DO
  END FUNCTION Trace2

END MODULE RayleighQuotientIteration

