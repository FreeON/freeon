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
!    PERFORM O(Lg N) WALKS ON THE POLETREE REPRESENTATION OF THE TOTAL 
!    ELECTRON DENSITY, USING THE K-D TREE DATA STRUCTURE TO APPLY THE 
!    PENETRATION ACCESABILITY CRITERIA
!
!    Author: Matt Challacombe
!------------------------------------------------------------------------------
MODULE TreeWalk
  USE DerivedTypes
  USE PoleNodeType
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  USE Globals
  USE QCTCThresholds
  USE MondoPoles 
  USE GammaFunctions
  USE PoleGlobals
  USE ERIGlobals
  USE Thresholding
  USE InvExp
  USE QCTCIndexing
  USE PoleTree
  USE Clock
  IMPLICIT NONE
#ifdef PAC_DEBUG
  REAL(DOUBLE) :: PC22,TPAC2
  REAL(DOUBLE),DIMENSION(1:HGLen) :: ErrBra
#endif
  REAL(DOUBLE),DIMENSION(1:HGLen) :: HGKet2, HGKet3, HGKet4
#ifdef MAC_DEBUG  
  REAL(DOUBLE) :: RR22,PC22,TPAC2,DeltaPrime
  REAL(DOUBLE),DIMENSION(3) :: POLECENTER
  LOGICAL PrintQ
  REAL(DOUBLE),DIMENSION(0:SPLen)         :: SPErrorKetC,SPErrorKetS
  REAL(DOUBLE),DIMENSION(0:SPLen)         :: SPERRORBRAC,SPERRORBRAS
#endif
!  REAL(DOUBLE),EXTERNAL :: MondoTimer
  REAL(DOUBLE) :: JWalk_Time_Start,Integral_Time_Start,Multipole_Time_Start
  REAL(DOUBLE) :: JWalk_Time,Integral_Time,Multipole_Time
  !---------------------------------------------------------------------
!  REAL(DOUBLE),EXTERNAL :: MONDOTIMER
  INTEGER :: NNearAv,NFarAv,NPrim,NInts
  TYPE PolePointer
     TYPE(PoleNode),POINTER :: P
  END TYPE PolePointer


  INTEGER,PARAMETER                     :: NumNodes=5000
  TYPE(PolePointer),DIMENSION(NumNodes) :: Near,Far
!  REAL(DOUBLE),DIMENSION(2*(1+MaxPoleEll*(MaxPoleEll+3)/2)*NumNodes) :: Wrk
  REAL(DOUBLE),DIMENSION(2000) :: O,W,V
CONTAINS  

  SUBROUTINE JWalk2(QC,P,Nucular_O)
!    USE InvExp
    TYPE(QCPrim)                     :: QC
    TYPE(PoleNode),POINTER           :: P
    LOGICAL,OPTIONAL                 :: Nucular_O
    TYPE(PoleNode),POINTER           :: Q
    LOGICAL                          :: Nucular    
    REAL(DOUBLE)                     :: PQ2
    REAL(DOUBLE)                     :: CoTan,OneOvPQ,OneOvPQxy,RS,SQ,PQToThMnsL,LocalThresh
    REAL(DOUBLE)                     :: TT,TwoC,COne,SOne,CTwo,STwo,RTE,RPE,T,PQHalf,JTau
    REAL(DOUBLE)                     :: SqrtW,RDist,LeftS,RightS,Err,TPAC,RMnsD,R,BigR,SmallR,Xpt,ET
    REAL(DOUBLE)                     :: PQX,PQY,PQZ,OmegaMin,MACError,MACError2,PACError
    INTEGER                          :: is,I,J,K,Ell,EllP,EllQ,LenP,LenQ,LCode
    INTEGER                          :: TotEll,PQEll,PQLen
    LOGICAL                          :: PAC,MAC,Leave
    REAL(DOUBLE),DIMENSION(3)        :: Ext,PQv
    INTEGER                          :: LP,MP,NP,LQ,MQ,NQ,PDex,QDex,MinBraEll,MinKetEll,LPQ
    INTEGER                          :: N,NNear,NFar

!!$    REAL(DOUBLE),DIMENSION(0:10,0:10):: MACTnsr
!!$    REAL(DOUBLE),DIMENSION(NumNodes,3)    :: PQ
!!$    REAL(DOUBLE),DIMENSION(0:2*HGEll,0:2*HGEll,0:2*HGEll,0:2*HGEll) :: MDR
    !-----------------------------------------------------------------------------------------------
    ! 

    IF(PRESENT(Nucular_O))THEN
       Nucular=Nucular_O
    ELSE
       Nucular=.FALSE.
    ENDIF
    !
    JWalk_Time_Start=MTimer()
    !
    EllP=QC%Prim%Ell
    LenP=LHGTF(EllP)
    !
    Q=>P
    NFar=0
    NNear=0
    !
    JTau=TauTwo
    !
    RTE=QC%Prim%Zeta*MinZab
    RPE=QC%Prim%Zeta+MinZab
    OmegaMin=RTE/RPE 
    !    
    DO 
       !
       MAC=.FALSE.
       PQx=QC%Prim%Pw(1)-Q%Box%Center(1)
       PQy=QC%Prim%Pw(2)-Q%Box%Center(2)
       PQz=QC%Prim%Pw(3)-Q%Box%Center(3)
       !
!!$       PQ2=Zero
!!$       IF(PQx<-Q%Box%Half(1))THEN
!!$          PQHalf=PQx+Q%Box%Half(1)
!!$          PQ2=PQ2+PQHalf*PQHalf          
!!$       ELSEIF(PQx>Q%Box%Half(1))THEN
!!$          PQHalf=PQx-Q%Box%Half(1)
!!$          PQ2=PQ2+PQHalf*PQHalf
!!$       ENDIF
!!$       IF(PQy<-Q%Box%Half(2))THEN
!!$          PQHalf=PQy+Q%Box%Half(2)
!!$          PQ2=PQ2+PQHalf*PQHalf
!!$       ELSEIF(PQy>Q%Box%Half(2))THEN
!!$          PQHalf=PQy-Q%Box%Half(2)
!!$          PQ2=PQ2+PQHalf*PQHalf
!!$       ENDIF
!!$       IF(PQz<-Q%Box%Half(3))THEN
!!$          PQHalf=PQz+Q%Box%Half(3)
!!$          PQ2=PQ2+PQHalf*PQHalf
!!$       ELSEIF(PQz>Q%Box%Half(3))THEN
!!$          PQHalf=PQz-Q%Box%Half(3)
!!$          PQ2=PQ2+PQHalf*PQHalf
!!$       ENDIF
!!$       T=OmegaMin*PQ2
!!$       IF(T>8D0)THEN
!!$
       IF(ABS(PQx)>QC%Box%Half(1)+Q%Box%Half(1).OR.  &
          ABS(PQy)>QC%Box%Half(2)+Q%Box%Half(2).OR.  &
          ABS(PQz)>QC%Box%Half(3)+Q%Box%Half(3))THEN 


#ifdef PAC_DEBUG
          CALL ErrCompare(QC,Q,1D-5)
#endif
          ! To here, we are outside the PAC distance and pennetration effects are 
          ! negligible.  Now check to see if the multipole translation error is acceptable; MAC=.TRUE.
          PQx=QC%Prim%Pw(1)-Q%Pole%Center(1)
          PQy=QC%Prim%Pw(2)-Q%Pole%Center(2)
          PQz=QC%Prim%Pw(3)-Q%Pole%Center(3)
          R=SQRT(PQx*PQx+PQy*PQy+PQz*PQz)
          ! Compute the putative multipole translation error
          LCode=100*QC%Prim%Ell+Q%Herm%Ell
          SELECT CASE(LCode)
             INCLUDE "MACErrBnd6.Inc"
          CASE DEFAULT
             CALL Halt('No explicit code for case LCode = '  &
                  //TRIM(IntToChar(LCode))//' in MACErrBnd6.Inc')
          END SELECT
          ! Is it acceptable?
          MAC=MACError<TauMAC      
#ifdef MAC_DEBUG
          CALL MACCompare(QC,Q,MACError)
#endif
       ENDIF
       !        
       IF(MAC)THEN
          NFar=NFar+1
          Far(NFar)%P=>Q
          IF(ASSOCIATED(Q%Next))THEN
             Q=>Q%Next
          ELSE
             EXIT             
          ENDIF
       ELSEIF(Q%Leaf)THEN
          NNear=NNear+1
          Near(NNear)%P=>Q
          IF(ASSOCIATED(Q%Next))THEN
             Q=>Q%Next
          ELSE
             EXIT
          ENDIF
       ELSEIF(ASSOCIATED(Q%Left))THEN
          Q=>Q%Left
       ELSEIF(ASSOCIATED(Q%Right))THEN 
          Q=>Q%Right
       ELSE
          EXIT
       ENDIF
    ENDDO
    !
    !---------------------------------------------------------------------------------------------
    Multipole_Time_Start=MTimer()
    JWalk_Time=JWalk_Time+(Multipole_Time_Start-JWalk_Time_Start)
     LCode=100*QC%Prim%Ell+MaxPoleEll
    SELECT CASE(LCode)	
    INCLUDE "CTraX.Inc"
     !---------------------------------------------------------------------------------------------
!!$    DO N=1,NFar
!!$       Q=>Far(N)%P
!!$       PQ=QC%Prim%Pw-Q%Pole%Center
!!$       CALL CTraXE(QC%Prim%Ell,Q%Pole%Ell,PQ(1),PQ(2),PQ(3), &
!!$            Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
!!$    ENDDO
!!
!!$    DO N=1,NFar
!!$       Q=>Far(N)%P
!!$       LP=QC%Prim%Ell
!!$       LQ=Q%Pole%Ell
!!$       LPQ=LP+LQ 
!!$       PQv=(QC%Prim%Pw-Q%Pole%Center)
!!$       CALL IrRegular(LPQ,PQv(1),PQv(2),PQv(3))
!!$       CALL CTraX77(LP,LQ,SPKetC(0),SPKetS(0),Cpq,Spq,Q%Pole%C(0),Q%Pole%S(0))
!!$    ENDDO
    !---------------------------------------------------------------------------------------------
    Integral_Time_Start=MTimer()
    Multipole_Time=Multipole_Time+(Integral_Time_Start-Multipole_Time_Start)
!!$
!!$    DO N=1,NNear
!!$      is=Near(N)%P%HERM%Stack
!!$      CALL RAhmadiJAlmlof95c(EllP,LenP,NuclearExpnt,JTau,QC%Prim%Zeta,QC%Prim%Pw,HGKet,HGStack(is)) 
!!$   ENDDO
!!$    DO N=1,NNear
!!$       Q=>Near(N)%P
!!$       EllP=QC%Prim%Ell
!!$       LenP=LHGTF(EllP)
!!$       DO EllQ=0,Q%HERM%Ell
!!$          Nq=Q%HERM%Nq(EllQ)
!!$          IF(Nq.NE.0)THEN 
!!$             LenQ=LHGTF(EllQ)
!!$             PQEll=EllP+EllQ
!!$             PQLen=LHGTF(PQEll)
!!$             CALL VMD3(QC%Prim%Ell,LHGTF(QC%Prim%Ell),EllQ,LHGTF(EllQ),          &
!!$                       QC%Prim%Ell+EllQ,LHGTF(QC%Prim%Ell+EllQ),Q%HERM%Nq(EllQ), &
!!$                       QC%Prim%Pw,QC%Prim%Zeta,Q%HERM%Cent(EllQ)%D,              &
!!$                       Q%HERM%Zeta(EllQ)%D,Q%HERM%Coef(EllQ)%D,HGKet)
!!$             NInts=NInts+Nq
!!$          ENDIF
!!$       ENDDO
!!$    ENDDO
!!$
!!$
    DO N=1,NNear
       Q=>Near(N)%P
       EllP=QC%Prim%Ell
       LenP=LHGTF(EllP)
       DO EllQ=0,Q%HERM%Ell
          Nq=Q%HERM%Nq(EllQ)
          IF(Nq.NE.0)THEN 
             LenQ=LHGTF(EllQ)
             PQEll=EllP+EllQ
             PQLen=LHGTF(PQEll)
             CALL RAhmadiJAlmlof95(HGEll4,Nq,EllP,EllQ,LenP,LenQ,PQEll,PQLen,         &
                  NuclearExpnt,JTau,Q%HERM%IHlf(EllQ)%D(1),QC%Prim%Zeta,              &
                  Q%HERM%Zeta(EllQ)%D(1),QC%Prim%Pw(1),Q%HERM%Cent(EllQ)%D(1,1),      &
                  Q%HERM%Coef(EllQ)%D(1,1),HGKet(1))
             NInts=NInts+Nq
          ENDIF
       ENDDO
    ENDDO
    !
    Integral_Time=Integral_Time+(MTimer()-Integral_Time_Start)
    NFarAv=NFarAv+NFar
    NNearAv=NNearAv+NNear

    RETURN
  END SUBROUTINE JWalk2

#ifdef MAC_DEBUG
  !
  RECURSIVE SUBROUTINE MACCompare(P,Q,Err)
    TYPE(QCPrim)                     :: P
    TYPE(PoleNode)                   :: Q
    INTEGER                          :: I,TOTELL
    REAL(DOUBLE)                     :: RDIST,SQRTW,LEFTS,RIGHTS,NEWMACERROR,TRUMACERROR,Err
    REAL(DOUBLE)                     :: PQ2,PQxy,PQx,PQy,PQx2,PQy2,PQz,Omega, &
         CoTan,OneOvPQ,OneOvPQxy,RS,SQ,PQToThMnsL, &
         TT,TwoC,COne,SOne,CTwo,STwo,RTE,RPE,T,Upq,PC2,Err1,Err2,Zeta,Eta
    INTEGER                          :: J,Ell,LCode
    LOGICAL                          :: True
    INTEGER                          :: LP,MP,NP,LQ,MQ,NQ,PDex,QDex,PZ,QZ,LPQ
    INTEGER         :: LM,LenSP
    REAL(DOUBLE),DIMENSION(3)        :: PQ

    TotEll = P%Prim%Ell+Q%Herm%Ell
    NewMACError=Err

    PQ=(P%Prim%Pw-Q%Pole%Center)

!    CALL CTraXE(P%Prim%Ell,Q%Pole%Ell,PQ(1),PQ(2),PQ(3), &
!         Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))

!    DeltaPrime=0D0
!    POLECENTER=Q%Pole%Center
!!$
    SPKetC=0D0
    SPKetS=0D0
    SPErrorKetC=0D0 
    SPErrorKetS=0D0 
!!$
    LP=P%Prim%Ell
    LQ=Q%Pole%Ell
    LPQ=LP+LQ        
    CALL IrRegular(LPQ,PQ(1),PQ(2),PQ(3))
    CALL CTraX77(LP,LQ,SPKetC(0),SPKetS(0),Cpq,Spq,Q%Pole%C(0),Q%Pole%S(0))

    CALL MACError(P,Q)
    ! Here is the differences between point charges and multipole expansions
    LenSP  = LSP(P%Prim%Ell)
    Err1=Zero

!    WRITE(*,*)' SPKETC = ',SPKetC(0:3)
!    WRITE(*,*)' SPErrC = ',SPErrorKetC(0:3)
    DO LM=0,LenSP
       Err1=Err1+ABS(SPERRORBRAC(LM)*(SPKetC(LM)-SPErrorKetC(LM))) &
                +ABS(SPERRORBRAS(LM)*(SPKetS(LM)-SPErrorKetS(LM))) 
    ENDDO
    TruMACError=Err1
    IF(ABS(TruMACError)<1D-20)RETURN
   
    Err1=Zero
    DO LM=0,LenSP
       Err1=Err1+SPERRORBRAC(LM)*SPErrorKetC(LM)+SPERRORBRAS(LM)*SPErrorKetS(LM)
    ENDDO

!    IF(ABS(ERR1).GT.1D30)THEN

    IF(ABS(TruMACError)>1D1*NewMACError.AND.NewMACError>1D-12)THEN

       
       WRITE(*,*)' Numb = ',Q%Box%Number
       WRITE(*,*)' Pw   = ',P%Prim%Pw(1)
       WRITE(*,*)' Zeta = ',P%Prim%Zeta
       WRITE(*,*)' - - - - -  - - - - -  - - - - -  - - - - - '
       WRITE(*,*)' Splt = ',Q%ChargeSplitOK
       WRITE(*,*)' Leaf = ',Q%Leaf
       WRITE(*,*)' TruE = ',TruMACError
       WRITE(*,*)' NewE = ',NewMACError
       WRITE(*,*)' Tier = ',Q%Box%Tier
       IF(Q%LEAF)WRITE(*,*)' Nq   = ',Q%HERM%Nq
       WRITE(*,*)' TotL = ',TotEll
       WRITE(*,*)' PEll = ',P%Prim%Ell       
       WRITE(*,*)' PoleL= ',Q%Pole%Ell
       WRITE(*,*)' Qdelta = ',Q%MAC%Delta
       WRITE(*,*)' DeltaP = ',DeltaPrime

       WRITE(*,*)' QMAC   = ',Q%MAC%O
       WRITE(*,*)' PMAC   = ',P%MAC%O

       WRITE(*,*)' Zero Order Mac Prod = ',Q%MAC%O(0)*P%MAC%O(0)
!       WRITE(*,*)' Zero Order Mac Tens = ', &        
!            (Q%MAC%Delta**(1.D0 + Q%Pole%Ell)*rr22**(-1.D0 - Q%Pole%Ell))/(-Q%MAC%Delta + rr22)

       WRITE(*,*)' PCnt = ',P%Prim%Pw(1)
       WRITE(*,*)' QCnt = ',Q%Box%Center(1)
       WRITE(*,*)' QHlf = ',Q%Box%Half(1)
       WRITE(*,11)'QCENTR = ',Q%POLE%CENTER(1)
       WRITE(*,11)' PQ    = ',PQ(1)
       WRITE(*,*)'-----------------------------------'
       WRITE(*,11)'SPKETC = ',SPKetC(0:LenSP)
       WRITE(*,11)'ERKETC = ',SPErrorKetC(0:LenSP)
       WRITE(*,*)'-----------------------------------'
       WRITE(*,11)'SPKETS = ',SPKetS(0:LenSP)
       WRITE(*,11)'ERKETS = ',SPErrorKetS(0:LenSP)
       WRITE(*,*)'-----------------------------------'
       WRITE(*,11)'QPOLEC = ',Q%Pole%C(1:4)

11     format(A10,12(2x,D12.6))


!    LP=0
    LQ=Q%Pole%Ell
    LPQ=LP+LQ 
    SPKetC=0D0
    SPKetS=0D0
    CALL IrRegular(LPQ,PQ(1),PQ(2),PQ(3))
!    CALL CTraX77P(LP,LQ,SPKetC(0),SPKetS(0),Cpq,Spq,Q%Pole%C(0),Q%Pole%S(0))
    CALL CTraX77(LP,LQ,SPKetC(0),SPKetS(0),Cpq,Spq,Q%Pole%C(0),Q%Pole%S(0))
    
    WRITE(*,11)'SPKETC = ',SPKetC(0:3)
    WRITE(*,11)'SPKETS = ',SPKetS(0:3)

    STOP
ENDIF
    WRITE(*,22)Q%Box%Tier,Q%Box%Number,TotEll,LOG10(TruMACError),LOG10(NewMACError)
22  FORMAT(3(2x,I4),2(2x,F9.3))

    IF(TruMACError>1D2*NewMACError.AND.TruMACError>1D-10.AND.P%Prim%Ell==0)THEN
       !    IF(PRINTQ)THEN
       WRITE(*,36)TruMACError
36     FORMAT(' TruMACError = ',D12.6)
!       STOP 'Well fucked '
    ENDIF

    !          ENDIF
  END SUBROUTINE MACCOMPARE

  RECURSIVE SUBROUTINE MACError(P,Q)
    TYPE(QCPrim)    :: P
    TYPE(PoleNode)  :: Q
    REAL(DOUBLE)    :: Error,LE,RE,PiZ
    REAL(DOUBLE)  :: ERR1,ERR2,PQX,PQY,PQZ,PQ2,RTE,RPE,OMEGA,T,UPQ
    INTEGER       :: I,PQEll,ELL,LP,MP,NP,LQ,MQ,NQ,PDEX,QDEX,LPQ

    REAL(DOUBLE),DIMENSION(3) :: PQ
    REAL(DOUBLE),DIMENSION(0:2*HGEll,0:2*HGEll,0:2*HGEll,0:2*HGEll) :: MDR
    REAL(DOUBLE)                    :: CramCo,MixMax,HGInEq,Err3,PCo,QCo
    REAL(DOUBLE),PARAMETER          :: K3=1.09D0**3
    REAL(DOUBLE),DIMENSION(0:12),PARAMETER :: Fact=(/1D0,1D0,2D0,6D0,24D0,120D0,      &
         720D0,5040D0,40320D0,362880D0,   &
         3628800D0,39916800D0,479001600D0/)

    REAL(DOUBLE), DIMENSION(0:100) :: AuxR,SPTempC,SPTempS
    INTEGER :: LenHG,LenSP

    IF(Q%Leaf)THEN
       Err1=Zero
       Err2=Zero
       DO Ell=0,Q%Herm%Ell
          LenHG=LHGTF(Ell)
          LenSP=LSP(Ell)
          DO I=1,Q%Herm%NQ(Ell)
             PiZ=(Pi/Q%Herm%Zeta(Ell)%D(I))**(1.5D0)

             CALL HGToSP_Direct(Ell,LenHG,LenSP,PiZ,Q%Herm%Coef(Ell)%D(1:LenHG,I), &
                                SPTempC(0:LenSP),SPTempS(0:LenSP))
             LP=P%Prim%Ell
             LQ=Ell
             LPQ=LP+LQ        

!!$             PQ=Q%Herm%Cent(Ell)%D(:,I)-POLECENTER
!!$             DeltaPrime=MAX(DeltaPrime,SQRT(DOT_PRODUCT(PQ,PQ)))

             PQ(1)=P%Prim%Pw(1)-Q%Herm%Cent(Ell)%D(1,I)
             PQ(2)=P%Prim%Pw(2)-Q%Herm%Cent(Ell)%D(2,I)
             PQ(3)=P%Prim%Pw(3)-Q%Herm%Cent(Ell)%D(3,I)
             CALL IrRegular(LPQ,PQ(1),PQ(2),PQ(3))
             CALL CTraX77(LP,LQ,SPErrorKetC(0),SPErrorKetS(0),Cpq,Spq,SPTempC(0:LenSP),SPTempS(0:LenSP))
          ENDDO
       ENDDO
    ELSEIF(ASSOCIATED(Q%Left).AND.ASSOCIATED(Q%Right))THEN
       CALL MACError(P,Q%Left)
       CALL MACError(P,Q%Right)
    ELSEIF(ASSOCIATED(Q%Left))THEN
       CALL MACError(P,Q%Left)
    ELSEIF(ASSOCIATED(Q%Right))THEN
       CALL MACError(P,Q%Right)
    ENDIF
  END SUBROUTINE MACError


  RECURSIVE SUBROUTINE MACError2(P,Q)
    TYPE(QCPrim)    :: P
    TYPE(PoleNode)  :: Q
    REAL(DOUBLE)    :: Error,LE,RE,PiZ
    REAL(DOUBLE)  :: ERR1,ERR2,PQX,PQY,PQZ,PQ2,RTE,RPE,OMEGA,T,UPQ
    INTEGER       :: I,PQEll,ELL,LP,MP,NP,LQ,MQ,NQ,PDEX,QDEX,LPQ

    REAL(DOUBLE),DIMENSION(3) :: PQ
    REAL(DOUBLE),DIMENSION(0:2*HGEll,0:2*HGEll,0:2*HGEll,0:2*HGEll) :: MDR
    REAL(DOUBLE)                    :: CramCo,MixMax,HGInEq,Err3,PCo,QCo
    REAL(DOUBLE),PARAMETER          :: K3=1.09D0**3
    REAL(DOUBLE),DIMENSION(0:12),PARAMETER :: Fact=(/1D0,1D0,2D0,6D0,24D0,120D0,      &
         720D0,5040D0,40320D0,362880D0,   &
         3628800D0,39916800D0,479001600D0/)

    REAL(DOUBLE), DIMENSION(0:100) :: AuxR,SPTempC,SPTempS
    INTEGER :: LenHG,LenSP

    IF(Q%Leaf)THEN
       Err1=Zero
       Err2=Zero
       DO Ell=0,Q%Herm%Ell
          LenHG=LHGTF(Ell)
          LenSP=LSP(Ell)
          DO I=1,Q%Herm%NQ(Ell)
             PiZ=(Pi/Q%Herm%Zeta(Ell)%D(I))**(1.5D0)

             CALL HGToSP_Direct(Ell,LenHG,LenSP,PiZ,Q%Herm%Coef(Ell)%D(1:LenHG,I), &
                                SPTempC(0:LenSP),SPTempS(0:LenSP))
             LP=P%Prim%Ell
             LQ=Ell
             LPQ=LP+LQ        

             PQ=Q%Herm%Cent(Ell)%D(:,I)-POLECENTER

             DeltaPrime=MAX(DeltaPrime,SQRT(DOT_PRODUCT(PQ,PQ)))

             PQ(1)=P%Prim%Pw(1)-Q%Herm%Cent(Ell)%D(1,I)
             PQ(2)=P%Prim%Pw(2)-Q%Herm%Cent(Ell)%D(2,I)
             PQ(3)=P%Prim%Pw(3)-Q%Herm%Cent(Ell)%D(3,I)

             SPErrorketC=0D0

             CALL IrRegular(LPQ,PQ(1),PQ(2),PQ(3))
             CALL CTraX77(LP,LQ,SPErrorKetC(0),SPErrorKetS(0),Cpq,Spq,SPTempC(0:LenSP),SPTempS(0:LenSP))
          ENDDO
       ENDDO
    ELSEIF(ASSOCIATED(Q%Left).AND.ASSOCIATED(Q%Right))THEN
       CALL MACError(P,Q%Left)
       CALL MACError(P,Q%Right)
    ELSEIF(ASSOCIATED(Q%Left))THEN
       CALL MACError(P,Q%Left)
    ELSEIF(ASSOCIATED(Q%Right))THEN
       CALL MACError(P,Q%Right)
    ENDIF
  END SUBROUTINE MACError2

#endif

#ifdef PAC_DEBUG
  !
  RECURSIVE SUBROUTINE ErrCompare(P,Q,Err)
    TYPE(QCPrim)                     :: P
    TYPE(PoleNode)                   :: Q
    INTEGER                          :: TOTELL
    REAL(DOUBLE)                     :: RDIST,SQRTW,LEFTS,RIGHTS,NEWPACERROR,TRUPACERROR,Err,ACh,Ch
    REAL(DOUBLE)                     :: PQ2,PQ,PQxy,PQx,PQy,PQx2,PQy2,PQz,Omega, &
         CoTan,OneOvPQ,OneOvPQxy,RS,SQ,PQToThMnsL, &
         TT,TwoC,COne,SOne,CTwo,STwo,RTE,RPE,T,Upq,PC2,Err1,Err2,Zeta,Eta,Ex
    INTEGER                          :: I,J,Ell,LCode
    LOGICAL                          :: True
    INTEGER                          :: LP,MP,NP,LQ,MQ,NQ,PDex,QDex,PZ,QZ
    TYPE(BBox)                       :: TmpBox

    IF(Q%Box%Number==141)THEN

    TotEll = P%Prim%Ell+Q%Herm%Ell
    NewPACError=Err
    ! Compute the actual PAC error
    IF(ASSOCIATED(Q%Left).AND.ASSOCIATED(Q%Right))THEN
       TruPACError=PACError(P,Q%Left)+PACError(P,Q%Right)
    ELSEIF(ASSOCIATED(Q%Left))THEN
       TruPACError=PACError(P,Q%Left)
    ELSEIF(ASSOCIATED(Q%Right))THEN
       TruPACError=PACError(P,Q%Right)
    ELSE
       TruPACError=PACError(P,Q)
    ENDIF
    IF(ABS(TruPACError)<1D-8)RETURN
    WRITE(*,22)Q%Box%Tier,Q%Box%Number,TotEll,LOG10(ABS(TruPACError)+1D-40),LOG10(NewPACError)
22  FORMAT(3(2x,I4),2(2x,F9.3))
!!$
    IF(ABS(TruPACError)>NewPACError)THEN

       WRITE(*,*)' Leaf = ',Q%Leaf
       WRITE(*,*)' TruE = ',TruPACError
       WRITE(*,*)' NewE = ',NewPACError
       WRITE(*,*)' Tier = ',Q%Box%Tier
       WRITE(*,*)' Numb = ',Q%Box%Number
       IF(Q%Leaf)&
       WRITE(*,*)' Nq   = ',Q%HERM%Nq
       WRITE(*,*)' TotL = ',TotEll
       WRITE(*,*)' PEll = ',P%Prim%Ell
       WRITE(*,*)' --------------------------------------------------------------------'
       WRITE(*,*)' PCnt = ',P%Prim%Pw
       WRITE(*,*)' P Box= ',P%Box%Half
       WRITE(*,*)' P%Zeta ',P%Prim%Zeta
       WRITE(*,*)' --------------------------------------------------------------------'
       WRITE(*,*)' QCnt = ',Q%Box%Center
       WRITE(*,*)' Q Box= ',Q%Box%Half
       WRITE(*,*)' --------------------------------------------------------------------'

       IF(Q%Leaf)THEN
          Ch=0D0
          ACh=0D0
          DO Ell=0,Q%Herm%Ell
             WRITE(*,*)' Ell = ',Ell
             DO I=1,Q%Herm%NQ(Ell)
                IF(Q%Herm%Zeta(Ell)%D(I)<1D10)THEN
                   PQx=P%Prim%Pw(1)-Q%Herm%Cent(Ell)%D(1,I)
                   PQy=P%Prim%Pw(2)-Q%Herm%Cent(Ell)%D(2,I)
                   PQz=P%Prim%Pw(3)-Q%Herm%Cent(Ell)%D(3,I)
                   PQ2=PQx*PQx+PQy*PQy+PQz*PQz
                   RTE=P%Prim%Zeta*Q%Herm%Zeta(Ell)%D(I)
                   RPE=P%Prim%Zeta+Q%Herm%Zeta(Ell)%D(I)
                   Omega=RTE/RPE 
                   T=Omega*PQ2
                   IF(T<1D0)THEN
                      TmpBox%BndBox(:,1)=Q%Herm%Cent(Ell)%D(:,I)
                      TmpBox%BndBox(:,2)=Q%Herm%Cent(Ell)%D(:,I)
                      Ex=Extent(Ell,Q%Herm%Zeta(Ell)%D(I),Q%Herm%Coef(Ell)%D(:,I),Tau_O=TauPAC,Potential_O=.TRUE.,ExtraEll_O=0)
                      TmpBox=ExpandBox(TmpBox,Ex)
                      
                      WRITE(*,*)'1 BoxOutSideBox?? ',BoxOutSideBox(Q%Box,P%Box) 
                      WRITE(*,*)'2 BoxOutSideBox?? ',BoxOutSideBox(TmpBox,P%Box) 
                      


                      WRITE(*,*)' Zeta = ',Q%Herm%Zeta(Ell)%D(I),' T = ',T                      

                      WRITE(*,*)' QBox = '
                      CALL PrintBBox(Q%BOX,6)
                      WRITE(*,*)' TmpBox = '
                      CALL PrintBBox(TmpBox,6)
                      CALL Halt(' In FillRhoLeaf: Multipole center outside of BBox ')
                   ENDIF

!                   WRITE(*,*)PQ2
                   Ch=Ch+Q%Herm%Coef(Ell)%D(1,I)*(Pi/Q%Herm%Zeta(Ell)%D(I))**(3D0/2D0)
                   ACh=ACh+ABS(Q%Herm%Coef(Ell)%D(1,I)*(Pi/Q%Herm%Zeta(Ell)%D(I))**(3D0/2D0))
                END IF
             ENDDO
          ENDDO          
          WRITE(*,*)' Total Charge = ',Ch
          WRITE(*,*)' Total Charge = ',ACh
       ENDIF
       STOP
    ENDIF


    ENDIF

  END SUBROUTINE ERRCOMPARE

  RECURSIVE FUNCTION PACError(P,Q) RESULT(Error)
    TYPE(QCPrim)    :: P
    TYPE(PoleNode)  :: Q
    REAL(DOUBLE)    :: Error,LE,RE
    REAL(DOUBLE)  :: ERR1,ERR2,PQX,PQY,PQZ,PQ2,RTE,RPE,OMEGA,T,UPQ
    INTEGER       :: I,PQEll,ELL,LP,MP,NP,LQ,MQ,NQ,PDEX,QDEX

    REAL(DOUBLE),DIMENSION(0:2*HGEll,0:2*HGEll,0:2*HGEll,0:2*HGEll) :: MDR
    REAL(DOUBLE)                    :: CramCo,MixMax,HGInEq,Err3,PCo,QCo
    REAL(DOUBLE),PARAMETER          :: K3=1.09D0**3
    REAL(DOUBLE),DIMENSION(0:12),PARAMETER :: Fact=(/1D0,1D0,2D0,6D0,24D0,120D0,      &
         720D0,5040D0,40320D0,362880D0,   &
         3628800D0,39916800D0,479001600D0/)

    REAL(DOUBLE), DIMENSION(0:100) :: AuxR

    IF(Q%Leaf)THEN
       Err1=Zero
       Err2=Zero
       DO Ell=0,Q%Herm%Ell
          DO I=1,Q%Herm%NQ(Ell)
             IF(Q%Herm%Zeta(Ell)%D(I)<1D10)THEN
                PQx=P%Prim%Pw(1)-Q%Herm%Cent(Ell)%D(1,I)
                PQy=P%Prim%Pw(2)-Q%Herm%Cent(Ell)%D(2,I)
                PQz=P%Prim%Pw(3)-Q%Herm%Cent(Ell)%D(3,I)
                PQ2=PQx*PQx+PQy*PQy+PQz*PQz
                RTE=P%Prim%Zeta*Q%Herm%Zeta(Ell)%D(I)
                RPE=P%Prim%Zeta+Q%Herm%Zeta(Ell)%D(I)
                Omega=RTE/RPE 
                T=Omega*PQ2
                PQx=-PQx; PQy=-PQy; PQz=-PQz
                Upq=TwoPi5x2/(RTE*SQRT(RPE))
                PQEll=P%Prim%Ell+Ell
                CALL ErrInts(2*HGEll,PQEll,AuxR,Omega,T)
                CALL MD3TRR(2*HGEll,PQEll,MDR,AuxR,Upq,PQx,PQy,PQz)
                DO LP=0,P%Prim%Ell
                   DO MP=0,P%Prim%Ell-LP
                      DO NP=0,P%Prim%Ell-LP-MP 
                         PDex=LMNDex(LP,MP,NP)                   
                         DO LQ=0,Ell
                            DO MQ=0,Ell-LQ
                               DO NQ=0,Ell-LQ-MQ
                                  QDex=LMNDex(LQ,MQ,NQ)
                                  Err1=Err1+ErrBra(PDex)*MDR(LP+LQ,MP+MQ,NP+NQ,0)*Q%Herm%Coef(Ell)%D(QDex,I)
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO
       Error=Err1


       WRITE(*,*)' Err1= ',Err1

    ELSEIF(ASSOCIATED(Q%Left).AND.ASSOCIATED(Q%Right))THEN
       Error=PACError(P,Q%Left)+PACError(P,Q%Right)
    ELSEIF(ASSOCIATED(Q%Left))THEN
       Error=PACError(P,Q%Left)
    ELSEIF(ASSOCIATED(Q%Right))THEN
       Error=PACError(P,Q%Right)
    ENDIF
  END FUNCTION PACError

#endif

  SUBROUTINE VMD3(EllP,LenP,EllQ,LenQ,EllPQ,LenPQ,Nc,P,ZetaP,Q,ZetaQ,QCoeff,Ket)

    USE GammaFunctions
    !
    INTEGER :: EllP,LenP,EllQ,LenQ,EllPQ,LenPQ,Nc
    REAL(DOUBLE),DIMENSION(Nc,0:LenPQ)    :: R
    REAL(DOUBLE),DIMENSION(Nc,0:EllPQ)   :: AuxR
    REAL(DOUBLE),DIMENSION(3)            :: P
    REAL(DOUBLE),DIMENSION(3,Nc)         :: Q
    REAL(DOUBLE),DIMENSION(Nc)           :: ZetaQ
    REAL(DOUBLE),DIMENSION(Nc)           :: PQx,PQy,PQz
    REAL(DOUBLE),DIMENSION(LenQ,Nc)      :: QCoeff
    REAL(DOUBLE),DIMENSION(LenP)         :: Ket
    !
    REAL(DOUBLE)                  :: Omega,T,PQ2
    INTEGER                       :: MaxL,LTot 
    REAL(DOUBLE),PARAMETER        :: Switch=26.0D0
    INTEGER,PARAMETER             :: LPlus= 300
    INTEGER,PARAMETER             :: L2=12+LPlus
    REAL(DOUBLE),DIMENSION(0:L2)  :: F
    REAL(DOUBLE)                  :: SqrtT,ET,OneOvT,FJ,TwoT,OmegaJ,TwoO
    INTEGER                       :: I,J,K,L,M,N,LP,MP,NP,LQ,MQ,NQ


    REAL(DOUBLE) :: ZetaP,RTE,RPE,Upq,FLT1 !,Mm1,Nm1,Lm1

    INTEGER :: IDX1,IDX2,IDX3,IP,IQ,IPQ


    !    INTEGER:: Jp0,Jp1,Np0_Jp0,Nm1_Jp0,Nm2_Jp0,Nm1_Jp1,Nm2_Jp1, &
    !         Meq1_Jp0,Meq0_Jp0,Meq0_Jp1,Mp0_Jp0,Mm1_Jp1,Mm2_Jp1, &
    !         Leq1_Jp0,Leq0_Jp1,Lp0_Jp0,Lm1_Jp1,Lm2_Jp1

    INTEGER :: Pdex,Qdex,PQdex


    REAL(DOUBLE),DIMENSION(LenP) :: Ket1



    !-------------------------------------------------------------------------
!!$    WRITE(*,*)' EllP = ',EllP,' EllQ = ',EllQ
!!$    WRITE(*,*)' EllPQ = ',EllPQ,' LenPQ = ',LenPQ

    DO I=1,Nc
       PQx(I)=-(P(1)-Q(1,I))
       PQy(I)=-(P(2)-Q(2,I))
       PQz(I)=-(P(3)-Q(3,I))
       PQ2=PQx(I)*PQx(I)+PQy(I)*PQy(I)+PQz(I)*PQz(I)

       IF(ZetaP==NuclearExpnt.AND.EllQ==0.AND.(PQ2.LT.1D-16.AND.ZetaQ(I).EQ.ZetaP))THEN
          AuxR(I,:)=0D0
       ELSE
          RTE=ZetaP*ZetaQ(I)
          RPE=ZetaP+ZetaQ(I)
          Omega=RTE/RPE
          Upq=TwoPi5x2/(RTE*SQRT(RPE))
          T=Omega*PQ2

          
!          WRITE(*,*)' T = ',T
!          WRITE(*,*)' Q = ',Q(1,I),' T = ',T

          IF(T==Zero)THEN
             OmegaJ=One
             TwoO=-Two*Omega
             DO J=0,EllPQ
                AuxR(I,J)=Upq*OmegaJ/DBLE(2*J+1)
                OmegaJ=TwoO*OmegaJ
             ENDDO
          ELSE
             IF(T.LT.Switch) THEN
                ! F_{j}(T)=(2*T*F_{j+1}+Exp[-T])/(2*j+1)
                !             WRITE(*,*)' T = ',T
                TwoT=Two*T
                ET=DEXP(-T)
                !             WRITE(*,*)' ET = ',ET
                FJ=Zero
                DO J=EllPQ+LPlus,0,-1                
                   F(J)=FJ
                   !                WRITE(*,*)' J = ',J,' FJ = ',F(J)
                   FJ=(TwoT*F(J)+ET)/(Two*DBLE(J)-One)
                ENDDO
             ELSE
                ! Multipole approx and upward recursion
                SqrtT=SQRT(T)
                OneOvT=One/T
                F(0)=SqrtPi/(Two*SqrtT) 
                DO J=1,EllPQ
                   F(J)=F(J-1)*(DBLE(J)-Half)*OneOvT
                ENDDO
             ENDIF
             OmegaJ=One
             TwoO=-Two*Omega
             DO J=0,EllPQ
                AuxR(I,J)=Upq*OmegaJ*F(J)
                OmegaJ=TwoO*OmegaJ
             ENDDO
          ENDIF
       ENDIF

!!$
!       WRITE(*,333)ZetaP,ZetaQ(I),T,F(0)
333    FORMAT('ZP = ',D12.6,' ZQ = ',D12.6,' T = ',D12.6,' F = ',D12.6)
    ENDDO


    !
    DO J=EllPQ,0,-1
       DO L=EllPQ-J,1,-1
          IDX1=LMNx(L,0,0)
          IDX2=LMNx(L-1,0,0)
          IDX3=LMNx(L-2,0,0)
          Flt1=DBLE(L-1)
          DO I=1,Nc
             R(I,IDX1)=PQx(I)*R(I,IDX2)+Flt1*R(I,IDX3)
          ENDDO
          DO M=EllPQ-J-L,1,-1
             IDX1=LMNx(M,L,0)
             IDX2=LMNx(M-1,L,0)
             IDX3=LMNx(M-2,L,0)
             FLT1=DBLE(M-1)
             DO I=1,Nc
                R(I,IDX1)=PQx(I)*R(I,IDX2)+Flt1*R(I,IDX3)
             ENDDO
             DO N=EllPQ-J-L-M,1,-1
                IDX1=LMNx(N,M,L)
                IDX2=LMNx(N-1,M,L)
                IDX3=LMNx(N-2,M,L)
                FLT1=DBLE(N-1)
                DO I=1,Nc
                   R(I,IDX1)=PQx(I)*R(I,IDX2)+Flt1*R(I,IDX3)
                ENDDO
             ENDDO
             IDX1=LMNx(0,M,L)
             IDX2=LMNx(0,M-1,L)
             IDX3=LMNx(0,M-2,L)
             FLT1=DBLE(M-1)
             DO I=1,Nc
                R(I,IDX1)=PQy(I)*R(I,IDX2)+Flt1*R(I,IDX3)
             ENDDO
             IDX1=LMNx(M,0,L)
             IDX2=LMNx(M-1,0,L)
             IDX3=LMNx(M-2,0,L)
             FLT1=DBLE(M-1)
             DO I=1,Nc
                R(I,IDX1)=PQx(I)*R(I,IDX2)+Flt1*R(I,IDX3)
             ENDDO
          ENDDO
          IDX1=LMNx(0,L,0)
          IDX2=LMNx(0,L-1,0)
          IDX3=LMNx(0,L-2,0)
          FLT1=DBLE(L-1)
          DO I=1,Nc
             R(I,IDX1)=PQy(I)*R(I,IDX2)+Flt1*R(I,IDX3)
          ENDDO
          IDX1=LMNx(0,0,L)
          IDX2=LMNx(0,0,L-1)
          IDX3=LMNx(0,0,L-2)
          FLT1=DBLE(L-1)
          DO I=1,Nc
             R(I,IDX1)=PQz(I)*R(I,IDX2)+Flt1*R(I,IDX3)
          ENDDO
       ENDDO
       DO I=1,Nc
          R(I,1)=AuxR(I,J)
       ENDDO
    ENDDO


    DO I=1,LenPQ
       WRITE(*,*)I,R(1,I)
    ENDDO

    STOP




!!$    CALL VMD3RR(MaxCluster,HGEll4,Nc,EllPQ,LMNx(-1,-1,-1),AuxR(1,0),PQx(1),PQy(1),PQz(1),R3T(1,1))
!!$    CALL KetKontract(MaxCluster,Nc,LenP*LenQ,PLMNx(EllP,EllQ)%I(1),QLMNx(EllP,EllQ)%I(1), &
!!$         PQLMNx(EllP,EllQ)%I(1),R3T(1,1),QCoeff(1,1),Ket(1))
!!$    LenPQ=LenP*LenQ
!!$    DO J=1,LenPQ
!!$       IP=PLMNx(EllP,EllQ)%I(J)
!!$       IQ=QLMNx(EllP,EllQ)%I(J)
!!$       IPQ=PQLMNx(EllP,EllQ)%I(J)
!!$       DO I=1,Nc
!!$            WRITE(*,11)J,IP,IQ,IPQ,R(I,IPQ),QCoeff(I,IQ) 
!!$ 11         FORMAT(4(I4,", "),2(D12.6,", "))
!!$!          Ket(IP)=Ket(IP)+R(I,IPQ)*QCoeff(I,IQ)
!!$       ENDDO
!!$    ENDDO



!!$    DO I=1,LenP
!!$       IF(ABS( Ket(I)-Ket1(I) )>1D-8 )THEN
!!$          WRITE(*,*)' I = ',I,Ket(I),KEt1(I)
!!$          STOP
!!$       ENDIF
!!$    ENDDO
    !

    !    DO I=1,Nc
    !       WRITE(*,*)I,DOT_PRODUCT(QCoeff(:,I),QCoeff(:,I))
    !    ENDDO

    DO LP=0,EllP
       DO MP=0,EllP-LP
          DO NP=0,EllP-LP-MP 
             Pdex=LMNx(LP,MP,NP)
             DO LQ=0,EllQ
                DO MQ=0,EllQ-LQ
                   DO NQ=0,EllQ-LQ-MQ
                      Qdex=LMNx(LQ,MQ,NQ)
                      PQdex=LMNx(LP+LQ,MP+MQ,NP+NQ)
                      DO I=1,Nc
!                         IF(QCoeff(QDex,I).NE.0D0.AND.Q(1,I)==0D0)THEN
                             Ket(PDex)=Ket(Pdex)+R(I,PQdex)*QCoeff(Qdex,I)
!                             WRITE(*,*)' Ket = ',Ket(1)
!                             WRITE(*,444)I,P(1),Q(1,I),R(I,PQdex)*QCoeff(Qdex,I)
444                          FORMAT(I4,' px = ',D12.6,', qx = ',D12.6,', R*Q = ',D12.6)
!                          ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE VMD3



END MODULE TreeWalk


