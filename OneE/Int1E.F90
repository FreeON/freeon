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
MODULE Int1E
!=================================================================================
!H MODULE Int1E
!H  Author:
!H  V. Weber
!H
!H This MODULE contains:
!H  PUBLIC:
!H  o SUB
!H
!H  PRIVATE:
!H  o SUB
!H
!H  OPTIONS:
!H
!H  THRESHOLD: Use -DINT1E_THRESH to threshold the primitives.
#define INT1E_THRESH
!H  DEBUGING : Use -DINT1E_DBUG to print some stuff.
!H  INFO     : Use -DINT1E_INFO to print some stuff.
!H
!H  Comments:
!H
!---------------------------------------------------------------------------------
  !
  USE DerivedTypes
  USE GlobalScalars
  USE AtomPairs
  IMPLICIT NONE
  !
  PRIVATE
!---------------------------------------------------------------------------------
! PUBLIC DECLARATIONS
!---------------------------------------------------------------------------------
  PUBLIC  :: I1E_Ovlap_AtBlk,I1E_Kinet_AtBlk,I1E_EDipl_AtBlk,I1E_EQuad_AtBlk
  PUBLIC  :: I1E_Ovlap1C_AtBlk,I1E_Kinet1C_AtBlk,I1E_EDipl1C_AtBlk
  PUBLIC  :: I1E_Ovlap2C_AtBlk,I1E_Kinet2C_AtBlk
  PUBLIC  :: I1E_KinMom_AtBlk,I1E_DipVel_AtBlk,I1E_rKinet_AtBlk
  PUBLIC  :: I1E_MssVel_AtBlk
  PUBLIC  :: I1E_NucAtt_AtBlk
  !
!---------------------------------------------------------------------------------
! PRIVATE DECLARATIONS
!---------------------------------------------------------------------------------
  PRIVATE :: I1E_Ovlap_OSRecurs,I1E_NucAtt_OSRecurs
  !
!---------------------------------------------------------------------------------
! SEMI-GLOBAL PARAMETER
!---------------------------------------------------------------------------------
  INTEGER, PARAMETER, PRIVATE :: I1E_MAXL=5 !MAXL is 1=S,2=P,3=D,4=F,5=G,...
  !
!---------------------------------------------------------------------------------
! SEMI-GLOBAL VARIABLES
!---------------------------------------------------------------------------------
  REAL(DOUBLE), DIMENSION(0:I1E_MAXL+4,0:I1E_MAXL+4), PRIVATE :: RRX,RRY,RRZ
  REAL(DOUBLE), DIMENSION(0:2*I1E_MAXL,0:(I1E_MAXL+1)**3,0:(I1E_MAXL+1)**3), PRIVATE :: RR0
  !
CONTAINS
  !
  !
  SUBROUTINE I1E_Ovlap_OSRecurs(PA,PB,Gamma,LMaxI,LMaxJ)
!---------------------------------------------------------------------------------
!H Compute OS recurrence relation:
!H $$
!H   (a+1_i|b)=PA_i(a|b)+\frac{1}{\zeta}N_i(a)(a-1_i|b)
!H    +\frac{1}{\zeta}N_i(b)(a|b-1_i)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    INTEGER                    :: LMaxI,LMaxJ
    REAL(DOUBLE)               :: Gamma
    REAL(DOUBLE), DIMENSION(3) :: PA,PB
    !-------------------------------------------------------------------
    INTEGER                    :: I,IM1,IP1,J,JM1,JP1
    REAL(DOUBLE)               :: G
    !-------------------------------------------------------------------
    !
    G=1.0D0/(2.0D0*Gamma)
    !
    ! Initialize.
    !
    RRX(0,0)=1.0D0
    RRY(0,0)=1.0D0
    RRZ(0,0)=1.0D0
    !
    ! Upward recursion in I for J=0.
    !
    RRX(1,0)=PA(1)
    RRY(1,0)=PA(2)
    RRZ(1,0)=PA(3)
    !
    DO I=1,LMaxI-1
       IP1=I+1
       IM1=I-1
       RRX(IP1,0)=DBLE(I)*G*RRX(IM1,0)+PA(1)*RRX(I,0)
       RRY(IP1,0)=DBLE(I)*G*RRY(IM1,0)+PA(2)*RRY(I,0)
       RRZ(IP1,0)=DBLE(I)*G*RRZ(IM1,0)+PA(3)*RRZ(I,0)
    ENDDO
    !
    ! Upward recursion in J for all I's .
    !
    RRX(0,1)=PB(1)
    RRY(0,1)=PB(2)
    RRZ(0,1)=PB(3)
    !
    DO I=1,LMaxI
       IM1=I-1
       RRX(I,1)=DBLE(I)*G*RRX(IM1,0)+PB(1)*RRX(I,0)
       RRY(I,1)=DBLE(I)*G*RRY(IM1,0)+PB(2)*RRY(I,0)
       RRZ(I,1)=DBLE(I)*G*RRZ(IM1,0)+PB(3)*RRZ(I,0)
    ENDDO
    !
    DO J=1,LMaxJ-1
       JP1=J+1
       JM1=J-1
       RRX(0,JP1)=DBLE(J)*G*RRX(0,JM1)+PB(1)*RRX(0,J)
       RRY(0,JP1)=DBLE(J)*G*RRY(0,JM1)+PB(2)*RRY(0,J)
       RRZ(0,JP1)=DBLE(J)*G*RRZ(0,JM1)+PB(3)*RRZ(0,J)
       DO I=1,LMaxI
          IM1=I-1
          RRX(I,JP1)=G*(DBLE(I)*RRX(IM1,J)+DBLE(J)*RRX(I,JM1))+PB(1)*RRX(I,J)
          RRY(I,JP1)=G*(DBLE(I)*RRY(IM1,J)+DBLE(J)*RRY(I,JM1))+PB(2)*RRY(I,J)
          RRZ(I,JP1)=G*(DBLE(I)*RRZ(IM1,J)+DBLE(J)*RRZ(I,JM1))+PB(3)*RRZ(I,J)
       ENDDO
    ENDDO
    !
  END SUBROUTINE I1E_Ovlap_OSRecurs
  !
  !
  SUBROUTINE I1E_NucAtt_OSRecurs(PA,PB,PC,Gamma,LMaxI,LMaxJ)
!---------------------------------------------------------------------------------
!H Compute OS recurrence relation:
!H $$
!H   (a+1_i|...|b)=PA_i(a|b)+\frac{1}{\zeta}N_i(a)(a-1_i|b)
!H    +\frac{1}{\zeta}N_i(b)(a|b-1_i)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    INTEGER                    :: LMaxI,LMaxJ
    REAL(DOUBLE)               :: Gamma
    REAL(DOUBLE), DIMENSION(3) :: PA,PB,PC
    !-------------------------------------------------------------------
    INTEGER                    :: A,B,M,OAX,OAY,OAZ,OBX,OBY,OBZ
    INTEGER                    :: MMax,IX,IY,IZ,JX,JY,JZ,IdA,IdB
    REAL(DOUBLE)               :: G,T,U
    REAL(DOUBLE), PARAMETER    :: InvSqrtPi=5.6418958354775628695D-1
    !-------------------------------------------------------------------
    !
    OAZ=1
    OAY=LMaxI+1
    OAX=OAY*OAY
    OBZ=1
    OBY=LMaxJ+1
    OBX=OBY*OBY
    MMax=LMaxI+LMaxJ
    G=1.0D0/(2.0D0*Gamma)
    T=2D0*SQRT(Gamma)*InvSqrtPi
    U=Gamma*(PC(1)*PC(1)+PC(2)*PC(2)+PC(3)*PC(3))
    !
    ! Initialize RR.
    !
    DO M=0,MMax
       RR0(M,0,0)=T*GammaF(M,U)
    ENDDO
    !
    ! Upward recursion in I with J=0.
    !
    DO A=1,LMaxI
       DO IX=0,A
       DO IY=0,A-IX
          IZ=A-IX-IY
          IdA=IX*OAX+IY*OAY+IZ*OAZ
          IF(IZ.GT.0) THEN
             DO M=0,MMax-A
                RR0(M,IdA,0)=PA(3)*RR0(M,IdA-OAZ,0)-PC(3)*RR0(M+1,IdA-OAZ,0)
             ENDDO
             IF(IZ.GT.1) THEN
                DO M=0,MMax-A
                   RR0(M,IdA,0)=RR0(M,IdA,0)+ &
                          G*DBLE(IZ-1)*(RR0(M,IdA-2*OAZ,0)-RR0(M+1,IdA-2*OAZ,0))
                ENDDO
             ENDIF
          ELSEIF(IY.GT.0) THEN
             DO M=0,MMax-A
                RR0(M,IdA,0)=PA(2)*RR0(M,IdA-OAY,0)-PC(2)*RR0(M+1,IdA-OAY,0)
             ENDDO
             IF(IY.GT.1) THEN
                DO M=0,MMax-A
                   RR0(M,IdA,0)=RR0(M,IdA,0)+ &
                          G*DBLE(IY-1)*(RR0(M,IdA-2*OAY,0)-RR0(M+1,IdA-2*OAY,0))
                ENDDO
             ENDIF
          ELSEIF(IX.GT.0) THEN
             DO M=0,MMax-A
                RR0(M,IdA,0)=PA(1)*RR0(M,IdA-OAX,0)-PC(1)*RR0(M+1,IdA-OAX,0)
             ENDDO
             IF(IX.GT.1) THEN
                DO M=0,MMax-A
                   RR0(M,IdA,0)=RR0(M,IdA,0)+ &
                          G*DBLE(IX-1)*(RR0(M,IdA-2*OAX,0)-RR0(M+1,IdA-2*OAX,0))
                ENDDO
             ENDIF
          ELSE
             STOP 9991
          ENDIF
       ENDDO
       ENDDO
    ENDDO
    !
    ! Upward recursion in J with all possible I's.
    !
    DO A=0,LMaxI
       DO IX=0,A
       DO IY=0,A-IX
          IZ=A-IX-IY
          IdA=IX*OAX+IY*OAY+IZ*OAZ
          DO B=1,LMaxJ
             DO JX=0,B
             DO JY=0,B-JX
                JZ=B-JX-JY
                IdB=JX*OBX+JY*OBY+JZ*OBZ
                IF(JZ.GT.0) THEN
                   DO M=0,MMax-A-B
                      RR0(M,IdA,IdB)=PB(3)*RR0(M,IdA,IdB-OBZ)-PC(3)*RR0(M+1,IdA,IdB-OBZ)
                   ENDDO
                   IF(JZ.GT.1) THEN
                      DO M=0,MMax-A-B
                         RR0(M,IdA,IdB)=RR0(M,IdA,IdB)+&
                                G*DBLE(JZ-1)*(RR0(M,IdA,IdB-2*OBZ)-RR0(M+1,IdA,IdB-2*OBZ))
                      ENDDO
                   ENDIF
                   IF(IZ.GT.0) THEN
                      DO M=0,MMax-A-B
                         RR0(M,IdA,IdB)=RR0(M,IdA,IdB)+&
                                G*DBLE(IZ)*(RR0(M,IdA-OAZ,IdB-OBZ)-RR0(M+1,IdA-OAZ,IdB-OBZ))
                      ENDDO
                   ENDIF
                ELSEIF(JY.GT.0) THEN
                   DO M=0,MMax-A-B
                      RR0(M,IdA,IdB)=PB(2)*RR0(M,IdA,IdB-OBY)-PC(2)*RR0(M+1,IdA,IdB-OBY)
                   ENDDO
                   IF(JY.GT.1) THEN
                      DO M=0,MMax-A-B
                         RR0(M,IdA,IdB)=RR0(M,IdA,IdB)+&
                                G*DBLE(JY-1)*(RR0(M,IdA,IdB-2*OBY)-RR0(M+1,IdA,IdB-2*OBY))
                      ENDDO
                   ENDIF
                   IF(IY.GT.0) THEN
                      DO M=0,MMax-A-B
                         RR0(M,IdA,IdB)=RR0(M,IdA,IdB)+&
                                G*DBLE(IY)*(RR0(M,IdA-OAY,IdB-OBY)-RR0(M+1,IdA-OAY,IdB-OBY))
                      ENDDO
                   ENDIF
                ELSEIF(JX.GT.0) THEN
                   DO M=0,MMax-A-B
                      RR0(M,IdA,IdB)=PB(1)*RR0(M,IdA,IdB-OBX)-PC(1)*RR0(M+1,IdA,IdB-OBX)
                   ENDDO
                   IF(JX.GT.1) THEN
                      DO M=0,MMax-A-B
                         RR0(M,IdA,IdB)=RR0(M,IdA,IdB)+&
                                G*DBLE(JX-1)*(RR0(M,IdA,IdB-2*OBX)-RR0(M+1,IdA,IdB-2*OBX))
                      ENDDO
                   ENDIF
                   IF(IX.GT.0) THEN
                      DO M=0,MMax-A-B
                         RR0(M,IdA,IdB)=RR0(M,IdA,IdB)+&
                                G*DBLE(IX)*(RR0(M,IdA-OAX,IdB-OBX)-RR0(M+1,IdA-OAX,IdB-OBX))
                      ENDDO
                   ENDIF
                ELSE
                   STOP 1234
                ENDIF
             ENDDO
             ENDDO
          ENDDO
       ENDDO
       ENDDO
    ENDDO
    !
  END SUBROUTINE I1E_NucAtt_OSRecurs
  !
  !
  SUBROUTINE I1E_Ovlap_AtBlk(BS,Pair,I1E)
!---------------------------------------------------------------------------------
!H Compute the overlap matrix:
!H $$
!H   (a|b)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BSET)                               :: BS
    TYPE(AtomPair)                           :: Pair
    REAL(DOUBLE), DIMENSION(Pair%NA,Pair%NB) :: I1E
    !-------------------------------------------------------------------
    INTEGER :: CFA,CFB,IndexA,IndexB,StartLA,StartLB
    INTEGER :: StopLA,StopLB,MaxLA,MaxLB,LMNA,LMNB
    INTEGER :: KA,KB,PFA,PFB,IA,IB,LA,LB,MA,MB,NA,NB
    REAL(DOUBLE), DIMENSION(3) :: P,A,B,PA,PB
    REAL(DOUBLE) :: AB2,ZA,ZB,CA,CB,Gamma,InvGa,Dum,Ov
    !-------------------------------------------------------------------
    !
    KA=Pair%KA
    KB=Pair%KB
    A(1)=Pair%A(1)
    A(2)=Pair%A(2)
    A(3)=Pair%A(3)
    B(1)=Pair%B(1)
    B(2)=Pair%B(2)
    B(3)=Pair%B(3)
    AB2=Pair%AB2
    !
    ! Loop over contracted functions A and B.
    !
    IndexA=0
    DO CFA=1,BS%NCFnc%I(KA)
       StartLA = BS%LStrt%I(CFA,KA)
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       IndexB=0
       DO CFB=1,BS%NCFnc%I(KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          !
          ! Begin loops over primitives.
          !
          DO PFA=BS%NPFnc%I(CFA,KA),1,-1
             ZA=BS%Expnt%D(PFA,CFA,KA)
          DO PFB=BS%NPFnc%I(CFB,KB),1,-1
             ZB=BS%Expnt%D(PFB,CFB,KB)
             !
             Gamma=ZA+ZB
             InvGa=1.0D0/Gamma
             Dum=ZA*ZB*AB2*InvGa
#ifdef INT1E_THRESH
             IF(Dum>PrimPairDistanceThreshold) EXIT
#endif
             Ov=EXP(-Dum)*SQRT(PI*InvGa)*PI*InvGa
             !
             P(1)=(ZA*A(1)+ZB*B(1))*InvGa
             P(2)=(ZA*A(2)+ZB*B(2))*InvGa
             P(3)=(ZA*A(3)+ZB*B(3))*InvGa
             PA(1)=P(1)-A(1)
             PA(2)=P(2)-A(2)
             PA(3)=P(3)-A(3)
             PB(1)=P(1)-B(1)
             PB(2)=P(2)-B(2)
             PB(3)=P(3)-B(3)
             !
             ! Compute OS RR.
             !
             CALL I1E_Ovlap_OSRecurs(PA(1),PB(1),Gamma,MaxLA,MaxLB)
             !
             ! Contract each buffer into appropriate location.
             !
             IB=IndexB
             DO LMNB=StartLB,StopLB
                LB=BS%LxDex%I(LMNB)
                MB=BS%LyDex%I(LMNB)
                NB=BS%LzDex%I(LMNB)
                CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                IB=IB+1
             IA=IndexA
             DO LMNA=StartLA,StopLA
                LA=BS%LxDex%I(LMNA)
                MA=BS%LyDex%I(LMNA)
                NA=BS%LzDex%I(LMNA)
                CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                IA=IA+1
                !
                I1E(IA,IB)=I1E(IA,IB)+RRX(LA,LB)*RRY(MA,MB)*RRZ(NA,NB)*Ov*CA*CB
                !
             ENDDO
             ENDDO
             !
          ENDDO
#ifdef INT1E_THRESH
          IF(PFB.EQ.BS%NPFnc%I(CFB,KB)) EXIT
#endif
          ENDDO
          !
          IndexB=IndexB+StopLB-StartLB+1
       ENDDO
       IndexA=IndexA+StopLA-StartLA+1
    ENDDO
    !
  END SUBROUTINE I1E_Ovlap_AtBlk
  !
  !
  REAL(DOUBLE) FUNCTION I1E_Ovlap_Max(A,ZA,MaxLA,B,ZB,MaxLB)
!---------------------------------------------------------------------------------
!H Compute the overlap matrix:
!H $$
!H   (a|b)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    INTEGER                    :: MaxLA,MaxLB
    REAL(DOUBLE)               :: ZA,ZB
    REAL(DOUBLE), DIMENSION(3) :: A,B
    !-------------------------------------------------------------------
    INTEGER                    :: IA,JA,IB,JB,LA,LB,MA,MB,NA,NB
    REAL(DOUBLE), DIMENSION(3) :: P,PA,PB
    REAL(DOUBLE)               :: AB2,Gamma,InvGa,Dum,Ov
    !-------------------------------------------------------------------
    !
    AB2=(A(1)-B(1))**2+(A(2)-B(2))**2+(A(3)-B(3))**2
    Gamma=ZA+ZB
    InvGa=1.0D0/Gamma
    Dum=ZA*ZB*AB2*InvGa
    Ov=EXP(-Dum)*SQRT(PI*InvGa)*PI*InvGa
    !
    P(1)=(ZA*A(1)+ZB*B(1))*InvGa
    P(2)=(ZA*A(2)+ZB*B(2))*InvGa
    P(3)=(ZA*A(3)+ZB*B(3))*InvGa
    PA(1)=P(1)-A(1)
    PA(2)=P(2)-A(2)
    PA(3)=P(3)-A(3)
    PB(1)=P(1)-B(1)
    PB(2)=P(2)-B(2)
    PB(3)=P(3)-B(3)
    !
    ! Compute OS RR.
    !
    CALL I1E_Ovlap_OSRecurs(PA(1),PB(1),Gamma,MaxLA,MaxLB)
    !
    I1E_Ovlap_Max=0D0
    !
    DO IA=0,MaxLA
       LA=MaxLA-IA
       DO JA=0,IA
          MA=IA-JA
          NA=JA
          DO IB=0,MaxLB
             LB=MaxLB-IB
             DO JB=0,IB
                MB=IB-JB
                NB=JB
                I1E_Ovlap_Max=MAX(I1E_Ovlap_Max, &
                       ABS(RRX(LA,LB)*RRY(MA,MB)*RRZ(NA,NB)*Ov)) !Add some normalization here!
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
  END FUNCTION I1E_Ovlap_Max
  !
  !
  SUBROUTINE I1E_Kinet_AtBlk(BS,Pair,I1E)
!---------------------------------------------------------------------------------
!H Compute the kinetic matrix:
!H $$
!H   -\frac{1}{2}(a|\nabla^2|b)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BSET)                               :: BS
    TYPE(AtomPair)                           :: Pair
    REAL(DOUBLE), DIMENSION(Pair%NA,Pair%NB) :: I1E
    !-------------------------------------------------------------------
    INTEGER :: CFA,CFB,IndexA,IndexB,StartLA,StartLB
    INTEGER :: StopLA,StopLB,MaxLA,MaxLB,LMNA,LMNB
    INTEGER :: KA,KB,PFA,PFB,IA,IB,LA,LB,MA,MB,NA,NB
    REAL(DOUBLE), DIMENSION(3) :: P,A,B,PA,PB
    REAL(DOUBLE) :: AB2,ZA,ZB,CA,CB,Gamma,InvGa,Dum,Ov,OvT
    REAL(DOUBLE) :: X0,Y0,Z0,TX,TY,TZ
    !-------------------------------------------------------------------
    !
    KA=Pair%KA
    KB=Pair%KB
    A(1)=Pair%A(1)
    A(2)=Pair%A(2)
    A(3)=Pair%A(3)
    B(1)=Pair%B(1)
    B(2)=Pair%B(2)
    B(3)=Pair%B(3)
    AB2=Pair%AB2
    !
    ! Loop over contracted functions A and B.
    !
    IndexA=0
    DO CFA=1,BS%NCFnc%I(KA)
       StartLA = BS%LStrt%I(CFA,KA)
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       IndexB=0
       DO CFB=1,BS%NCFnc%I(KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          !
          ! Begin loops over primitives.
          !
          DO PFA=BS%NPFnc%I(CFA,KA),1,-1
             ZA=BS%Expnt%D(PFA,CFA,KA)
          DO PFB=BS%NPFnc%I(CFB,KB),1,-1
             ZB=BS%Expnt%D(PFB,CFB,KB)
             !
             Gamma=ZA+ZB
             InvGa=1.0D0/Gamma
             Dum=ZA*ZB*AB2*InvGa
#ifdef INT1E_THRESH
             IF(Dum>PrimPairDistanceThreshold) EXIT
#endif
             Ov=EXP(-Dum)*SQRT(PI*InvGa)*PI*InvGa
             !
             P(1)=(ZA*A(1)+ZB*B(1))*InvGa
             P(2)=(ZA*A(2)+ZB*B(2))*InvGa
             P(3)=(ZA*A(3)+ZB*B(3))*InvGa
             PA(1)=P(1)-A(1)
             PA(2)=P(2)-A(2)
             PA(3)=P(3)-A(3)
             PB(1)=P(1)-B(1)
             PB(2)=P(2)-B(2)
             PB(3)=P(3)-B(3)
             !
             ! Compute OS RR.
             !
             CALL I1E_Ovlap_OSRecurs(PA(1),PB(1),Gamma,MaxLA+2,MaxLB)
             !
             ! Contract each buffer into appropriate location.
             !
             IB=IndexB
             DO LMNB=StartLB,StopLB
                LB=BS%LxDex%I(LMNB)
                MB=BS%LyDex%I(LMNB)
                NB=BS%LzDex%I(LMNB)
                CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                IB=IB+1
             IA=IndexA
             DO LMNA=StartLA,StopLA
                LA=BS%LxDex%I(LMNA)
                MA=BS%LyDex%I(LMNA)
                NA=BS%LzDex%I(LMNA)
                CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                IA=IA+1
                !
                X0=RRX(LA,LB)
                Y0=RRY(MA,MB)
                Z0=RRZ(NA,NB)
                OvT=Ov*CA*CB
                TX=ZA*DBLE(2*LA+1)*RRX(LA,LB)-2.0D0*ZA*ZA*RRX(LA+2,LB)
                IF(LA.GT.1) TX=TX-0.5D0*DBLE(LA*(LA-1))*RRX(LA-2,LB)
                TY=ZA*DBLE(2*MA+1)*RRY(MA,MB)-2.0D0*ZA*ZA*RRY(MA+2,MB)
                IF(MA.GT.1) TY=TY-0.5D0*DBLE(MA*(MA-1))*RRY(MA-2,MB)
                TZ=ZA*DBLE(2*NA+1)*RRZ(NA,NB)-2.0D0*ZA*ZA*RRZ(NA+2,NB)
                IF(NA.GT.1) TZ=TZ-0.5D0*DBLE(NA*(NA-1))*RRZ(NA-2,NB)
                I1E(IA,IB)=I1E(IA,IB)+OvT*(TX*Y0*Z0+X0*TY*Z0+X0*Y0*TZ)
                !
             ENDDO
             ENDDO
             !
          ENDDO
#ifdef INT1E_THRESH
          IF(PFB.EQ.BS%NPFnc%I(CFB,KB)) EXIT
#endif
          ENDDO
          !
          IndexB=IndexB+StopLB-StartLB+1
       ENDDO
       IndexA=IndexA+StopLA-StartLA+1
    ENDDO
    !
  END SUBROUTINE I1E_Kinet_AtBlk
  !
  !
  SUBROUTINE I1E_EDipl_AtBlk(BS,Pair,IXYZ,MOrig,I1E)
!---------------------------------------------------------------------------------
!H Compute the electric dipole matrix:
!H $$
!H   (a|({\bf r-O})|b)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BSET)                               :: BS
    TYPE(AtomPair)                           :: Pair
    INTEGER                                  :: IXYZ
    REAL(DOUBLE), DIMENSION(3)               :: MOrig
    REAL(DOUBLE), DIMENSION(Pair%NA,Pair%NB) :: I1E
    !-------------------------------------------------------------------
    INTEGER :: CFA,CFB,IndexA,IndexB,StartLA,StartLB
    INTEGER :: StopLA,StopLB,MaxLA,MaxLB,LMNA,LMNB
    INTEGER :: KA,KB,PFA,PFB,IA,IB,LA,LB,MA,MB,NA,NB
    REAL(DOUBLE), DIMENSION(3) :: P,A,B,PA,PB,AO
    REAL(DOUBLE) :: AB2,ZA,ZB,CA,CB,Gamma,InvGa,Dum,Ov,OvT
    REAL(DOUBLE) :: X0,Y0,Z0,X1,Y1,Z1
    !-------------------------------------------------------------------
    !
    KA=Pair%KA
    KB=Pair%KB
    A(1)=Pair%A(1)
    A(2)=Pair%A(2)
    A(3)=Pair%A(3)
    B(1)=Pair%B(1)
    B(2)=Pair%B(2)
    B(3)=Pair%B(3)
    AB2=Pair%AB2
    AO(1)=A(1)-MOrig(1)
    AO(2)=A(2)-MOrig(2)
    AO(3)=A(3)-MOrig(3)
    !
    ! Loop over contracted functions A and B.
    !
    IndexA=0
    DO CFA=1,BS%NCFnc%I(KA)
       StartLA = BS%LStrt%I(CFA,KA)
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       IndexB=0
       DO CFB=1,BS%NCFnc%I(KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          !
          ! Begin loops over primitives.
          !
          DO PFA=BS%NPFnc%I(CFA,KA),1,-1
             ZA=BS%Expnt%D(PFA,CFA,KA)
          DO PFB=BS%NPFnc%I(CFB,KB),1,-1
             ZB=BS%Expnt%D(PFB,CFB,KB)
             !
             Gamma=ZA+ZB
             InvGa=1.0D0/Gamma
             Dum=ZA*ZB*AB2*InvGa
#ifdef INT1E_THRESH
             IF(Dum>PrimPairDistanceThreshold) EXIT
#endif
             Ov=EXP(-Dum)*SQRT(PI*InvGa)*PI*InvGa
             !
             P(1)=(ZA*A(1)+ZB*B(1))*InvGa
             P(2)=(ZA*A(2)+ZB*B(2))*InvGa
             P(3)=(ZA*A(3)+ZB*B(3))*InvGa
             PA(1)=P(1)-A(1)
             PA(2)=P(2)-A(2)
             PA(3)=P(3)-A(3)
             PB(1)=P(1)-B(1)
             PB(2)=P(2)-B(2)
             PB(3)=P(3)-B(3)
             !
             ! Compute OS RR.
             !
             CALL I1E_Ovlap_OSRecurs(PA(1),PB(1),Gamma,MaxLA+1,MaxLB)
             !
             ! Contract each buffer into appropriate location.
             !
             IB=IndexB
             DO LMNB=StartLB,StopLB
                LB=BS%LxDex%I(LMNB)
                MB=BS%LyDex%I(LMNB)
                NB=BS%LzDex%I(LMNB)
                CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                IB=IB+1
             IA=IndexA
             DO LMNA=StartLA,StopLA
                LA=BS%LxDex%I(LMNA)
                MA=BS%LyDex%I(LMNA)
                NA=BS%LzDex%I(LMNA)
                CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                IA=IA+1
                !
                X0=RRX(LA,LB)
                Y0=RRY(MA,MB)
                Z0=RRZ(NA,NB)
                X1=RRX(LA+1,LB)
                Y1=RRY(MA+1,MB)
                Z1=RRZ(NA+1,NB)
                OvT=Ov*CA*CB
                IF(IXYZ.EQ.1)I1E(IA,IB)=I1E(IA,IB)+OvT*(X1+AO(1)*X0)*Y0*Z0
                IF(IXYZ.EQ.2)I1E(IA,IB)=I1E(IA,IB)+OvT*(Y1+AO(2)*Y0)*Z0*X0
                IF(IXYZ.EQ.3)I1E(IA,IB)=I1E(IA,IB)+OvT*(Z1+AO(3)*Z0)*X0*Y0
                !
             ENDDO
             ENDDO
             !
          ENDDO
#ifdef INT1E_THRESH
          IF(PFB.EQ.BS%NPFnc%I(CFB,KB)) EXIT
#endif
          ENDDO
          !
          IndexB=IndexB+StopLB-StartLB+1
       ENDDO
       IndexA=IndexA+StopLA-StartLA+1
    ENDDO
    !
  END SUBROUTINE I1E_EDipl_AtBlk
  !
  !
  SUBROUTINE I1E_EQuad_AtBlk(BS,Pair,IXYZ,MOrig,I1E)
!---------------------------------------------------------------------------------
!H Compute the electric quadrupole matrix:
!H $$
!H   (a|({\bf r-O})({\bf r-O})^T|b)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BSET)                               :: BS
    TYPE(AtomPair)                           :: Pair
    INTEGER                                  :: IXYZ
    REAL(DOUBLE), DIMENSION(3)               :: MOrig
    REAL(DOUBLE), DIMENSION(Pair%NA,Pair%NB) :: I1E
    !-------------------------------------------------------------------
    INTEGER :: CFA,CFB,IndexA,IndexB,StartLA,StartLB
    INTEGER :: StopLA,StopLB,MaxLA,MaxLB,LMNA,LMNB
    INTEGER :: KA,KB,PFA,PFB,IA,IB,LA,LB,MA,MB,NA,NB
    REAL(DOUBLE), DIMENSION(3) :: P,A,B,PA,PB,AO
    REAL(DOUBLE) :: AB2,ZA,ZB,CA,CB,Gamma,InvGa,Dum,Ov,OvT
    REAL(DOUBLE) :: X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2
    !-------------------------------------------------------------------
    !
    KA=Pair%KA
    KB=Pair%KB
    A(1)=Pair%A(1)
    A(2)=Pair%A(2)
    A(3)=Pair%A(3)
    B(1)=Pair%B(1)
    B(2)=Pair%B(2)
    B(3)=Pair%B(3)
    AB2=Pair%AB2
    AO(1)=A(1)-MOrig(1)
    AO(2)=A(2)-MOrig(2)
    AO(3)=A(3)-MOrig(3)
    !
    ! Loop over contracted functions A and B.
    !
    IndexA=0
    DO CFA=1,BS%NCFnc%I(KA)
       StartLA = BS%LStrt%I(CFA,KA)
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       IndexB=0
       DO CFB=1,BS%NCFnc%I(KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          !
          ! Begin loops over primitives.
          !
          DO PFA=BS%NPFnc%I(CFA,KA),1,-1
             ZA=BS%Expnt%D(PFA,CFA,KA)
          DO PFB=BS%NPFnc%I(CFB,KB),1,-1
             ZB=BS%Expnt%D(PFB,CFB,KB)
             !
             Gamma=ZA+ZB
             InvGa=1.0D0/Gamma
             Dum=ZA*ZB*AB2*InvGa
#ifdef INT1E_THRESH
             IF(Dum>PrimPairDistanceThreshold) EXIT
#endif
             Ov=EXP(-Dum)*SQRT(PI*InvGa)*PI*InvGa
             !
             P(1)=(ZA*A(1)+ZB*B(1))*InvGa
             P(2)=(ZA*A(2)+ZB*B(2))*InvGa
             P(3)=(ZA*A(3)+ZB*B(3))*InvGa
             PA(1)=P(1)-A(1)
             PA(2)=P(2)-A(2)
             PA(3)=P(3)-A(3)
             PB(1)=P(1)-B(1)
             PB(2)=P(2)-B(2)
             PB(3)=P(3)-B(3)
             !
             ! Compute OS RR.
             !
             CALL I1E_Ovlap_OSRecurs(PA(1),PB(1),Gamma,MaxLA+2,MaxLB)
             !
             ! Contract each buffer into appropriate location.
             !
             IB=IndexB
             DO LMNB=StartLB,StopLB
                LB=BS%LxDex%I(LMNB)
                MB=BS%LyDex%I(LMNB)
                NB=BS%LzDex%I(LMNB)
                CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                IB=IB+1
             IA=IndexA
             DO LMNA=StartLA,StopLA
                LA=BS%LxDex%I(LMNA)
                MA=BS%LyDex%I(LMNA)
                NA=BS%LzDex%I(LMNA)
                CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                IA=IA+1
                !
                X0=RRX(LA,LB)
                Y0=RRY(MA,MB)
                Z0=RRZ(NA,NB)
                X1=RRX(LA+1,LB)
                Y1=RRY(MA+1,MB)
                Z1=RRZ(NA+1,NB)
                X2=RRX(LA+2,LB)
                Y2=RRY(MA+2,MB)
                Z2=RRZ(NA+2,NB)
                OvT=Ov*CA*CB
                IF(IXYZ.EQ.1)I1E(IA,IB)=I1E(IA,IB)+OvT* &
                       (X2+2.0D0*AO(1)*X1+AO(1)*AO(1)*X0)*Y0*Z0
                IF(IXYZ.EQ.2)I1E(IA,IB)=I1E(IA,IB)+OvT* &
                       (X1+AO(1)*X0)*(Y1+AO(2)*Y0)*Z0
                IF(IXYZ.EQ.3)I1E(IA,IB)=I1E(IA,IB)+OvT* &
                       (X1+AO(1)*X0)*(Z1+AO(3)*Z0)*Y0
                IF(IXYZ.EQ.4)I1E(IA,IB)=I1E(IA,IB)+OvT* &
                       (Y2+2.0D0*AO(2)*Y1+AO(2)*AO(2)*Y0)*Z0*X0
                IF(IXYZ.EQ.5)I1E(IA,IB)=I1E(IA,IB)+OvT* &
                       (Y1+AO(2)*Y0)*(Z1+AO(3)*Z0)*X0
                IF(IXYZ.EQ.6)I1E(IA,IB)=I1E(IA,IB)+OvT* &
                       (Z2+2.0D0*AO(3)*Z1+AO(3)*AO(3)*Z0)*X0*Y0
                !
             ENDDO
             ENDDO
             !
          ENDDO
#ifdef INT1E_THRESH
          IF(PFB.EQ.BS%NPFnc%I(CFB,KB)) EXIT
#endif
          ENDDO
          !
          IndexB=IndexB+StopLB-StartLB+1
       ENDDO
       IndexA=IndexA+StopLA-StartLA+1
    ENDDO
    !
  END SUBROUTINE I1E_EQuad_AtBlk
  !
  !
  SUBROUTINE I1E_Ovlap1C_AtBlk(BS,Pair,I1E)
!---------------------------------------------------------------------------------
!H Compute the first differentiated overlap matrix with respect
!H to geometric distortions:
!H $$
!H   \frac{\partial}{\partial A_i}(a|b)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BSET)                                 :: BS
    TYPE(AtomPair)                             :: Pair
    REAL(DOUBLE), DIMENSION(Pair%NA,Pair%NB,3) :: I1E
    !-------------------------------------------------------------------
    INTEGER :: CFA,CFB,IndexA,IndexB,StartLA,StartLB
    INTEGER :: StopLA,StopLB,MaxLA,MaxLB,LMNA,LMNB
    INTEGER :: KA,KB,PFA,PFB,IA,IB,LA,LB,MA,MB,NA,NB
    REAL(DOUBLE), DIMENSION(3) :: P,A,B,PA,PB
    REAL(DOUBLE) :: AB2,ZA,ZB,CA,CB,Gamma,InvGa,Dum,Ov,OvT
    REAL(DOUBLE) :: X0,Y0,Z0,AX,AY,AZ
    !-------------------------------------------------------------------
    !
    KA=Pair%KA
    KB=Pair%KB
    A(1)=Pair%A(1)
    A(2)=Pair%A(2)
    A(3)=Pair%A(3)
    B(1)=Pair%B(1)
    B(2)=Pair%B(2)
    B(3)=Pair%B(3)
    AB2=Pair%AB2
    !
    ! Loop over contracted functions A and B.
    !
    IndexA=0
    DO CFA=1,BS%NCFnc%I(KA)
       StartLA = BS%LStrt%I(CFA,KA)
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       IndexB=0
       DO CFB=1,BS%NCFnc%I(KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          !
          ! Begin loops over primitives.
          !
          DO PFA=BS%NPFnc%I(CFA,KA),1,-1
             ZA=BS%Expnt%D(PFA,CFA,KA)
          DO PFB=BS%NPFnc%I(CFB,KB),1,-1
             ZB=BS%Expnt%D(PFB,CFB,KB)
             !
             Gamma=ZA+ZB
             InvGa=1.0D0/Gamma
             Dum=ZA*ZB*AB2*InvGa
#ifdef INT1E_THRESH
             IF(Dum>PrimPairDistanceThreshold) EXIT
#endif
             Ov=EXP(-Dum)*SQRT(PI*InvGa)*PI*InvGa
             !
             P(1)=(ZA*A(1)+ZB*B(1))*InvGa
             P(2)=(ZA*A(2)+ZB*B(2))*InvGa
             P(3)=(ZA*A(3)+ZB*B(3))*InvGa
             PA(1)=P(1)-A(1)
             PA(2)=P(2)-A(2)
             PA(3)=P(3)-A(3)
             PB(1)=P(1)-B(1)
             PB(2)=P(2)-B(2)
             PB(3)=P(3)-B(3)
             !
             ! Compute OS RR.
             !
             CALL I1E_Ovlap_OSRecurs(PA(1),PB(1),Gamma,MaxLA+1,MaxLB)
             !
             ! Contract each buffer into appropriate location.
             !
             IB=IndexB
             DO LMNB=StartLB,StopLB
                LB=BS%LxDex%I(LMNB)
                MB=BS%LyDex%I(LMNB)
                NB=BS%LzDex%I(LMNB)
                CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                IB=IB+1
             IA=IndexA
             DO LMNA=StartLA,StopLA
                LA=BS%LxDex%I(LMNA)
                MA=BS%LyDex%I(LMNA)
                NA=BS%LzDex%I(LMNA)
                CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                IA=IA+1
                !
                X0=RRX(LA,LB)
                Y0=RRY(MA,MB)
                Z0=RRZ(NA,NB)
                OvT=Ov*CA*CB
                AX=2.0D0*ZA*RRX(LA+1,LB)
                AY=2.0D0*ZA*RRY(MA+1,MB)
                AZ=2.0D0*ZA*RRZ(NA+1,NB)
                IF(LA.GT.0)AX=AX-DBLE(LA)*RRX(LA-1,LB)
                IF(MA.GT.0)AY=AY-DBLE(MA)*RRY(MA-1,MB)
                IF(NA.GT.0)AZ=AZ-DBLE(NA)*RRZ(NA-1,NB)
                I1E(IA,IB,1)=I1E(IA,IB,1)+OvT*AX*Y0*Z0
                I1E(IA,IB,2)=I1E(IA,IB,2)+OvT*X0*AY*Z0
                I1E(IA,IB,3)=I1E(IA,IB,3)+OvT*X0*Y0*AZ
                !
             ENDDO
             ENDDO
             !
          ENDDO
#ifdef INT1E_THRESH
          IF(PFB.EQ.BS%NPFnc%I(CFB,KB)) EXIT
#endif
          ENDDO
          !
          IndexB=IndexB+StopLB-StartLB+1
       ENDDO
       IndexA=IndexA+StopLA-StartLA+1
    ENDDO
    !
  END SUBROUTINE I1E_Ovlap1C_AtBlk
  !
  !
  SUBROUTINE I1E_Ovlap2C_AtBlk(BS,Pair,I1E)
!---------------------------------------------------------------------------------
!H Compute the second differentiated overlap matrix with respect
!H to geometric distortions:
!H $$
!H   \frac{\partial^2}{\partial A_i\partial A_j}(a|b)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BSET)                                 :: BS
    TYPE(AtomPair)                             :: Pair
    REAL(DOUBLE), DIMENSION(Pair%NA,Pair%NB,6) :: I1E
    !-------------------------------------------------------------------
    INTEGER :: CFA,CFB,IndexA,IndexB,StartLA,StartLB
    INTEGER :: StopLA,StopLB,MaxLA,MaxLB,LMNA,LMNB
    INTEGER :: KA,KB,PFA,PFB,IA,IB,LA,LB,MA,MB,NA,NB
    REAL(DOUBLE), DIMENSION(3) :: P,A,B,PA,PB
    REAL(DOUBLE) :: AB2,ZA,ZB,CA,CB,Gamma,InvGa,Dum,Ov,OvT,ZA2
    REAL(DOUBLE) :: X0,Y0,Z0,AX,AY,AZ,AX2,AY2,AZ2
    !-------------------------------------------------------------------
    !
    KA=Pair%KA
    KB=Pair%KB
    A(1)=Pair%A(1)
    A(2)=Pair%A(2)
    A(3)=Pair%A(3)
    B(1)=Pair%B(1)
    B(2)=Pair%B(2)
    B(3)=Pair%B(3)
    AB2=Pair%AB2
    !
    ! Loop over contracted functions A and B.
    !
    IndexA=0
    DO CFA=1,BS%NCFnc%I(KA)
       StartLA = BS%LStrt%I(CFA,KA)
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       IndexB=0
       DO CFB=1,BS%NCFnc%I(KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          !
          ! Begin loops over primitives.
          !
          DO PFA=BS%NPFnc%I(CFA,KA),1,-1
             ZA=BS%Expnt%D(PFA,CFA,KA)
          DO PFB=BS%NPFnc%I(CFB,KB),1,-1
             ZB=BS%Expnt%D(PFB,CFB,KB)
             !
             Gamma=ZA+ZB
             InvGa=1.0D0/Gamma
             Dum=ZA*ZB*AB2*InvGa
#ifdef INT1E_THRESH
             IF(Dum>PrimPairDistanceThreshold) EXIT
#endif
             Ov=EXP(-Dum)*SQRT(PI*InvGa)*PI*InvGa
             !
             P(1)=(ZA*A(1)+ZB*B(1))*InvGa
             P(2)=(ZA*A(2)+ZB*B(2))*InvGa
             P(3)=(ZA*A(3)+ZB*B(3))*InvGa
             PA(1)=P(1)-A(1)
             PA(2)=P(2)-A(2)
             PA(3)=P(3)-A(3)
             PB(1)=P(1)-B(1)
             PB(2)=P(2)-B(2)
             PB(3)=P(3)-B(3)
             !
             ! Compute OS RR.
             !
             CALL I1E_Ovlap_OSRecurs(PA(1),PB(1),Gamma,MaxLA+1,MaxLB)
             !
             ! Contract each buffer into appropriate location.
             !
             IB=IndexB
             DO LMNB=StartLB,StopLB
                LB=BS%LxDex%I(LMNB)
                MB=BS%LyDex%I(LMNB)
                NB=BS%LzDex%I(LMNB)
                CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                IB=IB+1
             IA=IndexA
             DO LMNA=StartLA,StopLA
                LA=BS%LxDex%I(LMNA)
                MA=BS%LyDex%I(LMNA)
                NA=BS%LzDex%I(LMNA)
                CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                IA=IA+1
                !
                X0=RRX(LA,LB)
                Y0=RRY(MA,MB)
                Z0=RRZ(NA,NB)
                OvT=Ov*CA*CB
                ZA2=ZA*ZA
                AX=2.0D0*ZA*RRX(LA+1,LB)
                AY=2.0D0*ZA*RRY(MA+1,MB)
                AZ=2.0D0*ZA*RRZ(NA+1,NB)
                IF(LA.GT.0)AX=AX-DBLE(LA)*RRX(LA-1,LB)
                IF(MA.GT.0)AY=AY-DBLE(MA)*RRY(MA-1,MB)
                IF(NA.GT.0)AZ=AZ-DBLE(NA)*RRZ(NA-1,NB)
                AX2=4D0*ZA2*RRX(LA+2,LB)-2D0*ZA*DBLE(2*LA+1)*X0
                AY2=4D0*ZA2*RRY(MA+2,MB)-2D0*ZA*DBLE(2*MA+1)*Y0
                AZ2=4D0*ZA2*RRZ(NA+2,NB)-2D0*ZA*DBLE(2*NA+1)*Z0
                IF(LA.GT.1)AX2=AX2+DBLE(LA*(LA-1))*RRX(LA-2,LB)
                IF(MA.GT.1)AY2=AY2+DBLE(MA*(MA-1))*RRY(MA-2,MB)
                IF(NA.GT.1)AZ2=AZ2+DBLE(NA*(NA-1))*RRZ(NA-2,NB)
                I1E(IA,IB,1)=I1E(IA,IB,1)+OvT*AX2*Y0 *Z0  !xx
                I1E(IA,IB,2)=I1E(IA,IB,2)+OvT*AX *AY *Z0  !xy
                I1E(IA,IB,3)=I1E(IA,IB,3)+OvT*AX *Y0 *AZ  !xz
                I1E(IA,IB,4)=I1E(IA,IB,4)+OvT*X0 *AY2*Z0  !yy
                I1E(IA,IB,5)=I1E(IA,IB,5)+OvT*X0 *AY *AZ  !yz
                I1E(IA,IB,6)=I1E(IA,IB,6)+OvT*X0 *Y0 *AZ2 !zz
                !
             ENDDO
             ENDDO
             !
          ENDDO
#ifdef INT1E_THRESH
          IF(PFB.EQ.BS%NPFnc%I(CFB,KB)) EXIT
#endif
          ENDDO
          !
          IndexB=IndexB+StopLB-StartLB+1
       ENDDO
       IndexA=IndexA+StopLA-StartLA+1
    ENDDO
    !
  END SUBROUTINE I1E_Ovlap2C_AtBlk
  !
  !
  SUBROUTINE I1E_Kinet1C_AtBlk(BS,Pair,I1E)
!---------------------------------------------------------------------------------
!H Compute the first differentiated kinetic matrix with respect
!H to geometric distortions:
!H $$
!H   -\frac{1}{2}\frac{\partial}{\partial A_i}(a|\nabla^2|b)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BSET)                                 :: BS
    TYPE(AtomPair)                             :: Pair
    REAL(DOUBLE), DIMENSION(Pair%NA,Pair%NB,3) :: I1E
    !-------------------------------------------------------------------
    INTEGER :: CFA,CFB,IndexA,IndexB,StartLA,StartLB
    INTEGER :: StopLA,StopLB,MaxLA,MaxLB,LMNA,LMNB
    INTEGER :: KA,KB,PFA,PFB,IA,IB,LA,LB,MA,MB,NA,NB
    REAL(DOUBLE), DIMENSION(3) :: P,A,B,PA,PB
    REAL(DOUBLE) :: AB2,ZA,ZB,CA,CB,Gamma,InvGa,Dum,Ov,OvT,ZA2,ZA3
    REAL(DOUBLE) :: X0,Y0,Z0,AX,AY,AZ,TX,TY,TZ,TAX,TAY,TAZ
    !-------------------------------------------------------------------
    !
    KA=Pair%KA
    KB=Pair%KB
    A(1)=Pair%A(1)
    A(2)=Pair%A(2)
    A(3)=Pair%A(3)
    B(1)=Pair%B(1)
    B(2)=Pair%B(2)
    B(3)=Pair%B(3)
    AB2=Pair%AB2
    !
    ! Loop over contracted functions A and B.
    !
    IndexA=0
    DO CFA=1,BS%NCFnc%I(KA)
       StartLA = BS%LStrt%I(CFA,KA)
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       IndexB=0
       DO CFB=1,BS%NCFnc%I(KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          !
          ! Begin loops over primitives.
          !
          DO PFA=BS%NPFnc%I(CFA,KA),1,-1
             ZA=BS%Expnt%D(PFA,CFA,KA)
          DO PFB=BS%NPFnc%I(CFB,KB),1,-1
             ZB=BS%Expnt%D(PFB,CFB,KB)
             !
             Gamma=ZA+ZB
             InvGa=1.0D0/Gamma
             Dum=ZA*ZB*AB2*InvGa
#ifdef INT1E_THRESH
             IF(Dum>PrimPairDistanceThreshold) EXIT
#endif
             Ov=EXP(-Dum)*SQRT(PI*InvGa)*PI*InvGa
             !
             P(1)=(ZA*A(1)+ZB*B(1))*InvGa
             P(2)=(ZA*A(2)+ZB*B(2))*InvGa
             P(3)=(ZA*A(3)+ZB*B(3))*InvGa
             PA(1)=P(1)-A(1)
             PA(2)=P(2)-A(2)
             PA(3)=P(3)-A(3)
             PB(1)=P(1)-B(1)
             PB(2)=P(2)-B(2)
             PB(3)=P(3)-B(3)
             !
             ! Compute OS RR.
             !
             CALL I1E_Ovlap_OSRecurs(PA(1),PB(1),Gamma,MaxLA+3,MaxLB)
             !
             ! Contract each buffer into appropriate location.
             !
             IB=IndexB
             DO LMNB=StartLB,StopLB
                LB=BS%LxDex%I(LMNB)
                MB=BS%LyDex%I(LMNB)
                NB=BS%LzDex%I(LMNB)
                CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                IB=IB+1
             IA=IndexA
             DO LMNA=StartLA,StopLA
                LA=BS%LxDex%I(LMNA)
                MA=BS%LyDex%I(LMNA)
                NA=BS%LzDex%I(LMNA)
                CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                IA=IA+1
                !
                X0=RRX(LA,LB)
                Y0=RRY(MA,MB)
                Z0=RRZ(NA,NB)
                OvT=Ov*CA*CB
                ZA2=ZA*ZA
                ZA3=ZA2*ZA
                AX =2D0*ZA*RRX(LA+1,LB)
                TX =ZA*DBLE(2*LA+1)*RRX(LA,LB)-2D0*ZA2*RRX(LA+2,LB)
                TAX=6D0*ZA2*DBLE(LA+1)*RRX(LA+1,LB)-4D0*ZA3*RRX(LA+3,LB)
                IF(LA.GT.0) THEN
                   AX =AX-DBLE(LA)*RRX(LA-1,LB)
                   TAX=TAX-3D0*ZA*DBLE(LA*LA)*RRX(LA-1,LB)
                   IF(LA.GT.1) THEN
                      TX=TX -0.5D0*DBLE(LA*(LA-1))*RRX(LA-2,LB)
                      IF(LA.GT.2)TAX=TAX+0.5D0*DBLE(LA*(LA-1)*(LA-2))*RRX(LA-3,LB)
                   ENDIF
                ENDIF
                AY =2D0*ZA*RRY(MA+1,MB)
                TY =ZA*DBLE(2*MA+1)*RRY(MA,MB)-2D0*ZA2*RRY(MA+2,MB)
                TAY=6D0*ZA2*DBLE(MA+1)*RRY(MA+1,MB)-4D0*ZA3*RRY(MA+3,MB)
                IF(MA.GT.0) THEN
                   AY =AY-DBLE(MA)*RRY(MA-1,MB)
                   TAY=TAY-3D0*ZA*DBLE(MA*MA)*RRY(MA-1,MB)
                   IF(MA.GT.1) THEN
                      TY =TY -0.5D0*DBLE(MA*(MA-1))*RRY(MA-2,MB)
                      IF(MA.GT.2)TAY=TAY+0.5D0*DBLE(MA*(MA-1)*(MA-2))*RRY(MA-3,MB)
                   ENDIF
                ENDIF
                AZ =2D0*ZA*RRZ(NA+1,NB)
                TZ =ZA*DBLE(2*NA+1)*RRZ(NA,NB)-2D0*ZA2*RRZ(NA+2,NB)
                TAZ=6D0*ZA2*DBLE(NA+1)*RRZ(NA+1,NB)-4D0*ZA3*RRZ(NA+3,NB)
                IF(NA.GT.0) THEN
                   AZ =AZ-DBLE(NA)*RRZ(NA-1,NB)
                   TAZ=TAZ-3D0*ZA*DBLE(NA*NA)*RRZ(NA-1,NB)
                   IF(NA.GT.1) THEN
                      TZ =TZ -0.5D0*DBLE(NA*(NA-1))*RRZ(NA-2,NB)
                      IF(NA.GT.2)TAZ=TAZ+0.5D0*DBLE(NA*(NA-1)*(NA-2))*RRZ(NA-3,NB)
                   ENDIF
                ENDIF
                I1E(IA,IB,1)=I1E(IA,IB,1)+OvT*(TAX*Y0*Z0+AX*TY *Z0+AX*Y0*TZ )
                I1E(IA,IB,2)=I1E(IA,IB,2)+OvT*(TX *AY*Z0+X0*TAY*Z0+X0*AY*TZ )
                I1E(IA,IB,3)=I1E(IA,IB,3)+OvT*(TX *Y0*AZ+X0*TY *AZ+X0*Y0*TAZ)
                !
             ENDDO
             ENDDO
             !
          ENDDO
#ifdef INT1E_THRESH
          IF(PFB.EQ.BS%NPFnc%I(CFB,KB)) EXIT
#endif
          ENDDO
          !
          IndexB=IndexB+StopLB-StartLB+1
       ENDDO
       IndexA=IndexA+StopLA-StartLA+1
    ENDDO
    !
  END SUBROUTINE I1E_Kinet1C_AtBlk
  !
  !
  SUBROUTINE I1E_Kinet2C_AtBlk(BS,Pair,I1E)
!---------------------------------------------------------------------------------
!H Compute the second differentiated kinetic matrix with respect
!H to geometric distortions:
!H $$
!H   -\frac{1}{2}\frac{\partial^2}{\partial A_i\partial A_j}(a|\nabla^2|b)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BSET)                                 :: BS
    TYPE(AtomPair)                             :: Pair
    REAL(DOUBLE), DIMENSION(Pair%NA,Pair%NB,6) :: I1E
    !-------------------------------------------------------------------
    INTEGER :: CFA,CFB,IndexA,IndexB,StartLA,StartLB
    INTEGER :: StopLA,StopLB,MaxLA,MaxLB,LMNA,LMNB
    INTEGER :: KA,KB,PFA,PFB,IA,IB,LA,LB,MA,MB,NA,NB
    REAL(DOUBLE), DIMENSION(3) :: P,A,B,PA,PB
    REAL(DOUBLE) :: AB2,ZA,ZB,CA,CB,Gamma,InvGa,Dum,Ov,OvT,ZA2,ZA3,ZA4
    REAL(DOUBLE) :: X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,AX,AY,AZ,AX2,AY2,AZ2
    REAL(DOUBLE) :: TX,TY,TZ,TAX,TAY,TAZ,TAX2,TAY2,TAZ2
    !-------------------------------------------------------------------
    !
    KA=Pair%KA
    KB=Pair%KB
    A(1)=Pair%A(1)
    A(2)=Pair%A(2)
    A(3)=Pair%A(3)
    B(1)=Pair%B(1)
    B(2)=Pair%B(2)
    B(3)=Pair%B(3)
    AB2=Pair%AB2
    !
    ! Loop over contracted functions A and B.
    !
    IndexA=0
    DO CFA=1,BS%NCFnc%I(KA)
       StartLA = BS%LStrt%I(CFA,KA)
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       IndexB=0
       DO CFB=1,BS%NCFnc%I(KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          !
          ! Begin loops over primitives.
          !
          DO PFA=BS%NPFnc%I(CFA,KA),1,-1
             ZA=BS%Expnt%D(PFA,CFA,KA)
             ZA2=ZA *ZA
             ZA3=ZA2*ZA
             ZA4=ZA3*ZA
          DO PFB=BS%NPFnc%I(CFB,KB),1,-1
             ZB=BS%Expnt%D(PFB,CFB,KB)
             !
             Gamma=ZA+ZB
             InvGa=1.0D0/Gamma
             Dum=ZA*ZB*AB2*InvGa
#ifdef INT1E_THRESH
             IF(Dum>PrimPairDistanceThreshold) EXIT
#endif
             Ov=EXP(-Dum)*SQRT(PI*InvGa)*PI*InvGa
             !
             P(1)=(ZA*A(1)+ZB*B(1))*InvGa
             P(2)=(ZA*A(2)+ZB*B(2))*InvGa
             P(3)=(ZA*A(3)+ZB*B(3))*InvGa
             PA(1)=P(1)-A(1)
             PA(2)=P(2)-A(2)
             PA(3)=P(3)-A(3)
             PB(1)=P(1)-B(1)
             PB(2)=P(2)-B(2)
             PB(3)=P(3)-B(3)
             !
             ! Compute OS RR.
             !
             CALL I1E_Ovlap_OSRecurs(PA(1),PB(1),Gamma,MaxLA+4,MaxLB)
             !
             ! Contract each buffer into appropriate location.
             !
             IB=IndexB
             DO LMNB=StartLB,StopLB
                LB=BS%LxDex%I(LMNB)
                MB=BS%LyDex%I(LMNB)
                NB=BS%LzDex%I(LMNB)
                CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                IB=IB+1
             IA=IndexA
             DO LMNA=StartLA,StopLA
                LA=BS%LxDex%I(LMNA)
                MA=BS%LyDex%I(LMNA)
                NA=BS%LzDex%I(LMNA)
                CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                IA=IA+1
                !
                X0=RRX(LA,LB)
                Y0=RRY(MA,MB)
                Z0=RRZ(NA,NB)
                X1=RRX(LA+1,LB)
                Y1=RRY(MA+1,MB)
                Z1=RRZ(NA+1,NB)
                X2=RRX(LA+2,LB)
                Y2=RRY(MA+2,MB)
                Z2=RRZ(NA+2,NB)
                OvT=Ov*CA*CB
                AX  =2D0*ZA*X1
                AX2 =4D0*ZA2*X2-2D0*ZA*DBLE(2*LA+1)*X0
                TX  =ZA*DBLE(2*LA+1)*X0-2D0*ZA2*X2
                TAX =6D0*ZA2*DBLE(LA+1)*X1-4D0*ZA3*RRX(LA+3,LB)
                TAX2=-8D0*ZA4*RRX(LA+4,LB)+8D0*ZA3*DBLE(2*LA+3)*X2 &
                       -6D0*ZA2*DBLE(2*LA*(LA+1)+1)*X0
                IF(LA.GT.0) THEN
                   AX =AX-DBLE(LA)*RRX(LA-1,LB)
                   TAX=TAX-3D0*ZA*DBLE(LA*LA)*RRX(LA-1,LB)
                   IF(LA.GT.1) THEN
                      AX2 =AX2+DBLE(LA*(LA-1))*RRX(LA-2,LB)
                      TX  =TX -0.5D0*DBLE(LA*(LA-1))*RRX(LA-2,LB)
                      TAX2=TAX2+2D0*ZA*DBLE(LA*(LA-1)*(2*LA-1))*RRX(LA-2,LB)
                      IF(LA.GT.2) THEN
                         TAX=TAX+0.5D0*DBLE(LA*(LA-1)*(LA-2))*RRX(LA-3,LB)
                         IF(LA.GT.3)TAX2=TAX2-0.5D0*DBLE(LA*(LA-1)*(LA-2)*(LA-3))*RRX(LA-4,LB)
                      ENDIF
                   ENDIF
                ENDIF
                AY  =2D0*ZA*Y1
                AY2 =4D0*ZA2*Y2-2D0*ZA*DBLE(2*MA+1)*Y0
                TY  =ZA*DBLE(2*MA+1)*Y0-2D0*ZA2*Y2
                TAY =6D0*ZA2*DBLE(MA+1)*Y1-4D0*ZA3*RRY(MA+3,MB)
                TAY2=-8D0*ZA4*RRY(MA+4,MB)+8D0*ZA3*DBLE(2*MA+3)*Y2 &
                       -6D0*ZA2*DBLE(2*MA*(MA+1)+1)*Y0
                IF(MA.GT.0) THEN
                   AY =AY-DBLE(MA)*RRY(MA-1,MB)
                   TAY=TAY-3D0*ZA*DBLE(MA*MA)*RRY(MA-1,MB)
                   IF(MA.GT.1) THEN
                      AY2 =AY2+DBLE(MA*(MA-1))*RRY(MA-2,MB)
                      TY  =TY -0.5D0*DBLE(MA*(MA-1))*RRY(MA-2,MB)
                      TAY2=TAY2+2D0*ZA*DBLE(MA*(MA-1)*(2*MA-1))*RRY(MA-2,MB)
                      IF(MA.GT.2) THEN
                         TAY=TAY+0.5D0*DBLE(MA*(MA-1)*(MA-2))*RRY(MA-3,MB)
                         IF(MA.GT.3)TAY2=TAY2-0.5D0*DBLE(MA*(MA-1)*(MA-2)*(MA-3))*RRY(MA-4,MB)
                      ENDIF
                   ENDIF
                ENDIF
                AZ  =2D0*ZA*Z1
                AZ2 =4D0*ZA2*Z2-2D0*ZA*DBLE(2*NA+1)*Z0
                TZ  =ZA*DBLE(2*NA+1)*Z0-2D0*ZA2*Z2
                TAZ =6D0*ZA2*DBLE(NA+1)*Z1-4D0*ZA3*RRZ(NA+3,NB)
                TAZ2=-8D0*ZA4*RRZ(NA+4,NB)+8D0*ZA3*DBLE(2*NA+3)*Z2 &
                       -6D0*ZA2*DBLE(2*NA*(NA+1)+1)*Z0
                IF(NA.GT.0) THEN
                   AZ =AZ-DBLE(NA)*RRZ(NA-1,NB)
                   TAZ=TAZ-3D0*ZA*DBLE(NA*NA)*RRZ(NA-1,NB)
                   IF(NA.GT.1) THEN
                      AZ2 =AZ2+DBLE(NA*(NA-1))*RRZ(NA-2,NB)
                      TZ  =TZ -0.5D0*DBLE(NA*(NA-1))*RRZ(NA-2,NB)
                      TAZ2=TAZ2+2D0*ZA*DBLE(NA*(NA-1)*(2*NA-1))*RRZ(NA-2,NB)
                      IF(NA.GT.2) THEN
                         TAZ=TAZ+0.5D0*DBLE(NA*(NA-1)*(NA-2))*RRZ(NA-3,NB)
                         IF(NA.GT.3)TAZ2=TAZ2-0.5D0*DBLE(NA*(NA-1)*(NA-2)*(NA-3))*RRZ(NA-4,NB)
                      ENDIF
                   ENDIF
                ENDIF
                I1E(IA,IB,1)=I1E(IA,IB,1)+OvT*(TAX2*Y0 *Z0 + AX2*TY  *Z0  + AX2*Y0 *TZ  ) !xx
                I1E(IA,IB,2)=I1E(IA,IB,2)+OvT*(TAX *AY *Z0 + AX *TAY *Z0  + AX *AY *TZ  ) !xy
                I1E(IA,IB,3)=I1E(IA,IB,3)+OvT*(TAX *Y0 *AZ + AX *TY  *AZ  + AX *Y0 *TAZ ) !xz
                I1E(IA,IB,4)=I1E(IA,IB,4)+OvT*(TX  *AY2*Z0 + X0 *TAY2*Z0  + X0 *AY2*TZ  ) !yy
                I1E(IA,IB,5)=I1E(IA,IB,5)+OvT*(TX  *AY *AZ + X0 *TAY *AZ  + X0 *AY *TAZ ) !yz
                I1E(IA,IB,6)=I1E(IA,IB,6)+OvT*(TX  *Y0*AZ2 + X0 *TY  *AZ2 + X0 *Y0 *TAZ2) !zz
                !
             ENDDO
             ENDDO
             !
          ENDDO
#ifdef INT1E_THRESH
          IF(PFB.EQ.BS%NPFnc%I(CFB,KB)) EXIT
#endif
          ENDDO
          !
          IndexB=IndexB+StopLB-StartLB+1
       ENDDO
       IndexA=IndexA+StopLA-StartLA+1
    ENDDO
    !
  END SUBROUTINE I1E_Kinet2C_AtBlk
  !
  !
  SUBROUTINE I1E_EDipl1C_AtBlk(BS,Pair,MOrig,I1E)
!---------------------------------------------------------------------------------
!H Compute the electric dipole matrix:
!H $$
!H   \frac{\partial}{\partial A_i}(a|({\bf r-O})|b)\quad
!H   \frac{\partial}{\partial O}(a|({\bf r-O})|b)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BSET)                                  :: BS
    TYPE(AtomPair)                              :: Pair
    REAL(DOUBLE), DIMENSION(3)                  :: MOrig
    REAL(DOUBLE), DIMENSION(Pair%NA,Pair%NB,10) :: I1E
    !-------------------------------------------------------------------
    INTEGER :: CFA,CFB,IndexA,IndexB,StartLA,StartLB
    INTEGER :: StopLA,StopLB,MaxLA,MaxLB,LMNA,LMNB
    INTEGER :: KA,KB,PFA,PFB,IA,IB,LA,LB,MA,MB,NA,NB
    REAL(DOUBLE), DIMENSION(3) :: P,A,B,PA,PB,AO
    REAL(DOUBLE) :: AB2,ZA,ZB,CA,CB,Gamma,InvGa,Dum,Ov,OvT
    REAL(DOUBLE) :: X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,AX0,AY0,AZ0,AX1,AY1,AZ1
    !-------------------------------------------------------------------
    !
    KA=Pair%KA
    KB=Pair%KB
    A(1)=Pair%A(1)
    A(2)=Pair%A(2)
    A(3)=Pair%A(3)
    B(1)=Pair%B(1)
    B(2)=Pair%B(2)
    B(3)=Pair%B(3)
    AB2=Pair%AB2
    AO(1)=A(1)-MOrig(1)
    AO(2)=A(2)-MOrig(2)
    AO(3)=A(3)-MOrig(3)
    !
    ! Loop over contracted functions A and B.
    !
    IndexA=0
    DO CFA=1,BS%NCFnc%I(KA)
       StartLA = BS%LStrt%I(CFA,KA)
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       IndexB=0
       DO CFB=1,BS%NCFnc%I(KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          !
          ! Begin loops over primitives.
          !
          DO PFA=BS%NPFnc%I(CFA,KA),1,-1
             ZA=BS%Expnt%D(PFA,CFA,KA)
          DO PFB=BS%NPFnc%I(CFB,KB),1,-1
             ZB=BS%Expnt%D(PFB,CFB,KB)
             !
             Gamma=ZA+ZB
             InvGa=1.0D0/Gamma
             Dum=ZA*ZB*AB2*InvGa
#ifdef INT1E_THRESH
             IF(Dum>PrimPairDistanceThreshold) EXIT
#endif
             Ov=EXP(-Dum)*SQRT(PI*InvGa)*PI*InvGa
             !
             P(1)=(ZA*A(1)+ZB*B(1))*InvGa
             P(2)=(ZA*A(2)+ZB*B(2))*InvGa
             P(3)=(ZA*A(3)+ZB*B(3))*InvGa
             PA(1)=P(1)-A(1)
             PA(2)=P(2)-A(2)
             PA(3)=P(3)-A(3)
             PB(1)=P(1)-B(1)
             PB(2)=P(2)-B(2)
             PB(3)=P(3)-B(3)
             !
             ! Compute OS RR.
             !
             CALL I1E_Ovlap_OSRecurs(PA(1),PB(1),Gamma,MaxLA+2,MaxLB)
             !
             ! Contract each buffer into appropriate location.
             !
             IB=IndexB
             DO LMNB=StartLB,StopLB
                LB=BS%LxDex%I(LMNB)
                MB=BS%LyDex%I(LMNB)
                NB=BS%LzDex%I(LMNB)
                CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                IB=IB+1
             IA=IndexA
             DO LMNA=StartLA,StopLA
                LA=BS%LxDex%I(LMNA)
                MA=BS%LyDex%I(LMNA)
                NA=BS%LzDex%I(LMNA)
                CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                IA=IA+1
                !
                X0=RRX(LA,LB)
                Y0=RRY(MA,MB)
                Z0=RRZ(NA,NB)
                X1=RRX(LA+1,LB)
                Y1=RRY(MA+1,MB)
                Z1=RRZ(NA+1,NB)
                X2=RRX(LA+2,LB)
                Y2=RRY(MA+2,MB)
                Z2=RRZ(NA+2,NB)
                OvT=Ov*CA*CB
                AX0=2.0D0*ZA*X1
                AY0=2.0D0*ZA*Y1
                AZ0=2.0D0*ZA*Z1
                IF(LA.GT.0)AX0=AX0-DBLE(LA)*RRX(LA-1,LB)
                IF(MA.GT.0)AY0=AY0-DBLE(MA)*RRY(MA-1,MB)
                IF(NA.GT.0)AZ0=AZ0-DBLE(NA)*RRZ(NA-1,NB)
                AX1=2.0D0*ZA*X2-DBLE(LA)*X0
                AY1=2.0D0*ZA*Y2-DBLE(MA)*Y0
                AZ1=2.0D0*ZA*Z2-DBLE(NA)*Z0
                !x
                I1E(IA,IB, 1)=I1E(IA,IB, 1)+OvT*(AX1+AO(1)*AX0)*Y0 *Z0
                I1E(IA,IB, 2)=I1E(IA,IB, 2)+OvT*(Y1 +AO(2)*Y0 )*Z0 *AX0
                I1E(IA,IB, 3)=I1E(IA,IB, 3)+OvT*(Z1 +AO(3)*Z0 )*AX0*Y0
                !y
                I1E(IA,IB, 4)=I1E(IA,IB, 4)+OvT*(X1 +AO(1)*X0 )*AY0*Z0
                I1E(IA,IB, 5)=I1E(IA,IB, 5)+OvT*(AY1+AO(2)*AY0)*Z0 *X0
                I1E(IA,IB, 6)=I1E(IA,IB, 6)+OvT*(Z1 +AO(3)*Z0 )*X0 *AY0
                !z
                I1E(IA,IB, 7)=I1E(IA,IB, 7)+OvT*(X1 +AO(1)*X0 )*Y0 *AZ0
                I1E(IA,IB, 8)=I1E(IA,IB, 8)+OvT*(Y1 +AO(2)*Y0 )*AZ0*X0
                I1E(IA,IB, 9)=I1E(IA,IB, 9)+OvT*(AZ1+AO(3)*AZ0)*X0 *Y0
                !O
                I1E(IA,IB,10)=I1E(IA,IB,10)-OvT*X0*Y0*Z0
                !
             ENDDO
             ENDDO
             !
          ENDDO
#ifdef INT1E_THRESH
          IF(PFB.EQ.BS%NPFnc%I(CFB,KB)) EXIT
#endif
          ENDDO
          !
          IndexB=IndexB+StopLB-StartLB+1
       ENDDO
       IndexA=IndexA+StopLA-StartLA+1
    ENDDO
    !
  END SUBROUTINE I1E_EDipl1C_AtBlk
  !
  !
  SUBROUTINE I1E_KinMom_AtBlk(BS,Pair,IXYZ,MOrig,I1E)
!---------------------------------------------------------------------------------
!H Compute the first differentiated kinetic moment matrix with respect
!H to geometric distortions:
!H $$
!H   (a|({\bf r-O})\times\nabla|b)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BSET)                               :: BS
    TYPE(AtomPair)                           :: Pair
    INTEGER                                  :: IXYZ
    REAL(DOUBLE), DIMENSION(3)               :: MOrig
    REAL(DOUBLE), DIMENSION(Pair%NA,Pair%NB) :: I1E
    !-------------------------------------------------------------------
    INTEGER :: CFA,CFB,IndexA,IndexB,StartLA,StartLB
    INTEGER :: StopLA,StopLB,MaxLA,MaxLB,LMNA,LMNB
    INTEGER :: KA,KB,PFA,PFB,IA,IB,LA,LB,MA,MB,NA,NB
    REAL(DOUBLE), DIMENSION(3) :: P,A,B,PA,PB,AO
    REAL(DOUBLE) :: AB2,ZA,ZB,CA,CB,Gamma,InvGa,Dum,Ov,OvT
    REAL(DOUBLE) :: X0,Y0,Z0,X1,Y1,Z1,DX,DY,DZ
    !-------------------------------------------------------------------
    !
    KA=Pair%KA
    KB=Pair%KB
    A(1)=Pair%A(1)
    A(2)=Pair%A(2)
    A(3)=Pair%A(3)
    B(1)=Pair%B(1)
    B(2)=Pair%B(2)
    B(3)=Pair%B(3)
    AB2=Pair%AB2
    AO(1)=A(1)-MOrig(1)
    AO(2)=A(2)-MOrig(2)
    AO(3)=A(3)-MOrig(3)
    !
    ! Loop over contracted functions A and B.
    !
    IndexA=0
    DO CFA=1,BS%NCFnc%I(KA)
       StartLA = BS%LStrt%I(CFA,KA)
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       IndexB=0
       DO CFB=1,BS%NCFnc%I(KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          !
          ! Begin loops over primitives.
          !
          DO PFA=BS%NPFnc%I(CFA,KA),1,-1
             ZA=BS%Expnt%D(PFA,CFA,KA)
          DO PFB=BS%NPFnc%I(CFB,KB),1,-1
             ZB=BS%Expnt%D(PFB,CFB,KB)
             !
             Gamma=ZA+ZB
             InvGa=1.0D0/Gamma
             Dum=ZA*ZB*AB2*InvGa
#ifdef INT1E_THRESH
             IF(Dum>PrimPairDistanceThreshold) EXIT
#endif
             Ov=EXP(-Dum)*SQRT(PI*InvGa)*PI*InvGa
             !
             P(1)=(ZA*A(1)+ZB*B(1))*InvGa
             P(2)=(ZA*A(2)+ZB*B(2))*InvGa
             P(3)=(ZA*A(3)+ZB*B(3))*InvGa
             PA(1)=P(1)-A(1)
             PA(2)=P(2)-A(2)
             PA(3)=P(3)-A(3)
             PB(1)=P(1)-B(1)
             PB(2)=P(2)-B(2)
             PB(3)=P(3)-B(3)
             !
             ! Compute OS RR.
             !
             CALL I1E_Ovlap_OSRecurs(PA(1),PB(1),Gamma,MaxLA+1,MaxLB)
             !
             ! Contract each buffer into appropriate location.
             !
             IB=IndexB
             DO LMNB=StartLB,StopLB
                LB=BS%LxDex%I(LMNB)
                MB=BS%LyDex%I(LMNB)
                NB=BS%LzDex%I(LMNB)
                CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                IB=IB+1
             IA=IndexA
             DO LMNA=StartLA,StopLA
                LA=BS%LxDex%I(LMNA)
                MA=BS%LyDex%I(LMNA)
                NA=BS%LzDex%I(LMNA)
                CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                IA=IA+1
                !
                X0=RRX(LA,LB)
                Y0=RRY(MA,MB)
                Z0=RRZ(NA,NB)
                X1=RRX(LA+1,LB)
                Y1=RRY(MA+1,MB)
                Z1=RRZ(NA+1,NB)
                OvT=Ov*CA*CB
                DX=-2D0*ZA*RRX(LA+1,LB)
                DY=-2D0*ZA*RRY(MA+1,MB)
                DZ=-2D0*ZA*RRZ(NA+1,NB)
                IF(LA.GT.0)DX=DX+DBLE(LA)*RRX(LA-1,LB)
                IF(MA.GT.0)DY=DY+DBLE(MA)*RRY(MA-1,MB)
                IF(NA.GT.0)DZ=DZ+DBLE(NA)*RRZ(NA-1,NB)
                IF(IXYZ.EQ.1)I1E(IA,IB)=I1E(IA,IB)+OvT*((Y1+AO(2)*Y0)*DZ-(Z1+AO(3)*Z0)*DY)*X0
                IF(IXYZ.EQ.2)I1E(IA,IB)=I1E(IA,IB)+OvT*((Z1+AO(3)*Z0)*DX-(X1+AO(1)*X0)*DZ)*Y0
                IF(IXYZ.EQ.3)I1E(IA,IB)=I1E(IA,IB)+OvT*((X1+AO(1)*X0)*DY-(Y1+AO(2)*Y0)*DX)*Z0
                !
             ENDDO
             ENDDO
             !
          ENDDO
#ifdef INT1E_THRESH
          IF(PFB.EQ.BS%NPFnc%I(CFB,KB)) EXIT
#endif
          ENDDO
          !
          IndexB=IndexB+StopLB-StartLB+1
       ENDDO
       IndexA=IndexA+StopLA-StartLA+1
    ENDDO
    !
  END SUBROUTINE I1E_KinMom_AtBlk
  !
  !
  SUBROUTINE I1E_DipVel_AtBlk(BS,Pair,IXYZ,I1E)
!---------------------------------------------------------------------------------
!H Compute the dipole velocity matrix:
!H $$
!H   (a|\nabla|b)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BSET)                               :: BS
    TYPE(AtomPair)                           :: Pair
    INTEGER                                  :: IXYZ
    REAL(DOUBLE), DIMENSION(Pair%NA,Pair%NB) :: I1E
    !-------------------------------------------------------------------
    INTEGER :: CFA,CFB,IndexA,IndexB,StartLA,StartLB
    INTEGER :: StopLA,StopLB,MaxLA,MaxLB,LMNA,LMNB
    INTEGER :: KA,KB,PFA,PFB,IA,IB,LA,LB,MA,MB,NA,NB
    REAL(DOUBLE), DIMENSION(3) :: P,A,B,PA,PB
    REAL(DOUBLE) :: AB2,ZA,ZB,CA,CB,Gamma,InvGa,Dum,Ov,OvT
    REAL(DOUBLE) :: X0,Y0,Z0,DX,DY,DZ
    !-------------------------------------------------------------------
    !
    KA=Pair%KA
    KB=Pair%KB
    A(1)=Pair%A(1)
    A(2)=Pair%A(2)
    A(3)=Pair%A(3)
    B(1)=Pair%B(1)
    B(2)=Pair%B(2)
    B(3)=Pair%B(3)
    AB2=Pair%AB2
    !
    ! Loop over contracted functions A and B.
    !
    IndexA=0
    DO CFA=1,BS%NCFnc%I(KA)
       StartLA = BS%LStrt%I(CFA,KA)
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       IndexB=0
       DO CFB=1,BS%NCFnc%I(KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          !
          ! Begin loops over primitives.
          !
          DO PFA=BS%NPFnc%I(CFA,KA),1,-1
             ZA=BS%Expnt%D(PFA,CFA,KA)
          DO PFB=BS%NPFnc%I(CFB,KB),1,-1
             ZB=BS%Expnt%D(PFB,CFB,KB)
             !
             Gamma=ZA+ZB
             InvGa=1.0D0/Gamma
             Dum=ZA*ZB*AB2*InvGa
#ifdef INT1E_THRESH
             IF(Dum>PrimPairDistanceThreshold) EXIT
#endif
             Ov=EXP(-Dum)*SQRT(PI*InvGa)*PI*InvGa
             !
             P(1)=(ZA*A(1)+ZB*B(1))*InvGa
             P(2)=(ZA*A(2)+ZB*B(2))*InvGa
             P(3)=(ZA*A(3)+ZB*B(3))*InvGa
             PA(1)=P(1)-A(1)
             PA(2)=P(2)-A(2)
             PA(3)=P(3)-A(3)
             PB(1)=P(1)-B(1)
             PB(2)=P(2)-B(2)
             PB(3)=P(3)-B(3)
             !
             ! Compute OS RR.
             !
             CALL I1E_Ovlap_OSRecurs(PA(1),PB(1),Gamma,MaxLA+1,MaxLB)
             !
             ! Contract each buffer into appropriate location.
             !
             IB=IndexB
             DO LMNB=StartLB,StopLB
                LB=BS%LxDex%I(LMNB)
                MB=BS%LyDex%I(LMNB)
                NB=BS%LzDex%I(LMNB)
                CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                IB=IB+1
             IA=IndexA
             DO LMNA=StartLA,StopLA
                LA=BS%LxDex%I(LMNA)
                MA=BS%LyDex%I(LMNA)
                NA=BS%LzDex%I(LMNA)
                CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                IA=IA+1
                !
                X0=RRX(LA,LB)
                Y0=RRY(MA,MB)
                Z0=RRZ(NA,NB)
                OvT=Ov*CA*CB
                DX=-2D0*ZA*RRX(LA+1,LB)
                DY=-2D0*ZA*RRY(MA+1,MB)
                DZ=-2D0*ZA*RRZ(NA+1,NB)
                IF(LA.GT.0)DX=DX+DBLE(LA)*RRX(LA-1,LB)
                IF(MA.GT.0)DY=DY+DBLE(MA)*RRY(MA-1,MB)
                IF(NA.GT.0)DZ=DZ+DBLE(NA)*RRZ(NA-1,NB)
                IF(IXYZ.EQ.1)I1E(IA,IB)=I1E(IA,IB)+OvT*DX*Y0*Z0
                IF(IXYZ.EQ.2)I1E(IA,IB)=I1E(IA,IB)+OvT*X0*DY*Z0
                IF(IXYZ.EQ.3)I1E(IA,IB)=I1E(IA,IB)+OvT*X0*Y0*DZ
                !
             ENDDO
             ENDDO
             !
          ENDDO
#ifdef INT1E_THRESH
          IF(PFB.EQ.BS%NPFnc%I(CFB,KB)) EXIT
#endif
          ENDDO
          !
          IndexB=IndexB+StopLB-StartLB+1
       ENDDO
       IndexA=IndexA+StopLA-StartLA+1
    ENDDO
    !
  END SUBROUTINE I1E_DipVel_AtBlk
  !
  !
  SUBROUTINE I1E_MssVel_AtBlk(BS,Pair,I1E)
!---------------------------------------------------------------------------------
!H Compute the mass velocity matrix:
!H $$
!H   (a|\nabla^2\cdot\nabla^2|b)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BSET)                               :: BS
    TYPE(AtomPair)                           :: Pair
    REAL(DOUBLE), DIMENSION(Pair%NA,Pair%NB) :: I1E
    !-------------------------------------------------------------------
    INTEGER :: CFA,CFB,IndexA,IndexB,StartLA,StartLB
    INTEGER :: StopLA,StopLB,MaxLA,MaxLB,LMNA,LMNB
    INTEGER :: KA,KB,PFA,PFB,IA,IB,LA,LB,MA,MB,NA,NB
    REAL(DOUBLE), DIMENSION(3) :: P,A,B,PA,PB
    REAL(DOUBLE) :: AB2,ZA,ZB,CA,CB,Gamma,InvGa,Dum,Ov,OvT,ZA2,ZA3,ZA4
    REAL(DOUBLE) :: X0,Y0,Z0,TX,TY,TZ
    !-------------------------------------------------------------------
    !
    KA=Pair%KA
    KB=Pair%KB
    A(1)=Pair%A(1)
    A(2)=Pair%A(2)
    A(3)=Pair%A(3)
    B(1)=Pair%B(1)
    B(2)=Pair%B(2)
    B(3)=Pair%B(3)
    AB2=Pair%AB2
    !
    ! Loop over contracted functions A and B.
    !
    IndexA=0
    DO CFA=1,BS%NCFnc%I(KA)
       StartLA = BS%LStrt%I(CFA,KA)
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       IndexB=0
       DO CFB=1,BS%NCFnc%I(KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          !
          ! Begin loops over primitives.
          !
          DO PFA=BS%NPFnc%I(CFA,KA),1,-1
             ZA=BS%Expnt%D(PFA,CFA,KA)
             ZA2=ZA *ZA
             ZA3=ZA2*ZA
             ZA4=ZA3*ZA
          DO PFB=BS%NPFnc%I(CFB,KB),1,-1
             ZB=BS%Expnt%D(PFB,CFB,KB)
             !
             Gamma=ZA+ZB
             InvGa=1.0D0/Gamma
             Dum=ZA*ZB*AB2*InvGa
#ifdef INT1E_THRESH
             IF(Dum>PrimPairDistanceThreshold) EXIT
#endif
             Ov=EXP(-Dum)*SQRT(PI*InvGa)*PI*InvGa
             !
             P(1)=(ZA*A(1)+ZB*B(1))*InvGa
             P(2)=(ZA*A(2)+ZB*B(2))*InvGa
             P(3)=(ZA*A(3)+ZB*B(3))*InvGa
             PA(1)=P(1)-A(1)
             PA(2)=P(2)-A(2)
             PA(3)=P(3)-A(3)
             PB(1)=P(1)-B(1)
             PB(2)=P(2)-B(2)
             PB(3)=P(3)-B(3)
             !
             ! Compute OS RR.
             !
             CALL I1E_Ovlap_OSRecurs(PA(1),PB(1),Gamma,MaxLA+4,MaxLB)
             !
             ! Contract each buffer into appropriate location.
             !
             IB=IndexB
             DO LMNB=StartLB,StopLB
                LB=BS%LxDex%I(LMNB)
                MB=BS%LyDex%I(LMNB)
                NB=BS%LzDex%I(LMNB)
                CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                IB=IB+1
             IA=IndexA
             DO LMNA=StartLA,StopLA
                LA=BS%LxDex%I(LMNA)
                MA=BS%LyDex%I(LMNA)
                NA=BS%LzDex%I(LMNA)
                CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                IA=IA+1
                !
                X0=RRX(LA,LB)
                Y0=RRY(MA,MB)
                Z0=RRZ(NA,NB)
                OvT=Ov*CA*CB
                TX=16D0*ZA4*RRX(LA+4,LB)-16D0*ZA3*DBLE(2*LA+3)*RRX(LA+2,LB) &
                       +12D0*ZA2*DBLE(LA**2+(LA+1)**2)*X0
                IF(LA.GT.1) THEN
                   TX=TX-4D0*ZA*DBLE(LA*(LA-1)*(2*LA-1))*RRX(LA-2,LB)
                   IF(LA.GT.3) TX=TX+DBLE(LA*(LA-1)*(LA-2)*(LA-3))*RRX(LA-4,LB)
                ENDIF
                TY=16D0*ZA4*RRY(MA+4,MB)-16D0*ZA3*DBLE(2*MA+3)*RRY(MA+2,MB) &
                       +12D0*ZA2*DBLE(MA**2+(MA+1)**2)*Y0
                IF(MA.GT.1) THEN
                   TY=TY-4D0*ZA*DBLE(MA*(MA-1)*(2*MA-1))*RRY(MA-2,MB)
                   IF(MA.GT.3) TY=TY+DBLE(MA*(MA-1)*(MA-2)*(MA-3))*RRY(MA-4,MB)
                ENDIF
                TZ=16D0*ZA4*RRZ(NA+4,NB)-16D0*ZA3*DBLE(2*NA+3)*RRZ(NA+2,NB) &
                       +12D0*ZA2*DBLE(NA**2+(NA+1)**2)*Z0
                IF(NA.GT.1) THEN
                   TZ=TZ-4D0*ZA*DBLE(NA*(NA-1)*(2*NA-1))*RRZ(NA-2,NB)
                   IF(NA.GT.3) TZ=TZ+DBLE(NA*(NA-1)*(NA-2)*(NA-3))*RRZ(NA-4,NB)
                ENDIF
                I1E(IA,IB)=I1E(IA,IB)+OvT*(TX*Y0*Z0+X0*TY*Z0+X0*Y0*TZ)
                !
             ENDDO
             ENDDO
             !
          ENDDO
#ifdef INT1E_THRESH
          IF(PFB.EQ.BS%NPFnc%I(CFB,KB)) EXIT
#endif
          ENDDO
          !
          IndexB=IndexB+StopLB-StartLB+1
       ENDDO
       IndexA=IndexA+StopLA-StartLA+1
    ENDDO
    !
  END SUBROUTINE I1E_MssVel_AtBlk
  !
  !
  SUBROUTINE I1E_NucAtt_AtBlk(GM,BS,Pair,I1E)
!---------------------------------------------------------------------------------
!H Compute the nuclear attraction matrix:
!H $$
!H   \sum_{k=1}^{N_{atm}}(a|\frac{Z_k}{|{\bf r-r}_k|}|b)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(CRDS)                               :: GM
    TYPE(BSET)                               :: BS
    TYPE(AtomPair)                           :: Pair
    REAL(DOUBLE), DIMENSION(Pair%NA,Pair%NB) :: I1E
    !-------------------------------------------------------------------
    INTEGER :: CFA,CFB,IndexA,IndexB,StartLA,StartLB
    INTEGER :: StopLA,StopLB,MaxLA,MaxLB,LMNA,LMNB
    INTEGER :: KA,KB,PFA,PFB,IA,IB,LA,LB,MA,MB,NA,NB,At
    INTEGER :: OAX,OAY,OAZ,OBX,OBY,OBZ,IdxA,IdxB
    REAL(DOUBLE), DIMENSION(3) :: P,A,B,PA,PB,PC
    REAL(DOUBLE) :: AB2,ZA,ZB,CA,CB,Gamma,InvGa,Dum,Ov,Z
    REAL(DOUBLE) :: X0,Y0,Z0,DX,DY,DZ
    !-------------------------------------------------------------------
    !
    KA=Pair%KA
    KB=Pair%KB
    A(1)=Pair%A(1)
    A(2)=Pair%A(2)
    A(3)=Pair%A(3)
    B(1)=Pair%B(1)
    B(2)=Pair%B(2)
    B(3)=Pair%B(3)
    AB2=Pair%AB2
    !
    ! Loop over contracted functions A and B.
    !
    IndexA=0
    DO CFA=1,BS%NCFnc%I(KA)
       StartLA = BS%LStrt%I(CFA,KA)
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       OAZ=1
       OAY=MaxLA+1
       OAX=OAY*OAY
       IndexB=0
       DO CFB=1,BS%NCFnc%I(KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          OBZ=1
          OBY=MaxLB+1
          OBX=OBY*OBY
          !
          ! Begin loops over primitives.
          !
          DO PFA=BS%NPFnc%I(CFA,KA),1,-1
             ZA=BS%Expnt%D(PFA,CFA,KA)
          DO PFB=BS%NPFnc%I(CFB,KB),1,-1
             ZB=BS%Expnt%D(PFB,CFB,KB)
             !
             Gamma=ZA+ZB
             InvGa=1.0D0/Gamma
             Dum=ZA*ZB*AB2*InvGa
#ifdef INT1E_THRESH
             IF(Dum>PrimPairDistanceThreshold) EXIT
#endif
             Ov=EXP(-Dum)*SQRT(PI*InvGa)*PI*InvGa
             !
             P(1)=(ZA*A(1)+ZB*B(1))*InvGa
             P(2)=(ZA*A(2)+ZB*B(2))*InvGa
             P(3)=(ZA*A(3)+ZB*B(3))*InvGa
             PA(1)=P(1)-A(1)
             PA(2)=P(2)-A(2)
             PA(3)=P(3)-A(3)
             PB(1)=P(1)-B(1)
             PB(2)=P(2)-B(2)
             PB(3)=P(3)-B(3)
             !
             ! Begin the loop over atoms.
             !
             DO At=1,NAtoms
                !
                Z=GM%AtNum%D(At)
                PC(1)=P(1)-GM%Carts%D(1,At)
                PC(2)=P(2)-GM%Carts%D(2,At)
                PC(3)=P(3)-GM%Carts%D(3,At)
                !
                ! Compute OS RR.
                !
                CALL I1E_NucAtt_OSRecurs(PA(1),PB(1),PC(1),Gamma,MaxLA,MaxLB)
                !
                ! Contract each buffer into appropriate location.
                !
                IB=IndexB
                DO LMNB=StartLB,StopLB
                   LB=BS%LxDex%I(LMNB)
                   MB=BS%LyDex%I(LMNB)
                   NB=BS%LzDex%I(LMNB)
                   CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                   IB=IB+1
                   IdxB=LB*OBX+MB*OBY+NB*OBZ
                IA=IndexA
                DO LMNA=StartLA,StopLA
                   LA=BS%LxDex%I(LMNA)
                   MA=BS%LyDex%I(LMNA)
                   NA=BS%LzDex%I(LMNA)
                   CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                   IA=IA+1
                   IdxA=LA*OAX+MA*OAY+NA*OAZ
                   !
                   I1E(IA,IB)=I1E(IA,IB)-Ov*CA*CB*Z*RR0(0,IdxA,IdxB)
                   !
                ENDDO
                ENDDO
                !
             ENDDO
             !
          ENDDO
#ifdef INT1E_THRESH
          IF(PFB.EQ.BS%NPFnc%I(CFB,KB)) EXIT
#endif
          ENDDO
          !
          IndexB=IndexB+StopLB-StartLB+1
       ENDDO
       IndexA=IndexA+StopLA-StartLA+1
    ENDDO
    !
  END SUBROUTINE I1E_NucAtt_AtBlk
  !
  !
  SUBROUTINE I1E_EField_AtBlk()
!---------------------------------------------------------------------------------
!H Compute the electric field matrix:
!H $$
!H   \sum_{k=1}^{N_{atm}}(a|\frac{{\bf r-r}_k}{|{\bf r-r}_k|^3}|b)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE

  END SUBROUTINE I1E_EField_AtBlk
  !
  !
  SUBROUTINE I1E_EFieldGrd_AtBlk()
    IMPLICIT NONE

  END SUBROUTINE I1E_EFieldGrd_AtBlk
  !
  !
  SUBROUTINE I1E_rKinet_AtBlk(BS,Pair,IXYZ,I1E)
!---------------------------------------------------------------------------------
!H Compute the modified kinetic matrix:
!H $$
!H   -\frac{1}{2}(a|{\bf r}\nabla^2|b)
!H $$
!---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BSET)                               :: BS
    TYPE(AtomPair)                           :: Pair
    INTEGER                                  :: IXYZ
    REAL(DOUBLE), DIMENSION(Pair%NA,Pair%NB) :: I1E
    !-------------------------------------------------------------------
    INTEGER :: CFA,CFB,IndexA,IndexB,StartLA,StartLB
    INTEGER :: StopLA,StopLB,MaxLA,MaxLB,LMNA,LMNB
    INTEGER :: KA,KB,PFA,PFB,IA,IB,LA,LB,MA,MB,NA,NB
    REAL(DOUBLE), DIMENSION(3) :: P,A,B,PA,PB
    REAL(DOUBLE) :: AB2,ZA,ZB,CA,CB,Gamma,InvGa,Dum,Ov,OvT,ZA2
    REAL(DOUBLE) :: X0,Y0,Z0,X1,Y1,Z1,TX0,TY0,TZ0,TX1,TY1,TZ1
    !-------------------------------------------------------------------
    !
    KA=Pair%KA
    KB=Pair%KB
    A(1)=Pair%A(1)
    A(2)=Pair%A(2)
    A(3)=Pair%A(3)
    B(1)=Pair%B(1)
    B(2)=Pair%B(2)
    B(3)=Pair%B(3)
    AB2=Pair%AB2
    !
    ! Loop over contracted functions A and B.
    !
    IndexA=0
    DO CFA=1,BS%NCFnc%I(KA)
       StartLA = BS%LStrt%I(CFA,KA)
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       IndexB=0
       DO CFB=1,BS%NCFnc%I(KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          !
          ! Begin loops over primitives.
          !
          DO PFA=BS%NPFnc%I(CFA,KA),1,-1
             ZA=BS%Expnt%D(PFA,CFA,KA)
             ZA2=ZA*ZA
          DO PFB=BS%NPFnc%I(CFB,KB),1,-1
             ZB=BS%Expnt%D(PFB,CFB,KB)
             !
             Gamma=ZA+ZB
             InvGa=1.0D0/Gamma
             Dum=ZA*ZB*AB2*InvGa
#ifdef INT1E_THRESH
             IF(Dum>PrimPairDistanceThreshold) EXIT
#endif
             Ov=EXP(-Dum)*SQRT(PI*InvGa)*PI*InvGa
             !
             P(1)=(ZA*A(1)+ZB*B(1))*InvGa
             P(2)=(ZA*A(2)+ZB*B(2))*InvGa
             P(3)=(ZA*A(3)+ZB*B(3))*InvGa
             PA(1)=P(1)-A(1)
             PA(2)=P(2)-A(2)
             PA(3)=P(3)-A(3)
             PB(1)=P(1)-B(1)
             PB(2)=P(2)-B(2)
             PB(3)=P(3)-B(3)
             !
             ! Compute OS RR.
             !
             CALL I1E_Ovlap_OSRecurs(PA(1),PB(1),Gamma,MaxLA+3,MaxLB)
             !
             ! Contract each buffer into appropriate location.
             !
             IB=IndexB
             DO LMNB=StartLB,StopLB
                LB=BS%LxDex%I(LMNB)
                MB=BS%LyDex%I(LMNB)
                NB=BS%LzDex%I(LMNB)
                CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                IB=IB+1
             IA=IndexA
             DO LMNA=StartLA,StopLA
                LA=BS%LxDex%I(LMNA)
                MA=BS%LyDex%I(LMNA)
                NA=BS%LzDex%I(LMNA)
                CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                IA=IA+1
                !
                X0=RRX(LA,LB)
                Y0=RRY(MA,MB)
                Z0=RRZ(NA,NB)
                X1=RRX(LA+1,LB)
                Y1=RRY(MA+1,MB)
                Z1=RRZ(NA+1,NB)
                OvT=Ov*CA*CB
                TX0=ZA*DBLE(2*LA+1)*RRX(LA,LB)-2D0*ZA*ZA*RRX(LA+2,LB)
                TY0=ZA*DBLE(2*MA+1)*RRY(MA,MB)-2D0*ZA*ZA*RRY(MA+2,MB)
                TZ0=ZA*DBLE(2*NA+1)*RRZ(NA,NB)-2D0*ZA*ZA*RRZ(NA+2,NB)
                IF(LA.GT.1) TX0=TX0-0.5D0*DBLE(LA*(LA-1))*RRX(LA-2,LB)
                IF(MA.GT.1) TY0=TY0-0.5D0*DBLE(MA*(MA-1))*RRY(MA-2,MB)
                IF(NA.GT.1) TZ0=TZ0-0.5D0*DBLE(NA*(NA-1))*RRZ(NA-2,NB)
                TX1=-2D0*ZA2*RRX(LA+3,LB)+ZA*DBLE(2*LA+3)*X1
                TY1=-2D0*ZA2*RRY(MA+3,MB)+ZA*DBLE(2*MA+3)*Y1
                TZ1=-2D0*ZA2*RRZ(NA+3,NB)+ZA*DBLE(2*NA+3)*Z1
                IF(LA.GT.0) TX1=TX1-0.5D0*DBLE(LA*(LA+1))*RRX(LA-1,LB)
                IF(MA.GT.0) TY1=TY1-0.5D0*DBLE(MA*(MA+1))*RRY(MA-1,MB)
                IF(NA.GT.0) TZ1=TZ1-0.5D0*DBLE(NA*(NA+1))*RRZ(NA-1,NB)
                IF(IXYZ.EQ.1)I1E(IA,IB)=I1E(IA,IB)+OvT*(TX1*Y0*Z0+X1*TY0*Z0+X1*Y0*TZ0+ &
                                                  A(1)*(TX0*Y0*Z0+X0*TY0*Z0+X0*Y0*TZ0))
                IF(IXYZ.EQ.2)I1E(IA,IB)=I1E(IA,IB)+OvT*(TX0*Y1*Z0+X0*TY1*Z0+X0*Y1*TZ0+ &
                                                  A(2)*(TX0*Y0*Z0+X0*TY0*Z0+X0*Y0*TZ0))
                IF(IXYZ.EQ.3)I1E(IA,IB)=I1E(IA,IB)+OvT*(TX0*Y0*Z1+X0*TY0*Z1+X0*Y0*TZ1+ &
                                                  A(3)*(TX0*Y0*Z0+X0*TY0*Z0+X0*Y0*TZ0))
                !
             ENDDO
             ENDDO
             !
          ENDDO
#ifdef INT1E_THRESH
          IF(PFB.EQ.BS%NPFnc%I(CFB,KB)) EXIT
#endif
          ENDDO
          !
          IndexB=IndexB+StopLB-StartLB+1
       ENDDO
       IndexA=IndexA+StopLA-StartLA+1
    ENDDO
    !
  END SUBROUTINE I1E_rKinet_AtBlk
  !
  !
  SUBROUTINE I1E_Init()
    IMPLICIT NONE

  END SUBROUTINE I1E_Init
  !
  !
  SUBROUTINE I1E_Finalize()
    IMPLICIT NONE

  END SUBROUTINE I1E_Finalize


END MODULE Int1E
