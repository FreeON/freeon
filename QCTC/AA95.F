C     -----------------------------------------------------------------
C     Reza Ahmadi and Jan Almlof
C     Chemical Physics Letters; 1 Dec. 1995; vol.246, no.4-5, p.364-70
C     Buffed by M. Challacombe
C     -----------------------------------------------------------------
      SUBROUTINE RAhmadiJAlmlof95(BigEll,Nc,EllP,EllQ,LenP,LenQ,
     >                            PQEll,PQLen,NukeX,Threshold,QEst,
     >                            ZetaP,ZetaQ,P,Q,QCo,Ket)
C
      IMPLICIT NONE
C
      INTEGER BigN,BigEll,Nc,EllP,EllQ,EllPQ,LenP,LenQ,LenPQ,PQEll,PQLen
      INTEGER I,J,K,L,NcDim
      REAL*8  P(3),Q(3,*),ZetaP,ZetaQ(*),NukeX,Threshold
      REAL*8  QEst(*),QCo(LenQ,*),Ket(*)
C
      INTEGER iPQx,iPQy,iPQz,iAux,iArr,iEnd
C
      INTEGER iClust,MxEll,MxLen,iWork
      PARAMETER (iClust=511)
      PARAMETER (MxEll=13)
      PARAMETER (MxLen=(MxEll+1)*(MxEll+2)*(MxEll+3)/6)
      PARAMETER (iWork=3*iClust+iClust*(MxLen+MxEll+1))
      REAL*8 Wrk(iWork)
C
      EllPQ=EllP+EllQ
      LenPQ=LenP*LenQ
C
C     FIND ONLY DISTRIBUTIONS THAT SATISFY THE INTEGRAL INEQUALITY
c      IF(QEst(1).LT.Threshold)THEN
c         RETURN
c      ELSEIF(QEst(Nc).LT.Threshold)THEN
c         K=0
c         L=Nc
c         Nc=Nc/2
c         DO WHILE(K.LT.Nc.AND.L.GT.Nc)
c            IF(QEst(Nc).LE.Threshold)THEN
c               L=Nc
c            ELSE
c               K=Nc
c            ENDIF
c            Nc=(K+L)/2
c         ENDDO
c      ENDIF
C
      iPQx=1
      iPQy=iPQx+Nc
      iPQz=iPQy+Nc
      iAux=iPQz+Nc
      iArr=iAux+Nc*(PQEll+1)
      iEnd=iArr+Nc*(PQLen+1)
C
      IF(iEnd>iWork)
     >   STOP 'Memory exceeded in AA95 '
C
C     CHECK TO SEE IF WE ARE DOING NUCLEAR-TOTAL OR ELECTRON-TOTAL INTEGRALS
      IF(ZetaP.EQ.NukeX.AND.EllQ==0)THEN
C        Special purpose (slow) code that avoids nuclear self interaction
         IF(EllP==0)THEN
            CALL ENuE0(Nc,P,Q,ZetaP,ZetaQ,QCo(1,1),Ket)
         ELSEIF(EllP==1)THEN
            CALL ENuE1(Nc,P,Q,ZetaP,ZetaQ,QCo,Ket)
         ELSE
            WRITE(*,*)' ZetaP = ',ZetaP
            WRITE(*,*)' STOPPED OUTWITH ELLP = ',ELLP,' ELLQ = ',ELLQ
            STOP 'No second or higher nuclear derivatives in AA.F'
         ENDIF
      ELSE
C
         SELECT CASE(EllPQ)
         CASE(0)
#ifdef PAPIEX
            CALL PAPIEX_START(1,"< AuxIGen_0 > ")
#endif
            CALL AuxIGen0(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
#ifdef PAPIEX
            CALL PAPIEX_STOP(1)
            CALL PAPIEX_START(2,"< MD3TRR_0 > ")
#endif
            CALL MD3TRR0(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))


#ifdef PAPIEX
            CALL PAPIEX_STOP(2)
            CALL PAPIEX_START(3,"< KTrax_0 > ")
#endif
            CALL KTrax_0_0(Nc,Wrk(iArr),QCo,Ket)
#ifdef PAPIEX
            CALL PAPIEX_STOP(3)
#endif
         CASE(1)
#ifdef PAPIEX
            CALL PAPIEX_START(4,"< AuxIGen_1 >")
#endif
            CALL AuxIGen1(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
#ifdef PAPIEX
            CALL PAPIEX_STOP(4)
            CALL PAPIEX_START(5,"< MD3TRR_1 > ")
#endif
            CALL MD3TRR1(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
#ifdef PAPIEX
            CALL PAPIEX_STOP(5)
            CALL PAPIEX_START(6,"< KTrax_1 > ")
#endif
            IF(EllP==0)THEN
               CALL KTrax_0_1(Nc,Wrk(iArr),QCo,Ket)
            ELSE
               CALL KTrax_1_0(Nc,Wrk(iArr),QCo,Ket)
            ENDIF
#ifdef PAPIEX
            CALL PAPIEX_STOP(6)
#endif
         CASE(2)
#ifdef PAPIEX
            CALL PAPIEX_START(7,"< AuxIGen_2 > ")
#endif
            CALL AuxIGen2(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
#ifdef PAPIEX
            CALL PAPIEX_STOP(7)
            CALL PAPIEX_START(8,"< MD3TRR_2 >")
#endif
            CALL MD3TRR2(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
#ifdef PAPIEX
            CALL PAPIEX_STOP(8)
            CALL PAPIEX_START(9,"< KTrax_2 > ")
#endif
            IF(EllP==0)THEN
               CALL KTrax_0_2(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_1(Nc,Wrk(iArr),QCo,Ket)
            ELSE
               CALL KTrax_2_0(Nc,Wrk(iArr),QCo,Ket)
            ENDIF
#ifdef PAPIEX
            CALL PAPIEX_STOP(9)
#endif
         CASE(3)
#ifdef PAPIEX
            CALL PAPIEX_START(10,"< AuxIGen_3 > ")
#endif
            CALL AuxIGen3(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
#ifdef PAPIEX
            CALL PAPIEX_STOP(10)
            CALL PAPIEX_START(11,"< MD3TRR_3 > ")
#endif
            CALL MD3TRR3(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
#ifdef PAPIEX
            CALL PAPIEX_STOP(11)
            CALL PAPIEX_START(12,"< KTrax_3 > ")
#endif
            IF(EllP==0)THEN
               CALL KTrax_0_3(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_2(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_1(Nc,Wrk(iArr),QCo,Ket)
            ELSE
               CALL KTrax_3_0(Nc,Wrk(iArr),QCo,Ket)
            ENDIF
#ifdef PAPIEX
            CALL PAPIEX_STOP(12)
#endif
         CASE(4)
#ifdef PAPIEX
            CALL PAPIEX_START(13,"< AuxIGen_4 > ")
#endif
            CALL AuxIGen4(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
#ifdef PAPIEX
            CALL PAPIEX_STOP(13)
            CALL PAPIEX_START(14,"< MD3TRR_4 > ")
#endif
            CALL MD3TRR4(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
#ifdef PAPIEX
            CALL PAPIEX_STOP(14)
            CALL PAPIEX_START(15,"< KTrax_4 > ")
#endif
            IF(EllP==0)THEN
               CALL KTrax_0_4(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_3(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_2(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==3)THEN
               CALL KTrax_3_1(Nc,Wrk(iArr),QCo,Ket)
            ELSE
               CALL KTrax_4_0(Nc,Wrk(iArr),QCo,Ket)
            ENDIF
#ifdef PAPIEX
            CALL PAPIEX_STOP(15)
#endif
         CASE(5)
#ifdef PAPIEX
            CALL PAPIEX_START(16,"< AuxIGen_5 > ")
#endif
            CALL AuxIGen5(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
#ifdef PAPIEX
            CALL PAPIEX_STOP(16)
            CALL PAPIEX_START(17,"< MD3TRR_5 > ")
#endif
            CALL MD3TRR5(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
#ifdef PAPIEX
            CALL PAPIEX_STOP(17)
            CALL PAPIEX_START(18,"< KTrax_5 > ")
#endif
            IF(EllP==0)THEN
               CALL KTrax_0_5(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_4(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_3(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==3)THEN
               CALL KTrax_3_2(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==4)THEN
               CALL KTrax_4_1(Nc,Wrk(iArr),QCo,Ket)
            ELSE
               CALL KTrax_5_0(Nc,Wrk(iArr),QCo,Ket)
            ENDIF
#ifdef PAPIEX
            CALL PAPIEX_STOP(18)
#endif
         CASE(6)
#ifdef PAPIEX
            CALL PAPIEX_START(19,"< AuxIGen_6 > ")
#endif
            CALL AuxIGen6(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
#ifdef PAPIEX
            CALL PAPIEX_STOP(19)
            CALL PAPIEX_START(20,"< MD3TRR_6 > ")
#endif
            CALL MD3TRR6(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
#ifdef PAPIEX
            CALL PAPIEX_STOP(20)
            CALL PAPIEX_START(21,"< KTrax_6 > ")
#endif
            IF(EllP==0)THEN
               CALL KTrax_0_6(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_5(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_4(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==3)THEN
               CALL KTrax_3_3(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==4)THEN
               CALL KTrax_4_2(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==5)THEN
               CALL KTrax_5_1(Nc,Wrk(iArr),QCo,Ket)
            ELSE
               CALL KTrax_6_0(Nc,Wrk(iArr),QCo,Ket)
            ENDIF
#ifdef PAPIEX
            CALL PAPIEX_STOP(21)
#endif
         CASE(7)
#ifdef PAPIEX
            CALL PAPIEX_START(22,"< AuxIGen_7 > ")
#endif
            CALL AuxIGen7(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),Wrk(iAux))
#ifdef PAPIEX
            CALL PAPIEX_STOP(22)
            CALL PAPIEX_START(23,"< MD3TRR_7 > ")
#endif
            CALL MD3TRR7(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
#ifdef PAPIEX
            CALL PAPIEX_STOP(23)
            CALL PAPIEX_START(24,"< KTrax_7 > ")
#endif
            IF(EllP==0)THEN
               STOP '(0|3+4) in AA95 '
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_6(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_5(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==3)THEN
               CALL KTrax_3_4(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==4)THEN
               CALL KTrax_4_3(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==5)THEN
               CALL KTrax_5_2(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==6)THEN
               CALL KTrax_6_1(Nc,Wrk(iArr),QCo,Ket)
            ELSE
               CALL KTrax_7_0(Nc,Wrk(iArr),QCo,Ket)
            ENDIF
#ifdef PAPIEX
            CALL PAPIEX_STOP(24)
#endif
         CASE(8)
#ifdef PAPIEX
            CALL PAPIEX_START(25,"< AuxIGen_8 > ")
#endif
            CALL AuxIGen8(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
#ifdef PAPIEX
            CALL PAPIEX_STOP(25)
            CALL PAPIEX_START(26,"< MD3TRR_8 > ")
#endif
            CALL MD3TRR8(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
#ifdef PAPIEX
            CALL PAPIEX_STOP(26)
            CALL PAPIEX_START(27,"< KTrax_8 > ")
#endif
            IF(EllP==0)THEN
               STOP '(0|8) in AA95 '
            ELSEIF(EllP==1)THEN
               STOP '(1|7) in AA95 '
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_6(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==3)THEN
               CALL KTrax_3_5(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==4)THEN
               CALL KTrax_4_4(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==5)THEN
               CALL KTrax_5_3(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==6)THEN
               CALL KTrax_6_2(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==7)THEN
               CALL KTrax_7_1(Nc,Wrk(iArr),QCo,Ket)
            ELSE
               STOP '(8+|0) in AA95'
            ENDIF
#ifdef PAPIEX
            CALL PAPIEX_STOP(27)
#endif
         CASE(9)
#ifdef PAPIEX
            CALL PAPIEX_START(25,"< AuxIGen_9 > ")
#endif
            CALL AuxIGen9(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
#ifdef PAPIEX
            CALL PAPIEX_STOP(25)
            CALL PAPIEX_START(26,"< MD3TRR_8 > ")
#endif
            CALL MD3TRR9(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
#ifdef PAPIEX
            CALL PAPIEX_STOP(26)
            CALL PAPIEX_START(27,"< KTrax_9 > ")
#endif
            IF(EllP==0)THEN
               STOP '(0|9) in AA95 '
            ELSEIF(EllP==1)THEN
               STOP '(1|8) in AA95 '
            ELSEIF(EllP==2)THEN
               STOP '(2|7) in AA95 '
            ELSEIF(EllP==3)THEN
               CALL KTrax_3_6(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==4)THEN
               CALL KTrax_4_5(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==5)THEN
               CALL KTrax_5_4(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==6)THEN
               CALL KTrax_6_3(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==7)THEN
               CALL KTrax_7_2(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==8)THEN
               STOP '(8|1) in AA95 '
            ELSE
               STOP '(9+|0) in AA95'
            ENDIF
#ifdef PAPIEX
            CALL PAPIEX_STOP(27)
#endif
         CASE DEFAULT

            WRITE(*,*)' STOPPED OUTWITH ELLP = ',ELLP,' ELLQ = ',ELLQ,
     >                ' ELLPQ = ',ELLPQ
            STOP 'AA95 '
         END SELECT
C
      ENDIF
C
      RETURN
      END
C
#ifdef OVERLAP_PROJECTION
      SUBROUTINE AAOverlap(BigEll,NcDim,EllP,EllQ,LenP,LenQ,
     >                            PQEll,PQLen,NukeX,Threshold,QEst,
     >                            ZetaP,ZetaQ,P,Q,QCo,Ket)
C
      IMPLICIT NONE
C
      INTEGER BigN,BigEll,Nc,EllP,EllQ,EllPQ,LenP,LenQ,LenPQ,PQEll,PQLen
      INTEGER I,J,K,L,NcDim
      REAL*8  P(3),Q(3,*),ZetaP,ZetaQ(*),NukeX,Threshold
      REAL*8  QEst(*),QCo(LenQ,*),Ket(*)
C
      INTEGER iPQx,iPQy,iPQz,iAux,iArr,iEnd
C
      INTEGER iClust,MxEll,MxLen,iWork
      PARAMETER (iClust=511)
      PARAMETER (MxEll=13)
      PARAMETER (MxLen=(MxEll+1)*(MxEll+2)*(MxEll+3)/6)
      PARAMETER (iWork=3*iClust+iClust*(MxLen+MxEll+1))
      REAL*8 Wrk(iWork)
C
      Nc=NcDim
      EllPQ=EllP+EllQ
      LenPQ=LenP*LenQ
C
      iPQx=1
      iPQy=iPQx+Nc
      iPQz=iPQy+Nc
      iAux=iPQz+Nc
      iArr=iAux+Nc*(PQEll+1)
      iEnd=iArr+Nc*(PQLen+1)
C
      IF(iEnd>iWork)
     >   STOP 'Memory exceeded in AA95 '


      SELECT CASE(EllPQ)
      CASE(0)
         CALL OvrIGen0(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),
     >        Wrk(iPQz),Wrk(iAux))
         CALL MD3TRR0(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >        Wrk(iAux),Wrk(iArr))
         CALL KTrax_0_0(Nc,Wrk(iArr),QCo,Ket)
      CASE(1)
            CALL OvrIGen1(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
            CALL MD3TRR1(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
            IF(EllP==0)THEN
               CALL KTrax_0_1(Nc,Wrk(iArr),QCo,Ket)
            ELSE
               CALL KTrax_1_0(Nc,Wrk(iArr),QCo,Ket)
            ENDIF
         CASE(2)
            CALL OvrIGen2(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
            CALL MD3TRR2(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
            IF(EllP==0)THEN
               CALL KTrax_0_2(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_1(Nc,Wrk(iArr),QCo,Ket)
            ELSE
               CALL KTrax_2_0(Nc,Wrk(iArr),QCo,Ket)
            ENDIF
         CASE(3)
            CALL OvrIGen3(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
            CALL MD3TRR3(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
            IF(EllP==0)THEN
               CALL KTrax_0_3(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_2(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_1(Nc,Wrk(iArr),QCo,Ket)
            ELSE
               CALL KTrax_3_0(Nc,Wrk(iArr),QCo,Ket)
            ENDIF
         CASE(4)
            CALL OvrIGen4(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
            CALL MD3TRR4(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
            IF(EllP==0)THEN
               CALL KTrax_0_4(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_3(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_2(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==3)THEN
               CALL KTrax_3_1(Nc,Wrk(iArr),QCo,Ket)
            ELSE
               CALL KTrax_4_0(Nc,Wrk(iArr),QCo,Ket)
            ENDIF
         CASE(5)
            CALL OvrIGen5(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
            CALL MD3TRR5(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
            IF(EllP==0)THEN
               CALL KTrax_0_5(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_4(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_3(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==3)THEN
               CALL KTrax_3_2(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==4)THEN
               CALL KTrax_4_1(Nc,Wrk(iArr),QCo,Ket)
            ELSE
               CALL KTrax_5_0(Nc,Wrk(iArr),QCo,Ket)
            ENDIF
         CASE(6)
            CALL OvrIGen6(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
            CALL MD3TRR6(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
            IF(EllP==0)THEN
               CALL KTrax_0_6(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_5(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_4(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==3)THEN
               CALL KTrax_3_3(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==4)THEN
               CALL KTrax_4_2(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==5)THEN
               CALL KTrax_5_1(Nc,Wrk(iArr),QCo,Ket)
            ELSE
               CALL KTrax_6_0(Nc,Wrk(iArr),QCo,Ket)
            ENDIF
         CASE(7)
            CALL OvrIGen7(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),Wrk(iAux))
            CALL MD3TRR7(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
            IF(EllP==0)THEN
               STOP '(0|3+4) in AA95 '
            ELSEIF(EllP==1)THEN
               CALL KTrax_1_6(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_5(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==3)THEN
               CALL KTrax_3_4(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==4)THEN
               CALL KTrax_4_3(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==5)THEN
               CALL KTrax_5_2(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==6)THEN
               CALL KTrax_6_1(Nc,Wrk(iArr),QCo,Ket)
            ELSE
               CALL KTrax_7_0(Nc,Wrk(iArr),QCo,Ket)
            ENDIF
         CASE(8)
            CALL OvrIGen8(Nc,P,Q,ZetaP,ZetaQ,Wrk(iPQx),Wrk(iPQy),
     >                    Wrk(iPQz),Wrk(iAux))
            CALL MD3TRR8(Nc,Wrk(iPQx),Wrk(iPQy),Wrk(iPQz),
     >                   Wrk(iAux),Wrk(iArr))
            IF(EllP==0)THEN
               STOP '(0|8) in AA95 '
            ELSEIF(EllP==1)THEN
               STOP '(1|8) in AA95 '
            ELSEIF(EllP==2)THEN
               CALL KTrax_2_6(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==3)THEN
               CALL KTrax_3_5(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==4)THEN
               CALL KTrax_4_4(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==5)THEN
               CALL KTrax_5_3(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==6)THEN
               CALL KTrax_6_2(Nc,Wrk(iArr),QCo,Ket)
            ELSEIF(EllP==7)THEN
               CALL KTrax_7_1(Nc,Wrk(iArr),QCo,Ket)
            ELSE
               STOP '(8+|0) in AA95'
            ENDIF
         CASE DEFAULT
            STOP 'AA95 '
         END SELECT
C
      RETURN
      END
#endif
C
      SUBROUTINE IBounds(Ell,Len,Nq,Zt,Co,Pndex,Qndex,PQndex,Est)
C
      INTEGER Ell,Len,Nq
C
      REAL*8 Zt,Co(*),Est(*)
      INTEGER Pndex(*),Qndex(*),PQndex(*)
      REAL*8 Omega,Upq,Phase,o1,o2
      INTEGER TwoEll,I,J,IP,VP,IQ,VQ,IPQ
      INTEGER L,M,N,LMN,I0,I1,I2
C
      INTEGER MxEll,MxLen
      PARAMETER (MxEll=13)
      PARAMETER (MxLen=(MxEll+1)*(MxEll+2)*(MxEll+3)/6)
C
      REAL*8 AuxR(0:MxEll),Sgn(0:MxLen),R(0:MxLen)
      REAL*8 R2(0:MxEll,0:MxEll,0:MxEll)
C
      INTEGER LBegin,LMNDex
C
      LBegin(L)=(L*(L+1)*(L+2))/6+1
      LMNDex(L,M,N)=LBegin(L+M+N)+N*(2*(L+M+N)-N+3)/2+M
C      LHGTF(L)=(L+1)*(L+2)*(L+3)/6
C
      Omega  = 5D-1*Zt
      Upq=2.473942945119315D1*Zt**(-2.5D0)
      TwoEll = 2*Ell
      o1     = 1D0
      o2     = -2D0*Omega
      DO L=0,TwoEll
         AuxR(L)=o1/DBLE(2*L+1)
         o1=o1*o2
      ENDDO
C
      DO J=TwoEll,0,-1
         DO I0=TwoEll-J,1,-1
            IF(I0-1.LE.0)THEN
               R2(I0,0,0)=0.0D0
            ELSE
               R2(I0,0,0)=DBLE(I0-1)*R2(I0-2,0,0)
            ENDIF
            DO I1=TwoEll-J-I0,1,-1
               IF(I1-1.LE.0)THEN
                  R2(I1,I0,0)=0.0D0
               ELSE
                  R2(I1,I0,0)=DBLE(I1-1)*R2(I1-2,I0,0)
               ENDIF
               DO I2=TwoEll-J-I0-I1,1,-1
                  IF(I2-1.LE.0)THEN
                     R2(I2,I1,I0)=0.0D0
                  ELSE
                     R2(I2,I1,I0)=DBLE(I2-1)*R2(I2-2,I1,I0)
                  ENDIF
               ENDDO
               IF(I1-1.LE.0)THEN
                  R2(0,I1,I0)=0.0D0
                  R2(I1,0,I0)=0.0D0
               ELSE
                  R2(0,I1,I0)=DBLE(I1-1)*R2(0,I1-2,I0)
                  R2(I1,0,I0)=DBLE(I1-1)*R2(I1-2,0,I0)
               ENDIF
            ENDDO
            IF(I0-1.LE.0)THEN
               R2(0,I0,0)=0.0D0
               R2(0,0,I0)=0.0D0
            ELSE
               R2(0,I0,0)=DBLE(I0-1)*R2(0,I0-2,0)
               R2(0,0,I0)=DBLE(I0-1)*R2(0,0,I0-2)
            ENDIF
         ENDDO
         R2(0,0,0)=AuxR(J)
      ENDDO
C
      DO L=0,TwoEll
         DO M=0,TwoEll-L
            DO N=0,TwoEll-L-M
               LMN=LMNdex(L,M,N)
               Sgn(LMN)=(-1D0)**(L+M+N)
               R(LMN)=R2(L,M,N)
            ENDDO
         ENDDO
      ENDDO
C
      DO I=1,Nq
         Est(I)=0D0
      ENDDO
C
      DO J=1,Len*Len
         IP=Pndex(J)
         IQ=Qndex(J)
         IPQ=PQndex(J)
         Phase=Sgn(IP)
         DO I=1,Nq
            VP=(I-1)*Len+IP
            VQ=(I-1)*Len+IQ
            Est(I)=Est(I)+Phase*Co(VP)*Co(VQ)*R(IPQ)
         ENDDO
      ENDDO
C
      DO I=1,Nq
         Est(I)=SQRT(Upq*Est(I))
      ENDDO
C
      RETURN
      END
C
      SUBROUTINE QSORTI (ORD,N,A)
C
C==============SORTS THE ARRAY A(I),I=1,2,...,N BY PUTTING THE
C   ASCENDING ORDER VECTOR IN ORD.  THAT IS ASCENDING ORDERED A
C   IS A(ORD(I)),I=1,2,...,N; DESCENDING ORDER A IS A(ORD(N-I+1)),
C   I=1,2,...,N .  THIS SORT RUNS IN TIME PROPORTIONAL TO N LOG N .
C
C
C     ACM QUICKSORT - ALGORITHM #402 - IMPLEMENTED IN FORTRAN 66 BY
C                                 WILLIAM H. VERITY, WHV@PSUVM.PSU.EDU
C                                 CENTER FOR ACADEMIC COMPUTING
C                                 THE PENNSYLVANIA STATE UNIVERSITY
C                                 UNIVERSITY PARK, PA.  16802
C
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION ORD(N),POPLST(2,20)
      INTEGER A(N),X,XX,Z,ZZ,Y
C
C     TO SORT DIFFERENT INPUT TYPES, CHANGE THE FOLLOWING
C     SPECIFICATION STATEMENTS; FOR EXAMPLE, FOR FORTRAN CHARACTER
C     USE THE FOLLOWING:  CHARACTER *(*) A(N)
C
c      REAL*8 A(N),X,Y,Z,XX,ZZ
C
      NDEEP=0
      U1=N
      L1=1
      DO 1  I=1,N
    1 ORD(I)=I
    2 IF (U1.LE.L1) RETURN
C
    3 L=L1
      U=U1
C
C PART
C
    4 P=L
      Q=U
C     FOR CHARACTER SORTS, THE FOLLOWING 3 STATEMENTS WOULD BECOME
C     X = ORD(P)
C     Z = ORD(Q)
C     IF (A(X) .LE. A(Z)) GO TO 2
C
C     WHERE "CLE" IS A LOGICAL FUNCTION WHICH RETURNS "TRUE" IF THE
C     FIRST ARGUMENT IS LESS THAN OR EQUAL TO THE SECOND, BASED ON "LEN"
C     CHARACTERS.
C
      X=A(ORD(P))
      Z=A(ORD(Q))
      IF (X.LE.Z) GO TO 5
      Y=X
      X=Z
      Z=Y
      YP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=YP
    5 IF (U-L.LE.1) GO TO 15
      XX=X
      IX=P
      ZZ=Z
      IZ=Q
C
C LEFT
C
    6 P=P+1
      IF (P.GE.Q) GO TO 7
      X=A(ORD(P))
      IF (X.GE.XX) GO TO 8
      GO TO 6
    7 P=Q-1
      GO TO 13
C
C RIGHT
C
    8 Q=Q-1
      IF (Q.LE.P) GO TO 9
      Z=A(ORD(Q))
      IF (Z.LE.ZZ) GO TO 10
      GO TO 8
    9 Q=P
      P=P-1
      Z=X
      X=A(ORD(P))
C
C DIST
C
   10 IF (X.LE.Z) GO TO 11
      Y=X
      X=Z
      Z=Y
      IP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=IP
   11 IF (X.LE.XX) GO TO 12
      XX=X
      IX=P
   12 IF (Z.GE.ZZ) GO TO 6
      ZZ=Z
      IZ=Q
      GO TO 6
C
C OUT
C
   13 CONTINUE
      IF (.NOT.(P.NE.IX.AND.X.NE.XX)) GO TO 14
      IP=ORD(P)
      ORD(P)=ORD(IX)
      ORD(IX)=IP
   14 CONTINUE
      IF (.NOT.(Q.NE.IZ.AND.Z.NE.ZZ)) GO TO 15
      IQ=ORD(Q)
      ORD(Q)=ORD(IZ)
      ORD(IZ)=IQ
   15 CONTINUE
      IF (U-Q.LE.P-L) GO TO 16
      L1=L
      U1=P-1
      L=Q+1
      GO TO 17
   16 U1=U
      L1=Q+1
      U=P-1
   17 CONTINUE
      IF (U1.LE.L1) GO TO 18
C
C START RECURSIVE CALL
C
      NDEEP=NDEEP+1
      POPLST(1,NDEEP)=U
      POPLST(2,NDEEP)=L
      GO TO 3
   18 IF (U.GT.L) GO TO 4
C
C POP BACK UP IN THE RECURSION LIST
C
      IF (NDEEP.EQ.0) GO TO 2
      U=POPLST(1,NDEEP)
      L=POPLST(2,NDEEP)
      NDEEP=NDEEP-1
      GO TO 18
C
C END SORT
C END QSORT
C
      END
