      SUBROUTINE VRRg(N,LBra,LKet,NUE0,NUF0,WR,WZ,W,Z)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER Fr,To
      REAL*8 W(N,NUE0,NUF0,2)
      REAL*8 Z(N,*)
      REAL*8 WR(12,*)
      REAL*8 WZ( 5,*)
      SAVE Two,Three,Four,Five
      DATA Two,Three,Four,Five / 2.0D0,3.0D0,4.0D0,5.0D0 /


      MxLBra = LBra
      MxLKet = LKet
      LTot   = LBra + LKet 

      To = MOD(LTot+1,2)+1
      DO 20 Knt = 0,LTot
         Fr = To
         To = 3 - To

      CALL DCOPY(N,Z(1,LTot-Knt+1),1,W(1,1,1,To),1)
      IF(MxLBra.le.0.or.Knt.eq.0) GO TO 10

      DO I = 1,N
         S000 = W(I,1,1,To)
         T000 = W(I,1,1,Fr)

         W(I,2,1,To) = WR(1,I)*S000 + WR( 7,I)*T000
         W(I,3,1,To) = WR(3,I)*S000 + WR( 9,I)*T000
         W(I,4,1,To) = WR(5,I)*S000 + WR(11,I)*T000

      END DO
      IF(MxLBra.le.1.or.Knt.eq.1) GO TO 10

      DO I = 1,N
         PAx  = WR( 1,I)
         PAy  = WR( 3,I)
         PAz  = WR( 5,I)
         SgP  = WZ( 2,I)
         SP2  = WZ( 4,I)
         WPx  = WR( 7,I)
         WPy  = WR( 9,I)
         WPz  = WR(11,I)
         U000 = SgP*(W(I,1,1,To) - W(I,1,1,Fr)*SP2)
         S100 = W(I,2,1,To)
         S010 = W(I,3,1,To)
         S001 = W(I,4,1,To)
         T100 = W(I,2,1,Fr)
         T010 = W(I,3,1,Fr)
         T001 = W(I,4,1,Fr)
         W(I, 5,1,To) = PAx*S100 + WPx*T100 + U000
         W(I, 6,1,To) = PAy*S100 + WPy*T100
         W(I, 7,1,To) = PAy*S010 + WPy*T010 + U000
         W(I, 8,1,To) = PAx*S001 + WPx*T001
         W(I, 9,1,To) = PAz*S010 + WPz*T010
         W(I,10,1,To) = PAz*S001 + WPz*T001 + U000
      END DO
      IF(MxLBra.le.2.or.Knt.eq.2) GO TO 10

      DO I = 1,N
         PAx  = WR( 1,I)
         PAy  = WR( 3,I)
         PAz  = WR( 5,I)
         SgP  = WZ( 2,I)
         SP2  = WZ( 4,I)
         WPx  = WR( 7,I)
         WPy  = WR( 9,I)
         WPz  = WR(11,I)
         U100 = SgP*(W(I, 2,1,To) - W(I,2,1,Fr)*SP2)
         U010 = SgP*(W(I, 3,1,To) - W(I,3,1,Fr)*SP2)
         U001 = SgP*(W(I, 4,1,To) - W(I,4,1,Fr)*SP2)
         S200 = W(I, 5,1,To)
         S020 = W(I, 7,1,To)
         S011 = W(I, 9,1,To)
         S002 = W(I,10,1,To)
         T200 = W(I, 5,1,Fr)
         T020 = W(I, 7,1,Fr)
         T011 = W(I, 9,1,Fr)
         T002 = W(I,10,1,Fr)
         W(I,11,1,To) = PAx*S200 + WPx*T200 + U100*Two
         W(I,12,1,To) = PAy*S200 + WPy*T200
         W(I,13,1,To) = PAx*S020 + WPx*T020
         W(I,14,1,To) = PAy*S020 + WPy*T020 + U010*Two
         W(I,15,1,To) = PAz*S200 + WPz*T200
         W(I,16,1,To) = PAx*S011 + WPx*T011
         W(I,17,1,To) = PAz*S020 + WPz*T020
         W(I,18,1,To) = PAx*S002 + WPx*T002
         W(I,19,1,To) = PAy*S002 + WPy*T002
         W(I,20,1,To) = PAz*S002 + WPz*T002 + U001*Two
      END DO
      IF(MxLBra.le.3.or.Knt.eq.3) GO TO 10

      DO I = 1,N
         PAx  = WR( 1,I)
         PAy  = WR( 3,I)
         PAz  = WR( 5,I)
         SgP  = WZ( 2,I)
         SP2  = WZ( 4,I)
         WPx  = WR( 7,I)
         WPy  = WR( 9,I)
         WPz  = WR(11,I)
         U200 = SgP*(W(I, 5,1,To) - W(I, 5,1,Fr)*SP2)
         U020 = SgP*(W(I, 7,1,To) - W(I, 7,1,Fr)*SP2)
         U002 = SgP*(W(I,10,1,To) - W(I,10,1,Fr)*SP2)
         S300 = W(I,11,1,To)
         S120 = W(I,13,1,To)
         S030 = W(I,14,1,To)
         S201 = W(I,15,1,To)
         S012 = W(I,19,1,To)
         S003 = W(I,20,1,To)
         T300 = W(I,11,1,Fr)
         T120 = W(I,13,1,Fr)
         T030 = W(I,14,1,Fr)
         T201 = W(I,15,1,Fr)
         T012 = W(I,19,1,Fr)
         T003 = W(I,20,1,Fr)
         W(I,21,1,To) = PAx*S300 + WPx*T300 + U200*Three
         W(I,22,1,To) = PAy*S300 + WPy*T300
         W(I,23,1,To) = PAx*S120 + WPx*T120 + U020
         W(I,24,1,To) = PAx*S030 + WPx*T030
         W(I,25,1,To) = PAy*S030 + WPy*T030 + U020*Three
         W(I,26,1,To) = PAz*S300 + WPz*T300
         W(I,27,1,To) = PAy*S201 + WPy*T201
         W(I,28,1,To) = PAz*S120 + WPz*T120
         W(I,29,1,To) = PAz*S030 + WPz*T030
         W(I,30,1,To) = PAz*S201 + WPz*T201 + U200
         W(I,31,1,To) = PAx*S012 + WPx*T012
         W(I,32,1,To) = PAy*S012 + WPy*T012 + U002
         W(I,33,1,To) = PAx*S003 + WPx*T003
         W(I,34,1,To) = PAy*S003 + WPy*T003
         W(I,35,1,To) = PAz*S003 + WPz*T003 + U002*Three
      END DO
      IF(MxLBra.le.4.or.Knt.eq.4) GO TO 10

      DO I = 1,N
         PAx  = WR( 1,I)
         PAy  = WR( 3,I)
         PAz  = WR( 5,I)
         SgP  = WZ( 2,I)
         SP2  = WZ( 4,I)
         WPx  = WR( 7,I)
         WPy  = WR( 9,I)
         WPz  = WR(11,I)
         U300 = SgP*(W(I,11,1,To) - W(I,11,1,Fr)*SP2)
         U120 = SgP*(W(I,13,1,To) - W(I,13,1,Fr)*SP2)
         U030 = SgP*(W(I,14,1,To) - W(I,14,1,Fr)*SP2)
         U003 = SgP*(W(I,20,1,To) - W(I,20,1,Fr)*SP2)
         S400 = W(I,21,1,To)
         S310 = W(I,22,1,To)
         S220 = W(I,23,1,To)
         S130 = W(I,24,1,To)
         S040 = W(I,25,1,To)
         S301 = W(I,26,1,To)
         S121 = W(I,28,1,To)
         S031 = W(I,29,1,To)
         S202 = W(I,30,1,To)
         S103 = W(I,33,1,To)
         S013 = W(I,34,1,To)
         S004 = W(I,35,1,To)
         T400 = W(I,21,1,Fr)
         T310 = W(I,22,1,Fr)
         T220 = W(I,23,1,Fr)
         T130 = W(I,24,1,Fr)
         T040 = W(I,25,1,Fr)
         T301 = W(I,26,1,Fr)
         T121 = W(I,28,1,Fr)
         T031 = W(I,29,1,Fr)
         T202 = W(I,30,1,Fr)
         T103 = W(I,33,1,Fr)
         T013 = W(I,34,1,Fr)
         T004 = W(I,35,1,Fr)
         W(I,36,1,To) = PAx*S400 + WPx*T400 + U300*Four
         W(I,37,1,To) = PAy*S400 + WPy*T400
         W(I,38,1,To) = PAy*S310 + WPy*t310 + U300
         W(I,39,1,To) = PAx*S130 + WPx*T130 + U030
         W(I,40,1,To) = PAx*S040 + WPx*T040
         W(I,41,1,To) = PAy*S040 + WPy*T040 + U030*Four
         W(I,42,1,To) = PAz*S400 + WPz*T400
         W(I,43,1,To) = PAy*S301 + WPy*T301
         W(I,44,1,To) = PAz*S220 + WPz*T220
         W(I,45,1,To) = PAz*S130 + WPz*T130
         W(I,46,1,To) = PAz*S040 + WPz*T040
         W(I,47,1,To) = PAz*S301 + WPz*T301 + U300
         W(I,48,1,To) = PAy*S202 + WPy*T202
         W(I,49,1,To) = PAz*S121 + WPz*T121 + U120
         W(I,50,1,To) = PAz*S031 + WPz*T031 + U030
         W(I,51,1,To) = PAx*S103 + WPx*T103 + U003
         W(I,52,1,To) = PAx*S013 + WPx*T013
         W(I,53,1,To) = PAy*S013 + WPy*T013 + U003
         W(I,54,1,To) = PAx*S004 + WPx*T004
         W(I,55,1,To) = PAy*S004 + WPy*T004
         W(I,56,1,To) = PAz*S004 + WPz*T004 + U003*Four
      END DO
      IF(MxLBra.le.5.or.Knt.eq.5) GO TO 10

      DO I = 1,N
         PAx  = WR( 1,I)
         PAy  = WR( 3,I)
         PAz  = WR( 5,I)
         SgP  = WZ( 2,I)
         SP2  = WZ( 4,I)
         WPx  = WR( 7,I)
         WPy  = WR( 9,I)
         WPz  = WR(11,I)
         U400 = SgP*(W(I,21,1,To) - W(I,21,1,Fr)*SP2)
         U220 = SgP*(W(I,23,1,To) - W(I,23,1,Fr)*SP2)
         U130 = SgP*(W(I,24,1,To) - W(I,24,1,Fr)*SP2)
         U040 = SgP*(W(I,25,1,To) - W(I,25,1,Fr)*SP2)
         U301 = SgP*(W(I,26,1,To) - W(I,26,1,Fr)*SP2)
         U013 = SgP*(W(I,34,1,To) - W(I,34,1,Fr)*SP2)
         U004 = SgP*(W(I,35,1,To) - W(I,35,1,Fr)*SP2)
         S500 = W(I,36,1,To)
         S410 = W(I,37,1,To)
         S230 = W(I,39,1,To)
         S140 = W(I,40,1,To)
         S050 = W(I,41,1,To)
         S401 = W(I,42,1,To)
         S311 = W(I,43,1,To)
         S221 = W(I,44,1,To)
         S131 = W(I,45,1,To)
         S041 = W(I,46,1,To)
         S302 = W(I,47,1,To)
         S113 = W(I,52,1,To)
         S023 = W(I,53,1,To)
         S104 = W(I,54,1,To)
         S014 = W(I,55,1,To)
         S005 = W(I,56,1,To)
         T500 = W(I,36,1,Fr)
         T410 = W(I,37,1,Fr)
         T230 = W(I,39,1,Fr)
         T140 = W(I,40,1,Fr)
         T050 = W(I,41,1,Fr)
         T401 = W(I,42,1,Fr)
         T311 = W(I,43,1,Fr)
         T221 = W(I,44,1,Fr)
         T131 = W(I,45,1,Fr)
         T041 = W(I,46,1,Fr)
         T302 = W(I,47,1,Fr)
         T113 = W(I,52,1,Fr)
         T023 = W(I,53,1,Fr)
         T104 = W(I,54,1,Fr)
         T014 = W(I,55,1,Fr)
         T005 = W(I,56,1,Fr)
         W(I,57,1,To) = PAx*S500 + WPx*T500 + U400*Five
         W(I,58,1,To) = PAy*S500 + WPy*T500
         W(I,59,1,To) = PAy*S410 + WPy*T410 + U400
         W(I,60,1,To) = PAx*S230 + WPx*T230 + U130*Two
         W(I,61,1,To) = PAx*S140 + WPx*T140 + U040    
         W(I,62,1,To) = PAx*S050 + WPx*T050
         W(I,63,1,To) = PAy*S050 + WPy*T050 + U040*Five
         W(I,64,1,To) = PAz*S500 + WPz*T500
         W(I,65,1,To) = PAy*S401 + WPy*T401
         W(I,66,1,To) = PAy*S311 + WPy*T311 + U301
         W(I,67,1,To) = PAz*S230 + WPz*T230
         W(I,68,1,To) = PAz*S140 + WPz*T140
         W(I,69,1,To) = PAz*S050 + WPz*T050
         W(I,70,1,To) = PAz*S401 + WPz*T401 + U400
         W(I,71,1,To) = PAy*S302 + WPy*T302
         W(I,72,1,To) = PAz*S221 + WPz*T221 + U220
         W(I,73,1,To) = PAz*S131 + WPz*T131 + U130
         W(I,74,1,To) = PAz*S041 + WPz*T041 + U040
         W(I,75,1,To) = PAz*S302 + WPz*T302 + U301*Two
         W(I,76,1,To) = PAx*S113 + WPx*T113 + U013    
         W(I,77,1,To) = PAx*S023 + WPx*T023
         W(I,78,1,To) = PAy*S023 + WPy*T023 + U013*Two
         W(I,79,1,To) = PAx*S104 + WPx*T104 + U004    
         W(I,80,1,To) = PAx*S014 + WPx*T014
         W(I,81,1,To) = PAy*S014 + WPy*T014 + U004
         W(I,82,1,To) = PAx*S005 + WPx*T005
         W(I,83,1,To) = PAy*S005 + WPy*T005
         W(I,84,1,To) = PAz*S005 + WPz*T005 + U004*Five
      END DO
      IF(MxLBra.le.6.or.Knt.eq.6) GO TO 10
      CALL MondoHalt('MxLBra > 6 in VRRG')

10    CONTINUE
      IF(MxLKet.le.LTot-Knt+0) GO TO 20

      NumL = LFuncC(0,MIN(MxLBra,Knt-1))
      DO L = 1,NumL
      DO I = 1,N
         QBx = WR(2,I)
         QBy = WR(4,I)
         QBz = WR(6,I)
         S000 = W(I,L,1,To)
         T000 = W(I,L,1,Fr)
         W(I,L,2,To) = QBx*S000 + WR( 8,I)*T000
         W(I,L,3,To) = QBy*S000 + WR(10,I)*T000
         W(I,L,4,To) = QBz*S000 + WR(12,I)*T000
      END DO
      END DO
      DO L = 2,NumL
         CALL KonK2L(L,Lx,Ly,Lz)
         IF(Lx.gt.0) THEN
            CALL KonL2K(LP,Lx-1,Ly,Lz)
            DO I = 1,N
               Y = WZ(5,I) * Lx
               W(I,L,2,To) = W(I,L,2,To)+W(I,LP,1,Fr)*Y
            END DO
         END IF
         IF(Ly.gt.0) THEN
            CALL KonL2K(LP,Lx,Ly-1,Lz)
            DO I = 1,N
               Y = WZ(5,I) * Ly
               W(I,L,3,To) = W(I,L,3,To)+W(I,LP,1,Fr)*Y
            END DO
         END IF
         IF(Lz.gt.0) THEN
            CALL KonL2K(LP,Lx,Ly,Lz-1)
            DO I = 1,N
               Y = WZ(5,I) * Lz
               W(I,L,4,To) = W(I,L,4,To)+W(I,LP,1,Fr)*Y
            END DO
         END IF
      END DO

      IF(MxLKet.le.LTot-Knt+1) GO TO 20

      NumL = LFuncC(0,MIN(MxLBra,Knt-2))
      DO L = 1,NumL
      DO I = 1,N
         QBx = WR(2,I)
         QBy = WR(4,I)
         QBz = WR(6,I)
         SgQ = WZ(1,I)
         SQ2 = WZ(3,I)
         WQx = WR( 8,I)
         WQy = WR(10,I)
         WQz = WR(12,I)
         U000 = SgQ*(W(I,L,1,To) - W(I,L,1,Fr)*SQ2)
         S100 = W(I,L,2,To)
         S010 = W(I,L,3,To)
         S001 = W(I,L,4,To)
         T100 = W(I,L,2,Fr)
         T010 = W(I,L,3,Fr)
         T001 = W(I,L,4,Fr)
         W(I,L, 5,To) = QBx*S100 + WQx*T100 + U000
         W(I,L, 6,To) = QBy*S100 + WQy*T100
         W(I,L, 7,To) = QBy*S010 + WQy*T010 + U000
         W(I,L, 8,To) = QBx*S001 + WQx*T001
         W(I,L, 9,To) = QBz*S010 + WQz*T010
         W(I,L,10,To) = QBz*S001 + WQz*T001 + U000
      END DO
      END DO
      DO L = 2,NumL
         CALL KonK2L(L,Lx,Ly,Lz)
         IF(Lx.gt.0) THEN
            CALL KonL2K(LP,Lx-1,Ly,Lz)
            DO I = 1,N
               Y = WZ(5,I) * Lx
               W(I,L,5,To) = W(I,L,5,To) + W(I,LP,2,Fr)*Y
               W(I,L,8,To) = W(I,L,8,To) + W(I,LP,4,Fr)*Y
            END DO
         END IF
         IF(Ly.gt.0) THEN
            CALL KonL2K(LP,Lx,Ly-1,Lz)
            DO I = 1,N
               Y = WZ(5,I) * Ly
               W(I,L,6,To) = W(I,L,6,To) + W(I,LP,2,Fr)*Y
               W(I,L,7,To) = W(I,L,7,To) + W(I,LP,3,Fr)*Y
            END DO
         END IF
         IF(Lz.gt.0) THEN
            CALL KonL2K(LP,Lx,Ly,Lz-1)
            DO I = 1,N
               Y = WZ(5,I) * Lz
               W(I,L, 9,To) = W(I,L, 9,To)+W(I,LP,3,Fr)*Y
               W(I,L,10,To) = W(I,L,10,To)+W(I,LP,4,Fr)*Y
            END DO
         END IF
      END DO
      IF(MxLKet.le.LTot-Knt+2) GO TO 20

      NumL = LFuncC(0,MIN(MxLBra,Knt-3))
      DO L = 1,NumL
      DO I = 1,N
         QBx = WR(2,I)
         QBy = WR(4,I)
         QBz = WR(6,I)
         SgQ = WZ(1,I)
         SQ2 = WZ(3,I)
         WQx = WR( 8,I)
         WQy = WR(10,I)
         WQz = WR(12,I)
         U100 = SgQ*(W(I,L, 2,To) - W(I,L,2,Fr)*SQ2)
         U010 = SgQ*(W(I,L, 3,To) - W(I,L,3,Fr)*SQ2)
         U001 = SgQ*(W(I,L, 4,To) - W(I,L,4,Fr)*SQ2)
         S200 = W(I,L, 5,To)
         S020 = W(I,L, 7,To)
         S011 = W(I,L, 9,To)
         S002 = W(I,L,10,To)
         T200 = W(I,L, 5,Fr)
         T020 = W(I,L, 7,Fr)
         T011 = W(I,L, 9,Fr)
         T002 = W(I,L,10,Fr)
         W(I,L,11,To) = QBx*S200 + WQx*T200 + U100*Two
         W(I,L,12,To) = QBy*S200 + WQy*T200
         W(I,L,13,To) = QBx*S020 + WQx*T020
         W(I,L,14,To) = QBy*S020 + WQy*T020 + U010*Two
         W(I,L,15,To) = QBz*S200 + WQz*T200
         W(I,L,16,To) = QBx*S011 + WQx*T011
         W(I,L,17,To) = QBz*S020 + WQz*T020
         W(I,L,18,To) = QBx*S002 + WQx*T002
         W(I,L,19,To) = QBy*S002 + WQy*T002
         W(I,L,20,To) = QBz*S002 + WQz*T002 + U001*Two
      END DO
      END DO
      DO L = 2,NumL
         CALL KonK2L(L,Lx,Ly,Lz)
         IF(Lx.gt.0) THEN
            CALL KonL2K(LP,Lx-1,Ly,Lz)
            DO I = 1,N
               Y = WZ(5,I) * Lx
               W(I,L,11,To) = W(I,L,11,To)+W(I,LP, 5,Fr)*Y
               W(I,L,13,To) = W(I,L,13,To)+W(I,LP, 7,Fr)*Y
               W(I,L,16,To) = W(I,L,16,To)+W(I,LP, 9,Fr)*Y
               W(I,L,18,To) = W(I,L,18,To)+W(I,LP,10,Fr)*Y
            END DO
         END IF
         IF(Ly.gt.0) THEN
            CALL KonL2K(LP,Lx,Ly-1,Lz)
            DO I = 1,N
               Y = WZ(5,I) * Ly
               W(I,L,12,To) = W(I,L,12,To)+W(I,LP, 5,Fr)*Y
               W(I,L,14,To) = W(I,L,14,To)+W(I,LP, 7,Fr)*Y
               W(I,L,19,To) = W(I,L,19,To)+W(I,LP,10,Fr)*Y
            END DO
         END IF
         IF(Lz.gt.0) THEN
            CALL KonL2K(LP,Lx,Ly,Lz-1)
            DO I = 1,N
               Y = WZ(5,I) * Lz
               W(I,L,15,To) = W(I,L,15,To)+W(I,LP, 5,Fr)*Y
               W(I,L,17,To) = W(I,L,17,To)+W(I,LP, 7,Fr)*Y
               W(I,L,20,To) = W(I,L,20,To)+W(I,LP,10,Fr)*Y
            END DO
         END IF
      END DO
      IF(MxLKet.le.LTot-Knt+3) GO TO 20

      NumL = LFuncC(0,MIN(MxLBra,Knt-4))
      DO L = 1,NumL
      DO I = 1,N
         QBx = WR(2,I)
         QBy = WR(4,I)
         QBz = WR(6,I)
         SgQ = WZ(1,I)
         SQ2 = WZ(3,I)
         WQx = WR( 8,I)
         WQy = WR(10,I)
         WQz = WR(12,I)
         U200 = SgQ*(W(I,L, 5,To) - W(I,L, 5,Fr)*SQ2)
         U020 = SgQ*(W(I,L, 7,To) - W(I,L, 7,Fr)*SQ2)
         U002 = SgQ*(W(I,L,10,To) - W(I,L,10,Fr)*SQ2)
         S300 = W(I,L,11,To)
         S120 = W(I,L,13,To)
         S030 = W(I,L,14,To)
         S201 = W(I,L,15,To)
         S012 = W(I,L,19,To)
         S003 = W(I,L,20,To)
         T300 = W(I,L,11,Fr)
         T120 = W(I,L,13,Fr)
         T030 = W(I,L,14,Fr)
         T201 = W(I,L,15,Fr)
         T012 = W(I,L,19,Fr)
         T003 = W(I,L,20,Fr)
         W(I,L,21,To) = QBx*S300 + WQx*T300 + U200*Three
         W(I,L,22,To) = QBy*S300 + WQy*T300
         W(I,L,23,To) = QBx*S120 + WQx*T120 + U020
         W(I,L,24,To) = QBx*S030 + WQx*T030
         W(I,L,25,To) = QBy*S030 + WQy*T030 + U020*Three
         W(I,L,26,To) = QBz*S300 + WQz*T300
         W(I,L,27,To) = QBy*S201 + WQy*T201
         W(I,L,28,To) = QBz*S120 + WQz*T120
         W(I,L,29,To) = QBz*S030 + WQz*T030
         W(I,L,30,To) = QBz*S201 + WQz*T201 + U200
         W(I,L,31,To) = QBx*S012 + WQx*T012
         W(I,L,32,To) = QBy*S012 + WQy*T012 + U002
         W(I,L,33,To) = QBx*S003 + WQx*T003
         W(I,L,34,To) = QBy*S003 + WQy*T003
         W(I,L,35,To) = QBz*S003 + WQz*T003 + U002*Three
      END DO
      END DO
      DO L = 2,NumL
         CALL KonK2L(L,Lx,Ly,Lz)
         IF(Lx.gt.0) THEN
            CALL KonL2K(LP,Lx-1,Ly,Lz)
            DO I = 1,N
               Y = WZ(5,I) * Lx
               W(I,L,21,To) = W(I,L,21,To)+W(I,LP,11,Fr)*Y
               W(I,L,23,To) = W(I,L,23,To)+W(I,LP,13,Fr)*Y
               W(I,L,24,To) = W(I,L,24,To)+W(I,LP,14,Fr)*Y
               W(I,L,31,To) = W(I,L,31,To)+W(I,LP,19,Fr)*Y
               W(I,L,33,To) = W(I,L,33,To)+W(I,LP,20,Fr)*Y
            END DO
         END IF
         IF(Ly.gt.0) THEN
            CALL KonL2K(LP,Lx,Ly-1,Lz)
            DO I = 1,N
               Y = WZ(5,I) * Ly
               W(I,L,22,To) = W(I,L,22,To)+W(I,LP,11,Fr)*Y
               W(I,L,25,To) = W(I,L,25,To)+W(I,LP,14,Fr)*Y
               W(I,L,27,To) = W(I,L,27,To)+W(I,LP,15,Fr)*Y
               W(I,L,32,To) = W(I,L,32,To)+W(I,LP,19,Fr)*Y
               W(I,L,34,To) = W(I,L,34,To)+W(I,LP,20,Fr)*Y
            END DO
         END IF
         IF(Lz.gt.0) THEN
            CALL KonL2K(LP,Lx,Ly,Lz-1)
            DO I = 1,N
               Y = WZ(5,I) * Lz
               W(I,L,26,To) = W(I,L,26,To)+W(I,LP,11,Fr)*Y
               W(I,L,28,To) = W(I,L,28,To)+W(I,LP,13,Fr)*Y
               W(I,L,29,To) = W(I,L,29,To)+W(I,LP,14,Fr)*Y
               W(I,L,30,To) = W(I,L,30,To)+W(I,LP,15,Fr)*Y
               W(I,L,35,To) = W(I,L,35,To)+W(I,LP,20,Fr)*Y
            END DO
         END IF
      END DO
      IF(MxLKet.le.LTot-Knt+4) GO TO 20

      NumL = LFuncC(0,MIN(MxLBra,Knt-5))
      DO L = 1,NumL
      DO I = 1,N
         QBx = WR(2,I)
         QBy = WR(4,I)
         QBz = WR(6,I)
         SgQ = WZ(1,I)
         SQ2 = WZ(3,I)
         WQx = WR( 8,I)
         WQy = WR(10,I)
         WQz = WR(12,I)
         U300 = SgQ*(W(I,L,11,To) - W(I,L,11,Fr)*SQ2)
         U120 = SgQ*(W(I,L,13,To) - W(I,L,13,Fr)*SQ2)
         U030 = SgQ*(W(I,L,14,To) - W(I,L,14,Fr)*SQ2)
         U003 = SgQ*(W(I,L,20,To) - W(I,L,20,Fr)*SQ2)
         S400 = W(I,L,21,To)
         S310 = W(I,L,22,To)
         S220 = W(I,L,23,To)
         S130 = W(I,L,24,To)
         S040 = W(I,L,25,To)
         S301 = W(I,L,26,To)
         S121 = W(I,L,28,To)
         S031 = W(I,L,29,To)
         S202 = W(I,L,30,To)
         S103 = W(I,L,33,To)
         S013 = W(I,L,34,To)
         S004 = W(I,L,35,To)
         T400 = W(I,L,21,Fr)
         T310 = W(I,L,22,Fr)
         T220 = W(I,L,23,Fr)
         T130 = W(I,L,24,Fr)
         T040 = W(I,L,25,Fr)
         T301 = W(I,L,26,Fr)
         T121 = W(I,L,28,Fr)
         T031 = W(I,L,29,Fr)
         T202 = W(I,L,30,Fr)
         T103 = W(I,L,33,Fr)
         T013 = W(I,L,34,Fr)
         T004 = W(I,L,35,Fr)
         W(I,L,36,To) = QBx*S400 + WQx*T400 + U300*Four
         W(I,L,37,To) = QBy*S400 + WQy*T400
         W(I,L,38,To) = QBy*S310 + WQy*T310 + U300
         W(I,L,39,To) = QBx*S130 + WQx*T130 + U030
         W(I,L,40,To) = QBx*S040 + WQx*T040
         W(I,L,41,To) = QBy*S040 + WQy*T040 + U030*Four
         W(I,L,42,To) = QBz*S400 + WQz*T400
         W(I,L,43,To) = QBy*S301 + WQy*T301
         W(I,L,44,To) = QBz*S220 + WQz*T220
         W(I,L,45,To) = QBz*S130 + WQz*T130
         W(I,L,46,To) = QBz*S040 + WQz*T040
         W(I,L,47,To) = QBz*S301 + WQz*T301 + U300
         W(I,L,48,To) = QBy*S202 + WQy*T202
         W(I,L,49,To) = QBz*S121 + WQz*T121 + U120
         W(I,L,50,To) = QBz*S031 + WQz*T031 + U030
         W(I,L,51,To) = QBx*S103 + WQx*T103 + U003
         W(I,L,52,To) = QBx*S013 + WQx*T013
         W(I,L,53,To) = QBy*S013 + WQy*T013 + U003
         W(I,L,54,To) = QBx*S004 + WQx*T004
         W(I,L,55,To) = QBy*S004 + WQy*T004
         W(I,L,56,To) = QBz*S004 + WQz*T004 + U003*Four
      END DO
      END DO
      DO L = 2,NumL
        CALL KonK2L(L,Lx,Ly,Lz)
         IF(Lx.gt.0) THEN
            CALL KonL2K(LP,Lx-1,Ly,Lz)
            DO I = 1,N
               Y = WZ(5,I) * Lx
               W(I,L,36,To) = W(I,L,36,To)+W(I,LP,21,Fr)*Y
               W(I,L,39,To) = W(I,L,39,To)+W(I,LP,24,Fr)*Y
               W(I,L,40,To) = W(I,L,40,To)+W(I,LP,25,Fr)*Y
               W(I,L,51,To) = W(I,L,51,To)+W(I,LP,33,Fr)*Y
               W(I,L,52,To) = W(I,L,52,To)+W(I,LP,34,Fr)*Y
               W(I,L,54,To) = W(I,L,54,To)+W(I,LP,35,Fr)*Y
            END DO
         END IF
         IF(Ly.gt.0) THEN
            CALL KonL2K(LP,Lx,Ly-1,Lz)
            DO I = 1,N
               Y = WZ(5,I) * Ly
               W(I,L,37,To) = W(I,L,37,To)+W(I,LP,21,Fr)*Y
               W(I,L,38,To) = W(I,L,38,To)+W(I,LP,22,Fr)*Y
               W(I,L,41,To) = W(I,L,41,To)+W(I,LP,25,Fr)*Y
               W(I,L,43,To) = W(I,L,43,To)+W(I,LP,26,Fr)*Y
               W(I,L,48,To) = W(I,L,48,To)+W(I,LP,30,Fr)*Y
               W(I,L,53,To) = W(I,L,53,To)+W(I,LP,34,Fr)*Y
               W(I,L,55,To) = W(I,L,55,To)+W(I,LP,35,Fr)*Y
            END DO
         END IF
         IF(Lz.gt.0) THEN
            CALL KonL2K(LP,Lx,Ly,Lz-1)
            DO I = 1,N
               Y = WZ(5,I) * Lz
               W(I,L,42,To) = W(I,L,42,To)+W(I,LP,21,Fr)*Y
               W(I,L,44,To) = W(I,L,44,To)+W(I,LP,23,Fr)*Y
               W(I,L,45,To) = W(I,L,45,To)+W(I,LP,24,Fr)*Y
               W(I,L,46,To) = W(I,L,46,To)+W(I,LP,25,Fr)*Y
               W(I,L,47,To) = W(I,L,47,To)+W(I,LP,26,Fr)*Y
               W(I,L,49,To) = W(I,L,49,To)+W(I,LP,28,Fr)*Y
               W(I,L,50,To) = W(I,L,50,To)+W(I,LP,29,Fr)*Y
               W(I,L,56,To) = W(I,L,56,To)+W(I,LP,35,Fr)*Y
            END DO
         END IF
      END DO
      IF(MxLKet.le.LTot-Knt+5) GO TO 20

      NumL = LFuncC(0,MIN(MxLBra,Knt-6))
      DO L = 1,NumL
      DO I = 1,N
         QBx = WR(2,I)
         QBy = WR(4,I)
         QBz = WR(6,I)
         SgQ = WZ(1,I)
         SQ2 = WZ(3,I)
         WQx = WR( 8,I)
         WQy = WR(10,I)
         WQz = WR(12,I)
         U400 = SgQ*(W(I,L,21,To) - W(I,L,21,Fr)*SQ2)
         U220 = SgQ*(W(I,L,23,To) - W(I,L,23,Fr)*SQ2)
         U130 = SgQ*(W(I,L,24,To) - W(I,L,24,Fr)*SQ2)
         U040 = SgQ*(W(I,L,25,To) - W(I,L,25,Fr)*SQ2)
         U301 = SgQ*(W(I,L,26,To) - W(I,L,26,Fr)*SQ2)
         U013 = SgQ*(W(I,L,34,To) - W(I,L,34,Fr)*SQ2)
         U004 = SgQ*(W(I,L,35,To) - W(I,L,35,Fr)*SQ2)
         S500 = W(I,L,36,To)
         S410 = W(I,L,37,To)
         S230 = W(I,L,39,To)
         S140 = W(I,L,40,To)
         S050 = W(I,L,41,To)
         S401 = W(I,L,42,To)
         S311 = W(I,L,43,To)
         S221 = W(I,L,44,To)
         S131 = W(I,L,45,To)
         S041 = W(I,L,46,To)
         S302 = W(I,L,47,To)
         S113 = W(I,L,52,To)
         S023 = W(I,L,53,To)
         S104 = W(I,L,54,To)
         S014 = W(I,L,55,To)
         S005 = W(I,L,56,To)
         T500 = W(I,L,36,Fr)
         T410 = W(I,L,37,Fr)
         T230 = W(I,L,39,Fr)
         T140 = W(I,L,40,Fr)
         T050 = W(I,L,41,Fr)
         T401 = W(I,L,42,Fr)
         T311 = W(I,L,43,Fr)
         T221 = W(I,L,44,Fr)
         T131 = W(I,L,45,Fr)
         T041 = W(I,L,46,Fr)
         T302 = W(I,L,47,Fr)
         T113 = W(I,L,52,Fr)
         T023 = W(I,L,53,Fr)
         T104 = W(I,L,54,Fr)
         T014 = W(I,L,55,Fr)
         T005 = W(I,L,56,Fr)
         W(I,L,57,To) = QBx*S500 + WQx*T500 + U400*Five
         W(I,L,58,To) = QBy*S500 + WQy*T500             
         W(I,L,59,To) = QBy*S410 + WQy*T410 + U400     
         W(I,L,60,To) = QBx*S230 + WQx*T230 + U130*Two
         W(I,L,61,To) = QBx*S140 + WQx*T140 + U040    
         W(I,L,62,To) = QBx*S050 + WQx*T050            
         W(I,L,63,To) = QBy*S050 + WQy*T050 + U040*Five
         W(I,L,64,To) = QBz*S500 + WQz*T500             
         W(I,L,65,To) = QBy*S401 + WQy*T401             
         W(I,L,66,To) = QBy*S311 + WQy*T311 + U301
         W(I,L,67,To) = QBz*S230 + WQz*T230        
         W(I,L,68,To) = QBz*S140 + WQz*T140        
         W(I,L,69,To) = QBz*S050 + WQz*T050
         W(I,L,70,To) = QBz*S401 + WQz*T401 + U400
         W(I,L,71,To) = QBy*S302 + WQy*T302        
         W(I,L,72,To) = QBz*S221 + WQz*T221 + U220
         W(I,L,73,To) = QBz*S131 + WQz*T131 + U130
         W(I,L,74,To) = QBz*S041 + WQz*T041 + U040
         W(I,L,75,To) = QBz*S302 + WQz*T302 + U301*Two
         W(I,L,76,To) = QBx*S113 + WQx*T113 + U013    
         W(I,L,77,To) = QBx*S023 + WQx*T023            
         W(I,L,78,To) = QBy*S023 + WQy*T023 + U013*Two
         W(I,L,79,To) = QBx*S104 + WQx*T104 + U004    
         W(I,L,80,To) = QBx*S014 + WQx*T014            
         W(I,L,81,To) = QBy*S014 + WQy*T014 + U004
         W(I,L,82,To) = QBx*S005 + WQx*T005        
         W(I,L,83,To) = QBy*S005 + WQy*T005        
         W(I,L,84,To) = QBz*S005 + WQz*T005 + U004*Five
      END DO
      END DO
      DO L = 2,NumL
         CALL KonK2L(L,Lx,Ly,Lz)
         IF(Lx.gt.0) THEN
            CALL KonL2K(LP,Lx-1,Ly,Lz)
            DO I = 1,N
               Y = WZ(5,I) * Lx
               W(I,L,57,To) = W(I,L,57,To)+W(I,LP,36,Fr)*Y
               W(I,L,60,To) = W(I,L,60,To)+W(I,LP,39,Fr)*Y
               W(I,L,61,To) = W(I,L,61,To)+W(I,LP,40,Fr)*Y
               W(I,L,62,To) = W(I,L,62,To)+W(I,LP,41,Fr)*Y
               W(I,L,76,To) = W(I,L,76,To)+W(I,LP,52,Fr)*Y
               W(I,L,77,To) = W(I,L,77,To)+W(I,LP,53,Fr)*Y
               W(I,L,79,To) = W(I,L,79,To)+W(I,LP,54,Fr)*Y
               W(I,L,80,To) = W(I,L,80,To)+W(I,LP,55,Fr)*Y
               W(I,L,82,To) = W(I,L,82,To)+W(I,LP,56,Fr)*Y
            END DO
         END IF
         IF(Ly.gt.0) THEN
            CALL KonL2K(LP,Lx,Ly-1,Lz)
            DO I = 1,N
               Y = WZ(5,I) * Ly
               W(I,L,58,To) = W(I,L,58,To)+W(I,LP,36,Fr)*Y
               W(I,L,59,To) = W(I,L,59,To)+W(I,LP,37,Fr)*Y
               W(I,L,63,To) = W(I,L,63,To)+W(I,LP,41,Fr)*Y
               W(I,L,65,To) = W(I,L,65,To)+W(I,LP,42,Fr)*Y
               W(I,L,66,To) = W(I,L,66,To)+W(I,LP,43,Fr)*Y
               W(I,L,71,To) = W(I,L,71,To)+W(I,LP,47,Fr)*Y
               W(I,L,78,To) = W(I,L,78,To)+W(I,LP,53,Fr)*Y
               W(I,L,81,To) = W(I,L,81,To)+W(I,LP,55,Fr)*Y
               W(I,L,83,To) = W(I,L,83,To)+W(I,LP,56,Fr)*Y
            END DO
         END IF
         IF(Lz.gt.0) THEN
            CALL KonL2K(LP,Lx,Ly,Lz-1)
            DO I = 1,N
               Y = WZ(5,I) * Lz
               W(I,L,64,To) = W(I,L,64,To)+W(I,LP,36,Fr)*Y
               W(I,L,67,To) = W(I,L,67,To)+W(I,LP,39,Fr)*Y
               W(I,L,68,To) = W(I,L,68,To)+W(I,LP,40,Fr)*Y
               W(I,L,69,To) = W(I,L,69,To)+W(I,LP,41,Fr)*Y
               W(I,L,70,To) = W(I,L,70,To)+W(I,LP,42,Fr)*Y
               W(I,L,72,To) = W(I,L,72,To)+W(I,LP,44,Fr)*Y
               W(I,L,73,To) = W(I,L,73,To)+W(I,LP,45,Fr)*Y
               W(I,L,74,To) = W(I,L,74,To)+W(I,LP,46,Fr)*Y
               W(I,L,75,To) = W(I,L,75,To)+W(I,LP,47,Fr)*Y
               W(I,L,84,To) = W(I,L,84,To)+W(I,LP,56,Fr)*Y
            END DO
         END IF
      END DO
      IF(MxLKet.le.LTot-Knt+6) GO TO 20
      CALL MondoHalt('MxLKet > 6 in VRRG')

20    CONTINUE
      RETURN
      END

      FUNCTION LFuncC(Lmin,Lmax)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER LFuncC,Lmin,Lmax
      IF(Lmax.GE.Lmin) THEN
         LFuncC=((Lmax+1)*(Lmax+2)*(Lmax+3) -
     >            Lmin*(Lmin+1)*(Lmin+2))/6      
      ELSE
         LFuncC=0
      ENDIF
      RETURN
      END

      SUBROUTINE KonL2K(L,Lx,Ly,Lz)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER L,Lx,Ly,Lz
      L=LMNDex(Lx,Ly,Lz)
      RETURN
      END
