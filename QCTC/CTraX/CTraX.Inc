      CASE(4)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX4(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(104)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX104(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(204)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX204(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(6)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX6(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(106)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX106(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(206)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX206(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(306)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX306(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(406)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX406(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(8)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX8(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(108)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX108(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(208)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX208(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(308)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX308(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(408)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX408(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(508)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX508(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(608)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX608(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(10)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX10(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(110)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX110(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(210)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX210(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(310)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX310(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(410)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX410(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(510)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX510(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(610)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX610(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(710)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX710(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(12)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX12(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(112)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX112(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(212)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX212(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(312)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX312(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(412)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX412(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(512)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX512(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(612)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX612(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(712)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX712(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(14)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX14(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(114)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX114(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(214)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX214(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(314)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX314(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(414)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX414(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(514)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX514(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(614)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX614(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(714)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX714(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(16)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX16(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(116)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX116(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(216)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX216(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(316)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX316(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(416)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX416(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(516)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX516(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(616)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX616(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(716)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX716(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(18)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX18(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(118)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX118(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(218)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX218(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(318)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX318(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(418)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX418(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(518)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX518(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(618)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX618(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(718)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX718(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(20)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX20(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(120)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX120(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(220)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX220(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(320)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX320(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(420)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX420(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(520)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX520(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(620)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX620(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE(720)
         DO N=1,NFar
            Q=>Far(N)%P
            PQx=QC%Prim%P(1)-Q%Pole%Center(1)
            PQy=QC%Prim%P(2)-Q%Pole%Center(2)
            PQz=QC%Prim%P(3)-Q%Pole%Center(3)
            CALL CTraX720(PQx,PQy,PQz,Q%Pole%C(0),Q%Pole%S(0),SPKetC(0),SPKetS(0))
         ENDDO
      CASE DEFAULT
         STOP ' Unsuported expansion orders in TreeWalk.F90'
      END SELECT
