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
      CASE DEFAULT
         STOP ' Unsuported expansion orders in TreeWalk.F90'
      END SELECT
