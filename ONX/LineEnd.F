C
C     Matt Challacombe
C     Los Alamos National Laboratory
C     Copywrite 1997, The University of California
C
      INTEGER FUNCTION LineEnd(String)
C
      CHARACTER*(*) String
C
      L=LEN(String)
      DO I=1,L
         IF(String(I:I).EQ.' ')THEN
            LineEnd=I-1
            RETURN
         ENDIF
      ENDDO
C
      LineEnd=L
C
      RETURN
      END
C
