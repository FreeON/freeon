      SUBROUTINE KonK2L(K,Lx,Ly,Lz)
      IMPLICIT INTEGER (a-z)
      PARAMETER (MaxK=84)
      INTEGER Li(3,MaxK)
      SAVE   Li
      DATA ((Li(i,j),i=1,3),j=1,35)
     $        / 0,0,0,
     $          1,0,0,  0,1,0,
     $                  0,0,1,
     $          2,0,0,  1,1,0,  0,2,0,
     $                  1,0,1,  0,1,1,
     $                          0,0,2,
     $          3,0,0,  2,1,0,  1,2,0,  0,3,0,
     $                  2,0,1,  1,1,1,  0,2,1,
     $                          1,0,2,  0,1,2,
     $                                  0,0,3,
     $          4,0,0,  3,1,0,  2,2,0,  1,3,0,  0,4,0,
     $                  3,0,1,  2,1,1,  1,2,1,  0,3,1,
     $                          2,0,2,  1,1,2,  0,2,2,
     $                                  1,0,3,  0,1,3,
     $                                          0,0,4 /
      DATA ((Li(i,j),i=1,3),j=36,MaxK)
     $        / 5,0,0,  4,1,0,  3,2,0,  2,3,0,  1,4,0,  0,5,0,
     $                  4,0,1,  3,1,1,  2,2,1,  1,3,1,  0,4,1,
     $                          3,0,2,  2,1,2,  1,2,2,  0,3,2,
     $                                  2,0,3,  1,1,3,  0,2,3,
     $                                          1,0,4,  0,1,4,
     $                                                  0,0,5,
     $          6,0,0,  5,1,0,  4,2,0,  3,3,0,  2,4,0,  1,5,0,  0,6,0,
     $                  5,0,1,  4,1,1,  3,2,1,  2,3,1,  1,4,1,  0,5,1,
     $                          4,0,2,  3,1,2,  2,2,2,  1,3,2,  0,4,2,
     $                                  3,0,3,  2,1,3,  1,2,3,  0,3,3,
     $                                          2,0,4,  1,1,4,  0,2,4,
     $                                                  1,0,5,  0,1,5,
     $                                                          0,0,6 /
      IF(K.le.0) THEN

         CALL MondoHalt('Non-positive K in KonK2L')

      ELSE IF(K.le.MaxK) THEN

         Lx = Li(1,K)
         Ly = Li(2,K)
         Lz = Li(3,K)

      ELSE

c        We are outside the range of our lookup table so we first find
c        the total angular momentum L and then increment from there...

         DO 30 L = 0,10000
            Num = (L+1)*(L+2)*(L+3)/6 + 1
            IF(K.lt.Num) THEN

c              Now count up until we find the Kth orbital...

               Knt = L*(L+1)*(L+2)/6
               DO 20 Lz = 0,L
                  DO 10 Ly = 0,L-Lz
                     Lx  = L - Ly - Lz
                     Knt = Knt + 1
                     IF(K.eq.Knt) RETURN
10                CONTINUE
20             CONTINUE
               CALL MondoHalt('Yuk #1 in KonK2L')
            END IF
30       CONTINUE
         CALL MondoHalt('Yuk #2 in KonK2L')
      END IF
      END

