      SUBROUTINE VRRL(N,NVRR,NR,NS,SLOC,VLOC,W,R,WR,WZ)
      IMPLICIT NONE
      INTEGER N,NR,NS,NT,I,J,NVRR
      INTEGER I1,I2,I3,I4,I5,I6,I7
      INTEGER I8,I9,I10,I11,I12,I13
      INTEGER SLOC(NR),VLOC(10,NS)
      REAL*8  F9,F10
      REAL*8  W(N,NVRR),R(N,NR)
      REAL*8  WR(12,N),WZ(5,N)

      DO J=1,NR
        I1=SLOC(J)
        DO I=1,N
          W(I,I1)=R(I,J)
        ENDDO
      ENDDO
      
      DO J=1,NS
        NT=VLOC(1,J)
        IF(NT.EQ.1) THEN
          I2=VLOC(2,J)          
          I3=VLOC(3,J)          
          I4=VLOC(4,J)          
          I5=VLOC(5,J)          
          I11=I2+6              
          DO I=1,N
            W(I,I3)=WR(I2,I)*W(I,I4)+WR(I11,I)*W(I,I5)
          ENDDO
        ELSEIF(NT.EQ.2) THEN
          I2=VLOC(2,J)         
          I3=VLOC(3,J)         
          I4=VLOC(4,J)         
          I5=VLOC(5,J)         
          I7=VLOC(7,J)        
          I8=VLOC(8,J)
          F10=FLOAT(VLOC(10,J))  
          I11=I2+6             
          I12=MOD(I3,2)+1      
          DO I=1,N
            W(I,I3)=WR(I2,I)*W(I,I4)+WR(I11,I)*W(I,I5)
     >              +F10*WZ(5,I)*W(I,I8)
          ENDDO    
        ELSEIF(NT.EQ.3) THEN
          I2=VLOC(2,J)        
          I3=VLOC(3,J)        
          I4=VLOC(4,J)        
          I5=VLOC(5,J)       
          I6=VLOC(6,J)
          I7=VLOC(7,J)         
          F9=FLOAT(VLOC(9,J))
          I11=I2+6             
          I12=MOD(I2,2)+1
          I13=I12+2            
          DO I=1,N
            W(I,I3)=WR(I2,I)*W(I,I4)+WR(I11,I)*W(I,I5)
     >              +F9*WZ(I12,I)*(W(I,I6)-WZ(I13,I)*W(I,I7))
          ENDDO  
        ELSEIF(NT.EQ.4) THEN
          I2=VLOC(2,J)          ! IPA
          I3=VLOC(3,J)          ! ILOC1
          I4=VLOC(4,J)          ! ILOC2
          I5=VLOC(5,J)          ! ILOC3
          I6=VLOC(6,J)          ! ILOC4
          I7=VLOC(7,J)          ! ILOC5
          I8=VLOC(8,J)          ! ILOC6
          F9=FLOAT(VLOC(9,J))   ! F1
          F10=FLOAT(VLOC(10,J)) ! F2
          I11=I2+6              ! IWP
          I12=MOD(I2,2)+1       ! IZE
          I13=I12+2             ! IZE2
          DO I=1,N

            W(I,I3)=WR(I2,I)*W(I,I4)+WR(I11,I)*W(I,I5)
     >              +F10*WZ(5,I)*W(I,I8)
     >              +F9*WZ(I12,I)*(W(I,I6)-WZ(I13,I)*W(I,I7))
          ENDDO   
        ELSE
          WRITE(*,*) "Illegal NT in VRR: ",NT
          STOP
        ENDIF
      ENDDO
      RETURN
      END
