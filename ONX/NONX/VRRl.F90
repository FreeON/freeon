SUBROUTINE VRRl(N,NVRR,NR,NS,SLOC,VLOC,IB,W,R)
  USE DerivedTypes
  USE GlobalScalars
  USE ONXParameters
  USE PrettyPrint
  IMPLICIT NONE

  INTEGER            :: N,NVRR,NR,NS
  INTEGER            :: SLOC(NR),VLOC(10,NS)
  TYPE(IBuf)         :: IB
  REAL(DOUBLE)       :: W(N,NVRR)
  REAL(DOUBLE)       :: R(N,NR)
  INTEGER            :: I,J,NT
  INTEGER            :: I1,I2,I3,I4,I5,I6,I7
  INTEGER            :: I8,I9,I10,I11,I12,I13
  REAL(DOUBLE)       :: F9,F10

  IF (N.EQ.0) RETURN
  
  DO J=1,NR
    I1=SLOC(J)
    DO I=1,N
      W(I,I1)=R(I,J)
    ENDDO
  ENDDO

  DO J=1,NS
    NT=VLOC(1,J)
    SELECT CASE (NT)
      CASE (1)
        I2=VLOC(2,J)
        I3=VLOC(3,J)
        I4=VLOC(4,J)
        I5=VLOC(5,J)
        I11=I2+6
        DO I=1,N
          W(I,I3)=IB%WR%D(I,I2)*W(I,I4)+IB%WR%D(I,I11)*W(I,I5)
        ENDDO
      CASE (2)
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
          W(I,I3)=IB%WR%D(I,I2)*W(I,I4)+IB%WR%D(I,I11)*W(I,I5)          &
                  +F10*IB%WZ%D(I,5)*W(I,I8)
        ENDDO   
      CASE (3)
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
          W(I,I3)=IB%WR%D(I,I2)*W(I,I4)+IB%WR%D(I,I11)*W(I,I5)         &
                  +F9*IB%WZ%D(I,I12)*(W(I,I6)-IB%WZ%D(I,I13)*W(I,I7))
        ENDDO
      CASE (4)
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
          W(I,I3)=IB%WR%D(I,I2)*W(I,I4)+IB%WR%D(I,I11)*W(I,I5)         &
                  +F10*IB%WZ%D(I,5)*W(I,I8)                            &
                  +F9*IB%WZ%D(I,I12)*(W(I,I6)-IB%WZ%D(I,I13)*W(I,I7))
        ENDDO  
      CASE DEFAULT
        WRITE(*,*) "NT=",NT
        CALL Halt(' Illegal NT in ONX:VRRl')
    END SELECT
  ENDDO 
END SUBROUTINE VRRl
