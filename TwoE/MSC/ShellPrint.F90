SUBROUTINE ShellPrint(NBFA,NBFB,NBFC,NBFD,CFA,CFB,CFC,CFD, &
                      AtAO,AtBO,AtCO,AtDO,OFFA,OFFB,OFFC,OFFD,IntType,I)
  USE DerivedTypes
  INTEGER :: NBFA,NBFB,NBFC,NBFD,CFA,CFB,CFC,CFD,AtAO,AtBO,AtCO,AtDO,OFFA,OFFB,OFFC,OFFD
  INTEGER :: LDA,LDB,LDC,LDD,IA,IB,IC,ID,OffSet,IntType,OA,OB,OC,OD
  REAL(DOUBLE) :: I(*)
!  WRITE(*,*)'(* ============= IntType = ',IntType,'======================*)'
  LDA=1
  LDB=NBFA
  LDC=NBFA*NBFB 
  LDD=NBFA*NBFB*NBFC
  OA=OFFA-1
  OB=OFFB-1
  OC=OFFC-1
  OD=OFFD-1

  DO IA=1,CFA
     DO IC=1,CFC
        DO IB=1,CFB
           DO ID=1,CFD
              OffSet=(OA+IA-1)*LDA+(OB+IB-1)*LDB+(OC+IC-1)*LDC+(OD+ID-1)*LDD+1
              
!!$              IF(AtAO+IA+OA==23.AND.AtCO+IC+OC==3.AND.AtBO+IB+OB==10.AND.AtDO+ID+OD==2)THEN
!!$                 WRITE(*,*)'IntType = ',IntType
!!$                 WRITE(*,*)' OFF = ',OA,OB,OC,OD
!!$                 WRITE(*,*)'  CF = ',CFA,CFC,CFB,CFD
!!$                 WRITE(*,*)' III = ',IA,IC,IB,ID
!!$                 WRITE(*,101)OffSet,AtAO+IA+OA,AtCO+IC+OC,AtBO+IB+OB,AtDO+ID+OD,I(OffSet)
!!$                 STOP
!!$              ENDIF

              IF(ABS(I(OffSet))>1D-12)THEN 
!                 WRITE(*,101)OffSet,IA+OA,IC+OC,IB+OB,ID+OD,I(OffSet)
                 WRITE(44,301)AtAO+IA+OA,AtCO+IC+OC,AtBO+IB+OB,AtDO+ID+OD, &
                          FRACTION(I(OffSet)),EXPONENT(I(OffSet))
              ENDIF
              IF(ABS(I(OffSet))>1D10)THEN
                 WRITE(*,*)'IntType = ',IntType
                 WRITE(*,101)OffSet,AtAO+IA+OA,AtCO+IC+OC,AtBO+IB+OB,AtDO+ID+OD,I(OffSet)
                 STOP ' BAD INT !!!'
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO
101 FORMAT(I10,' (',I3,',',I3,'|',I3,',',I3,') = ',F20.16)
201 FORMAT(' (',I3,',',I3,'|',I3,',',I3,') = ',D22.16)
301 FORMAT(1x,'Int2[[',I3,',',I3,',',I3,',',I3,']] = ',F19.16,'*2^(',I3,');')
END SUBROUTINE ShellPrint
