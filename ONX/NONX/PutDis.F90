  SUBROUTINE PutDis(N,iDis,iPrm,I,J,IndexA,KonAC,DB)
    USE DerivedTypes
    IMPLICIT NONE
    TYPE(DBuf),INTENT(INOUT) :: DB
    INTEGER,INTENT(IN)       :: N,I,J,IndexA,KonAC
    INTEGER,INTENT(INOUT)    :: iDis,iPrm
    INTEGER                  :: I0,I1,I2,I3,ILoc

    DB%DisPtr%I(1,IndexA,J,I)=N
    DB%DisPtr%I(2,IndexA,J,I)=iDis
    DB%DisPtr%I(3,IndexA,J,I)=iPrm
    DO I1=1,N
      ILoc=DB%BufT%I(I1,J,I)
      DO I2=1,DB%MAXC
        I0=I2+(I1-1)*DB%MAXC
        DB%DisBuf%D(I0) = DB%TBufC%D(I2,ILoc)
      END DO  ! I2
      DO I2=1,KonAC
        DO I3=1,DB%MAXP
          I0=I3+(I2-1)*DB%MAXP+(I-1)*DB%MAXP*(KonAC+DB%MInfo)
          DB%PrmBuf%D(I0) = DB%TBufP%D(I3,I2,ILoc)
        END DO  ! I3
      END DO  ! I2
    END DO  ! I1
    iDis = iDis + N*DB%MAXC
    iPrm = iPrm + N*(KonAC+DB%MInfo)*DB%MAXP
  END SUBROUTINE PutDis
