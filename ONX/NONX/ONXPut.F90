MODULE ONXPut
!H=================================================================================
!H MODULE ONXPut
!H This MODULE contains:
!H  o SUB PutDis
!H  o SUB PutSubBlk
!H---------------------------------------------------------------------------------
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  !
  IMPLICIT NONE
  PRIVATE
  !
!--------------------------------------------------------------------------------- 
! PUBLIC DECLARATIONS
!--------------------------------------------------------------------------------- 
  PUBLIC :: PutDis
  PUBLIC :: PutSubBlk
  !
CONTAINS
  !
  SUBROUTINE PutDis(N,iDis,iPrm,I,J,IndexA,KonAC,BufT,DB)
!H--------------------------------------------------------------------------------- 
!H SUBROUTINE PutDis(N,iDis,iPrm,I,J,IndexA,KonAC,BufT,DB)
!H
!H---------------------------------------------------------------------------------
   IMPLICIT NONE
    TYPE(DBuf),INTENT(INOUT)  :: DB
    TYPE(INT_RNK3),INTENT(IN) :: BufT
    INTEGER,INTENT(IN)        :: N,I,J,IndexA,KonAC
    INTEGER,INTENT(INOUT)     :: iDis,iPrm
    INTEGER                   :: I0,I1,I2,I3,ILoc

    DB%DisPtr%I(1,IndexA,J,I)=N
    DB%DisPtr%I(2,IndexA,J,I)=iDis
    DB%DisPtr%I(3,IndexA,J,I)=iPrm
    DO I1=1,N
      ILoc=BufT%I(I1,J,I)
      DO I2=1,DB%MAXC
        I0=I2+(I1-1)*DB%MAXC+iDis-1
        DB%DisBuf%D(I0) = DB%TBufC%D(I2,ILoc)
      END DO  ! I2
      DO I2=1,KonAC
        DO I3=1,DB%MAXP
          I0=I3+(I2-1)*DB%MAXP+(I1-1)*DB%MAXP*(KonAC+DB%MInfo)+iPrm-1
          DB%PrmBuf%D(I0) = DB%TBufP%D(I3,I2,ILoc)
        END DO  ! I3
      END DO  ! I2
    END DO  ! I1
    iDis = iDis + N*DB%MAXC
    iPrm = iPrm + N*(KonAC+DB%MInfo)*DB%MAXP
  END SUBROUTINE PutDis
  !
  !
  SUBROUTINE PutSubBlk(I,N,A1,A2,B1,B2,RS,CS,A,B)
!H--------------------------------------------------------------------------------- 
!H SUBROUTINE PutSubBlk(I,N,A1,A2,B1,B2,RS,CS,A,B)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN)         :: I,N,A1,A2,B1,B2,RS,CS
    REAL(DOUBLE),INTENT(INOUT) :: A(A1,A2)
    REAL(DOUBLE),INTENT(IN)    :: B(N,B1,B2)
    INTEGER                    :: IC,IR,RS1,CS1
    !
    RS1=RS-1
    CS1=CS-1
    DO IC=1,B2
       DO IR=1,B1
          A(IR+RS1,IC+CS1)=A(IR+RS1,IC+CS1)+B(I,IR,IC)
       END DO
    END DO
    !
  END SUBROUTINE PutSubBlk
  !
END MODULE ONXPut
