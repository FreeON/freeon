SUBROUTINE RangeOfExchange(BSc,GMc,BSp,GMp,D,NameBuf)
  USE DerivedTypes
  USE GlobalScalars
  USE Macros
  USE MemMan
  USE MondoMPI
  USE Stats
  USE ONXParameters
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! Basis set, coordinates, ect.
!--------------------------------------------------------------------------------
  TYPE(BSET),INTENT(IN)    :: BSc,BSp
  TYPE(CRDS),INTENT(IN)    :: GMc,GMp
#ifdef PARALLEL
  TYPE(DBCSR)              :: D
#else
  TYPE(BCSR)               :: D
#endif
  INTEGER                  :: AtA,AtB,AtC,AtD
  INTEGER                  :: ri,ci
  INTEGER                  :: KA,KB,NBFA,NBFB
  REAL(DOUBLE)             :: Ax,Ay,Az,Bx,By,Bz
  REAL(DOUBLE)             :: Cx,Cy,Cz,Dx,Dy,Dz
  REAL(DOUBLE)             :: AB2,AC2,CD2
  TYPE(INT_VECT)           :: NameBuf
  TYPE(DBL_VECT)           :: MR


  NRows=0
  NCols=0
  NElem=0
  NameBuf%I=0
  DenRange=0.0D0
#ifdef PARALLEL
  DO AtC=Beg%I(MyID),End%I(MyID)
    ri=AtC-Beg%I(MyID)+1
#else
  DO AtC=1,NAtoms
    ri=AtC
#endif
    Cx=GMp%Carts%D(1,AtC)
    Cy=GMp%Carts%D(2,AtC)
    Cz=GMp%Carts%D(3,AtC)
    DO ci=D%RowPt%I(ri),D%RowPt%I(ri+1)-1
      AtD=D%ColPt%I(ci)
      Dx=GMp%Carts%D(1,AtD)
      Dy=GMp%Carts%D(2,AtD)
      Dz=GMp%Carts%D(3,AtD)
      CD2=(Cx-Dx)*(Cx-Dx) + &
          (Cy-Dy)*(Cy-Dy) + &
          (Cz-Dz)*(Cz-Dz)
      DenRange=MAX(DenRange,SQRT(CD2)*1.1D0)
    END DO
  END DO

  ONXRange=DenRange+2.0D0*DisRange

  DO AtA=1,NAtoms
    Ax=GMc%Carts%D(1,AtA)
    Ay=GMc%Carts%D(2,AtA)
    Az=GMc%Carts%D(3,AtA)
#ifdef PARALLEL
    DO AtC=Beg%I(MyID),End%I(MyID)
#else
    DO ATC=1,NAtoms
#endif
      Cx=GMp%Carts%D(1,AtC)
      Cy=GMp%Carts%D(2,AtC)
      Cz=GMp%Carts%D(3,AtC)
      AC2=(Ax-Cx)*(Ax-Cx) + &
          (Ay-Cy)*(Ay-Cy) + &
          (Az-Cz)*(Az-Cz)
      IF (SQRT(AC2).LE.DisRange) THEN
        NRows=NRows+1
        NameBuf%I(NRows)=AtA
        EXIT
      END IF
    END DO ! AtC
  END DO ! AtA

!
! Loop over elements of the exchange matrix
!
  DO ri=1,NRows
    AtA=NameBuf%I(ri)
    KA=GMc%AtTyp%I(AtA)
    NBFA=BSc%BfKnd%I(KA)
    Ax=GMc%Carts%D(1,AtA)
    Ay=GMc%Carts%D(2,AtA)
    Az=GMc%Carts%D(3,AtA)
    DO ci=1,NRows
      AtB=NameBuf%I(ci)
      KB=GMc%AtTyp%I(AtB)
      NBFB=BSc%BfKnd%I(KB)
      Bx=GMc%Carts%D(1,AtB)
      By=GMc%Carts%D(2,AtB)
      Bz=GMc%Carts%D(3,AtB)
      AB2=(Ax-Bx)*(Ax-Bx) + &
          (Ay-By)*(Ay-By) + &
          (Az-Bz)*(Az-Bz)
      IF (SQRT(AB2).LE.ONXRange) THEN
        NCols=NCols+1
        NElem=NElem+NBFA*NBFB
      END IF
    END DO ! ci
  END DO ! ri

#ifdef PARALLEL
  IF (MyID==ROOT) THEN
    CALL New(MR,NPrc,0)
  END IF
  CALL Gather(ONXRange,MR)
  IF (MyID==ROOT) THEN
    MatRange=GetMax(NPrc,MR,0)
    CALL Delete(MR)
  END IF
#else
  MatRange=ONXRange
#endif

END SUBROUTINE RangeOfExchange

