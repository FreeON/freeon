SUBROUTINE Digest(N,NA,NB,NC,ND,L1,L2,L3,L4,IntSwitch,K,W,D)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE ONXParameters
  USE ONXMemory
  IMPLICIT NONE
  INTEGER,INTENT(IN)        :: N,IntSwitch
  INTEGER,INTENT(IN)        :: NA,NB,NC,ND
  INTEGER,INTENT(IN)        :: L1,L2,L3,L4
  REAL(DOUBLE),INTENT(OUT)  :: K(N,NA,NB)
  REAL(DOUBLE),INTENT(IN)   :: W(N,L1,L2,L3,L4)
  REAL(DOUBLE),INTENT(IN)   :: D(NC,ND)
  INTEGER                   :: I,IA,IB,IC,ID
  REAL(DOUBLE)              :: Dcd

  K=0.0D0
  SELECT CASE (IntSwitch)
  CASE (11)
    DO IB=1,L4
    DO ID=1,L3
    DO IA=1,L2
    DO IC=1,L1
      Dcd=D(IC,ID)
      DO I=1,N
        K(I,IA,IB)=K(I,IA,IB)-W(I,IC,IA,ID,IB)*Dcd
      END DO
    END DO
    END DO
    END DO
    END DO
  CASE (01)
    DO IB=1,L4
    DO ID=1,L3
    DO IC=1,L2
      Dcd=D(IC,ID)
    DO IA=1,L1
      DO I=1,N
        K(I,IA,IB)=K(I,IA,IB)-W(I,IA,IC,ID,IB)*Dcd
      END DO
    END DO
    END DO
    END DO
    END DO
  CASE (10)
    DO ID=1,L4
    DO IB=1,L3
    DO IA=1,L2
    DO IC=1,L1
      Dcd=D(IC,ID)
      DO I=1,N
        K(I,IA,IB)=K(I,IA,IB)-W(I,IC,IA,IB,ID)*Dcd
      END DO
    END DO
    END DO
    END DO
    END DO
  CASE (00)
    DO ID=1,L4
    DO IB=1,L3
    DO IC=1,L2
      Dcd=D(IC,ID)
    DO IA=1,L1
      DO I=1,N
        K(I,IA,IB)=K(I,IA,IB)-W(I,IA,IC,IB,ID)*Dcd
      END DO
    END DO
    END DO
    END DO
    END DO
  CASE DEFAULT
    CALL Halt(' Illegal IntSwitch in ONX:Digest')
  END SELECT
END SUBROUTINE Digest
