SUBROUTINE CompareBlock(AtAOff,AtBOff,AtCOff,AtDOff,NA,NB,NC,ND,I,dI)
  USE DerivedTypes
  INTEGER :: AtAOff,AtBOff,AtCOff,AtDOff,NA,NB,NC,ND,IA,IB,IC,ID
  REAL(DOUBLE) :: Div
  REAL(DOUBLE), DIMENSION (1:NA,1:NC,1:NB,1:ND) :: I,dI
  DO IA=1,NA
     DO IC=1,NC
        DO IB=1,NB
           DO ID=1,ND
              Div=1.0d0
              IF(ABS(dI(IA,IC,IB,ID)).GT.1d-16)Div=dI(IA,IC,IB,ID)
              IF(ABS((I(IA,IC,IB,ID)-dI(IA,IC,IB,ID))/Div)> 1D-4)THEN 
                 WRITE(*,101)IA+AtAOff,IC+AtCOff,IB+AtBOff,ID+AtDOff,I(IA,IC,IB,ID),dI(IA,IC,IB,ID)
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO
!101 FORMAT(' (',I3,',',I3,'|',I3,',',I3,') = ',F20.16,', ',F20.16)
101 FORMAT(' (',I3,',',I3,'|',I3,',',I3,') = ',E24.16,', ',E24.16)
201 FORMAT(' (',I3,',',I3,'|',I3,',',I3,') = ',D22.16)
301 FORMAT(1x,'Int2[[',I3,',',I3,',',I3,',',I3,']] = ',F19.16,'*2^(',I3,');')
  
END SUBROUTINE CompareBlock

PROGRAM dN4KTest
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE Macros
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE ShellPairStruct
  USE GetAtomPairMod
  IMPLICIT NONE
  !-------------------------------------------------------------------

  !  TYPE(CList2) , DIMENSION(:), POINTER  :: List
  TYPE(CRDS)              :: GMc
  TYPE(BSET)              :: BS
  TYPE(BCSR)              :: D
  !  TYPE(CellSet)           :: CS_OUT
  !  TYPE(ANode2) , POINTER                :: AtAList,AtAListTmp,NodeA
  TYPE(AtomInfo)                        :: ACAtmInfo,BDAtmInfo
  INTEGER                               :: AtA,AtC,AtB,AtD,KA,KC,KB,KD, &
       CFA,CFC,CFB,CFD,iFAC,iFBD,NBFA,NBFC,NBFB,NBFD,AtAOff,AtBOff,AtCOff,AtDOff, & 
       AALen,BBLen,CCLen,DDLen,NIntBlk
  INTEGER                               :: NCell,IntType,LocNInt,Off,iDir,iDiff,I
  INTEGER                               :: indab,indcd,iptrdab,iptrdcd,IXYZ
  REAL(DOUBLE)                          :: Tmp1,Grad
  TYPE(AtomPr), DIMENSION(100)          :: ACAtmPair,BDAtmPair
  TYPE(ONX2OffSt) :: OffSet
  INTEGER, PARAMETER                     :: MaxBK=1000
  REAL(DOUBLE) , DIMENSION(MaxBK)       :: C,Cm,Cp
  REAL(DOUBLE) , DIMENSION(MaxBK,12)    :: dC
  REAL(DOUBLE) , DIMENSION(:,:), ALLOCATABLE :: GradKx
  TYPE(INT_VECT)                        :: BColIdx,DColIdx
  TYPE(DBL_RNK2)                             :: GradAux
  REAL(DOUBLE) , EXTERNAL               :: DGetAbsMax
  TYPE(ARGMT)                    :: Args
  CHARACTER(LEN=4) :: Prog="dN4K"

  CALL StartUp(Args,Prog)

  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GMc,Tag_O=CurGeom)
  
  CALL Get(D,TrixFile('D',Args,0))

  AtAOff=0
  DO AtA=1,NAtoms ! Run over AtA
     KA=GMc%AtTyp%I(AtA)
     ACAtmInfo%K1=KA
     NBFA=BS%BfKnd%I(KA)
     ACAtmInfo%Atm1X=GMc%Carts%D(1,AtA)
     ACAtmInfo%Atm1Y=GMc%Carts%D(2,AtA)
     ACAtmInfo%Atm1Z=GMc%Carts%D(3,AtA)
     AtCOff=0
     DO AtC=1,NAtoms ! Run over AtC
        KC=GMc%AtTyp%I(AtC)
        ACAtmInfo%K2=KC
        NBFC=BS%BfKnd%I(KC)
        ACAtmInfo%Atm2X=GMc%Carts%D(1,AtC)
        ACAtmInfo%Atm2Y=GMc%Carts%D(2,AtC)
        ACAtmInfo%Atm2Z=GMc%Carts%D(3,AtC)
        AtBOff=0
        DO AtB=1,NAtoms ! Run over AtB
           KB=GMc%AtTyp%I(AtB)
           NBFB=BS%BfKnd%I(KB)
           BDAtmInfo%K1=KB
           BDAtmInfo%Atm1X=GMc%Carts%D(1,AtB)
           BDAtmInfo%Atm1Y=GMc%Carts%D(2,AtB)
           BDAtmInfo%Atm1Z=GMc%Carts%D(3,AtB)
           AtDOff=0
           DO AtD=1,NAtoms ! Run over AtD
              KD=GMc%AtTyp%I(AtD)
              NBFD=BS%BfKnd%I(KD)
              BDAtmInfo%K2=KD
              BDAtmInfo%Atm2X=GMc%Carts%D(1,AtD)
              BDAtmInfo%Atm2Y=GMc%Carts%D(2,AtD)
              BDAtmInfo%Atm2Z=GMc%Carts%D(3,AtD)
              DO iDir=3,12,3              
                 IF(iDir==3)THEN
                    Tmp1=ACAtmInfo%Atm1Z
                 ELSEIF(iDir==6)THEN
                    Tmp1=ACAtmInfo%Atm2Z
                 ELSEIF(iDir==9)THEN
                    Tmp1=BDAtmInfo%Atm1Z
                 ELSEIF(iDir==12)THEN
                    Tmp1=BDAtmInfo%Atm2Z
                 ENDIF
                 !-------------------------------------------------------------------------------------------
                 ! Differencing !!!
                 DO iDiff=1,2
                    IF(iDiff==1)THEN
                       IF(iDir==3)THEN
                          ACAtmInfo%Atm1Z=Tmp1-1D-4
                       ELSEIF(iDir==6)THEN
                          ACAtmInfo%Atm2Z=Tmp1-1D-4
                       ELSEIF(iDir==9)THEN
                          BDAtmInfo%Atm1Z=Tmp1-1D-4
                       ELSEIF(iDir==12)THEN
                          BDAtmInfo%Atm2Z=Tmp1-1D-4
                       ENDIF
                    ELSE
                       IF(iDir==3)THEN
                          ACAtmInfo%Atm1Z=Tmp1+1D-4
                       ELSEIF(iDir==6)THEN
                          ACAtmInfo%Atm2Z=Tmp1+1D-4
                       ELSEIF(iDir==9)THEN
                          BDAtmInfo%Atm1Z=Tmp1+1D-4
                       ELSEIF(iDir==12)THEN
                          BDAtmInfo%Atm2Z=Tmp1+1D-4
                       ENDIF
                    ENDIF
                    CALL GetAtomPair2(ACAtmInfo,ACAtmPair,BS,CS_OUT%CellCarts%D(1,1))
                    CALL GetAtomPair2(BDAtmInfo,BDAtmPair,BS,CS_OUT%CellCarts%D(1,1))
                    IF(iDir==3)THEN
                       ACAtmInfo%Atm1Z=Tmp1
                    ELSEIF(iDir==6)THEN
                       ACAtmInfo%Atm2Z=Tmp1
                    ELSEIF(iDir==9)THEN
                       BDAtmInfo%Atm1Z=Tmp1
                    ELSEIF(iDir==12)THEN
                       BDAtmInfo%Atm2Z=Tmp1
                    ENDIF
                    CALL DBL_VECT_EQ_DBL_SCLR(NBFA*NBFB*NBFC*NBFD,C(1),0.0d0)
                    iFAC=0
                    OffSet%A=1
                    DO CFA=1,BS%NCFnc%I(KA)
                       AALen=BS%LStop%I(CFA,KA)-BS%LStrt%I(CFA,KA)+1
                       OffSet%C=1
                       DO CFC=1,BS%NCFnc%I(KC)
                          CCLen=BS%LStop%I(CFC,KC)-BS%LStrt%I(CFC,KC)+1
                          iFBD=0
                          OffSet%B=1
                          iFAC=iFAC+1
                          DO CFB=1,BS%NCFnc%I(KB)
                             BBLen=BS%LStop%I(CFB,KB)-BS%LStrt%I(CFB,KB)+1
                             OffSet%D=1
                             DO CFD=1,BS%NCFnc%I(KD)
                                iFBD=iFBD+1
                                DDLen=BS%LStop%I(CFD,KD)-BS%LStrt%I(CFD,KD)+1
                                IntType=ACAtmPair(iFAC)%SP%IntType*10000+BDAtmPair(iFBD)%SP%IntType
                                INCLUDE 'ERIInterfaceB.Inc'
                                OffSet%D=OffSet%D+DDLen
                             ENDDO ! End blkfunc on D
                             OffSet%B=OffSet%B+BBLen
                          ENDDO ! End blkfunc on B
                          OffSet%C=OffSet%C+CCLen
                       ENDDO ! End blkfunc on C
                       OffSet%A=OffSet%A+AALen                 !
                    ENDDO ! End blkfunc on A
                    IF(iDiff==1)THEN
                       Cm(1:NBFA*NBFB*NBFC*NBFD)=C(1:NBFA*NBFB*NBFC*NBFD)
                    ELSE
                       Cp(1:NBFA*NBFB*NBFC*NBFD)=C(1:NBFA*NBFB*NBFC*NBFD)
                    ENDIF
                 ENDDO
                 C(1:NBFA*NBFB*NBFC*NBFD)=(Half*1D4)*(Cp(1:NBFA*NBFB*NBFC*NBFD)-Cm(1:NBFA*NBFB*NBFC*NBFD))
                 !-------------------------------------------------------------------------------------------
                 ! Analytics !!! 
                 CALL GetAtomPair2(ACAtmInfo,ACAtmPair,BS,CS_OUT%CellCarts%D(1,1))
                 CALL GetAtomPair2(BDAtmInfo,BDAtmPair,BS,CS_OUT%CellCarts%D(1,1))
                 DO I=1,12
                    CALL DBL_VECT_EQ_DBL_SCLR(NBFA*NBFB*NBFC*NBFD,dC(1,I),BIG_DBL)
                    !CALL DBL_VECT_EQ_DBL_SCLR(NBFA*NBFB*NBFC*NBFD,dC(1,I),0D0)
                 ENDDO
                 NIntBlk=MaxBK ! NBFA*NBFB*NBFC*NBFD
                 iFAC=0
                 OffSet%A=1
                 DO CFA=1,BS%NCFnc%I(KA)
                    AALen=BS%LStop%I(CFA,KA)-BS%LStrt%I(CFA,KA)+1
                    OffSet%C=1
                    DO CFC=1,BS%NCFnc%I(KC)
                       CCLen=BS%LStop%I(CFC,KC)-BS%LStrt%I(CFC,KC)+1
                       iFBD=0
                       OffSet%B=1
                       iFAC=iFAC+1
                       DO CFB=1,BS%NCFnc%I(KB)
                          BBLen=BS%LStop%I(CFB,KB)-BS%LStrt%I(CFB,KB)+1
                          OffSet%D=1
                          DO CFD=1,BS%NCFnc%I(KD)
                             iFBD=iFBD+1
                             DDLen=BS%LStop%I(CFD,KD)-BS%LStrt%I(CFD,KD)+1
                             IntType=ACAtmPair(iFAC)%SP%IntType*10000+BDAtmPair(iFBD)%SP%IntType
                             INCLUDE 'dERIInterfaceB.Inc'
                             OffSet%D=OffSet%D+DDLen
                          ENDDO ! End blkfunc on D
                          OffSet%B=OffSet%B+BBLen
                       ENDDO ! End blkfunc on B
                       OffSet%C=OffSet%C+CCLen
                    ENDDO ! End blkfunc on C
                    OffSet%A=OffSet%A+AALen                 
                 ENDDO ! End blkfunc on A
                 CALL CompareBlock(AtAOff,AtBOff,AtCOff,AtDOff,NBFA,NBFB,NBFC,NBFD,C(1),dC(1,iDir))
              ENDDO ! iDir
              AtDOff=AtDOff+NBFD
           ENDDO ! D
           AtBOff=AtBOff+NBFB
        ENDDO ! B
        AtCOff=AtCOff+NBFC
     ENDDO ! C
     AtAOff=AtAOff+NBFA
  ENDDO ! A

END PROGRAM dN4KTest

