SUBROUTINE CompareBlock(AtAOff,AtBOff,AtCOff,AtDOff,NA,NB,NC,ND,I,dI)
  USE DerivedTypes
  INTEGER :: AtAOff,AtBOff,AtCOff,AtDOff,NA,NB,NC,ND,IA,IB,IC,ID
  REAL(DOUBLE), DIMENSION (1:NA,1:NC,1:NB,1:ND) :: I,dI
  DO IA=1,NA
     DO IC=1,NC
        DO IB=1,NB
           DO ID=1,ND
              IF(ABS(I(IA,IC,IB,ID)-dI(IA,IC,IB,ID))> 1D-4)THEN 
              WRITE(*,101)IA+AtAOff,IC+AtCOff,IB+AtBOff,ID+AtDOff,I(IA,IC,IB,ID),dI(IA,IC,IB,ID)
!              WRITE(33,301)IA+AtAOff,IC+AtCOff,IB+AtBOff,ID+AtDOff, &
!                          FRACTION(I(IA,IC,IB,ID)),EXPONENT(I(IA,IC,IB,ID))
              !WRITE(33,301)IA+AtAOff,IC+AtCOff,IB+AtBOff,ID+AtDOff,I(IA,IC,IB,ID)
              !WRITE(iOut,301)I,J,K,L,FRACTION(Int),EXPONENT(Int)
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO
101 FORMAT(' (',I3,',',I3,'|',I3,',',I3,') = ',F20.16,', ',F20.16)
201 FORMAT(' (',I3,',',I3,'|',I3,',',I3,') = ',D22.16)
301 FORMAT(1x,'Int2[[',I3,',',I3,',',I3,',',I3,']] = ',F19.16,'*2^(',I3,');')
  
END SUBROUTINE CompareBlock


PROGRAM N4KTest
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
  TYPE(CRDS)              :: GM
  TYPE(BSET)              :: BS
  !  TYPE(CellSet)           :: CS_OUT
  !  TYPE(ANode2) , POINTER                :: AtAList,AtAListTmp,NodeA
  TYPE(AtomInfo)                        :: ACAtmInfo,BDAtmInfo
  INTEGER                               :: AtA,AtC,AtB,AtD,KA,KC,KB,KD, &
       CFA,CFC,CFB,CFD,CFAC,CFBD,NBFA,NBFC,NBFB,NBFD,AtAOff,AtBOff,AtCOff,AtDOff, & 
       AALen,BBLen,CCLen,DDLen,NIntBlk
  INTEGER                               :: NCell,IntType,LocNInt,Off,iDir,iDiff,I
  REAL(DOUBLE)                          :: Tmp1
  TYPE(AtomPr), DIMENSION(100)          :: ACAtmPair,BDAtmPair
  TYPE(ONX2OffSt) :: OffSet
  INTEGER, PARAMETER                     :: MaxBK=1000
  REAL(DOUBLE) , DIMENSION(MaxBK)       :: C,Cm,Cp
  REAL(DOUBLE) , DIMENSION(MaxBK,12)    :: dC
  REAL(DOUBLE) , EXTERNAL               :: DGetAbsMax
  TYPE(ARGMT)                    :: Args
  CHARACTER(LEN=3) :: Prog="N4K"

  CALL StartUp(Args,Prog)

  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
  AtAOff=0
  DO AtA=1,NAtoms ! Run over AtA
     KA=GM%AtTyp%I(AtA)
     ACAtmInfo%K1=KA
     NBFA=BS%BfKnd%I(KA)
     ACAtmInfo%Atm1X=GM%Carts%D(1,AtA)
     ACAtmInfo%Atm1Y=GM%Carts%D(2,AtA)
     ACAtmInfo%Atm1Z=GM%Carts%D(3,AtA)
     AtCOff=0
     DO AtC=1,NAtoms ! Run over AtC
        KC=GM%AtTyp%I(AtC)
        ACAtmInfo%K2=KC
        NBFC=BS%BfKnd%I(KC)
        ACAtmInfo%Atm2X=GM%Carts%D(1,AtC)
        ACAtmInfo%Atm2Y=GM%Carts%D(2,AtC)
        ACAtmInfo%Atm2Z=GM%Carts%D(3,AtC)
        AtBOff=0
        DO AtB=1,NAtoms ! Run over AtB
           KB=GM%AtTyp%I(AtB)
           NBFB=BS%BfKnd%I(KB)
           BDAtmInfo%K1=KB
           BDAtmInfo%Atm1X=GM%Carts%D(1,AtB)
           BDAtmInfo%Atm1Y=GM%Carts%D(2,AtB)
           BDAtmInfo%Atm1Z=GM%Carts%D(3,AtB)
           AtDOff=0
           DO AtD=1,NAtoms ! Run over AtD
              KD=GM%AtTyp%I(AtD)
              NBFD=BS%BfKnd%I(KD)
              BDAtmInfo%K2=KD
              BDAtmInfo%Atm2X=GM%Carts%D(1,AtD)
              BDAtmInfo%Atm2Y=GM%Carts%D(2,AtD)
              BDAtmInfo%Atm2Z=GM%Carts%D(3,AtD)
!-------------------------------------------------------------------------------------------
              ! Differencing !!!
!
!              iDir=3
!              iDir=6
!              iDir=9
              iDir=12
!              Tmp1=ACAtmInfo%Atm1Z
!              Tmp1=ACAtmInfo%Atm2Z
!              Tmp1=BDAtmInfo%Atm1Z
              Tmp1=BDAtmInfo%Atm2Z
              DO iDiff=1,2
                 IF(iDiff==1)THEN
!                    ACAtmInfo%Atm1Z=Tmp1-1D-4
!                    ACAtmInfo%Atm2Z=Tmp1-1D-4
!                    BDAtmInfo%Atm1Z=Tmp1-1D-4
                    BDAtmInfo%Atm2Z=Tmp1-1D-4
                 ELSE
!                    ACAtmInfo%Atm1Z=Tmp1+1D-4
!                    ACAtmInfo%Atm2Z=Tmp1+1D-4
!                    BDAtmInfo%Atm1Z=Tmp1+1D-4
                    BDAtmInfo%Atm2Z=Tmp1+1D-4
                 ENDIF
                 CALL GetAtomPair(ACAtmInfo,ACAtmPair,BS,CS_OUT%CellCarts%D(1,1))
                 CALL GetAtomPair(BDAtmInfo,BDAtmPair,BS,CS_OUT%CellCarts%D(1,1))
!                 ACAtmInfo%Atm1Z=Tmp1
!                 ACAtmInfo%Atm2Z=Tmp1
!                 BDAtmInfo%Atm1Z=Tmp1
                 BDAtmInfo%Atm2Z=Tmp1

                 CALL DBL_VECT_EQ_DBL_SCLR(NBFA*NBFB*NBFC*NBFD,C(1),0.0d0)
                 CFAC=0
                 OffSet%A=1
                 DO CFA=1,BS%NCFnc%I(KA)
                    AALen=BS%LStop%I(CFA,KA)-BS%LStrt%I(CFA,KA)+1
                    OffSet%C=1
                    DO CFC=1,BS%NCFnc%I(KC)
                       CCLen=BS%LStop%I(CFC,KC)-BS%LStrt%I(CFC,KC)+1
                       CFBD=0
                       OffSet%B=1
                       CFAC=CFAC+1
                       DO CFB=1,BS%NCFnc%I(KB)
                          BBLen=BS%LStop%I(CFB,KB)-BS%LStrt%I(CFB,KB)+1
                          OffSet%D=1
                          DO CFD=1,BS%NCFnc%I(KD)
                             CFBD=CFBD+1
                             DDLen=BS%LStop%I(CFD,KD)-BS%LStrt%I(CFD,KD)+1
                             IntType=ACAtmPair(CFAC)%SP%IntType*10000+BDAtmPair(CFBD)%SP%IntType
!                             WRITE(*,*)' IntType = ',ACAtmPair(CFAC)%SP%IntType,BDAtmPair(CFBD)%SP%IntType
                             INCLUDE 'ERIInclude.Inc'
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
              CALL GetAtomPair(ACAtmInfo,ACAtmPair,BS,CS_OUT%CellCarts%D(1,1))
              CALL GetAtomPair(BDAtmInfo,BDAtmPair,BS,CS_OUT%CellCarts%D(1,1))
              DO I=1,12
                 CALL DBL_VECT_EQ_DBL_SCLR(NBFA*NBFB*NBFC*NBFD,dC(1,I),0.0d0)
              ENDDO
              NIntBlk=MaxBK ! NBFA*NBFB*NBFC*NBFD
              CFAC=0
              OffSet%A=1
              DO CFA=1,BS%NCFnc%I(KA)
                 AALen=BS%LStop%I(CFA,KA)-BS%LStrt%I(CFA,KA)+1
                 OffSet%C=1
                 DO CFC=1,BS%NCFnc%I(KC)
                    CCLen=BS%LStop%I(CFC,KC)-BS%LStrt%I(CFC,KC)+1
                    CFBD=0
                    OffSet%B=1
                    CFAC=CFAC+1
                    DO CFB=1,BS%NCFnc%I(KB)
                       BBLen=BS%LStop%I(CFB,KB)-BS%LStrt%I(CFB,KB)+1
                       OffSet%D=1
                       DO CFD=1,BS%NCFnc%I(KD)
                          CFBD=CFBD+1
                          DDLen=BS%LStop%I(CFD,KD)-BS%LStrt%I(CFD,KD)+1
                          IntType=ACAtmPair(CFAC)%SP%IntType*10000+BDAtmPair(CFBD)%SP%IntType
!                          WRITE(*,*)' IntType = ',ACAtmPair(CFAC)%SP%IntType,BDAtmPair(CFBD)%SP%IntType
                          INCLUDE 'dERIInclude.Inc'
                          OffSet%D=OffSet%D+DDLen
                       ENDDO ! End blkfunc on D
                       OffSet%B=OffSet%B+BBLen
                    ENDDO ! End blkfunc on B
                    OffSet%C=OffSet%C+CCLen
                 ENDDO ! End blkfunc on C
                 OffSet%A=OffSet%A+AALen                 
              ENDDO ! End blkfunc on A
              CALL CompareBlock(AtAOff,AtBOff,AtCOff,AtDOff,NBFA,NBFB,NBFC,NBFD,C(1),dC(1,iDir))
              AtDOff=AtDOff+NBFD
           ENDDO ! D
           AtBOff=AtBOff+NBFB
        ENDDO ! B
        AtCOff=AtCOff+NBFC
     ENDDO ! C
     AtAOff=AtAOff+NBFA
  ENDDO ! A
END PROGRAM N4KTest

