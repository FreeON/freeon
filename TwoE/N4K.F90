SUBROUTINE SetBigIBlock(AtAOff,AtBOff,AtCOff,AtDOff,NA,NB,NC,ND,I)
  USE DerivedTypes
  INTEGER :: AtAOff,AtBOff,AtCOff,AtDOff,NA,NB,NC,ND,IA,IB,IC,ID
  REAL(DOUBLE) :: I(1:NA,1:NC,1:NB,1:ND)
  DO IA=1,NA
     DO IC=1,NC
        DO IB=1,NB
           DO ID=1,ND
              IF(ABS(I(IA,IC,IB,ID))>1D-12)THEN 
!              WRITE(*,101)IA+AtAOff,IC+AtCOff,IB+AtBOff,ID+AtDOff,I(IA,IC,IB,ID)
              WRITE(33,301)IA+AtAOff,IC+AtCOff,IB+AtBOff,ID+AtDOff, &
                          FRACTION(I(IA,IC,IB,ID)),EXPONENT(I(IA,IC,IB,ID))
              !WRITE(33,301)IA+AtAOff,IC+AtCOff,IB+AtBOff,ID+AtDOff,I(IA,IC,IB,ID)
              !WRITE(iOut,301)I,J,K,L,FRACTION(Int),EXPONENT(Int)
              ENDIF

           ENDDO
        ENDDO
     ENDDO
  ENDDO
101 FORMAT(' (',I3,',',I3,'|',I3,',',I3,') = ',F20.16)
201 FORMAT(' (',I3,',',I3,'|',I3,',',I3,') = ',D22.16)
301 FORMAT(1x,'Int2[[',I3,',',I3,',',I3,',',I3,']] = ',F19.16,'*2^(',I3,');')
  
END SUBROUTINE SetBigIBlock

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
    TYPE(CRDS)              :: GMc
    TYPE(BSET)              :: BS
    !  TYPE(CellSet)           :: CS_OUT
    !  TYPE(ANode2) , POINTER                :: AtAList,AtAListTmp,NodeA
    TYPE(AtomInfo)                        :: ACAtmInfo,BDAtmInfo
    INTEGER                               :: AtA,AtC,AtB,AtD,KA,KC,KB,KD, &
         CFA,CFC,CFB,CFD,iFAC,iFBD,NBFA,NBFC,NBFB,NBFD,AtAOff,AtBOff,AtCOff,AtDOff, & 
         AALen,BBLen,CCLen,DDLen,NCFA,NCFB,NCFC,NCFD
    INTEGER                               :: NCell,IntType,LocNInt,Off
    REAL(DOUBLE)                          :: RInt,AC2,BD2,NInts
    TYPE(AtomPr), DIMENSION(100)          :: ACAtmPair,BDAtmPair
    TYPE(ONX2OffSt) :: OffSet

    REAL(DOUBLE) , DIMENSION(50)          :: RIntCell  ! this should be declared somewhere
    REAL(DOUBLE) , DIMENSION(30**4)       :: C         ! this should be declared somewhere
    REAL(DOUBLE) , DIMENSION(50,50,50,50) :: I
    REAL(DOUBLE) , EXTERNAL               :: DGetAbsMax
    TYPE(ARGMT)                    :: Args
    CHARACTER(LEN=3) :: Prog="N4K"

    CALL StartUp(Args,Prog)

    CALL Get(BS,Tag_O=CurBase)
    CALL Get(GMc,Tag_O=CurGeom)
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
          ACAtmInfo%Atm12X=ACAtmInfo%Atm1X-ACAtmInfo%Atm2X
          ACAtmInfo%Atm12Y=ACAtmInfo%Atm1Y-ACAtmInfo%Atm2Y
          ACAtmInfo%Atm12Z=ACAtmInfo%Atm1Z-ACAtmInfo%Atm2Z
          AC2=(ACAtmInfo%Atm12X)**2+(ACAtmInfo%Atm12Y)**2+(ACAtmInfo%Atm12Z)**2
          CALL GetAtomPair2(ACAtmInfo,ACAtmPair,BS,CS_OUT%CellCarts%D(1,1))
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
                BDAtmInfo%Atm12X=BDAtmInfo%Atm1X-BDAtmInfo%Atm2X
                BDAtmInfo%Atm12Y=BDAtmInfo%Atm1Y-BDAtmInfo%Atm2Y
                BDAtmInfo%Atm12Z=BDAtmInfo%Atm1Z-BDAtmInfo%Atm2Z
                BD2=(BDAtmInfo%Atm12X)**2+(BDAtmInfo%Atm12Y)**2+(BDAtmInfo%Atm12Z)**2
                CALL GetAtomPair2(BDAtmInfo,BDAtmPair,BS,CS_OUT%CellCarts%D(1,1))
                CALL DBL_VECT_EQ_DBL_SCLR(NBFA*NBFB*NBFC*NBFD,C(1),0.0d0)
!                CALL DBL_VECT_EQ_DBL_SCLR(NBFA*NBFB*NBFC*NBFD,C(1),BIG_DBL)
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

#ifdef LDSJFLSDFJSLDFJ
   SELECT CASE(IntType)
   CASE(3030303)
     CALL IntB3030303(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
                  OffSet%A  ,1             , & 
                  OffSet%C-1,NBFA*NBFB     , & 
                  OffSet%B-1,NBFA          , & 
                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))

!     CALL Int3333(ACAtmPair(iFAC)%SP%Cst(1,1),ACAtmPair(iFAC)%SP%L, & 
!                  BDAtmPair(iFBD)%SP%Cst(1,1),BDAtmPair(iFBD)%SP%L, & 
!                  ACAtmPair(iFAC)%SP%AtmInfo,BDAtmPair(iFBD)%SP%AtmInfo, & 
!                  OffSet%A  ,1             , & 
!                  OffSet%C-1,NBFA*NBFB     , & 
!                  OffSet%B-1,NBFA          , & 
!                  OffSet%D-1,NBFA*NBFB*NBFC,GMc%PBC,C(1))
     LocNInt=81
   END SELECT
#endif

                            INCLUDE 'ERIInterfaceB.Inc'

                            CALL ShellPrint(NBFA,NBFB,NBFC,NBFD,AALen,BBLen,CCLen,DDLen,  &
                                            AtAOff,AtBOff,AtCOff,AtDOff,                  &
                                            OffSet%A,OffSet%B,OffSet%C,OffSet%D,IntType,C(1))

                            OffSet%D=OffSet%D+DDLen
                         ENDDO ! End blkfunc on D
                         OffSet%B=OffSet%B+BBLen
                      ENDDO ! End blkfunc on B
                      OffSet%C=OffSet%C+CCLen
                   ENDDO ! End blkfunc on C
                   OffSet%A=OffSet%A+AALen                 !
                ENDDO ! End blkfunc on A
                CALL SetBigIBlock(AtAOff,AtBOff,AtCOff,AtDOff,NBFA,NBFB,NBFC,NBFD,C(1))
                AtDOff=AtDOff+NBFD
             ENDDO ! D
             AtBOff=AtBOff+NBFB
          ENDDO ! B
          AtCOff=AtCOff+NBFC
       ENDDO ! C
       AtAOff=AtAOff+NBFA
    ENDDO ! A
  END PROGRAM N4KTest
