!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe and  C. J. Tymczak
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!    COMPUTE THE Kintic Energy Matrix
!
PROGRAM GradT
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE AtomPairs
  USE GradTBlock
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
!-------------------------------------------------------- 
!     FUNCTION TBlok(BS,MD,Pair)
!        USE DerivedTypes
!        USE GlobalScalars
!        IMPLICIT NONE
!        TYPE(BSET),                INTENT(IN)    :: BS
!        TYPE(DBL_RNK4),            INTENT(IN)    :: MD
!        TYPE(Pair),                INTENT(IN)    :: Pair
!        REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB)  :: TVck 
!     END SUBROUTINE TBlok
!--------------------------------------------------------
! 
#ifdef PARALLEL
  TYPE(DBCSR)         :: P
#else
  TYPE(BCSR)          :: P
#endif
#ifdef PERIODIC 
  LOGICAL             :: NotMakeBlock,InitBlock
  INTEGER             :: NC,NCells,NX,NY,NZ,I
  REAL(DOUBLE)        :: Bx,By,Bz
#endif
  TYPE(AtomPair)      :: Pair
!
  TYPE(BSET)          :: BS
  TYPE(CRDS)          :: GM
  TYPE(DBL_RNK4)      :: MD
  TYPE(DBL_RNK2)      :: TrPT
!
  TYPE(ARGMT)         :: Args
  INTEGER             :: PP,R,AtA,AtB,NN,IStrtP,IStopP,LP,JP,MB,MA
                         
  CHARACTER(LEN=5),PARAMETER :: Prog='GradT'
!--------------------------------------- 
! Start up macro
!
  CALL StartUp(Args,Prog)
!----------------------------------------------
! Get basis set and geometry
!
  CALL Get(BS)
  CALL Get(GM)
!---------------------------------------------- 
! Allocations 
!
  CALL New(MD,(/3,BS%NASym+4,BS%NASym+4,2*BS%NASym+4/),(/1,-2,-2,0/))
  CALL New(P)
  CALL Get(P,TrixFile('D',Args,1))
  CALL New(TrPT,(/3,NAtoms/))
!-----------------------------------------------
! Set the Distance Overlap Threshhold
!-----------------------------------------------
  CALL SetOverLapDist(Thresholds%Dist)
!-----------------------------------------------
! Initialize the matrix and associated indecies
!
  P%RowPt%I(1)=1
#ifdef PERIODIC
!-----------------------------------------------
! Calculate the Number of Cells
!
  CALL SetCellNumber(GM)
#endif
!-----------------------------------------------
! Main loops
!
  TrPT%D=Zero
#ifdef PARALLEL
  DO AtA=Beg%I(MyId),End%I(MyId)
#else
  DO AtA=1,NAtoms
#endif
     MA=BSiz%I(AtA)
     IStrtP=P%RowPt%I(AtA)
     IStopP=P%RowPt%I(AtA+1)-1
     DO JP=IStrtP,IStopP
        AtB=P%ColPt%I(JP)
        PP=P%BlkPt%I(JP)
        MB=BSiz%I(AtB)
        IF(SetAtomPair(GM,BS,AtB,AtA,Pair)) THEN
           TrPT%D(:,AtA)=TrPT%D(:,AtA) + &
                GradTBlok(BS,MD,AtB,AtA,AtA,P%MTrix%D(PP:PP+MA*MB-1),Pair)
           IF (AtA.NE.AtB) THEN
              TrPT%D(:,AtB)=TrPT%D(:,AtB) + &
                   GradTBlok(BS,MD,AtB,AtA,AtB,P%MTrix%D(PP:PP+MA*MB-1),Pair)
           ENDIF
        ENDIF
     ENDDO
  ENDDO
!------------------------------------------------------------
! Tidy up
! 
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(MD)
  CALL ShutDown(Prog)
!
END PROGRAM GradT
