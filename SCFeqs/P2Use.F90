!    COMPUTES AN ORTHOGONAL GUESS DENSITY MATRIX EITHER FROM
!    A SUPERPOSITION OF DIAGONAL, ATOMIC LEWIS STRUCTURE BLOCKS 
!    OR FROM A ORTHOGONAL DENSITY MATRIX COORESPONDING TO A DIFFERENT 
!    (HOPEFULLY CLOSE) GEOMETRY OR LEVEL OF ACCURACY 
!    Author: Matt Challacombe
!-------------------------------------------------------------------------
PROGRAM P2Use
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE PBlokGuess
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
  TYPE(ARGMT)                   :: Args
  TYPE(BSET)                    :: BS
  TYPE(CRDS)                    :: GM
#ifdef PARALLEL
  TYPE(DBCSR) &
#else
  TYPE(BCSR)  & 
#endif
                                :: S,P,X,T0,T1,T2
  TYPE(DBL_RNK2)                :: BlkP
  REAL(DOUBLE)                  :: Scale,TrP,Fact,ECount, &
       DeltaP,OldDeltaP,DensityDev
  INTEGER                       :: I,J,JP,AtA,Q,R,T,KA,NBFA, &
       NPur,PcntPNon0,OldFileID
  CHARACTER(LEN=2)              :: Cycl
  LOGICAL                       :: Present
  CHARACTER(LEN=DEFAULT_CHR_LEN):: Mssg,BName
  CHARACTER(LEN=5),PARAMETER    :: Prog='P2Use'
!------------------------------------------------------------------------------- 
  ! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  Cycl=IntToChar(Args%I%I(1))
  ! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
  IF(SCFActn=='Extrapolate')THEN
     CALL Halt(' Extrapolation turned off, need non-orthogonal SP2 or TS4... ')
  ELSEIF(SCFActn=='Restart')THEN  
     ! Close current group and HDF
     CALL CloseHDFGroup(H5GroupID)
     CALL CloseHDF(HDFFileID)
     ! Open old group and HDF
     OldFileID=OpenHDF(Restart)
     HDF_CurrentID=OpenHDFGroup(OldFileID,"Clone #"//TRIM(IntToChar(MyClone)))
     ! Find the current density matrix
     CALL Get(P,'CurrentDM',CheckPoint_O=.TRUE.)
     ! Close it up 
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(OldFileID)
     ! Reopen current group and HDF
     HDFFileID=OpenHDF(H5File)
     H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
     HDF_CurrentID=H5GroupID
  ELSEIF(SCFActn=='Project')THEN
     ! Allocations 
     CALL New(P)
     CALL New(T0)
     CALL New(T1)
     CALL New(T2)
     ! Get previous geometries orthogonal density matrix 
     CALL Get(P,TrixFile('OrthoD',Args,-1))     
#ifdef PARALLEL
     IF(MyId==ROOT)THEN
#endif
        INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
#ifdef PARALLEL
     ENDIF
     CALL BCast(Present)
#endif
     IF(Present)THEN     
        CALL Get(T1,TrixFile('X',Args))   ! T1=S_new^(-1/2)
        CALL Multiply(T1,P,T2)            ! T2=S_new^(-1/2).P_old
        CALL Multiply(T2,T1,T0)           ! P_new_AO=S_new^(-1/2).P_old.S_new^(-1/2)
        CALL Filter(P,T0)                 ! T1=Filter[P_new_AO]
     ELSE
        CALL Get(T1,TrixFile('Z',Args))   ! T1=Z_new
        CALL Multiply(T1,P,T2)            ! T2=Z.P_old
        CALL Get(T1,TrixFile('ZT',Args))  ! T1=Z^T
        CALL Multiply(T2,T1,T0)           ! P_new_AO=Z.P_old.Z^T
        CALL Filter(P,T0)                 ! T1=Filter[P_new_AO]
     ENDIF
     CALL Delete(T0)
     CALL Delete(T1)
     CALL Delete(T2)
  ELSEIF(SCFActn=='DensitySuperposition')THEN
     ! Make a diagonal guess
     CALL Get(BName,'bsetname',CurBase)
     IF(INDEX(BName,'STO')/=0)THEN
        ! Compute a diagonal guess as the superposition of atomic lewis 
        ! structure occupancies--works only for minimal (STO) basis sets
        CALL New(BlkP,(/MaxBlkSize**2,NAtoms/))
        DO I=1,NAtoms
           CALL FillPBlok(BSiz%I(I),INT(GM%AtNum%D(I)),BlkP%D(:,I))
        ENDDO
        CALL SetToI(P,BlkP)
!        CALL SetEq(P,P2)
        ! Check for the correct elctron count
        TrP=Trace(P)
        IF(ABS(TrP-DBLE(NEl/Two))>1.D-10) &
             CALL Warn(' In P2Use, TrP = '//TRIM(DblToChar(TrP)))
        CALL Delete(BlkP)
#ifdef NEW_NOTVERYGOOD_BLOCKDIAGONAL
        CALL Get(S,TrixFile('S',Args,Stats_O=Current))
        ! Set global workspace for FunkOnSqMat
        CALL SetDSYEVWork(MaxBlkSize)
        CALL New(X)
        T=1; R=1; X%RowPt%I(1)=1
        DO I=1,NAtoms
           DO JP=S%RowPt%I(I),S%RowPt%I(I+1)-1
              J=S%ColPt%I(JP)
              IF(I==J)THEN
                 CALL FunkOnSqMat(BSiz%I(I),InvSqRt,S%MTrix%D(S%BlkPt%I(JP)),X%MTrix%D(R))
                                 !,PrintValues_O=.TRUE.,Unit_O=6)
                 X%ColPt%I(T)=I
                 X%BlkPt%I(T)=R
                 R=R+BSiz%I(I)**2
                 T=T+1 
                 X%RowPt%I(I+1)=T        
              ENDIF
           ENDDO
        ENDDO
        X%NBlks=T-1
        X%NNon0=R-1
        CALL New(T1)
        CALL Multiply(X,P,T1)   ! T1=S^(-1/2).DiagonalGuess
        CALL Multiply(T1,X,P)   ! P=S^(-1/2).DiagonalGuess.S^(-1/2)
        CALL Delete(X)
        CALL Delete(T1)
#else
#ifdef PARALLEL
        IF(MyId==ROOT)THEN
#endif
           INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
#ifdef PARALLEL
        ENDIF
        CALL BCast(Present)
#endif
        IF(Present)THEN     
           CALL Get(T1,TrixFile('X',Args))   ! T1=S_new^(-1/2)
           CALL Multiply(T1,P,T2)            ! T2=S_new^(-1/2).P_old
           CALL Multiply(T2,T1,T0)           ! P_new_AO=S_new^(-1/2).P_old.S_new^(-1/2)
           CALL Filter(P,T0)                 ! T1=Filter[P_new_AO]
        ELSE
           CALL Get(T1,TrixFile('Z',Args))   ! T1=Z_new
           CALL Multiply(T1,P,T2)            ! T2=Z.P_old
           CALL Get(T1,TrixFile('ZT',Args))  ! T1=Z^T
           CALL Multiply(T2,T1,T0)           ! P_new_AO=Z.P_old.Z^T
           CALL Filter(P,T0)                 ! T1=Filter[P_new_AO]
        ENDIF
        CALL Delete(T0)
        CALL Delete(T1)
        CALL Delete(T2)
#endif
     ENDIF
  ELSE
     ! Guess == Core
     CALL New(BlkP,(/MaxBlkSize**2,NAtoms/))
     DO I=1,NAtoms
        BlkP%D(:,I)=Zero
     ENDDO
     CALL SetToI(P,BlkP)
     CALL Delete(BlkP)
  ENDIF
  ! IO for the non-orthogonal P 
  CALL Put(P,TrixFile('D',Args,0))
  CALL PChkSum(P,'P['//TRIM(Cycl)//']',Prog)
  CALL PPrint( P,'P['//TRIM(Cycl)//']')
  CALL Plot(   P,'P_'//TRIM(Cycl))
  ! Tidy up ...
  CALL Delete(GM)
  CALL Delete(BS)
  CALL Delete(P)
  CALL ShutDown(Prog)   
END PROGRAM P2Use

      
