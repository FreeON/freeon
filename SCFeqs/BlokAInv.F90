  !    COMPUTE THE INCOMPLETE INVERSE CHOLESKY FACTOR 
  !    OF THE OVERLAP MATRIX Z=S^(-L) 
  !    Author: Matt Challacombe
  !-----------------------------------------------------------------
  PROGRAM BlokAInv
    USE DerivedTypes
    USE GlobalScalars
    USE GlobalCharacters
    USE DenMatMethods
    USE InOut
    USE PrettyPrint
    USE MemMan
    USE Parse
    USE Macros
    USE LinAlg
    USE AInv
#ifdef NAG
    USE F90_UNIX_ENV
#endif
    IMPLICIT NONE
    TYPE(BCSR)          :: A,Z,Zt,DiagD
    !#ifdef SPATIAL_THRESHOLDING
    TYPE(BCSR)          :: T1,T2
    !#endif 
#ifdef FIND_CONDA
    TYPE(DBL_RNK2)      :: B,C
#endif
    TYPE(BSET)          :: BS
    TYPE(CRDS)          :: GM
    TYPE(INT_VECT)      :: AiFlg,ZiFlg,ColPt,BlkPt
    TYPE(ARGMT)         :: Args
    INTEGER             :: I,J,Q,R,IDex,JDex,ZDex,ZBlk,NIJ, &
         n,ni,msiz,strtai,stopai,strtaj,stopaj, &
         strtzi,stopzi,nj,strtzj,stopzj,jcol,k,kdex, &
         aiblk,ajblk,zjblk,m,ziblk,icol,zrowpt,zcolpt, & 
         zblkpt,NewBloks,EndBloks,IRow,JRow,ZBlksPreFilter,ZBlksPostFilter
    TYPE(DBL_VECT)      :: Blk1,Blk2
    TYPE(DBL_RNK2)      :: P,DA
    REAL(DOUBLE)        :: Op,Mx0,B2Norm

    REAL(DOUBLE)        :: AInvDistanceThresh,IRowX,IRowY,IRowZ
    TYPE(AtomPair)      :: Pair
    REAL(DOUBLE), &
         EXTERNAL :: DDOT

    LOGICAL             :: TEST_AINV
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
    CHARACTER(LEN=8),&
         PARAMETER  :: Prog='BlokAInv'
!-----------------------------------------------------------------------------------------------------------
    ! Start up macro
    CALL StartUp(Args,Prog)
    ! Get basis set and geometry
    CALL Get(BS,Tag_O=CurBase)
    CALL Get(GM,Tag_O=CurGeom)
    TEST_AINV=.TRUE.

    CALL SussTrix('AINVThreshold',Prog)

#ifdef SPATIAL_THRESHOLDING
    IF(GM%Ordrd==SFC_HILBERT.OR.GM%Ordrd==SFC_PEANO)THEN
       CALL OpenASCII(InpFile,Inp)
       IF(OptDblQ(Inp,Prog,AInvDistanceThresh))THEN
          WRITE(*,*)'From input, AInvDistanceThresh = ',AInvDistanceThresh
       ELSE
          AInvDistanceThresh=1.D3
          WRITE(*,*)' Using default AInvDistanceThreshold=',AInvDistanceThresh
       ENDIF
       IF(OptKeyQ(Inp,Prog,'NoTest'))THEN
          TEST_AINV=.FALSE.
       ENDIF
       CLOSE(UNIT=Inp)
    ELSE
       AInvDistanceThresh=1.D24
       TEST_AINV=.FALSE.
       WRITE(*,*)' Using O(N^2) AInv due to non-local ordering. Try HOrder or ZOrder '
    ENDIF
#endif
    ! Allocations 
    CALL New(A)
    CALL Get(A,TrixFile('S',Args))
    ! Blocked AINV
    CALL BlockedAInv(A,Thresholds%Trix,GM,AInvDistanceThresh,Z,Zt,PerfMon)

    IF(PrintFlags%Key>=DEBUG_MEDIUM)THEN
       CALL PPrint(PerfMon,Prog)  
       CALL PPrint(PerfMon,Prog,Unit_O=6)  
    ENDIF
    ! Consistency check
    IF(TEST_AINV)THEN
       CALL New(T1)
       CALL New(T2)
       CALL Multiply(Zt,A,T1)
       CALL Multiply(T1,Z,T2)
       !
       CALL SetToI(T1)
       CALL Multiply(T1,-One)
       CALL Add(T1,T2,A)
       Mx0=Max(A)
       !
       IF(AInvDistanceThresh/=Zero)THEN     
          Mssg='Max(Z^t.A.Z-I)='//TRIM(DblToShrtChar(Mx0))                 &
               //', DistThrsh='//TRIM(DblToShrtChar(AInvDistanceThresh))  &
               //', TrixThrsh='//TRIM(DblToShrtChar(Thresholds%Trix))
       ELSE
          Mssg='Max(Z^t.A.Z-I)='//TRIM(DblToShrtChar(Mx0))                 &
               //', TrixThrsh='//TRIM(DblToShrtChar(Thresholds%Trix))
       ENDIF
       IF(Mx0>1.D2*Thresholds%Trix)THEN
          CALL Warn('In BlokAInv, failed test: '//TRIM(Mssg))
       ENDIF
       IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
          Mssg=ProcessName(Prog)//TRIM(Mssg)
          WRITE(*,*)TRIM(Mssg)
          CALL OpenASCII(OutFile,Out)
          CALL PrintProtectL(Out)
          WRITE(Out,*)TRIM(Mssg)
          CALL PrintProtectR(Out)
          CLOSE(UNIT=Out,STATUS='KEEP')
       ENDIF
    ENDIF
    !  Put Z and ZT to disk
    CALL Put(Z,TrixFile('Z',Args))
    CALL Put(Zt,TrixFile('ZT',Args))
    !  Debug
    CALL PChkSum(Z,'Z',Proc_O=Prog)
    CALL PPrint(Z,'Z')
    CALL Plot(Z,'Z')
    !  Tidy up
    CALL Delete(Z)
    CALL Delete(ZT)
    CALL ShutDown(Prog) 
  END PROGRAM BlokAInv







