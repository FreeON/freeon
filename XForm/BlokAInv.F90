!    COMPUTE THE INCOMPLETE INVERSE CHOLESKY FACTOR 
!    OF THE OVERLAP MATRIX Z=S^(-L) USING BENZI AND TUMAS 
!    BLOCKED AINV
!    Author: Matt Challacombe
!-----------------------------------------------------------------
PROGRAM BlokAInv
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE MatFunk
  USE AtomPairs
#ifdef NAG
   USE F90_UNIX_ENV
#endif
  IMPLICIT NONE
  TYPE(BCSR)          :: A,Z,Zt,DiagD
#ifdef SPATIAL_THRESHOLDING
  TYPE(BCSR)          :: T1,T2
#endif 
#ifdef EXTREME_DEBUG
  TYPE(DBL_RNK2)      :: DnsD1,DnsD2,DnsD3,DnsZ1,DnsZ2,DnsZ3,B,C
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

#ifdef EXTREME_DEBUG
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: AInvFile
#endif
!-----------------------------------------------------------------------------------------------------------
! Start up macro
!
  CALL StartUp(Args,Prog)

#ifdef EXTREME_DEBUG

#endif
!
! Get basis set and geometry
!
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
!
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
     ELSE
        TEST_AINV=.TRUE.
        WRITE(*,*)' Testing AInv ... '
     ENDIF
     CLOSE(UNIT=Inp)
  ELSE
     AInvDistanceThresh=1.D24
     TEST_AINV=.FALSE.
     WRITE(*,*)' Using O(N^2) AInv due to non-local ordering. Try HOrder or ZOrder '
  ENDIF
  
#endif
!
! Allocations 
!
  CALL New(A)
  CALL Get(A,TrixFile('S',Args))
#ifdef FIND_CONDA
  CALL New(B,(/NBasF,NBasF/))
  CALL New(C,(/NBasF,NBasF/))
  CALL SetEq(B,A)
  CALL SetDSYEVWork(NBasF)
  CALL FunkOnSqMat(NBasF,Inverse,B%D,C%D,PrintCond_O=.TRUE.)
  CALL Delete(B)
  CALL Delete(C)
  CALL UnSetDSYEVWork()
#endif
  CALL New(Z)
!
!  PrintFlags%Mat=PLOT_MATRICES
#ifdef USE_METIS
  CALL MetisReorder(A)
#endif
! Set global workspace for FunkOnSqMat
  CALL SetDSYEVWork(MaxBlkSize)
! Allocate intermediate blocks 
  CALL New(Blk1,MaxBlkSize*MaxBlkSize)
  CALL New(Blk2,MaxBlkSize*MaxBlkSize)
! Allocate diagonal "pivot" blocks
  CALL New(P,(/MaxBlkSize*MaxBlkSize,NAtoms/))
! Allocate coloumn flags
  CALL New(AiFlg,NAtoms); AiFlg%I=0
  CALL New(ZiFlg,NAtoms); ZiFlg%I=0
! Allocate temporaries for symbolic Z_I
  CALL New(ColPt,MaxBlks)
  CALL New(BlkPt,MaxBlks)
!
! Start with the identity; Z=I
!
  CALL SetToI(Z)
! Index for new Z bloks
  ZBlk=Z%NNon0+1
!
! Main outer loop down rows 
!
  DO IRow=1,NAtoms
!
!    Set IRow coordinates for distance based screening
     IRowX=GM%Carts%D(1,IRow)
     IRowY=GM%Carts%D(2,IRow)
     IRowZ=GM%Carts%D(3,IRow)
!
     NI=BSiz%I(IRow)
     StrtAI=A%RowPt%I(IRow);StopAI=A%RowPt%I(IRow+1)-1
     StrtZI=Z%RowPt%I(IRow);StopZI=Z%RowPt%I(IRow+1)-1
     DO I=StrtAI,StopAI; AiFlg%I(A%ColPt%I(I))=I; ENDDO
!    Store symbolic structure of Z_I
     ZDex=StrtZI
     DO I=StrtZI,StopZI
        IDex=Z%ColPt%I(I)
        ZiFlg%I(IDex)=ZDex
        ColPt%I(ZDex)=IDex
        BlkPt%I(ZDex)=Z%BlkPt%I(I)
        ZDex=ZDex+1
     ENDDO
!
!    Inner loop down rows ...
!
     DO JRow=1,IRow-1
!
#ifdef SPATIAL_THRESHOLDING
        IF(((IRowX-GM%Carts%D(1,JRow))**2+ &
            (IRowY-GM%Carts%D(2,JRow))**2+ &
            (IRowZ-GM%Carts%D(3,JRow))**2)<AInvDistanceThresh)THEN
#endif
           NJ=BSiz%I(JRow)
           StrtAJ=A%RowPt%I(JRow);StopAJ=A%RowPt%I(JRow+1)-1
           StrtZJ=Z%RowPt%I(JRow);StopZJ=Z%RowPt%I(JRow+1)-1
!
!          Blk1=P^(j-1)_i=[A^t_j].[Z^(j-1)_i]; Going down rows over N: (NxNJ)^T.(NxNI)
!
           NIJ=NI*NJ
           Blk1%D(1:NIJ)=Zero
           DO J=StrtAJ,StopAJ
              JDex=A%ColPt%I(J)
#ifdef SPATIAL_THRESHOLDING
              IF(((IRowX-GM%Carts%D(1,JDex))**2+ &
                  (IRowY-GM%Carts%D(2,JDex))**2+ &
                  (IRowZ-GM%Carts%D(3,JDex))**2)<AInvDistanceThresh)THEN
#endif
!
                 IDex=ZiFlg%I(JDex)
                 IF(IDex/=0)THEN                 
                    ZiBlk=BlkPt%I(IDex) 
                    AjBlk=A%BlkPt%I(J)
                    M=BSiz%I(JDex)
                    CALL DGEMM_NN(NJ,M,NI,One,A%MTrix%D(AjBlk),Z%MTrix%D(ZiBlk),Blk1%D)
                    PerfMon%FLOP=PerfMon%FLOP+DBLE(NIJ*M)
                 ENDIF
#ifdef SPATIAL_THRESHOLDING
              ENDIF
#endif
           ENDDO
!
!          Blk2=[P^(j-1)_j]^(-1).[P^(j-1)_i]
!
           CALL DGEMM_NNc(NJ,NJ,NI,One,Zero,P%D(:,JRow),Blk1%D,Blk2%D)
           PerfMon%FLOP=PerfMon%FLOP+DBLE(NIJ*NJ)
! 
!          Check the magintude of Blk2.  Update Z_I only if Blk2 is "large" enough.       
!        
           B2Norm=SQRT(DDOT(NI*NJ,Blk2%D,1,Blk2%D,1))
           PerfMon%FLOP=PerfMon%FLOP+DBLE(NIJ)
!        
           IF(B2Norm>Thresholds%Trix*1.D-1)THEN
!
!             Z^j_i=Z^(j-1)_i-[Z^(j-1)_j].{[P^(j-1)_j]^(-1).[P^(j-1)_i]}
!             Update going down rows:(NxNI)=(NxNI)+(NxNJ).(NJxNI)
!
              DO JDex=StrtZJ,StopZJ
                 JCol=ColPt%I(JDex) 
#ifdef SPATIAL_THRESHOLDING
                 IF(((IRowX-GM%Carts%D(1,JCol))**2+ &
                     (IRowY-GM%Carts%D(2,JCol))**2+ &
                     (IRowZ-GM%Carts%D(3,JCol))**2)<AInvDistanceThresh)THEN
#endif
                    ZjBlk=BlkPt%I(JDex) 
                    IDex =ZiFlg%I(JCol)
                    M=BSiz%I(JCol)
                    PerfMon%FLOP=PerfMon%FLOP+DBLE(NIJ*M)
                    IF(IDex/=0)THEN                  
                      ZiBlk=BlkPt%I(IDex) 
                      CALL DGEMM_NNc(M,NJ,NI,-One,One,Z%MTrix%D(ZjBlk),Blk2%D,Z%MTrix%D(ZiBlk))
                    ELSE
                      ZiBlk=ZBlk
                      CALL DGEMM_NNc(M,NJ,NI,-One,Zero,Z%MTrix%D(ZjBlk),Blk2%D,Z%MTrix%D(ZiBlk))
                      ZiFlg%I(JCol)=ZDex
                      ColPt%I(ZDex)=JCol
                      BlkPt%I(ZDex)=ZiBlk
                      ZDex=ZDex+1 
                      ZBlk=ZBlk+M*NI
                    ENDIF
#ifdef SPATIAL_THRESHOLDING
                 ENDIF
#endif
              ENDDO
#ifdef SPATIAL_THRESHOLDING
           ENDIF
#endif
!
       ENDIF
!
     ENDDO ! end inner loop over JRow
!
!    Reup symbolic structure of Z
!
     IF(IRow>1)THEN
        NewBloks=ZDex-StopZI-1
        EndBloks=Z%NBlks-StopZJ-1
        Z%RowPt%I(IRow+1:NAtoms+1)=Z%RowPt%I(IRow+1:NAtoms+1)+NewBloks
        ColPt%I(ZDex:ZDex+EndBloks-1)=Z%ColPt%I(StopZI+1:Z%NBlks)
        BlkPt%I(ZDex:ZDex+EndBloks-1)=Z%BlkPt%I(StopZI+1:Z%NBlks)
        Z%NBlks=Z%NBlks+NewBloks
        Z%NNon0=ZBlk
        Z%ColPt%I(StrtZI:Z%NBlks)=ColPt%I(StrtZI:Z%NBlks)
        Z%BlkPt%I(StrtZI:Z%NBlks)=BlkPt%I(StrtZI:Z%NBlks)
     ENDIF
!
!    Blk1=P^(i-1)_i=[A^t_i].[Z^(i-1)_i]; Going down rows:(NxNI)^T.(NxNI)
!
     Blk1%D=Zero 
     StrtZI=Z%RowPt%I(IRow);StopZI=Z%RowPt%I(IRow+1)-1
     DO IDex=StrtZI,StopZI
        ICol=Z%ColPt%I(IDex)
        KDex=AiFlg%I(ICol)
        IF(KDex/=0)THEN                  
           ZiBlk=Z%BlkPt%I(IDex)
           AiBlk=A%BlkPt%I(KDex)
           M=BSiz%I(ICol)
           CALL DGEMM_NN(NI,M,NI,One,A%MTrix%D(AiBlk),Z%MTrix%D(ZiBlk),Blk1%D)
           PerfMon%FLOP=PerfMon%FLOP+DBLE(NI*M*NI)
        ENDIF
     ENDDO
!
!    P(I)=[P^(i-1)_i]^(-1)
!
     CALL FunkOnSqMat(NI,Inverse,Blk1%D,P%D(:,IRow))
!
!    Estimated performance; 2 DGEMMS+1 DSYEV
!
     PerfMon%FLOP=PerfMon%FLOP+DBLE((2+6)*NI**3)
!
!    Reset flags for column flags   
!
     DO I=StrtZI,ZDex-1;ZiFlg%I(ColPt%I(I))=0; ENDDO
     DO J=StrtAI,StopAI; AiFlg%I(A%ColPt%I(J))=0; ENDDO
!
  ENDDO ! end main loop over IRow
!
! Finishing touches on Z 
!
  Z%NAtms=NAtoms
  Z%NBlks=ZDex-1
  Z%NNon0=ZBlk-1
!
! Free some memory
!
  CALL Delete(ColPt)
  CALL Delete(BlkPt)
#ifdef SPATIAL_THRESHOLDING
#else 
  IF(TEST_AINV) &
     CALL Delete(A)
#endif
!
! Compute dimensions of DiagD & allocate it
!
  DiagD%NBlks=NAtoms
  DiagD%NNon0=0
  DO I=1,NAtoms
     DiagD%NNon0=DiagD%NNon0+BSiz%I(I)**2
  ENDDO
  CALL New(DiagD,(/NAtoms,DiagD%NBlks,DiagD%NNon0/))
!
! DiagD=P^(-1/2) in BCSR format
!
  DO I=1,NAtoms
     N=BSiz%I(I)
     CALL FunkOnSqMat(N,SqRoot,P%D(:,I),Blk1%D)
!    Estimated performance; 2 DGEMMS+1 DSYEV
     PerfMon%FLOP=PerfMon%FLOP+DBLE((2+6)*NI**3)
     P%D(1:N*N,I)=Blk1%D(1:N*N)     
  ENDDO
  CALL SetToI(DiagD,P)
!
! Free some more memory
!
  CALL Delete(P)
  CALL Delete(Blk1)
  CALL Delete(Blk2)
  CALL UnSetDSYEVWork()
!
! This is workspace for Z^t
!
  CALL New(Zt)!,(/NAtoms,Z%NBlks,Z%NNon0/))
!
! Symbolic transpose only, bloks in place 
!
  CALL XPose(Z)
! 
  ZBlksPreFilter=Z%NBlks
!
! Final Z=P^(-1/2).Z
!
  CALL Multiply(Z,DiagD,Zt)
  CALL Filter(Z,Zt)
  ZBlksPostFilter=Z%NBlks
!
! Full transpose
!
  CALL XPose(Z,Zt)
!  
! Account for multiplies AND adds in DGEMMs
!
  PerfMon%FLOP=PerfMon%Flop*Two
!
  CALL Elapsed_TIME(PerfMon,'Accum')
  CALL PPrint(PerfMon,Prog,Unit_O=6)  
!
#ifdef SPATIAL_THRESHOLDING 
  IF(TEST_AINV)THEN
!
!    Consistency check
!    
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
     Mssg=     'Max(Z^t.A.Z-I)='//TRIM(DblToShrtChar(Mx0))                 &
        //', DistanceThreshold='//TRIM(DblToShrtChar(AInvDistanceThresh))  &
        //', Thresholds%Trix='//TRIM(DblToShrtChar(Thresholds%Trix))
!
     IF(Mx0>1.D1*Thresholds%Trix)THEN
        CALL Halt('In BlokAInv, failed test: '//TRIM(Mssg))
     ELSE
        Mssg=Prog//' :: '//TRIM(Mssg)
        WRITE(*,*)TRIM(Mssg)
        CALL OpenASCII(OutFile,Out)
        CALL PrintProtectL(Out)
        WRITE(Out,*)TRIM(Mssg)
        CALL PrintProtectR(Out)
        CLOSE(UNIT=Out,STATUS='KEEP')
!
!        CALL GetEnv('MONDO_WORK',AInvFile)
!        AInvFile=TRIM(AInvFile)//'/AInv.dat'
!        CALL OpenASCII(AInvFile,35)
!        WRITE(35,*)NAtoms,ZBlksPreFilter,ZBlksPostFilter,Mx0
!        CLOSE(35)
      ENDIF
   ENDIF
#endif
!
!  Put Z and ZT to disk
!  
   CALL Put(Z,TrixFile('Z',Args))
   CALL Put(Zt,TrixFile('ZT',Args))
!
!  Debug
!
   CALL PChkSum(Z,'Z',Prog)
   CALL PPrint(Z,'Z')
   CALL Plot(Z,'Z')
!
!  Tidy up
!
   CALL Delete(Z)
   CALL Delete(ZT)
   CALL Delete(DiagD)
   CALL ShutDown(Prog) 
END PROGRAM 







