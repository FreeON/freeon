!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!    COMPUTE THE INCOMPLETE INVERSE CHOLESKY FACTOR 
!    OF THE OVERLAP MATRIX Z=S^(-L) USING BENZI AND TUMAS AINV
!
!
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
  IMPLICIT NONE
  TYPE(BCSR)          :: A,Z,Zt,DiagD
#ifdef EXTREME_DEBUG
  TYPE(BCSR)          :: T1,T2
  TYPE(DBL_RNK2)      :: DnsD1,DnsD2,DnsD3,DnsZ1,DnsZ2,DnsZ3
#endif
  TYPE(BSET)          :: BS
  TYPE(CRDS)          :: GM
  TYPE(INT_VECT)      :: AiFlg,ZiFlg,ColPt,BlkPt
  TYPE(ARGMT)         :: Args
  INTEGER :: I,J,Q,R,IDex,JDex,ZDex,ZBlk,NIJ, &
             n,ni,msiz,strtai,stopai,strtaj,stopaj,strtzi,stopzi,nj,strtzj,stopzj,jcol,k,kdex, &
             aiblk,ajblk,zjblk,m,ziblk,icol,zrowpt,zcolpt,zblkpt,NewBloks,EndBloks,IRow,JRow
  TYPE(DBL_VECT) :: Blk1,Blk2
  TYPE(DBL_RNK2) :: P
  REAL(DOUBLE) :: Op,Mx0,B2Norm
  REAL(DOUBLE), EXTERNAL :: DDOT
  REAL(DOUBLE), DIMENSION(10,10) :: T
  CHARACTER(LEN=8),PARAMETER                :: Prog='BlokAInv'
#ifdef EXTREME_PRINT_DEBUG
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: ZIChar,ZJChar,AIChar, AJChar
#endif
!----------------------------------------------------------------------------------------------------------- 
! Start up macro
!
  CALL StartUp(Args,Prog)
!
! Get basis set and geometry
!
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
!
! Allocations 
!
  CALL New(Z)
  CALL New(A)
  CALL Get(A,TrixFile('S',Args))
! Set global workspace for FunkOnSqMat
  CALL SetDSYEVWork(MaxBlkSize)
! Allocate intermediate blocks 
  CALL New(Blk1,MaxBlkSize*MaxBlkSize)
  CALL New(Blk2,MaxBlkSize*MaxBlkSize)
! Allocate diagonal "pivot" blocks
  WRITE(*,*)' MaxBlkSize = ',MaxBlkSize,' NAtoms = ',NAtoms
  CALL New(P,(/MaxBlkSize*MaxBlkSize,NAtoms/))
  WRITE(*,*)' MaxBlkSize = ',MaxBlkSize,' NAtoms = ',NAtoms
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
        NJ=BSiz%I(JRow)
        StrtAJ=A%RowPt%I(JRow);StopAJ=A%RowPt%I(JRow+1)-1
        StrtZJ=Z%RowPt%I(JRow);StopZJ=Z%RowPt%I(JRow+1)-1
!
!       Blk1=P^(j-1)_i=[A^t_j].[Z^(j-1)_i]; Going down rows over N: (NxNJ)^T.(NxNI)
!
        NIJ=NI*NJ
        Blk1%D(1:NIJ)=Zero
        DO J=StrtAJ,StopAJ
           JDex=A%ColPt%I(J)
           IDex=ZiFlg%I(JDex)
           IF(IDex/=0)THEN                 
              ZiBlk=BlkPt%I(IDex) 
              AjBlk=A%BlkPt%I(J)
              M=BSiz%I(JDex)
              CALL DGEMM_NN(NJ,M,NI,One,A%MTrix%D(AjBlk),Z%MTrix%D(ZiBlk),Blk1%D)
              PerfMon%FLOP=PerfMon%FLOP+DBLE(NIJ*M)
           ENDIF
        ENDDO
!
!       Blk2=[P^(j-1)_j]^(-1).[P^(j-1)_i]
!
        CALL DGEMM_NNc(NJ,NJ,NI,One,Zero,P%D(:,JRow),Blk1%D,Blk2%D)
        PerfMon%FLOP=PerfMon%FLOP+DBLE(NIJ*NJ)
! 
!       Check the magintude of Blk2.  Update Z_I only if Blk2 is "large" enough.       
!        
        B2Norm=SQRT(DDOT(NI*NJ,Blk2%D,1,Blk2%D,1))
        PerfMon%FLOP=PerfMon%FLOP+DBLE(NIJ)
!        
        IF(B2Norm>Thresholds%Trix*1.D-1)THEN

!
!          Z^j_i=Z^(j-1)_i-[Z^(j-1)_j].{[P^(j-1)_j]^(-1).[P^(j-1)_i]}
!          Update going down rows:(NxNI)=(NxNI)+(NxNJ).(NJxNI)
!
           DO JDex=StrtZJ,StopZJ
              JCol =ColPt%I(JDex) 
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

           ENDDO
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
#ifdef EXTREME_DEBUG
#else 
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
!
     PerfMon%FLOP=PerfMon%FLOP+DBLE((2+6)*NI**3)
!
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
! Final Z=P^(-1/2).Z
!
  CALL Multiply(Z,DiagD,Zt)
  CALL Filter(Z,Zt)
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
#ifdef EXTREME_DEBUG
!
! Consistency check
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
  WRITE(*,*)' Mx0 = ',Mx0
! CALL ShutDown(Prog)
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







