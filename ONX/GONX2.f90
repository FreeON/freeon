!------------------------------------------------------------------------------
!    This code is part of the FreeON suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the FreeON group for review, and possible
!    dissemination in future releases.
!------------------------------------------------------------------------------

#include "MondoConfig.h"

PROGRAM GONX2

#ifndef PARALLEL
#undef ONX2_PARALLEL
#endif

  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE Functionals

  USE ONX2DataType
  USE ONXCtrSclg   , ONLY: TrnMatBlk
  USE ONX2List     , ONLY: AllocList,DeAllocList,MakeGList,PrintList
  USE GONX2ComputDK, ONLY: ComputDK

#ifdef ONX2_PARALLEL
  USE MondoMPI
  USE FastMatrices
  USE ONXGet       , ONLY: Get_Essential_RowCol,GetOffArr,Set_DFASTMAT_EQ_DBCSR2,GetDab
  USE PartDrv      , ONLY: PDrv_Initialize,PDrv_Finalize
#else
  USE ONXGet       , ONLY: GetOffArr
#endif
  !
  IMPLICIT NONE
  !
#ifdef ONX2_PARALLEL
  TYPE(FASTMAT), POINTER     :: DFMcd,DFMab,KxFM
#else
  TYPE(BCSR)                 :: D
#endif
  TYPE(BSET)                 :: BSc
  TYPE(CRDS)                 :: GMc
  TYPE(ARGMT)                :: Args
!--------------------------------------------------------------------------------
#ifdef ONX2_PARALLEL
  TYPE(DBL_RNK2)             :: GradTmp
  TYPE(INT_VECT)             :: APt,BPt,CPt,DPt
#endif
  TYPE(DBL_RNK2)             :: GradX,GradAux,BoxX
  TYPE(DBL_VECT)             :: GTmp
  REAL(DOUBLE)               :: KScale
!--------------------------------------------------------------------------------
#ifdef ONX2_PARALLEL
  TYPE(DBL_VECT)             :: TmGxArr,TmMLArr,TmTMArr,TmALArr,TmDLAr,TmDLArr
!  REAL(DOUBLE), EXTERNAL     :: MondoTimer
  INTEGER                    :: CMin,CMax,DMin,DMax,IErr
  INTEGER                    :: ANbr,BNbr,CNbr,DNbr
#endif
  REAL(DOUBLE)               :: Time1,Time2
  REAL(DOUBLE)               :: TmTM,TmML,TmGx,TmAL,TmDL
  CHARACTER(LEN=*),PARAMETER :: Prog='GONX2'
  LOGICAL                    :: DoStrs

  TYPE(INT_RNK2) :: OffArr
#ifdef ONX2_PARALLEL
  TYPE(CList), DIMENSION(:), POINTER :: ListC,ListD
#else
  TYPE(CList), DIMENSION(:), POINTER :: ListC
#endif

  INTEGER :: ixyz,jxyz,A1,A2
  DoStrs=.TRUE.!.FALSE.

#ifdef ONX2_PARALLEL
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
#else
  CALL StartUp(Args,Prog)
#endif

  ! [FIXME] Fix stack (set it to unlimited). This should really be fixed so that
  ! we don't have to screw with the stack.
  CALL UnlimitStack()

  CALL Get(BSiz ,'atsiz',Tag_O=PrvBase)
  CALL Get(OffS ,'atoff',Tag_O=PrvBase)
  CALL Get(NBasF,'nbasf',Tag_O=PrvBase)

  CALL Get(BSc,Tag_O=CurBase)
  CALL Get(GMc,Tag_O=CurGeom)

  !------------------------------------------------
  ! Initialization and allocations.
  !
#ifdef ONX2_PARALLEL
  NULLIFY(ListC,ListD)
#else
  NULLIFY(ListC)
#endif

  TmTM=Zero;TmML=Zero;TmGx=Zero;TmAL=Zero;TmDL=Zero
  CALL New(GradX,(/3,NAtoms/))
  CALL DBL_VECT_EQ_DBL_SCLR(3*NAtoms,GradX%D(1,1),0.0D0)

  CALL New(GradAux,(/3,NAtoms/))
  CALL DBL_VECT_EQ_DBL_SCLR(3*NAtoms,GradAux%D(1,1),0.0D0)

  CALL New(OffArr,(/BSc%NCtrt,BSc%NKind/))

  !STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS 
  CALL New(BoxX,(/3,3/))

  CALL DBL_VECT_EQ_DBL_SCLR(9,BoxX%D(1,1),0.0D0)
  !STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS 

  CALL GetBufferSize(GMc,BSc,GMc,BSc)

  CALL GetOffArr(OffArr,BSc)

  !------------------------------------------------
  ! Get denstiy matrix.
  !
  SELECT CASE(SCFActn)
  CASE('ForceEvaluation')
    ! We are getting NSMat now from hdf via StartUp.
#ifdef ONX2_PARALLEL
     CALL PDrv_Initialize(DFMcd,TrixFile('D',Args,1),'GONXPart',Args,WhenPartS_O='GEO')
     !NSMat=DFMcd%NSMat
#else
     CALL Get(D,TrixFile('D',Args,1))
     !NSMat=D%NSMat
#endif
  CASE DEFAULT
     CALL Halt('GONX2: Do not recognize this action <'//TRIM(SCFActn)//'>.')
  END SELECT

  !------------------------------------------------
  ! Normalize the density matrix
  !
#ifdef ONX2_PARALLEL
  Time1 = MondoTimer()
  CALL TrnMatBlk(BSc,GMc,DFMcd)
  Time2 = MondoTimer()
#else
  CALL CPU_TIME(Time1)
  CALL TrnMatBlk(BSc,GMc,D)
  CALL CPU_TIME(Time2)
#endif
  TmTM = Time2-Time1
  !
  !------------------------------------------------
  ! Allocate the list(s) and get the buffer sizes.
  !
#ifdef ONX2_PARALLEL
  Time1 = MondoTimer()
  CALL Get_Essential_RowCol(DFMcd,CPt,CNbr,CMin,CMax,DPt,DNbr,DMin,DMax)
  !write(*,*) 'CMin',CMin,'CMax',CMax,'MyID',MyID
  !write(*,*) 'DMin',DMin,'DMax',DMax,'MyID',MyID
  !
  CALL AllocList(ListC,CMin,CMax)
  CALL AllocList(ListD,DMin,DMax)
  Time2 = MondoTimer()
#else
  CALL CPU_TIME(Time1)
  CALL AllocList(ListC,1,NAtoms)
  CALL CPU_TIME(Time2)
#endif
  TmAL = Time2-Time1
  !
  !------------------------------------------------
  ! Make the distribution list(s).
  !
#ifdef ONX2_PARALLEL
  Time1 = MondoTimer()
  CALL MakeGList(ListC,GMc,BSc,CS_OUT,CPt,CNbr,APt,ANbr)
  CALL MPI_Barrier(MONDO_COMM,IErr)
  CALL MakeGList(ListD,GMc,BSc,CS_OUT,DPt,DNbr,BPt,BNbr)
  CALL MPI_Barrier(MONDO_COMM,IErr)
  Time2 = MondoTimer()
#else
  CALL CPU_TIME(Time1)
  CALL MakeGList(ListC,GMc,BSc,CS_OUT)
  CALL CPU_TIME(Time2)
#endif
  TmML = Time2-Time1
  !
  !------------------------------------------------
  ! Print list.
  !
#ifdef ONX2_DBUG
#ifdef ONX2_PARALLEL
  CALL PrintList(ListC)
  CALL PrintList(ListD)
#else
  CALL PrintList(ListC)
#endif
#endif
  !
  !------------------------------------------------
  ! Get the second density matrix and normalization.
  !
#ifdef ONX2_PARALLEL
  CALL GetDab(DFMab,APt,ANbr,BPt,BNbr,Args)
  Time1 = MondoTimer()
  CALL TrnMatBlk(BSc,GMc,DFMab)
  Time2 = MondoTimer()
  TmTM = TmTM+Time2-Time1
#endif

  !------------------------------------------------
  ! Compute Exchange Forces.
  !
#ifdef ONX2_PARALLEL
  Time1 = MondoTimer()
  CALL ComputDK(DFMcd,DFMab,GradX,BoxX,DoStrs,ListC,ListD,OffArr,GMc,BSc,CS_OUT)
  Time2 = MondoTimer()
#else
  CALL CPU_TIME(Time1)
  ! Nick, December 5 2007:
  !
  ! With the intel ifort compiler (and maybe others), and a small stacksize
  ! limit, the call to ComputDK() might fail with a segmentation fault. Raising
  ! the stacksize limit fixes this. We should fix this call to prevent this from
  ! happening even with small stacksize limit.
  CALL ComputDK(D,GradX,BoxX,DoStrs,ListC,ListC,OffArr,GMc,BSc,CS_OUT)
  CALL CPU_TIME(Time2)
#endif
  TmGx = Time2-Time1

  !------------------------------------------------
  ! Free up some space. Deallocate the list(s).
  !
#ifdef ONX2_PARALLEL
  Time1 = MondoTimer()
  CALL DeAllocList(ListC)
  CALL DeAllocList(ListD)
  Time2 = MondoTimer()
#else
  CALL CPU_TIME(Time1)
  CALL DeAllocList(ListC)
  CALL CPU_TIME(Time2)
#endif
  TmDL = Time2-Time1
  !
  !------------------------------------------------
  ! Redistribute partition informations.
  !
#ifdef ONX2_PARALLEL
  CALL PDrv_Finalize(DFMcd,CollectInPar_O=.TRUE.)
  CALL Delete_FastMat1(DFMcd)
  CALL Delete_FastMat1(DFMab)
#else
  CALL Delete(D)
#endif
  !
  !------------------------------------------------
  ! Add Exchange Gradient.
  !
  CALL Get(GradAux,'Gradients',Tag_O=CurGeom)
  KScale=ExactXScale(ModelChem)
  IF(NSMat.GT.1) CALL DSCAL(3*NAtoms,0.5D0,GradX%D(1,1),1)
#ifdef ONX2_PARALLEL
  CALL New(GradTmp,(/3,NAtoms/))
  CALL DBL_VECT_EQ_DBL_SCLR(3*NAtoms,GradTmp%D(1,1),0.0d0)
  CALL MPI_REDUCE(GradX%D(1,1),GradTmp%D(1,1),3*NAtoms,MPI_DOUBLE_PRECISION, &
       &          MPI_SUM,ROOT,MONDO_COMM,IErr)
! Print out.
  CALL New(GTmp,3*NAtoms)
  DO IXYZ=1,GMc%NAtms
     A1=3*(IXYZ-1)+1
     A2=3*IXYZ
     GTmp%D(A1:A2) = GradTmp%D(1:3,IXYZ)
  ENDDO
  CALL PChkSum(GTmp,'dKx/dR['//TRIM(CurGeom)//']',Proc_O=Prog)  
  CALL DAXPY(3*NAtoms,KScale,GradTmp%D(1,1),1,GradAux%D(1,1),1)
! Print Out Forces
  CALL Print_Force(GMc,GTmp,'X Force')
  CALL Print_Force(GMc,GTmp,'X Force',Unit_O=6)
  CALL Delete(GTmp)
  CALL Delete(GradTmp)
! STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS 
  CALL New(GradTmp,(/3,3/))
  CALL DBL_VECT_EQ_DBL_SCLR(9,GradTmp%D(1,1),0.0d0)
  IF(DoStrs) THEN
     IF(NSMat.GT.1) CALL DSCAL(9,0.5D0,BoxX%D(1,1),1)
     CALL MPI_REDUCE(BoxX%D(1,1),GradTmp%D(1,1),9,MPI_DOUBLE_PRECISION, &
          &          MPI_SUM,ROOT,MONDO_COMM,IErr)
!    Zero the Lower Triange
     DO IXYZ=1,3
        DO JXYZ=1,IXYZ-1
           GradTmp%D(IXYZ,JXYZ)   = 1.D8
        ENDDO
     ENDDO
     CALL DAXPY(9,KScale,GradTmp%D(1,1),1,GMc%PBC%LatFrc%D(1,1),1)
  ENDIF
! Print Out the Lattice Forces
  CALL Print_LatForce(GMc,GradTmp%D,'X Lattice Force')
  CALL Print_LatForce(GMc,GradTmp%D,'X Lattice Force',Unit_O=6)
  CALL Delete(GradTmp)
#else
  ! Print out.
  CALL New(GTmp,3*NAtoms)
  DO IXYZ=1,GMc%NAtms
     A1=3*(IXYZ-1)+1
     A2=3*IXYZ
     GTmp%D(A1:A2) = GradX%D(1:3,IXYZ)
  ENDDO
!
  CALL PChkSum(GTmp,'dKx/dR['//TRIM(CurGeom)//']',Proc_O=Prog) 
! Print Out Forces
  CALL Print_Force(GMc,GTmp,'X Force')
  CALL Print_Force(GMc,GTmp,'X Force',Unit_O=6)
  CALL Delete(GTmp)
  CALL DAXPY(3*NAtoms,KScale,GradX%D(1,1),1,GradAux%D(1,1),1)
!
! STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS 
  IF(DoStrs) THEN
     IF(NSMat.GT.1) CALL DSCAL(9,0.5D0,BoxX%D(1,1),1)
!    Zero the Lower Triange
     DO IXYZ=1,3
        DO JXYZ=1,IXYZ-1
           BoxX%D(IXYZ,JXYZ) = Zero
        ENDDO
     ENDDO
     CALL DAXPY(9,KScale,BoxX%D(1,1),1,GMc%PBC%LatFrc%D(1,1),1)
  ENDIF
! Print Out the Lattice Forces
  CALL Print_LatForce(GMc,BoxX%D,'X Lattice Force')
  CALL Print_LatForce(GMc,BoxX%D,'X Lattice Force',Unit_O=6)
#endif
  !------------------------------------------------
  ! Save Exchange Gradients and Stress.
  !
#if defined(PARALLEL_CLONES)
  IF(MRank(MPI_COMM_WORLD) == ROOT) THEN
    CALL MondoLog(DEBUG_NONE, Prog, "writing gradients to hdf", "Clone "//TRIM(IntToChar(MyClone)))
    CALL Put(GradAux,'Gradients',Tag_O=CurGeom)
  ELSE
    CALL MondoLog(DEBUG_NONE, Prog, "sending gradients to clone 1", "Clone "//TRIM(IntToChar(MyClone)))
  ENDIF
#else
  CALL Put(GradAux,'Gradients',Tag_O=CurGeom)
#endif
  !
  !STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS 
  IF(DoStrs) THEN
#if defined(PARALLEL_CLONES)
    IF(MRank(MPI_COMM_WORLD) == ROOT) THEN
      CALL MondoLog(DEBUG_NONE, Prog, "writing stress to hdf", "Clone "//TRIM(IntToChar(MyClone)))
      CALL Put(GMc%PBC%LatFrc,'latfrc',Tag_O=CurGeom)
    ELSE
      CALL MondoLog(DEBUG_NONE, Prog, "sending stress to clone 1", "Clone "//TRIM(IntToChar(MyClone)))
    ENDIF
#else
    CALL Put(GMc%PBC%LatFrc,'latfrc',Tag_O=CurGeom)
#endif
  ENDIF
  !STRESS STRESS STRESS STRESS STRESS STRESS STRESS STRESS 
  !
  !------------------------------------------------
  ! Timing.
  !
#ifdef GONX2_INFO
#ifdef ONX2_PARALLEL
  !
  ! End Total Timing
  CALL New(TmGxArr,NPrc)
  CALL New(TmMLArr,NPrc)
  CALL New(TmTMArr,NPrc)
  CALL New(TmALArr,NPrc)
  CALL New(TmDLArr,NPrc)
  !
  CALL MPI_Gather(TmGx,1,MPI_DOUBLE_PRECISION,TmGxArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmML,1,MPI_DOUBLE_PRECISION,TmMLArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmTM,1,MPI_DOUBLE_PRECISION,TmTMArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmAL,1,MPI_DOUBLE_PRECISION,TmALArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmDL,1,MPI_DOUBLE_PRECISION,TmDLArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  !
  IF(MyID.EQ.ROOT) THEN
     !
     ! Imbalance stuff.
     CALL PImbalance(TmGxArr ,NPrc,Prog_O='ComputeGK')
     !
     WRITE(*,1001) SUM(TmALArr%D )/DBLE(NPrc),MINVAL(TmALArr%D ),MAXVAL(TmALArr%D )
     WRITE(*,1002) SUM(TmMLArr%D )/DBLE(NPrc),MINVAL(TmMLArr%D ),MAXVAL(TmMLArr%D )
     WRITE(*,1003) SUM(TmTMArr%D )/DBLE(NPrc),MINVAL(TmTMArr%D ),MAXVAL(TmTMArr%D )
     WRITE(*,1004) SUM(TmGxArr%D )/DBLE(NPrc),MINVAL(TmGxArr%D ),MAXVAL(TmGxArr%D )
     WRITE(*,1005) SUM(TmDLArr%D )/DBLE(NPrc),MINVAL(TmDLArr%D ),MAXVAL(TmDLArr%D )
     !WRITE(*,1006) SUM(NERIsArr%D)           ,MINVAL(NERIsArr%D),MAXVAL(NERIsArr%D)
     !
1001 FORMAT(' GONX: Ave TmAL = ',F15.2,', Min TmAL = ',F15.2,', Max TmAL = ',F15.2)
1002 FORMAT(' GONX: Ave TmML = ',F15.2,', Min TmML = ',F15.2,', Max TmML = ',F15.2)
1003 FORMAT(' GONX: Ave TmTM = ',F15.2,', Min TmTM = ',F15.2,', Max TmTM = ',F15.2)
1004 FORMAT(' GONX: Ave TmGx = ',F15.2,', Min TmGx = ',F15.2,', Max TmGx = ',F15.2)
1005 FORMAT(' GONX: Ave TmDL = ',F15.2,', Min TmDL = ',F15.2,', Max TmDL = ',F15.2)
     !1006 FORMAT(' ONX: Tot ERI  = ',F15.2,', Min ERI  = ',F15.2,', Max ERI  = ',F15.2)
  ENDIF
  !
  CALL Delete(TmALArr)
  CALL Delete(TmMLArr)
  CALL Delete(TmTMArr)
  CALL Delete(TmGxArr)
  CALL Delete(TmDLArr)
  !
#else
  !
  WRITE(*,1001) TmAL
  WRITE(*,1002) TmML
  WRITE(*,1003) TmTM
  WRITE(*,1004) TmGx
  WRITE(*,1005) TmDL
  !WRITE(*,1006) ....
  !
1001 FORMAT(' GONX: Ave TmAL = ',F15.2)
1002 FORMAT(' GONX: Ave TmML = ',F15.2)
1003 FORMAT(' GONX: Ave TmTM = ',F15.2)
1004 FORMAT(' GONX: Ave TmGx = ',F15.2)
1005 FORMAT(' GONX: Ave TmDL = ',F15.2)
  !1006 FORMAT(' ONX: Tot ERI  = ',F15.2)
  !
#endif
#endif
  !
  !------------------------------------------------
  ! Deallocate Arrays.
  !
  CALL Delete(GradX)
  CALL Delete(GradAux)
  CALL Delete(BoxX)
  CALL Delete(OffArr)
  !
  CALL Delete(BSc)
  CALL Delete(GMc)
  !
  CALL ShutDown(Prog)
  !
END PROGRAM GONX2

