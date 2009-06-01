!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
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
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
!    GENERIC IO ROUTINES FOR MONDOSCF TYPES
!    Author: Matt Challacombe and CK Gan
!-------------------------------------------------------------------------------

#include "MondoConfig.h"

MODULE InOut
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE ProcessControl
  USE MemMan
  USE Indexing
  USE Parse
  USE MondoLogger
  USE Utilities

#ifdef PARALLEL
  USE MondoMPI
#endif
  USE SetXYZ

  IMPLICIT NONE

  INTERFACE Get
    MODULE PROCEDURE Get_INT_SCLR, Get_INT_VECT, Get_INT_RNK2, &
         Get_INT_RNK3, Get_INT_RNK4, Get_DBL_SCLR, &
         Get_DBL_VECT, Get_DBL_RNK2, Get_DBL_RNK3, &
         Get_DBL_RNK4, Get_CHR_SCLR, Get_CHR10_VECT,&
         Get_Sp1x1,                                &
         Get_BondD,    Get_AtmB,                   &
         Get_CHR_VECT, Get_LOG_VECT, Get_INTC,     &
         Get_LOG_SCLR, Get_BSET,     Get_CRDS,     &
         Get_TOLS,     Get_BCSR,     Get_PBCFit,   &
#ifdef PARALLEL
         Get_DBCSR,                                &
#endif
         Get_PBCInfo,  Get_CellSet,                &
         Get_ARGMT,    Get_HGRho,                  &
         Get_CMPoles

  END INTERFACE
  INTERFACE Put
    MODULE PROCEDURE Put_INT_SCLR, Put_INT_VECT, Put_INT_RNK2, &
         Put_INT_RNK3, Put_INT_RNK4, Put_DBL_SCLR, &
         Put_DBL_VECT, Put_DBL_RNK2, Put_DBL_RNK3, &
         Put_DBL_RNK4, Put_DBL_RNK6, Put_CHR_SCLR, &
         Put_CHR_VECT, Put_LOG_VECT, Put_CHR10_VECT,&
         Put_Sp1x1,    Put_BondD,                  &
         Put_AtmB,     Put_INTC,     Put_PBCFit,   &
         Put_LOG_SCLR, Put_BSET,     Put_CRDS,     &
#ifdef PARALLEL
         Put_DBCSR,                                &
#endif
         Put_PBCInfo,  Put_CellSet,                &
         Put_TOLS,     Put_BCSR,     Put_HGRho,    &
         Put_CMPoles
  END INTERFACE

  INTEGER, EXTERNAL  :: HDF5CreateFile,HDF5OpenFile,HDF5CloseFile,         &
       HDF5CreateGroup,HDF5OpenGroup,HDF5CloseGroup,      &
       HDF5CloseData,HDF5ReadIntegerVector,               &
       HDF5ReadDoubleVector,HDF5WriteIntegerVector,       &
       HDF5WriteDoubleVector,HDF5SizeOfData

  INTEGER, PARAMETER :: NATIVE_DOUBLE=6, NATIVE_INT32=24

  TYPE META_DATA
    INTEGER            :: Dimension
    INTEGER            :: DataType
    INTEGER            :: DataId
    INTEGER            :: DataSpc
    INTEGER            :: Status
    INTEGER            :: Unlimited
    CHARACTER(LEN=DCL) :: VarName
  END TYPE META_DATA

  INTEGER, SAVE      :: HDF_CurrentID

CONTAINS
  !===============================================================================
  !  INITIALIZE AN HDF FILE
  !===============================================================================
  FUNCTION InitHDF(FileName) RESULT(FileID)
    CHARACTER(LEN=*), INTENT(IN) :: FileName
    INTEGER                      :: NC,STATUS,FileID
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      NC=StringLen(FileName)
      FileID=HDF5CreateFile(NC,Char2Ints(NC,FileName))
      IF(FileID==FAIL) &
           CALL Halt(' Failed to create the HDF file <'//TRIM(FileName)//'>.')
#ifdef PARALLEL
    ENDIF
#endif
  END FUNCTION InitHDF
  !===============================================================================
  !    Open a HDF file
  !===============================================================================
  FUNCTION OpenHDF(FileName) RESULT(FileID)
    CHARACTER(LEN=*),INTENT(IN) :: FileName
    INTEGER                     :: NC,STATUS
    INTEGER                     :: FileID
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      NC=StringLen(FileName)
      FileID=HDF5OpenFile(NC,Char2Ints(NC,FileName))
      IF(FileID==FAIL) THEN
        CALL Halt(' Failed to open the HDF file <'//TRIM(FileName)//'>.')
      ENDIF
#ifdef PARALLEL
    ENDIF
#endif
  END FUNCTION OpenHDF
  !===============================================================================
  ! CLOSE AN HDF FILE
  !===============================================================================
  SUBROUTINE CloseHDF(FileID)
    INTEGER :: Status,FileID
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      Status=HDF5CloseFile(FileID)
      IF(Status==FAIL) THEN
        CALL Halt(' Failed to close an HDF file with HDF_CurrentID=' &
             //'<'//TRIM(IntToChar(FileID))//'>. ')
      ENDIF
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE CloseHDF
  !===============================================================================
  ! INITIALIZE A HDF GROUP
  !===============================================================================
  FUNCTION InitHDFGroup(FileID,GroupName) RESULT(GroupID)
    CHARACTER(LEN=*), INTENT(IN) :: GroupName
    INTEGER                      :: FileID,NC,STATUS,GroupID
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      NC=StringLen(GroupName)
      GroupID=HDF5CreateGroup(FileID,NC,Char2Ints(NC,GroupName))
      IF(GroupID==FAIL) &
           CALL Halt(' Failed to create the Group file <'//TRIM(GroupName)//'>.')
#ifdef PARALLEL
    ENDIF
#endif
  END FUNCTION InitHDFGroup
  !===============================================================================
  ! OPEN A HDF GROUP
  !===============================================================================
  FUNCTION OpenHDFGroup(FileID,GroupName) RESULT(GroupID)
    CHARACTER(LEN=*),INTENT(IN) :: GroupName
    INTEGER                     :: NC,STATUS
    INTEGER                     :: FileID,GroupID
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      NC=StringLen(GroupName)
      GroupID=HDF5OpenGroup(FileID,NC,Char2Ints(NC,GroupName))
      IF(GroupID==FAIL) &
           CALL Halt(' Failed to open the Group file <'//TRIM(GroupName)//'>.')
#ifdef PARALLEL
    ENDIF
#endif
  END FUNCTION OpenHDFGroup
  !===============================================================================
  ! CLOSE A HDF GROUP
  !===============================================================================
  SUBROUTINE CloseHDFGroup(GroupID)
    INTEGER :: Status,GroupID
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      Status=HDF5CloseGroup(GroupID)
      IF(Status==FAIL) &
           CALL Halt(' Failed to close a Group with GroupID=' &
           //'<'//TRIM(IntToChar(GroupID))//'>. ')
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE CloseHDFGroup
  !===============================================================================
  !
  !===============================================================================
  SUBROUTINE CreateData(Meta)
    TYPE(META_DATA) :: Meta
    INTEGER         :: NC
    NC=StringLen(Meta%VarName)

    !    IF(Meta%Unlimited==1)CALL Halt(' UNLIMIT! ')
    CALL HDF5CreateData(HDF_CurrentID,Meta%DataType,Meta%Dimension,          &
         NC,Char2Ints(NC,Meta%VarName),Meta%Unlimited, &
         Meta%DataId,Meta%DataSpc)
    IF(Meta%DataId==FAIL) &
         CALL Halt(' Failed in CreateData:'//TRIM(MetaChar(Meta)))
    Meta%Status=SUCCEED
  END SUBROUTINE CreateData
  !===============================================================================
  !
  !===============================================================================
  SUBROUTINE OpenData(Meta,Put_O)
    TYPE(META_DATA)  :: Meta
    LOGICAL,OPTIONAL :: Put_O
    INTEGER          :: I,NC,SizeOf,SIZE_OF_HDF5_DATA
    INTEGER, DIMENSION(:), ALLOCATABLE :: VarNameInts

    Meta%Status=FAIL
    NC=StringLen(Meta%VarName)
    ALLOCATE(VarNameInts(NC))
    VarNameInts = Char2Ints(NC,Meta%VarName)
    CALL HDF5OpenData(HDF_CurrentID,NC,VarNameInts,Meta%DataId,Meta%DataSpc)
    DEALLOCATE(VarNameInts)
    IF(Meta%DataId==FAIL)THEN
      IF(PRESENT(Put_O))THEN
        CALL CreateData(Meta)
      ELSE
        CALL Halt('Failed in OpenData:'//TRIM(MetaChar(Meta)))
      ENDIF
    ELSE
      SizeOf=HDF5SizeOfData(Meta%DataSpc)
      IF(PRESENT(Put_O))THEN
        IF(SizeOf<Meta%Dimension.AND.Meta%Unlimited==1)THEN
          CALL HDF5ExtendData(Meta%DataId,Meta%DataSpc,Meta%Dimension)
        ELSEIF(SizeOf>Meta%Dimension)THEN
          CALL HDF5SelectData(Meta%DataId,Meta%DataSpc,Meta%Dimension)
        ELSEIF(SizeOf<Meta%Dimension)THEN
          CALL Halt('Failed in OpenData, need to redimension fixed form data: '//TRIM(MetaChar(Meta)))
        ENDIF
      ELSEIF(SizeOf>Meta%Dimension)THEN
        CALL HDF5SelectData(Meta%DataId,Meta%DataSpc,Meta%Dimension)
      ENDIF
    ENDIF
  END SUBROUTINE OpenData
  !===============================================================================
  !
  !===============================================================================
  SUBROUTINE CloseData(Meta)
    TYPE(META_DATA) :: Meta
    IF(Meta%Status==FAIL) THEN
      CALL Halt(' HDF R/W error for '//TRIM(MetaChar(Meta)))
    ENDIF
    Meta%Status=HDF5CloseData(Meta%DataId,Meta%DataSpc)
    IF(Meta%Status==FAIL) THEN
      CALL Halt(' HDF5CloseData error for '//TRIM(MetaChar(Meta)))
    ENDIF
  END SUBROUTINE CloseData
  !===============================================================================
  !
  !===============================================================================
  FUNCTION NameTag(VarName,Tag_O,Stats_O) RESULT(FullName)
    CHARACTER(LEN=*),         INTENT(IN)    :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
    INTEGER,DIMENSION(3),OPTIONAL,INTENT(IN):: Stats_O
    INTEGER,DIMENSION(3)                    :: ReCycStats
    CHARACTER(LEN=DEFAULT_CHR_LEN)          :: FullName
    INTEGER                                 :: Tag_INT_old, Tag_INT_new
    IF(PRESENT(Tag_O).AND.PRESENT(Stats_O))THEN
      CALL Halt(' Logic error in NameTag with both tag and stats passed through! ')
    ELSEIF(PRESENT(Tag_O))THEN
      ! Apply recycling.
      Tag_INT_old = CharToInt(Tag_O)
      Tag_INT_new = MOD(Tag_INT_old, RecycleHDF)
      FullName=TRIM(VarName)//TRIM(IntToChar(Tag_INT_new))
    ELSEIF(PRESENT(Stats_O))THEN
      ReCycStats=Stats_O
      ReCycStats(3)=MOD(ReCycStats(3), RecycleHDF)
      FullName=TRIM(VarName)//TRIM(StatsToChar(ReCycStats))
    ELSE
      FullName=TRIM(VarName)
    ENDIF
    CALL LowCase(FullName)
  END FUNCTION NameTag
  !===============================================================================
  !
  !===============================================================================
  FUNCTION SetMeta(VarName,DataType,SizeOf,Unlimit_O) RESULT(Meta)
    TYPE(META_DATA)              :: Meta
    CHARACTER(LEN=*), INTENT(IN) :: VarName
    INTEGER,          INTENT(IN) :: DataType,SizeOf
    LOGICAL,OPTIONAL, INTENT(IN) :: Unlimit_O
    IF(PRESENT(Unlimit_O))THEN
      IF(Unlimit_O)THEN
        Meta%Unlimited=1
      ELSE
        Meta%Unlimited=0
      ENDIF
    ELSE
      Meta%Unlimited=0
    ENDIF
    Meta%DataType=DataType
    Meta%Dimension=SizeOf
    Meta%VarName=VarName
    Meta%DataId=FAIL
    Meta%Status=FAIL
  END FUNCTION SetMeta
  !===============================================================================
  !
  !===============================================================================
  FUNCTION MetaChar(Meta) RESULT(CharMeta)
    TYPE(META_DATA)                :: Meta
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: CharMeta
    CHARACTER(LEN=1)               :: TF
    IF(Meta%Unlimited==1)THEN
      TF='T'
    ELSE
      TF='F'
    ENDIF
    CharMeta=Rtrn//'Meta%VarName  = <'//TRIM(Meta%VarName)//'>'       //Rtrn &
         //'Meta%DataType = '//TRIM(IntToChar(Meta%DataType)) //Rtrn &
         //'Meta%Dimension= '//TRIM(IntToChar(Meta%Dimension))//Rtrn &
         //'Meta%DataId   = '//TRIM(IntToChar(Meta%DataId))   //Rtrn &
         //'Meta%Status   = '//TRIM(IntToChar(Meta%Status))
  END FUNCTION MetaChar
  !===============================================================================
  !
  !===============================================================================
  FUNCTION StringLen(String)
    CHARACTER(LEN=*),INTENT(IN) :: String
    INTEGER                     :: StringLen
    StringLen=LEN(TRIM(String))
  END FUNCTION StringLen
  !===============================================================================
  !
  !===============================================================================
  FUNCTION Char2Ints(NC,String)
    CHARACTER(LEN=*),INTENT(IN) :: String
    INTEGER,         INTENT(IN) :: NC
    INTEGER,DIMENSION(NC)       :: Char2Ints
    INTEGER                     :: I
    DO I=1,NC
      Char2Ints(I)=ICHAR(String(I:I))
    ENDDO
  END FUNCTION Char2Ints
  !===============================================================================
  !
  !===============================================================================
  SUBROUTINE WriteIntegerVector(Meta,A)
    TYPE(META_DATA)                    :: Meta
    INTEGER, DIMENSION(Meta%Dimension) :: A
    Meta%Status=HDF5WriteIntegerVector(Meta%DataId,Meta%DataSpc,A)
  END SUBROUTINE WriteIntegerVector

  SUBROUTINE WriteDoubleVector(Meta,A)
    TYPE(META_DATA)                         :: Meta
    REAL(DOUBLE), DIMENSION(Meta%Dimension) :: A
    Meta%Status=HDF5WriteDoubleVector(Meta%DataId,Meta%DataSpc,A)
  END SUBROUTINE WriteDoubleVector

  SUBROUTINE ReadIntegerVector(Meta,A)
    TYPE(META_DATA)                    :: Meta
    INTEGER, DIMENSION(Meta%Dimension) :: A
    Meta%Status=HDF5ReadIntegerVector(Meta%DataId,Meta%DataSpc,A)
  END SUBROUTINE ReadIntegerVector

  SUBROUTINE ReadDoubleVector(Meta,A)
    TYPE(META_DATA)                         :: Meta
    REAL(DOUBLE), DIMENSION(Meta%Dimension) :: A
    Meta%Status=HDF5ReadDoubleVector(Meta%DataId,Meta%DataSpc,A)
  END SUBROUTINE ReadDoubleVector
  !===============================================================================
  !
  !===============================================================================
  SUBROUTINE Get_INT_SCLR(A,VarName,Tag_O)
    INTEGER,                  INTENT(INOUT) :: A
    CHARACTER(LEN=*),         INTENT(IN)    :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
    INTEGER,DIMENSION(1)                    :: B
    TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_INT_SCLR", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_INT_SCLR", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,1,.FALSE.)
      CALL OpenData(Meta)
      CALL ReadIntegerVector(Meta,B)
      CALL CloseData(Meta)
      A=B(1)
#ifdef PARALLEL
    ENDIF
    IF(InParallel)CALL Bcast(A)
#endif
  END SUBROUTINE Get_INT_SCLR
  !-------------------------------------------------------------------------------


  SUBROUTINE Get_INT_VECT(A,VarName,Tag_O)
    TYPE(INT_VECT),           INTENT(INOUT) :: A
    CHARACTER(LEN=*),         INTENT(IN)    :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
    TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT) THEN
#endif
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_INT_VECT", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_INT_VECT", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,SIZE(A%I,1),.FALSE.)
      CALL OpenData(Meta)
      CALL ReadIntegerVector(Meta,A%I)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
    IF(InParallel)CALL Bcast(A)
#endif
  END SUBROUTINE Get_INT_VECT
  !-------------------------------------------------------------------------------

  SUBROUTINE Get_INTC(A,VarName,Tag_O)
    TYPE(INTC),               INTENT(INOUT) :: A
    CHARACTER(LEN=*),         INTENT(IN)    :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
    !
    CALL Get(A%N,                'NIntC'//TRIM(VarName),Tag_O=Tag_O)
    IF(A%N==0) RETURN
    IF(AllocQ(A%Alloc)) CALL Delete(A)
    CALL New(A,A%N)
    CALL Get(A%Atoms,            'Atoms'//TRIM(VarName),Tag_O=Tag_O)
    CALL Get(A%Cells,            'Cells'//TRIM(VarName),Tag_O=Tag_O)
    CALL Get(A%Def,                'Def'//TRIM(VarName),Tag_O=Tag_O)
    CALL Get(A%Value,            'Value'//TRIM(VarName),Tag_O=Tag_O)
    CALL Get(A%ConstrValue,'ConstrValue'//TRIM(VarName),Tag_O=Tag_O)
    CALL Get(A%Active,          'Active'//TRIM(VarName),Tag_O=Tag_O)
    CALL Get(A%PredVal,        'PredVal'//TRIM(VarName),Tag_O=Tag_O)
    CALL Get(A%PredGrad,      'PredGrad'//TRIM(VarName),Tag_O=Tag_O)
    CALL Get(A%InvHess,       'PredGrad'//TRIM(VarName),Tag_O=Tag_O)
  END SUBROUTINE Get_INTC

  !-------------------------------------------------------------------------------

  SUBROUTINE Put_INTC(A,VarName,Tag_O)
    TYPE(INTC),               INTENT(IN) :: A
    CHARACTER(LEN=*),         INTENT(IN) :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    !
    !#ifdef PARALLEL
    !       IF(MyId==ROOT)THEN
    !#endif
    CALL Put(A%N,                'NIntC'//TRIM(VarName),Tag_O=Tag_O)
    IF(A%N<=0) RETURN
    CALL Put(A%Def,                'Def'//TRIM(VarName),Tag_O=Tag_O)
    CALL Put(A%Atoms,            'Atoms'//TRIM(VarName),Tag_O=Tag_O)
    CALL Put(A%Cells,            'Cells'//TRIM(VarName),Tag_O=Tag_O)
    CALL Put(A%Value,            'Value'//TRIM(VarName),Tag_O=Tag_O)
    CALL Put(A%ConstrValue,'ConstrValue'//TRIM(VarName),Tag_O=Tag_O)
    CALL Put(A%Active,          'Active'//TRIM(VarName),Tag_O=Tag_O)
    CALL Put(A%PredVal,        'PredVal'//TRIM(VarName),Tag_O=Tag_O)
    CALL Put(A%PredGrad,      'PredGrad'//TRIM(VarName),Tag_O=Tag_O)
    CALL Put(A%InvHess,       'PredGrad'//TRIM(VarName),Tag_O=Tag_O)
    !#ifdef PARALLEL
    !       ENDIF
    !#endif
  END SUBROUTINE Put_INTC
  !
  !---------------------------------------------------------------------
  !

  SUBROUTINE Get_INT_RNK2(A,VarName,Tag_O)
    TYPE(INT_RNK2),           INTENT(INOUT) :: A
    CHARACTER(LEN=*),         INTENT(IN)    :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
    TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT) THEN
#endif
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_INT_RNK2", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_INT_RNK2", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,SIZE(A%I,1)*SIZE(A%I,2),.FALSE.)
      CALL OpenData(Meta)
      CALL ReadIntegerVector(Meta,A%I)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
    IF(InParallel)CALL Bcast(A)
#endif
  END SUBROUTINE Get_INT_RNK2
  !-------------------------------------------------------------------------------


  SUBROUTINE Get_INT_RNK3(A,VarName,Tag_O)
    TYPE(INT_RNK3),           INTENT(INOUT) :: A
    CHARACTER(LEN=*),         INTENT(IN)    :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
    TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT) THEN
#endif
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_INT_RNK3", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_INT_RNK3", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,SIZE(A%I,1)*SIZE(A%I,2)*SIZE(A%I,3),.FALSE.)
      CALL OpenData(Meta)
      CALL ReadIntegerVector(Meta,A%I)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
    IF(InParallel)CALL Bcast(A)
#endif
  END SUBROUTINE Get_INT_RNK3
  !-------------------------------------------------------------------------------


  SUBROUTINE Get_INT_RNK4(A,VarName,Tag_O)
    TYPE(INT_RNK4),           INTENT(INOUT) :: A
    CHARACTER(LEN=*),         INTENT(IN)    :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
    TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT) THEN
#endif
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_INT_RNK4", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_INT_RNK4", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,SIZE(A%I,1)*SIZE(A%I,2)*SIZE(A%I,3)*SIZE(A%I,4),.FALSE.)
      CALL OpenData(Meta)
      CALL ReadIntegerVector(Meta,A%I)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
    IF(InParallel)CALL Bcast(A)
#endif
  END SUBROUTINE Get_INT_RNK4

  SUBROUTINE Get_DBL_SCLR(A,VarName,Stats_O,Tag_O)
    REAL(DOUBLE),             INTENT(INOUT)  :: A
    CHARACTER(LEN=*),         INTENT(IN)     :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)     :: Tag_O
    CHARACTER(LEN=DEFAULT_CHR_LEN)           :: Tag
    INTEGER,DIMENSION(3),OPTIONAL,INTENT(IN) :: Stats_O
    REAL(DOUBLE),DIMENSION(1)                :: B
    TYPE(META_DATA)                          :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      IF(PRESENT(Stats_O).AND.PRESENT(Tag_O)) THEN
        CALL Halt("[Get_DBL_SCLR] either Stats_O or Tag_O")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O=Tag_O,Stats_O=Stats_O),NATIVE_DOUBLE,1,.FALSE.)
      CALL OpenData(Meta)
      CALL ReadDoubleVector(Meta,B)
      CALL CloseData(Meta)
      A=B(1)
#ifdef PARALLEL
    ENDIF
    IF(InParallel)CALL Bcast(A)
#endif
  END SUBROUTINE Get_DBL_SCLR

  SUBROUTINE Get_DBL_VECT(A,VarName,Tag_O)
    TYPE(DBL_VECT),           INTENT(INOUT) :: A
    CHARACTER(LEN=*),         INTENT(IN)    :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
    TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT) THEN
#endif
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_DBL_VECT", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_DBL_VECT", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,SIZE(A%D,1),.FALSE.)
      CALL OpenData(Meta)
      CALL ReadDoubleVector(Meta,A%D)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
    IF(InParallel)CALL Bcast(A)
#endif
  END SUBROUTINE Get_DBL_VECT

  SUBROUTINE Get_DBL_RNK2(A,VarName,Tag_O)
    TYPE(DBL_RNK2),           INTENT(INOUT) :: A
    CHARACTER(LEN=*),         INTENT(IN)    :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
    TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT) THEN
#endif
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_DBL_RNK2", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_DBL_RNK2", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,SIZE(A%D,1)*SIZE(A%D,2),.FALSE.)
      CALL OpenData(Meta)
      CALL ReadDoubleVector(Meta,A%D)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
    IF(InParallel)CALL Bcast(A)
#endif
  END SUBROUTINE Get_DBL_RNK2

  SUBROUTINE Get_DBL_RNK3(A,VarName,Tag_O)
    TYPE(DBL_RNK3),           INTENT(INOUT) :: A
    CHARACTER(LEN=*),         INTENT(IN)    :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
    TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT) THEN
#endif
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_DBL_RNK3", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_DBL_RNK3", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3),.FALSE.)
      CALL OpenData(Meta)
      CALL ReadDoubleVector(Meta,A%D)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
    IF(InParallel)CALL Bcast(A)
#endif
  END SUBROUTINE Get_DBL_RNK3

  SUBROUTINE Get_DBL_RNK4(A,VarName,Tag_O)
    TYPE(DBL_RNK4),           INTENT(INOUT) :: A
    CHARACTER(LEN=*),         INTENT(IN)    :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
    TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT) THEN
#endif
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_DBL_RNK4", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_DBL_RNK4", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3)*SIZE(A%D,4),.FALSE.)
      CALL OpenData(Meta)
      CALL ReadDoubleVector(Meta,A%D)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
    IF(InParallel)CALL Bcast(A)
#endif
  END SUBROUTINE Get_DBL_RNK4

  SUBROUTINE Get_DBL_RNK6(A,VarName,Tag_O)
    TYPE(DBL_RNK6),           INTENT(INOUT) :: A
    CHARACTER(LEN=*),         INTENT(IN)    :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
    TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT) THEN
#endif
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_DBL_RNK6", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_DBL_RNK6", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3)*SIZE(A%D,4)*SIZE(A%D,5)*SIZE(A%D,6),.FALSE.)
      CALL OpenData(Meta)
      CALL ReadDoubleVector(Meta,A%D)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
    IF(InParallel)CALL Bcast(A)
#endif
  END SUBROUTINE Get_DBL_RNK6
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE Get_CHR_SCLR(A,VarName,Tag_O)
    CHARACTER(LEN=*),         INTENT(INOUT) :: A
    CHARACTER(LEN=*),         INTENT(IN)    :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
    INTEGER,DIMENSION(DEFAULT_CHR_LEN)      :: B=ICHAR(' ')
    INTEGER                                 :: I,N
    TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      N=LEN(A)
      IF(N>DEFAULT_CHR_LEN) THEN
        CALL Halt('Static strings overrun in Get_CHR_SCLR')
      ENDIF
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_CHR_SCLR", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_CHR_SCLR", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,DEFAULT_CHR_LEN,.FALSE.)
      CALL OpenData(Meta)
      CALL ReadIntegerVector(Meta,B)
      DO I=1,N
        A(I:I)=CHAR(B(I))
      ENDDO
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
    IF(InParallel)CALL Bcast(A)
#endif
  END SUBROUTINE Get_CHR_SCLR
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE Get_LOG_SCLR(A,VarName,Tag_O)
    LOGICAL,                  INTENT(INOUT) :: A
    CHARACTER(LEN=*),         INTENT(IN)    :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
    INTEGER,DIMENSION(1)                    :: ILog
    TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_LOG_SCLR", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_LOG_SCLR", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,1,.FALSE.)
      CALL OpenData(Meta)
      CALL ReadIntegerVector(Meta,ILog)
      CALL CloseData(Meta)
      IF(ILog(1)==0)THEN
        A=.FALSE.
      ELSEIF(ILog(1)==1)THEN
        A=.TRUE.
      ELSE
        CALL Halt('Error in Get_LOG_SCLR.')
      ENDIF
#ifdef PARALLEL
    ENDIF
    IF(InParallel)CALL Bcast(A)
#endif
  END SUBROUTINE Get_LOG_SCLR
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE Put_INT_SCLR(A,VarName,Tag_O)
    INTEGER,                  INTENT(IN) :: A
    CHARACTER(LEN=*),         INTENT(IN) :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    INTEGER,DIMENSION(1)                 :: B
    TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_INT_SCLR", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_INT_SCLR", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,1,.FALSE.)
      CALL OpenData(Meta,.TRUE.)
      B(1)=A
      CALL WriteIntegerVector(Meta,B)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Put_INT_SCLR

  SUBROUTINE Put_INT_VECT(A,VarName,N_O,Tag_O,Unlimit_O)
    TYPE(INT_VECT),           INTENT(IN) :: A
    CHARACTER(LEN=*),         INTENT(IN) :: VarName
    INTEGER,         OPTIONAL,INTENT(IN) :: N_O
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    LOGICAL,         OPTIONAL,INTENT(IN) :: Unlimit_O
    INTEGER                              :: N
    TYPE(META_DATA)                      :: Meta

#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      IF(PRESENT(N_O))THEN
        N=N_O
      ELSE
        N=SIZE(A%I,1)
      ENDIF
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_INT_VECT", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_INT_VECT", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,N,Unlimit_O)
      CALL OpenData(Meta,.TRUE.)
      CALL WriteIntegerVector(Meta,A%I)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
#endif

  END SUBROUTINE Put_INT_VECT

  SUBROUTINE Put_INT_RNK2(A,VarName,N_O,Tag_O,Unlimit_O)
    TYPE(INT_RNK2),           INTENT(IN) :: A
    CHARACTER(LEN=*),         INTENT(IN) :: VarName
    INTEGER,         OPTIONAL,INTENT(IN) :: N_O
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    LOGICAL,         OPTIONAL,INTENT(IN) :: Unlimit_O
    INTEGER                              :: N
    TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      IF(PRESENT(N_O))THEN
        N=N_O
      ELSE
        N=SIZE(A%I,1)*SIZE(A%I,2)
      ENDIF
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_INT_RNK2", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_INT_RNK2", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,N,Unlimit_O)
      CALL OpenData(Meta,.TRUE.)
      CALL WriteIntegerVector(Meta,A%I)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Put_INT_RNK2

  SUBROUTINE Put_INT_RNK3(A,VarName,N_O,Tag_O,Unlimit_O)
    TYPE(INT_RNK3),           INTENT(IN) :: A
    CHARACTER(LEN=*),         INTENT(IN) :: VarName
    INTEGER,         OPTIONAL,INTENT(IN) :: N_O
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    LOGICAL,         OPTIONAL,INTENT(IN) :: Unlimit_O
    INTEGER                              :: N
    TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      IF(PRESENT(N_O))THEN
        N=N_O
      ELSE
        N=SIZE(A%I,1)*SIZE(A%I,2)*SIZE(A%I,3)
      ENDIF
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_INT_RNK3", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_INT_RNK3", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,N,Unlimit_O)
      CALL OpenData(Meta,.TRUE.)
      CALL WriteIntegerVector(Meta,A%I)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Put_INT_RNK3

  SUBROUTINE Put_INT_RNK4(A,VarName,N_O,Tag_O,Unlimit_O)
    TYPE(INT_RNK4),           INTENT(IN) :: A
    CHARACTER(LEN=*),         INTENT(IN) :: VarName
    INTEGER,         OPTIONAL,INTENT(IN) :: N_O
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    LOGICAL,         OPTIONAL,INTENT(IN) :: Unlimit_O
    INTEGER                              :: N
    TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      IF(PRESENT(N_O))THEN
        N=N_O
      ELSE
        N=SIZE(A%I,1)*SIZE(A%I,2)*SIZE(A%I,3)*SIZE(A%I,4)
      ENDIF
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_INT_RNK4", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_INT_RNK4", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,N,Unlimit_O)
      CALL OpenData(Meta,.TRUE.)
      CALL WriteIntegerVector(Meta,A%I)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Put_INT_RNK4
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE Put_DBL_SCLR(A,VarName,Stats_O,Tag_O)
    REAL(DOUBLE),             INTENT(IN)     :: A
    CHARACTER(LEN=*),         INTENT(IN)     :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)     :: Tag_O
    CHARACTER(LEN=DEFAULT_CHR_LEN)           :: Tag
    INTEGER,DIMENSION(3),OPTIONAL,INTENT(IN) :: Stats_O
    REAL(DOUBLE),DIMENSION(1)                :: B
    TYPE(META_DATA)                          :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      Meta=SetMeta(NameTag(VarName,Tag_O=Tag_O,Stats_O=Stats_O),NATIVE_DOUBLE,1,.FALSE.)
      CALL OpenData(Meta,.TRUE.)
      B(1)=A
      CALL WriteDoubleVector(Meta,B)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Put_DBL_SCLR

  SUBROUTINE Put_DBL_VECT(A,VarName,N_O,Tag_O,Unlimit_O)
    TYPE(DBL_VECT),           INTENT(IN) :: A
    CHARACTER(LEN=*),         INTENT(IN) :: VarName
    INTEGER,         OPTIONAL,INTENT(IN) :: N_O
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    LOGICAL,         OPTIONAL,INTENT(IN) :: Unlimit_O
    INTEGER                              :: N
    TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      IF(PRESENT(N_O))THEN
        N=N_O
      ELSE
        N=SIZE(A%D,1)
      ENDIF
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_DBL_VECT", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_DBL_VECT", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,N,Unlimit_O)
      CALL OpenData(Meta,.TRUE.)
      CALL WriteDoubleVector(Meta,A%D)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Put_DBL_VECT
  !
  SUBROUTINE Put_DBL_RNK2(A,VarName,N_O,Tag_O,Unlimit_O)
    TYPE(DBL_RNK2),           INTENT(IN) :: A
    CHARACTER(LEN=*),         INTENT(IN) :: VarName
    INTEGER,         OPTIONAL,INTENT(IN) :: N_O
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    LOGICAL,         OPTIONAL,INTENT(IN) :: Unlimit_O
    INTEGER                              :: N
    TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      IF(PRESENT(N_O))THEN
        N=N_O
      ELSE
        N=SIZE(A%D,1)*SIZE(A%D,2)
      ENDIF
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_DBL_RNK2", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_DBL_RNK2", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,N,Unlimit_O)
      CALL OpenData(Meta,.TRUE.)
      CALL WriteDoubleVector(Meta,A%D)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Put_DBL_RNK2
  !
  SUBROUTINE Put_DBL_RNK3(A,VarName,N_O,Tag_O,Unlimit_O)
    TYPE(DBL_RNK3),           INTENT(IN) :: A
    CHARACTER(LEN=*),         INTENT(IN) :: VarName
    INTEGER,         OPTIONAL,INTENT(IN) :: N_O
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    LOGICAL,         OPTIONAL,INTENT(IN) :: Unlimit_O
    INTEGER                              :: N
    TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      IF(PRESENT(N_O))THEN
        N=N_O
      ELSE
        N=SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3)
      ENDIF
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_DBL_RNK3", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_DBL_RNK3", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,N,Unlimit_O)
      CALL OpenData(Meta,.TRUE.)
      CALL WriteDoubleVector(Meta,A%D)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Put_DBL_RNK3
  !
  SUBROUTINE Put_DBL_RNK4(A,VarName,N_O,Tag_O,Unlimit_O)
    TYPE(DBL_RNK4),           INTENT(IN) :: A
    CHARACTER(LEN=*),         INTENT(IN) :: VarName
    INTEGER,         OPTIONAL,INTENT(IN) :: N_O
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    LOGICAL,         OPTIONAL,INTENT(IN) :: Unlimit_O
    INTEGER                              :: N
    TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      IF(PRESENT(N_O))THEN
        N=N_O
      ELSE
        N=SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3)*SIZE(A%D,4)
      ENDIF
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_DBL_RNK4", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_DBL_RNK4", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,N,Unlimit_O)
      CALL OpenData(Meta,.TRUE.)
      CALL WriteDoubleVector(Meta,A%D)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Put_DBL_RNK4
  !
  SUBROUTINE Put_DBL_RNK6(A,VarName,N_O,Tag_O,Unlimit_O)
    TYPE(DBL_RNK6),           INTENT(IN) :: A
    CHARACTER(LEN=*),         INTENT(IN) :: VarName
    INTEGER,         OPTIONAL,INTENT(IN) :: N_O
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    LOGICAL,         OPTIONAL,INTENT(IN) :: Unlimit_O
    INTEGER                              :: N
    TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      IF(PRESENT(N_O))THEN
        N=N_O
      ELSE
        N=SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3) &
             *SIZE(A%D,4)*SIZE(A%D,5)*SIZE(A%D,6)
      ENDIF
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_DBL_RNK6", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_DBL_RNK6", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,N,Unlimit_O)
      CALL OpenData(Meta,.TRUE.)
      CALL WriteDoubleVector(Meta,A%D)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Put_DBL_RNK6
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE Put_CHR_SCLR(A,VarName,Tag_O)
    CHARACTER(LEN=*),         INTENT(IN) :: A
    CHARACTER(LEN=*),         INTENT(IN) :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    INTEGER,DIMENSION(DEFAULT_CHR_LEN)   :: B=ICHAR(' ')
    INTEGER                              :: I,N
    TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      N=LEN(A)
      IF(N>DEFAULT_CHR_LEN)THEN
        CALL Halt('Static strings overrun in Put_CHR_SCLR')
      ENDIF
      DO I=1,N
        B(I)=ICHAR(A(I:I))
      ENDDO
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_CHR_SCLR", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_CHR_SCLR", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32, &
           DEFAULT_CHR_LEN,.FALSE.)
      CALL OpenData(Meta,.TRUE.)
      CALL WriteIntegerVector(Meta,B)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Put_CHR_SCLR
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE Put_LOG_SCLR(A,VarName,Tag_O)
    LOGICAL,                  INTENT(IN) :: A
    CHARACTER(LEN=*),         INTENT(IN) :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    INTEGER,DIMENSION(1)                 :: ILog
    TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      ILog(1)=0
      IF(A)ILog(1)=1
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_LOG_SCLR", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_LOG_SCLR", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,1,.FALSE.)
      CALL OpenData(Meta,.TRUE.)
      CALL WriteIntegerVector(Meta,ILog)
      CALL CloseData(Meta)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Put_LOG_SCLR
  !-------------------------------------------------------------------------------
  !     Get a basis set
  !-------------------------------------------------------------------------------
  SUBROUTINE Get_BSET(BS,Tag_O)
    TYPE(BSET)                           :: BS
    CHARACTER(LEN=*),OPTIONAL            :: Tag_O
    IF(AllocQ(BS%Alloc))CALL Delete(BS)
    CALL Get(BS%NAtms,'natoms')
    CALL Get(BS%BName,'bsetname',Tag_O=Tag_O)
    CALL Get(BS%NBasF,'nbasf',Tag_O=Tag_O)
    CALL Get(BS%NKind,'nkind',Tag_O=Tag_O)
    CALL Get(BS%NCtrt,'nctrt',Tag_O=Tag_O)
    CALL Get(BS%NPrim,'nprim',Tag_O=Tag_O)
    CALL Get(BS%NASym,'nasym',Tag_O=Tag_O)
    CALL Get(BS%LMNLen,'lmnlen',Tag_O=Tag_O)
    CALL Get(BS%HasECPs,'hasecps',Tag_O=Tag_O)
    CALL Get(BS%MxProjL,'maxprojectorell',Tag_O=Tag_O)
    CALL Get(BS%Typ1Fnk,'numberoftype1primitives',Tag_O=Tag_O)
    CALL Get(BS%Typ2Fnk,'numberoftype2primitives',Tag_O=Tag_O)
    CALL New(BS)
    CALL Get(BS%Kinds,'kind',Tag_O=Tag_O)
    CALL Get(BS%NCFnc,'ncfuncs',Tag_O=Tag_O)
    CALL Get(BS%NPFnc,'npfuncs',Tag_O=Tag_O)
    CALL Get(BS%ASymm,'angsym',Tag_O=Tag_O)
    CALL Get(BS%Expnt,'exponent',Tag_O=Tag_O)
    CALL Get(BS%CCoef,'ccoefcnt',Tag_O=Tag_O)
    CALL Get(BS%BFKnd,'basfperkind',Tag_O=Tag_O)
    CALL Get(BS%LStrt,'lstart',Tag_O=Tag_O)
    CALL Get(BS%LStop,'lstop',Tag_O=Tag_O)
    CALL Get(BS%LxDex,'lxdex',Tag_O=Tag_O)
    CALL Get(BS%LyDex,'lydex',Tag_O=Tag_O)
    CALL Get(BS%LzDex,'lzdex',Tag_O=Tag_O)
    IF(BS%HasECPs)THEN
      CALL Get(BS%NTyp1PF,'typeoneprimitives',Tag_O=Tag_O)
      CALL Get(BS%NTyp2PF,'typetwoprimitives',Tag_O=Tag_O)
      CALL Get(BS%ProjEll,'projectorell',Tag_O=Tag_O)
      CALL Get(BS%Typ1Ell,'typeoneell',Tag_O=Tag_O)
      CALL Get(BS%Typ2Ell,'typetwoell',Tag_O=Tag_O)
      CALL Get(BS%Typ1Exp,'typeoneexponents',Tag_O=Tag_O)
      CALL Get(BS%Typ2Exp,'typetwoexponents',Tag_O=Tag_O)
      CALL Get(BS%Typ1CCo,'typeonecoefficients',Tag_O=Tag_O)
      CALL Get(BS%Typ2CCo,'typetwocoefficients',Tag_O=Tag_O)
    ENDIF
  END SUBROUTINE Get_BSET
  !-------------------------------------------------------------------------------
  !     Put a  basis set
  !
  SUBROUTINE Put_BSET(BS,Tag_O)
    TYPE(BSET),               INTENT(IN) :: BS
    CHARACTER(LEN=*),OPTIONAL            :: Tag_O
    CALL Put(BS%NAtms,'natoms')
    CALL Put(BS%BName,'bsetname',Tag_O=Tag_O)
    CALL Put(BS%NBasF,'nbasf',Tag_O=Tag_O)
    CALL Put(BS%NKind,'nkind',Tag_O=Tag_O)
    CALL Put(BS%NCtrt,'nctrt',Tag_O=Tag_O)
    CALL Put(BS%NPrim,'nprim',Tag_O=Tag_O)
    CALL Put(BS%NASym,'nasym',Tag_O=Tag_O)
    CALL Put(BS%LMNLen,'lmnlen',Tag_O=Tag_O)
    CALL Put(BS%Kinds,'kind',Tag_O=Tag_O)
    CALL Put(BS%BFKnd,'basfperkind',Tag_O=Tag_O)
    CALL Put(BS%NCFnc,'ncfuncs',Tag_O=Tag_O)
    CALL Put(BS%NPFnc,'npfuncs',Tag_O=Tag_O)
    CALL Put(BS%ASymm,'angsym',Tag_O=Tag_O)
    CALL Put(BS%Expnt,'exponent',Tag_O=Tag_O)
    CALL Put(BS%CCoef,'ccoefcnt',Tag_O=Tag_O)
    CALL Put(BS%LStrt,'lstart',Tag_O=Tag_O)
    CALL Put(BS%LStop,'lstop',Tag_O=Tag_O)
    CALL Put(BS%LxDex,'lxdex',Tag_O=Tag_O)
    CALL Put(BS%LyDex,'lydex',Tag_O=Tag_O)
    CALL Put(BS%LzDex,'lzdex',Tag_O=Tag_O)
    CALL Put(BS%HasECPs,'hasecps',Tag_O=Tag_O)
    CALL Put(BS%MxProjL,'maxprojectorell',Tag_O=Tag_O)
    CALL Put(BS%Typ1Fnk,'numberoftype1primitives',Tag_O=Tag_O)
    CALL Put(BS%Typ2Fnk,'numberoftype2primitives',Tag_O=Tag_O)
    IF(BS%HasECPs)THEN
      CALL Put(BS%NTyp1PF,'typeoneprimitives',Tag_O=Tag_O)
      CALL Put(BS%NTyp2PF,'typetwoprimitives',Tag_O=Tag_O)
      CALL Put(BS%ProjEll,'projectorell',Tag_O=Tag_O)
      CALL Put(BS%Typ1Ell,'typeoneell',Tag_O=Tag_O)
      CALL Put(BS%Typ2Ell,'typetwoell',Tag_O=Tag_O)
      CALL Put(BS%Typ1Exp,'typeoneexponents',Tag_O=Tag_O)
      CALL Put(BS%Typ2Exp,'typetwoexponents',Tag_O=Tag_O)
      CALL Put(BS%Typ1CCo,'typeonecoefficients',Tag_O=Tag_O)
      CALL Put(BS%Typ2CCo,'typetwocoefficients',Tag_O=Tag_O)
    ENDIF
  END SUBROUTINE Put_BSET
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  !             Get the Periodic Info
  SUBROUTINE Get_PBCInfo(PBC,Tag_O)
    TYPE(PBCInfo)                         :: PBC
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)  :: Tag_O
    CALL Get(PBC%Dimen     ,'Dimension' ,Tag_O=Tag_O)
    CALL Get(PBC%PFFMaxEll ,'PFFMaxEll' ,Tag_O=Tag_O)
    CALL Get(PBC%PFFWelSep ,'PFFWelSep' ,Tag_O=Tag_O)
    !! Depricated:
    !!                CALL Get(PBC%AtomW     ,'AtomWrap'  ,Tag_O=Tag_O)
    CALL Get(PBC%SuperCell ,'SuperCell' ,Tag_O=Tag_O)
    CALL Get(PBC%InVecForm ,'VectorForm',Tag_O=Tag_O)
    CALL Get(PBC%InAtomCrd ,'AtomicCrd' ,Tag_O=Tag_O)
    CALL Get(PBC%Translate ,'Translate' ,Tag_O=Tag_O)
    CALL Get(PBC%CellVolume,'CellVolume',Tag_O=Tag_O)
    CALL Get(PBC%Epsilon   ,'Epsilon'   ,Tag_O=Tag_O)
    CALL Get(PBC%DipoleFAC ,'DPoleFAC'  ,Tag_O=Tag_O)
    CALL Get(PBC%QupoleFAC ,'QPoleFAC'  ,Tag_O=Tag_O)
    CALL Get(PBC%AutoW     ,'AutoWrap'  ,Tag_O=Tag_O)
    CALL Get(PBC%CellCenter,'CellCenter',Tag_O=Tag_O)
    CALL Get(PBC%TransVec  ,'TransVec'  ,Tag_O=Tag_O)
    CALL Get(PBC%BoxShape  ,'BoxShape'  ,Tag_O=Tag_O)
    CALL Get(PBC%InvBoxSh  ,'InvBoxSh'  ,Tag_O=Tag_O)
    CALL Get(PBC%LatFrc    ,'LatFrc'    ,Tag_O=Tag_O)
  END SUBROUTINE Get_PBCInfo
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  !             Put the Periodic Info

  SUBROUTINE Put_PBCInfo(PBC,Tag_O)
    TYPE(PBCInfo),            INTENT(IN)  :: PBC
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)  :: Tag_O
    CALL Put(PBC%Dimen     ,'Dimension' ,Tag_O=Tag_O)
    CALL Put(PBC%PFFMaxEll ,'PFFMaxEll' ,Tag_O=Tag_O)
    CALL Put(PBC%PFFWelSep ,'PFFWelSep' ,Tag_O=Tag_O)
    !!                Depricated:
    !!                CALL Put(PBC%AtomW     ,'AtomWrap'  ,Tag_O=Tag_O)
    CALL Put(PBC%SuperCell ,'SuperCell' ,Tag_O=Tag_O)
    CALL Put(PBC%InVecForm ,'VectorForm',Tag_O=Tag_O)
    CALL Put(PBC%InAtomCrd ,'AtomicCrd' ,Tag_O=Tag_O)
    CALL Put(PBC%Translate ,'Translate' ,Tag_O=Tag_O)
    CALL Put(PBC%CellVolume,'CellVolume',Tag_O=Tag_O)
    CALL Put(PBC%Epsilon   ,'Epsilon'   ,Tag_O=Tag_O)
    CALL Put(PBC%DipoleFAC ,'DPoleFAC'  ,Tag_O=Tag_O)
    CALL Put(PBC%QupoleFAC ,'QPoleFAC'  ,Tag_O=Tag_O)
    !                WRITE(*,*)' DIPOLEFAC  = ',PBC%DipoleFac

    CALL Put(PBC%AutoW     ,'AutoWrap'  ,Tag_O=Tag_O)
    CALL Put(PBC%CellCenter,'CellCenter',Tag_O=Tag_O)
    !                WRITE(*,*)' CellCenter = ',PBC%CellCenter%D
    CALL Put(PBC%TransVec  ,'TransVec'  ,Tag_O=Tag_O)
    !                WRITE(*,*)' TransVec = ',PBC%TransVec%D
    CALL Put(PBC%BoxShape  ,'BoxShape'  ,Tag_O=Tag_O)
    !                WRITE(*,*)' BoxShape = ',PBC%BoxShape%D
    CALL Put(PBC%InvBoxSh  ,'InvBoxSh'  ,Tag_O=Tag_O)
    !                WRITE(*,*)' InvBoxSh = ',PBC%InvBoxSh%D
    CALL Put(PBC%LatFrc    ,'LatFrc'    ,Tag_O=Tag_O)
  END SUBROUTINE Put_PBCInfo
  !-------------------------------------------------------------------------------
  !     Get some coordinates
  SUBROUTINE Get_CRDS(GM,Tag_O)
    TYPE(CRDS)                           :: GM
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O

    !CALL MondoLog(DEBUG_NONE, "Get_CRDS", "getting CRDS...")

    ! If already allocated, delete old object.
    IF(AllocQ(GM%Alloc)) CALL Delete(GM)

    !-------------------------------------------------------------------------------
    !        Items that should not change with geometry...
    CALL Get(GM%NAtms,'natoms'       ,Tag_O=Tag_O)
    CALL Get(GM%Confg,'configuration',Tag_O=Tag_O)
    CALL Get(GM%NElec,'nel'          ,Tag_O=Tag_O)
    CALL Get(GM%NAlph,'nelalpha'     ,Tag_O=Tag_O)
    CALL Get(GM%NBeta,'nelbeta'      ,Tag_O=Tag_O)
    CALL Get(GM%TotCh,'charge'       ,Tag_O=Tag_O)
    CALL Get(GM%NKind,'nkind'        ,Tag_O=Tag_O)
    CALL Get(GM%InAU, 'inau'         ,Tag_O=Tag_O)

    ! Allocate new object (with the information we just read from file...).
    CALL New(GM)

    !-------------------------------------------------------------------------------
    !        Items that can change with geometry ...

    CALL Get(GM%ETotal    ,'gm_etot'      ,Tag_O=Tag_O)
    CALL Get(GM%ETotalPerSCF, "gm_etotalperscf",  Tag_O=Tag_O)
    CALL Get(GM%Ordrd     ,'reordered'    ,Tag_O=Tag_O)

    CALL Get(GM%AtTyp     ,'atomtype'     ,Tag_O=Tag_O)
    CALL Get(GM%AtNum     ,'atomicnumbers',Tag_O=Tag_O)
    CALL Get(GM%AtNam     ,'atomname'     ,Tag_O=Tag_O)
!!!!                CALL Get(GM%AtMMTyp   ,'mmtype'       ,Tag_O=Tag_O)
    CALL Get(GM%AtMss     ,'atomicmass'   ,Tag_O=Tag_O)
    CALL Get(GM%PBC                       ,Tag_O=Tag_O)
    CALL Get(GM%InCells   ,'incells'      ,Tag_O=Tag_O)
    CALL Get(GM%OvCells   ,'ovcells'      ,Tag_O=Tag_O)


    CALL Get(GM%BndBox    ,'boundingbox'  ,Tag_O=Tag_O)
    CALL Get(GM%CConstrain,'constraints'  ,Tag_O=Tag_O)
    CALL Get(GM%Carts     ,'cartesians'   ,Tag_O=Tag_O)
    CALL Get(GM%Velocity  ,'Velocity'     ,Tag_O=Tag_O)
    CALL Get(GM%Gradients ,'Gradients'    ,Tag_O=Tag_O)
    CALL Get(GM%Fext      ,'Fext'         ,Tag_O=Tag_O)

    !---------Variables we REALLY want to get rid of-------
    CALL Get(GM%BoxCarts  ,'LatticeCoord' ,Tag_O=Tag_O)
    CALL Get(GM%Displ     ,'Displ'        ,Tag_O=Tag_O)

    ! THIS PBC DISPL MAY BE UNCESSESARY.  LETS SEE IF WE CAN MAKEIT DIE!
    !                IF(PRESENT(Tag_O))THEN
    !                   CALL Get(GM%PBCDispl,Tag_O='KAROLYS_PBCDispl'//TRIM(Tag_O))
    !                ELSE
    !                   CALL Get(GM%PBCDispl,Tag_O='KAROLYS_PBCDispl')
    !                ENDIF

    !--------- to here ----------------------
    CALL Get(GM%LatticeOnly,'LatticeOnly' ,Tag_O=Tag_O)
    CALL Get(GM%AltCount  ,'AltCount'     ,Tag_O=Tag_O)

    !CALL MondoLog(DEBUG_NONE, "Get_CRDS", "done getting CRDS.")

  END SUBROUTINE Get_CRDS
  !-------------------------------------------------------------------------------
  !     Put a coordinate set

  SUBROUTINE Put_CRDS(GM,Tag_O)
    TYPE(CRDS),               INTENT(IN) :: GM
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O

    !-------------------------------------------------------------------------------
    !        Items that should not change with geometry...
!!$                IF(PRESENT(Tag_O))THEN
!!$                   WRITE(*,*)'==================================================================='
!!$                   WRITE(*,*)'==================================================================='
!!$                   WRITE(*,*)' PUTTING COORDINATES WITH TAG = ',TAG_O
!!$                   WRITE(*,*)'==================================================================='
!!$                   WRITE(*,*)'==================================================================='
!!$                ENDIF

    !CALL MondoLog(DEBUG_NONE, "Put_CRDS", "putting...")

    CALL Put(GM%NAtms,'natoms'       ,Tag_O=Tag_O)
    CALL Put(GM%Confg,'configuration',Tag_O=Tag_O)
    CALL Put(GM%NElec,'nel'          ,Tag_O=Tag_O)
    CALL Put(GM%NAlph,'nelalpha'     ,Tag_O=Tag_O)
    CALL Put(GM%NBeta,'nelbeta'      ,Tag_O=Tag_O)
    CALL Put(GM%TotCh,'charge'       ,Tag_O=Tag_O)
    CALL Put(GM%NKind,'nkind'        ,Tag_O=Tag_O)
    CALL Put(GM%InAU, 'inau'         ,Tag_O=Tag_O)
    !-------------------------------------------------------------------------------
    !        Items that can change with geometry ...
    CALL Put(GM%ETotal    ,'gm_etot'      ,Tag_O=Tag_O)
    CALL Put(GM%ETotalPerSCF, "gm_etotalperscf",  Tag_O=Tag_O)
    CALL Put(GM%Ordrd     ,'reordered'    ,Tag_O=Tag_O)
    CALL Put(GM%AtTyp     ,'atomtype'     ,Tag_O=Tag_O)
    CALL Put(GM%AtNum     ,'atomicnumbers',Tag_O=Tag_O)
    CALL Put(GM%AtNam     ,'atomname'     ,Tag_O=Tag_O)
!!!!                CALL Put(GM%AtMMTyp   ,'mmtype'       ,Tag_O=Tag_O)
    CALL Put(GM%AtMss     ,'atomicmass'   ,Tag_O=Tag_O)
    CALL Put(GM%BndBox    ,'boundingbox'  ,Tag_O=Tag_O)

    !                WRITE(*,*)'BndBox = ',GM%BndBox%D

    CALL Put(GM%CConstrain,'constraints'  ,Tag_O=Tag_O)
    CALL Put(GM%Carts     ,'cartesians'   ,Tag_O=Tag_O)

    !                WRITE(*,*)' Carts = ',GM%Carts%D

    CALL Put(GM%Velocity  ,'Velocity'     ,Tag_O=Tag_O)
    CALL Put(GM%Gradients ,'Gradients'    ,Tag_O=Tag_O)
    CALL Put(GM%Fext      ,'Fext'         ,Tag_O=Tag_O)
    ! PBC INFO ONE -- THE ONE WE SHOULD KEEP
    CALL Put(GM%PBC                       ,Tag_O=Tag_O)

    ! In case we don't have InCells allocated yet, we allocate something
    ! (admittedly a hack, but for now the best I have...).
    IF(.NOT. AllocQ(GM%InCells%Alloc)) THEN
      CALL New(GM%InCells, 0)
    ENDIF
    CALL Put(GM%InCells   ,'incells'      ,Tag_O=Tag_O)
    IF(.NOT. AllocQ(GM%OvCells%Alloc)) THEN
      CALL New(GM%OvCells, 0)
    ENDIF
    CALL Put(GM%OvCells   ,'ovcells'      ,Tag_O=Tag_O)
    !---------Variables we REALLY want to get rid of-------
    CALL Put(GM%BoxCarts  ,'LatticeCoord' ,Tag_O=Tag_O)

    !                WRITE(*,*)' BoxCarts = ',GM%BoxCarts%D

    CALL Put(GM%Displ     ,'Displ'        ,Tag_O=Tag_O)

    ! LETS SEE IF WE CAN GET RID OF THIS ONE!!
    ! PBC INFO TWO -- THE ONE WE SHOULD GET RID OF:
    !                IF(PRESENT(Tag_O))THEN
    !                   CALL Put(GM%PBCDispl,Tag_O='KAROLYS_PBCDispl'//TRIM(Tag_O))
    !                ELSE
    !                   CALL Put(GM%PBCDispl,Tag_O='KAROLYS_PBCDispl')
    !                ENDIF

    !--------- to here ----------------------
    CALL Put(GM%LatticeOnly,'LatticeOnly' ,Tag_O=Tag_O)
    CALL Put(GM%AltCount  ,'AltCount'     ,Tag_O=Tag_O)

    !                WRITE(*,*)'==================================================================='
    !                WRITE(*,*)'==================================================================='
    !                WRITE(*,*)' DONE DONE DONE PUTTING COORDINATES '
    !                WRITE(*,*)'==================================================================='
    !                WRITE(*,*)'==================================================================='

    !CALL MondoLog(DEBUG_NONE, "Put_CRDS", "done putting.")

  END SUBROUTINE Put_CRDS
  !-------------------------------------------------------------------------------
  !     Get a BCSR matrix

  SUBROUTINE Get_BCSR(A,Name,PFix_O,CheckPoint_O,BCast_O)
    TYPE(BCSR),               INTENT(INOUT) :: A
    CHARACTER(LEN=*),         INTENT(IN)    :: Name
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: PFix_O
    LOGICAL,         OPTIONAL,INTENT(IN)    :: CheckPoint_O
    LOGICAL,         OPTIONAL,INTENT(IN)    :: BCast_O
    REAL(DOUBLE)                            :: Dummy,Chk
    CHARACTER(LEN=DEFAULT_CHR_LEN)          :: FileName
    INTEGER                                 :: I,NSMat,NAtms,NBlks,NNon0,IOS
    LOGICAL                                 :: Exists,LimitsQ
    LOGICAL                                 :: Bcast

    !CALL MondoLog(DEBUG_MAXIMUM, "Get_BCSR", "getting BCSR from "//TRIM(Name))

    IF(PRESENT(BCast_O)) THEN
      Bcast = BCast_O
    ELSE
      Bcast = .FALSE.
    ENDIF
#ifdef PARALLEL
    IF(MyId==0)THEN
#endif
      IF(PRESENT(CheckPoint_O))THEN
        IF(CheckPoint_O)THEN
          NSMat=A%NSMat
          NAtms=A%NAtms
          NNon0=A%NNon0
          NBlks=A%NBlks

          CALL Get(NSMat,TRIM(Name)//'%NSMat')
          CALL Get(NAtms,TRIM(Name)//'%NAtms')
          CALL Get(NBlks,TRIM(Name)//'%NBlks')
          CALL Get(NNon0,TRIM(Name)//'%NNon0')

          IF(AllocQ(A%Alloc))THEN
            LimitsQ=                            &
                 (NAtms.GT.SIZE(A%RowPt%I)).OR. &
                 (NBlks.GT.SIZE(A%ColPt%I)).OR. &
                 (NBlks.GT.SIZE(A%BlkPt%I)).OR. &
                 (NNon0.GT.SIZE(A%MTrix%D))

            !LimitsQ=.NOT.                   &
            !     (A%NAtms<=SIZE(A%RowPt%I)).AND. &
            !     (A%NBlks<=SIZE(A%ColPt%I)).AND. &
            !     (A%NBlks<=SIZE(A%BlkPt%I)).AND. &
            !     (A%NNon0<=SIZE(A%MTrix%D))
            IF(LimitsQ)THEN
              CALL Delete(A)
              CALL New(A,(/NAtms,NBlks,NNon0/),NSMat_O=NSMat)
            ENDIF
          ELSE
            CALL New(A,(/NAtms,NBlks,NNon0/),NSMat_O=NSMat)
          ENDIF

          CALL Get(A%RowPt,TRIM(Name)//'%RowPt')
          CALL Get(A%ColPt,TRIM(Name)//'%ColPt')
          CALL Get(A%BlkPt,TRIM(Name)//'%BlkPt')
          CALL Get(A%MTrix,TRIM(Name)//'%MTrix')
#ifdef PARALLEL
          IF(Bcast) THEN
            CALL BcastBCSR(A)
          ENDIF
#endif
          RETURN
        ENDIF
      ENDIF

      IF(PRESENT(PFix_O))THEN
        FileName=TRIM(Name)//TRIM(PFix_O)
      ELSE
        FileName=Name
      ENDIF
      INQUIRE(FILE=TRIM(FileName),EXIST=Exists)
      IF(Exists)THEN
#ifdef FORMATTED
        OPEN(UNIT=Seq,FILE=FileName,STATUS='OLD', &
             FORM='FORMATTED',ACCESS='SEQUENTIAL')
#else
        OPEN(UNIT=Seq,FILE=FileName,STATUS='OLD', &
             FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
#endif
      ELSE
        CALL Halt('Get_BCSR could not find '//TRIM(FileName))
      ENDIF

#ifdef FORMATTED
      READ(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)NSMat,NAtms,NNon0,NBlks
      INCLUDE 'Formats.Inc'
#else
      READ(UNIT=Seq,Err=1,IOSTAT=IOS)NSMat,NAtms,NNon0,NBlks
#endif
      IF(AllocQ(A%Alloc))THEN
        IF(NSMat.GT.A%NSMat) THEN
          CALL Delete(A)
          CALL New(A,NSMat_O=NSMat)
        ENDIF
        LimitsQ=                            &
             (NAtms.GT.SIZE(A%RowPt%I)).OR. &
             (NBlks.GT.SIZE(A%ColPt%I)).OR. &
             (NBlks.GT.SIZE(A%BlkPt%I)).OR. &
             (NNon0.GT.SIZE(A%MTrix%D))
        IF(LimitsQ)THEN
          CALL MondoLog(DEBUG_MAXIMUM, "Get_BCSR", &
               'reallocate the matrix A%NSMat.EQ.NSMat = ' &
               //TRIM(LogicalToChar(A%NSMat.EQ.NSMat)))
          CALL Delete(A)
          CALL New(A,(/NAtoms,NBlks,NNon0/),NSMat_O=NSMat)
        ELSE
          A%NSMat=NSMat
          A%NAtms=NAtms
          A%NBlks=NBlks
          A%NNon0=NNon0
        ENDIF
      ELSE
        CALL New(A,(/NAtms,NBlks,NNon0/),NSMat_O=NSMat)
      ENDIF

#ifdef FORMATTED
      READ(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%RowPt%I(I),I=1,NAtoms+1)
      READ(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%ColPt%I(I),I=1,A%NBlks)
      READ(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%BlkPt%I(I),I=1,A%NBlks)
      READ(UNIT=Seq,FMT=66,Err=1,IOSTAT=IOS)(Dummy,i=1,A%NBlks)
      READ(UNIT=Seq,FMT=66,Err=1,IOSTAT=IOS)(A%MTrix%D(I),I=1,A%NNon0)
#else
      READ(UNIT=Seq,Err=2,IOSTAT=IOS)(A%RowPt%I(I),I=1,NAtoms+1)
      READ(UNIT=Seq,Err=3,IOSTAT=IOS)(A%ColPt%I(I),I=1,A%NBlks)
      READ(UNIT=Seq,Err=4,IOSTAT=IOS)(A%BlkPt%I(I),I=1,A%NBlks)
      READ(UNIT=Seq,Err=5,IOSTAT=IOS)(Dummy,I=1,A%NBlks)
      READ(UNIT=Seq,Err=6,IOSTAT=IOS)(A%MTrix%D(I),I=1,A%NNon0)
#endif

      CLOSE(UNIT=Seq,STATUS='KEEP')
#ifdef PARALLEL
    ENDIF
    IF(Bcast) THEN
      CALL BcastBCSR(A)
    ENDIF
#endif

!!$    Chk=Zero
!!$    DO I=1,MIN(NBasF**2,A%NNon0)
!!$       Chk=Chk+A%MTrix%D(I)*A%Mtrix%D(I)
!!$    ENDDO
!!$    Chk=SQRT(Chk)
!!$    CALL MondoLog(DEBUG_MAXIMUM, "Get_BCSR", "getting BCSR from "//TRIM(Name)//" Check = "//TRIM(DblToChar(Chk)))


    RETURN
1   CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Get_BCSR:1340')
2   CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Get_BCSR reading RowPt')
3   CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Get_BCSR reading ColPt')
4   CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Get_BCSR reading BlkPt')
5   CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Get_BCSR reading Dummy')
6   CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Get_BCSR reading MTrix')
  END SUBROUTINE Get_BCSR


  !-------------------------------------------------------------------------------
#ifdef PARALLEL
  SUBROUTINE BcastBCSR(A)
    TYPE(BCSR) :: A
    INTEGER :: IErr,NSMat,NAtms,NBlks,NNon0,i
    LOGICAL :: LimitsQ
    NSMat = A%NSMat
    NAtms = A%NAtms
    NBlks = A%NBlks
    NNon0 = A%NNon0

    CALL Bcast(NSMat)
    CALL Bcast(NAtms)
    CALL Bcast(NBlks)
    CALL Bcast(NNon0)

    IF(AllocQ(A%Alloc))THEN
      LimitsQ=                            &
           (NAtms.GT.SIZE(A%RowPt%I)).OR. &
           (NBlks.GT.SIZE(A%ColPt%I)).OR. &
           (NBlks.GT.SIZE(A%BlkPt%I)).OR. &
           (NNon0.GT.SIZE(A%MTrix%D))
      IF(LimitsQ.AND.MyID.EQ.0) CALL Halt('BcastBCSR: Something wrong there 1')
      IF(LimitsQ)THEN
        CALL Delete(A,OnAll_O=.TRUE.)
        CALL New(A,(/NAtms,NBlks,NNon0/),NSMat_O=NSMat,OnAll_O=.TRUE.)
      ENDIF
    ELSE
      IF(MyID.EQ.0) CALL Halt('BcastBCSR: Something wrong there 2')
      CALL New(A,(/NAtms,NBlks,NNon0/),NSMat_O=NSMat,OnAll_O=.TRUE.)
    ENDIF
    CALL Bcast(A%RowPt,N_O=NAtoms+1)
    CALL Bcast(A%ColPt,N_O=NBlks)
    CALL Bcast(A%BlkPt,N_O=NBlks)
    CALL Bcast(A%MTrix,N_O=NNon0)
  END SUBROUTINE BcastBCSR
#endif
  !-------------------------------------------------------------------------------
  !     Put a BCSR matrix

  SUBROUTINE Put_BCSR(A,Name,PFix_O,CheckPoint_O)
    TYPE(BCSR),               INTENT(IN) :: A
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: PFix_O
    CHARACTER(LEN=*),         INTENT(IN) :: Name
    LOGICAL,         OPTIONAL,INTENT(IN) :: CheckPoint_O
    CHARACTER(LEN=DEFAULT_CHR_LEN)       :: FileName
    LOGICAL                              :: Exists
    INTEGER                              :: I,IOS

!!$    REAL(DOUBLE) :: Chk
!!$    Chk=Zero
!!$    DO I=1,A%NNon0
!!$       Chk=Chk+A%MTrix%D(I)*A%Mtrix%D(I)
!!$    ENDDO
!!$    Chk=SQRT(Chk)
!!$    CALL MondoLog(DEBUG_MAXIMUM, "Put_BCSR", "getting BCSR from "//TRIM(Name)//" Check = "//TRIM(DblToChar(Chk)))

#ifdef PARALLEL
    IF(MyId==0)THEN
#endif
      IF(PRESENT(PFix_O))THEN
        FileName=TRIM(Name)//TRIM(PFix_O)
      ELSE
        FileName=Name
      ENDIF

      IF(PRESENT(CheckPoint_O))THEN
        IF(CheckPoint_O)THEN
          !write(*,*) 'put',A%NSMat
          CALL Put(A%NSMat,TRIM(Name)//'%NSMat')
          !call get(A%NSMat,TRIM(Name)//'%NSMat')
          !write(*,*) 'get',A%NSMat
          CALL Put(A%NAtms,TRIM(Name)//'%NAtms')
          CALL Put(A%NBlks,TRIM(Name)//'%NBlks')
          CALL Put(A%NNon0,TRIM(Name)//'%NNon0')
          CALL Put(A%RowPt,TRIM(Name)//'%RowPt',A%NAtms+1)
#ifdef PARALLEL_CLONES
          CALL Put(A%ColPt,TRIM(Name)//'%ColPt',A%NBlks)!,Unlimit_O=.TRUE.)
          !CALL Put(A%BlkPt,TRIM(Name)//'%BlkPt',A%NBlks*A%NSMat)!,Unlimit_O=.TRUE.)
          CALL Put(A%BlkPt,TRIM(Name)//'%BlkPt',A%NBlks)!,Unlimit_O=.TRUE.)
          CALL Put(A%MTrix,TRIM(Name)//'%MTrix',A%NNon0)!,Unlimit_O=.TRUE.)
#else
          CALL Put(A%ColPt,TRIM(Name)//'%ColPt',A%NBlks,Unlimit_O=.TRUE.)
          !CALL Put(A%BlkPt,TRIM(Name)//'%BlkPt',A%NBlks*A%NSMat,Unlimit_O=.TRUE.)
          CALL Put(A%BlkPt,TRIM(Name)//'%BlkPt',A%NBlks,Unlimit_O=.TRUE.)
          CALL Put(A%MTrix,TRIM(Name)//'%MTrix',A%NNon0,Unlimit_O=.TRUE.)
#endif
          RETURN
        ENDIF
      ENDIF

      ! WRITE TO BINARY FILE
      INQUIRE(FILE=FileName,EXIST=Exists)
      IF(Exists)THEN
#if FORMATTED
        OPEN(UNIT=Seq,FILE=FileName,STATUS='REPLACE', &
             FORM='FORMATTED',ACCESS='SEQUENTIAL')
#else
        OPEN(UNIT=Seq,FILE=FileName,STATUS='REPLACE', &
             FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
#endif
      ELSE
#if FORMATTED
        OPEN(UNIT=Seq,FILE=FileName,STATUS='NEW', &
             FORM='FORMATTED',ACCESS='SEQUENTIAL')
#else
        OPEN(UNIT=Seq,FILE=FileName,STATUS='NEW', &
             FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
#endif
      ENDIF
#ifdef FORMATTED
      WRITE(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)A%NSMat,NAtoms,A%NNon0,A%NBlks
      WRITE(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%RowPt%I(i),i=1,NAtoms+1)
      WRITE(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%ColPt%I(i),i=1,A%NBlks)
      !WRITE(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%BlkPt%I(i),i=1,A%NBlks*A%NSMat)
      WRITE(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%BlkPt%I(i),i=1,A%NBlks)
      WRITE(UNIT=Seq,FMT=66,Err=1,IOSTAT=IOS)(BIG_DBL,i=1,A%NBlks)
      WRITE(UNIT=Seq,FMT=66,Err=1,IOSTAT=IOS)(A%MTrix%D(i),i=1,A%NNon0)
      INCLUDE 'Formats.Inc'
#else
      !write(*,*) 'put',A%NSMat
      WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)A%NSMat,NAtoms,A%NNon0,A%NBlks
      WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%RowPt%I(I),I=1,NAtoms+1)
      WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%ColPt%I(I),I=1,A%NBlks)
      !WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%BlkPt%I(I),I=1,A%NBlks*A%NSMat)
      WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%BlkPt%I(I),I=1,A%NBlks)
      WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(BIG_DBL,I=1,A%NBlks)
      WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%MTrix%D(I),I=1,A%NNon0)
#endif
      CLOSE(UNIT=Seq,STATUS='KEEP')
#ifdef PARALLEL
    ENDIF
#endif

    RETURN
1   CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Put_BCSR.')
  END SUBROUTINE Put_BCSR
#ifdef PARALLEL
  !#ifdef MPIIO
  !
  !***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***
  !***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***
  !***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***
  !
  ! Get a DBCSR matrix
  SUBROUTINE Get_DBCSR_MPIIO_N_1(A,Name,PFix_O)
    TYPE(DBCSR),              INTENT(INOUT) :: A
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: PFix_O
    CHARACTER(LEN=*),         INTENT(IN)    :: Name
    CHARACTER(LEN=DEFAULT_CHR_LEN)          :: FileName
    TYPE(INT_VECT)                          :: INTArr
    INTEGER,DIMENSION(MPI_STATUS_SIZE)      :: Status
    INTEGER,PARAMETER                       :: INTSIZE=4,DBLSIZE=8
    INTEGER                                 :: IErr,Id,Header,Count,NPrcTmp,NInt,Dims(4)
    INTEGER(KIND=MPI_OFFSET_KIND)           :: OffSet,LOffSet,ArrGOff(NPrc)
    LOGICAL                                 :: Exists,LimitsQ
    !-------------------------------------------------------------------------------
    !
    IF(PRESENT(PFix_O))THEN
      FileName=TRIM(Name)//TRIM(PFix_O)
    ELSE
      FileName=Name
    ENDIF
    INQUIRE(FILE=TRIM(FileName),EXIST=Exists)
    IF(.NOT.Exists) CALL Halt(' Get_DBCSR_MPI_IO could not find '//TRIM(FileName))
    !
    CALL MPI_FILE_OPEN(MONDO_COMM,FileName,MPI_MODE_RDONLY,MPI_INFO_NULL,Id,IErr)
    !
    ! Read NPrc to the file and bcast
    IF(MyID.EQ.0) THEN
      CALL MPI_FILE_READ_AT(Id,0_MPI_OFFSET_KIND,NPrcTmp,1,MPI_INTEGER,Status,IErr)
    ENDIF
    CALL MPI_BCAST(NPrcTmp,1,MPI_INTEGER,0,MONDO_COMM,IErr)
    !
    ! ADD SOMETHING HERE IF WE WANT RESTART WITH DIFF NPrc
    IF(NPrcTmp.NE.NPrc) CALL Halt('Get_DBCSR_MPI_IO: current Nprc not consistent with the one in the file!')
    !
    ! Read the local offset
    CALL MPI_FILE_READ_AT_ALL(Id,INT(INTSIZE,MPI_OFFSET_KIND),ArrGOff(1), &
         &                NPrc*MPI_OFFSET_KIND,MPI_BYTE,Status,IErr)
    !
    ! Read the data
    OffSet=ArrGOff(MyID+1)
    CALL MPI_FILE_READ_AT_ALL(Id,OffSet,Dims(1),4,MPI_INTEGER,Status,IErr)
    !
    ! Check if right size
    CALL CheckAlloc_DBCSR(A,Dims)
    !
    NInt=2*(A%NAtms+1)+2*A%NBlks
    CALL New(INTArr,NInt)
    !
    ! Read the data
    OffSet=OffSet+4*INTSIZE
    CALL MPI_FILE_READ_AT_ALL(Id,OffSet,INTArr%I(1),NInt,MPI_INTEGER,Status,IErr)
    !CALL MPI_GET_COUNT(Status,MPI_INTEGER,Count,IErr)
    !write(*,'(A,I3,A,I6,A,I6)') 'process', MyID,' read 1',Count,' integers, expected ',NInt
    OffSet=OffSet+NInt*INTSIZE
    CALL MPI_FILE_READ_AT_ALL(Id,OffSet,A%MTrix%D(1),A%NNon0  ,MPI_DOUBLE_PRECISION,Status,IErr)
    !CALL MPI_GET_COUNT(Status,MPI_DOUBLE_PRECISION,Count,IErr)
    !write(*,'(A,I3,A,I6,A,I6)') 'process', MyID,' read 2',Count,' realss, expected ',A%NNon0
    !
    ! Unpack the data
    CALL INT_VECT_EQ_INT_VECT(A%NAtms+1,A%RowPt%I(1),INTArr%I(                      1))
    CALL INT_VECT_EQ_INT_VECT(A%NAtms+1,A%GRwPt%I(1),INTArr%I(  (A%NAtms+1)        +1))
    CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,A%ColPt%I(1),INTArr%I(2*(A%NAtms+1)        +1))
    CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,A%BlkPt%I(1),INTArr%I(2*(A%NAtms+1)+A%NBlks+1))
    !
    A%Node=MyId
    !
    CALL Delete(INTArr)
    CALL MPI_FILE_CLOSE(Id,IErr)
    !
  END SUBROUTINE Get_DBCSR_MPIIO_N_1
  !
  ! Put a DBCSR matrix
  SUBROUTINE Put_DBCSR_MPIIO_N_1(A,Name,PFix_O)
    TYPE(DBCSR), INTENT(INOUT)           :: A
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: PFix_O
    CHARACTER(LEN=*),         INTENT(IN) :: Name
    CHARACTER(LEN=DEFAULT_CHR_LEN)       :: FileName
    TYPE(INT_VECT)                       :: INTArr
    INTEGER,DIMENSION(MPI_STATUS_SIZE)   :: Status
    INTEGER,PARAMETER                    :: INTSIZE=4,DBLSIZE=8
    INTEGER                              :: IErr,Id,Dims(4),Count,Info,NInt
    INTEGER(KIND=MPI_OFFSET_KIND)        :: OffSet,Header,ArrGOff(NPrc),LOffSet
    !-------------------------------------------------------------------------------
    !
    IF(PRESENT(PFix_O))THEN
      FileName=TRIM(Name)//TRIM(PFix_O)
    ELSE
      FileName=Name
    ENDIF
    !
    CALL MPI_INFO_CREATE(Info,IErr)
    CALL MPI_INFO_SET(Info,'serialize_open' ,'yes'    ,IErr)
    CALL MPI_INFO_SET(Info,'striping_factor','2'      ,IErr) !2
    CALL MPI_INFO_SET(Info,'striping_unit'  ,'8388608',IErr) !33554432 !8388608!1048576
    CALL MPI_INFO_SET(Info,'start_iodevice' ,'1'      ,IErr) !3
    !
    CALL MPI_FILE_OPEN(MONDO_COMM,FileName,MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                                !&             MPI_INFO_NULL,Id,IErr)
         &             Info,Id,IErr)
    !
    ! Compute header
    Header=NPrc*MPI_OFFSET_KIND+INTSIZE
    CALL GetOffSet_DBCSR(A,ArrGOff,Header)
    !
    IF(MyID.EQ.0) THEN
      ! Write NPrc to the file
      CALL MPI_FILE_WRITE_AT(Id,0_MPI_OFFSET_KIND,NPrc,1,MPI_INTEGER,Status,IErr)
      ! Write the global offset
      CALL MPI_FILE_WRITE_AT(Id,INT(INTSIZE,MPI_OFFSET_KIND),ArrGOff(1), &
           &                 NPrc*MPI_OFFSET_KIND,MPI_BYTE,Status,IErr)
    ENDIF
    !
    ! Pack the data
    NInt=4+2*(A%NAtms+1)+2*A%NBlks
    CALL New(INTArr,NInt)
    INTArr%I(1)=A%NSMat
    INTArr%I(2)=A%NAtms
    INTArr%I(3)=A%NBlks
    INTArr%I(4)=A%NNon0
    CALL INT_VECT_EQ_INT_VECT(A%NAtms+1,INTArr%I(                      5),A%RowPt%I(1))
    CALL INT_VECT_EQ_INT_VECT(A%NAtms+1,INTArr%I(  (A%NAtms+1)        +5),A%GRwPt%I(1))
    CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,INTArr%I(2*(A%NAtms+1)        +5),A%ColPt%I(1))
    CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,INTArr%I(2*(A%NAtms+1)+A%NBlks+5),A%BlkPt%I(1))
    !
    ! Write the data
    OffSet=ArrGOff(MyID+1)
    CALL MPI_FILE_WRITE_AT_ALL(Id,OffSet,INTArr%I(1),NInt,MPI_INTEGER,Status,IErr)
    CALL MPI_GET_COUNT(Status,MPI_INTEGER,Count,IErr)
    IF(Count.NE.NInt) THEN
      CALL Halt('process '//TRIM(IntToChar(MyID))//' write '//TRIM(IntToChar(Count))// &
           &    ' integers, expected '//TRIM(IntToChar(NInt)))
    ENDIF
    !
    OffSet=OffSet+NInt*INTSIZE
    CALL MPI_FILE_WRITE_AT_ALL(Id,OffSet,A%MTrix%D(1),A%NNon0,MPI_DOUBLE_PRECISION,Status,IErr)
    CALL MPI_GET_COUNT(Status,MPI_DOUBLE_PRECISION,Count,IErr)
    IF(Count.NE.A%NNon0) THEN
      CALL Halt('process '//TRIM(IntToChar(MyID))//' write '//TRIM(IntToChar(Count))// &
           &    ' reals, expected '//TRIM(IntToChar(A%NNon0)))
    ENDIF
    !
    CALL Delete(INTArr)
    CALL MPI_FILE_CLOSE(Id,IErr)
    !
  END SUBROUTINE Put_DBCSR_MPIIO_N_1

  SUBROUTINE GetOffSet_DBCSR(A,ArrGOff,Header)
    TYPE(DBCSR)                  :: A
    INTEGER(KIND=MPI_OFFSET_KIND):: ArrGOff(NPrc),ArrLOff(NPrc)
    INTEGER(KIND=MPI_OFFSET_KIND):: Header
    INTEGER,PARAMETER            :: INTSIZE=4,DBLSIZE=8
    INTEGER                      :: THeader,IErr
    !
    ArrGOff=0
    ! Local
    THeader=                   4 *INTSIZE ! A%NSMat,NAtoms,A%NBlks,A%NNon0
    THeader=THeader+2*(A%NAtms+1)*INTSIZE ! A%RowPt%I,A%GRwPt%I
    THeader=THeader+2*(A%NBlks  )*INTSIZE ! A%ColPt%I,A%BlkPt%I
    THeader=THeader+  (A%NNon0  )*DBLSIZE ! A%MTrix%D
    !B%GUpDate=A%GUpDate
    !A%GClPt! ??
    ArrGOff(MyID+1)=THeader
    CALL MPI_ALLGATHER(ArrGOff(MyID+1),MPI_OFFSET_KIND,MPI_BYTE,ArrGOff(1), &
         &             MPI_OFFSET_KIND,MPI_BYTE,MONDO_COMM,IErr)
    CALL CalcGOffSet(ArrGOff,Header)
  END SUBROUTINE GetOffSet_DBCSR

  SUBROUTINE CalcGOffSet(ArrGOff,Header)
    INTEGER(KIND=MPI_OFFSET_KIND):: ArrGOff(NPrc),ArrLOff(NPrc),Header
    INTEGER                      :: I
    ArrLOff=ArrGOff
    ArrGOff(1)=Header
    DO I=1,NPrc-1
      ArrGOff(I+1)=ArrGOff(I)+ArrLOff(I)
    ENDDO
  END SUBROUTINE CalcGOffSet

  SUBROUTINE CheckAlloc_DBCSR(A,Dims)
    TYPE(DBCSR) :: A
    INTEGER     :: Dims(4)
    LOGICAL     :: LimitsQ
    IF(AllocQ(A%Alloc))THEN
      IF(Dims(1).GT.A%NSMat) THEN
        CALL Delete(A)
        CALL New(A,NSMat_O=Dims(1))
      ENDIF
      LimitsQ=                             &
           (Dims(2).GT.SIZE(A%RowPt%I)).OR. &
           (Dims(3).GT.SIZE(A%ColPt%I)).OR. &
           (Dims(3).GT.SIZE(A%BlkPt%I)).OR. &
           (Dims(4).GT.SIZE(A%MTrix%D))
      IF(LimitsQ)THEN
        write(*,*)'In CheckAlloc_DBCSR Reallocate the matrix A%NSMat.EQ.NSMat=',A%NSMat.EQ.Dims(1)
        CALL Delete(A)
        CALL New(A,(/Dims(2),Dims(3),Dims(4)/),NSMat_O=Dims(1))
      ELSE
        A%NSMat=Dims(1)
        A%NAtms=Dims(2)
        A%NBlks=Dims(3)
        A%NNon0=Dims(4)
      ENDIF
    ELSE
      CALL New(A,(/Dims(2),Dims(3),Dims(4)/),NSMat_O=Dims(1))
      A%NSMat=Dims(1)
      A%NAtms=Dims(2)
      A%NBlks=Dims(3)
      A%NNon0=Dims(4)
    ENDIF
  END SUBROUTINE CheckAlloc_DBCSR

  SUBROUTINE Put_DBCSR_N_1(A,Name,PFix_O)
    TYPE(DBCSR), INTENT(INOUT)           :: A
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: PFix_O
    CHARACTER(LEN=*),         INTENT(IN) :: Name
    CHARACTER(LEN=DEFAULT_CHR_LEN)       :: FileName
    TYPE(INT_VECT)                       :: INTArr
    TYPE(DBL_VECT)                       :: DBLArr
    INTEGER,PARAMETER                    :: INTSIZE=4,DBLSIZE=8
    INTEGER                              :: IErr,IOS,Dims(4),Count,Info,I,IPrc,NInt,NDbl
    INTEGER(KIND=MPI_OFFSET_KIND)        :: OffSet,Header,ArrGOff(NPrc),LOffSet
    LOGICAL                              :: Exists
    !-------------------------------------------------------------------------------
    !
    IF(PRESENT(PFix_O))THEN
      FileName=TRIM(Name)//TRIM(PFix_O)
    ELSE
      FileName=Name
    ENDIF
    !
    IF(MyID.EQ.0) THEN
      INQUIRE(FILE=FileName,EXIST=Exists)
      IF(Exists)THEN
        OPEN(UNIT=Seq,FILE=FileName,STATUS='REPLACE', &
             FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      ELSE
        OPEN(UNIT=Seq,FILE=FileName,STATUS='NEW', &
             FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      ENDIF

      !Write NPrc
      WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)NPrc

      !Write the ROOT bcsr matrix
      WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)A%NSMat,A%NAtms,A%NBlks,A%NNon0
      WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%RowPt%I(I),I=1,A%NAtms+1)
      WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%GRwPt%I(I),I=1,A%NAtms+1)
      WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%ColPt%I(I),I=1,A%NBlks)
      WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%BlkPt%I(I),I=1,A%NBlks)
      WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%MTrix%D(I),I=1,A%NNon0)
      CALL New(INTArr,2*(A%NAtms+1)+2*A%NBlks)
      CALL New(DBLArr,A%NNon0)
      !
      DO IPrc=1,NPrc-1
        CALL MPI_RECV(Dims(1),4,MPI_INTEGER,IPrc,101,MONDO_COMM,MPI_STATUS_IGNORE,IErr)
        NInt=2*(Dims(2)+1)+2*Dims(3)
        NDbl=Dims(4)
        IF(SIZE(INTArr%I).LT.NInt) THEN
          CALL Delete(INTArr)
          CALL New(INTArr,NInt)
        ENDIF
        IF(SIZE(DBLArr%D).LT.NDbl) THEN
          CALL Delete(DBLArr)
          CALL New(DBLArr,NDbl)
        ENDIF
        CALL MPI_RECV(INTArr%I(1),NInt,MPI_INTEGER         ,IPrc,201,MONDO_COMM,MPI_STATUS_IGNORE,IErr)
        CALL MPI_RECV(DBLArr%D(1),NDbl,MPI_DOUBLE_PRECISION,IPrc,301,MONDO_COMM,MPI_STATUS_IGNORE,IErr)
        WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(Dims(I),I=1,4)
        WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(INTArr%I(I),I=1,NInt)
        WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(DBLArr%D(I),I=1,NDbl)
      ENDDO
      CALL Delete(INTArr)
      CALL Delete(DBLArr)
      CLOSE(UNIT=Seq,STATUS='KEEP')
    ELSE
      Dims(1)=A%NSMat;Dims(2)=A%NAtms;Dims(3)=A%NBlks;Dims(4)=A%NNon0
      NInt=2*(Dims(2)+1)+2*Dims(3)
      NDbl=Dims(4)
      CALL New(INTArr,NInt)
      !CALL INT_VECT_EQ_INT_VECT(A%NAtms+1,INTArr%I(                1),A%RowPt%I(1))
      !CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,INTArr%I(A%NAtms        +2),A%ColPt%I(1))
      !CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,INTArr%I(A%NAtms+A%NBlks+2),A%BlkPt%I(1))
      CALL INT_VECT_EQ_INT_VECT(A%NAtms+1,INTArr%I(                      1),A%RowPt%I(1))
      CALL INT_VECT_EQ_INT_VECT(A%NAtms+1,INTArr%I(  (A%NAtms+1)        +1),A%GRwPt%I(1))
      CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,INTArr%I(2*(A%NAtms+1)        +1),A%ColPt%I(1))
      CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,INTArr%I(2*(A%NAtms+1)+A%NBlks+1),A%BlkPt%I(1))
      CALL MPI_SEND(     Dims(1),   4,MPI_INTEGER         ,0,101,MONDO_COMM,IErr)
      CALL MPI_SEND( INTArr%I(1),NInt,MPI_INTEGER         ,0,201,MONDO_COMM,IErr)
      CALL MPI_SEND(A%MTrix%D(1),NDbl,MPI_DOUBLE_PRECISION,0,301,MONDO_COMM,IErr)
      CALL Delete(INTArr)
    ENDIF

    RETURN
1   CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Put_DBCSR_N_1.')
  END SUBROUTINE Put_DBCSR_N_1

  SUBROUTINE Get_DBCSR_N_1(A,Name,PFix_O)
    TYPE(DBCSR),              INTENT(INOUT) :: A
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: PFix_O
    CHARACTER(LEN=*),         INTENT(IN)    :: Name
    CHARACTER(LEN=DEFAULT_CHR_LEN)          :: FileName
    TYPE(INT_VECT)                          :: INTArr
    TYPE(DBL_VECT)                          :: DBLArr
    INTEGER,DIMENSION(MPI_STATUS_SIZE)      :: Status
    INTEGER,PARAMETER                       :: INTSIZE=4,DBLSIZE=8
    INTEGER                                 :: IErr,IOS,Id,NInt,NDbl,Count,I,IPrc,NPrcTmp,Dims(4)
    LOGICAL                                 :: Exists
    !-------------------------------------------------------------------------------
    !
    IF(PRESENT(PFix_O))THEN
      FileName=TRIM(Name)//TRIM(PFix_O)
    ELSE
      FileName=Name
    ENDIF
    !
    IF(MyID.EQ.0) THEN
      INQUIRE(FILE=TRIM(FileName),EXIST=Exists)
      IF(.NOT.Exists) CALL Halt(' Get_DBCSR_N_1 could not find '//TRIM(FileName))
      !
      !
      OPEN(UNIT=Seq,FILE=FileName,STATUS='OLD', &
           FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      !
      ! Read NPrc to the file and bcast
      READ(UNIT=Seq,Err=1,IOSTAT=IOS)NPrcTmp
      READ(UNIT=Seq,Err=1,IOSTAT=IOS)(Dims(I),I=1,4)
      IF(NPrcTmp.NE.NPrc) CALL Halt('Get_DBCSR_N_1: current Nprc not consistent with the one in the file!')
      !
      ! Check if right size
      CALL CheckAlloc_DBCSR(A,Dims)
      !
      ! Read the ROOT bcsr matrix
      READ(UNIT=Seq,Err=1,IOSTAT=IOS)(A%RowPt%I(I),I=1,A%NAtms+1)
      READ(UNIT=Seq,Err=1,IOSTAT=IOS)(A%GRwPt%I(I),I=1,A%NAtms+1)
      READ(UNIT=Seq,Err=1,IOSTAT=IOS)(A%ColPt%I(I),I=1,A%NBlks)
      READ(UNIT=Seq,Err=1,IOSTAT=IOS)(A%BlkPt%I(I),I=1,A%NBlks)
      READ(UNIT=Seq,Err=1,IOSTAT=IOS)(A%MTrix%D(I),I=1,A%NNon0)
      CALL New(INTArr,2*(A%NAtms+1)+2*A%NBlks)
      CALL New(DBLArr,A%NNon0)
      !
      DO IPrc=1,NPrc-1
        READ(UNIT=Seq,Err=1,IOSTAT=IOS)(Dims(I),I=1,4)
        NInt=2*(Dims(2)+1)+2*Dims(3)
        NDbl=Dims(4)
        IF(SIZE(INTArr%I).LT.NInt) THEN
          CALL Delete(INTArr)
          CALL New(INTArr,NInt)
        ENDIF
        IF(SIZE(DBLArr%D).LT.NDbl) THEN
          CALL Delete(DBLArr)
          CALL New(DBLArr,NDbl)
        ENDIF
        READ(UNIT=Seq,Err=1,IOSTAT=IOS)(INTArr%I(I),I=1,NInt)
        READ(UNIT=Seq,Err=1,IOSTAT=IOS)(DBLArr%D(I),I=1,NDbl)
        !
        CALL MPI_SEND(    Dims(1),   4,MPI_INTEGER         ,IPrc,101,MONDO_COMM,IErr)
        CALL MPI_SEND(INTArr%I(1),NInt,MPI_INTEGER         ,IPrc,201,MONDO_COMM,IErr)
        CALL MPI_SEND(DBLArr%D(1),NDbl,MPI_DOUBLE_PRECISION,IPrc,301,MONDO_COMM,IErr)
      ENDDO
      CALL Delete(INTArr)
      CALL Delete(DBLArr)
      CLOSE(UNIT=Seq,STATUS='KEEP')
    ELSE
      CALL MPI_RECV(Dims(1),4,MPI_INTEGER,0,101,MONDO_COMM,MPI_STATUS_IGNORE,IErr)
      CALL CheckAlloc_DBCSR(A,Dims)
      NInt=2*(Dims(2)+1)+2*Dims(3)
      NDbl=Dims(4)
      CALL New(INTArr,NInt)
      CALL MPI_RECV( INTArr%I(1),NInt,MPI_INTEGER         ,0,201,MONDO_COMM,MPI_STATUS_IGNORE,IErr)
      CALL MPI_RECV(A%MTrix%D(1),NDbl,MPI_DOUBLE_PRECISION,0,301,MONDO_COMM,MPI_STATUS_IGNORE,IErr)
      CALL INT_VECT_EQ_INT_VECT(A%NAtms+1,A%RowPt%I(1),INTArr%I(                      1))
      CALL INT_VECT_EQ_INT_VECT(A%NAtms+1,A%GRwPt%I(1),INTArr%I(  (A%NAtms+1)        +1))
      CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,A%ColPt%I(1),INTArr%I(2*(A%NAtms+1)        +1))
      CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,A%BlkPt%I(1),INTArr%I(2*(A%NAtms+1)+A%NBlks+1))
      !CALL INT_VECT_EQ_INT_VECT(A%NAtms+1,A%RowPt%I(1),INTArr%I(                1))
      !CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,A%ColPt%I(1),INTArr%I(A%NAtms        +2))
      !CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,A%BlkPt%I(1),INTArr%I(A%NAtms+A%NBlks+2))
      CALL Delete(INTArr)
    ENDIF
    A%Node=MyId
    !
    RETURN
1   CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Get_DBCSR_N_1.')
  END SUBROUTINE Get_DBCSR_N_1


  !***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***
  !***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***
  !***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***MPI-IO***
  !
  !-------------------------------------------------------------------------------
  !     Get a DBCSR matrix

  SUBROUTINE Get_DBCSR(A,Name,PFix_O,CheckPoint_O)
    TYPE(DBCSR),              INTENT(INOUT) :: A
    LOGICAL,         OPTIONAL,INTENT(IN)    :: CheckPoint_O
    CHARACTER(LEN=*),         INTENT(IN)    :: Name
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: PFix_O
    TYPE(BCSR)                              :: B
#ifdef PARALLEL
    LOGICAL                                 :: InParTemp
#endif
    !-------------------------------------------------------------------------------
#ifdef PARALLEL
    IF(PRESENT(Checkpoint_O))THEN
      InParTemp=InParallel
      ! We must turn off the parallel broadcast at this point
      ! so that we only gather from HDF to the ROOT node
      InParallel=.FALSE.
    ENDIF
#endif
    CALL Get_BCSR(B,Name,PFix_O,CheckPoint_O)
#ifdef PARALLEL
    IF(PRESENT(Checkpoint_O))THEN
      InParallel=InParTemp
    ENDIF
#endif
    CALL SetEq(A,B)
    CALL Delete(B)
    A%Node=MyId
  END SUBROUTINE Get_DBCSR
  !-------------------------------------------------------------------------------
  !     Put a DBCSR matrix

  SUBROUTINE Put_DBCSR(A,Name,PFix_O,CheckPoint_O)
    TYPE(DBCSR), INTENT(INOUT)           :: A
    LOGICAL,         OPTIONAL,INTENT(IN) :: CheckPoint_O
    TYPE(BCSR)                           :: B
    CHARACTER(LEN=*),         INTENT(IN) :: Name
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: PFix_O
    !-------------------------------------------------------------------------------
    CALL SetEq(B,A)
    CALL Put_BCSR(B,Name,PFix_O,CheckPoint_O)
    CALL Delete(B)
    !-------------------------------------------------------------------------------
  END SUBROUTINE Put_DBCSR
#endif
  !-------------------------------------------------------------------------------
  !     Put thresholds

  SUBROUTINE Put_TOLS(NGLCT,Tag_O)
    TYPE(TOLS),               INTENT(IN) :: NGLCT
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    CHARACTER(LEN=DEFAULT_CHR_LEN)       :: Tag
    CALL Put(NGLCT%Cube,'cubeneglect',Tag_O=Tag_O)
    CALL Put(NGLCT%Trix,'trixneglect',Tag_O=Tag_O)
    CALL Put(NGLCT%Dist,'distneglect',Tag_O=Tag_O)
    CALL Put(NGLCT%TwoE,'twoeneglect',Tag_O=Tag_O)
    CALL Put(NGLCT%ETol,'enregyneglect',Tag_O=Tag_O)
    CALL Put(NGLCT%DTol,'densityneglect',Tag_O=Tag_O)
  END SUBROUTINE Put_TOLS

  !-------------------------------------------------------------------------------
  !     Get thresholds

  SUBROUTINE Get_TOLS(NGLCT,Tag_O)
    IMPLICIT NONE
    TYPE(TOLS)                           :: NGLCT
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    CALL Get(NGLCT%Cube,'cubeneglect',Tag_O=Tag_O)
    CALL Get(NGLCT%Trix,'trixneglect',Tag_O=Tag_O)
    CALL Get(NGLCT%Dist,'distneglect',Tag_O=Tag_O)
    CALL Get(NGLCT%TwoE,'twoeneglect',Tag_O=Tag_O)
    CALL Get(NGLCT%ETol,'enregyneglect',Tag_O=Tag_O)
    CALL Get(NGLCT%DTol,'densityneglect',Tag_O=Tag_O)
  END SUBROUTINE Get_TOLS

  !-------------------------------------------------------------------------------
  !     Get arguments from the command line

  SUBROUTINE Get_ARGMT(A)
#ifdef NAG
    USE F90_UNIX
    IMPLICIT NONE
#else
    IMPLICIT NONE
#endif

#if ! defined (GFORTRAN)
    INTEGER,EXTERNAL               :: IARGC
#endif
    TYPE(ARGMT)                    :: A
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: Tmp1,Tmp2
    INTEGER                        :: I,NArg,NChar,NInts
    INTEGER,PARAMETER              :: MaxArg=10
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      NArg=IARGC()
      NChar=0
      NInts=0
      DO I=1,NArg
        CALL GetArg(I,Tmp1)
        CALL LowCase(Tmp1)
        IF(SCAN(Tmp1,Lower)==0)THEN
          NInts=NInts+1
        ELSE
          NChar=NChar+1
        ENDIF
      ENDDO
#ifdef PARALLEL
    ENDIF
    IF(InParallel)THEN
      CALL Bcast(NInts)
      CALL Bcast(NChar)
    ENDIF
#endif
    CALL New(A,(/NChar,NInts/))
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      NChar=0
      NInts=0
      DO I=1,NArg
        CALL GetArg(I,Tmp1)
        Tmp2=Tmp1
        CALL LowCase(Tmp2)
        IF(SCAN(Tmp2,Lower)==0)THEN
          NInts=NInts+1
          A%I%I(NInts)=CharToInt(TRIM(Tmp1))
        ELSE
          NChar=NChar+1
          A%C%C(NChar)=TRIM(Tmp1)
        ENDIF
      ENDDO
#ifdef PARALLEL
    ENDIF
    IF(InParallel)THEN
      CALL Bcast(A%I)
      CALL Bcast_CHR_VECT(A%C)
    ENDIF
#endif
  END SUBROUTINE Get_ARGMT



  SUBROUTINE Get_HGRho(A,Name,Args,SCFCycle,BCast_O)
    TYPE(HGRho)                      :: A
    TYPE(ARGMT)                      :: Args
    INTEGER                          :: SCFCycle,IOS,I,NExpt,NDist,NCoef,NSDen
    CHARACTER(LEN=*)                 :: Name
    CHARACTER(LEN=DEFAULT_CHR_LEN)   :: FileName
    LOGICAL                          :: Exists
    LOGICAL,OPTIONAL,INTENT(IN)      :: BCast_O
    LOGICAL                          :: BcastQ

#ifdef PARALLEL
    IF(PRESENT(BCast_O)) THEN
      BcastQ = BCast_O
    ELSE
      BcastQ = .FALSE.
    ENDIF
#endif

#ifdef PARALLEL
    IF(MyID == ROOT) THEN
#endif
      FileName=TrixFile(Name,Args,SCFCycle)
      INQUIRE(FILE=FileName,EXIST=Exists)
      IF(Exists) THEN
        OPEN(UNIT=Seq,FILE=FileName,STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      ELSE
        CALL Halt(' Get_HGRho could not find '//TRIM(FileName))
      ENDIF

      CALL MondoLog(DEBUG_MAXIMUM, "Get_HGRho", "getting rho from "//TRIM(FileName))

      ! Allocate Memory
      READ(UNIT=Seq,Err=100,IOSTAT=IOS) NSDen,NExpt,NDist,NCoef
#ifdef PARALLEL
    ENDIF
    IF(BcastQ) THEN
      CALL Bcast(NSDen)
      CALL Bcast(NExpt)
      CALL Bcast(NDist)
      CALL Bcast(NCoef)
      CALL Bcast(NSDen)
    ENDIF
#endif
    CALL New_HGRho(A,(/NExpt,NDist,NCoef,NSDen/))

    !CALL MondoLog(DEBUG_MAXIMUM, "Get_HGRho", &
    !    "NSDen = "//TRIM(IntToChar(A%NSDen))//", " &
    !  //"NExpt = "//TRIM(IntToChar(A%NExpt))//", " &
    !  //"NDist = "//TRIM(IntToChar(A%NDist))//", " &
    !  //"NCoef = "//TRIM(IntToChar(A%NCoef))//", " &
    !  //"MyID = "//TRIM(IntToChar(MyID)))

    ! Read In the Density
#ifdef PARALLEL
    IF(MyID == ROOT) THEN
#endif
      READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%NQ%I(I)    ,I=1,A%NExpt)
      READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%OffQ%I(I)  ,I=1,A%NExpt)
      READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%OffR%I(I)  ,I=1,A%NExpt)
      READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Lndx%I(I)  ,I=1,A%NExpt)
      READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Expt%D(I)  ,I=1,A%NExpt)
      READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qx%D(I)    ,I=1,A%NDist)
      READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qy%D(I)    ,I=1,A%NDist)
      READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qz%D(I)    ,I=1,A%NDist)
      READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Co%D(I)    ,I=1,A%NCoef*A%NSDen)
#ifdef PARALLEL
    ENDIF
    IF(BcastQ) THEN
      Call Bcast(A%NQ,N_O=A%NExpt)
      Call Bcast(A%OffQ,N_O=A%NExpt)
      Call Bcast(A%OffR,N_O=A%NExpt)
      Call Bcast(A%Lndx,N_O=A%NExpt)
      Call Bcast(A%Expt,N_O=A%NExpt)
      Call Bcast(A%Qx,N_O=A%NDist)
      Call Bcast(A%Qy,N_O=A%NDist)
      Call Bcast(A%Qz,N_O=A%NDist)
      Call Bcast(A%Co,N_O=A%NCoef*A%NSDen)
    ENDIF
    IF(MyID == ROOT) THEN
#endif

      Close(UNIT=Seq,STATUS='KEEP')
#ifdef PARALLEL
    ENDIF
#endif
    RETURN
100 CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Get_HGRho.')
  END SUBROUTINE Get_HGRho
  !===============================================================================
  ! Write  the density to disk
  !===============================================================================
  SUBROUTINE Put_HGRho(A,Name,Args,SCFCycle)
    TYPE(HGRho)                      :: A
    TYPE(ARGMT)                      :: Args
    INTEGER                          :: I,SCFCycle,IOS
    REAL(DOUBLE)                     :: Dummy
    CHARACTER(LEN=*)                 :: Name
    CHARACTER(LEN=DEFAULT_CHR_LEN)   :: FileName
    LOGICAL                          :: Exists

    FileName=TrixFile(Name,Args,SCFCycle)
    INQUIRE(FILE=FileName,EXIST=Exists)
    IF(Exists) THEN
      OPEN(UNIT=Seq,FILE=FileName,STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    ELSE
      OPEN(UNIT=Seq,FILE=FileName,STATUS='NEW',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    ENDIF

    CALL MondoLog(DEBUG_MAXIMUM, "Put_HGRho", "putting rho to "//TRIM(FileName))

    ! Write density to disk
    !CALL MondoLog(DEBUG_MAXIMUM, "Put_HGRho", &
    !    "NSDen = "//TRIM(IntToChar(A%NSDen))//", " &
    !  //"NExpt = "//TRIM(IntToChar(A%NExpt))//", " &
    !  //"NDist = "//TRIM(IntToChar(A%NDist))//", " &
    !  //"NCoef = "//TRIM(IntToChar(A%NCoef))//", " &
    !  //"MyID = "//TRIM(IntToChar(MyID)))

    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS) A%NSDen,A%NExpt,A%NDist,A%NCoef
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%NQ%I(I)    ,I=1,A%NExpt)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%OffQ%I(I)  ,I=1,A%NExpt)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%OffR%I(I)  ,I=1,A%NExpt)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Lndx%I(I)  ,I=1,A%NExpt)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Expt%D(I)  ,I=1,A%NExpt)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qx%D(I)    ,I=1,A%NDist)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qy%D(I)    ,I=1,A%NDist)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qz%D(I)    ,I=1,A%NDist)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Co%D(I)    ,I=1,A%NCoef*A%NSDen)

    CLOSE(UNIT=Seq,STATUS='KEEP')
    RETURN
100 CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Put_HGRho.')
  END SUBROUTINE Put_HGRho

  !===============================================================================
  ! Get Cartesian multipoles
  !===============================================================================
  SUBROUTINE Get_CMPoles(A,Tag_O)
    TYPE(CMPoles)                    :: A
    CHARACTER(LEN=*),OPTIONAL        :: Tag_O

    IF(.NOT.AllocQ(A%Alloc))CALL New_CMPoles(A)
    ! Get will broadcast DPole and QPole automatically.
    CALL Get(A%DPole,'dipole',Tag_O=Tag_O)
    CALL Get(A%QPole,'quadrupole',Tag_O=Tag_O)
  END SUBROUTINE Get_CMPoles
  !===============================================================================
  ! Put Cartesian multipoles
  !===============================================================================
  SUBROUTINE Put_CMPoles(A,Tag_O)
    TYPE(CMPoles)                    :: A
    CHARACTER(LEN=*),OPTIONAL        :: Tag_O
    CALL Put(A%DPole,'dipole',Tag_O=Tag_O)
    CALL Put(A%QPole,'quadrupole',Tag_O=Tag_O)
  END SUBROUTINE Put_CMPoles

  !-------------------------------------------------------------------------------
  !     Open an ASCII file
  !
  SUBROUTINE OpenASCII(FileName,Unit,NewFile_O,OldFileQ_O,Rewind_O)
    CHARACTER(LEN=*), INTENT(IN) :: FileName
    INTEGER,          INTENT(IN) :: Unit
    LOGICAL, OPTIONAL,INTENT(IN) :: NewFile_O,OldFileQ_O,Rewind_O
    INTEGER                      :: IOS
    LOGICAL                      :: Opened, Exists

    !-------------------------------------------------------------------------------
    !        Does the file exist?

    INQUIRE(FILE=FileName,OPENED=Opened, &
         EXIST=Exists,ERR=11,IOSTAT=IOS)
    IF(PRESENT(OldFileQ_O))THEN
      IF(OldFileQ_O.AND.(.NOT.Exists)) &
           CALL HALT(' File '//TRIM(FileName)//' does not exist! ')
    ENDIF
    !-------------------------------------------------------------------------------
    IF(PRESENT(NewFile_O))THEN
      IF(NewFile_O.AND.Exists)THEN
        ! Open replace if already exists
        OPEN(UNIT=Unit,FILE=FileName, &
             ACCESS='SEQUENTIAL', FORM='FORMATTED', &
             ERR=11,IOSTAT=IOS,STATUS='REPLACE')
      ELSEIF(NewFile_O)THEN
        ! Open a new file
        OPEN(UNIT=Unit,FILE=FileName, &
             ACCESS='SEQUENTIAL',FORM='FORMATTED', &
             ERR=11,IOSTAT=IOS,STATUS='NEW')
      ELSEIF(Exists)THEN
        ! Open an old file (NewFile_O=.FALSE.)
        OPEN(UNIT=Unit,FILE=FileName, &
             ACCESS='SEQUENTIAL', FORM='FORMATTED', &
             POSITION='APPEND',ERR=11,IOSTAT=IOS,STATUS='OLD')
      ELSE
        CALL Halt(' Bad logic in OpenASCI' )
      ENDIF
    ELSEIF(PRESENT(Rewind_O))THEN
      ! Open existing file and position at the top
      IF(Rewind_O.AND.Exists)THEN
        OPEN(UNIT=Unit,FILE=FileName, &
             ACCESS='SEQUENTIAL', FORM='FORMATTED', &
             POSITION='REWIND',ERR=11,IOSTAT=IOS,STATUS='OLD')
      ELSEIF(Rewind_O)THEN
        ! Just open a new file
        OPEN(UNIT=Unit,FILE=FileName, &
             ACCESS='SEQUENTIAL',FORM='FORMATTED', &
             ERR=11,IOSTAT=IOS,STATUS='NEW')
      ELSE
        CALL Halt(' Bad logic in OpenASCI' )
      ENDIF
      !-------------------------------------------------------------------------------
      !        Open existing file and position at the bottom (default)

    ELSEIF(Exists.AND.(.NOT.Opened))THEN
      OPEN(UNIT=Unit,FILE=FileName, &
           ACCESS='SEQUENTIAL', FORM='FORMATTED', &
           POSITION='APPEND',ERR=11,IOSTAT=IOS,STATUS='OLD')
      !-------------------------------------------------------------------------------
      !        Create a new file and open it

    ELSEIF(Exists.AND.Opened)THEN
      CALL Warn(' File '//TRIM(FileName)//' already open')
    ELSE
      OPEN(UNIT=Unit,FILE=FileName, &
           ACCESS='SEQUENTIAL',FORM='FORMATTED', &
           ERR=11,IOSTAT=IOS,STATUS='NEW')
    ENDIF
    RETURN
11  CALL Halt(' OpenASCII ERROR: IOS='//TRIM(IntToChar(IOS))// &
         ' on file '//TRIM(FileName)//'.')
  END SUBROUTINE OpenASCII

  !-------------------------------------------------------------------------------
  SUBROUTINE Put_CHR_VECT(A,VarName,Tag_O)
    INTEGER                              :: I,N,II,NN
    TYPE(CHR_VECT) :: A
    CHARACTER(LEN=*),         INTENT(IN) :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    TYPE(META_DATA)                      :: Meta

#ifdef OLD_CHR_VECT
    INTEGER,DIMENSION(DCL)   :: B  !=ICHAR(' ')
#else
    INTEGER,ALLOCATABLE :: B(:)
    INTEGER :: RunInd,BufSize
#endif

#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif

#ifdef OLD_CHR_VECT
      NN=SIZE(A%C)
      DO II = 1, NN
        N=LEN(A%C(II))
        IF(N>DCL) CALL Halt('Static strings overrun in Put_CHR_VECT')
        DO I=1,N
          B(I)=ICHAR(A%C(II)(I:I))
        ENDDO
        Meta=SetMeta(NameTag(VarName,TRIM(IntToChar(II))// &
             TRIM(Tag_O)),NATIVE_INT32,N,.FALSE.)
        CALL OpenData(Meta,.TRUE.)
        CALL WriteIntegerVector(Meta,B)
        CALL CloseData(Meta)
      ENDDO
#else
      NN = SIZE(A%C)
      BufSize = 1
      DO II = 1, NN
        N = LEN(TRIM(A%C(II)))
        IF(N>DCL) CALL Halt('Static strings overrun in Put_CHR_VECT')
        BufSize = BufSize + N + 1
      ENDDO

      ALLOCATE(B(BufSize))
      RunInd = 1
      B(RunInd) = NN
      RunInd = RunInd + 1
      DO II = 1, NN
        N = LEN(TRIM(A%C(II)))
        DO I = 1, N
          B(RunInd) = ICHAR(A%C(II)(I:I))
          RunInd = RunInd + 1
        ENDDO
        B(RunInd) = -1000
        RunInd = RunInd + 1
      ENDDO
      IF(BufSize /= RunInd - 1) STOP 'ERR: Index problem in Put_CHR_VECT !'
      CALL Put(BufSize,NameTag(VarName,TRIM(IntToChar(0))//TRIM(Tag_O)))
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_CHR_VECT", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_CHR_VECT", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,TRIM(IntToChar(1))// &
           TRIM(Tag_O)),NATIVE_INT32,BufSize,.FALSE.)
      CALL OpenData(Meta,.TRUE.)
      CALL WriteIntegerVector(Meta,B)
      CALL CloseData(Meta)
      DEALLOCATE(B)
#endif
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Put_CHR_VECT

  !-------------------------------------------------------------------------------
  SUBROUTINE Get_CHR_VECT(A,VarName,Tag_O)
    INTEGER                                 :: I,N,II,NN
    TYPE(CHR_VECT),           INTENT(INOUT) :: A
    CHARACTER(LEN=*),         INTENT(IN)    :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
#ifdef OLD_CHR_VECT
    INTEGER,DIMENSION(DEFAULT_CHR_LEN)      :: B !=ICHAR(' ')
#else
    INTEGER,ALLOCATABLE :: B(:)
    INTEGER :: RunInd,StrInd,StrLen,BufSize
    CHARACTER(LEN=DCL) :: TEMP
#endif
    TYPE(META_DATA)                         :: Meta

#ifdef OLD_CHR_VECT
#else
    CALL Get(BufSize,NameTag(VarName,TRIM(IntToChar(0))//TRIM(Tag_O)))
#endif

#ifdef PARALLEL

    IF(MyId==ROOT)THEN
#endif

#ifdef OLD_CHR_VECT
      NN=SIZE(A%C)
      DO II = 1, NN
        N=LEN(A%C(II))
        IF(N>DEFAULT_CHR_LEN) CALL Halt('Static strings overrun in Get_CHR_VECT')
        Meta=SetMeta(NameTag(VarName,TRIM(IntToChar(II))//TRIM(Tag_O)),NATIVE_INT32,N,.FALSE.)

        CALL OpenData(Meta)
        CALL ReadIntegerVector(Meta,B(1))
        CALL CloseData(Meta)
        DO I=1,N
          A%C(II)(I:I)=CHAR(B(I))
        ENDDO
      ENDDO
#else
      NN = SIZE(A%C)
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_CHR_VECT", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_CHR_VECT", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,TRIM(IntToChar(1))//TRIM(Tag_O)),NATIVE_INT32,BufSize,.FALSE.)
      CALL OpenData(Meta)
      ALLOCATE(B(BufSize))
      CALL ReadIntegerVector(Meta,B(1))
      IF(NN /= B(1)) STOP 'ERR: SIZE problem in Get_CHR_VECT'
      RunInd = 2
      DO II = 1, NN
        TEMP = ' '
        StrInd = 1
        DO
          IF(B(RunInd) == -1000) THEN
            RunInd = RunInd + 1
            EXIT
          ELSE
            TEMP(StrInd:StrInd) = CHAR(B(RunInd))
            StrInd = StrInd + 1
            RunInd = RunInd+1
          ENDIF
        ENDDO
        A%C(II) = TEMP
      ENDDO
      IF(BufSize /= RunInd -1 ) STOP 'ERR: Index problem in Get_CHR_VECT'
      DEALLOCATE(B)
#endif

#ifdef PARALLEL
    ENDIF
    ! not supported yet. IF(InParallel)CALL Bcast(A)
    !! STOP 'ERROR : Bcast in Get_Chr_vect (InOut.F90) not supported!'
#endif
  END SUBROUTINE Get_CHR_VECT
  !-------------------------------------------------------------------------------
  SUBROUTINE Put_CHR10_VECT(A,VarName,Tag_O)
    INTEGER                               :: I,N,II,NN
    TYPE(CHR10_VECT)                      :: A
    CHARACTER(LEN=*),         INTENT(IN)  :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)  :: Tag_O
    TYPE(META_DATA)                       :: Meta
    INTEGER,ALLOCATABLE                   :: B(:)
    INTEGER                               :: RunInd,BufSize
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      NN = SIZE(A%C)
      BufSize = 1
      DO II = 1, NN
        N = LEN(TRIM(A%C(II)))
        IF(N>DCL) CALL Halt('Static strings overrun in Put_CHR10_VECT')
        BufSize = BufSize + N + 1
      ENDDO

      ALLOCATE(B(BufSize))
      RunInd = 1
      B(RunInd) = NN
      RunInd = RunInd + 1
      DO II = 1, NN
        N = LEN(TRIM(A%C(II)))
        DO I = 1, N
          B(RunInd) = ICHAR(A%C(II)(I:I))
          RunInd = RunInd + 1
        ENDDO
        B(RunInd) = -1000
        RunInd = RunInd + 1
      ENDDO
      IF(BufSize /= RunInd - 1) STOP 'ERR: Index problem in Put_CHR10_VECT !'
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_CHR10_VECT", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_CHR10_VECT", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      IF(PRESENT(Tag_O))THEN
        CALL Put(BufSize,NameTag(VarName//TRIM(IntToChar(0)), Tag_O))
        Meta=SetMeta(NameTag(VarName//TRIM(IntToChar(1)), Tag_O),NATIVE_INT32,BufSize,.FALSE.)
      ELSE
        CALL Put(BufSize,NameTag(VarName//TRIM(IntToChar(0))))
        Meta=SetMeta(NameTag(VarName//TRIM(IntToChar(1))),NATIVE_INT32,BufSize,.FALSE.)
      ENDIF
      CALL OpenData(Meta,.TRUE.)
      CALL WriteIntegerVector(Meta,B)
      CALL CloseData(Meta)
      DEALLOCATE(B)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Put_CHR10_VECT

  !-------------------------------------------------------------------------------
  SUBROUTINE Get_CHR10_VECT(A,VarName,Tag_O)
    INTEGER                                 :: I,N,II,NN
    TYPE(CHR10_VECT),           INTENT(INOUT) :: A
    CHARACTER(LEN=*),         INTENT(IN)    :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
    INTEGER,ALLOCATABLE :: B(:)
    INTEGER :: RunInd,StrInd,StrLen,BufSize
    CHARACTER(LEN=DCL) :: TEMP
    TYPE(META_DATA)                         :: Meta
    IF(PRESENT(Tag_O))THEN
      CALL Get(BufSize,NameTag(VarName//TRIM(IntToChar(0)), Tag_O))
    ELSE
      CALL Get(BufSize,NameTag(VarName//TRIM(IntToChar(0))))
    ENDIF
#ifdef PARALLEL

    IF(MyId==ROOT)THEN
#endif
      NN = SIZE(A%C)
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_CHR10_VECT", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_CHR10_VECT", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      IF(PRESENT(Tag_O))THEN
        Meta=SetMeta(NameTag(VarName//TRIM(IntToChar(1)), TRIM(Tag_O)),NATIVE_INT32,BufSize,.FALSE.)
      ELSE
        Meta=SetMeta(NameTag(VarName//TRIM(IntToChar(1))),NATIVE_INT32,BufSize,.FALSE.)
      ENDIF
      CALL OpenData(Meta)
      ALLOCATE(B(BufSize))
      CALL ReadIntegerVector(Meta,B(1))
      IF(NN /= B(1)) STOP 'ERR: SIZE problem in Get_CHR10_VECT'
      RunInd = 2
      DO II = 1, NN
        TEMP = ' '
        StrInd = 1
        DO
          IF(B(RunInd) == -1000) THEN
            RunInd = RunInd + 1
            EXIT
          ELSE
            TEMP(StrInd:StrInd) = CHAR(B(RunInd))
            StrInd = StrInd + 1
            RunInd = RunInd+1
          ENDIF
        ENDDO
        A%C(II) = TEMP
      ENDDO
      IF(BufSize /= RunInd -1 ) STOP 'ERR: Index problem in Get_CHR10_VECT'
      DEALLOCATE(B)
#ifdef PARALLEL
    ENDIF
    ! not supported yet. IF(InParallel)CALL Bcast(A)
    !! STOP 'ERROR : Bcast in Get_Chr10_vect (InOut.F90) not supported!'
#endif
  END SUBROUTINE Get_CHR10_VECT
  !-------------------------------------------------------------------------------
  SUBROUTINE Put_LOG_VECT(A,VarName,Tag_O)
    INTEGER                              :: I,N,II,NN
    TYPE(LOG_VECT)                       :: A
    CHARACTER(LEN=*),         INTENT(IN) :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    TYPE(INT_VECT)                       :: ILog
    TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      !
      NN=SIZE(A%L)
      CALL New(ILog,NN)
      DO I = 1, NN
        ILog%I(I)=0
        IF(A%L(I)) ILog%I(I)=1
      ENDDO
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_LOG_VECT", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Put_LOG_VECT", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,NN)
      CALL OpenData(Meta,.TRUE.)
      CALL WriteIntegerVector(Meta,ILog%I)
      CALL CloseData(Meta)
      CALL Delete(ILog)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Put_LOG_VECT
  !-------------------------------------------------------------------------------
  SUBROUTINE Get_LOG_VECT(A,VarName,Tag_O)
    INTEGER                              :: I,NN
    TYPE(LOG_VECT)                       :: A
    CHARACTER(LEN=*),         INTENT(IN) :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    TYPE(INT_VECT)                       :: ILog
    TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      NN=SIZE(A%L)
      CALL New(ILog,NN)
      IF(PRESENT(Tag_O)) THEN
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_LOG_VECT", "VarName = "//TRIM(VarName)//", Tag_O = "//TRIM(Tag_O))
      ELSE
        !CALL MondoLog(DEBUG_MAXIMUM, "Get_LOG_VECT", "VarName = "//TRIM(VarName)//", Tag_O not set")
      ENDIF
      Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,SIZE(ILog%I,1),.FALSE.)
      CALL OpenData(Meta)
      CALL ReadIntegerVector(Meta,ILog%I)
      CALL CloseData(Meta)
      DO I = 1, NN
        IF(ILog%I(I)==1) Then
          A%L(I) = .TRUE.
        ELSE
          A%L(I) = .FALSE.
        ENDIF
      ENDDO
      CALL Delete(ILog)
#ifdef PARALLEL
    ENDIF
    IF(InParallel) CALL Bcast(A)
#endif
  END SUBROUTINE Get_LOG_VECT
  !-------------------------------------------------------------------------------
  SUBROUTINE Put_CellSet(CS,Name_O,Tag_O,Unlimit_O)
    TYPE(CellSet)                  :: CS
    CHARACTER(Len=*),Optional            :: Name_O
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    LOGICAL,         OPTIONAL,INTENT(IN) :: Unlimit_O
    !
    IF(PRESENT(Name_O))THEN
      CALL Put(CS%Radius   ,TRIM(Name_O)//'_cell_radius',Tag_O=Tag_O)
      CALL Put(CS%NCells   ,TRIM(Name_O)//'_cell_number',Tag_O=Tag_O)
      CALL Put(CS%CellCarts,TRIM(Name_O)//'_cell_vectors',Tag_O=Tag_O,Unlimit_O=.TRUE.)
    ELSE
      CALL Put(CS%Radius   ,'cell_radius',Tag_O=Tag_O)
      CALL Put(CS%NCells   ,'cell_number',Tag_O=Tag_O)
      CALL Put(CS%CellCarts,TRIM(Name_O)//'_cell_vectors',Tag_O=Tag_O,Unlimit_O=.TRUE.)
    ENDIF
  END SUBROUTINE Put_CellSet

  SUBROUTINE Get_CellSet(CS,Name_O,Tag_O)
    TYPE(CellSet)                  :: CS
    CHARACTER(Len=*),OPTIONAL      :: Name_O
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    INTEGER                        :: NC

    IF(PRESENT(Name_O))THEN
      CALL Get(CS%Radius   ,TRIM(Name_O)//'_cell_radius',Tag_O=Tag_O)
      CALL Get(CS%NCells   ,TRIM(Name_O)//'_cell_number',Tag_O=Tag_O)
      ! Allocate a new CellSet with the information just read.
      CALL New_CellSet(CS,CS%NCells)
      CALL Get(CS%CellCarts,TRIM(Name_O)//'_cell_vectors',Tag_O=Tag_O)
    ELSE
      CALL Get(CS%Radius   ,'cell_radius',Tag_O=Tag_O)
      CALL Get(CS%NCells   ,'cell_number',Tag_O=Tag_O)
      ! Allocate a new CellSet with the information just read.
      CALL New_CellSet(CS,CS%NCells)
      CALL Get(CS%CellCarts,'cell_vectors',Tag_O=Tag_O)
    ENDIF
  END SUBROUTINE Get_CellSet
  !
  SUBROUTINE Get_Sp1x1(A,Tag,Symb_O)
    TYPE(Sp1x1)               :: A
    CHARACTER(LEN=*)          :: Tag
    LOGICAL,OPTIONAL          :: Symb_O
    !
    CALL Get(A%NRow,'NRow',Tag_O=Tag)
    CALL Get(A%NZ,'NZ',Tag_O=Tag)
    IF(PRESENT(Symb_O)) THEN
      CALL New(A,A%NRow,A%NZ,Symb_O=Symb_O)
    ELSE
      CALL New(A,A%NRow,A%NZ)
    ENDIF
    CALL Get(A%IA,'IA',Tag_O=Tag)
    CALL Get(A%JA,'JA',Tag_O=Tag)
    IF(PRESENT(Symb_O)) THEN
      IF(Symb_O) THEN
        RETURN
      ELSE
        CALL Get(A%AN,'AN',Tag_O=Tag)
      ENDIF
    ELSE
      CALL Get(A%AN,'AN',Tag_O=Tag)
    ENDIF
  END SUBROUTINE Get_Sp1x1
  !
  !----------------------------------------------------
  !
  SUBROUTINE Put_Sp1x1(A,Tag,Symb_O)
    TYPE(Sp1x1)      :: A
    CHARACTER(LEN=*) :: Tag
    LOGICAL,OPTIONAL :: Symb_O
    INTEGER          :: NRow,NZ
    !
    CALL Put(A%NRow,'NRow',Tag_O=Tag)
    CALL Put(A%NZ,'NZ',Tag_O=Tag)
    CALL Put(A%IA,'IA',Tag_O=Tag)
    CALL Put(A%JA,'JA',Tag_O=Tag)
    IF(PRESENT(Symb_O)) THEN
      IF(Symb_O) THEN
        RETURN
      ELSE
        CALL Put(A%AN,'AN',Tag_O=Tag)
      ENDIF
    ELSE
      CALL Put(A%AN,'AN',Tag_O=Tag)
    ENDIF
  END SUBROUTINE Put_Sp1x1
  !
  !----------------------------------------------------------------
  !
  SUBROUTINE Put_BondD(A,Name,Tag_O)
    TYPE(BONDDATA)   :: A
    CHARACTER(LEN=*) :: Name
    CHARACTER(LEN=*),OPTIONAL :: Tag_O
    !
    IF(.NOT.AllocQ(A%Alloc)) THEN
      A%N=0
    ENDIF
    CALL Put(A%N,TRIM(Name)//'N',Tag_O=Tag_O)
    IF(A%N/=0) THEN
      CALL Put(A%IJ,TRIM(Name)//'IJ',Tag_O=Tag_O)
      CALL Put(A%Length,TRIM(Name)//'Length',Tag_O=Tag_O)
      CALL Put(A%Type,TRIM(Name)//'Type',Tag_O=Tag_O)
      CALL Put(A%HBExtraSN,TRIM(Name)//'HBExtraSN',Tag_O=Tag_O)
      CALL Put(A%HBExtraNC,TRIM(Name)//'HBExtraNC',Tag_O=Tag_O)
      CALL Put(A%LonelyAtom,TRIM(Name)//'LonelyAtom',Tag_O=Tag_O)
    ENDIF
  END SUBROUTINE Put_BondD
  !
  !---------------------------------------------------------------
  !
  SUBROUTINE Get_BondD(A,Name,Tag_O)
    TYPE(BONDDATA)   :: A
    CHARACTER(LEN=*) :: Name
    CHARACTER(LEN=*),OPTIONAL :: Tag_O
    !
    CALL Get(A%N,TRIM(Name)//'N',Tag_O=Tag_O)
    IF(A%N/=0) THEN
      CALL New(A,A%N)
      CALL Get(A%IJ,TRIM(Name)//'IJ',Tag_O=Tag_O)
      CALL Get(A%Length,TRIM(Name)//'Length',Tag_O=Tag_O)
      CALL Get(A%Type,TRIM(Name)//'Type',Tag_O=Tag_O)
      CALL Get(A%HBExtraSN,TRIM(Name)//'HBExtraSN',Tag_O=Tag_O)
      CALL Get(A%HBExtraNC,TRIM(Name)//'HBExtraNC',Tag_O=Tag_O)
      CALL Get(A%LonelyAtom,TRIM(Name)//'LonelyAtom',Tag_O=Tag_O)
    ENDIF
  END SUBROUTINE Get_BondD
  !
  !---------------------------------------------------------------
  !
  SUBROUTINE Put_AtmB(A,Name,Tag_O)
    TYPE(ATOMBONDS)  :: A
    CHARACTER(LEN=*) :: Name
    CHARACTER(LEN=*),OPTIONAL :: Tag_O
    INTEGER          :: N1,N2
    !
    IF(.NOT.AllocQ(A%Alloc)) THEN
      A%N1=0
      A%N2=0
    ENDIF
    CALL Put(A%N1,TRIM(Name)//'N1',Tag_O=Tag_O)
    CALL Put(A%N2,TRIM(Name)//'N2',Tag_O=Tag_O)
    IF(A%N1==0) RETURN
    IF(A%N1/=0) THEN
      CALL Put(A%Count,TRIM(Name)//'Count',Tag_O=Tag_O)
      IF(A%N2/=0) THEN
        CALL Put(A%Bonds,TRIM(Name)//'Bonds',Tag_O=Tag_O)
        CALL Put(A%Atoms,TRIM(Name)//'Atoms',Tag_O=Tag_O)
      ENDIF
    ENDIF
  END SUBROUTINE Put_AtmB
  !
  !---------------------------------------------------------------
  !
  SUBROUTINE Get_AtmB(A,Name,Tag_O)
    TYPE(ATOMBONDS)  :: A
    CHARACTER(LEN=*) :: Name
    CHARACTER(LEN=*),OPTIONAL :: Tag_O
    INTEGER          :: N1,N2
    !
    CALL Get(A%N1,TRIM(Name)//'N1',Tag_O=Tag_O)
    CALL Get(A%N2,TRIM(Name)//'N2',Tag_O=Tag_O)
    CALL New(A,A%N1,A%N2)
    IF(A%N1/=0) THEN
      CALL Get(A%Count,TRIM(Name)//'Count',Tag_O=Tag_O)
      IF(A%N2/=0) THEN
        CALL Get(A%Bonds,TRIM(Name)//'Bonds',Tag_O=Tag_O)
        CALL Get(A%Atoms,TRIM(Name)//'Atoms',Tag_O=Tag_O)
      ENDIF
    ENDIF
  END SUBROUTINE Get_AtmB
  !
  !---------------------------------------------------------------
  !
  SUBROUTINE Get_PBCFit(A,Name,Tag_O)
    TYPE(PBCFits) :: A
    CHARACTER(LEN=*) :: Name
    CHARACTER(LEN=*),OPTIONAL :: Tag_O
    !
    CALL Get(A%MaxMem,TRIM(Name)//'MaxMem',Tag_O=Tag_O)
    CALL Get(A%ActMem,TRIM(Name)//'ActMem',Tag_O=Tag_O)
    IF(A%MaxMem>LattMaxMem) CALL Halt('A%MaxMem>LattMaxMem in Get_PBCFit')
    IF(.NOT.AllocQ(A%Alloc)) THEN
      CALL New(A,LattMaxMem)
    ELSE
      IF(A%MaxMem/=LattMaxMem) CALL Halt('A%MaxMem/=LattMaxMem in Get_PBCFit')
    ENDIF
    CALL Get(A%AWeights,TRIM(Name)//'AWeights',Tag_O=Tag_O)
    CALL Get(A%PBCValues,TRIM(Name)//'PBCValues',Tag_O=Tag_O)
    CALL Get(A%PBCGrads,TRIM(Name)//'PBCGrads',Tag_O=Tag_O)
  END SUBROUTINE Get_PBCFit
  !
  !---------------------------------------------------------------
  !
  SUBROUTINE Put_PBCFit(A,Name,Tag_O)
    TYPE(PBCFits)    :: A
    CHARACTER(LEN=*) :: Name
    CHARACTER(LEN=*),OPTIONAL :: Tag_O
    !
    CALL Put(A%MaxMem,TRIM(Name)//'MaxMem',Tag_O=Tag_O)
    CALL Put(A%ActMem,TRIM(Name)//'ActMem',Tag_O=Tag_O)
    CALL Put(A%AWeights,TRIM(Name)//'AWeights',Tag_O=Tag_O)
    CALL Put(A%PBCValues,TRIM(Name)//'PBCValues',Tag_O=Tag_O)
    CALL Put(A%PBCGrads,TRIM(Name)//'PBCGrads',Tag_O=Tag_O)
  END SUBROUTINE Put_PBCFit

  SUBROUTINE HDF5DeleteObject(ID, groupName, objectName)
    INTEGER, INTENT(IN)           :: ID
    CHARACTER(LEN=*), INTENT(IN)  :: groupName, objectName

    CALL HDF5Delete(ID, &
         LEN(TRIM(groupName)), Char2Ints(LEN(TRIM(groupName)), groupName), &
         LEN(TRIM(objectName)), Char2Ints(LEN(TRIM(objectName)), objectName))

  END SUBROUTINE HDF5DeleteObject

END MODULE InOut
