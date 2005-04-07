!    GENERIC IO ROUTINES FOR MONDOSCF TYPES
!    Author: Matt Challacombe and CK Gan
!-------------------------------------------------------------------------------
MODULE InOut
  USE DerivedTypes
  USE GlobalScalars   
  USE GlobalCharacters
  USE ProcessControl    
  USE MemMan
  USE Indexing
  USE Parse
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
       IF(FileID==FAIL) &
            CALL Halt(' Failed to open the HDF file <'//TRIM(FileName)//'>.') 
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
       STATUS=HDF5CloseFile(FileID)
       IF(STATUS==FAIL) &
            CALL Halt(' Failed to close an HDF file with HDF_CurrentID=' &
            //'<'//TRIM(IntToChar(FileID))//'>. ') 
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
       STATUS=HDF5CloseGroup(GroupID)
       IF(STATUS==FAIL) &
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
    Meta%Status=FAIL
    NC=StringLen(Meta%VarName)
    CALL HDF5OpenData(HDF_CurrentID,NC,Char2Ints(NC,Meta%VarName), &
         Meta%DataId,Meta%DataSpc)
    IF(Meta%DataId==FAIL)THEN
       IF(PRESENT(Put_O))THEN
          CALL CreateData(Meta)
       ELSE
          CALL Halt(' Failed in OpenData:'//TRIM(MetaChar(Meta)))
       ENDIF
    ELSE
       SizeOf=HDF5SizeOfData(Meta%DataSpc)
       IF(PRESENT(Put_O))THEN
          IF(SizeOf<Meta%Dimension.AND.Meta%UnLimited==1)THEN
             CALL HDF5ExtendData(Meta%DataId,Meta%DataSpc,Meta%Dimension)
          ELSEIF(SizeOf>Meta%Dimension)THEN
             CALL HDF5SelectData(Meta%DataId,Meta%DataSpc,Meta%Dimension)
          ELSEIF(SizeOf<Meta%Dimension)THEN
             CALL Halt(' Failed in OpenData, need to redimension fixed form data:' &
                  //TRIM(MetaChar(Meta)))
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
    IF(Meta%Status==FAIL) &
         CALL Halt(' HDF R/W error for '//TRIM(MetaChar(Meta)))
    Meta%Status=HDF5CloseData(Meta%DataId,Meta%DataSpc)
    IF(Meta%Status==FAIL) &
         CALL Halt(' HDF5CloseData error for '//TRIM(MetaChar(Meta)))
  END SUBROUTINE CloseData
  !===============================================================================
  !
  !===============================================================================
  FUNCTION NameTag(VarName,Tag_O) RESULT(FullName)
    CHARACTER(LEN=*),         INTENT(IN)    :: VarName
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
    CHARACTER(LEN=DEFAULT_CHR_LEN)          :: FullName
    IF(PRESENT(Tag_O))THEN
       FullName=TRIM(VarName)//TRIM(Tag_O)
    ELSE
       FullName=TRIM(VarName)
    ENDIF
    CALL LowCase(FullName)
  END FUNCTION NameTag
  !===============================================================================
  !
  !===============================================================================
  FUNCTION SetMeta(VarName,DataType,SizeOf,UnLimit_O) RESULT(Meta)
    TYPE(META_DATA)              :: Meta
    CHARACTER(LEN=*), INTENT(IN) :: VarName
    INTEGER,          INTENT(IN) :: DataType,SizeOf
    LOGICAL,OPTIONAL, INTENT(IN) :: UnLimit_O
    IF(PRESENT(UnLimit_O))THEN
       IF(UnLimit_O)THEN
          Meta%UnLimited=1
       ELSE
          Meta%UnLimited=0
       ENDIF
    ELSE
       Meta%UnLimited=0
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
    IF(Meta%UnLimited==1)THEN
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
    CHARACTER(LEN=*),INTENT(IN)          :: String
    INTEGER,         INTENT(IN)          :: NC 
    INTEGER,DIMENSION(NC)                :: Char2Ints
    INTEGER                              :: I 
    DO I=1,NC
       Char2Ints(I)=ICHAR(String(I:I))
    ENDDO
  END FUNCTION Char2Ints
  !===============================================================================
  !
  !===============================================================================
  SUBROUTINE WriteIntegerVector(Meta,A)       
    TYPE(META_DATA)                :: Meta
    INTEGER,                   &
         DIMENSION(Meta%Dimension)   :: A
    Meta%Status=HDF5WriteIntegerVector(Meta%DataId,Meta%DataSpc,A)
  END SUBROUTINE WriteIntegerVector

  SUBROUTINE WriteDoubleVector(Meta,A)       
    TYPE(META_DATA)                :: Meta
    REAL(DOUBLE),              &
         DIMENSION(Meta%Dimension)   :: A
    Meta%Status=HDF5WriteDoubleVector(Meta%DataId,Meta%DataSpc,A)
  END SUBROUTINE WriteDoubleVector

  SUBROUTINE ReadIntegerVector(Meta,A)       
    TYPE(META_DATA)               :: Meta
    INTEGER,                    &
         DIMENSION(Meta%Dimension)   :: A
    Meta%Status=HDF5ReadIntegerVector(Meta%DataId,Meta%DataSpc,A)
  END SUBROUTINE ReadIntegerVector

  SUBROUTINE ReadDoubleVector(Meta,A)       
    TYPE(META_DATA)               :: Meta
    REAL(DOUBLE),              &
         DIMENSION(Meta%Dimension)  :: A
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
       Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32, &
            SIZE(A%I,1),.FALSE.)
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
             Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32, &
                  SIZE(A%I,1)*SIZE(A%I,2),.FALSE.)
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
             Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32, &
                  SIZE(A%I,1)*SIZE(A%I,2)*SIZE(A%I,3), &
                  .FALSE.)
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
             Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,             &
                  SIZE(A%I,1)*SIZE(A%I,2)*SIZE(A%I,3)*SIZE(A%I,4), &
                  .FALSE.)
             CALL OpenData(Meta)
             CALL ReadIntegerVector(Meta,A%I)
             CALL CloseData(Meta)
#ifdef PARALLEL 
          ENDIF
          IF(InParallel)CALL Bcast(A)
#endif 
        END SUBROUTINE Get_INT_RNK4

        SUBROUTINE Get_DBL_SCLR(A,VarName,Tag_O)
          REAL(DOUBLE),             INTENT(INOUT)   :: A
          CHARACTER(LEN=*),         INTENT(IN)      :: VarName
          CHARACTER(LEN=*),OPTIONAL,INTENT(IN)      :: Tag_O
          REAL(DOUBLE),DIMENSION(1)                 :: B
          TYPE(META_DATA)                           :: Meta
#ifdef PARALLEL 
          IF(MyId==ROOT)THEN
#endif 
             Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,1,.FALSE.)
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
             Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE, &
                  SIZE(A%D,1),.FALSE.)
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
             Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE, &
                  SIZE(A%D,1)*SIZE(A%D,2),.FALSE.)
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
             Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE, &
                  SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3),  &
                  .FALSE.)
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
             Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,            &
                  SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3)*SIZE(A%D,4), &
                  .FALSE.)
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
             Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,            &
                  SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3)*SIZE(A%D,4)  &
                  *SIZE(A%D,5)*SIZE(A%D,6),.FALSE.)
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
             IF(N>DEFAULT_CHR_LEN) &
                  CALL Halt('Static strings overrun in Get_CHR_SCLR')
             Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32, &
                  DEFAULT_CHR_LEN,.FALSE.)
             CALL OpenData(Meta)
             CALL ReadIntegerVector(Meta,B)
             DO I=1,N; A(I:I)=CHAR(B(I)); ENDDO
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
                Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,1,.FALSE.)
                CALL OpenData(Meta,.TRUE.)
                B(1)=A
                CALL WriteIntegerVector(Meta,B)
                CALL CloseData(Meta)
#ifdef PARALLEL 
             ENDIF
#endif 
           END SUBROUTINE Put_INT_SCLR

           SUBROUTINE Put_INT_VECT(A,VarName,N_O,Tag_O,UnLimit_O)
             TYPE(INT_VECT),           INTENT(IN) :: A
             CHARACTER(LEN=*),         INTENT(IN) :: VarName
             INTEGER,         OPTIONAL,INTENT(IN) :: N_O
             CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
             LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
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
                Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,N,UnLimit_O)
                CALL OpenData(Meta,.TRUE.)
                CALL WriteIntegerVector(Meta,A%I)
                CALL CloseData(Meta)
#ifdef PARALLEL 
             ENDIF
#endif 
           END SUBROUTINE Put_INT_VECT

           SUBROUTINE Put_INT_RNK2(A,VarName,N_O,Tag_O,UnLimit_O)
             TYPE(INT_RNK2),           INTENT(IN) :: A
             CHARACTER(LEN=*),         INTENT(IN) :: VarName
             INTEGER,         OPTIONAL,INTENT(IN) :: N_O
             CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
             LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
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
                Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,N,UnLimit_O)
                CALL OpenData(Meta,.TRUE.)
                CALL WriteIntegerVector(Meta,A%I)
                CALL CloseData(Meta)
#ifdef PARALLEL 
             ENDIF
#endif 
           END SUBROUTINE Put_INT_RNK2

           SUBROUTINE Put_INT_RNK3(A,VarName,N_O,Tag_O,UnLimit_O)
             TYPE(INT_RNK3),           INTENT(IN) :: A
             CHARACTER(LEN=*),         INTENT(IN) :: VarName
             INTEGER,         OPTIONAL,INTENT(IN) :: N_O
             CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
             LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
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
                Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,N,UnLimit_O)
                CALL OpenData(Meta,.TRUE.)
                CALL WriteIntegerVector(Meta,A%I)
                CALL CloseData(Meta)
#ifdef PARALLEL 
             ENDIF
#endif 
           END SUBROUTINE Put_INT_RNK3

           SUBROUTINE Put_INT_RNK4(A,VarName,N_O,Tag_O,UnLimit_O)
             TYPE(INT_RNK4),           INTENT(IN) :: A
             CHARACTER(LEN=*),         INTENT(IN) :: VarName
             INTEGER,         OPTIONAL,INTENT(IN) :: N_O
             CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
             LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
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
                Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,N,UnLimit_O)
                CALL OpenData(Meta,.TRUE.)
                CALL WriteIntegerVector(Meta,A%I)
                CALL CloseData(Meta)
#ifdef PARALLEL 
             ENDIF
#endif 
           END SUBROUTINE Put_INT_RNK4
           !-------------------------------------------------------------------------------

           !-------------------------------------------------------------------------------
           SUBROUTINE Put_DBL_SCLR(A,VarName,Tag_O)
             REAL(DOUBLE),             INTENT(IN) :: A
             CHARACTER(LEN=*),         INTENT(IN) :: VarName
             CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
             REAL(DOUBLE),DIMENSION(1)            :: B
             TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
             IF(MyId==ROOT)THEN
#endif 
                Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,1,.FALSE.)
                CALL OpenData(Meta,.TRUE.)
                B(1)=A
                CALL WriteDoubleVector(Meta,B)
                CALL CloseData(Meta)
#ifdef PARALLEL 
             ENDIF
#endif 
           END SUBROUTINE Put_DBL_SCLR

           SUBROUTINE Put_DBL_VECT(A,VarName,N_O,Tag_O,UnLimit_O)
             TYPE(DBL_VECT),           INTENT(IN) :: A
             CHARACTER(LEN=*),         INTENT(IN) :: VarName
             INTEGER,         OPTIONAL,INTENT(IN) :: N_O
             CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
             LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
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
                Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,N,UnLimit_O)
                CALL OpenData(Meta,.TRUE.)
                CALL WriteDoubleVector(Meta,A%D)
                CALL CloseData(Meta)
#ifdef PARALLEL 
             ENDIF
#endif 
           END SUBROUTINE Put_DBL_VECT
           ! 
           SUBROUTINE Put_DBL_RNK2(A,VarName,N_O,Tag_O,UnLimit_O)
             TYPE(DBL_RNK2),           INTENT(IN) :: A
             CHARACTER(LEN=*),         INTENT(IN) :: VarName
             INTEGER,         OPTIONAL,INTENT(IN) :: N_O
             CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
             LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
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
                Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,N,UnLimit_O)
                CALL OpenData(Meta,.TRUE.)
                CALL WriteDoubleVector(Meta,A%D)
                CALL CloseData(Meta)
#ifdef PARALLEL 
             ENDIF
#endif 
           END SUBROUTINE Put_DBL_RNK2
           ! 
           SUBROUTINE Put_DBL_RNK3(A,VarName,N_O,Tag_O,UnLimit_O)
             TYPE(DBL_RNK3),           INTENT(IN) :: A
             CHARACTER(LEN=*),         INTENT(IN) :: VarName
             INTEGER,         OPTIONAL,INTENT(IN) :: N_O
             CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
             LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
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
                Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,N,UnLimit_O)
                CALL OpenData(Meta,.TRUE.)
                CALL WriteDoubleVector(Meta,A%D)
                CALL CloseData(Meta)
#ifdef PARALLEL 
             ENDIF
#endif 
           END SUBROUTINE Put_DBL_RNK3
           ! 
           SUBROUTINE Put_DBL_RNK4(A,VarName,N_O,Tag_O,UnLimit_O)
             TYPE(DBL_RNK4),           INTENT(IN) :: A
             CHARACTER(LEN=*),         INTENT(IN) :: VarName
             INTEGER,         OPTIONAL,INTENT(IN) :: N_O
             CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
             LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
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
                Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,N,UnLimit_O)
                CALL OpenData(Meta,.TRUE.)
                CALL WriteDoubleVector(Meta,A%D)
                CALL CloseData(Meta)
#ifdef PARALLEL 
             ENDIF
#endif 
           END SUBROUTINE Put_DBL_RNK4
           ! 
           SUBROUTINE Put_DBL_RNK6(A,VarName,N_O,Tag_O,UnLimit_O)
             TYPE(DBL_RNK6),           INTENT(IN) :: A
             CHARACTER(LEN=*),         INTENT(IN) :: VarName
             INTEGER,         OPTIONAL,INTENT(IN) :: N_O
             CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
             LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
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
                Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,N,UnLimit_O)
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
                DO I=1,N; B(I)=ICHAR(A(I:I)); ENDDO
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
!             Get the Periodic Info
              SUBROUTINE Get_PBCInfo(PBC,Tag_O)
                TYPE(PBCInfo)                         :: PBC
                CHARACTER(LEN=*),OPTIONAL,INTENT(IN)  :: Tag_O
                CALL Get(PBC%Dimen     ,'Dimension' ,Tag_O=Tag_O)
                CALL Get(PBC%PFFMaxEll ,'PFFMaxEll' ,Tag_O=Tag_O)
                CALL Get(PBC%PFFMaxLay ,'PFFMaxLay' ,Tag_O=Tag_O)
                CALL Get(PBC%PFFOvRide ,'PFFOvRide' ,Tag_O=Tag_O)
                CALL Get(PBC%AtomW     ,'AtomWrap'  ,Tag_O=Tag_O)
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
                CALL Put(PBC%PFFMaxLay ,'PFFMaxLay' ,Tag_O=Tag_O)
                CALL Put(PBC%PFFOvRide ,'PFFOvRide' ,Tag_O=Tag_O)
                CALL Put(PBC%AtomW     ,'AtomWrap'  ,Tag_O=Tag_O)
                CALL Put(PBC%SuperCell ,'SuperCell' ,Tag_O=Tag_O)
                CALL Put(PBC%InVecForm ,'VectorForm',Tag_O=Tag_O)
                CALL Put(PBC%InAtomCrd ,'AtomicCrd' ,Tag_O=Tag_O)
                CALL Put(PBC%Translate ,'Translate' ,Tag_O=Tag_O)
                CALL Put(PBC%CellVolume,'CellVolume',Tag_O=Tag_O)
                CALL Put(PBC%Epsilon   ,'Epsilon'   ,Tag_O=Tag_O)
                CALL Put(PBC%DipoleFAC ,'DPoleFAC'  ,Tag_O=Tag_O)
                CALL Put(PBC%QupoleFAC ,'QPoleFAC'  ,Tag_O=Tag_O)
                CALL Put(PBC%AutoW     ,'AutoWrap'  ,Tag_O=Tag_O)
                CALL Put(PBC%CellCenter,'CellCenter',Tag_O=Tag_O)
                CALL Put(PBC%TransVec  ,'TransVec'  ,Tag_O=Tag_O)
                CALL Put(PBC%BoxShape  ,'BoxShape'  ,Tag_O=Tag_O)
                CALL Put(PBC%InvBoxSh  ,'InvBoxSh'  ,Tag_O=Tag_O)
                CALL Put(PBC%LatFrc    ,'LatFrc'    ,Tag_O=Tag_O)
              END SUBROUTINE Put_PBCInfo
              !-------------------------------------------------------------------------------
              !     Get some coordinates
              SUBROUTINE Get_CRDS(GM,Tag_O)
                TYPE(CRDS)                           :: GM
                CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
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
                CALL Get(GM%InAu, 'inau'         ,Tag_O=Tag_O)
                CALL New(GM)
                !-------------------------------------------------------------------------------
                !        Items that can change with geometry ...       
                CALL Get(GM%ETotal    ,'gm_etot'      ,Tag_O=Tag_O)
                CALL Get(GM%Ordrd     ,'reordered'    ,Tag_O=Tag_O)
                CALL Get(GM%AtTyp     ,'atomtype'     ,Tag_O=Tag_O)
                CALL Get(GM%AtNum     ,'atomicnumbers',Tag_O=Tag_O)
                CALL Get(GM%AtNam     ,'atomname'     ,Tag_O=Tag_O)
                CALL Get(GM%AtMMTyp   ,'mmtype'       ,Tag_O=Tag_O)
                CALL Get(GM%AtMss     ,'atomicmass'   ,Tag_O=Tag_O)
                CALL Get(GM%PBC,Tag_O=Tag_O)
                CALL Get(GM%BndBox    ,'boundingbox'  ,Tag_O=Tag_O)
                CALL Get(GM%CConstrain,'constraints'  ,Tag_O=Tag_O)
                CALL Get(GM%Carts     ,'cartesians'   ,Tag_O=Tag_O)
                CALL Get(GM%BoxCarts  ,'LatticeCoord' ,Tag_O=Tag_O)
                CALL Get(GM%Velocity  ,'Velocity'     ,Tag_O=Tag_O)
                CALL Get(GM%Gradients ,'Gradients'    ,Tag_O=Tag_O)
                CALL Get(GM%Displ     ,'Displ'        ,Tag_O=Tag_O)
                CALL Get(GM%PBCDispl  ,Tag_O='PBCDispl'//TRIM(Tag_O))
                CALL Get(GM%LatticeOnly,'LatticeOnly' ,Tag_O=Tag_O)
                CALL Get(GM%AltCount  ,'AltCount'     ,Tag_O=Tag_O)
              END SUBROUTINE Get_CRDS
              !-------------------------------------------------------------------------------
              !     Put a coordinate set

              SUBROUTINE Put_CRDS(GM,Tag_O)
                TYPE(CRDS),               INTENT(IN) :: GM
                CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
                !-------------------------------------------------------------------------------
                !        Items that should not change with geometry...
                CALL Put(GM%NAtms,'natoms'       ,Tag_O=Tag_O)
                CALL Put(GM%Confg,'configuration',Tag_O=Tag_O)
                CALL Put(GM%NElec,'nel'          ,Tag_O=Tag_O)
                CALL Put(GM%NAlph,'nelalpha'     ,Tag_O=Tag_O)
                CALL Put(GM%NBeta,'nelbeta'      ,Tag_O=Tag_O)
                CALL Put(GM%TotCh,'charge'       ,Tag_O=Tag_O)
                CALL Put(GM%NKind,'nkind'        ,Tag_O=Tag_O)
                CALL Put(GM%InAu, 'inau'         ,Tag_O=Tag_O)
                !-------------------------------------------------------------------------------
                !        Items that can change with geometry ...       
                CALL Put(GM%ETotal    ,'gm_etot'      ,Tag_O=Tag_O)
                CALL Put(GM%Ordrd     ,'reordered'    ,Tag_O=Tag_O)
                CALL Put(GM%AtTyp     ,'atomtype'     ,Tag_O=Tag_O)
                CALL Put(GM%AtNum     ,'atomicnumbers',Tag_O=Tag_O)
                CALL Put(GM%AtNam     ,'atomname'     ,Tag_O=Tag_O)
                CALL Put(GM%AtMMTyp   ,'mmtype'       ,Tag_O=Tag_O)
                CALL Put(GM%AtMss     ,'atomicmass'   ,Tag_O=Tag_O)
                CALL Put(GM%PBC,Tag_O=Tag_O)
                CALL Put(GM%BndBox    ,'boundingbox'  ,Tag_O=Tag_O)
                CALL Put(GM%CConstrain,'constraints'  ,Tag_O=Tag_O)
                CALL Put(GM%Carts     ,'cartesians'   ,Tag_O=Tag_O)
                CALL Put(GM%BoxCarts  ,'LatticeCoord' ,Tag_O=Tag_O)
                CALL Put(GM%Velocity  ,'Velocity'     ,Tag_O=Tag_O)
                CALL Put(GM%Gradients ,'Gradients'    ,Tag_O=Tag_O)
                CALL Put(GM%Displ     ,'Displ'        ,Tag_O=Tag_O)
                CALL Put(GM%PBCDispl  ,Tag_O='PBCDispl'//TRIM(Tag_O))
                CALL Put(GM%LatticeOnly,'LatticeOnly' ,Tag_O=Tag_O)
                CALL Put(GM%AltCount  ,'AltCount'     ,Tag_O=Tag_O)
              END SUBROUTINE Put_CRDS
              !-------------------------------------------------------------------------------
              !     Get a BCSR matrix

              SUBROUTINE Get_BCSR(A,Name,PFix_O,CheckPoint_O,Bcast_O)
                TYPE(BCSR),               INTENT(INOUT) :: A     
                CHARACTER(LEN=*),         INTENT(IN)    :: Name
                CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: PFix_O
                LOGICAL,         OPTIONAL,INTENT(IN)    :: CheckPoint_O
                LOGICAL,         OPTIONAL,INTENT(IN)    :: Bcast_O
                REAL(DOUBLE)                            :: Dummy
                CHARACTER(LEN=DEFAULT_CHR_LEN)          :: FileName          
                INTEGER                                 :: I,NAtms,NBlks,NNon0,IOS
                LOGICAL                                 :: Exists,LimitsQ
                LOGICAL                                 :: Bcast



                IF(PRESENT(Bcast_O)) THEN
                   Bcast = Bcast_O
                ELSE
                   Bcast = .FALSE.
                ENDIF
#ifdef PARALLEL
                IF(MyId==0)THEN
#endif
                   IF(PRESENT(CheckPoint_O))THEN
                      IF(CheckPoint_O)THEN
                         CALL Get(A%NAtms,TRIM(Name)//'%NAtms')
                         CALL Get(A%NBlks,TRIM(Name)//'%NBlks')
                         CALL Get(A%NNon0,TRIM(Name)//'%NNon0')

                         IF(AllocQ(A%Alloc))THEN
                            LimitsQ=.NOT.                   &
                                 (A%NAtms<=SIZE(A%RowPt%I)).AND. &
                                 (A%NBlks<=SIZE(A%ColPt%I)).AND. &
                                 (A%NBlks<=SIZE(A%BlkPt%I)).AND. &
                                 (A%NNon0<=SIZE(A%MTrix%D))
                            IF(LimitsQ)THEN
                               CALL Delete(A)
                               CALL New(A,(/A%NAtms,A%NBlks,A%NNon0/))
                            ENDIF
                         ELSE
                            CALL New(A,(/A%NAtms,A%NBlks,A%NNon0/))
                         ENDIF
!
                         CALL Get(A%RowPt,TRIM(Name)//'%RowPt')
                         CALL Get(A%ColPt,TRIM(Name)//'%ColPt')
                         CALL Get(A%BlkPt,TRIM(Name)//'%BlkPt')
                         CALL Get(A%MTrix,TRIM(Name)//'%MTrix')
#ifdef PARALLEL
                         IF(Bcast) THEN
                            CALL BcastBCSR(A)
                         ENDIF
#endif
!!$
!!$                         WRITE(*,*)' IN GET, ATMS = ',A%NAtms,' Blks = ',A%NBlks,' NNon0 = ',A%NNon0
!!$                         WRITE(*,*)A%RowPt%I
!!$                         WRITE(*,*)A%ColPt%I
!!$                         WRITE(*,*) A%MTrix%D(1:A%NNon0)

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
                      CALL Halt(' Get_BCSR could not find '//TRIM(FileName))
                   ENDIF

#ifdef FORMATTED
                   READ(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)NAtms,NNon0,NBlks
                   INCLUDE 'Formats.Inc'
#else
                   READ(UNIT=Seq,Err=1,IOSTAT=IOS)NAtms,NNon0,NBlks
#endif
                   IF(AllocQ(A%Alloc))THEN
                      LimitsQ=.NOT.                         &
                           (NAtms<=SIZE(A%RowPt%I)).AND. &
                           (NBlks<=SIZE(A%ColPt%I)).AND. &
                           (NBlks<=SIZE(A%BlkPt%I)).AND. &
                           (NNon0<=SIZE(A%MTrix%D))
                      IF(LimitsQ)THEN
                         CALL Delete(A)
                         CALL New(A,(/NAtoms,NBlks,NNon0/))
                      ELSE
                         A%NAtms=NAtms
                         A%NBlks=NBlks
                         A%NNon0=NNon0
                      ENDIF
                   ELSE
                      CALL New(A,(/NAtms,NBlks,NNon0/))
                   ENDIF
#ifdef FORMATTED
                   READ(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%RowPt%I(I),I=1,NAtoms+1)
                   READ(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%ColPt%I(I),I=1,A%NBlks)
                   READ(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%BlkPt%I(I),I=1,A%NBlks)
                   READ(UNIT=Seq,FMT=66,Err=1,IOSTAT=IOS)(Dummy,i=1,A%NBlks)
                   READ(UNIT=Seq,FMT=66,Err=1,IOSTAT=IOS)(A%MTrix%D(I),I=1,A%NNon0)
#else
                   READ(UNIT=Seq,Err=1,IOSTAT=IOS)(A%RowPt%I(I),I=1,NAtoms+1)
                   READ(UNIT=Seq,Err=1,IOSTAT=IOS)(A%ColPt%I(I),I=1,A%NBlks)
                   READ(UNIT=Seq,Err=1,IOSTAT=IOS)(A%BlkPt%I(I),I=1,A%NBlks)
                   READ(UNIT=Seq,Err=1,IOSTAT=IOS)(Dummy,i=1,A%NBlks)
                   READ(UNIT=Seq,Err=1,IOSTAT=IOS)(A%MTrix%D(I),I=1,A%NNon0)
#endif

                   CLOSE(UNIT=Seq,STATUS='KEEP')
#ifdef PARALLEL
                ENDIF
                IF(Bcast) THEN
                   CALL BcastBCSR(A)
                ENDIF
#endif
                RETURN
1               CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Get_BCSR.')
              END SUBROUTINE Get_BCSR


              !-------------------------------------------------------------------------------
#ifdef PARALLEL
              SUBROUTINE BcastBCSR(A)
                TYPE(BCSR) :: A
                INTEGER :: IErr,NAtms,NBlks,NNon0
                LOGICAL :: LimitsQ
                NAtms = A%NAtms
                NBlks = A%NBlks
                NNon0 = A%NNon0

                CALL Bcast(NAtms)
                CALL Bcast(NBlks)
                CALL Bcast(NNon0)


                IF(AllocQ(A%Alloc))THEN
                   LimitsQ=.NOT.                         &
                        (NAtms<=SIZE(A%RowPt%I)).AND. &
                        (NBlks<=SIZE(A%ColPt%I)).AND. &
                        (NBlks<=SIZE(A%BlkPt%I)).AND. &
                        (NNon0<=SIZE(A%MTrix%D))
                   IF(LimitsQ)THEN
                      CALL Delete(A)
                      CALL New(A,(/NAtms,NBlks,NNon0/),OnAll_O=.TRUE.)
                   ENDIF
                ELSE
                   CALL New(A,(/NAtms,NBlks,NNon0/),OnAll_O=.TRUE.)
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
                         CALL Put(A%NAtms,TRIM(Name)//'%NAtms')
                         CALL Put(A%NBlks,TRIM(Name)//'%NBlks')
                         CALL Put(A%NNon0,TRIM(Name)//'%NNon0')
                         CALL Put(A%RowPt,TRIM(Name)//'%RowPt',A%NAtms+1)
#ifdef PARALLEL_CLONES
                         CALL Put(A%ColPt,TRIM(Name)//'%ColPt',A%NBlks)!,UnLimit_O=.TRUE.)
                         CALL Put(A%BlkPt,TRIM(Name)//'%BlkPt',A%NBlks)!,UnLimit_O=.TRUE.)
                         CALL Put(A%MTrix,TRIM(Name)//'%MTrix',A%NNon0)!,UnLimit_O=.TRUE.)
#else
                         CALL Put(A%ColPt,TRIM(Name)//'%ColPt',A%NBlks,UnLimit_O=.TRUE.)
                         CALL Put(A%BlkPt,TRIM(Name)//'%BlkPt',A%NBlks,UnLimit_O=.TRUE.)
                         CALL Put(A%MTrix,TRIM(Name)//'%MTrix',A%NNon0,UnLimit_O=.TRUE.)
#endif
!!$
!!$                         WRITE(*,*)' IN PUT, ATMS = ',A%NAtms,' Blks = ',A%NBlks,' NNon0 = ',A%NNon0
!!$                         WRITE(*,*) A%RowPt%I
!!$                         WRITE(*,*) A%ColPt%I
!!$                         WRITE(*,*) A%MTrix%D(1:A%NNon0)

                         RETURN
                      ENDIF
                   ENDIF

                   !           WRITE TO BINARY FILE
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
                   WRITE(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)NAtoms,A%NNon0,A%NBlks
                   WRITE(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%RowPt%I(i),i=1,NAtoms+1)
                   WRITE(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%ColPt%I(i),i=1,A%NBlks)
                   WRITE(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%BlkPt%I(i),i=1,A%NBlks)
                   WRITE(UNIT=Seq,FMT=66,Err=1,IOSTAT=IOS)(BIG_DBL,i=1,A%NBlks)
                   WRITE(UNIT=Seq,FMT=66,Err=1,IOSTAT=IOS)(A%MTrix%D(i),i=1,A%NNon0)
                   INCLUDE 'Formats.Inc'            
#else
                   WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)NAtoms,A%NNon0,A%NBlks
                   WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%RowPt%I(i),i=1,NAtoms+1)
                   WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%ColPt%I(i),i=1,A%NBlks)
                   WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%BlkPt%I(i),i=1,A%NBlks)
                   WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(BIG_DBL,i=1,A%NBlks)
                   WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%MTrix%D(i),i=1,A%NNon0)
#endif
                   CLOSE(UNIT=Seq,STATUS='KEEP')
#ifdef PARALLEL
                ENDIF
#endif
                RETURN
1               CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Put_BCSR.')
              END SUBROUTINE Put_BCSR
#ifdef PARALLEL
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
                INTEGER                          :: SCFCycle,IOS,I,NExpt,NDist,NCoef
                CHARACTER(LEN=*)                 :: Name 
                CHARACTER(LEN=DEFAULT_CHR_LEN)   :: FileName
                LOGICAL                          :: Exists
                LOGICAL,OPTIONAL,INTENT(IN)      :: Bcast_O
                LOGICAL                          :: BcastQ

#ifdef PARALLEL
                IF(PRESENT(Bcast_O)) THEN
                   BcastQ = Bcast_O
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

                   ! Allocate Memory

                   READ(UNIT=Seq,Err=100,IOSTAT=IOS) NExpt,NDist,NCoef
#ifdef PARALLEL
                ENDIF
                IF(BcastQ) THEN
                   CALL Bcast(NExpt)
                   CALL Bcast(NDist)
                   CALL Bcast(NCoef)
                ENDIF
#endif
                CALL New_HGRho(A,(/NExpt,NDist,NCoef/))

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
                   READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Co%D(I)    ,I=1,A%NCoef)
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
                   Call Bcast(A%Co,N_O=A%NCoef)
                ENDIF
                IF(MyID == ROOT) THEN
#endif

                   Close(UNIT=Seq,STATUS='KEEP')
#ifdef PARALLEL
                ENDIF
#endif
                RETURN
100             CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Get_HGRho.')
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
                !   
                FileName=TrixFile(Name,Args,SCFCycle)
                INQUIRE(FILE=FileName,EXIST=Exists)
                IF(Exists) THEN
                   OPEN(UNIT=Seq,FILE=FileName,STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
                ELSE
                   OPEN(UNIT=Seq,FILE=FileName,STATUS='NEW',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
                ENDIF

                !   Write density to disk

                WRITE(UNIT=Seq,Err=100,IOSTAT=IOS) A%NExpt,A%NDist,A%NCoef
                WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%NQ%I(I)    ,I=1,A%NExpt)
                WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%OffQ%I(I)  ,I=1,A%NExpt)
                WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%OffR%I(I)  ,I=1,A%NExpt)
                WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Lndx%I(I)  ,I=1,A%NExpt)
                WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Expt%D(I)  ,I=1,A%NExpt)
                WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qx%D(I)    ,I=1,A%NDist)
                WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qy%D(I)    ,I=1,A%NDist)
                WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qz%D(I)    ,I=1,A%NDist)
                WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Co%D(I)    ,I=1,A%NCoef)

                CLOSE(UNIT=Seq,STATUS='KEEP')
                RETURN
100             CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Put_HGRho.')
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
11              CALL Halt(' OpenASCII ERROR: IOS='//TRIM(IntToChar(IOS))// &
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
                      DO I=1,N; B(I)=ICHAR(A%C(II)(I:I)); ENDDO
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
                        DO I=1,N; A%C(II)(I:I)=CHAR(B(I)); ENDDO
                      ENDDO
#else
                      NN = SIZE(A%C)
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
                INTEGER                              :: I,N,II,NN
                TYPE(CHR10_VECT) :: A
                CHARACTER(LEN=*),         INTENT(IN) :: VarName
                CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
                TYPE(META_DATA)                      :: Meta

#ifdef OLD_CHR10_VECT
                INTEGER,DIMENSION(DCL)   :: B  !=ICHAR(' ')
#else
                INTEGER,ALLOCATABLE :: B(:)
                INTEGER :: RunInd,BufSize
#endif

#ifdef PARALLEL
                IF(MyId==ROOT)THEN
#endif

#ifdef OLD_CHR10_VECT
                   NN=SIZE(A%C)
                   DO II = 1, NN
                      N=LEN(A%C(II))
                      IF(N>DCL) CALL Halt('Static strings overrun in Put_CHR10_VECT')
                      DO I=1,N; B(I)=ICHAR(A%C(II)(I:I)); ENDDO
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
                   CALL Put(BufSize,NameTag(VarName,TRIM(IntToChar(0))//TRIM(Tag_O)))
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
                 END SUBROUTINE Put_CHR10_VECT

                 !-------------------------------------------------------------------------------
                 SUBROUTINE Get_CHR10_VECT(A,VarName,Tag_O)
                   INTEGER                                 :: I,N,II,NN
                   TYPE(CHR10_VECT),           INTENT(INOUT) :: A
                   CHARACTER(LEN=*),         INTENT(IN)    :: VarName
                   CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
#ifdef OLD_CHR10_VECT
                   INTEGER,DIMENSION(DEFAULT_CHR_LEN)      :: B !=ICHAR(' ')
#else
                   INTEGER,ALLOCATABLE :: B(:)
                   INTEGER :: RunInd,StrInd,StrLen,BufSize
                   CHARACTER(LEN=DCL) :: TEMP
#endif
                   TYPE(META_DATA)                         :: Meta

#ifdef OLD_CHR10_VECT
#else
                   CALL Get(BufSize,NameTag(VarName,TRIM(IntToChar(0))//TRIM(Tag_O)))
#endif

#ifdef PARALLEL 

                   IF(MyId==ROOT)THEN
#endif 

#ifdef OLD_CHR10_VECT
                      NN=SIZE(A%C)
                      DO II = 1, NN
                        N=LEN(A%C(II))
                        IF(N>DEFAULT_CHR_LEN) CALL Halt('Static strings overrun in Get_CHR10_VECT')
                        Meta=SetMeta(NameTag(VarName,TRIM(IntToChar(II))//TRIM(Tag_O)),NATIVE_INT32,N,.FALSE.)

                        CALL OpenData(Meta)
                        CALL ReadIntegerVector(Meta,B(1))
                        CALL CloseData(Meta)
                        DO I=1,N; A%C(II)(I:I)=CHAR(B(I)); ENDDO
                      ENDDO
#else
                      NN = SIZE(A%C)
                      Meta=SetMeta(NameTag(VarName,TRIM(IntToChar(1))//TRIM(Tag_O)),NATIVE_INT32,BufSize,.FALSE.)
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
#endif

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
                         Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32, &
                              SIZE(ILog%I,1),.FALSE.)
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
    LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
!
    IF(PRESENT(Name_O))THEN
       CALL Put(CS%Radius   ,TRIM(Name_O)//'_cell_radius',Tag_O=Tag_O)
       CALL Put(CS%NCells   ,TRIM(Name_O)//'_cell_number',Tag_O=Tag_O)
       CALL Put(CS%CellCarts,TRIM(Name_O)//'_cell_vectors',Tag_O=Tag_O,Unlimit_O=Unlimit_O)
    ELSE
       CALL Put(CS%Radius   ,'cell_radius',Tag_O=Tag_O)
       CALL Put(CS%NCells   ,'cell_number',Tag_O=Tag_O)
       CALL Put(CS%CellCarts,'cell_vectors',Tag_O=Tag_O,Unlimit_O=Unlimit_O)
    ENDIF
  END SUBROUTINE Put_CellSet

  SUBROUTINE Get_CellSet(CS,Name_O,Tag_O)
    TYPE(CellSet)                  :: CS
    CHARACTER(Len=*),OPTIONAL      :: Name_O
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
    INTEGER                        :: NC
!
    IF(PRESENT(Name_O))THEN
       CALL Get(CS%Radius   ,TRIM(Name_O)//'_cell_radius',Tag_O=Tag_O)
       CALL Get(CS%NCells   ,TRIM(Name_O)//'_cell_number',Tag_O=Tag_O)
       CALL New_CellSet(CS,CS%NCells)
       CALL Get(CS%CellCarts,TRIM(Name_O)//'_cell_vectors',Tag_O=Tag_O)
    ELSE
       CALL Get(CS%Radius   ,'cell_radius',Tag_O=Tag_O)
       CALL Get(CS%NCells   ,'cell_number',Tag_O=Tag_O)
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
END MODULE

