!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 1999, The University of California
!
!
!    GENERIC IO ROUTINES FOR MONDOSCF TYPES
!
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
                       Get_DBL_RNK4, Get_CHR_SCLR,  &
                       Get_LOG_SCLR, Get_BSET,     Get_CRDS,     &
                       Get_TOLS,     Get_BCSR,                   & 
#ifdef PARALLEL
                       Get_DBCSR,                                &
#endif
                       Get_ARGMT,    Get_HGRho

   END INTERFACE     
   INTERFACE Put
      MODULE PROCEDURE Put_INT_SCLR, Put_INT_VECT, Put_INT_RNK2, &
                       Put_INT_RNK3, Put_INT_RNK4, Put_DBL_SCLR, &
                       Put_DBL_VECT, Put_DBL_RNK2, Put_DBL_RNK3, &
                       Put_DBL_RNK4, Put_DBL_RNK6, Put_CHR_SCLR, &
                       Put_LOG_SCLR, Put_BSET,     Put_CRDS,     &
#ifdef PARALLEL
                       Put_DBCSR,                                &
#endif
                       Put_TOLS,     Put_BCSR,     Put_HGRho
   END INTERFACE     
   TYPE META_DATA
      INTEGER            :: Dimension
      INTEGER            :: DataType
      INTEGER            :: DataId
      INTEGER            :: DataSpc
      INTEGER            :: Status
      INTEGER            :: Unlimited
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: VarName 
   END TYPE
   INTEGER, SAVE      :: FileID
   INTEGER, EXTERNAL  :: HDF5CreateFile,HDF5OpenFile,HDF5CloseFile,         &
                         HDF5CloseData,HDF5ReadIntegerVector,               &
                         HDF5ReadDoubleVector,HDF5WriteIntegerVector,       &
                         HDF5WriteDoubleVector,HDF5SizeOfData
   INTEGER, PARAMETER :: NATIVE_DOUBLE=6, NATIVE_INT32=24
   CONTAINS 
!================================================================
!
!  FILE MANIPULATION ROUTINES
!
!================================================================
!
!     Initialize a HDF file
!
      SUBROUTINE InitHDF(FileName)
         CHARACTER(LEN=*), INTENT(IN) :: FileName
         INTEGER                      :: NC,STATUS
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            NC=StringLen(FileName)
            FileID=HDF5CreateFile(NC,Char2Ints(NC,FileName))
            IF(FileID==FAIL) &
            CALL Halt(' Failed to create the HDF file' &
                      //'<'//TRIM(FileName)//'>.') 
#ifdef PARALLEL
         ENDIF
#endif  
      END SUBROUTINE InitHDF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    Open a HDF file
! 
      SUBROUTINE OpenHDF(FileName)
         CHARACTER(LEN=*),INTENT(IN) :: FileName
         INTEGER                     :: NC,STATUS
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            NC=StringLen(FileName)
            FileID=HDF5OpenFile(NC,Char2Ints(NC,FileName))
            IF(FileID==FAIL) &
            CALL Halt(' Failed to open the HDF file' &
                      //'<'//TRIM(FileName)//'>.') 
#ifdef PARALLEL
         ENDIF
#endif  
      END SUBROUTINE OpenHDF
!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
      SUBROUTINE CloseHDF()
         INTEGER :: Status
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            STATUS=HDF5CloseFile(FileID)
            IF(FileID==FAIL) &
            CALL Halt(' Failed to close an HDF file with FileID=' &
                      //'<'//TRIM(IntToChar(FileID))//'>.') 
#ifdef PARALLEL
         ENDIF
#endif  
      END SUBROUTINE CloseHDF
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------  
      SUBROUTINE CreateData(Meta)
         TYPE(META_DATA) :: Meta
         INTEGER         :: NC
         NC=StringLen(Meta%VarName)
         CALL HDF5CreateData(FileID,Meta%DataType,Meta%Dimension,          &
                             NC,Char2Ints(NC,Meta%VarName),Meta%Unlimited, &
                             Meta%DataId,Meta%DataSpc)
         IF(Meta%DataId==FAIL) &
         CALL Halt(' Failed in CreateData:'//TRIM(MetaChar(Meta)))
         Meta%Status=SUCCEED
      END SUBROUTINE CreateData
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      SUBROUTINE OpenData(Meta,Put_O)
         TYPE(META_DATA)  :: Meta
         LOGICAL,OPTIONAL :: Put_O
         INTEGER          :: I,NC,SizeOf,SIZE_OF_HDF5_DATA
         Meta%Status=FAIL
         NC=StringLen(Meta%VarName)
         CALL HDF5OpenData(FileID,NC,Char2Ints(NC,Meta%VarName), &
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
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
      SUBROUTINE CloseData(Meta)
         TYPE(META_DATA) :: Meta
         IF(Meta%Status==FAIL) &
         CALL Halt(' HDF R/W error for '//TRIM(MetaChar(Meta)))
         Meta%Status=HDF5CloseData(Meta%DataId,Meta%DataSpc)
         IF(Meta%Status==FAIL) &
         CALL Halt(' HDF5CloseData error for '//TRIM(MetaChar(Meta)))
      END SUBROUTINE CloseData
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
! 
!----------------------------------------------------------------------
      FUNCTION StringLen(String)
         CHARACTER(LEN=*),INTENT(IN) :: String
         INTEGER                     :: StringLen
         StringLen=LEN(TRIM(String))
      END FUNCTION StringLen
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
      FUNCTION Char2Ints(NC,String)
         CHARACTER(LEN=*),INTENT(IN)          :: String
         INTEGER,         INTENT(IN)          :: NC 
         INTEGER,DIMENSION(NC)                :: Char2Ints
         INTEGER                              :: I 
         DO I=1,NC
            Char2Ints(I)=ICHAR(String(I:I))
         ENDDO
      END FUNCTION Char2Ints
!-----------------------------------------------------------------------
!  
!  
!  
      SUBROUTINE WriteIntegerVector(Meta,A)       
         TYPE(META_DATA)                :: Meta
         INTEGER,                   &
            DIMENSION(Meta%Dimension)   :: A
         Meta%Status=HDF5WriteIntegerVector(Meta%DataId,Meta%DataSpc,A)
      END SUBROUTINE WriteIntegerVector
!
      SUBROUTINE WriteDoubleVector(Meta,A)       
         TYPE(META_DATA)                :: Meta
         REAL(DOUBLE),              &
            DIMENSION(Meta%Dimension)   :: A
         Meta%Status=HDF5WriteDoubleVector(Meta%DataId,Meta%DataSpc,A)
      END SUBROUTINE WriteDoubleVector
!
      SUBROUTINE ReadIntegerVector(Meta,A)       
         TYPE(META_DATA)               :: Meta
         INTEGER,                    &
           DIMENSION(Meta%Dimension)   :: A
         Meta%Status=HDF5ReadIntegerVector(Meta%DataId,Meta%DataSpc,A)
      END SUBROUTINE ReadIntegerVector
!
      SUBROUTINE ReadDoubleVector(Meta,A)       
         TYPE(META_DATA)               :: Meta
         REAL(DOUBLE),              &
            DIMENSION(Meta%Dimension)  :: A
         Meta%Status=HDF5ReadDoubleVector(Meta%DataId,Meta%DataSpc,A)
      END SUBROUTINE ReadDoubleVector
!--------------------------------------------------------------------
!
!
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
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_INT_SCLR
!--------------------------------------------------------------------
!
!
    SUBROUTINE Get_INT_VECT(A,VarName,Tag_O)
       TYPE(INT_VECT),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) &
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32, &
                       SIZE(A%I,1),.FALSE.)
          CALL OpenData(Meta)
          CALL ReadIntegerVector(Meta,A%I)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_INT_VECT
!--------------------------------------------------------------------
!
!
    SUBROUTINE Get_INT_RNK2(A,VarName,Tag_O)
       TYPE(INT_RNK2),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) &
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32, &
                       SIZE(A%I,1)*SIZE(A%I,2),.FALSE.)
          CALL OpenData(Meta)
          CALL ReadIntegerVector(Meta,A%I)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_INT_RNK2
!--------------------------------------------------------------------
!
!
    SUBROUTINE Get_INT_RNK3(A,VarName,Tag_O)
       TYPE(INT_RNK3),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) &
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32, &
                       SIZE(A%I,1)*SIZE(A%I,2)*SIZE(A%I,3), &
                       .FALSE.)
          CALL OpenData(Meta)
          CALL ReadIntegerVector(Meta,A%I)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_INT_RNK3
!--------------------------------------------------------------------
!
!
    SUBROUTINE Get_INT_RNK4(A,VarName,Tag_O)
       TYPE(INT_RNK4),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) &
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,             &
                       SIZE(A%I,1)*SIZE(A%I,2)*SIZE(A%I,3)*SIZE(A%I,4), &
                       .FALSE.)
          CALL OpenData(Meta)
          CALL ReadIntegerVector(Meta,A%I)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_INT_RNK4
!
    SUBROUTINE Get_DBL_SCLR(A,VarName,Tag_O)
       REAL(DOUBLE),             INTENT(INOUT)   :: A
       CHARACTER(LEN=*),         INTENT(IN)      :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)      :: Tag_O
       REAL(DOUBLE),DIMENSION(1)                 :: B
       TYPE(META_DATA)                           :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,1,.FALSE.)
          CALL OpenData(Meta)
          CALL ReadDoubleVector(Meta,B)
          CALL CloseData(Meta)
          A=B(1)
#ifdef PARALLEL 
       ENDIF       
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_DBL_SCLR
!
    SUBROUTINE Get_DBL_VECT(A,VarName,Tag_O)
       TYPE(DBL_VECT),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) &
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE, &
                       SIZE(A%D,1),.FALSE.)
          CALL OpenData(Meta)
          CALL ReadDoubleVector(Meta,A%D)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_DBL_VECT
!
    SUBROUTINE Get_DBL_RNK2(A,VarName,Tag_O)
       TYPE(DBL_RNK2),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) &
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE, &
                       SIZE(A%D,1)*SIZE(A%D,2),.FALSE.)
          CALL OpenData(Meta)
          CALL ReadDoubleVector(Meta,A%D)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_DBL_RNK2
!
    SUBROUTINE Get_DBL_RNK3(A,VarName,Tag_O)
       TYPE(DBL_RNK3),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) &
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE, &
                       SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3),  &
                       .FALSE.)
          CALL OpenData(Meta)
          CALL ReadDoubleVector(Meta,A%D)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_DBL_RNK3
!
    SUBROUTINE Get_DBL_RNK4(A,VarName,Tag_O)
       TYPE(DBL_RNK4),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) &
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,            &
                       SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3)*SIZE(A%D,4), &
                       .FALSE.)
          CALL OpenData(Meta)
          CALL ReadDoubleVector(Meta,A%D)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_DBL_RNK4
!
    SUBROUTINE Get_DBL_RNK6(A,VarName,Tag_O)
       TYPE(DBL_RNK6),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) &
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,            &
                       SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3)*SIZE(A%D,4)  &
                      *SIZE(A%D,5)*SIZE(A%D,6),.FALSE.)
          CALL OpenData(Meta)
          CALL ReadDoubleVector(Meta,A%D)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_DBL_RNK6
!---------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------------
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
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_CHR_SCLR
!---------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------------
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
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_LOG_SCLR
!---------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------------
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
!
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
!
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
!
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
!
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
!---------------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------------------
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
!
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
!----------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
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
!----------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
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
!----------------------------------------------------------------------------------------
!     Get a basis set
!----------------------------------------------------------------------------------------
      SUBROUTINE Get_BSET(BS,Tag_O)
         TYPE(BSET),              INTENT(OUT) :: BS
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
      END SUBROUTINE Get_BSET
!------------------------------------------------------------------
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
      END SUBROUTINE Put_BSET
!-------------------------------------------------------------------
!     Get some coordinates
!
      SUBROUTINE Get_CRDS(GM,Tag_O)
         TYPE(CRDS),           INTENT(OUT)    :: GM
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
         IF(AllocQ(GM%Alloc))CALL Delete(GM)
!-------------------------------------------------------------------
!        Items that should not change with geometry...
!
         CALL Get(GM%NAtms,'natoms')
         CALL Get(GM%Confg,'configuration')
         CALL Get(GM%NElec,'nel')
         CALL Get(GM%NAlph,'nelalpha')
         CALL Get(GM%NBeta,'nelbeta')
         CALL Get(GM%TotCh,'charge')
         CALL Get(GM%NKind,'nkind')
         CALL Get(GM%InAu, 'inau')
#ifdef PERIODIC
         CALL Get(GM%AtomW     ,'AtomWrap')
         CALL Get(GM%AutoW(1)  ,'AutoWrapX')
         CALL Get(GM%AutoW(2)  ,'AutoWrapY')
         CALL Get(GM%AutoW(3)  , 'AutoWrapZ')
         CALL Get(GM%InVecForm ,'VectorForm')
         CALL Get(GM%InAtomCrd ,'AtomicCrd')
         CALL Get(GM%NoTransVec,'NoTranVec')
#endif
         CALL New(GM)
!----------------------------------------------------------------
!        Items that can change with geometry ...       
!
         CALL Get(GM%Ordrd,'reordered',Tag_O=Tag_O)
         CALL Get(GM%AtTyp,'atomtype',Tag_O=Tag_O)
         CALL Get(GM%AtNum,'atomicnumbers',Tag_O=Tag_O)
         CALL Get(GM%Carts,'cartesians',Tag_O=Tag_O)
         CALL Get(GM%BndBox,'boundingbox',Tag_O=Tag_O)
#ifdef PERIODIC
         CALL Get(GM%TransVec,'Originvector',Tag_O=Tag_O)
         CALL Get(GM%BoxShape,'Boxshape',Tag_O=Tag_O)
         CALL Get(GM%InvBoxSh,'InverseBoxshape',Tag_O=Tag_O)
         CALL Get(GM%BoxCarts,'LatticeCoord',Tag_O=Tag_O)
#endif
      END SUBROUTINE Get_CRDS
!---------------------------------------------------------------------
!     Put a coordinate set
!
      SUBROUTINE Put_CRDS(GM,Tag_O)
         TYPE(CRDS),               INTENT(IN) :: GM
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
!-------------------------------------------------------------------
!        Items that should not change with geometry...
!
         CALL Put(GM%NAtms,'natoms')
         CALL Put(GM%Confg,'configuration')
         CALL Put(GM%NElec,'nel')
         CALL Put(GM%NAlph,'nelalpha')
         CALL Put(GM%NBeta,'nelbeta')
         CALL Put(GM%TotCh,'charge')
         CALL Put(GM%NKind,'nkind')
         CALL Put(GM%InAu, 'inau')
#ifdef PERIODIC
         CALL Put(GM%AtomW     ,'AtomWrap')
         CALL Put(GM%AutoW(1)  ,'AutoWrapX')
         CALL Put(GM%AutoW(2)  ,'AutoWrapY')
         CALL Put(GM%AutoW(3)  ,'AutoWrapZ')
         CALL Put(GM%InVecForm ,'VectorForm')
         CALL Put(GM%InAtomCrd ,'AtomicCrd')
         CALL Put(GM%NoTransVec,'NoTranVec')
#endif
!----------------------------------------------------------------
!        Items that can change with geometry ...       
!
         CALL Put(GM%Ordrd,'reordered',Tag_O=Tag_O)
         CALL Put(GM%AtNum,'atomicnumbers',Tag_O=Tag_O)
         CALL Put(GM%AtTyp,'atomtype',Tag_O=Tag_O)
         CALL Put(GM%Carts,'cartesians',Tag_O=Tag_O)
         CALL Put(GM%BndBox,'boundingbox',Tag_O=Tag_O)
#ifdef PERIODIC
         CALL Put(GM%TransVec,'Originvector',Tag_O=Tag_O)
         CALL Put(GM%BoxShape,'Boxshape',Tag_O=Tag_O)
         CALL Put(GM%InvBoxSh,'InverseBoxshape',Tag_O=Tag_O)
         CALL Put(GM%BoxCarts,'LatticeCoord',Tag_O=Tag_O)
#endif
      END SUBROUTINE Put_CRDS
!---------------------------------------------------------------------
!     Get a BCSR matrix
!
      SUBROUTINE Get_BCSR(A,Name,PFix_O,CheckPoint_O)
         TYPE(BCSR),               INTENT(INOUT) :: A     
         CHARACTER(LEN=*),         INTENT(IN)    :: Name
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: PFix_O
         LOGICAL,         OPTIONAL,INTENT(IN)    :: CheckPoint_O
         REAL(DOUBLE)                            :: Dummy
         CHARACTER(LEN=DEFAULT_CHR_LEN)          :: FileName          
         INTEGER                                 :: I,NAtms,NBlks,NNon0,IOS
         LOGICAL                                 :: Exists,LimitsQ
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
                  CALL Get(A%RowPt,TRIM(Name)//'%RowPt')
                  CALL Get(A%ColPt,TRIM(Name)//'%ColPt')
                  CALL Get(A%BlkPt,TRIM(Name)//'%BlkPt')
                  CALL Get(A%MTrix,TRIM(Name)//'%MTrix')
                  RETURN
               ENDIF
            ENDIF
!
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
#endif
         RETURN
       1 CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Get_BCSR.')
      END SUBROUTINE Get_BCSR
!---------------------------------------------------------------------
!     Put a BCSR matrix
!
      SUBROUTINE Put_BCSR(A,Name,PFix_O,BlksName_O,Non0Name_O,CheckPoint_O)
         TYPE(BCSR),               INTENT(IN) :: A             
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: PFix_O,BlksName_O,Non0Name_O
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
!
            IF(PRESENT(CheckPoint_O))THEN
               IF(CheckPoint_O)THEN
                  CALL Put(A%NAtms,TRIM(Name)//'%NAtms')
                  CALL Put(A%NBlks,TRIM(Name)//'%NBlks')
                  CALL Put(A%NNon0,TRIM(Name)//'%NNon0')
                  CALL Put(A%RowPt,TRIM(Name)//'%RowPt',A%NAtms+1)
                  CALL Put(A%ColPt,TRIM(Name)//'%ColPt',A%NBlks,UnLimit_O=.TRUE.)
                  CALL Put(A%BlkPt,TRIM(Name)//'%BlkPt',A%NBlks,UnLimit_O=.TRUE.)
                  CALL Put(A%MTrix,TRIM(Name)//'%MTrix',A%NNon0,UnLimit_O=.TRUE.)
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
!--------------------------------------------
!           KLUGE To maintain complience with old (F77) version of MONDO
            IF(PRESENT(BlksName_O))CALL Put(A%NBlks,BlksName_O)
            IF(PRESENT(Non0Name_O))CALL Put(A%NNon0,Non0Name_O)
!           End KLUGE 
!--------------------------------------------
#ifdef PARALLEL
          ENDIF
#endif
          RETURN
        1 CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Put_BCSR.')
      END SUBROUTINE Put_BCSR
#ifdef PARALLEL
!---------------------------------------------------------------------
!     Get a DBCSR matrix
!
      SUBROUTINE Get_DBCSR(A,Name,PFix_O,FromHDF_O)
         TYPE(DBCSR),              INTENT(INOUT) :: A     
         TYPE(BCSR)                              :: B                          
         CHARACTER(LEN=*),         INTENT(IN)    :: Name
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: PFix_O
         LOGICAL,         OPTIONAL,INTENT(IN)    :: FromHDF_O
!----------------------------------------------------------------------------
!        KLUGE To maintain compliance with non-parallel f77 version of MONDO
!        Ultimately, should do distributed IO
         CALL Get_BCSR(B,Name,PFix_O,FromHDF_O)
         CALL SetEq(A,B) 
         CALL Delete(B)
         A%Node=MyId
!        End KLUGE 
!--------------------------------------------
      END SUBROUTINE Get_DBCSR
!---------------------------------------------------------------------
!     Put a DBCSR matrix
!
      SUBROUTINE Put_DBCSR(A,Name,PFix_O,ToHDF_O,BlksName_O,Non0Name_O)
         TYPE(DBCSR), INTENT(INOUT)           :: A     
         TYPE(BCSR)                           :: B                          
         CHARACTER(LEN=*),         INTENT(IN) :: Name
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: PFix_O,BlksName_O,Non0Name_O
         LOGICAL, OPTIONAL,        INTENT(IN) :: ToHDF_O
!----------------------------------------------------------------------------
!        KLUGE To maintain compliance with non-parallel f77 version of MONDO
!        Ultimately, should do distributed IO
!
         CALL SetEq(B,A)
         CALL Put_BCSR(B,Name,PFix_O,ToHDF_O,BlksName_O,Non0Name_O)
         CALL Delete(B)
!        End KLUGE 
!--------------------------------------------
      END SUBROUTINE Put_DBCSR
#endif
!------------------------------------------------------------------
!     Put thresholds 
!
      SUBROUTINE Put_TOLS(NGLCT,Tag_O)
         TYPE(TOLS),               INTENT(IN) :: NGLCT
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
         CHARACTER(LEN=DEFAULT_CHR_LEN)       :: Tag
         IF(PRESENT(Tag_O))THEN
            Tag='_'//TRIM(Tag_O)
         ELSE
            Tag=''
         ENDIF 
         CALL Put(NGLCT%Cube,'cubeneglect',Tag_O=Tag_O)
         CALL Put(NGLCT%Trix,'trixneglect',Tag_O=Tag_O)
         CALL Put(NGLCT%Dist,'distneglect',Tag_O=Tag_O)
         CALL Put(NGLCT%TwoE,'twoeneglect',Tag_O=Tag_O)
         CALL Put(NGLCT%ETol,'enregyneglect',Tag_O=Tag_O)
         CALL Put(NGLCT%DTol,'densityneglect',Tag_O=Tag_O)
      END SUBROUTINE Put_TOLS

!------------------------------------------------------------------
!     Get thresholds 
!
      SUBROUTINE Get_TOLS(NGLCT,Tag_O)
         TYPE(TOLS),              INTENT(OUT) :: NGLCT
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
!        CHARACTER(LEN=DEFAULT_CHR_LEN)       :: Tag
!        IF(PRESENT(Tag_O))THEN
!           Tag='_'//TRIM(Tag_O)
!        ELSE
!           Tag=''
!        ENDIF 
         CALL Get(NGLCT%Cube,'cubeneglect',Tag_O=Tag_O)
         CALL Get(NGLCT%Trix,'trixneglect',Tag_O=Tag_O)
         CALL Get(NGLCT%Dist,'distneglect',Tag_O=Tag_O)
         CALL Get(NGLCT%TwoE,'twoeneglect',Tag_O=Tag_O)
         CALL Get(NGLCT%ETol,'enregyneglect',Tag_O=Tag_O)
         CALL Get(NGLCT%DTol,'densityneglect',Tag_O=Tag_O)
      END SUBROUTINE Get_TOLS
!--------------------------------------------------------------------
!     Get arguments from the command line
!
      SUBROUTINE Get_ARGMT(A)
#ifdef NAG
         USE F90_UNIX
#else
         INTEGER,EXTERNAL               :: IARGC
#endif
         TYPE(ARGMT),INTENT(OUT)        :: A
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
            CALL BCast(NInts)
            CALL BCast(NChar)
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
            CALL BCast(A%I)
            CALL BCast_CHR_VECT(A%C) ! Stupid PGI compiler ...
         ENDIF
#endif
      END SUBROUTINE Get_ARGMT

  SUBROUTINE Get_HGRho(A,Name,Args,SCFCycle,REst_O)
    TYPE(HGRho)                      :: A
    TYPE(ARGMT)                      :: Args
    INTEGER                          :: SCFCycle,IOS,I,NExpt,NDist,NCoef
    REAL(DOUBLE)                     :: Dummy
    CHARACTER(LEN=*)                 :: Name 
    CHARACTER(LEN=DEFAULT_CHR_LEN)   :: FileName
    LOGICAL,OPTIONAL                 :: REst_O
    LOGICAL                          :: Exists
!   
    FileName=TrixFile(Name,Args,SCFCycle)
    INQUIRE(FILE=FileName,EXIST=Exists)
    IF(Exists) THEN
       OPEN(UNIT=Seq,FILE=FileName,STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    ELSE
       CALL Halt(' Get_HGRho could not find '//TRIM(FileName))
    ENDIF
!
!   Allocate Memory
!
    READ(UNIT=Seq,Err=100,IOSTAT=IOS) NExpt,NDist,NCoef
    CALL New_HGRho(A,(/NExpt,NDist,NCoef/),REst_O)
!
!   Read In the Density
!
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%NQ%I(I)    ,I=1,A%NExpt)
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%OffQ%I(I)  ,I=1,A%NExpt)
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%OffR%I(I)  ,I=1,A%NExpt)
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Lndx%I(I)  ,I=1,A%NExpt)
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Expt%D(I)  ,I=1,A%NExpt)
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qx%D(I)    ,I=1,A%NDist)
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qy%D(I)    ,I=1,A%NDist)
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qz%D(I)    ,I=1,A%NDist)
    IF(PRESENT(REst_O)) THEN
       IF(REst_O) THEN
          READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Est%D(I),I=1,A%NDist)
       ENDIF
    ENDIF
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Co%D(I) ,I=1,A%NCoef)
!
    Close(UNIT=Seq,STATUS='KEEP')
    RETURN
100 CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Get_HGRho.')
  END SUBROUTINE Get_HGRho
!========================================================================================
! Write  the density to disk 
!========================================================================================
  SUBROUTINE Put_HGRho(A,Name,Args,SCFCycle,REst_O)
    TYPE(HGRho)                      :: A
    TYPE(ARGMT)                      :: Args
    INTEGER                          :: I,SCFCycle,IOS
    REAL(DOUBLE)                     :: Dummy
    CHARACTER(LEN=*)                 :: Name 
    CHARACTER(LEN=DEFAULT_CHR_LEN)   :: FileName
    LOGICAL,OPTIONAL                 :: REst_O
    LOGICAL                          :: Exists
!   
    FileName=TrixFile(Name,Args,SCFCycle)
    INQUIRE(FILE=FileName,EXIST=Exists)
    IF(Exists) THEN
       OPEN(UNIT=Seq,FILE=FileName,STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    ELSE
       OPEN(UNIT=Seq,FILE=FileName,STATUS='NEW',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    ENDIF
!
!   Write density to disk
!
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS) A%NExpt,A%NDist,A%NCoef
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%NQ%I(I)    ,I=1,A%NExpt)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%OffQ%I(I)  ,I=1,A%NExpt)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%OffR%I(I)  ,I=1,A%NExpt)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Lndx%I(I)  ,I=1,A%NExpt)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Expt%D(I)  ,I=1,A%NExpt)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qx%D(I)    ,I=1,A%NDist)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qy%D(I)    ,I=1,A%NDist)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qz%D(I)    ,I=1,A%NDist)
    IF(PRESENT(REst_O)) THEN
       IF(REst_O) THEN
          WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Est%D(I),I=1,A%NDist)
       ENDIF
    ENDIF
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Co%D(I) ,I=1,A%NCoef)
!
    CLOSE(UNIT=Seq,STATUS='KEEP')
    RETURN
100 CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Put_HGRho.')
  END SUBROUTINE Put_HGRho
!------------------------------------------------------------------
!     Open an ASCII file   
! 
      SUBROUTINE OpenASCII(FileName,Unit,NewFile_O,Rewind_O)
         CHARACTER(LEN=*), INTENT(IN) :: FileName
         INTEGER,          INTENT(IN) :: Unit
         LOGICAL, OPTIONAL,INTENT(IN) :: NewFile_O,Rewind_O
         INTEGER                      :: IOS
         LOGICAL                      :: Opened, Exists
!------------------------------------------------------------------
!        Does the file exist?
!
         INQUIRE(FILE=FileName,OPENED=Opened, &
                 EXIST=Exists,ERR=11,IOSTAT=IOS)
!------------------------------------------------------------------
!        Open a new file
!
         IF(PRESENT(NewFile_O))THEN
            IF(NewFile_O.AND.Exists)THEN
!              Open replace if already exists
               OPEN(UNIT=Unit,FILE=FileName, &
                    ACCESS='SEQUENTIAL', FORM='FORMATTED', &
                    ERR=11,IOSTAT=IOS,STATUS='REPLACE')
            ELSE
!              Open a new file
               OPEN(UNIT=Unit,FILE=FileName, &
                    ACCESS='SEQUENTIAL',FORM='FORMATTED', &
                    ERR=11,IOSTAT=IOS,STATUS='NEW')         
            ENDIF
!------------------------------------------------------------------
!        Open existing file and position at the top
!
         ELSEIF(PRESENT(Rewind_O))THEN
            IF(Rewind_O.AND.Exists)THEN
               OPEN(UNIT=Unit,FILE=FileName, &
                    ACCESS='SEQUENTIAL', FORM='FORMATTED', &
                    POSITION='REWIND',ERR=11,IOSTAT=IOS,STATUS='OLD')
            ELSE
               CALL Halt(' Logic error 2 in OpenASCII ')
            ENDIF
!------------------------------------------------------------------
!        Open existing file and position at the bottom (default)
!
         ELSEIF(Exists.AND.(.NOT.Opened))THEN
            OPEN(UNIT=Unit,FILE=FileName, &
                 ACCESS='SEQUENTIAL', FORM='FORMATTED', &
                 POSITION='APPEND',ERR=11,IOSTAT=IOS,STATUS='OLD')
!------------------------------------------------------------------
!        Create a new file and open it
!
         ELSEIF(Exists.AND.Opened)THEN
           CALL Warn(' File '//TRIM(FileName)//' already open')
         ELSE
           OPEN(UNIT=Unit,FILE=FileName, &
                ACCESS='SEQUENTIAL',FORM='FORMATTED', &
                ERR=11,IOSTAT=IOS,STATUS='NEW')         
         ENDIF
         RETURN
      11 CALL Halt(' OpenASCII ERROR: IOS='//TRIM(IntToChar(IOS))// &
                   ' on file '//TRIM(FileName)//'.') 
      END SUBROUTINE OpenASCII

END MODULE

