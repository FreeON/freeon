!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
! Read in and output intermediates for the development
! of new methods for solving the SCF equations
!
PROGRAM ToDev
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE SetXYZ
  USE LinAlg
#ifdef NAG
   USE F90_UNIX_ENV
#endif
  IMPLICIT NONE
  TYPE(ARGMT)                    :: Args
  TYPE(BCSR)                     :: S
  CHARACTER(LEN=5),PARAMETER     :: Prog='ToDev'
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: DevFile
!--------------------------------------------------------------------
!
  CALL StartUp(Args,Prog)
!--------------------------------------------------------------------
  CALL GetEnv('MONDO_HOME',DevFile)
  DevFile=TRIM(DevFile)//'/SCFeqs/MMA/SCFDev/scf.in'
  CALL OpenASCII(DevFile,35,NewFile_O=.TRUE.)
  PrintFlags%Fmt=DEBUG_MMASTYLE
  PrintFlags%Mat=DEBUG_MATRICES
  WRITE(35,*)' (* All matrices are in an orthogonal rep *)'
  WRITE(35,*)' NBasF  = ',NBasF,';'
  WRITE(35,*)' NEl    = ',NEl,';'
  CALL New(S)
  CALL Get(s,TrixFile('OrthoF',Args,0))
  CALL PPrint(s,'f',FileName_O=DevFile,Unit_O=35)   
  CALL get(s,TrixFile('OrthoD',Args,1))
  CALL PPrint(s,'pp',FileName_O=DevFile,Unit_O=35)   
  CLOSE(UNIT=35,STATUS='KEEP')
!----------------------------------------------
  CALL Delete(s)
  CALL ShutDown(Prog)
END PROGRAM 


