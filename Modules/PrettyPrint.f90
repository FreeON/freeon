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
!    GENERIC PRETTY PRINTING FOR MONDOSCF TYPES
!    Author: Matt Challacombe
!---------------------------------------------------------

#include "MondoConfig.h"

MODULE PrettyPrint
   USE DerivedTypes
   USE GlobalCharacters
   USE GlobalScalars
   USE GlobalObjects
   USE MondoLogger
   USE BasisSetParameters
   USE Clock
   USE Parse
   USE SetXYZ
   USE InOut

   USE ParsingConstants
!   USE Order
#ifdef PARALLEL
   USE MondoMPI
#endif

   IMPLICIT NONE

   INTERFACE PPrint
      MODULE PROCEDURE Print_INT_SCLR,   Print_DBL_SCLR,  &
                       Print_CHR_SCLR,   Print_INT_VECT,  &
                       Print_DBL_VECT,   Print_DBL_RNK2,  &
                       Print_DBL_Rank2A, Print_BSET,      &
                       Print_CRDS,       Print_BCSR,      &
#ifdef PARALLEL
                       Print_DBCSR,                       &
#endif
                       Print_PBCInfo,    Print_CellSet,   &
                       Print_MEMS,       Print_TIME,      &
                       Print_HGRho
   END INTERFACE

   INTERFACE PChkSum
      MODULE PROCEDURE Print_CheckSum_DBL_VECT
      MODULE PROCEDURE Print_CheckSum_DBL_RNK2
      MODULE PROCEDURE Print_CheckSum_BCSR
      MODULE PROCEDURE Print_CheckSum_HGRho
#ifdef PARALLEL
      MODULE PROCEDURE Print_CheckSum_DBCSR
#endif
   END INTERFACE

   INTERFACE CheckSum
     MODULE PROCEDURE CheckSum_BCSR
   END INTERFACE CheckSum

   CHARACTER(LEN=DEFAULT_CHR_LEN) :: String
   CHARACTER(LEN=*),PARAMETER     :: CheckEq = ' CheckSum = '
   CONTAINS
      SUBROUTINE TimeStamp(Mssg,Enter_O)
         CHARACTER(LEN=*),INTENT(IN) :: Mssg
         LOGICAL,OPTIONAL,INTENT(IN) :: Enter_O
         CHARACTER(LEN=8)            :: DDate,DDay,HMS
         CHARACTER(LEN=10)           :: TTime
         CHARACTER(LEN=5)            :: Zone
         INTEGER, DIMENSION(8)       :: Values
         LOGICAL                     :: Enter
#ifdef PARALLEL
         IF(MyID==0)THEN
#endif
         Enter=.TRUE.; IF(PRESENT(Enter_O))Enter=Enter_O
         CALL DATE_AND_TIME(DDate,TTime,Zone,Values)
         DDay=DDate(5:6)//'/'//DDate(7:8)//'/'//DDate(3:4)
         HMS=TTime(1:2)//':'//TTime(3:4)//':'//TTime(5:6)
         IF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
            CALL MondoLog(DEBUG_NONE, "TimeStamp", '(*'//TRIM(Mssg)//' '//TRIM(DDay)//' @ '//TRIM(HMS)//'*)')
         ELSEIF(Enter)THEN
            CALL MondoLog(DEBUG_NONE, "TimeStamp", '<<'//TRIM(Mssg)//' '//TRIM(DDay)//' @ '//TRIM(HMS)//'>>')
         ELSE
            CALL MondoLog(DEBUG_NONE, "TimeStamp", '>-'//TRIM(Mssg)//' '//TRIM(DDay)//' @ '//TRIM(HMS)//'-<')
         ENDIF
#ifdef PARALLEL
         ENDIF
#endif
      END SUBROUTINE TimeStamp

      FUNCTION OpenPU(FileName_O,Unit_O,NewFile_O)
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: FileName_O
         INTEGER,         OPTIONAL,INTENT(IN) :: Unit_O
         LOGICAL,         OPTIONAL,INTENT(IN) :: NewFile_O
         INTEGER                              :: OpenPU,PU
         IF(PRESENT(FileName_O).AND.PRESENT(Unit_O))THEN
            PU=Unit_O
            CALL OpenASCII(FileName_O,PU,NewFile_O=NewFile_O)
         ELSEIF(PRESENT(FileName_O).AND.(.NOT.PRESENT(Unit_O)))THEN
            PU=Tmp
            CALL OpenASCII(FileName_O,PU,NewFile_O=NewFile_O)
         ELSEIF(PRESENT(Unit_O))THEN
            PU=Unit_O
         ELSE
            PU=Out
            CALL OpenASCII(OutFile,PU,NewFile_O=NewFile_O)
         ENDIF
         OpenPU=PU
      END FUNCTION OpenPU

      SUBROUTINE ClosePU(U)
         INTEGER,INTENT(IN) :: U
         IF(U/=6)CLOSE(UNIT=U,STATUS='KEEP')
      END SUBROUTINE ClosePU

      SUBROUTINE Print_CHR_SCLR(X,Name_O,FileName_O,Unit_O)
         CHARACTER(LEN=*),INTENT(IN)          :: X
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Name_O,FileName_O
         INTEGER, OPTIONAL,INTENT(IN)         :: Unit_O
         INTEGER                              :: PU
         PU=OpenPU(FileName_O,Unit_O)
         IF(PrintFlags%Fmt==DEBUG_MMASTYLE.AND.PRESENT(Name_O))THEN
            WRITE(PU,11)TRIM(Name_O),TRIM(X)
         ELSEIF(PRESENT(Name_O))THEN
            WRITE(PU,12)TRIM(Name_O),TRIM(X)
         ELSEIF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
            WRITE(PU,13)TRIM(X)
         ELSE
            WRITE(PU,14)TRIM(X)
         ENDIF
         CALL ClosePU(PU)
      11 FORMAT(1x,A,' = ',A,';')
      12 FORMAT(1x,A,' = ',A)
      13 FORMAT('(* ',A,' *)')
      14 FORMAT(1x,A)
      END SUBROUTINE Print_CHR_SCLR
!----------------------------------------------------------------PRINT AN INTEGER
      SUBROUTINE Print_INT_SCLR(X,Name,FileName_O,Unit_O,Protect_O)
         INTEGER,                   INTENT(IN) :: X
         CHARACTER(LEN=*),          INTENT(IN) :: Name
         CHARACTER(LEN=*), OPTIONAL,INTENT(IN) :: FileName_O
         INTEGER,          OPTIONAL,INTENT(IN) :: Unit_O
         LOGICAL,          OPTIONAL,INTENT(IN) :: Protect_O
         INTEGER                               :: PU
         CHARACTER(LEN=INTERNAL_INT_LEN)       :: CTmp,Id
         INTEGER                               :: I,J,M,N
         LOGICAL                               :: Protect
!--------------------------------------------------------
         CTmp=IntToChar(X)
         IF(PRESENT(Protect_O))THEN
            Protect=Protect_O
         ELSE
            Protect=.TRUE.
         ENDIF
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            PU=OpenPU(FileName_O,Unit_O)
            IF(Protect) &
               CALL PrintProtectL(PU)
            CALL ClosePU(PU)
#ifdef PARALLEL
         ENDIF
         IF(InParallel)THEN
            DO I=0,NPrc-1
               CALL AlignNodes()
               IF(MyId==I)THEN
                  PU=OpenPU(FileName_O,Unit_O)
                  Id=IntToChar(MyId)
                  IF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
                     WRITE(PU,*)TRIM(Name),'[',TRIM(Id),'] = ',TRIM(CTmp),';'
                  ELSE
                     WRITE(PU,*)TRIM(Name),'[',TRIM(Id),'] := ',TRIM(CTmp)
                  ENDIF
                  CALL ClosePU(PU)
               ENDIF
            ENDDO
         ELSE
#endif
            PU=OpenPU(FileName_O,Unit_O)
            IF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
               WRITE(PU,*)TRIM(Name),' = ',TRIM(CTmp),';'
            ELSE
               WRITE(PU,*)TRIM(Name),' := ',TRIM(CTmp)
            ENDIF
            CALL ClosePU(PU)
#ifdef PARALLEL
         ENDIF
         IF(MyId==ROOT)THEN
#endif
            PU=OpenPU(FileName_O,Unit_O)
            IF(Protect) &
               CALL PrintProtectR(PU)
            CALL ClosePU(PU)
#ifdef PARALLEL
         ENDIF
#endif
      END SUBROUTINE Print_INT_SCLR
!----------------------------------------------------------------
!     PRINT A DOUBLE
!
      SUBROUTINE Print_DBL_SCLR(X,Name,FileName_O,Unit_O)
         REAL(DOUBLE),              INTENT(IN) :: X
         CHARACTER(LEN=*),          INTENT(IN) :: Name
         CHARACTER(LEN=*), OPTIONAL,INTENT(IN) :: FileName_O
         INTEGER,          OPTIONAL,INTENT(IN) :: Unit_O
         INTEGER                               :: PU
         PU=OpenPU(FileName_O,Unit_O)
         IF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
            WRITE(PU,11)Name,FRACTION(X),EXPONENT(X)
         ELSE
            String=TRIM(Name)//' = '//TRIM(DblToChar(X))
            WRITE(PU,"(A)")TRIM(String)
         ENDIF
         CALL ClosePU(PU)
      11 FORMAT(A,' = ',F19.16,'*2^(',I4,');')
      END SUBROUTINE Print_DBL_SCLR
!----------------------------------------------------------------
!     PRINT AN INT_VECT
!
      SUBROUTINE Print_INT_VECT(A,Name,N_O,M_O,FileName_O,Unit_O)
         TYPE(INT_VECT),            INTENT(IN) :: A
         CHARACTER(LEN=*),          INTENT(IN) :: Name
         CHARACTER(LEN=*), OPTIONAL,INTENT(IN) :: FileName_O
         INTEGER,          OPTIONAL,INTENT(IN) :: Unit_O,M_O,N_O
         TYPE(CHR_VECT)                        :: CA
         CHARACTER(LEN=INTERNAL_INT_LEN)       :: Id
         INTEGER                               :: PU
         INTEGER                               :: I,J,M,N
         M=1;         IF(PRESENT(M_O))M=M_O
         N=SIZE(A%I); IF(PRESENT(N_O))N=N_O
         CALL New(CA,N,M)
         DO J=M,N; CA%C(J)=IntToChar(A%I(J)); ENDDO
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            PU=OpenPU(FileName_O,Unit_O)
            IF(PU/=6)THEN
               CALL PrintProtectL(PU)
               CLOSE(PU)
            ENDIF
#ifdef PARALLEL
         ENDIF
         DO I=0,NPrc-1
            CALL AlignNodes()
            IF(MyId==I)THEN
               PU=OpenPU(FileName_O,Unit_O)
               Id=IntToChar(MyId)
               WRITE(PU,*)TRIM(Name),'[',TRIM(Id),'] := ',(TRIM(CA%C(J)),', ',J=M,N)
               CALL ClosePU(PU)
            ENDIF
         ENDDO
#else
         PU=OpenPU(FileName_O,Unit_O)
         WRITE(PU,*)TRIM(Name),' := ',(TRIM(CA%C(J)),', ',J=1,N)
         CALL ClosePU(PU)
#endif
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            PU=OpenPU(FileName_O,Unit_O)
            IF(PU/=6)THEN
               CALL PrintProtectR(PU)
               CLOSE(PU)
            ENDIF
#ifdef PARALLEL
         ENDIF
#endif
         CALL Delete(CA)
      END SUBROUTINE Print_INT_VECT
!----------------------------------------------------------------
!     PRINT AN INT_VECT
!
      SUBROUTINE Print_DBL_VECT(A,Name,N_O,M_O,FileName_O,Unit_O)
         TYPE(DBL_VECT),            INTENT(IN) :: A
         CHARACTER(LEN=*),          INTENT(IN) :: Name
         CHARACTER(LEN=*), OPTIONAL,INTENT(IN) :: FileName_O
         CHARACTER(LEN=DEFAULT_CHR_LEN)        :: outputString
         INTEGER,          OPTIONAL,INTENT(IN) :: Unit_O,M_O,N_O
         TYPE(CHR_VECT)                        :: CA
         CHARACTER(LEN=INTERNAL_INT_LEN)       :: Id
         INTEGER                               :: PU
         INTEGER                               :: I,J,M,N
         M=1;         IF(PRESENT(M_O))M=M_O
         N=SIZE(A%D); IF(PRESENT(N_O))N=N_O
         CALL New(CA,N,M)
         DO J=M,N; CA%C(J)=DblToMedmChar(A%D(J)); ENDDO
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            PU=OpenPU(FileName_O,Unit_O)
            IF(PU/=6)THEN
               CALL PrintProtectL(PU)
               CLOSE(PU)
            ENDIF
#ifdef PARALLEL
         ENDIF
         DO I=0,NPrc-1
            CALL AlignNodes()
            IF(MyId==I)THEN
               PU=OpenPU(FileName_O,Unit_O)
               Id=IntToChar(MyId)
               WRITE(PU,*)TRIM(Name),'[',TRIM(Id),'] := ',(TRIM(CA%C(J)),', ',J=M,N)
               CALL ClosePU(PU)
            ENDIF
         ENDDO
#else
         WRITE(outputString, *) TRIM(Name), ' := ', (TRIM(CA%C(J)), ', ', J=1, N)
         CALL MondoLog(DEBUG_NONE, "", TRIM(ADJUSTL(outputString)))
#endif
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            PU=OpenPU(FileName_O,Unit_O)
            IF(PU/=6)THEN
               CALL PrintProtectR(PU)
               CLOSE(PU)
            ENDIF
#ifdef PARALLEL
         ENDIF
#endif
         CALL Delete(CA)
      END SUBROUTINE Print_DBL_VECT

!----------------------------------------------------------------PRINT BASIS SET
      SUBROUTINE Print_BSET(BS,Unit_O)
        TYPE(BSET)                  :: BS
        INTEGER,OPTIONAL,INTENT(IN) :: Unit_O
        INTEGER                     :: NC,NP,MinL,MaxL
        INTEGER                     :: I,J,K,L,M,PU

        !IF(PrintFlags%Set/=DEBUG_BASISSET)RETURN
        PU=OpenPU(Unit_O=Unit_O)
        CALL PrintProtectL(PU)
        WRITE(PU,*)'Internal representation of the ', &
             TRIM(BS%BName),' basis set: '//Rtrn
        DO I=1,BS%NKind
           WRITE(PU,1001)
           NC=BS%NCFnc%I(I)
           WRITE(PU,1002)Ats(BS%Kinds%I(I)),NC
           DO J=1,NC
              NP=BS%NPFnc%I(J,I)
              MinL=LBegin(BS%ASymm%I(1,J,I))
              MaxL=LEnd(BS%ASymm%I(2,J,I))
              IF(NP==1)THEN
                 WRITE(PU,1103)J,NP
                 IF(MaxL==1)THEN
                    WRITE(PU,1104)
                 ELSE
                    WRITE(PU,1004)
                 ENDIF
              ELSE
                 WRITE(PU,1003)J,NP
                 WRITE(PU,1004)
              ENDIF
              WRITE(PU,1005)(ASymmTyps(L),L=MinL,MaxL)
              DO K=1,NP
                 WRITE(PU,1006)K,BS%Expnt%D(K,J,I), &
                      (BS%CCoef%D(M,K,J,I),M=MinL,MaxL)
              ENDDO
!              WRITE(PU,*)Rtrn
           ENDDO
           IF(BS%HasECPs)THEN
              IF(BS%NTyp1PF%I(I)>0)THEN
                 K=BS%NCoreEl%D(I)
                 WRITE(PU,2000)Ats(BS%Kinds%I(I)),K
                 WRITE(PU,2001)
                 WRITE(PU,2002)
                 DO J=1,BS%NTyp1PF%I(I)
                    WRITE(PU,2006)BS%Typ1Ell%I(J,I),BS%Typ1Exp%D(J,I),BS%Typ1CCo%D(J,I)
                 ENDDO
                 DO K=0,BS%ProjEll%I(I)
                    WRITE(PU,2003)K
                    WRITE(PU,2002)
                    DO J=1,BS%NTyp2PF%I(K,I)
                       WRITE(PU,2006)BS%Typ2Ell%I(J,K,I),BS%Typ2Exp%D(J,K,I),BS%Typ2CCo%D(J,K,I)
                    ENDDO
                 ENDDO
              ENDIF
           ENDIF
        ENDDO
        WRITE(PU,1001)
        CALL PrintProtectR(PU)
        CALL ClosePU(PU)

2000    FORMAT(1x,A2,' has an ECP replacing ',I2,' electrons')
2001    FORMAT(1x,'Type 1 ECP contraction')
2002    FORMAT(1x,' Ell      Exponent       CCoefficient ')
2003    FORMAT(1x,'Type 2 ECP contraction, projector symmetry:',I2)
2006    FORMAT(1x,I2,6x,2(1x,D14.8))

1001    FORMAT(72('='))
1002    FORMAT(1x,A2,' has ',I2,' associated contractions,')
1003    FORMAT(1x,'Contraction ',I2,' involves', &
             I2,' primitives : ')
1004    FORMAT(1x,'Primitive  Exponent',7x,'Normalized Coeficients')
1103    FORMAT(1x,'Contraction ',I2,' involves', &
             I2,' primitive : ')
1104    FORMAT(1x,'Primitive  Exponent',7x,'Normalized Coeficient')

1005    FORMAT(18x,20(12x,A3))
1006    FORMAT(1x,I2,6x,8(1x,D14.8))
1007    FORMAT(60('='))
      END SUBROUTINE Print_BSET
!----------------------------------------------------------------PRINT PBC
!
!----------------------------------------------------------------PRINT PBC
      SUBROUTINE Print_PBCInfo(PBC,FileName_O,Unit_O)
        TYPE(PBCInfo)                        :: PBC
        INTEGER                              :: I,K,PU
        LOGICAL                              :: Opened
        INTEGER,         OPTIONAL,INTENT(IN) :: Unit_O
        CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: FileName_O
        CHARACTER(LEN=DEFAULT_CHR_LEN)       :: Mssg
!
        PU=OpenPU(FileName_O=FileName_O,Unit_O=Unit_O)

        Mssg='PFFMaxEll = '//TRIM(IntToChar(PBC%PFFMaxEll))
        WRITE(PU,*)TRIM(Mssg)
        Mssg='PFFWellSep = '//TRIM(IntToChar(PBC%PFFWelSep))
        WRITE(PU,*)TRIM(Mssg)

        Mssg='PBCs in '//TRIM(IntToChar(PBC%Dimen))//' dimensions'
        WRITE(PU,*)TRIM(Mssg)

        Mssg="Cell Volume = "//TRIM(DblToMedmChar(PBC%CellVolume))
        WRITE(PU,*)TRIM(Mssg)

        Mssg='PBC         = ('//TRIM(IntToChar(PBC%AutoW%I(1)))//',' &
                              //TRIM(IntToChar(PBC%AutoW%I(2)))//',' &
                              //TRIM(IntToChar(PBC%AutoW%I(3)))//')'
        WRITE(PU,*)TRIM(Mssg)

        Mssg=' TransVec   = ('//TRIM(DblToMedmChar(PBC%TransVec%D(1)))//', ' &
                              //TRIM(DblToMedmChar(PBC%TransVec%D(2)))//', ' &
                              //TRIM(DblToMedmChar(PBC%TransVec%D(3)))//') '
        WRITE(PU,*)TRIM(Mssg)
        Mssg=' CellCenter = ('//TRIM(DblToMedmChar(PBC%CellCenter%D(1)))//', ' &
                              //TRIM(DblToMedmChar(PBC%CellCenter%D(2)))//', ' &
                              //TRIM(DblToMedmChar(PBC%CellCenter%D(3)))//') '
        WRITE(PU,*)TRIM(Mssg)
        WRITE(PU,*)' Lattice Vectors: '
        Mssg=' a =  ('//TRIM(DblToMedmChar(PBC%BoxShape%D(1,1)))//', ' &
                      //TRIM(DblToMedmChar(PBC%BoxShape%D(2,1)))//', ' &
                      //TRIM(DblToMedmChar(PBC%BoxShape%D(3,1)))//') '
        WRITE(PU,*)TRIM(Mssg)
        Mssg=' b =  ('//TRIM(DblToMedmChar(PBC%BoxShape%D(1,2)))//', ' &
                      //TRIM(DblToMedmChar(PBC%BoxShape%D(2,2)))//', ' &
                      //TRIM(DblToMedmChar(PBC%BoxShape%D(3,2)))//') '
        WRITE(PU,*)TRIM(Mssg)
        Mssg=' c =  ('//TRIM(DblToMedmChar(PBC%BoxShape%D(1,3)))//', ' &
                      //TRIM(DblToMedmChar(PBC%BoxShape%D(2,3)))//', ' &
                      //TRIM(DblToMedmChar(PBC%BoxShape%D(3,3)))//') '
        WRITE(PU,*)TRIM(Mssg)
        WRITE(PU,*)' Inverse Lattice Vectors: '
        Mssg=' 1/a = ('//TRIM(DblToMedmChar(PBC%InvBoxSh%D(1,1)))//', ' &
                       //TRIM(DblToMedmChar(PBC%InvBoxSh%D(2,1)))//', ' &
                       //TRIM(DblToMedmChar(PBC%InvBoxSh%D(3,1)))//') '
        WRITE(PU,*)TRIM(Mssg)
        Mssg=' 1/b = ('//TRIM(DblToMedmChar(PBC%InvBoxSh%D(1,2)))//', ' &
                       //TRIM(DblToMedmChar(PBC%InvBoxSh%D(2,2)))//', ' &
                       //TRIM(DblToMedmChar(PBC%InvBoxSh%D(3,2)))//') '
        WRITE(PU,*)TRIM(Mssg)
        Mssg=' 1/c = ('//TRIM(DblToMedmChar(PBC%InvBoxSh%D(1,3)))//', ' &
                       //TRIM(DblToMedmChar(PBC%InvBoxSh%D(2,3)))//', ' &
                       //TRIM(DblToMedmChar(PBC%InvBoxSh%D(3,3)))//') '
        WRITE(PU,*)TRIM(Mssg)
        CALL ClosePU(PU)
!
     END SUBROUTINE Print_PBCInfo
!--------------------------------------------------------------------------
! Print the CellSet
!--------------------------------------------------------------------------
  SUBROUTINE Print_CellSet(CS,Name,Proc_O,Unit_O)
     TYPE(CellSet)                    :: CS
     CHARACTER(LEN=*),OPTIONAL        :: Proc_O
     CHARACTER(LEN=*)                 :: Name
     INTEGER,OPTIONAL                 :: Unit_O
     INTEGER                          :: PU
     CHARACTER(LEN=2*DEFAULT_CHR_LEN) :: CellStr
!-----------------------------------------------------------------------------------------
     IF(PrintFlags%Key<DEBUG_MAXIMUM.OR.CS%NCells==1)RETURN
     IF(.NOT. AllocQ(CS%Alloc))CALL Halt(' Cells not allocated in Print_CellSet')
#ifdef PARALLEL
   IF(MyID == ROOT) THEN
#endif
     PU=OpenPU(Unit_O=Unit_O)
     IF(PRESENT(Proc_O)) THEN
        CellStr=ProcessName(Proc_O)//'Cells in '//TRIM(Name)//' = '//TRIM(IntToChar(CS%NCells))
     ELSE
        CellStr='Cells in '//TRIM(Name)//' = '//TRIM(IntToChar(CS%NCells))
     ENDIF
     WRITE(PU,*)TRIM(CellStr)
     CALL ClosePU(PU)
#ifdef PARALLEL
   ENDIF
#endif
  END SUBROUTINE Print_CellSet
!----------------------------------------------------------------PRINT COORDINATES
!
!----------------------------------------------------------------PRINT COORDINATES

     SUBROUTINE XSFPreamble(Steps,FileName,Unit)
       INTEGER                     :: Steps
       INTEGER,         INTENT(IN) :: Unit
       CHARACTER(LEN=*),INTENT(IN) :: FileName
       CHARACTER(LEN=DEFAULT_CHR_LEN)       :: Mssg
       CALL OpenASCII(FileName,Unit,Rewind_O=.TRUE.)
       Mssg='ANIMSTEPS  '//IntToChar(Steps)
       WRITE(Unit,*)TRIM(Mssg)
       CLOSE(Unit)
     END SUBROUTINE XSFPreamble

     SUBROUTINE Print_CRDS(GM, FileName_O, Unit_O, PrintGeom_O, NewFile_O, Clone_O, CrdInAng_O, Remark_O, Gradients_O)
       TYPE(CRDS)                           :: GM
       INTEGER                              :: K
       LOGICAL                              :: Opened
       INTEGER,         OPTIONAL,INTENT(IN) :: Unit_O,Clone_O
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: FileName_O,PrintGeom_O,Remark_O,Gradients_O
       LOGICAL,         OPTIONAL,INTENT(IN) :: NewFile_O,CrdInAng_O
       LOGICAL                              :: CrdInAng
       INTEGER                              :: I,PU
       CHARACTER(LEN=2*DCL)                 :: Mssg
       CHARACTER(LEN=DCL)                   :: AuxChar
       CHARACTER(LEN=10)                    :: Gradients
       CHARACTER(LEN=2)                     :: Atom
       CHARACTER(LEN=4)                     :: PDBAtom
       REAL(DOUBLE)                         :: AA,GAA
       REAL(DOUBLE), DIMENSION(3)           :: VTmp
       REAL(DOUBLE)                         :: A,B,C,Alpha,Beta,Gamma

       ! Are we printing out gradients or forces?
       IF(PRESENT(Gradients_O))THEN
          Gradients=TRIM(Gradients_O)
       ELSE
          Gradients=""
       ENDIF
       !Check if the coordinates are in Angstrom.
       IF(PRESENT(CrdInAng_O))THEN
          CrdInAng=CrdInAng_O
       ELSE
          CrdInAng=.FALSE.
       ENDIF
       IF(CrdInAng)THEN
          AA=One
          GAA=One
       ELSEIF(TRIM(Gradients)=='Velocities')THEN
          AA=AUToAngstroms
          GAA=AUToAngstroms/InternalTimeToFemtoseconds
       ELSEIF(TRIM(Gradients)=='Gradients')THEN
          AA=AUToAngstroms
          GAA=au2eV/AUToAngstroms
       ELSE
          AA=AUToAngstroms
       ENDIF
       !
#ifdef PARALLEL
       IF(MyId==ROOT)THEN
#endif
          PU=OpenPU(FileName_O=FileName_O,Unit_O=Unit_O,NewFile_O=NewFile_O)
          IF(PRESENT(PrintGeom_O))THEN
             IF(PrintGeom_O=='XYZ'.OR.PrintGeom_O=='XSF') THEN
!============================================================================
                IF(PrintGeom_O=='XYZ')THEN
!============================================================================
                   ! Print XYZ format
                   Mssg=IntToChar(GM%NAtms)
                   WRITE(PU,*)TRIM(Mssg)
                   IF(PRESENT(Remark_O))THEN
                      Mssg=TRIM(Remark_O)//', <SCF> = '//TRIM(DblToChar(GM%ETotal))//' Hartree, '//TRIM(DblToChar(GM%ETotal*au2eV))//' eV'
                   ELSE
                      Mssg='Geom #'//TRIM(IntToChar(GM%Confg)) &
                           //', <SCF> = '//TRIM(DblToChar(GM%ETotal))//' Hartree, '//TRIM(DblToChar(GM%ETotal*au2eV))//' eV'
                   ENDIF
                   IF(PRESENT(Clone_O)) &
                        Mssg='Clone # '//TRIM(IntToChar(Clone_O))//" / "//TRIM(Mssg)
                   IF(GM%PBC%Dimen/=0)THEN
                      CALL VecToAng(GM%PBC,a,b,c,alpha,beta,gamma)
                      Mssg=TRIM(Mssg)//", PBC = " &
                           //TRIM(FltToMedmChar(a*AA))//" " &
                           //TRIM(FltToMedmChar(b*AA))//" " &
                           //TRIM(FltToMedmChar(c*AA))//" " &
                           //TRIM(FltToMedmChar(alpha))//" " &
                           //TRIM(FltToMedmChar(beta))//" " &
                           //TRIM(FltToMedmChar(gamma))
                   ENDIF
                   WRITE(PU,"(A)") TRIM(Mssg)
!============================================================================
                ELSE  !! XSF
!============================================================================
                   IF(GM%PBC%Dimen==0)THEN
                      Mssg='MOLECULE'
                   ELSEIF(GM%PBC%Dimen==1)THEN
                      Mssg='WIRE'
                   ELSEIF(GM%PBC%Dimen==2)THEN
                      Mssg='SLAB'
                   ELSE
                      Mssg='CRYSTAL'
                   ENDIF
                   WRITE(PU,*)TRIM(Mssg)
                   Mssg='PRIMVEC '//TRIM(IntToChar(GM%Confg+1))
                   WRITE(PU,*)TRIM(Mssg)
                   ! Cell Vectors
                   GM%PBC%BoxShape%D=GM%PBC%BoxShape%D*AA
                   WRITE(PU,111)GM%PBC%BoxShape%D(1,:)
                   WRITE(PU,111)GM%PBC%BoxShape%D(2,:)
                   WRITE(PU,111)GM%PBC%BoxShape%D(3,:)
111                FORMAT(1X,3(F14.5,' '))
                   GM%PBC%BoxShape%D=GM%PBC%BoxShape%D/AA
                   Mssg='PRIMCOORD '//TRIM(IntToChar(GM%Confg+1))
                   WRITE(PU,*)TRIM(Mssg)
                   WRITE(PU,*)GM%NAtms,1
                ENDIF
                ! XYZ style formating
                DO I=1,GM%NAtms
                   IF(GM%CConstrain%I(I)==1) THEN
                      AuxChar='C'
                   ELSE IF(GM%CConstrain%I(I)==2) THEN
                      AuxChar='R'
                   ELSE
                      AuxChar=' '
                   ENDIF
                   Atom=GM%AtNam%C(I)
                   CALL UpCase(Atom)
                   IF(TRIM(Gradients)=="Gradients")THEN
                      WRITE(PU,222)Atom,(GM%Carts%D(K,I)*AA,K=1,3),(-GM%Gradients%D(K,I)*GAA,K=1,3),TRIM(ADJUSTL(AuxChar))
                   ELSEIF(TRIM(Gradients)=="Velocities")THEN
                      WRITE(PU,222)Atom,(GM%Carts%D(K,I)*AA,K=1,3),(GM%Velocity%D(K,I)*GAA,K=1,3),TRIM(ADJUSTL(AuxChar))
                   ELSE
                      WRITE(PU,223)Atom,(GM%Carts%D(K,I)*AA,K=1,3),TRIM(ADJUSTL(AuxChar))
                   ENDIF
                ENDDO
222             FORMAT(1X,A2,6(F22.16,' '),A1)
223             FORMAT(1X,A2,3(F22.16,' '),A1)

!!$                DO I=1,GM%NAtms
!!$                   IF(GM%CConstrain%I(I)==1) THEN
!!$                      AuxChar='C'
!!$                   ELSE IF(GM%CConstrain%I(I)==2) THEN
!!$                      AuxChar='R'
!!$                   ELSE
!!$                      AuxChar=' '
!!$                   ENDIF
!!$                   Atom=GM%AtNam%C(I)
!!$                   CALL UpCase(Atom)
!!$                   WRITE(PU,223)Atom,(GM%BoxCarts%D(K,I),K=1,3),TRIM(ADJUSTL(AuxChar))
!!$223                FORMAT(1X,A2,3(F22.16,' '),A1,' << FRACTIONALS ')
!!$                ENDDO
!============================================================================
             ELSEIF(PrintGeom_O=='PDB')THEN
!============================================================================
                AA=One/AngstromsToAU
                WRITE(PU,11)GM%Confg
11              FORMAT('MODEL  ',I6)
                IF(PRESENT(Remark_O))THEN
                   Mssg=TRIM(Remark_O)//', <SCF> = '//TRIM(DblToChar(GM%ETotal))//' Hartree, '//TRIM(DblToChar(GM%ETotal*au2eV))//' eV'
                ELSE
                   Mssg=' <SCF> = '//TRIM(DblToChar(GM%ETotal))//' Hartree, '//TRIM(DblToChar(GM%ETotal*au2eV))//' eV'
                ENDIF
                WRITE(PU,22) TRIM(Mssg)
22              FORMAT('REMARK   1  ', A)
                IF(GM%PBC%Dimen/=0)THEN
                   CALL VecToAng(GM%PBC,a,b,c,alpha,beta,gamma)
                   WRITE(PU,33)a*AA,b*AA,c*AA,alpha,beta,gamma
33                 FORMAT('CRYST1',3F9.3,3F7.2,1x,11A1,I4)
                ENDIF
                DO I=1,GM%NAtms
                   Atom=GM%AtNam%C(I)
                   IF(GM%CConstrain%I(I)==1) THEN
                      AuxChar='C'
                   ELSE IF(GM%CConstrain%I(I)==2) THEN
                      AuxChar='R'
                   ELSE
                      AuxChar=' '
                   ENDIF
                   CALL UpCase(Atom)
                   IF(TRIM(Gradients)=='Gradients')THEN
                      WRITE(PU,44)'ATOM  ',  I,     Atom//'  ', ' ',  'UNK',   ' ', 1, ' ' ,   (GM%Carts%D(K,I)*AA,K=1,3), &
                                                                                              (-GM%Gradients%D(K,I)*GAA,K=1,3),AuxChar
                   ELSEIF(TRIM(Gradients)=='Velocities')THEN
                      WRITE(PU,44)'ATOM  ',  I,     Atom//'  ', ' ',  'UNK',   ' ', 1, ' ' ,   (GM%Carts%D(K,I)*AA,K=1,3), &
                                                                                               (GM%Velocity%D(K,I)*GAA,K=1,3),AuxChar
                   ELSE
                      WRITE(PU,45)'ATOM  ',  I,     Atom//'  ', ' ',  'UNK',   ' ', 1, ' ' ,   (GM%Carts%D(K,I)*AA,K=1,3),One,Zero,AuxChar
                   ENDIF
                ENDDO
                WRITE(PU,55)
44              FORMAT(     A6,      I5, 1X,A4,    A1,   A3,  1X, A1 ,I4,  A1,3X,6F8.3,1x,A2)
45              FORMAT(     A6,      I5, 1X,A4,    A1,   A3,  1X, A1 ,I4,  A1,3X,3F8.3,2F6.2,10X,A2,A2 )
55              FORMAT('ENDMDL')
!============================================================================
             ELSEIF(PrintGeom_O=='CIF')THEN
!============================================================================
                ! See http://www.iucr.ac.uk/iucr-top/cif/cifdic_html/1/cif_core.dic/index.html
                ! for the definition of the entries.
                !
                CALL VecToAng(GM%PBC,A,B,C,Alpha,Beta,Gamma)
                !
                IF(PRESENT(Remark_O))THEN
                   Mssg='data_'//TRIM(Remark_O)//'_<SCF>_=_'//TRIM(DblToChar(GM%ETotal*au2eV))//' eV'
                ELSE
                   Mssg='data_'//'Geom#'//TRIM(IntToChar(GM%Confg))//'_<SCF>_=_'//TRIM(DblToChar(GM%ETotal*au2eV))//'_eV'
                ENDIF
                Mssg=Squish(Mssg)
                WRITE(PU,500) Mssg
                !Cell
                WRITE(PU,501)
                WRITE(PU,502) A
                WRITE(PU,503) B
                WRITE(PU,504) C
                WRITE(PU,505) Alpha
                WRITE(PU,506) Beta
                WRITE(PU,507) Gamma
                WRITE(PU,508) GM%PBC%CellVolume
                !Atoms
                WRITE(PU,550)
                WRITE(PU,551)
                WRITE(PU,552)
                WRITE(PU,553)
                WRITE(PU,554)
                DO I=1,GM%NAtms
                   Atom=GM%AtNam%C(I)
                   Vtmp(1:3)=MATMUL(GM%PBC%InvBoxSh%D,GM%Carts%D(1:3,I))
                   CALL UpCase(Atom)
                   WRITE(PU,555)Atom,(VTmp(K),K=1,3)
                ENDDO
                !
                !Cell
500             FORMAT(A60)
501             FORMAT('_space_group_symop_operation_xyz ''x, y, z''')
502             FORMAT('_cell_length_a    ',F14.5)
503             FORMAT('_cell_length_b    ',F14.5)
504             FORMAT('_cell_length_c    ',F14.5)
505             FORMAT('_cell_angle_alpha ',F14.5)
506             FORMAT('_cell_angle_beta  ',F14.5)
507             FORMAT('_cell_angle_gamma ',F14.5)
508             FORMAT('_cell_volume      ',F14.5)
                !Atoms
550             FORMAT('loop_')
551             FORMAT('_atom_site_type_symbol')
552             FORMAT('_atom_site_fract_x')
553             FORMAT('_atom_site_fract_y')
554             FORMAT('_atom_site_fract_z')
555             FORMAT(1X,A,1X,3F14.5)
             ENDIF
!============================================================================
          ELSEIF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
!============================================================================
             Mssg=' NAtoms = '//IntToChar(GM%NAtms)//';'
             WRITE(PU,*)TRIM(Mssg)
             WRITE(PU,*)'(* Coordinates are in AU *)'
             DO I=1,GM%NAtms
                Mssg='R['//TRIM(IntToChar(I))//']={'           &
                     //DblToMMAChar(GM%Carts%D(1,I))//','  &
                     //DblToMMAChar(GM%Carts%D(2,I))//','  &
                     //DblToMMAChar(GM%Carts%D(3,I))//'};'
                WRITE(PU,*)TRIM(Mssg)
             ENDDO
          ELSE
             !  Default is the full dump of the internal representation
             CALL PrintProtectL(PU)
             WRITE(PU,3)
             WRITE(PU,*)'Internal representation of the geometry:'
             IF(GM%InAU)THEN
                WRITE(PU,*)'Geometry originally in AU has not been rescaled.'
             ELSE
                WRITE(PU,*)'Geometry originally in Angstroms has been converted to AU.'
             ENDIF
             IF(GM%Ordrd==SFC_HILBERT)THEN
                WRITE(PU,*)'Geometry has been reordered using the Hilbert curve.'
             ELSEIF(GM%Ordrd==SFC_PEANO)THEN
                WRITE(PU,*)'Geometry has been reordered using the Peano curve.'
             ELSEIF(GM%Ordrd==SFC_RANDOM)THEN
                WRITE(PU,*)'Geometry has been randomly reordered.'
             ELSE
                WRITE(PU,*)'Geometry has not been reordered.'
             ENDIF
             Mssg='Number of electrons = '//TRIM(IntToChar(GM%NElec))
             WRITE(PU,*)TRIM(Mssg)
             WRITE(PU,2)
             CALL ClosePU(PU)
             CALL Print_PBCInfo(GM%PBC,FileName_O,Unit_O)
             PU=OpenPU(FileName_O=FileName_O,Unit_O=Unit_O)
             WRITE(PU,*)' Cartesian coordinates in AU:'
             DO I=1,GM%NAtms
                Mssg=TRIM(IntToChar(I))//'   '//Ats(INT(GM%AtNum%D(I))) &  !!!! correct only for integer charged QM atoms
                     //'   '//DblToMedmChar(GM%Carts%D(1,I))          &
                     //'   '//DblToMedmChar(GM%Carts%D(2,I))          &
                     //'   '//DblToMedmChar(GM%Carts%D(3,I))
                WRITE(PU,*)TRIM(Mssg)
             ENDDO
             WRITE(PU,3)
             CALL PrintProtectR(PU)
          ENDIF
          CALL ClosePU(PU)
#ifdef PARALLEL
       ENDIF
#endif
2      FORMAT(72('-'))
3      FORMAT(72('='))
     END SUBROUTINE Print_CRDS
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!    Print a BCSR matrix
!
     SUBROUTINE Print_BCSR(A,Name,FileName_O,Unit_O)
       TYPE(BCSR)                           :: A
       TYPE(DBL_RNK2)                       :: B
       CHARACTER(LEN=*),INTENT(IN)          :: Name
       INTEGER,         OPTIONAL,INTENT(IN) :: Unit_O
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: FileName_O

       IF(PrintFlags%Mat/=DEBUG_MATRICES.AND.(.NOT.PRESENT(Unit_O)))RETURN
#ifdef PARALLEL
       IF(MyId==ROOT)THEN
#endif
         CALL SetEq(B,A)
         CALL Print_DBL_RNK2(B,Name,FileName_O,Unit_O)
         CALL Delete(B)
#ifdef PARALLEL
       ENDIF
#endif
     END SUBROUTINE Print_BCSR
#ifdef PARALLEL
!-----------------------------------------------------------------------------
!    Print a DBCSR matrix
!
     SUBROUTINE Print_DBCSR(A,Name,FileName_O,Node_O,Distrib_O)
        TYPE(DBCSR), INTENT(INOUT)           :: A
        TYPE(DBL_RNK2)                       :: B
        TYPE(BCSR)                           :: C
        CHARACTER(LEN=*),INTENT(IN)          :: Name
        CHARACTER(LEN=DEFAULT_CHR_LEN)       :: Name2
        CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: FileName_O
        INTEGER,OPTIONAL                     :: Node_O
        LOGICAL,OPTIONAL                     :: Distrib_O
        INTEGER                              :: I
        IF(PrintFlags%Mat/=DEBUG_MATRICES)RETURN
        IF(PRESENT(Distrib_O))THEN
           IF(Distrib_O)THEN
              CALL SetEq(B,A)
              DO I=0,NPrc-1
                 IF(InParallel)CALL AlignNodes()
                 IF(MyId==I)THEN
                    Name2=TRIM(Name)//'['//TRIM(IntToChar(I))//']'
                    CALL Print_DBL_RNK2(B,Name2,OutFile)
                 ENDIF
              ENDDO
              CALL Delete(B)
           ELSE
              CALL Halt(' Logic error 1 in Print_DBCSR ')
           ENDIF
        ELSEIF(PRESENT(Node_O))THEN
           IF(MyId==Node_O)THEN
              CALL SetEq(B,A)
              Name2=TRIM(Name)//'['//TRIM(IntToChar(Node_O))//']'
              CALL Print_DBL_RNK2(B,Name2,OutFile)
              CALL Delete(B)
           ENDIF
        ELSE
           CALL SetEq(C,A)
           CALL Print_BCSR(C,Name,Filename_O)
           CALL Delete(C)
        ENDIF
     END SUBROUTINE Print_DBCSR
#endif

     SUBROUTINE Print_DBL_RNK2(A,Name,FileName_O,Unit_O)
        TYPE(DBL_RNK2)             :: A
        CHARACTER(LEN=*)           :: Name
        INTEGER, OPTIONAL          :: Unit_O
        CHARACTER(LEN=*), OPTIONAL :: FileName_O
        IF(.NOT.AllocQ(A%Alloc)) &
           CALL Halt(' Matrix not allocated in Print_DBL_RNK2.')
        CALL Print_DBL_Rank2A(A%D,Name,FileName_O,Unit_O)
     END SUBROUTINE Print_DBL_RNK2

     SUBROUTINE Print_DBL_Rank2A(A,Name,FileName_O,Unit_O)
        REAL(DOUBLE),DIMENSION(:,:),INTENT(IN) :: A
        INTEGER                                :: I,J,K,L,M,N,Unit
        CHARACTER(LEN=*)                       :: Name
        INTEGER, OPTIONAL                      :: Unit_O
        CHARACTER(LEN=*),OPTIONAL,INTENT(IN)   :: FileName_O
        Unit=Out; IF(PRESENT(Unit_O))Unit=Unit_O
        IF(PRESENT(FileName_O).AND.Unit/=6)THEN
           CALL OpenASCII(FileName_O,Unit)
        ELSEIF(Unit/=6)THEN
           CALL OpenASCII(OutFile,Unit)
        ENDIF
        M=SIZE(A,1); N=SIZE(A,2)
        IF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
           IF(M/=N) &
              CALL Halt(' Print_DBL_Rank2A does not currently'//   &
                        ' support rectangular matrices with'//     &
                        ' PrintFlags%Fmt==DEBUG_MMASTYLE:'//Rtrn// &
                        ' attempting to print '//TRIM(NAME))
           WRITE(Unit,100)
           WRITE(Unit,101)
           WRITE(Unit,102)N
           DO I=1,N,2
              K=MIN(I+1,N)
              IF(K-I+1.EQ.2)THEN
                 WRITE(Unit,202)(J,J=I,K)
              ELSEIF(K-I+1.EQ.1)THEN
                 WRITE(Unit,201)(J,J=I,K)
              ENDIF
              DO L=1,N
                  IF(K-I+1.EQ.2)THEN
                     WRITE(Unit,302)L,L, &
                       (FRACTION(A(L,J)),EXPONENT(A(L,J)),J=I,K)
                  ELSEIF(K-I+1.EQ.1)THEN
                     WRITE(Unit,301)L,L, &
                       (FRACTION(A(L,J)),EXPONENT(A(L,J)),J=I,K)
                  ENDIF
              ENDDO
           ENDDO
           WRITE(Unit,103)N
           WRITE(Unit,104)TRIM(Name)
           WRITE(Unit,99)
        ELSEIF(PrintFlags%Fmt == DEBUG_MMSTYLE) THEN
           WRITE(Unit, 401) TRIM(Name)
           WRITE(Unit, "(A)") "%%MatrixMarket matrix coordinate real general"
           K = 0
           DO I = 1, N
              DO J = 1, M
                 IF(A(I, J) /= 0.0D0) THEN
                   K = K+1
                 ENDIF
              ENDDO
           ENDDO
           WRITE(Unit, "(I5,x,I5,x,I10)") M, N, K
           DO I = 1, N
              DO J = 1, M
                 IF(A(I, J) /= 0.0D0) THEN
                    WRITE(Unit, "(I5,1x,I5,1x,F16.8)") I, J, A(I, J)
                 ENDIF
              ENDDO
           ENDDO
        ELSE
           WRITE(Unit,401)TRIM(Name)

           DO I=1,N,10
              K=MIN(I+9,N)
              WRITE(Unit,501)(J,J=I,K)
              DO L=1,M
                WRITE(Unit,701) L,(A(L,J),J=I,K)
              ENDDO
           ENDDO

        ENDIF
        IF(Out/=6) CLOSE(Out)
        RETURN

   100  FORMAT(' (*',65('='),'*)')
    99  FORMAT(' mv[m_List]:=MatrixForm[ N[ Chop[m,0.001],4]]; ')
   101  FORMAT(1x,'Ap[x_List,y_List]:=Append[x,y];')
   102  FORMAT(1x,'Tmp=Table[{},{i,',I4,'}];')
   103  FORMAT(1x,'Do[Tmp[[i]]=Flatten[Tmp[[i]]],{i,1,',I4,'}];')
   104  FORMAT(1x,A,' = Tmp ; ')
   201  FORMAT('(*           ',1(14x,i4),'   *)')
   202  FORMAT('(*           ',2(14x,i4),'   *)')
   301  FORMAT(1x,'Tmp[[',I4,']]=Ap[Tmp[[',I4,']],', &
                         '{',F19.16,'*2^(',I4,')}];')
   302  FORMAT(1x,'Tmp[[',I4,']]=Ap[Tmp[[',I4,']],', &
                       '{',1(F19.16,'*2^(',I4,'),'), &
                             F19.16,'*2^(',I4,')}];')
   401  FORMAT(2x,A,':')
   501  FORMAT(T2,10I16)
   601  FORMAT(I5,10F10.5)
   701  FORMAT(I5,10D16.8)
!
   END SUBROUTINE Print_DBL_Rank2A
!==================================================================
!
!    Print Check Sums
!
!==================================================================
   SUBROUTINE Print_CheckSum_DBL_VECT(A,Name,Unit_O,Proc_O)
        TYPE(DBL_VECT), INTENT(IN)           :: A
        REAL(DOUBLE)                         :: Chk
        CHARACTER(LEN=*)                     :: Name
        INTEGER,         OPTIONAL,INTENT(IN) :: Unit_O
        CHARACTER(LEN=*),OPTIONAL            :: Proc_O
        INTEGER                              :: I,PU
        CHARACTER(LEN=DEFAULT_CHR_LEN)       :: ChkStr
!----------------------------------------------------------------------------------------
        IF(PrintFlags%Key/=DEBUG_MAXIMUM.AND. &
           PrintFlags%Chk/=DEBUG_CHKSUMS)RETURN
!---------------------------------------------------------------------------------------
!       Compute check sum
        Chk=Zero
        DO I=LBOUND(A%D,1),UBOUND(A%D,1)
           Chk=Chk+A%D(I)*A%D(I)
        ENDDO
        Chk=SQRT(Chk)
#ifdef PARALLEL
        IF(MyID==ROOT)THEN
#endif
           ChkStr=CheckSumString(Chk,Name,Proc_O)
!          Write check string
           PU=OpenPU(Unit_O=Unit_O)
           WRITE(PU,'(A)')TRIM(ChkStr)
           CALL ClosePU(PU)
#ifdef PARALLEL
        ENDIF
#endif
      END SUBROUTINE Print_CheckSum_DBL_VECT
!----------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
   SUBROUTINE Print_CheckSum_DBL_RNK2(A,Name,Unit_O,Proc_O)
        TYPE(DBL_RNK2), INTENT(IN)           :: A
        REAL(DOUBLE)                         :: Chk
        CHARACTER(LEN=*)                     :: Name
        INTEGER,         OPTIONAL,INTENT(IN) :: Unit_O
        CHARACTER(LEN=*),OPTIONAL            :: Proc_O
        INTEGER                              :: I,J,PU
        CHARACTER(LEN=DEFAULT_CHR_LEN)       :: ChkStr
!----------------------------------------------------------------------------------------
        IF(PrintFlags%Key/=DEBUG_MAXIMUM.AND. &
           PrintFlags%Chk/=DEBUG_CHKSUMS)RETURN
!---------------------------------------------------------------------------------------
!       Compute check sum
        Chk=Zero
        DO I=LBOUND(A%D,1),UBOUND(A%D,1)
           DO J=LBOUND(A%D,2),UBOUND(A%D,2)
              Chk=Chk+A%D(I,J)*A%D(I,J)
           ENDDO
        ENDDO
        Chk=SQRT(Chk)
#ifdef PARALLEL
        IF(MyID==ROOT)THEN
#endif
           ChkStr=CheckSumString(Chk,Name,Proc_O)
!          Write check string
           PU=OpenPU(Unit_O=Unit_O)
           WRITE(PU,'(A)')TRIM(ChkStr)
           CALL ClosePU(PU)
#ifdef PARALLEL
        ENDIF
#endif
      END SUBROUTINE Print_CheckSum_DBL_RNK2

      FUNCTION CheckSum_BCSR(A) RESULT(Chk)
        TYPE(BCSR), INTENT(IN) :: A
        REAL(DOUBLE)           :: Chk
        INTEGER                :: I

        Chk=Zero
        DO I=1,A%NNon0
           Chk=Chk+A%MTrix%D(I)*A%Mtrix%D(I)
        ENDDO
        Chk=SQRT(Chk)
      END FUNCTION CheckSum_BCSR

      SUBROUTINE Print_CheckSum_BCSR(A,Name,Proc_O,Unit_O)
        TYPE(BCSR), INTENT(IN)               :: A
        REAL(DOUBLE)                         :: Chk
        CHARACTER(LEN=*)                     :: Name
        INTEGER,         OPTIONAL,INTENT(IN) :: Unit_O
        CHARACTER(LEN=*),OPTIONAL            :: Proc_O
        INTEGER                              :: I,J,PU,AtA,AtB,P,R,RN,NA,K1,K2
        CHARACTER(LEN=DEFAULT_CHR_LEN)       :: ChkStr
        REAL(DOUBLE),DIMENSION(5,5)          :: Temp
!----------------------------------------------------------------------------------------
        IF(PrintFlags%Key/=DEBUG_MAXIMUM .AND. &
             PrintFlags%Chk/=DEBUG_CHKSUMS)RETURN
!---------------------------------------------------------------------------------------
!       Compute check sum
#ifdef PARALLEL
        IF(MyID==ROOT)THEN
#endif
           Chk = CheckSum(A)
           ! Create check string
           ChkStr=CheckSumString(Chk,Name,Proc_O)
           ! Write check string
           CALL MondoLog(DEBUG_MAXIMUM, "", ChkStr)
#ifdef PARALLEL
        ENDIF
#endif
      END SUBROUTINE Print_CheckSum_BCSR
!----------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
#ifdef PARALLEL
   SUBROUTINE Print_CheckSum_DBCSR(A,Name,Proc_O,Unit_O)
      TYPE(DBCSR), INTENT(IN)   :: A
      CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proc_O
      INTEGER,         OPTIONAL,INTENT(IN) :: Unit_O
      REAL(DOUBLE)              :: Chk
      REAL(DOUBLE), EXTERNAL    :: DDot
      REAL(DOUBLE)              :: DotPrd
      CHARACTER(LEN=*)          :: Name
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: ChkStr
      INTEGER :: I,PU
!----------------------------------------------------------------------------------------
    IF(PrintFlags%Key/=DEBUG_MAXIMUM.AND. &
       PrintFlags%Chk/=DEBUG_CHKSUMS)RETURN
!---------------------------------------------------------------------------------------
      Chk=Zero
      DO I=1,A%NNon0
         Chk=Chk+A%MTrix%D(I)*A%Mtrix%D(I)
      ENDDO
      DotPrd=Reduce(Chk)
      IF(MyID==ROOT)THEN
         Chk=SQRT(DotPrd)
         ! Create check string
         ChkStr=CheckSumString(Chk,Name,Proc_O)
         PU=OpenPU(Unit_O=Unit_O)
         WRITE(PU,'(A)')ChkStr
         CALL ClosePU(PU)
      ENDIF
   END SUBROUTINE Print_CheckSum_DBCSR
#endif
!----------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
  SUBROUTINE Print_CheckSum_HGRho(A,Name,DistRho,Proc_O,Unit_O)
    TYPE(HGRho)                          :: A
    CHARACTER(LEN=*)                     :: Name
    INTEGER,         OPTIONAL,INTENT(IN) :: Unit_O
    CHARACTER(LEN=*),OPTIONAL            :: Proc_O
    LOGICAL                              :: DistRho
    INTEGER                              :: PU,I,L,M,N,LMN,jadd,zq,iq,oq,orr,Ell,LenKet,NQ
    REAL(DOUBLE)                         :: ChkCoef,ChkExp,ChkLMN,Zt
    CHARACTER(LEN=2*DEFAULT_CHR_LEN)     :: ChkStr
#ifdef PARALLEL
    REAL(DOUBLE)                         :: SChkCoef,SChkExp,SChkLMN
#endif
!----------------------------------------------------------------------------------------
    IF(PrintFlags%Key/=DEBUG_MAXIMUM.AND. &
       PrintFlags%Chk/=DEBUG_CHKSUMS)RETURN
!---------------------------------------------------------------------------------------
    IF(.NOT. AllocQ(A%Alloc)) THEN
       CALL Halt(' Density not allocated in Print_CheckSum_HGRho')
    ENDIF
!
    ChkCoef=Zero
    ChkExp =Zero
    ChkLMN =Zero
    DO zq=1,A%NExpt-1
       Zt =A%Expt%D(zq)
       oq =A%OffQ%I(zq)
       orr=A%OffR%I(zq)
       Ell    = A%Lndx%I(zq)
       LenKet = LHGTF(Ell)
       DO iq=1,A%NQ%I(zq)
          jadd=orr+(iq-1)*LenKet
          DO LMN=1,LenKet
             ChkCoef=ChkCoef+ ABS(A%Co%D(jadd+LMN))
             ChkExp =ChkExp + Zt*ABS(A%Co%D(jadd+LMN))
             ChkLMN =ChkLMN + DBLE(LMN)*ABS(A%Co%D(jadd+LMN))
          ENDDO
       ENDDO
    ENDDO
#ifdef PARALLEL
    SChkCoef = Reduce(ChkCoef)
    SChkExp  = Reduce(ChkExp)
    SChkLMN  = Reduce(ChkLMN)
    IF(MyID == ROOT) THEN
       IF(DistRho) THEN
          ChkCoef = SChkCoef
          ChkExp  = SChkExp
          ChkLMN  = SChkLMN
       ENDIF
#endif
!      Create check string
       ChkStr=CheckSumString(ChkCoef,Name//TRIM("-Coef"),Proc_O)
       PU=OpenPU(Unit_O=Unit_O)
       WRITE(PU,'(A)') TRIM(ChkStr)
       CALL ClosePU(PU)
       ChkStr=CheckSumString(ChkExp,Name//TRIM("-Expt"),Proc_O)
       PU=OpenPU(Unit_O=Unit_O)
       WRITE(PU,'(A)') TRIM(ChkStr)
       CALL ClosePU(PU)
       ChkStr=CheckSumString(ChkLMN,Name//TRIM("-LMNs"),Proc_O)
       PU=OpenPU(Unit_O=Unit_O)
       WRITE(PU,'(A)') TRIM(ChkStr)
       CALL ClosePU(PU)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Print_CheckSum_HGRho
!----------------------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------------------
  FUNCTION CheckSumString(Chk,Name,Proc_O) RESULT(ChkStr)
    REAL(DOUBLE)                         :: Chk
    CHARACTER(LEN=*)                     :: Name
    CHARACTER(LEN=*),OPTIONAL            :: Proc_O
    CHARACTER(LEN=DEFAULT_CHR_LEN)       :: ChkStr
    IF(PRESENT(Proc_O).AND.PrintFlags%Fmt/=DEBUG_MMASTYLE)THEN
       ChkStr=ProcessName(Proc_O)//TRIM(Name)//CheckEq//TRIM(DblToChar(Chk))
    ELSEIF(PrintFlags%Fmt/=DEBUG_MMASTYLE)THEN
       ChkStr=TRIM(Name)//CheckEq//TRIM(DblToChar(Chk))
    ELSEIF(PRESENT(Proc_O).AND.PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
       ChkStr='(* '//TRIM(Proc_O)//' *)'//'ChkSum'//TRIM(Name)       &
            //' = '//TRIM(FltToChar(FRACTION(Chk))) &
            //'*2^('//TRIM(IntToChar(EXPONENT(Chk)))//');'
    ELSEIF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
       ChkStr='ChkSum'//TRIM(Name)//' = '//TRIM(FltToChar(FRACTION(Chk))) &
            //'*2^('//TRIM(IntToChar(EXPONENT(Chk)))//');'
    ELSE
       CALL Halt('Logic error in Print_CheckSum_BCSR')
    ENDIF
  END FUNCTION CheckSumString
!========================================================================================
!     Print Out the Density
!========================================================================================
  SUBROUTINE Print_HGRho(A,Name,FileName_O,Unit_O)
    TYPE(HGRho)                      :: A
    CHARACTER(LEN=*)                 :: Name
    CHARACTER(LEN=*),OPTIONAL        :: FileName_O
    INTEGER,OPTIONAL                 :: Unit_O
    INTEGER                          :: PU
    INTEGER                          :: L,M,N,LMN,NPrim,iadd,jadd,zq,iq,oq,orr,Ell,LenKet
    CHARACTER(LEN=DEFAULT_CHR_LEN)   :: Strng

    IF(.NOT. AllocQ(A%Alloc)) THEN
       CALL Halt('Density not allocated in PPrint_HGRho')
    ENDIF

    IF(PRESENT(Unit_O) .AND. PRESENT(FileName_O)) THEN
      PU = OpenPU(Unit_O = Unit_O, FileName_O = FileName_O)
    ELSE IF(PRESENT(Unit_O)) THEN
      PU = OpenPU(Unit_O = Unit_O)
    ELSE IF(PRESENT(FileName_O)) THEN
      PU = OpenPU(FileName_O = FileName_O)
    ELSE
      PU = OpenPU()
    ENDIF

    IF(PrintFlags%Fmt==DEBUG_MMASTYLE) THEN
       NPrim=0
       WRITE(PU,*)'ClearAll[Rho];'
       WRITE(PU,*)'Rho[R_List]:=Module[{zeta,zs,Q,RQ,ExpRQ,Lx,Ly,Lz,RhoSum},'
       WRITE(PU,*)'RhoSum=0;'
       DO zq=1,A%NExpt-1
          oq =A%OffQ%I(zq)
          orr=A%OffR%I(zq)
          Ell=A%Lndx%I(zq)
          LenKet=LHGTF(Ell)
          Strng=Squish('zeta='//DblToMMAChar(A%Expt%D(zq))//';')
          WRITE(PU,*)TRIM(Strng)
          WRITE(PU,*)'zs=Sqrt[zeta];'
          IF(A%NQ%I(zq).NE.0) THEN
             DO iq=1,A%NQ%I(zq)
                NPrim=NPrim+1
                iadd=oq+iq
                jadd=orr+(iq-1)*LenKet
                Strng=Squish('Q={'//DblToMMAChar(A%Qx%D(iadd))   &
                             //','//DblToMMAChar(A%Qy%D(iadd))   &
                             //','//DblToMMAChar(A%Qz%D(iadd))//'};')
                WRITE(PU,*)TRIM(Strng)
                WRITE(PU,*)'RQ=R-Q;'
                WRITE(PU,*)'ExpRQ=Exp[-zeta*RQ.RQ];'
                Strng='Do[Lx[l]=zs^l*HermiteH[l,zs*RQ[[1]]];'//Rtrn &
                   //'    Ly[l]=zs^l*HermiteH[l,zs*RQ[[2]]];'//Rtrn &
                   //'    Lz[l]=zs^l*HermiteH[l,zs*RQ[[3]]];'//Rtrn &
                   //',{l,0,'//TRIM(IntToChar(Ell))//'}]; '
                WRITE(PU,*)TRIM(Strng)
                DO L=0,Ell
                DO M=0,Ell-L
                DO N=0,Ell-L-M
                   LMN=LMNDex(L,M,N)+jadd
                   Strng=Squish('RhoSum=RhoSum+Lx['//IntToChar(L)//']*Ly['//IntToChar(M) &
                              //']*Lz['//IntToChar(N)//']'//'*ExpRQ*'//DblToMMAChar(A%Co%D(LMN))//';')
                   WRITE(PU,*)TRIM(Strng)
                ENDDO
                ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDDO
       WRITE(PU,*)'RhoSum];'
    ELSEIF(PrintFlags%Key==DEBUG_MEDIUM .OR. PrintFlags%Key==DEBUG_MAXIMUM .OR. PU==6) THEN
       NPrim=0
       WRITE(PU,30) Name
       WRITE(PU,31)
       WRITE(PU,32) A%NExpt
       WRITE(PU,33) A%NDist
       WRITE(PU,34) A%NCoef
       WRITE(PU,31)
       IF(PrintFlags%Key==DEBUG_MAXIMUM .AND. A%NDist .LT. 100) THEN
          DO zq=1,A%NExpt
             oq =A%OffQ%I(zq)
             orr=A%OffR%I(zq)
             Ell    = A%Lndx%I(zq)
             LenKet = LHGTF(Ell)
             IF(A%NQ%I(zq).NE.0) THEN
                DO iq=1,A%NQ%I(zq)
                   NPrim=NPrim+1
                   iadd=oq+iq
                   jadd=orr+(iq-1)*LenKet
                   WRITE(PU,20) NPrim,Ell,A%Expt%D(zq)
                   WRITE(PU,21) A%Qx%D(iadd),A%Qy%D(iadd),A%Qz%D(iadd)
                   WRITE(PU,22)
                   DO L=0,Ell
                      DO M=0,Ell-L
                         DO N=0,Ell-L-M
                            LMN=LMNDex(L,M,N)+jadd
                            WRITE(PU,24) L,M,N
                            WRITE(PU,25) A%Co%D(LMN)
                         ENDDO
                      ENDDO
                   ENDDO
                   WRITE(PU,31)
                ENDDO
             ENDIF
          ENDDO
       ENDIF
    ENDIF
!
    CALL ClosePU(PU)
!
    RETURN
!
10  FORMAT(' zeta[',I4,']=SetPrecision[',F19.16,'*2^(' ,I4, ')},50];')
11  FORMAT(' Ell[',I4,']=',I2,';')
12  FORMAT(' Q[',I4,']=SetPrecision[{',F19.16,'*2^(' ,I4, '),', &
                                       F19.16,'*2^(' ,I4, '),', &
                                       F19.16,'*2^(' ,I4, ')},50];')
14  FORMAT(' Co[',I4,', ',I2,', ',I2,', ',I2,']=SetPrecision[',F19.16,'*2^(',I4,'),50];')
15  FORMAT(' NPrim=',I6)
!
20  FORMAT(2x,'Primative #',I4,2x,' Max L = ',I4,2x,' EXPONENT = ',D14.8)
21  FORMAT(2x,'QR     = (',D14.8,', ',D14.8,', ',D14.8,')')
22  FORMAT(2x,'=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=')
24  FORMAT(3x,'L=',I4,' M=',I4,' N=',I4)
25  FORMAT(3x,'RhoCo(L,M,N) = ',D14.8)
30  FORMAT(1x,A,' in a Hermite Gaussian basis: ')
31  FORMAT(1x,'=========================================================================')
32  FORMAT(1x,' Number of Exponents    = ',I8,2x)
33  FORMAT(1x,' Number of Distribution = ',I8,2x)
34  FORMAT(1x,' Number of Coeffients   = ',I8,2x)
  END SUBROUTINE Print_HGRho
!========================================================================================
!     Print Out the Cartesian Multipoles
!========================================================================================
  SUBROUTINE Print_CMPoles(A,Proc_O,Unit_O)
    TYPE(CMPoles)                    :: A
    CHARACTER(LEN=*),OPTIONAL        :: Proc_O
    INTEGER,OPTIONAL                 :: Unit_O
    INTEGER                          :: PU
    INTEGER                          :: I,J
    REAL(DOUBLE)                     :: TrQ
    CHARACTER(LEN=DEFAULT_CHR_LEN)   :: Mssg
!-------------------------------------------------------------------------------------
    IF(.NOT. AllocQ(A%Alloc)) THEN
       CALL Halt(' Multipoles not allocated in PPrint_CMPoles')
    ENDIF
    TrQ=A%QPole%D(1)+A%QPole%D(2)+A%QPole%D(3)
    Mssg=ProcessName(Proc_O,'Moments')                  &
        //'<r> = ('//TRIM(DblToShrtChar(A%DPole%D(1)))  &
        //', '//TRIM(DblToShrtChar(A%DPole%D(2)))       &
        //', '//TRIM(DblToShrtChar(A%DPole%D(3)))       &
        //'), <r^2> = '//TRIM(DblToShrtChar(TrQ))
    PU=OpenPU(Unit_O=Unit_O)
    WRITE(PU,*)TRIM(Mssg)
    CALL ClosePU(PU)
!
#ifdef CJsStuff
    WRITE(PU,10)
    WRITE(PU,11) Strng
    WRITE(PU,12)
    WRITE(PU,20) A%DPole%D(1),A%DPole%D(2),A%DPole%D(3)
    WRITE(PU,13)
    WRITE(PU,20) A%QPole%D(1),A%QPole%D(2),A%QPole%D(3)
    WRITE(PU,14)
    WRITE(PU,20) A%QPole%D(4),A%QPole%D(5),A%QPole%D(6)
    WRITE(PU,10)
10  FORMAT(54('='))
11  FORMAT(1x,A54)
12  FORMAT(6x,' p(x) ',12x,' p(y) ',12x,' p(z) ')
13  FORMAT(6x,'q(x*x)',12x,'q(y*y)',12x,'q(z*z)')
14  FORMAT(6x,'q(x*y)',12x,'q(x*z)',12x,'q(y*z)')
20  FORMAT(1x,D15.8,3x,D15.8,3x,D15.8)
#endif
!
  END SUBROUTINE Print_CMPoles
!========================================================================================
! Print The Forces
!========================================================================================
  SUBROUTINE Print_Force(GM,Frc,Name_O,Unit_O)
    TYPE(CRDS)                     :: GM
    TYPE(DBL_VECT)                 :: Frc
    CHARACTER(LEN=*),OPTIONAL      :: Name_O
    INTEGER,OPTIONAL               :: Unit_O
    INTEGER                        :: AtA,A1,PU,I,J

    IF(PrintFlags%Key /= DEBUG_MAXIMUM .OR. PrintFlags%MM /= DEBUG_FRC) RETURN
#ifdef PARALLEL
    IF(MyID /= ROOT) RETURN
#endif
    IF(PRESENT(Unit_O)) THEN
      CALL MondoLog(DEBUG_NONE, "Print_Force", "Unit_O = "//TRIM(IntToChar(Unit_O)))
      CALL Halt("[FIXME]")
    ENDIF

    IF(PRESENT(Name_O)) THEN
      CALL MondoLog(DEBUG_NONE, "Force", Name_O)
    ENDIF

    CALL MondoLog(DEBUG_NONE, "Force", "Atom  Z  Forces (eV/A)")
    DO AtA = 1,GM%Natms
       A1=3*(AtA-1)+1
       CALL MondoLog(DEBUG_NONE, "Force", TRIM(IntToChar(AtA))//" "// &
         TRIM(IntToChar(INT(GM%AtNum%D(AtA))))//" "// &
         TRIM(DblToChar((au2eV/AUToAngstroms)*Frc%D(A1)))//" "// &
         TRIM(DblToChar((au2eV/AUToAngstroms)*Frc%D(A1+1)))//" "// &
         TRIM(DblToChar((au2eV/AUToAngstroms)*Frc%D(A1+2))))
    ENDDO

  END SUBROUTINE Print_Force

!========================================================================================
! Print The Lattice Forces
!========================================================================================
  SUBROUTINE Print_LatForce(GM,LFrc,Name_O,Unit_O)
    TYPE(CRDS)                     :: GM
    REAL(DOUBLE),DIMENSION(3,3)    :: LFrc
    CHARACTER(LEN=*),OPTIONAL      :: Name_O
    INTEGER,OPTIONAL               :: Unit_O
    INTEGER                        :: AtA,A1,A2,PU,I,J
!----------------------------------------------------------------------------------------
    IF(PrintFlags%Key /= DEBUG_MAXIMUM .OR. PrintFlags%MM /= DEBUG_FRC) RETURN
    IF(GM%PBC%Dimen == 0) RETURN
#ifdef PARALLEL
    IF(MyID /= ROOT) RETURN
#endif
    PU=OpenPU(Unit_O=Unit_O)
    IF(PRESENT(Name_O)) WRITE(PU,32) TRIM(Name_O)
    WRITE(PU,33)
    IF(GM%PBC%Dimen > 0) THEN
       WRITE(PU,32) TRIM('        div(a)            div(b)           div(c)')
       WRITE(PU,33)
       DO I=1,3
          WRITE(PU,40) LFrc(I,1:3)
       ENDDO
       WRITE(PU,33)
    ENDIF
    CALL ClosePU(PU)
!
32  FORMAT(A)
33  FORMAT(58('-'))
35  FORMAT(I6,'  ',I3,'  ',3(F14.10))
40  FORMAT(3(1X,E18.10))
!
  END SUBROUTINE Print_LatForce
#ifdef EXTREME_DEBUG
!========================================================================================
! Print The Gradients
!========================================================================================
  SUBROUTINE Print_CheckSum_Force(Frc,Name,Proc_O,Unit_O)
    TYPE(DBL_VECT)                 :: Frc
    CHARACTER(LEN=*)               :: Name
    CHARACTER(LEN=*),OPTIONAL      :: Proc_O
    INTEGER,OPTIONAL               :: Unit_O
    INTEGER                        :: AtA,A1,PU
    REAL(DOUBLE)                   :: FX,FY,FZ,ChkF
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: ChkStr
!
    IF(PrintFlags%Key/=DEBUG_MAXIMUM .AND. PrintFlags%Chk/=DEBUG_CHKSUMS) RETURN
!
    FX = Zero
    FY = Zero
    FZ = Zero
    DO AtA = 1,NAtoms
       A1 = 3*(AtA-1)+1
       FX = FX+Frc%D(A1)
       FY = FY+Frc%D(A1+1)
       FZ = FZ+Frc%D(A1+2)
    ENDDO
    ChkF = SQRT(FX*FX+FY*FY+FZ*FZ)
!
    IF(PRESENT(Proc_O)) THEN
       ChkStr=ProcessName(Proc_O)//TRIM(Name)// " ForceSum = " //TRIM(DblToChar(ChkF))
    ELSE
       ChkStr=TRIM(Name)// " ForceSum = " //TRIM(DblToChar(ChkF))
    ENDIF
!
    PU=OpenPU(Unit_O=Unit_O)
    WRITE(PU,'(A)') ChkStr
    CALL ClosePU(PU)
!
  END SUBROUTINE Print_CheckSum_Force
#endif
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
  FUNCTION BlockIndexString(Name,I,K,J,Id)
    INTEGER                         :: I,K
    INTEGER, OPTIONAL               :: J,Id
    CHARACTER(LEN=*)                :: Name
    CHARACTER(LEN=INTERNAL_INT_LEN) :: cI,cJ,cK,cId
    CHARACTER(LEN=64)               :: BlockIndexString
    cI=IntToChar(I)
    cK=IntToChar(K)
    IF(PRESENT(J).AND.PRESENT(Id))THEN
       cJ=IntToChar(J)
       cId=IntToChar(Id)
       BlockIndexString=TRIM(Name)//'(Id='//TRIM(cId)//',I='// &
            TRIM(cI)//',J='//TRIM(cK)//',K='//TRIM(cJ)//')'
    ELSEIF(PRESENT(J))THEN
       cJ=IntToChar(J)
       BlockIndexString=TRIM(Name)//'('// &
            TRIM(cI)//',J='//TRIM(cJ)//','//TRIM(cK)//')'
    ELSEIF(PRESENT(Id))THEN
       cId=IntToChar(Id)
       BlockIndexString=TRIM(Name)//'(Id='//TRIM(cId)//','// &
            TRIM(cI)//','//TRIM(cK)//')'
    ELSE
       BlockIndexString=TRIM(Name)//'('//TRIM(cI)//','//TRIM(cK)//')'
    ENDIF
  END FUNCTION BlockIndexString
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
  SUBROUTINE Elapsed_TIME(T,Init_O,Proc_O)
    TYPE(TIME),           INTENT(INOUT)  :: T
    REAL(DOUBLE)                         :: WallTm,CPUSTm,TimeTot,TimeMax, &
         TimeMin,TimeAve,Imb
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Init_O,Proc_O
    CHARACTER(LEN=2*DEFAULT_CHR_LEN)     :: Mssg
    CHARACTER(LEN=10)                    :: Mssg10
    INTEGER                              :: I,PU
    IF(PRESENT(Init_O))THEN
       IF(Init_O=='Init')THEN
          T%FLOP=0
          T%Wall=Zero
          T%CPUS=Zero
          T%WStrt=WallSec()
          T%CStrt=CPUSec()
       ELSEIF(Init_O=='Start')THEN
          T%WStrt=WallSec()
          T%CStrt=CPUSec()
       ELSEIF(Init_O=='Accum')THEN
          WallTm=WallSec()
          CPUSTm=CPUSec()
          T%Wall=T%Wall+(WallTm-T%WStrt)
          T%CPUS=T%CPUS+(CPUSTm-T%CStrt)
          T%WStrt=WallTm
          T%CStrt=CPUSTm
       ELSE
          CALL Halt('Unknown option '//TRIM(Init_O)//' in Print_Elapsed_TIME')
       ENDIF
    ENDIF
#ifdef PARALLEL

#ifdef NewPrint
#else
    IF(PRESENT(Proc_O))THEN
       IF(InParallel)THEN
          TimeTot=Reduce(T%Wall,MPI_SUM)
          TimeMax=Reduce(T%Wall,MPI_MAX)
          TimeMin=Reduce(T%Wall,MPI_MIN)
       ENDIF
       IF(MyId==ROOT.AND.InParallel.AND.PrintFlags%Key>=DEBUG_MEDIUM)THEN
          ! Compute relative imbalance
          TimeAve=TimeTot/DBLE(NPrc)
          IF(TimeTot/=Zero)THEN
             Imb=ABS(TimeMax-TimeAve)/TimeAve
          ELSE
             Imb=Zero
          ENDIF
          Mssg=ProcessName(Proc_O)//'Parallel Statistics '
          PU=OpenPU()
          WRITE(PU,*)TRIM(Mssg)
          Mssg=       '                    Imblnce = '//TRIM(DblToShrtChar(Imb))      &
               //Rtrn//'                     MinTime = '//TRIM(DblToShrtChar(TimeMin)) &
               //Rtrn//'                     AvgTime = '//TRIM(DblToShrtChar(TimeAve)) &
               //Rtrn//'                     MaxTime = '//TRIM(DblToShrtChar(TimeMax))
          WRITE(PU,*)TRIM(Mssg)
          CLOSE(Out)
       ENDIF
    ENDIF
#endif

#endif
  END SUBROUTINE Elapsed_TIME
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
  SUBROUTINE Print_TIME(T,Proc_O,FileName_O,Unit_O,BareBones_O)
    TYPE(TIME),               INTENT(INOUT) :: T
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Proc_O,FileName_O
    INTEGER, OPTIONAL                       :: Unit_O
    REAL(DOUBLE)                            :: Elapsed_CPUS,Elapsed_WALL,FLOPS
    LOGICAL, OPTIONAL                       :: BareBones_O
    CHARACTER(LEN=DEFAULT_CHR_LEN)          :: Proc
    CHARACTER(LEN=DEFAULT_CHR_LEN*2)        :: Mssg
    INTEGER                                 :: PU
!------------------------------------------------------------------------------------
#ifdef PARALLEL
    IF(InParallel)  &
         CALL AlignNodes()
    IF(MyId==ROOT)THEN
#endif
       Elapsed_CPUS=T%CPUS
       Elapsed_Wall=T%Wall
       PU=OpenPU(FileName_O,Unit_O)
       CALL PrintProtectL(PU)
#ifdef PARALLEL
    ENDIF
    IF(InParallel)THEN
       FLOPS=Reduce(T%FLOP)
    ELSE
#endif
       FLOPS=T%FLOP
#ifdef PARALLEL
    ENDIF
#endif
!----------------------------------------------------
    IF(PRESENT(Proc_O))THEN
       Proc=Proc_O
    ELSE
       Proc=Blnk
    ENDIF
!-------------------------------------------------------------------------
    IF(PRESENT(BareBones_O))THEN
       IF(.NOT.BareBones_O)  &
            CALL Halt(' Logic error in Print_TIME')
#ifdef PARALLEL
       IF(MyId==ROOT)THEN
          WRITE(PU,10)NPrc,Elapsed_Wall,MFlops(FLOPS,Elapsed_Wall)
10        FORMAT(I4,' ',D12.6,' ',I10)
       ENDIF
#else
       WRITE(PU,10) Proc,Elapsed_CPUS,Elapsed_Wall
10     FORMAT(' ',a22,' :: (CPU,WALL) = ',F16.2,' ',F16.2)
#endif
       CALL PrintProtectR(PU)
       CLOSE(Out)
       RETURN
    ENDIF
!-------------------------------------------------------------------------
#ifdef PARALLEL
    IF(MyID==ROOT)THEN
#endif
       IF(MFlops(FLOPS,Elapsed_Wall)>Zero)THEN
          IF(Elapsed_CPUS/=Zero)THEN
#ifdef PARALLEL
             Mssg=ProcessName(TRIM(Proc))//'CPU (Sec,MFLOPS) = (' &
                  //TRIM(DblToShrtChar(Elapsed_CPUS))//', '     &
                  //TRIM(IntToChar(MFlops(FLOPS,Elapsed_CPUS))) &
                  //'), WALL(Sec,MFLOPS) = ('                  &
                  //TRIM(DblToShrtChar(Elapsed_Wall))//', '     &
                  //TRIM(IntToChar(MFlops(FLOPS,Elapsed_Wall))) &
                  //'), NProc = '//TRIM(IntToChar(NPrc))
#else
             Mssg=ProcessName(TRIM(Proc))//'CPU (Sec,MFLOPS) = (' &
                  //TRIM(DblToShrtChar(Elapsed_CPUS))//', '     &
                  //TRIM(IntToChar(MFlops(FLOPS,Elapsed_CPUS))) &
                  //'), WALL(Sec,MFLOPS) = ('                  &
                  //TRIM(DblToShrtChar(Elapsed_Wall))//', '     &
                  //TRIM(IntToChar(MFlops(FLOPS,Elapsed_Wall))) &
                  //')'
#endif
          ELSE
#ifdef PARALLEL
             Mssg=ProcessName(TRIM(Proc))//'WALL (Sec,MFLOPS) = (' &
                  //TRIM(DblToShrtChar(Elapsed_Wall))//', '      &
                  //TRIM(IntToChar(MFlops(FLOPS,Elapsed_Wall)))  &
                  //'), NProc = '//TRIM(IntToChar(NPrc))
#else
             Mssg=ProcessName(TRIM(Proc))//'WALL (Sec,MFLOPS) = (' &
                  //TRIM(DblToShrtChar(Elapsed_Wall))//', '      &
                  //TRIM(IntToChar(MFlops(FLOPS,Elapsed_Wall)))  &
                  //')'
#endif
          ENDIF
       ELSE
          IF(Elapsed_CPUS>Zero)THEN
#ifdef PARALLEL
             Mssg=ProcessName(TRIM(Proc))//'CPU Sec = '   &
                  //TRIM(DblToShrtChar(Elapsed_CPUS))   &
                  //', WALL (Sec) = '                   &
                  //TRIM(DblToShrtChar(Elapsed_Wall))   &
                  //', NProc = '//TRIM(IntToChar(NPrc))
#else
             Mssg=ProcessName(TRIM(Proc))//'CPU Sec = '  &
                  //TRIM(DblToShrtChar(Elapsed_CPUS))  &
                  //', WALL Sec = '                    &
                  //TRIM(DblToShrtChar(Elapsed_Wall))
#endif
          ELSE
#ifdef PARALLEL
             Mssg=ProcessName(TRIM(Proc))//'WALL Sec = '  &
                  //TRIM(DblToShrtChar(Elapsed_Wall))      &
                  //', NProc = '//TRIM(IntToChar(NPrc))
#else
             Mssg=ProcessName(TRIM(Proc))//'WALL (Sec) = '  &
                  //TRIM(DblToShrtChar(Elapsed_Wall))
#endif
          ENDIF
       ENDIF
       WRITE(PU,"(A)")TRIM(Mssg)
       CALL PrintProtectR(PU)
       CLOSE(Out)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Print_TIME
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
  FUNCTION MFlops(Flops,Sec)
    REAL(DOUBLE), INTENT(IN) :: Flops,Sec
    INTEGER                  :: MFlops
    IF(Sec<=Zero.OR.Flops<=Zero)THEN
       MFlops=0
    ELSE
       MFlops=INT(Flops*1.0D-6/Sec)
    ENDIF
  END FUNCTION MFlops
!---------------------------------------------------------------------
!    PRINT MEMORY STATISTICS
!---------------------------------------------------------------------
  SUBROUTINE Print_MEMS(A,Proc,Unit_O)
    TYPE(MEMS),INTENT(IN)            :: A
    CHARACTER(LEN=*),INTENT(IN)      :: Proc
    INTEGER,OPTIONAL                 :: Unit_O
    INTEGER                          :: I,L,PU
    CHARACTER(LEN=2*DEFAULT_CHR_LEN) :: Mssg
#ifdef PARALLEL
    INTEGER :: MaxAllocs,MaxDeAllocs,IErr,MaxMemTab,MaxMaxMem
    ! This is way to much info in parallel!!
    RETURN
#endif

    IF(PrintFlags%Key/=DEBUG_MAXIMUM)RETURN
    L=LEN(TRIM(Proc))
#ifdef PARALLEL
    IF(InParallel)THEN
       IF(MyId==ROOT)THEN
          PU=OpenPU(Unit_O=Unit_O)
          CALL PrintProtectL(PU)
          CLOSE(PU)
       ENDIF

#ifdef NewPrint
       CALL MPI_Reduce(A%Allocs,MaxAllocs,1,MPI_INTEGER,MPI_MAX,ROOT,MPI_COMM_WORLD,IErr)
       CALL MPI_Reduce(A%DeAllocs,MaxDeAllocs,1,MPI_INTEGER,MPI_MAX,ROOT,MPI_COMM_WORLD,IErr)
       CALL MPI_Reduce(A%MemTab,MaxMemTab,1,MPI_INTEGER,MPI_MAX,ROOT,MPI_COMM_WORLD,IErr)
       CALL MPI_Reduce(A%MaxMem,MaxMaxMem,1,MPI_INTEGER,MPI_MAX,ROOT,MPI_COMM_WORLD,IErr)
       IF(MyID == ROOT) THEN
         Mssg=TRIM(Proc)//' :: MaxAllocs='//TRIM(IntToChar(MaxAllocs)) &
              //', MaxDeAllocs='//TRIM(IntToChar(MaxDeAllocs)) &
             //', MaxMemTab='//TRIM(IntToChar(MaxMemTab)) &
            //', MaxMaxMem='//TRIM(IntToChar(MaxMaxMem))
         PU=OpenPU(Unit_O=Unit_O)
         WRITE(PU,"(A)") TRIM(Mssg)
         CLOSE(PU)
       ENDIF
#else
       DO I=0,NPrc-1
          IF(InParallel)CALL AlignNodes()
          IF(MyId==I)THEN
             Mssg=TRIM(Proc)//'#'//TRIM(IntToChar(I))           &
                  //' :: Allocs='//TRIM(IntToChar(A%Allocs))    &
                  //',  DeAllocs='//TRIM(IntToChar(A%DeAllocs)) &
                  //', '//TRIM(IntToChar(A%MemTab ))            &
                  //' bytes are presently allocated. '          &
                  //'A max of '//TRIM(IntToChar(A%MaxMem))      &
                  //' bytes were allocated.'
             PU=OpenPU(Unit_O=Unit_O)
             WRITE(PU,"(A)")Mssg
             CLOSE(PU)
          ENDIF
       ENDDO
       CALL AlignNodes()
#endif
       IF(MyId==ROOT)THEN
          PU=OpenPU()
          CALL PrintProtectR(PU)
          CLOSE(PU)
       ENDIF
    ELSE
#endif
       PU=OpenPU()
       CALL PrintProtectL(PU)
       Mssg=ProcessName(Proc)                             &
            //'Allocs='//TRIM(IntToChar(A%Allocs))        &
            //',  DeAllocs='//TRIM(IntToChar(A%DeAllocs)) &
            //', '//TRIM(IntToChar(A%MemTab ))            &
            //' bytes are presently allocated. '          &
            //'A max of '//TRIM(IntToChar(A%MaxMem))      &
            //' bytes were allocated.'
       WRITE(PU,"(A)")TRIM(Mssg)
       CALL PrintProtectR(PU)
       CLOSE(PU)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Print_MEMS

!--------------------------------------------------------------------------
! Print the CellSet
!--------------------------------------------------------------------------
  SUBROUTINE PPrint_CellSet(CS,Name,FileName_O,Unit_O)
    TYPE(CellSet)                    :: CS
    CHARACTER(LEN=*)                 :: Name
    CHARACTER(LEN=*),OPTIONAL        :: FileName_O
    INTEGER,OPTIONAL                 :: Unit_O
    INTEGER                          :: OutU
    INTEGER                          :: NC
    REAL(DOUBLE)                     :: RMax,R2

    IF(.NOT. AllocQ(CS%Alloc)) THEN
       CALL Halt(' Cells are  not allocated in PPrint_CellSet')
    ENDIF
    IF(PRESENT(Unit_O)) THEN
       OutU=Unit_O
    ELSE
       OutU=Out
    ENDIF
    IF(PRESENT(FileName_O) .AND. OutU /= 6) THEN
       CALL OpenASCII(FileName_O,OutU)
    ELSEIF(OutU /= 6) THEN
       CALL OpenASCII(OutFile,OutU)
    ENDIF
!
    RMax = Zero
    DO NC = 1,CS%NCells
       R2 = CS%CellCarts%D(1,NC)**2+CS%CellCarts%D(2,NC)**2+CS%CellCarts%D(3,NC)**2
       IF(R2 > RMax) RMax = R2
    ENDDO
    RMax = SQRT(RMax)
!
!!$    IF(PrintFlags%Key==DEBUG_MEDIUM) THEN
!!$       WRITE(OutU,30)
!!$       WRITE(OutU,10) Name
!!$       WRITE(OutU,11) CS%NCells,RMax
!!$       WRITE(OutU,30)
!!$    ELSEIF(PrintFlags%Key==DEBUG_MAXIMUM) THEN
       WRITE(OutU,30)
       WRITE(OutU,10) Name
       WRITE(OutU,11) CS%NCells,RMax
       IF(CS%NCells /= 0) THEN
          WRITE(OutU,20)
          WRITE(OutU,31)
          DO NC = 1,CS%NCells
             WRITE(OutU,21) NC,CS%CellCarts%D(1,NC),CS%CellCarts%D(2,NC),CS%CellCarts%D(3,NC)
          ENDDO
       ENDIF
       WRITE(OutU,30)
!!$    ENDIF
    RETURN
10  FORMAT(1x,A)
11  FORMAT(1x,'Number of Cells    = ',I8,2x,'   Maxium Range   = ',D15.8)
12  FORMAT(1x,'Distance Threshold = ',D10.4,'   Pair Threshold = ',D15.8)
20  FORMAT(1x,'<Coordinates of the Cells>')
21  FORMAT(2x,'Cell #',I3,2x,'Q = (',D15.8,', ',D15.8,', ',D15.8,' )')
30  FORMAT(1x,'=========================================================================')
31  FORMAT(1x,'=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=')
  END SUBROUTINE PPrint_CellSet

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
  SUBROUTINE PrintProtectL(Unit)
    INTEGER :: Unit
    IF(PrintFlags%Fmt==DEBUG_MMASTYLE.AND.Unit/=6) &
         WRITE(Unit,*)LeftParenStar
  END SUBROUTINE PrintProtectL
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
  SUBROUTINE PrintProtectR(Unit)
    INTEGER :: Unit
    IF(PrintFlags%Fmt==DEBUG_MMASTYLE.AND.Unit/=6) &
         WRITE(Unit,*)RightParenStar
  END SUBROUTINE PrintProtectR

    SUBROUTINE PImbalance(A,NPrc,Prog_O)
    TYPE(DBL_VECT) :: A
    INTEGER :: NPrc,I
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Prog_O
    CHARACTER(LEN=100) :: Prog
    REAL(DOUBLE) :: TmMin,TmMax,Imbalance,DevSum

    IF(Present(Prog_O)) THEN
      Prog = Prog_O
    ELSE
      Prog = ''
    ENDIF
    TmMax = -100.0D0
    TmMin = 1.D+10
    DO I = 1, NPrc
      TmMax = Max(TmMax,A%D(I))
      TmMin = Min(TmMin,A%D(I))
    ENDDO
    CALL OpenASCII(OutFile,Out)
    WRITE(*,*) TRIM(Prog)//': TmMax= ',TmMax,', TmMin=',TmMin
    WRITE(Out,*) TRIM(Prog)//': TmMax= ',TmMax,', TmMin=',TmMin

    DevSum = 0.0D0
    DO I = 1, NPrc
      DevSum = DevSum + (TmMax-A%D(I))
    ENDDO
    Imbalance = DevSum/(NPrc*TmMax)
    WRITE(*,*) TRIM(Prog)//': Imbalance = ',Imbalance
    WRITE(Out,*) TRIM(Prog)//': Imbalance = ',Imbalance
    CLOSE(Out,STATUS='KEEP')

  END SUBROUTINE PImbalance

!
!==================================================================
!
!    PLOT A BCSR MATRIX
!
!==================================================================
!      FUNCTION Dot(N,V1,V2)
!         REAL(DOUBLE) :: Dot
!         REAL(DOUBLE), DIMENSION(:)  :: V1,V2
!         INTEGER :: I,N
!         Dot=Zero
!         DO I=1,N
!            Dot=Dot+V1(I)*V2(I)
!         ENDDO
!      END FUNCTION Dot
!
END MODULE
